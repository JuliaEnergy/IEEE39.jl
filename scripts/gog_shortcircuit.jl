#=
Short-circuit with line-tripping simulation on a 2×2 IEEE39 Grid of Grids.

Builds a 2×2 lattice of IEEE39 copies (with mild power-flow perturbation and
random generator omissions), introduces a mid-line short-circuit on the chosen
line, trips the faulted line after 150 ms, and simulates 10 s of post-fault
transient dynamics.

Outputs
  gog_sc_angles.png   – per-copy rotor angles δ (shows inter-area transients)
  gog_sc_voltages.png – per-copy bus-voltage angles ∠V (shows fault propagation
                        across tie lines into neighbouring copies)
  gog_nodes.json / gog_lines.json – full network description for downstream use
=#

using Pkg
Pkg.activate(@__DIR__)

using IEEE39
using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using SciMLBase
using Plots
using Random
using Logging
Logging.disable_logging(Logging.Warn)

## ── Simulation parameters ─────────────────────────────────────────────────────

const N          = 2     # lattice size → N×N = 4 copies
const N_VERT     = 39    # buses per single-copy IEEE39 network
const N_LINES_39 = 46    # branches per single-copy IEEE39 network

fault_line = 1           # line index to fault (within copy 0, bus 1–2)
t_fault    = 0.1         # fault onset [s]
t_trip     = 0.25        # line trip   [s]  (150 ms fault duration)
fault_pos  = 0.5         # position on line (0 = src end, 1 = dst end)
T_sim      = 10.0        # simulation horizon [s]
saveat     = 0.0:0.01:T_sim

## ── Build 2×2 GoG ─────────────────────────────────────────────────────────────

rng = MersenneTwister(42)
println("Building $(N)×$(N) Grid of Grids  (pert=0.05, p_omit=0.1)…")
nw, nws = build_IEEE39_grid_of_grids(N; pert=0.05, p_omit=0.1, rng, max_tries=30)

n_tie      = ne(nw) - N^2 * N_LINES_39
center_sni = (N ÷ 2) * N + (N ÷ 2)          # 0-indexed copy that keeps slack

println("  $(nv(nw)) buses,  $(ne(nw)) branches  ($(n_tie) inter-copy tie lines)")
println("  Centre/slack copy: sni=$(center_sni)  " *
        "(buses $(N_VERT*center_sni+1)–$(N_VERT*(center_sni+1)))")
let ge = get_graphelement(nw[EIndex(fault_line)])
    println("  Fault: sc on line $(fault_line) at t=$(t_fault)s, " *
            "trip at t=$(t_trip)s  (bus $(ge.src) → bus $(ge.dst))")
end

## ── Export network description ────────────────────────────────────────────────

export_nodes_json(nw, nws, joinpath(@__DIR__, "gog_nodes.json"); include_nf=false)
export_lines_json(nw, nws, joinpath(@__DIR__, "gog_lines.json"))
println("  Network JSON → scripts/gog_nodes.json  and  gog_lines.json")

## ── Sparsify Jacobian (significant speedup for large networks) ────────────────

NetworkDynamics.set_jac_prototype!(nw; remove_conditions=true)

## ── Simulate ──────────────────────────────────────────────────────────────────

cb   = sc_and_trip(fault_line, t_fault, t_trip, fault_pos)
prob = ODEProblem(nw, nws, (0.0, T_sim); add_nw_cb=cb)

println("\nSimulating…")
t0  = time()
sol = solve(prob, FBDF();
            saveat,
            abstol=1e-4, reltol=1e-4,
            maxiters=200_000, dtmin=1e-6, force_dtmin=true)
elapsed = round(time() - t0; digits=1)
println("  Retcode: $(sol.retcode)  |  $(length(sol.t)) saved time points  |  $(elapsed)s")

SciMLBase.successful_retcode(sol.retcode) ||
    @warn "Solver did not succeed — plots may be incomplete"

## ── Per-copy index helpers ────────────────────────────────────────────────────

# Global bus index for 0-indexed copy k, local bus b
gv(k, b) = N_VERT * k + b

machine_locals = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39]
load_locals    = [3, 4, 7, 8, 12, 15, 16, 18, 20, 21, 23, 25, 26, 28, 29]

# Human-readable copy labels; centre copy is labelled explicitly
copy_label(k) = "Copy $k" * (k == center_sni ? "  [centre / slack]" : "")

## ── Plot 1 – per-copy machine rotor angles ────────────────────────────────────
# Each subplot shows the generator angles δ in one copy of the IEEE39 grid.
# The fault copy (0) reacts immediately; dynamics in neighbouring copies arrive
# via the tie lines a few cycles later.

p1 = plot(layout=(N^2, 1), size=(960, 200*N^2 + 60),
          left_margin=8Plots.mm, bottom_margin=2Plots.mm)

for k in 0:N^2-1
    idxs = vidxs(nw, [gv(k, b) for b in machine_locals], r"machine₊δ$")
    isempty(idxs) && continue
    plot!(p1, sol; idxs, subplot=k+1, legend=false,
          ylabel="δ [rad]",
          title=copy_label(k), titleloc=:left, titlefont=8,
          xlabel=k == N^2-1 ? "t [s]" : "")
    vline!(p1, [t_fault]; subplot=k+1, color=:red,    lw=1.5, ls=:solid, label="fault")
    vline!(p1, [t_trip];  subplot=k+1, color=:darkorange, lw=1.5, ls=:dash,  label="trip")
end

## ── Plot 2 – per-copy bus voltage angles ──────────────────────────────────────
# Voltage angle perturbations at load buses; a non-zero response in copies 1–3
# confirms that the fault propagates through the tie-line connections.

p2 = plot(layout=(N^2, 1), size=(960, 200*N^2 + 60),
          left_margin=8Plots.mm, bottom_margin=2Plots.mm)

for k in 0:N^2-1
    idxs = vidxs(nw, [gv(k, b) for b in load_locals], :busbar₊u_arg)
    isempty(idxs) && continue
    plot!(p2, sol; idxs, subplot=k+1, legend=false,
          ylabel="∠V [rad]",
          title=copy_label(k), titleloc=:left, titlefont=8,
          xlabel=k == N^2-1 ? "t [s]" : "")
    vline!(p2, [t_fault]; subplot=k+1, color=:red,       lw=1.5, ls=:solid, label="fault")
    vline!(p2, [t_trip];  subplot=k+1, color=:darkorange, lw=1.5, ls=:dash,  label="trip")
end

## ── Save figures ──────────────────────────────────────────────────────────────

out_angles   = joinpath(@__DIR__, "gog_sc_angles.png")
out_voltages = joinpath(@__DIR__, "gog_sc_voltages.png")
savefig(p1, out_angles)
savefig(p2, out_voltages)
println("\nFigures saved:")
println("  $out_angles")
println("  $out_voltages")
