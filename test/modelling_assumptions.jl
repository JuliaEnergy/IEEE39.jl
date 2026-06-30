# ─────────────────────────────────────────────────────────────────────────────
# Tests that document the modelling assumptions of the IEEE39 package.
# Each testset pins a behaviour that the rest of the code silently relies on, so
# that a future change which breaks the assumption fails loudly here.
# (Included from runtests.jl; the `using`s below also let it run standalone.)
# ─────────────────────────────────────────────────────────────────────────────
using Test
using IEEE39
using IEEE39: DATA_DIR, node_nf_linearization
using PowerDynamics
using NetworkDynamics
using OrdinaryDiffEq
using SciMLBase
using Random
using CSV, DataFrames
using LinearAlgebra: norm

@testset "Modelling assumptions" begin

    # ── Parameter assignment ────────────────────────────────────────────────
    # `apply_csv_params!` (and the branch loader) assign CSV columns to model
    # parameters via a `Regex(col * "\$")` match. A column that matches *zero*
    # symbols is silently ignored (its value never reaches the model); one that
    # matches *several* overwrites unrelated parameters. Both are bugs, so pin
    # every data column to exactly one parameter on a representative bus/edge.
    @testset "CSV columns map to exactly one model parameter" begin
        nw = set_IEEE39_PF_init(get_IEEE39_base())
        nmatch(syms, col) = count(s -> occursin(Regex(col * "\$"), string(s)), syms)
        function check(csv, keycols, syms)
            df = CSV.read(joinpath(DATA_DIR, csv), DataFrame)
            for col in names(df)
                col in keycols && continue
                @test nmatch(syms, col) == 1
            end
        end
        check("load.csv",    ["bus"], psym(nw[VIndex(3)]))    # pure load
        check("machine.csv", ["bus"], psym(nw[VIndex(30)]))   # controlled machine
        check("avr.csv",     ["bus"], psym(nw[VIndex(30)]))
        check("gov.csv",     ["bus"], psym(nw[VIndex(30)]))
        check("machine.csv", ["bus"], psym(nw[VIndex(39)]))   # uncontrolled machine+load
        check("branch.csv",  ["src_bus", "dst_bus", "transformer"], psym(nw[EIndex(1)]))
    end

    # ── Normal-form linearization ───────────────────────────────────────────
    # `nf_linearization` assumes every dynamic bus has a 2-in/2-out (current ->
    # voltage) interface and a (near) steady-state operating point. After the
    # residual gate was relaxed to a warning it must succeed — returning finite
    # matrices — for every bus category, not just well-damped machines.
    @testset "nf_linearization works for every bus category" begin
        nw  = set_IEEE39_PF_init(get_IEEE39_base())
        nws = initialize_from_pf(nw; verbose=false)
        for i in (1, 3, 30, 31, 39)   # junction, load, machine, slack m+l, unctrld m+l
            lti = node_nf_linearization(nw, nws, i)
            for key in ("M", "A", "B", "C", "D")
                @test all(isfinite, reduce(vcat, lti[key]))
            end
            @test length(lti["C"]) == 2          # output is (ln|V|, arg V)
            @test length(lti["B"]) == length(lti["A"])
        end
    end

    # ── Killswitch semantics ────────────────────────────────────────────────
    # `add_killswitch_to_vm` zeros the derivatives of a bus's *differential*
    # states when its `killswitch` parameter is set. This FREEZES the machine's
    # internal dynamics; it does NOT disconnect the bus, which keeps enforcing
    # its current-injection (algebraic) coupling to the network.
    @testset "Killswitch freezes dynamics but keeps the bus connected" begin
        nw  = set_IEEE39_PF_init(get_IEEE39_base(; add_killswitch=true))
        nws = initialize_from_pf(nw; verbose=false)
        nws_kill = deepcopy(nws)
        nws_kill.p.v[30, :killswitch] = 1.0      # freeze the bus-30 machine
        cb   = sc_and_clear(1, 0.1, 0.15, 0.5)   # disturbance on a remote line
        prob = ODEProblem(nw, nws_kill, (0.0, 3.0); add_nw_cb=cb)
        sol  = solve(prob, Rodas4(); abstol=1e-7, reltol=1e-7, maxiters=10^6,
                     dtmin=1e-10, force_dtmin=true)
        @test SciMLBase.successful_retcode(sol.retcode)

        ω30  = sol[VIndex(30, :ctrld_gen₊machine₊ω)]   # frozen machine
        ω32  = sol[VIndex(32, :ctrld_gen₊machine₊ω)]   # free machine
        ur30 = sol[VIndex(30, :busbar₊u_r)]
        ui30 = sol[VIndex(30, :busbar₊u_i)]
        v30  = @. sqrt(ur30^2 + ui30^2)

        @test maximum(abs, ω30 .- ω30[1]) < 1e-9     # speed truly frozen
        @test maximum(abs, ω32 .- ω32[1]) > 1e-5     # neighbouring machine reacts
        @test maximum(abs, v30 .- v30[1]) > 1e-3     # bus voltage still responds
    end

    # ── Distributed-slack initialization & stability ────────────────────────
    # The distributed-slack power-flow/init must produce a genuine equilibrium
    # (not merely "finite"): the network RHS residual should be ~machine-eps.
    # And that equilibrium must be dynamically stable — a small perturbation of
    # the generator speeds (a consistent perturbation of differential states)
    # must decay back to *a* steady state rather than growing without bound.
    @testset "Distributed-slack init is a stable equilibrium" begin
        nw   = set_IEEE39_PF_init(get_IEEE39_base(; distributed_slack=true))
        pfnw = powerflow_model(nw)
        pfs0 = solve_powerflow(nw; pfnw, pfs0=NWState(pfnw), verbose=false)
        nws  = initialize_from_pf(nw; pfs=pfs0,
                   default_overrides=interface_values(pfs0),
                   tol=1e-3, nwtol=1e-3, verbose=false)
        u0 = uflat(nws); p0 = pflat(nws)

        # 1. genuine steady state (much tighter than the 1e-2 base-test bound)
        du = similar(u0); nw(du, u0, p0, NaN)
        @test norm(du, Inf) < 1e-8

        # 2. small perturbation of the generator speeds (consistent IC)
        speed_idx = [VIndex(i, s) for i in 1:39 for s in sym(nw[VIndex(i)])
                     if endswith(string(s), "₊ω")]
        @test !isempty(speed_idx)
        rng = MersenneTwister(TEST_RNG_SEED)
        nws_p = deepcopy(nws)
        for idx in speed_idx
            nws_p[idx] = nws[idx] + 1e-3 * randn(rng)
        end

        prob = ODEProblem(nw, uflat(nws_p), (0.0, 60.0), p0)
        sol  = solve(prob, Rodas4(); abstol=1e-8, reltol=1e-8, maxiters=10^6,
                     dtmin=1e-12, force_dtmin=true)
        @test SciMLBase.successful_retcode(sol.retcode)
        @test all(isfinite, sol.u[end])                 # bounded, no blow-up
        duT = similar(u0); nw(duT, sol.u[end], p0, NaN)
        @test norm(duT, Inf) < 1e-3                      # decayed to a steady state
    end

    # ── Grid-of-grids distributed-slack variation (reconciled semantics) ─────
    # `build_IEEE39_grid_of_grids` now varies the bus-31 lead/follower generation
    # with the same subtract-load / scale-gen / add-back semantics as
    # `generate_gog_powerflow_variation`. Exercise that path (omissions + noise
    # in distributed-slack mode) and confirm it initializes to a steady state.
    @testset "Grid-of-grids distributed slack with omissions and perturbations" begin
        nw, nws = build_IEEE39_grid_of_grids(2;
            pert=0.05, p_omit=0.1, rng=MersenneTwister(TEST_RNG_SEED),
            distributed_slack=true, max_tries=30)
        u0 = uflat(nws); p0 = pflat(nws)
        @test all(isfinite, u0)
        @test all(isfinite, p0)
        du = similar(u0); nw(du, u0, p0, NaN)
        @test norm(du, Inf) < 1e-2
    end

end  # Modelling assumptions
