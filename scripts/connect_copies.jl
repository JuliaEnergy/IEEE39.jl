#=

Build a large network by tiling the IEEE39 grid into an N×N lattice of copies,
wired together with random inter-grid lines, with varied power flow and randomly
omitted generators (see `build_IEEE39_grid_of_grids`), then run a short-circuit
simulation on it.

=#

using Pkg
Pkg.activate(@__DIR__)

using PowerDynamics
using Graphs
using OrdinaryDiffEq
using Plots
using IEEE39
using NetworkDynamics
using SparseConnectivityTracer
using Random
using Logging
# Logging.disable_logging(Logging.Warn)

##

nw_large, nws_init = build_IEEE39_grid_of_grids(2; pert=0.1, p_omit=0.2, rng=MersenneTwister(1))

## Export the node types and parameters of the generated network

export_nodes_json(nw_large, nws_init, joinpath(@__DIR__, "nodes.json"))

##

NetworkDynamics.set_jac_prototype!(nw_large; remove_conditions=true)

##

prob = ODEProblem(nw_large, nws_init, (0.0, 5.0); add_nw_cb=sc_and_trip(1, 0.1, 0.15, 0.1))
tol = 1e-3

t_start = time()
sol = solve(prob, FBDF(), maxiters=10000, abstol=tol, reltol=tol, dtmin=1e-4, force_dtmin=true)
println("solve time: ", time() - t_start)

##

plot(sol, idxs=vidxs(nw_large, 1:39, r"machine₊δ$"))
