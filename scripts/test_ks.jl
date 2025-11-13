
using Pkg
Pkg.activate(@__DIR__)

using PowerDynamics
using Graphs
using Plots
using Random
using OrdinaryDiffEq
using IEEE39
using StatsBase: sample
using SciMLBase
using DSP
using Tables
using Arrow
using SparseConnectivityTracer
using Logging
Logging.disable_logging(Logging.Warn)

##

function simulate_sc_and_kill(nw, nws; t_sc_to_trip, sc_bus, pos, kill_list, solver=Rodas4(), tol=1e-3, maxiters=10000, abstol=tol, reltol=tol, dtmin = 1e-4, force_dtmin = true)
    prob = ODEProblem(nw, copy(uflat(nws)), (0.0, 5.0), copy(pflat(nws)); callback=sc_and_trip_and_kill(sc_bus, 0.1, 0.1 + t_sc_to_trip, pos, kill_list))
    solve(prob, solver; saveat=0.:0.01:5.0, maxiters, abstol, reltol, dtmin, force_dtmin)
end

##
nw_base = get_IEEE39_base(;add_killswitch=true)
nw_39 = set_IEEE39_PF_init(nw_base)
pfnw_39 = powerflow_model(nw_39)
pfs0_39 = NWState(pfnw_39)

# NetworkDynamics.set_jac_prototype!(nw_39; dense=true)

##

nws = initialize_from_pf(nw_39)

##

sol = simulate_sc_and_kill(nw_39, nws; t_sc_to_trip = 0.1, pos=0.1, sc_bus = 1, kill_list=[39])

plot(sol, idxs=vidxs(nw_39, 39, "u_mag"))

##