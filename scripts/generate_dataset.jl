#=

This script uses the IEEE39 Bus system provided by this package
to generate a dataset of short circuit simulations.

The dataset consist of random variations of the power flow.
In 40% of cases we also deactivte one line, in a further 20% we deactivate two lines.

For every simulation, two variations of the transient stability index are calcualted and saved.

The script also shows how to rerun individual simulations from the dataset, and plot their behavior.

=#

##

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

function simulate_sc(nw, nws; t_sc_to_trip, sc_bus, pos, solver=Rodas4(), tol=1e-3, maxiters=10000, abstol=tol, reltol=tol, dtmin = 1e-4, force_dtmin = true)
    prob = ODEProblem(nw, nws, (0.0, 5.0); add_nw_cb=sc_and_trip(sc_bus, 0.1, 0.1 + t_sc_to_trip, pos))
    solve(prob, solver; saveat=0.:0.01:5.0, maxiters, abstol, reltol, dtmin, force_dtmin)
end


##

_tsi(delta_a) = (pi - delta_a) / (pi + delta_a)

function tsi(angles::Matrix)
    mean_angle = angle(sum(exp.(1.0im .* angles[:, 1])))
    angles .-= mean_angle .- π
    angles = mod2pi.(angles)
    angles = unwrap(angles; dims=2)
    max_a = maximum(angles; dims=1) |> vec
    min_a = minimum(angles; dims=1) |> vec
    delta_a = max_a .- min_a
    minimum(_tsi.(delta_a))
end

function analyze_sol(nw, nw_state, sol)
    if SciMLBase.successful_retcode(sol)
        machine_angles = hcat(sol(sol.t, idxs=vidxs(nw, :, r"machine₊δ$")).u...)
        tsi_machine_delta = tsi(machine_angles)
        bus_angles = hcat(sol(sol.t, idxs=vidxs(nw, :, :busbar₊u_arg)).u...)
        tsi_bus_delta = tsi(bus_angles)
        return (; successful_retcode = true, tsi_machine_delta, tsi_bus_delta)
    else
        return (; successful_retcode = true, tsi_machine_delta = NaN, tsi_bus_delta = NaN)
    end
end

##

function gen_NNF_and_run_short_circuits(nw, pfnw, pfs0; debug=false, verbose=false)

    pfs = deepcopy(pfs0)
    
    # Initial distribution:
    # 40% N0
    # 40% N1
    # 20% N2

    edges = collect(1:46)
    e_fail = Int[]
    thresh = rand()
    if 0.4 < thresh < 0.8
        global e_fail = sample(edges, 1, replace=false)
        pfs.e[e_fail[1]].p["active"] = 0.
    elseif thresh >= 0.8
        global e_fail = sample(edges, 2, replace=false)
        pfs.e[e_fail[1]].p["active"] = 0.
        pfs.e[e_fail[2]].p["active"] = 0.
    end

    sc_buses = setdiff(edges, e_fail)
    
    nws_varied, pfs_varied, (;P_31, Q_31, P_39, Q_39) = try
        generate_powerflow_variation(nw, pfnw, pfs)
    catch e
        verbose && println(e)
        return nothing
    end

    trip_times = 0.05:0.05:0.2
    sc_positions = [0.1, 0.9]

    # For convenience we package all possible combinations of trip_times, positions and sc_bus into an array of named tuples.
    experiments = [(; t_sc_to_trip, sc_bus, pos) for (t_sc_to_trip, sc_bus, pos) in Base.Iterators.product(trip_times, sc_buses, sc_positions)]

    res = []

    if debug
        experiments = sample(experiments, 10, replace=false)
    end 

    for exp_paras in experiments
        verbose && println(exp_paras)
        sol = simulate_sc(nw, nws_varied; exp_paras...)
        push!(res, (; exp_paras... , analyze_sol(nw, nws_varied, sol)...))
    end

    sim_results = columntable(res)
    pg_state = (; e_fail, P_31, Q_31, P_39, Q_39, u_pf = uflat(pfs_varied), p_pf = pflat(pfs_varied), u_nw = uflat(nws_varied), p_nw = pflat(nws_varied))    
    
    return (; sim_results..., pg_state...)
end

##

function generate_dataset(N, filename, nw, pfnw, pfs0; debug=false, verbose=false)

    verbose && println("Running simulation 1/$N")
    data = nothing
    while isnothing(data)
        data = gen_NNF_and_run_short_circuits(nw, pfnw, pfs0; verbose, debug)
        verbose && println("No data generated")
    end

    results = (; [k => [data[k]] for k in keys(data)]...)
    for i in 2:N

        verbose && println("Running simulation $i/$N")
        data = nothing
        while isnothing(data)
            data = gen_NNF_and_run_short_circuits(nw, pfnw, pfs0; verbose, debug)
            verbose && println("No data generated")
        end

        for k in keys(data)
            push!(results[k], data[k])
        end
    end
    Arrow.write(filename, results)
end

##
nw_base = get_IEEE39_base()
nw_39 = set_IEEE39_PF_init(nw_base)
pfnw_39 = powerflow_model(nw_39)
pfs0_39 = NWState(pfnw_39)

NetworkDynamics.set_jac_prototype!(nw_39; remove_conditions=true)


##
data_file = joinpath(@__DIR__, "debug.arrow")
Random.seed!(42)

t_start = time()
wd = generate_dataset(10, data_file, nw_39, pfnw_39, pfs0_39; verbose=true, debug=true)
t_end = time()
println(t_end - t_start)

##

data = Arrow.Table(data_file)

##

NNF = 3
exp_idx = 5

t_sc_to_trip = data[:t_sc_to_trip][NNF][exp_idx]
sc_bus = data[:sc_bus][NNF][exp_idx]
pos = data[:pos][NNF][exp_idx]
tsi_machine_delta = data[:tsi_machine_delta][NNF][exp_idx]
tsi_bus_delta = data[:tsi_bus_delta][NNF][exp_idx]


nw_initial_state = NWState(nw_39, copy(data[:u_nw][NNF]), copy(data[:p_nw][NNF]))
sol_test = simulate_sc(nw_39, nw_initial_state; t_sc_to_trip, sc_bus, pos)
results = analyze_sol(nw_39, nw_initial_state, sol_test)

println("Results match: $(isapprox(tsi_machine_delta, results.tsi_machine_delta)), $(isapprox(tsi_bus_delta, results.tsi_bus_delta))")

plot(sol_test, idxs=vidxs(nw_39, :, r"machine₊δ$"))

##

plot(sol_test, idxs=vidxs(nw_39, :, :busbar₊u_arg), legend=false)

