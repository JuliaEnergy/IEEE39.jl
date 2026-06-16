
using Graphs
using Random
using CSV
using DataFrames
using Dates



DATA_DIR = joinpath(pkgdir(@__MODULE__), "ieee39data")

function generate_powerflow_variation_RES(nw_IEEE39_with_PF, pfnw, pfs0; pert=0.1)
    datapath = joinpath("C:\\", "Users", "lbitcher", "Downloads", "input-files", "Input files", "RT", "Solar")
    load_pos = [3, 4, 7, 8, 12, 15, 16, 18, 20, 21, 23, 24, 25, 26, 27, 28, 29, 31, 39]

    nw = nw_IEEE39_with_PF
    pfs = deepcopy(pfs0)
    # Maybe define copy?
    new_gen = generate_gen_data(datapath)
    random_state = trunc(Int,rand()*366)

    # 31 is the slack bus in the power flow, but we still can vary the load located at the bus
    # initialize_from_pf takes parameters from the defaults specified when the component is constructed
    # We can obtain these defaults, and override them when calling initialize_from_pf 

    P_31 = get_default(nw, VIndex(31, :ZIPLoad₊Pset)) * (1. + pert .* randn())
    Q_31 = get_default(nw, VIndex(31, :ZIPLoad₊Qset)) * (1. + pert .* randn())

    # 39 is another large generator and load, and modeled as PV on the powerflow side
    Q_39 = get_default(nw, VIndex(39, :ZIPLoad₊Qset)) * (1. + pert .* randn())

    # For the power we want to independently vary the power supplied by the generator and the load at bus 39.

    P_39 = get_default(nw, VIndex(39, :ZIPLoad₊Pset))

    pfs.v[39].p["P"] -= P_39 # Thus we substract the power supplied by the load from the PV node parameters.
    P_39 = (1. + pert .* randn()) * P_39 # Vary the power of the load at bus 39

    pfs.v[:].p["P"] .*= (1. .+ pert .* randn(length(pfs.v[:].p["P"])))
    init_sum = sum(pfs0.v[load_pos[1:end-1]].p["P"]) + P_39
    gen_sum = sum(new_gen[1:9,random_state])
    scaling_factor = -init_sum / gen_sum

    pfs.v[[30,32,33,34,35,36,37,38,39]].p["P"] .= new_gen[1:9,random_state] .* scaling_factor
    pfs.v[:].p["Q"] .*= (1. .+ pert .* randn(length(pfs.v[:].p["Q"])))
    println(pfs.v[:].p["Q"])

    
    pfs.v[39].p["P"] += P_39 # And add the varied load power back to the PV node parameters.
    println(pfs.v[:].p["P"])

    pfs = solve_powerflow(nw; pfnw, pfs0=pfs, verbose=false)
    println("after Powerflow")
    # Provide the overrides for the varied power and reactive power of the loads at bus 31 and 39 
    internal_variations = Dict(
        VIndex(31, :ZIPLoad₊Pset) => P_31,
        VIndex(31, :ZIPLoad₊Qset) => Q_31,
        VIndex(39, :ZIPLoad₊Pset) => P_39,
        VIndex(39, :ZIPLoad₊Qset) => Q_39
    )

    default_overrides = merge(interface_values(pfs), internal_variations)

    nw_state = initialize_from_pf(nw; pfs = pfs, tol=1e-3, nwtol=1e-3, default_overrides, verbose=false)
    return nw_state, pfs, (;P_31, Q_31, P_39, Q_39)
end

function generate_powerflow_variation_loads(nw_IEEE39_with_PF, pfnw, pfs0; pert=0.1)
    datapath = joinpath("C:\\", "Users", "lbitcher", "Coding", "PDmodels", "IEEE39.jl", "ieee39data", "Loads.csv")
    load_pos = [3, 4, 7, 8, 12, 15, 16, 18, 20, 21, 23, 24, 25, 26, 27, 28, 29]#, 31, 39]
    other_pos = setdiff(1:39, load_pos)

    nw = nw_IEEE39_with_PF
    pfs = deepcopy(pfs0)
    # Maybe define copy?
    ts, new_loads = generate_load_data(datapath)
    random_state = trunc(Int,rand()*length(ts))
    # println(random_state)

    # 31 is the slack bus in the power flow, but we still can vary the load located at the bus
    # initialize_from_pf takes parameters from the defaults specified when the component is constructed
    # We can obtain these defaults, and override them when calling initialize_from_pf 

    P_31 = get_default(nw, VIndex(31, :ZIPLoad₊Pset)) * (1. + pert .* randn())
    Q_31 = get_default(nw, VIndex(31, :ZIPLoad₊Qset)) * (1. + pert .* randn())

    # 39 is another large generator and load, and modeled as PV on the powerflow side
    Q_39 = get_default(nw, VIndex(39, :ZIPLoad₊Qset)) * (1. + pert .* randn())
    # For the power we want to independently vary the power supplied by the generator and the load at bus 39.
    P_39 = get_default(nw, VIndex(39, :ZIPLoad₊Pset))

    pfs.v[39].p["P"] -= P_39 # Thus we substract the power supplied by the load from the PV node parameters.
    # println(length(collect(new_loads[random_state, 1:end-2])))
    # println(pfs.v[load_pos].p["P"])
    pfs.v[load_pos].p["P"] .= collect(new_loads[random_state, 1:end-2])
    # init_sum = sum(pfs0.v[1:29].p["P"])
    # gen_sum = sum(new_gen[1:8,random_state])
    # scaling_factor = -init_sum / gen_sum

    pfs.v[other_pos].p["P"] .*= (1. .+ pert .* randn(length(pfs.v[other_pos].p["P"])))
    pfs.v[:].p["Q"] .*= (1. .+ pert .* randn(length(pfs.v[:].p["Q"])))
    # println(pfs.v[:].p["Q"])

    P_39 = (1. + pert .* randn()) * P_39 # Vary the power of the load at bus 39
    
    pfs.v[39].p["P"] += P_39 # And add the varied load power back to the PV node parameters.
    # println(pfs.v[:].p["P"])

    pfs = solve_powerflow(nw; pfnw, pfs0=pfs, verbose=false)
    println("after Powerflow")
    # Provide the overrides for the varied power and reactive power of the loads at bus 31 and 39 
    internal_variations = Dict(
        VIndex(31, :ZIPLoad₊Pset) => P_31,
        VIndex(31, :ZIPLoad₊Qset) => Q_31,
        VIndex(39, :ZIPLoad₊Pset) => P_39,
        VIndex(39, :ZIPLoad₊Qset) => Q_39
    )

    default_overrides = merge(interface_values(pfs), internal_variations)

    nw_state = initialize_from_pf(nw; pfs = pfs, tol=1e-3, nwtol=1e-3, default_overrides, verbose=false)
    return nw_state, pfs, (;P_31, Q_31, P_39, Q_39)
end

function generate_load_data(input::String)
    data = CSV.read(input, DataFrame, header=1)
    # println(data[2:5,1])
    res = map(x -> DateTime(x,"yyyy.mm.dd HH.MMM.SS"), data[2:end,1])
    # idxs = findall(x -> Dates.hour(x) .== 13, res)
    output = parse.(Float64, data[2:end,2:end])
    return res,output
end

generate_load_data(datapath)


function generate_gen_data(input::String)
    directory = readdir(input)
    # data = CSV.read(input, DataFrame)
    output = zeros(9,366)
    i=1
    for file in directory
        !endswith(file,".csv") && continue
        data = CSV.read(joinpath(input, file), DataFrame)
        # println(data[9:17,1])
        res = map(x -> DateTime(x,"mm/dd/yy HH:MMM"), data[:,1])
        idxs = findall(x -> Dates.hour(x) .== 13, res)
        output[i,:] = data[idxs,2]
        # println(res)
        i+=1
        i>9 && break
    end
    # if endswith(input, ".csv")
    #     data = CSV.read(input, DataFrame)
    #     return data
    # end
    return output
    
end

# generate_gen_data(joinpath("C:\\", "Users", "lbitcher", "Downloads", "input-files", "Input files", "RT", "Solar"));
