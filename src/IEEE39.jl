module IEEE39

using Graphs
using Random
using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using NetworkDynamics
using DataFrames
using CSV
using DiffEqCallbacks
using PrecompileTools
using SciMLBase

include("ieee39_base.jl")
export get_IEEE39_base, set_IEEE39_PF_init

include("ieee39_random_pf.jl")
export generate_powerflow_variation

include("utilities.jl")
export sc_and_trip, sc_and_trip_and_kill

function pre()
    nw_base = get_IEEE39_base()
    nw = set_IEEE39_PF_init(nw_base)
    pfnw = powerflow_model(nw)
    pfs0 = NWState(pfnw)
    pfs_init = initialize_from_pf(nw);
end

precompile(pre, ())

end # module IEEE39
