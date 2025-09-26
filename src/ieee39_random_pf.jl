
using Graphs
using Random

function generate_powerflow_variation(nw_IEEE39_with_PF, pfnw, pfs0; pert=0.1)

    nw = nw_IEEE39_with_PF
    pfs = deepcopy(pfs0)
    # Maybe define copy?

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

    pfs.v[:].p["P"] .*= (1. .+ pert .* randn(length(pfs.v[:].p["P"]))) # Vary the remaining power
    pfs.v[:].p["Q"] .*= (1. .+ pert .* randn(length(pfs.v[:].p["Q"])))

    P_39 = (1. + pert .* randn()) * P_39 # Vary the power of the load at bus 39
    
    pfs.v[39].p["P"] += P_39 # And add the varied load power back to the PV node parameters.

    pfs = solve_powerflow(nw; pfnw, pfs0=pfs, verbose=false)

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
