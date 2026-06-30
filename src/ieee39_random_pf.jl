
using Graphs
using Random

# Internal helper: return dict[key] if present, else perturb base randomly.
_pfv_resolve(d, key, base, pert, rng) =
    isnothing(d) || !haskey(d, key) ? base * (1 + pert * randn(rng)) : d[key]

"""
    generate_powerflow_variation(nw, pfnw, pfs0; pert=0.1,
                                  load_P=nothing, load_Q=nothing,
                                  rng=Random.default_rng())

Generate a power-flow variation of the IEEE39 network `nw` (standard or
distributed-slack) and return `(nws, pfs, loads)`, where `nws` is the
initialised `NWState`, `pfs` is the solved power-flow state, and `loads` is
a named tuple `(; P_31, Q_31, P_39, Q_39)` recording the load set-points used.

## How the variation is applied

For every bus that exposes a `"P"` power-flow parameter the parameter is
multiplied by `1 + pert·ε` (ε ~ N(0,1)).  This covers:
- PQ load buses (standard and DS modes): their `P` set-point.
- PV generator buses (standard mode): their net-injection `P`.
- Distributed-slack lead and follower buses (DS mode): their `Pbase`.
  NetworkDynamics exposes `"P"` as an alias for `"Pbase"`, so the same
  batch operation scales generator set-points in both modes.

Bus 39 (mixed generator + load) is handled specially: its load component is
subtracted from the pf `P` before the batch scale (so only the generator side
is scaled), then a separately-varied load value is added back.  In DS mode the
same code works transparently because `"P"` aliases `"Pbase"` for the follower.

## Keywords
- `pert` – multiplicative noise standard deviation (default 0.1).
- `load_P` – `Dict{Int,Float64}` mapping bus index to an absolute active-power
  set-point.  Overrides the random perturbation for that bus.
- `load_Q` – same for reactive power.
- `rng` – random number generator (default: task-local).
"""
function generate_powerflow_variation(nw, pfnw, pfs0;
                                       pert=0.1,
                                       load_P=nothing,
                                       load_Q=nothing,
                                       rng=Random.default_rng())
    pfs    = deepcopy(pfs0)
    n_vert = nv(nw)

    # Load defaults for the two mixed generator+load buses (31 and 39)
    P31_base = get_default(nw, VIndex(31, :ZIPLoad₊Pset))
    Q31_base = get_default(nw, VIndex(31, :ZIPLoad₊Qset))
    P39_base = get_default(nw, VIndex(39, :ZIPLoad₊Pset))
    Q39_base = get_default(nw, VIndex(39, :ZIPLoad₊Qset))

    # Resolve the final load set-points for gen+load buses
    P_31 = _pfv_resolve(load_P, 31, P31_base, pert, rng)
    Q_31 = _pfv_resolve(load_Q, 31, Q31_base, pert, rng)
    P_39 = _pfv_resolve(load_P, 39, P39_base, pert, rng)
    Q_39 = _pfv_resolve(load_Q, 39, Q39_base, pert, rng)

    # For bus 39: remove the load contribution from the net-injection "P" so
    # that the batch scale below only touches the generator side.
    # In DS mode "P" is an alias for the follower's "Pbase", so this works
    # identically for both standard and distributed-slack networks.
    pfs.v[39].p["P"] -= P39_base

    # Scale every bus's "P" (and "Q") power-flow parameter:
    #   standard mode  – PQ load buses and PV generator buses;
    #                    the slack bus (31) has no "P" and is skipped.
    #   distributed-slack – PQ buses AND all lead/follower "Pbase" values
    #                    (NetworkDynamics exposes them as "P").
    np = length(pfs.v[:].p["P"])
    nq = length(pfs.v[:].p["Q"])
    pfs.v[:].p["P"] .*= (1 .+ pert .* randn(rng, np))
    pfs.v[:].p["Q"] .*= (1 .+ pert .* randn(rng, nq))

    # Add back the independently-varied bus-39 load
    pfs.v[39].p["P"] += P_39

    # Override any directly-specified pure PQ load buses
    pq_load_buses = [b for b in 1:n_vert
                     if :ZIPLoad₊Pset ∈ nw[VIndex(b)].psym && b ∉ (31, 39)]
    if !isnothing(load_P)
        for b in pq_load_buses
            haskey(load_P, b) && (pfs.v[b].p["P"] = load_P[b])
        end
    end
    if !isnothing(load_Q)
        for b in pq_load_buses
            haskey(load_Q, b) && (pfs.v[b].p["Q"] = load_Q[b])
        end
    end

    pfs = solve_powerflow(nw; pfnw, pfs0=pfs, verbose=false)

    internal_variations = Dict(
        VIndex(31, :ZIPLoad₊Pset) => P_31,
        VIndex(31, :ZIPLoad₊Qset) => Q_31,
        VIndex(39, :ZIPLoad₊Pset) => P_39,
        VIndex(39, :ZIPLoad₊Qset) => Q_39,
    )
    default_overrides = merge(interface_values(pfs), internal_variations)
    nws = initialize_from_pf(nw; pfs, tol=1e-3, nwtol=1e-3,
                             default_overrides, verbose=false)
    return nws, pfs, (; P_31, Q_31, P_39, Q_39)
end

"""
    generate_gog_powerflow_variation(nw_gog, pfnw, pfs0; N, n_vert=39,
                                      pert=0.1, load_P=nothing, load_Q=nothing,
                                      distributed_slack=false,
                                      rng=Random.default_rng())

Generate a power-flow variation for a grid-of-grids network built by
`build_IEEE39_grid_of_grids(N; ...)` and return `(nws, pfs)`.

Global bus indices follow the convention `n_vert * sni + b` (copy index `sni`
in `0:N²-1`, local bus `b` in `1:n_vert`).

## Keywords
- `N` – grid side length (total copies `= N²`).
- `n_vert` – local vertices per copy (default 39).
- `pert`, `rng` – as in `generate_powerflow_variation`.
- `load_P` – `Dict{Int,Float64}` mapping **global** bus index to an absolute
  active-power set-point.  Buses absent from the dict are perturbed by `pert`.
- `load_Q` – same for reactive power.
- `distributed_slack` – must match the value used when the network was built.
  Needed to determine whether the centre copy's bus 31 (lead in DS mode,
  slack in standard mode) exposes a `"P"` power-flow parameter.
"""
function generate_gog_powerflow_variation(nw_gog, pfnw, pfs0;
                                           N,
                                           n_vert=39,
                                           pert=0.1,
                                           load_P=nothing,
                                           load_Q=nothing,
                                           distributed_slack=false,
                                           rng=Random.default_rng())
    pfs        = deepcopy(pfs0)
    n_copies   = N * N
    center_sni = (N ÷ 2) * N + (N ÷ 2)
    gv(sni, b) = n_vert * sni + b

    # Step 1: subtract load from gen+load buses before the batch scale so that
    #         only the generator side is scaled.
    #
    #  bus 39 copies – pfPV when intact (has a default for ZIPLoad₊Pset);
    #                  pure-load bus when omitted (no default → skip, let Step 2 scale it).
    #  bus 31 non-centre – pfPV (standard) or DS follower ("P"="Pbase" alias).
    #                       Either way "P" is accessible; subtract load.
    #  bus 31 centre    – pfSlack in standard (no "P" → skip);
    #                       DS lead in DS mode ("P"="Pbase" → subtract).
    for sni in 0:(n_copies - 1)
        gidx39 = gv(sni, 39)
        if :ZIPLoad₊Pset ∈ nw_gog[VIndex(gidx39)].psym &&
                has_default(nw_gog, VIndex(gidx39, :ZIPLoad₊Pset))
            pfs.v[gidx39].p["P"] -= get_default(nw_gog, VIndex(gidx39, :ZIPLoad₊Pset))
        end
        # Centre standard = slack (no P). All other cases have "P".
        if sni != center_sni || distributed_slack
            gidx31 = gv(sni, 31)
            pfs.v[gidx31].p["P"] -= get_default(nw_gog, VIndex(gidx31, :ZIPLoad₊Pset))
        end
    end

    # Step 2: scale all "P" (and "Q") pf parameters in one batch.
    #   - PQ buses: their load P set-point.
    #   - PV buses: their net-injection P (generator side, load already subtracted).
    #   - DS followers and lead: their Pbase (via the "P"→"Pbase" alias).
    #   Standard-mode slack (centre bus 31) has no "P" and is skipped automatically.
    np = length(pfs.v[:].p["P"])
    nq = length(pfs.v[:].p["Q"])
    pfs.v[:].p["P"] .*= (1 .+ pert .* randn(rng, np))
    pfs.v[:].p["Q"] .*= (1 .+ pert .* randn(rng, nq))

    # Step 3: add back the independently-varied loads and collect internal_variations
    internal_variations = Dict()
    for sni in 0:(n_copies - 1)
        # bus 39
        gidx39 = gv(sni, 39)
        if :ZIPLoad₊Pset ∈ nw_gog[VIndex(gidx39)].psym &&
                has_default(nw_gog, VIndex(gidx39, :ZIPLoad₊Pset))
            P39_base = get_default(nw_gog, VIndex(gidx39, :ZIPLoad₊Pset))
            Q39_base = get_default(nw_gog, VIndex(gidx39, :ZIPLoad₊Qset))
            P_39 = _pfv_resolve(load_P, gidx39, P39_base, pert, rng)
            Q_39 = _pfv_resolve(load_Q, gidx39, Q39_base, pert, rng)
            pfs.v[gidx39].p["P"] += P_39
            internal_variations[VIndex(gidx39, :ZIPLoad₊Pset)] = P_39
            internal_variations[VIndex(gidx39, :ZIPLoad₊Qset)] = Q_39
        end
        # bus 31 (all copies)
        gidx31 = gv(sni, 31)
        P31_base = get_default(nw_gog, VIndex(gidx31, :ZIPLoad₊Pset))
        Q31_base = get_default(nw_gog, VIndex(gidx31, :ZIPLoad₊Qset))
        P_31 = _pfv_resolve(load_P, gidx31, P31_base, pert, rng)
        Q_31 = _pfv_resolve(load_Q, gidx31, Q31_base, pert, rng)
        internal_variations[VIndex(gidx31, :ZIPLoad₊Pset)] = P_31
        internal_variations[VIndex(gidx31, :ZIPLoad₊Qset)] = Q_31
        if sni != center_sni || distributed_slack
            pfs.v[gidx31].p["P"] += P_31
        end
    end

    # Step 4: apply any directly-specified load_P / load_Q to pure PQ buses
    if !isnothing(load_P) || !isnothing(load_Q)
        for global_b in 1:nv(nw_gog)
            sni     = (global_b - 1) ÷ n_vert
            local_b = global_b - n_vert * sni
            local_b ∈ (31, 39) && continue
            :ZIPLoad₊Pset ∉ nw_gog[VIndex(global_b)].psym && continue
            !isnothing(load_P) && haskey(load_P, global_b) &&
                (pfs.v[global_b].p["P"] = load_P[global_b])
            !isnothing(load_Q) && haskey(load_Q, global_b) &&
                (pfs.v[global_b].p["Q"] = load_Q[global_b])
        end
    end

    pfs = solve_powerflow(nw_gog; pfnw, pfs0=pfs, verbose=false)
    default_overrides = merge(interface_values(pfs), internal_variations)
    nws = initialize_from_pf(nw_gog; pfs, tol=1e-3, nwtol=1e-3,
                             default_overrides, verbose=false)
    return nws, pfs
end
