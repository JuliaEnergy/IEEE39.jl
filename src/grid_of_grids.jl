
using Graphs
using Random
using PowerDynamics
using CSV
using DataFrames

"""
    build_IEEE39_grid_of_grids(N; pert=0.1, p_omit=0.2, rng=Random.default_rng(), max_tries=20)

Build an `N×N` lattice of IEEE39 grid copies wired together into a single large
`Network` and return `(nw_large, nws_init)`. One-off helper for large-network
experiments.

  - Copies are placed on an `N×N` lattice (copy index `r*N + c`, with `r,c` in
    `0:N-1`). Each pair of right/down neighbours is joined by one extra line
    whose parameters are copied from a random existing branch, connecting two
    random buses (one in each neighbouring copy).
  - Only the copy in the centre of the lattice keeps its slack bus (bus 31); in
    every other copy bus 31 is turned into a PV bus, so the whole grid has a
    single slack.
  - Each non-slack generator bus (30, 32–39) in each copy is independently
    dropped with probability `p_omit`: machine-only buses become junction
    (Kirchhoff) buses, the machine+load bus (39) becomes a load bus. To keep each
    copy roughly self-balanced (so the lone central slack and the tie lines are
    not overloaded), the copy's PQ loads are scaled down by the dropped
    generation.
  - The power flow is then varied systematically following
    [`generate_powerflow_variation`](@ref): every generator power and every load
    power is scaled by `(1 + pert·randn())`. At the generator+load buses (39, and
    31 where it is PV) the generator and load powers are varied independently; at
    the central slack only its load is varied.

Because the random draws occasionally produce a power flow that does not
converge, the whole construction is retried with fresh draws (advancing `rng`)
up to `max_tries` times.
"""
function build_IEEE39_grid_of_grids(N; pert=0.1, p_omit=0.2, rng=Random.default_rng(), max_tries=20)
    nw_39 = set_IEEE39_PF_init(get_IEEE39_base())
    n_vert = nv(nw_39)
    n_copies = N * N
    center_sni = (N ÷ 2) * N + (N ÷ 2)   # the single copy that keeps its slack

    vms = nw_39[VIndex(:)]
    ems = nw_39[EIndex(:)]

    # Reuse existing, tested bus models as templates when a generator is dropped.
    junction_template = nw_39[VIndex(1)]   # pure Kirchhoff bus
    load_template = nw_39[VIndex(20)]      # pure load bus
    gen_buses = [30, 32, 33, 34, 35, 36, 37, 38, 39]  # droppable (31 kept as gen)

    # Load setpoints of the two generator+load buses.
    load_P31 = get_default(nw_39[VIndex(31)], :ZIPLoad₊Pset)
    load_Q31 = get_default(nw_39[VIndex(31)], :ZIPLoad₊Qset)
    load_P39 = get_default(nw_39[VIndex(39)], :ZIPLoad₊Pset)
    load_Q39 = get_default(nw_39[VIndex(39)], :ZIPLoad₊Qset)

    # Nominal power-flow setpoints, read straight from the base pf model.
    base_pf = NWState(powerflow_model(nw_39))
    genP = Dict(b => base_pf.v[b].p["P"] for b in gen_buses)        # scheduled gen (bus 39: net)
    pq_load_buses = [b for b in 1:n_vert if (:ZIPLoad₊Pset in nw_39[VIndex(b)].psym) && b ∉ (31, 39)]
    L_pq = sum(abs(base_pf.v[b].p["P"]) for b in pq_load_buses)     # total PQ load magnitude

    # PV setpoint for the bus-31 slacks we convert to PV: the slack's reference
    # net injection and voltage from the bus table.
    bus_df = CSV.read(joinpath(DATA_DIR, "bus.csv"), DataFrame)
    row31 = bus_df[findfirst(==(31), bus_df.bus), :]
    P31_pv, V31 = row31.P, row31.V

    subname(sni, name) = Symbol("net$(sni)_$(name)")
    gv(sni, b) = n_vert * sni + b

    function attempt()
        nw_vertices = []
        nw_edges = []
        load_factor = ones(n_copies)
        omit_per_copy = [Int[] for _ in 1:n_copies]
        for sni in 0:(n_copies - 1)
            offset = n_vert * sni

            # Decide omissions for this copy and how much load to shed to compensate.
            omit = [b for b in gen_buses if rand(rng) < p_omit]
            omit_per_copy[sni + 1] = omit
            # Generation removed: PV buses lose their setpoint; bus 39 turns from a
            # net injector into its full load (so it sheds genP[39] - load_P39).
            P_lost = sum(b == 39 ? genP[39] - load_P39 : genP[b] for b in omit; init=0.0)
            f = clamp((L_pq - P_lost) / L_pq, 0.05, 1.0)
            load_factor[sni + 1] = f

            for vm in vms
                b = get_graphelement(vm)
                gvidx = offset + b
                if b in omit
                    if b == 39  # generator+load bus -> load bus
                        rep = VertexModel(copy(load_template); vidx=gvidx, name=subname(sni, "bus$b"))
                        set_pfmodel!(rep, pfPQ(P=load_P39, Q=load_Q39))
                        set_guess!(rep, :ZIPLoad₊Pset, load_P39)
                        set_guess!(rep, :ZIPLoad₊Qset, load_Q39)
                    else        # machine-only bus -> junction (Kirchhoff) bus
                        rep = VertexModel(copy(junction_template); vidx=gvidx, name=subname(sni, "bus$b"))
                    end
                    push!(nw_vertices, rep)
                else
                    vmc = VertexModel(copy(vm); vidx=gvidx, name=subname(sni, vm.name))
                    # only the centre copy keeps its slack; elsewhere bus 31 -> PV
                    if b == 31 && sni != center_sni
                        set_pfmodel!(vmc, pfPV(P=P31_pv, V=V31))
                    end
                    push!(nw_vertices, vmc)
                end
            end
            for em in ems
                ge = get_graphelement(em)
                push!(nw_edges, EdgeModel(copy(em); src=ge.src + offset, dst=ge.dst + offset))
            end
        end

        # Wire neighbouring copies on the lattice (right + down neighbours).
        for r in 0:(N - 1), c in 0:(N - 1)
            sni = r * N + c
            for nb in (c < N - 1 ? sni + 1 : nothing, r < N - 1 ? sni + N : nothing)
                isnothing(nb) && continue
                em = rand(rng, ems)
                src = rand(rng, 1:n_vert) + n_vert * sni
                dst = rand(rng, 1:n_vert) + n_vert * nb
                push!(nw_edges, EdgeModel(copy(em); src, dst))
            end
        end

        nw_large = Network(nw_vertices, nw_edges)

        # --- systematic power-flow variation, cf. generate_powerflow_variation ---
        pfnw = powerflow_model(nw_large)
        pfs = NWState(pfnw)
        internal_variations = Dict()

        # 1. At the generator+load PV buses, remove the load from the net injection
        #    so that the global scaling below varies the generator power alone.
        for sni in 0:(n_copies - 1)
            39 ∉ omit_per_copy[sni + 1] && (pfs.v[gv(sni, 39)].p["P"] -= load_P39)
            sni != center_sni && (pfs.v[gv(sni, 31)].p["P"] -= load_P31)
        end
        # 2. Shed each copy's PQ loads to compensate the generation it dropped.
        for sni in 0:(n_copies - 1)
            f = load_factor[sni + 1]
            f == 1 && continue
            for b in pq_load_buses
                pfs.v[gv(sni, b)].p["P"] *= f
                pfs.v[gv(sni, b)].p["Q"] *= f
            end
        end
        # 3. Vary every generator and load power.
        pfs.v[:].p["P"] .*= (1 .+ pert .* randn(rng, length(pfs.v[:].p["P"])))
        pfs.v[:].p["Q"] .*= (1 .+ pert .* randn(rng, length(pfs.v[:].p["Q"])))
        # 4. Vary the loads of the generator+load buses independently and add them
        #    back to the net injection; at the central slack the load is varied via
        #    an initialization override only (it has no power-flow setpoint).
        for sni in 0:(n_copies - 1)
            if 39 ∉ omit_per_copy[sni + 1]
                P39 = (1 + pert * randn(rng)) * load_P39
                Q39 = (1 + pert * randn(rng)) * load_Q39
                pfs.v[gv(sni, 39)].p["P"] += P39
                internal_variations[VIndex(gv(sni, 39), :ZIPLoad₊Pset)] = P39
                internal_variations[VIndex(gv(sni, 39), :ZIPLoad₊Qset)] = Q39
            end
            P31 = (1 + pert * randn(rng)) * load_P31
            Q31 = (1 + pert * randn(rng)) * load_Q31
            sni != center_sni && (pfs.v[gv(sni, 31)].p["P"] += P31)
            internal_variations[VIndex(gv(sni, 31), :ZIPLoad₊Pset)] = P31
            internal_variations[VIndex(gv(sni, 31), :ZIPLoad₊Qset)] = Q31
        end

        pfs = solve_powerflow(nw_large; pfnw, pfs0=pfs, verbose=false)
        nws = initialize_from_pf(nw_large; pfs, tol=1e-3, nwtol=1e-3,
            default_overrides=merge(interface_values(pfs), internal_variations), verbose=false)
        return nw_large, nws
    end

    for attempt_i in 1:max_tries
        try
            return attempt()
        catch
            attempt_i == max_tries && rethrow()
        end
    end
end
