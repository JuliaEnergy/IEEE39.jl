using Test
using IEEE39
using PowerDynamics
using NetworkDynamics
using OrdinaryDiffEq
using Graphs
using SciMLBase
using Random
using LinearAlgebra: norm
using JSON

const TEST_RNG_SEED = 42

# ─────────────────────────────────────────────────────────────────────────────
# Helpers shared across testsets
# ─────────────────────────────────────────────────────────────────────────────

function base_pf_setup(; distributed_slack=false)
    nw   = set_IEEE39_PF_init(get_IEEE39_base(; distributed_slack))
    pfnw = powerflow_model(nw)
    pfs0 = NWState(pfnw)
    if distributed_slack
        pfs0 = solve_powerflow(nw; pfnw, pfs0, verbose=false)
    end
    nw, pfnw, pfs0
end

function run_mild_perturbation(nw, nws; T=3.0)
    u0 = uflat(nws); p0 = pflat(nws)
    speed_idx = [VIndex(i, s) for i in 1:nv(nw)
                 for s in sym(nw[VIndex(i)]) if endswith(string(s), "₊ω")]
    rng_local = MersenneTwister(TEST_RNG_SEED)
    nws_p = deepcopy(nws)
    for idx in speed_idx
        nws_p[idx] = nws[idx] + 1e-3 * randn(rng_local)
    end
    prob = ODEProblem(nw, uflat(nws_p), (0.0, T), p0)
    sol  = solve(prob, Rodas4(); abstol=1e-4, reltol=1e-4,
                 maxiters=100_000, dtmin=1e-8, force_dtmin=true)
    @test SciMLBase.successful_retcode(sol.retcode)
    @test all(isfinite, sol.u[end])
end

# ─────────────────────────────────────────────────────────────────────────────
@testset "IEEE39 Base Model" begin

    @testset "Standard initialization" begin
        nw  = set_IEEE39_PF_init(get_IEEE39_base())
        nws = initialize_from_pf(nw; verbose=false)
        u0  = uflat(nws)
        p0  = pflat(nws)
        @test all(isfinite, u0)
        @test all(isfinite, p0)
        du = similar(u0)
        nw(du, u0, p0, NaN)
        @test norm(du, Inf) < 1e-4
    end

    @testset "Distributed slack initialization" begin
        nw, pfnw, pfs0 = base_pf_setup(; distributed_slack=true)
        nws = initialize_from_pf(nw; pfs=pfs0,
                  default_overrides=interface_values(pfs0),
                  tol=1e-3, nwtol=1e-3, verbose=false)
        u0 = uflat(nws); p0 = pflat(nws)
        @test all(isfinite, u0)
        @test all(isfinite, p0)
        du = similar(u0)
        nw(du, u0, p0, NaN)
        @test norm(du, Inf) < 1e-2
    end

    @testset "export_nodes_json viability" begin
        nw  = set_IEEE39_PF_init(get_IEEE39_base())
        nws = initialize_from_pf(nw; verbose=false)
        tmp = tempname() * ".json"
        export_nodes_json(nw, nws, tmp; include_nf=false)
        nodes = JSON.parsefile(tmp)

        @test length(nodes) == nv(nw)
        @test all(n -> haskey(n, "index"),      nodes)
        @test all(n -> haskey(n, "type"),       nodes)
        @test all(n -> haskey(n, "pf_type"),    nodes)
        @test all(n -> haskey(n, "parameters"), nodes)

        types = Dict(n["index"] => n["type"] for n in nodes)
        @test types[1]  == "junction"
        @test types[3]  == "load"
        @test types[30] == "machine"
        @test types[31] == "machine_load"
        @test types[39] == "machine_load"

        vm31 = nw[VIndex(31)]
        n31  = nodes[findfirst(n -> n["index"] == 31, nodes)]
        ep   = n31["parameters"]
        for s in vm31.psym
            @test haskey(ep, string(s))
            @test ep[string(s)] ≈ nws.p.v[31, s]
        end
        rm(tmp)
    end

    @testset "export_lines_json viability" begin
        nw  = set_IEEE39_PF_init(get_IEEE39_base())
        nws = initialize_from_pf(nw; verbose=false)
        tmp = tempname() * ".json"
        export_lines_json(nw, nws, tmp)
        lines = JSON.parsefile(tmp)

        @test length(lines) == ne(nw)
        @test all(l -> haskey(l, "index"),      lines)
        @test all(l -> haskey(l, "src"),        lines)
        @test all(l -> haskey(l, "dst"),        lines)
        @test all(l -> haskey(l, "parameters"), lines)

        l1  = lines[findfirst(l -> l["index"] == 1, lines)]
        @test l1["src"] == 1
        @test l1["dst"] == 2

        em1 = nw[EIndex(1)]
        ep  = l1["parameters"]
        for s in em1.psym
            @test haskey(ep, string(s))
            @test ep[string(s)] ≈ nws.p.e[1, s]
        end
        rm(tmp)
    end

    @testset "Short-circuit-and-clear plausibility" begin
        nw  = set_IEEE39_PF_init(get_IEEE39_base())
        nws = initialize_from_pf(nw; verbose=false)
        u0  = uflat(nws)
        p0  = pflat(nws)

        cb   = sc_and_clear(1, 1.0, 1.05, 0.5)
        prob = ODEProblem(nw, nws, (0.0, 120.0); add_nw_cb=cb)
        sol  = solve(prob, Rodas4(); abstol=1e-6, reltol=1e-6,
                     maxiters=500_000, dtmin=1e-8, force_dtmin=true)

        @test SciMLBase.successful_retcode(sol.retcode)
        @test all(isfinite, sol.u[end])
        du = zeros(length(u0))
        nw(du, sol.u[end], p0, NaN)
        @test norm(du, Inf) < 1e-4
    end

    @testset "Distributed slack – short-circuit-and-clear plausibility" begin
        nw, pfnw, pfs0 = base_pf_setup(; distributed_slack=true)
        nws = initialize_from_pf(nw; pfs=pfs0,
                  default_overrides=interface_values(pfs0),
                  tol=1e-3, nwtol=1e-3, verbose=false)
        u0 = uflat(nws); p0 = pflat(nws)

        cb   = sc_and_clear(1, 1.0, 1.05, 0.5)
        prob = ODEProblem(nw, nws, (0.0, 120.0); add_nw_cb=cb)
        sol  = solve(prob, Rodas4(); abstol=1e-6, reltol=1e-6,
                     maxiters=500_000, dtmin=1e-8, force_dtmin=true)

        @test SciMLBase.successful_retcode(sol.retcode)
        @test all(isfinite, sol.u[end])
        du = zeros(length(u0))
        nw(du, sol.u[end], p0, NaN)
        @test norm(du, Inf) < 1e-4
    end

end  # IEEE39 Base Model

# ─────────────────────────────────────────────────────────────────────────────
@testset "generate_powerflow_variation" begin

    @testset "Standard – random perturbation" begin
        nw, pfnw, pfs0 = base_pf_setup()
        nws, pfs, loads = generate_powerflow_variation(nw, pfnw, pfs0;
                              pert=0.05, rng=MersenneTwister(TEST_RNG_SEED))
        @test all(isfinite, uflat(nws))
        @test all(isfinite, pflat(nws))
        du = similar(uflat(nws))
        nw(du, uflat(nws), pflat(nws), NaN)
        @test norm(du, Inf) < 1e-2
    end

    @testset "Standard – specified load P and Q" begin
        nw, pfnw, pfs0 = base_pf_setup()
        # Set load at bus 3 and bus 20 to fixed values
        P_spec = Dict(3 => -3.0, 20 => -5.5)
        Q_spec = Dict(3 => -0.02, 20 => -0.9)
        nws, pfs, loads = generate_powerflow_variation(nw, pfnw, pfs0;
                              pert=0.05, load_P=P_spec, load_Q=Q_spec,
                              rng=MersenneTwister(TEST_RNG_SEED))
        @test all(isfinite, uflat(nws))
        # Verify the specified values appear in the power-flow solution
        @test pfs.v[3].p["P"]  ≈ P_spec[3]
        @test pfs.v[3].p["Q"]  ≈ Q_spec[3]
        @test pfs.v[20].p["P"] ≈ P_spec[20]
    end

    @testset "Standard – specified load P for gen+load bus 39" begin
        nw, pfnw, pfs0 = base_pf_setup()
        P39_spec = -10.0
        nws, pfs, loads = generate_powerflow_variation(nw, pfnw, pfs0;
                              load_P=Dict(39 => P39_spec),
                              rng=MersenneTwister(TEST_RNG_SEED))
        @test all(isfinite, uflat(nws))
        # The load component of bus 39 must match the specification
        @test loads.P39 ≈ P39_spec
    end

    @testset "Distributed slack – random perturbation" begin
        nw, pfnw, pfs0 = base_pf_setup(; distributed_slack=true)
        pfs0_before = deepcopy(pfs0)
        nws, pfs, loads = generate_powerflow_variation(nw, pfnw, pfs0;
                              pert=0.05, rng=MersenneTwister(TEST_RNG_SEED))
        @test all(isfinite, uflat(nws))
        @test all(isfinite, pflat(nws))
        du = similar(uflat(nws))
        nw(du, uflat(nws), pflat(nws), NaN)
        @test norm(du, Inf) < 1e-2
        # Verify generator Pbase values are varied (not fixed)
        # NetworkDynamics exposes Pbase as "P" for DS followers
        @test pfs.v[30].p["P"] ≠ pfs0_before.v[30].p["P"]   # pure gen bus
        @test pfs.v[32].p["P"] ≠ pfs0_before.v[32].p["P"]   # pure gen bus
    end

    @testset "Distributed slack – specified load P and Q" begin
        nw, pfnw, pfs0 = base_pf_setup(; distributed_slack=true)
        P_spec = Dict(3 => -3.0, 20 => -5.5)
        Q_spec = Dict(3 => -0.02, 20 => -0.9)
        nws, pfs, loads = generate_powerflow_variation(nw, pfnw, pfs0;
                              pert=0.05, load_P=P_spec, load_Q=Q_spec,
                              rng=MersenneTwister(TEST_RNG_SEED))
        @test all(isfinite, uflat(nws))
        @test pfs.v[3].p["P"]  ≈ P_spec[3]
        @test pfs.v[3].p["Q"]  ≈ Q_spec[3]
        @test pfs.v[20].p["P"] ≈ P_spec[20]
    end

end  # generate_powerflow_variation

# ─────────────────────────────────────────────────────────────────────────────
@testset "generate_gog_powerflow_variation" begin

    @testset "Single-slack 2x2 – random perturbation" begin
        nw, nws0 = build_IEEE39_grid_of_grids(2;
            pert=0.0, p_omit=0.0, rng=MersenneTwister(TEST_RNG_SEED))
        pfnw = powerflow_model(nw)
        pfs0 = NWState(pfnw)

        nws, pfs = generate_gog_powerflow_variation(nw, pfnw, pfs0;
                       N=2, pert=0.05, rng=MersenneTwister(TEST_RNG_SEED))
        @test all(isfinite, uflat(nws))
        @test all(isfinite, pflat(nws))
        du = similar(uflat(nws))
        nw(du, uflat(nws), pflat(nws), NaN)
        @test norm(du, Inf) < 1e-2
    end

    @testset "Single-slack 2x2 – specified load P" begin
        nw, nws0 = build_IEEE39_grid_of_grids(2;
            pert=0.0, p_omit=0.0, rng=MersenneTwister(TEST_RNG_SEED))
        pfnw = powerflow_model(nw)
        pfs0 = NWState(pfnw)

        # Fix bus 3 of copy 0 (global index 3) and bus 3 of copy 1 (global 42)
        P_spec = Dict(3 => -3.5, 42 => -3.2)
        nws, pfs = generate_gog_powerflow_variation(nw, pfnw, pfs0;
                       N=2, pert=0.05, load_P=P_spec,
                       rng=MersenneTwister(TEST_RNG_SEED))
        @test all(isfinite, uflat(nws))
        @test pfs.v[3].p["P"]  ≈ P_spec[3]
        @test pfs.v[42].p["P"] ≈ P_spec[42]
    end

    @testset "Distributed-slack 2x2 – random perturbation" begin
        nw, nws0 = build_IEEE39_grid_of_grids(2;
            pert=0.0, p_omit=0.0, rng=MersenneTwister(TEST_RNG_SEED),
            distributed_slack=true)
        pfnw = powerflow_model(nw)
        pfs0 = NWState(pfnw)

        nws, pfs = generate_gog_powerflow_variation(nw, pfnw, pfs0;
                       N=2, pert=0.05, distributed_slack=true,
                       rng=MersenneTwister(TEST_RNG_SEED))
        @test all(isfinite, uflat(nws))
        @test all(isfinite, pflat(nws))
        du = similar(uflat(nws))
        nw(du, uflat(nws), pflat(nws), NaN)
        @test norm(du, Inf) < 1e-2
    end

end  # generate_gog_powerflow_variation

# ─────────────────────────────────────────────────────────────────────────────
@testset "IEEE39 Grid of Grids" begin

    @testset "Single-slack 2x2 grid - initialization" begin
        nw, nws = build_IEEE39_grid_of_grids(2;
            pert=0.0, p_omit=0.0, rng=MersenneTwister(TEST_RNG_SEED))
        u0 = uflat(nws); p0 = pflat(nws)
        @test all(isfinite, u0)
        @test all(isfinite, p0)
        du = similar(u0); nw(du, u0, p0, NaN)
        @test norm(du, Inf) < 1e-2
    end

    @testset "Distributed-slack 2x2 grid - initialization" begin
        nw, nws = build_IEEE39_grid_of_grids(2;
            pert=0.0, p_omit=0.0, rng=MersenneTwister(TEST_RNG_SEED),
            distributed_slack=true)
        u0 = uflat(nws); p0 = pflat(nws)
        @test all(isfinite, u0)
        @test all(isfinite, p0)
        du = similar(u0); nw(du, u0, p0, NaN)
        @test norm(du, Inf) < 1e-2
    end

    @testset "Single-slack 2x2 grid - short-circuit simulation" begin
        nw, nws = build_IEEE39_grid_of_grids(2;
            pert=0.0, p_omit=0.0, rng=MersenneTwister(TEST_RNG_SEED))
        cb   = sc_and_trip(1, 1.0, 1.1, 0.5)
        prob = ODEProblem(nw, nws, (0.0, 2.0); add_nw_cb=cb)
        sol  = solve(prob, Rodas4(); abstol=1e-4, reltol=1e-4, maxiters=10_000,
                     dtmin=1e-6, force_dtmin=true)
        @test SciMLBase.successful_retcode(sol.retcode)
        @test all(isfinite, sol.u[end])
    end

    @testset "Distributed-slack 2x2 grid - short-circuit simulation" begin
        nw, nws = build_IEEE39_grid_of_grids(2;
            pert=0.0, p_omit=0.0, rng=MersenneTwister(TEST_RNG_SEED),
            distributed_slack=true)
        cb   = sc_and_trip(1, 1.0, 1.1, 0.5)
        prob = ODEProblem(nw, nws, (0.0, 2.0); add_nw_cb=cb)
        sol  = solve(prob, Rodas4(); abstol=1e-4, reltol=1e-4, maxiters=10_000,
                     dtmin=1e-6, force_dtmin=true)
        @test SciMLBase.successful_retcode(sol.retcode)
        @test all(isfinite, sol.u[end])
    end

    @testset "Single-slack with omissions and perturbations" begin
        nw, nws = build_IEEE39_grid_of_grids(2;
            pert=0.05, p_omit=0.1, rng=MersenneTwister(TEST_RNG_SEED),
            max_tries=20)
        @test all(isfinite, uflat(nws))
        @test all(isfinite, pflat(nws))
    end

end  # IEEE39 Grid of Grids

# ─────────────────────────────────────────────────────────────────────────────
@testset "Network reconstruction from JSON" begin

    function reconstruction_test(nw, nws; T=3.0)
        nodes_f = tempname() * ".json"
        lines_f = tempname() * ".json"
        export_nodes_json(nw, nws, nodes_f; include_nf=false)
        export_lines_json(nw, nws, lines_f)
        nw_r = import_network_from_json(nodes_f, lines_f)
        rm(nodes_f); rm(lines_f)

        @test nv(nw_r) == nv(nw)
        @test ne(nw_r) == ne(nw)

        u0 = uflat(nws); p0 = pflat(nws)
        sol1 = solve(ODEProblem(nw,   u0, (0.0, T), p0), Rodas4();
                     abstol=1e-8, reltol=1e-8, maxiters=100_000, dtmin=1e-10,
                     force_dtmin=true)
        sol2 = solve(ODEProblem(nw_r, u0, (0.0, T), p0), Rodas4();
                     abstol=1e-8, reltol=1e-8, maxiters=100_000, dtmin=1e-10,
                     force_dtmin=true)
        @test SciMLBase.successful_retcode(sol1.retcode)
        @test SciMLBase.successful_retcode(sol2.retcode)
        @test sol1.u[end] ≈ sol2.u[end]
    end

    @testset "Base – Slack" begin
        nw  = set_IEEE39_PF_init(get_IEEE39_base())
        nws = initialize_from_pf(nw; verbose=false)
        reconstruction_test(nw, nws)
    end

    @testset "Base – DS" begin
        nw, pfnw, pfs0 = base_pf_setup(; distributed_slack=true)
        nws = initialize_from_pf(nw; pfs=pfs0,
                  default_overrides=interface_values(pfs0),
                  tol=1e-3, nwtol=1e-3, verbose=false)
        reconstruction_test(nw, nws)
    end

    @testset "GoG 2x2 – Slack – omit=0.0" begin
        nw, nws = build_IEEE39_grid_of_grids(2;
            pert=0.0, p_omit=0.0, rng=MersenneTwister(TEST_RNG_SEED))
        reconstruction_test(nw, nws)
    end

    @testset "GoG 2x2 – Slack – omit=0.1" begin
        nw, nws = build_IEEE39_grid_of_grids(2;
            pert=0.0, p_omit=0.1, rng=MersenneTwister(TEST_RNG_SEED), max_tries=20)
        reconstruction_test(nw, nws)
    end

end  # Network reconstruction from JSON

# ─────────────────────────────────────────────────────────────────────────────
@testset "Mild perturbation – all valid combinations" begin

    @testset "Base – Slack – omit=0.0" begin
        nw  = set_IEEE39_PF_init(get_IEEE39_base())
        nws = initialize_from_pf(nw; verbose=false)
        run_mild_perturbation(nw, nws)
    end

    @testset "Base – DS – omit=0.0" begin
        nw   = set_IEEE39_PF_init(get_IEEE39_base(; distributed_slack=true))
        pfnw = powerflow_model(nw)
        pfs0 = solve_powerflow(nw; pfnw, pfs0=NWState(pfnw), verbose=false)
        nws  = initialize_from_pf(nw; pfs=pfs0,
                   default_overrides=interface_values(pfs0),
                   tol=1e-3, nwtol=1e-3, verbose=false)
        run_mild_perturbation(nw, nws)
    end

    for ds in (false, true), omit in (0.0, 0.1)
        ds_label   = ds    ? "DS"    : "Slack"
        omit_label = omit == 0.0 ? "omit=0.0" : "omit=0.1"
        @testset "GoG – $ds_label – $omit_label" begin
            nw, nws = build_IEEE39_grid_of_grids(2;
                pert=0.0, p_omit=omit,
                rng=MersenneTwister(TEST_RNG_SEED),
                distributed_slack=ds, max_tries=20)
            run_mild_perturbation(nw, nws)
        end
    end

end  # Mild perturbation

# ─────────────────────────────────────────────────────────────────────────────
# Modelling-assumption / regression tests (documents the assumptions above).
include("modelling_assumptions.jl")
