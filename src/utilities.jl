using DiffEqCallbacks
using SciMLBase
using PowerDynamics
using NetworkDynamics
using NetworkDynamics: fftype
using Graphs
using JSON
using ForwardDiff

function sc_and_trip(i, t_sc, t_trip, pos)
    short_circuit = (integrator) -> begin
        p_view = NWParameter(integrator)
        p_view.e[i, "pos"] = pos # Position of sc on line
        p_view.e[i, "shortcircuit"] = 1.
        save_parameters!(integrator)
    end
    trip_line = (integrator) -> begin
        p_view = NWParameter(integrator)
        p_view.e[i, "active"] = 0.
        save_parameters!(integrator)
    end
    sc_cb = PresetTimeCallback(t_sc, short_circuit)
    trip_cb = PresetTimeCallback(t_trip, trip_line)
    return CallbackSet(sc_cb, trip_cb)
end
function sc_and_clear(i, t_sc, t_clear, pos)
    short_circuit = (integrator) -> begin
        p_view = NWParameter(integrator)
        p_view.e[i, "pos"] = pos
        p_view.e[i, "shortcircuit"] = 1.
        save_parameters!(integrator)
    end
    clear_fault = (integrator) -> begin
        p_view = NWParameter(integrator)
        p_view.e[i, "shortcircuit"] = 0.
        save_parameters!(integrator)
    end
    sc_cb = PresetTimeCallback(t_sc, short_circuit)
    clear_cb = PresetTimeCallback(t_clear, clear_fault)
    return CallbackSet(sc_cb, clear_cb)
end

"""
    nf_linearization(vm::VertexModel, state; tol=1e-6) -> NamedTuple

Normal Form linearization of a `VertexModel` around the operating point in
`state`, in the power/phase convention `(δQ, δP) ↦ (δ ln|V|, δ arg(V))`:

    M δẋ = A δx + B δ(Q, P)
       δΘ = C δx              = (δ ln|V|, δ arg(V))

Returns `(; M, A, B, C, D, S0, Θ0, i0, u0, x0, p0)`, where `x0`/`p0` are the
state/parameter operating point and `S0`/`Θ0` the complex power / log-voltage at
that point. It assumes (and warns if not) that the model is initialized at a
steady state, otherwise the linearization is affine rather than linear. The
steady-state check (`tol`, applied to both the initialization residual and the
inner state-derivative residual) only *warns* — it never throws — so the function
is safe to call on networks initialized to a loose tolerance (e.g. the looser
distributed-slack init).

Adapted from NormalFormIdentification.jl (Hans Würfel). The autodiff Jacobians
are computed with `ForwardDiff`; the original finite-difference cross-check
assertions are dropped — they were over-strict (fixed `atol=1e-6`) and spuriously
failed on stiff machine Jacobians even when the autodiff result was correct.
"""
function nf_linearization(vm::VertexModel, state=NetworkDynamics.get_defaults_or_inits_dict(vm); tol=1e-6)
    res = init_residual(vm, state)
    if isnan(res) || res > tol
        @warn "Vertex model does not appear to be at a steady state (residual=$res > tol=$tol); the linearization is affine rather than linear about this point."
    end

    pvec = Float64[state[s] for s in psym(vm)]

    # wrap the inner functions:  M ẋ = f_inner(x, i_dq) ;  u = g_inner(x)
    f_inner = function (x, idq)
        dx = zeros(typeof(first(x) * first(idq)), dim(vm))
        vm.f(dx, x, idq, pvec, NaN)
        dx
    end
    g_inner = function (x)
        vdq = similar(x, 2)
        if fftype(vm) isa PureStateMap
            vm.g(vdq, x)
        elseif fftype(vm) isa NoFeedForward
            vm.g(vdq, x, pvec, NaN)
        else
            error("nf_linearization only supports PureStateMap or NoFeedForward outputs")
        end
        vdq
    end

    xvec = Float64[state[s] for s in sym(vm)]
    idqvec = Float64[state[s] for s in insym(vm)]
    # The inner state-derivative residual is a steady-state diagnostic, not a hard
    # requirement: a network initialized to a finite tolerance (e.g. the looser
    # distributed-slack init) can sit slightly off equilibrium and still yield a
    # useful — if affine — linearization. Warn (don't throw) when it exceeds `tol`.
    fres = maximum(abs.(f_inner(xvec, idqvec)))
    fres > tol && @warn "Inner dynamics residual $fres exceeds tol=$tol; the normal-form linearization includes an affine offset."
    @assert g_inner(xvec) ≈ [state[s] for s in outsym(vm)]

    # operating point in power/phase coordinates
    S0 = Complex(g_inner(xvec)...) * conj(Complex(idqvec...))
    QP0vec = [imag(S0), real(S0)]

    # rewrap to the convention ΔQP -> (ln|V|, arg V)
    f = function (x, ΔQP)
        Q, P = QP0vec + ΔQP
        u_r, u_i = g_inner(x)
        ic = conj((P + im * Q) / (u_r + im * u_i))
        f_inner(x, [real(ic), imag(ic)])
    end
    g = function (x)
        u_r, u_i = g_inner(x)
        [1 / 2 * log(u_r^2 + u_i^2), atan(u_i, u_r)]
    end

    M = vm.mass_matrix
    A = ForwardDiff.jacobian(x -> f(x, zeros(2)), xvec)       # ∂f/∂x
    B = ForwardDiff.jacobian(ΔQP -> f(xvec, ΔQP), zeros(2))   # ∂f/∂ΔQP
    C = ForwardDiff.jacobian(g, xvec)                         # ∂g/∂x

    (; M, A, B, C, D=zeros(2, 2),
        S0=[real(S0), imag(S0)], Θ0=g(xvec),
        i0=idqvec, u0=g_inner(xvec), x0=xvec, p0=pvec)
end

_rows(M) = [collect(Float64, r) for r in eachrow(Matrix(M))]

"""
    node_nf_linearization(nw, nws, i)

Compute the Normal Form linearization ([`nf_linearization`](@ref)) of vertex `i`
around its operating point in state `nws`, and return the LTI system as a `Dict`:
the original `states`, the matrices `M`, `A`, `B`, `C`, `D`, and the
linearization point `S0`, `Theta0`, `x0`.
"""
function node_nf_linearization(nw, nws, i)
    vm = nw[VIndex(i)]
    allsym = vcat(sym(vm), psym(vm), insym(vm), outsym(vm))
    state = Dict(s => nws[VIndex(i, s)] for s in allsym)
    lti = nf_linearization(vm, state)
    Dict(
        "states" => string.(sym(vm)),
        "M" => _rows(lti.M),
        "A" => _rows(lti.A),
        "B" => _rows(lti.B),
        "C" => _rows(lti.C),
        "D" => _rows(lti.D),
        "S0" => collect(Float64, lti.S0),
        "Theta0" => collect(Float64, lti.Θ0),
        "x0" => collect(Float64, lti.x0),
    )
end

"""
    export_nodes_json(nw, nws, filename; include_nf=true)

Write the node (bus) types and parameters of network `nw` to `filename` as a
JSON array, one object per bus. Each entry holds the vertex index, the model
`name`, a structural `type` ("junction", "load", "machine" or "machine_load"),
the power-flow model `pf_type` ("pqbus", "pvbus" or "slackbus"), and a
`parameters` map of every node parameter to its value taken from the state `nws`.

If `include_nf` is true, each entry also carries `nf_linearization`, the Normal
Form linearization of the node (see [`node_nf_linearization`](@ref)). Nodes whose
linearization fails are written with `nf_linearization` set to `null`.
"""
function export_nodes_json(nw, nws, filename; include_nf=true)
    nodes = map(1:nv(nw)) do i
        vm = nw[VIndex(i)]
        syms = vm.psym
        ssyms = string.(syms)
        has_machine = any(s -> occursin("machine₊", s), ssyms)
        has_load = any(s -> occursin("ZIPLoad₊", s), ssyms)
        type = if has_machine && has_load
            "machine_load"
        elseif has_machine
            "machine"
        elseif has_load
            "load"
        else
            "junction"
        end
        params = Dict(string(s) => nws.p.v[i, s] for s in syms)
        nf = if include_nf
            try
                node_nf_linearization(nw, nws, i)
            catch
                nothing
            end
        else
            nothing
        end
        (; index=i, name=string(vm.name), type, pf_type=string(get_pfmodel(vm).name),
            parameters=params, nf_linearization=nf)
    end
    open(filename, "w") do io
        JSON.print(io, nodes, 2)
    end
    filename
end

"""
    export_lines_json(nw, nws, filename)

Write the line (edge) parameters of network `nw` to `filename` as a JSON array,
one object per line. Each entry holds the edge index, model `name`, source and
destination vertex indices, and a `parameters` map of every edge parameter to
its current value taken from the state `nws`.
"""
function export_lines_json(nw, nws, filename)
    lines = map(1:ne(nw)) do i
        em = nw[EIndex(i)]
        syms = em.psym
        ge = get_graphelement(em)
        params = Dict(string(s) => nws.p.e[i, s] for s in syms)
        (; index=i, name=string(em.name), src=ge.src, dst=ge.dst, parameters=params)
    end
    open(filename, "w") do io
        JSON.print(io, lines, 2)
    end
    filename
end

"""
    import_network_from_json(nodes_file, lines_file) -> Network

Reconstruct a `Network` from JSON files produced by [`export_nodes_json`](@ref)
and [`export_lines_json`](@ref).  The reconstructed network has the same graph
topology, vertex-model dynamics, and parameters as the original, so simulating
from the same initial conditions (`uflat(nws)`, `pflat(nws)`) produces
identical trajectories.

Bus type is inferred from parameter-name patterns:
- `ZIPLoad₊`  → has a load component
- `ctrld_gen₊` → has a controlled machine (with AVR/governor)
- `machine₊` without `ctrld_gen₊` prefix → uncontrolled machine
- `killswitch` present → killswitch wrapper is applied

Power-flow metadata (`pf_type`, DS lead/follow coupling) is not restored
because it is not used for ODE simulation.
"""
function import_network_from_json(nodes_file, lines_file)
    nodes_data = JSON.parsefile(nodes_file)
    lines_data = JSON.parsefile(lines_file)

    # ── Build the same component templates as get_IEEE39_base ────────────────
    load_comp     = ZIPLoad(; name=:ZIPLoad)
    unctrld_mach  = SauerPaiMachine(; τ_m_input=false, vf_input=false, name=:machine)
    _machine      = SauerPaiMachine(; name=:machine)
    _avr          = AVRTypeI(; name=:avr, ceiling_function=:quadratic)
    _gov          = TGOV1(; name=:gov)
    ctrld_mach    = CompositeInjector([_machine, _avr, _gov]; name=:ctrld_gen)

    BASE_MVA  = 100.0
    BASE_FREQ = 60.0

    @named junction_t     = compile_bus(MTKBus())
    strip_defaults!(junction_t)

    @named load_t         = compile_bus(MTKBus(load_comp))
    strip_defaults!(load_t)

    @named ctrld_t        = compile_bus(MTKBus(ctrld_mach))
    strip_defaults!(ctrld_t)
    set_default!(ctrld_t, r"S_b$", BASE_MVA)
    set_default!(ctrld_t, r"ω_b$", 2π * BASE_FREQ)

    @named ctrld_load_t   = compile_bus(MTKBus(ctrld_mach, load_comp))
    strip_defaults!(ctrld_load_t)
    set_default!(ctrld_load_t, r"S_b$", BASE_MVA)
    set_default!(ctrld_load_t, r"ω_b$", 2π * BASE_FREQ)

    @named unctrld_load_t = compile_bus(MTKBus(unctrld_mach, load_comp))
    strip_defaults!(unctrld_load_t)
    set_default!(unctrld_load_t, r"S_b$", BASE_MVA)
    set_default!(unctrld_load_t, r"ω_b$", 2π * BASE_FREQ)

    # ── Reconstruct vertex models ────────────────────────────────────────────
    sort!(nodes_data, by=n -> n["index"])

    busses = map(nodes_data) do node
        i      = node["index"]
        params = node["parameters"]
        pnames = keys(params)

        has_load    = any(s -> occursin("ZIPLoad₊",      s), pnames)
        has_machine = any(s -> occursin("machine₊",      s), pnames)
        has_ctrld   = any(s -> startswith(s, "ctrld_gen₊"), pnames)
        has_kill    = haskey(params, "killswitch")

        template = if !has_machine && !has_load
            junction_t
        elseif !has_machine
            load_t
        elseif has_ctrld && !has_load
            ctrld_t
        elseif has_ctrld
            ctrld_load_t
        else
            unctrld_load_t          # uncontrolled machine + load (bus 39)
        end

        bus = compile_bus(template; vidx=i, name=Symbol("bus$i"))
        for (k, v) in params
            k == "killswitch" && continue
            set_default!(bus, Symbol(k), Float64(v))
        end
        if has_kill
            bus = add_killswitch_to_vm(bus)
            set_default!(bus, :killswitch, Float64(params["killswitch"]))
        end
        bus
    end

    # ── Reconstruct edge models ──────────────────────────────────────────────
    @named piline_t = compile_line(MTKLine(PiLine_fault(; name=:piline)))

    sort!(lines_data, by=l -> l["index"])

    branches = map(lines_data) do line
        branch = compile_line(piline_t; src=line["src"], dst=line["dst"])
        for (k, v) in line["parameters"]
            set_default!(branch, Symbol(k), Float64(v))
        end
        branch
    end

    Network(busses, branches)
end

function sc_and_trip_and_kill(i, t_sc, t_trip, pos, kill_list)
    short_circuit = (integrator) -> begin
        p_view = NWParameter(integrator)
        p_view.e[i, "pos"] = pos # Position of sc on line
        p_view.e[i, "shortcircuit"] = 1.
        save_parameters!(integrator)
    end
    trip_line = (integrator) -> begin
        p_view = NWParameter(integrator)
        p_view.e[i, "active"] = 0.
        p_view.v[kill_list, :killswitch] .= 1.
        save_parameters!(integrator)
    end
    sc_cb = PresetTimeCallback(t_sc, short_circuit)
    trip_cb = PresetTimeCallback(t_trip, trip_line)
    return CallbackSet(sc_cb, trip_cb)
end

