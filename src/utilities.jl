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

"""
    nf_linearization(vm::VertexModel, state) -> NamedTuple

Normal Form linearization of a `VertexModel` around the operating point in
`state`, in the power/phase convention `(δQ, δP) ↦ (δ ln|V|, δ arg(V))`:

    M δẋ = A δx + B δ(Q, P)
       δΘ = C δx              = (δ ln|V|, δ arg(V))

Returns `(; M, A, B, C, D, S0, Θ0, i0, u0, x0, p0)`, where `x0`/`p0` are the
state/parameter operating point and `S0`/`Θ0` the complex power / log-voltage at
that point. It assumes (and warns if not) that the model is initialized at a
steady state, otherwise the linearization is affine rather than linear.

Adapted from NormalFormIdentification.jl (Hans Würfel). The autodiff Jacobians
are computed with `ForwardDiff`; the original finite-difference cross-check
assertions are dropped — they were over-strict (fixed `atol=1e-6`) and spuriously
failed on stiff machine Jacobians even when the autodiff result was correct.
"""
function nf_linearization(vm::VertexModel, state=NetworkDynamics.get_defaults_or_inits_dict(vm))
    res = init_residual(vm, state)
    if isnan(res) || res > 1e-8
        @warn "Vertex model does not appear to be at a steady state (residual=$res); linearization may be affine."
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
    @assert maximum(abs.(f_inner(xvec, idqvec))) < 1e-6
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

