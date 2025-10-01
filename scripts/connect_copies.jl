
using Pkg
Pkg.activate(@__DIR__)

using PowerDynamics
using Graphs
using OrdinaryDiffEq
using Plots
using IEEE39
using SparseConnectivityTracer
using Logging
# Logging.disable_logging(Logging.Warn)

##
nw_base = get_IEEE39_base()
nw_39 = set_IEEE39_PF_init(nw_base)


##

add_subnet_name(i, s) = Symbol("$(i)_$s")

make_subnet_VM(vm, sni, n_vert) = VertexModel(copy(vm), vidx=(n_vert * sni) + get_graphelement(vm), name=add_subnet_name("net$sni", vm.name))

function make_subnet_EM(em, sni, n_vert)
    EM_edge = get_graphelement(em)
    EdgeModel(copy(em), src=EM_edge.src + (n_vert * sni), dst = EM_edge.dst + (n_vert * sni))
end

##

n_vert = length(nw_39[VIndex(:)])

nw_vertices = []
nw_edges = []
for sni in 0:30
    println(n_vert * sni)
    append!(nw_vertices, [make_subnet_VM(vm, sni, n_vert) for vm in nw_39[VIndex(:)]])
    append!(nw_edges, [make_subnet_EM(em, sni, n_vert) for em in nw_39[EIndex(:)]])
end

for sni in 0:29
    em = rand(nw_39[EIndex(:)])
    em_edge = get_graphelement(em)
    append!(nw_edges, [EdgeModel(copy(em), src=em_edge.src + (n_vert * sni), dst = em_edge.dst + (n_vert * (sni + 1)))])
end

##

nw_large = Network(nw_vertices, nw_edges)

##

nws_init = initialize_from_pf(nw_large)

##

NetworkDynamics.set_jac_prototype!(nw_large; remove_conditions=true)

##

prob = ODEProblem(nw_large, copy(uflat(nws_init)), (0.0, 5.0), copy(pflat(nws_init)); callback=sc_and_trip(1, 0.1, 0.15, 0.1))
tol=1e-3

##

t_start=time()
sol = solve(prob, FBDF(), maxiters=10000, abstol=tol, reltol=tol, dtmin = 1e-4, force_dtmin = true)
t_end = time()
println(t_end-t_start)

##
using CUDA
using LinearSolve
t_start=time()
sol = solve(prob, KenCarp47(; linsolve = KLUFactorization()), maxiters=10000, abstol=tol, reltol=tol, dtmin = 1e-4, force_dtmin = true)
t_end = time()
println(t_end-t_start)


##

t_start=time()
sol = solve(prob, Rodas4P2(), maxiters=10000, abstol=tol, reltol=tol, dtmin = 1e-4, force_dtmin = true)
t_end = time()
println(t_end-t_start)

##

t_start=time()
sol = solve(prob, Rosenbrock23(), maxiters=10000, abstol=tol, reltol=tol, dtmin = 1e-4, force_dtmin = true)
t_end = time()
println(t_end-t_start)

##

t_start=time()
sol = solve(prob, QNDF(), maxiters=10000, abstol=tol, reltol=tol, dtmin = 1e-4, force_dtmin = true)
t_end = time()
println(t_end-t_start)

##

t_start=time()
sol = solve(prob, QNDF(; linsolve = KLUFactorization()), maxiters=10000, abstol=tol, reltol=tol, dtmin = 1e-4, force_dtmin = true)
t_end = time()
println(t_end-t_start)

##

t_start=time()
sol = solve(prob, QNDF(; linsolve = CudaOffloadLUFactorization()), maxiters=10000, abstol=tol, reltol=tol, dtmin = 1e-4, force_dtmin = true)
t_end = time()
println(t_end-t_start)

##
using CUDA
using Sparspak
using CUSOLVERRF

nw_cu = CUDA.adapt(CuArray{Float32}, Network(nw_large; execution=PolyesterExecution{true}(), aggregator=PolyesterAggregator(+)))
prob = ODEProblem(nw_cu, copy(uflat(nws_init)) |> cu, (0.0, 5.0), copy(pflat(nws_init)) |> cu; callback=sc_and_trip(1, 0.1, 0.15, 0.1))

##
CUDA.allowscalar()
t_start=time()
sol = solve(prob, QNDF(; linsolve=CUSOLVERRFFactorization()), maxiters=10000, abstol=tol, reltol=tol, dtmin = 1e-4, force_dtmin = true)
t_end = time()
println(t_end-t_start)

##

plot(sol, idxs=vidxs(nw_large, 1:39, r"machine₊δ$"))