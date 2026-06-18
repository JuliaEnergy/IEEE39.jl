using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PowerDynamics
using IEEE39
using DataFrames
using CSV

# Load the base network and initialize power flow
nw_base = get_IEEE39_base()
nw_39 = set_IEEE39_PF_init(nw_base)
pfnw_39 = powerflow_model(nw_39)
pfs0_39 = NWState(pfnw_39)

println("="^70)
println("POWER FLOW DIAGNOSTICS FOR IEEE 39-BUS SYSTEM")
println("="^70)

# Get bus and load data
DATA_DIR = joinpath(pkgdir(IEEE39), "ieee39data")
bus_df = CSV.read(joinpath(DATA_DIR, "bus.csv"), DataFrame)
load_df = CSV.read(joinpath(DATA_DIR, "load.csv"), DataFrame)

# Calculate total power balance from bus.csv
println("\n1. INITIAL POWER BALANCE (from bus.csv):")
println("-"^70)

total_load_P = sum(bus_df[bus_df.category .== "load", :P])
total_load_Q = sum(bus_df[bus_df.category .== "load", :Q])
total_gen_P = sum(bus_df[bus_df.category .!= "load" .&& bus_df.category .!= "junction", :P])

# Account for buses with both load and generation
mixed_P = bus_df[(bus_df.has_gen .== true) .&& (bus_df.has_load .== true), :P]
mixed_Q = bus_df[(bus_df.has_gen .== true) .&& (bus_df.has_load .== true), :Q]

println("Load Buses P:                 $(round(total_load_P; digits=2)) MW")
println("Load Buses Q:                 $(round(total_load_Q; digits=2)) MVAr")
println("Generator Buses P:            $(round(total_gen_P; digits=2)) MW")
println("Mixed Gen+Load Buses:         $(length(mixed_P)) buses")
if length(mixed_P) > 0
    println("  Mixed P:                  $(round(sum(mixed_P); digits=2)) MW")
    println("  Mixed Q:                  $(round(sum(mixed_Q); digits=2)) MVAr")
end

# Calculate imbalance
imbalance_P = total_gen_P + sum(mixed_P) - total_load_P
imbalance_Q = sum(bus_df.Q) - total_load_Q

println("\nPower Balance Summary:")
println("  Total Generation:           $(round(total_gen_P + sum(mixed_P); digits=2)) MW")
println("  Total Load:                 $(round(total_load_P; digits=2)) MW")
println("  ⚠ IMBALANCE (ΔP):          $(round(imbalance_P; digits=2)) MW")

if abs(imbalance_P) > 1.0
    println("\n  *** CRITICAL: Active power mismatch of $(round(abs(imbalance_P); digits=2)) MW ***")
    println("  This MUST be fixed before power flow can converge!")
end

println("\n2. ATTEMPTING POWER FLOW SOLUTION:")
println("-"^70)

try
    pfs = solve_powerflow(nw_39; pfnw=pfnw_39, pfs0=pfs0_39, verbose=true, maxiter=100)
    println("\n✓ Power flow converged successfully!")

    println("\nBus Voltage Summary (first 10 buses):")
    for i in 1:min(10, length(pfs.v))
        V = abs(pfs.v[i].u[1] + 1im*pfs.v[i].u[2])
        δ = atan(pfs.v[i].u[2], pfs.v[i].u[1])
        bus_type = bus_df[i, :bus_type]
        println("  Bus $i ($bus_type): V = $(round(V; digits=4)) pu, δ = $(round(rad2deg(δ); digits=2))°")
    end

catch e
    println("✗ Power flow FAILED to converge")
    println("Error: $e")

    println("\n3. DIAGNOSIS & RECOMMENDATIONS:")
    println("-"^70)

    if abs(imbalance_P) > 1.0
        println("\n ★ PRIMARY ISSUE: Active Power Imbalance")
        println("━"^70)
        println("The power flow cannot converge because generation ≠ load.")
        println("Current imbalance: $(round(imbalance_P; digits=2)) MW")
        println("\n FIX OPTIONS:")
        if imbalance_P > 0
            println("  Excess generation of $(round(imbalance_P; digits=2)) MW needs to be absorbed:")
            println("  1. Increase loads by $(round(imbalance_P; digits=2)) MW")
            println("  2. Reduce generation by $(round(imbalance_P; digits=2)) MW")
            println("  3. Reduce generator P values at buses 30-38")
        else
            deficit = abs(imbalance_P)
            println("  Generation deficit of $deficit MW needs additional supply:")
            println("  1. Increase generation by $deficit MW")
            println("  2. Reduce loads by $deficit MW")
            println("  3. Increase Slack bus (bus 31) generation")
        end
    else
        println("\n ★ SECONDARY ISSUE: Reactive Power or Convergence Problem")
        println("━"^70)
        println("Power balance is OK, but other issues prevent convergence:")
        println("  • Reactive power may be out of generator limits")
        println("  • Bus voltages may be violating physical constraints")
        println("  • Network may have numerical issues")
    end

    println("\n ★ RECOMMENDED FIXES (in order):")
    println("━"^70)
    println("1. BALANCE ACTIVE POWER FIRST")
    println("   - Make sure ΣP_gen = ΣP_load")
    println("   - Start without random perturbations")
    println("")
    println("2. VERIFY REACTIVE POWER LIMITS")
    println("   - Each generator has Q_min and Q_max")
    println("   - Check that Q values are within bounds")
    println("")
    println("3. REDUCE PERTURBATIONS")
    println("   - Set perturbation factor to 0% initially")
    println("   - Verify baseline solution works")
    println("   - Then gradually increase perturbation (e.g., 5% → 10%)")
    println("")
    println("4. CHECK VOLTAGE BOUNDS")
    println("   - PV buses should have V targets in range [0.9-1.1] pu")
    println("   - Slack bus voltage should be realistic")
end

println("\n" * "="^70)
