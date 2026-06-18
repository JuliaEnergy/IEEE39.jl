using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DataFrames
using CSV

DATA_DIR = joinpath(dirname(@__DIR__), "ieee39data")

# Read data files
bus_df = CSV.read(joinpath(DATA_DIR, "bus.csv"), DataFrame)
load_df = CSV.read(joinpath(DATA_DIR, "load.csv"), DataFrame)

println("="^70)
println("IEEE 39-BUS SYSTEM - POWER BALANCE ANALYSIS")
println("="^70)

println("\n1. BUS DATA SUMMARY:")
println("-"^70)
println("Bus Type Distribution:")
for bus_type in unique(bus_df.bus_type)
    count = sum(bus_df.bus_type .== bus_type)
    println("  $bus_type: $count buses")
end

println("\n2. POWER BALANCE:")
println("-"^70)

# Sum all P values (loads are negative, generators positive)
total_P = sum(bus_df.P)
total_Q = sum(bus_df.Q)

println("Sum of all bus P values:     $(round(total_P; digits=3)) MW")
println("Sum of all bus Q values:     $(round(total_Q; digits=3)) MVAr")

# Breakdown by bus type
println("\nDetailed Breakdown:")
for cat in ["load", "junction", "ctrld_machine", "ctrld_machine_load", "unctrld_machine_load"]
    mask = bus_df.category .== cat
    if any(mask)
        P_sum = sum(bus_df[mask, :P])
        Q_sum = sum(bus_df[mask, :Q])
        count = sum(mask)
        println("  $cat ($count buses):")
        println("    P = $(round(P_sum; digits=3)) MW, Q = $(round(Q_sum; digits=3)) MVAr")
        # Show bus numbers
        bus_nums = bus_df[mask, :bus]
        println("    Buses: $(join(bus_nums, ", "))")
    end
end

# Analysis
println("\n3. ANALYSIS:")
println("-"^70)

if abs(total_P) < 0.1
    println("✓ Active power is BALANCED (ΣP ≈ 0 MW)")
else
    println("✗ Active power IMBALANCE of $(round(abs(total_P); digits=2)) MW")
    println("  This is why the power flow cannot converge!")
end

# Load breakdown
load_buses = bus_df[bus_df.category .== "load", :]
if nrow(load_buses) > 0
    total_load_P = sum(load_buses.P)
    println("\nLoad Buses ($(nrow(load_buses)) buses):")
    for row in eachrow(load_buses)
        println("  Bus $(Int(row.bus)): P = $(round(row.P; digits=2)) MW, Q = $(round(row.Q; digits=2)) MVAr")
    end
    println("  TOTAL: P = $(round(total_load_P; digits=2)) MW")
end

# Generator breakdown
gen_buses = bus_df[(bus_df.has_gen .== true) .&& (bus_df.category .!= "load"), :]
if nrow(gen_buses) > 0
    total_gen_P = sum(gen_buses.P)
    println("\nGenerator Buses ($(nrow(gen_buses)) buses):")
    for row in eachrow(gen_buses)
        println("  Bus $(Int(row.bus)) ($(row.bus_type)): P = $(round(row.P; digits=2)) MW, Q = $(round(row.Q; digits=2)) MVAr")
    end
    println("  TOTAL: P = $(round(total_gen_P; digits=2)) MW")
end

# Mixed buses
mixed_buses = bus_df[(bus_df.has_gen .== true) .&& (bus_df.has_load .== true), :]
if nrow(mixed_buses) > 0
    total_mixed_P = sum(mixed_buses.P)
    println("\nMixed Gen+Load Buses ($(nrow(mixed_buses)) buses):")
    for row in eachrow(mixed_buses)
        println("  Bus $(Int(row.bus)): P = $(round(row.P; digits=2)) MW (net = gen - load)")
    end
    println("  TOTAL: P = $(round(total_mixed_P; digits=2)) MW")
end

println("\n4. RECOMMENDATIONS:")
println("-"^70)
if abs(total_P) > 0.1
    println("\n★ FIX THE POWER IMBALANCE FIRST ★\n")
    if total_P > 0
        excess = total_P
        println("Current excess generation: $excess MW")
        println("\nOptions to balance:")
        println("  1. Reduce generator P values at PV buses (30-39)")
        println("  2. Increase load P values at load buses")
        println("  3. Adjust slack bus (31) settings")
    else
        deficit = abs(total_P)
        println("Current generation deficit: $deficit MW")
        println("\nOptions to balance:")
        println("  1. Increase generator P values at PV buses (30-39)")
        println("  2. Reduce load P values at load buses")
        println("  3. Increase slack bus (31) capacity")
    end
else
    println("✓ Power balance is OK!")
    println("\nIf power flow still fails to converge, investigate:")
    println("  • Reactive power violations (Q limits)")
    println("  • Voltage instability at certain buses")
    println("  • Reactive power perturbations too large")
end

println("\n" * "="^70)
