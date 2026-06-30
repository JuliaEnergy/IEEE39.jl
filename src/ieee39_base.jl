
using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using NetworkDynamics
using DataFrames
using CSV
using LinearAlgebra: diag, isdiag

DATA_DIR = joinpath(pkgdir(@__MODULE__), "ieee39data")

function apply_csv_params!(bus, table, bus_index)
    row_idx = findfirst(table.bus .== bus_index)

    ## Apply all parameters except "bus" column
    row = table[row_idx, :]
    for col_name in names(table)
        if col_name != "bus"
            set_default!(bus, Regex(col_name*"\$"), row[col_name])
        end
    end
end

struct KillswitchWrapper{F}
    f::F
    odes::Vector{Int}
end

function (kw::KillswitchWrapper{F})(du, u, ein, p, t) where {F}
    kw.f(du, u, ein, p, t)
    killswitch = p[end]
    if !iszero(killswitch)
        for i in kw.odes
            du[i] = 0
        end
    end
    nothing
end

function add_killswitch_to_vm(vm)
    @assert isdiag(vm.mass_matrix)
    odes = findall(!iszero, diag(vm.mass_matrix))
    kwf = KillswitchWrapper(vm.f, odes)
    newp = vcat(psym(vm), :killswitch)

    vm_kill = VertexModel(vm, f=kwf, psym=newp)
    set_default!(vm_kill, :killswitch, 0)

    vm_kill
end

function get_MATPOWER_LineRatings_MVA()
    mp_branch_df = CSV.read(joinpath(DATA_DIR, "mp_branch.csv"), DataFrame)
    mp_branch_df.rateA
end

function get_MATPOWER_LineRatings_pu(BASE_MVA = 100.0)
    mp_branch_df = CSV.read(joinpath(DATA_DIR, "mp_branch.csv"), DataFrame)
    mp_branch_df.rateA ./ BASE_MVA
end

function get_IEEE39_base(; add_killswitch=false, distributed_slack=false)


    branch_df = CSV.read(joinpath(DATA_DIR, "branch.csv"), DataFrame)
    bus_df = CSV.read(joinpath(DATA_DIR, "bus.csv"), DataFrame)
    load_df = CSV.read(joinpath(DATA_DIR, "load.csv"), DataFrame)
    machine_df = CSV.read(joinpath(DATA_DIR, "machine.csv"), DataFrame)
    avr_df = CSV.read(joinpath(DATA_DIR, "avr.csv"), DataFrame)
    gov_df = CSV.read(joinpath(DATA_DIR, "gov.csv"), DataFrame)

    BASE_MVA = 100.0
    BASE_FREQ = 60.0


    load = ZIPLoad(;name=:ZIPLoad)

    uncontrolled_machine = SauerPaiMachine(;
        τ_m_input=false,  ## No external mechanical torque input
        vf_input=false,   ## No external field voltage input
        name=:machine,
    )

    _machine = SauerPaiMachine(;
        name=:machine,
    )
    _avr = AVRTypeI(;
        name=:avr,
        ceiling_function=:quadratic,
    )
    _gov = TGOV1(; name=:gov,)

    controlled_machine = CompositeInjector(
        [_machine, _avr, _gov],
        name=:ctrld_gen
    )

    @named junction_bus_template = compile_bus(MTKBus())
    strip_defaults!(junction_bus_template)  ## Clear default parameters for manual setting

    @named load_bus_template = compile_bus(MTKBus(load))
    strip_defaults!(load_bus_template)

    @named ctrld_machine_bus_template = compile_bus(
        MTKBus(controlled_machine);
    )
    strip_defaults!(ctrld_machine_bus_template)

    set_default!(ctrld_machine_bus_template, r"S_b$", BASE_MVA)
    set_default!(ctrld_machine_bus_template, r"ω_b$", 2π*BASE_FREQ)

    @named ctrld_machine_load_bus_template = compile_bus(
        MTKBus(controlled_machine, load);
    )
    strip_defaults!(ctrld_machine_load_bus_template)
    set_default!(ctrld_machine_load_bus_template, r"S_b$", BASE_MVA)
    set_default!(ctrld_machine_load_bus_template, r"ω_b$", 2π*BASE_FREQ)

    @named unctrld_machine_load_bus_template = compile_bus(
        MTKBus(uncontrolled_machine, load);
    )
    strip_defaults!(unctrld_machine_load_bus_template)
    set_default!(unctrld_machine_load_bus_template, r"S_b$", BASE_MVA)
    set_default!(unctrld_machine_load_bus_template, r"ω_b$", 2π*BASE_FREQ)


    busses = []
    for row in eachrow(bus_df)
        i = row.bus

        ## Select template based on bus category
        bus = if row.category == "junction"
            compile_bus(junction_bus_template; vidx=i, name=Symbol("bus$i"))
        elseif row.category == "load"
            compile_bus(load_bus_template; vidx=i, name=Symbol("bus$i"))
        elseif row.category == "ctrld_machine"
            compile_bus(ctrld_machine_bus_template; vidx=i, name=Symbol("bus$i"))
        elseif row.category == "ctrld_machine_load"
            compile_bus(ctrld_machine_load_bus_template; vidx=i, name=Symbol("bus$i"))
        elseif row.category == "unctrld_machine_load"
            compile_bus(unctrld_machine_load_bus_template; vidx=i, name=Symbol("bus$i"))
        end

        ## Apply component parameters from CSV files
        row.has_load && apply_csv_params!(bus, load_df, i)
        row.has_gen && apply_csv_params!(bus, machine_df, i)
        row.has_avr && apply_csv_params!(bus, avr_df, i)
        row.has_gov && apply_csv_params!(bus, gov_df, i)

        ## Set power flow model based on bus type
        pf_model = if row.bus_type == "PQ"
            pfPQ(P=row.P, Q=row.Q)
        elseif row.bus_type == "PV"
            if distributed_slack
                dist_slack_follow_current_source(P=row.P, V=row.V, idx=i, lead_vidx=31)
            else
                pfPV(P=row.P, V=row.V)
            end
        elseif row.bus_type == "Slack"
            if distributed_slack
                dist_slack_lead_current_source(P=row.P, V=row.V, δ=0, idx=i)
            else
                pfSlack(V=row.V, δ=0)
            end
        end
        set_pfmodel!(bus, pf_model)

        push!(busses, bus)
    end

    @named piline_template = compile_line(MTKLine(PiLine_fault(;name=:piline)))

    branches = []
    for row in eachrow(branch_df)
        ## Create line instance with topology
        line = compile_line(piline_template; src=row.src_bus, dst=row.dst_bus)

        ## Apply electrical parameters from CSV data
        for col_name in names(branch_df)
            if col_name ∉ ["src_bus", "dst_bus", "transformer"]
                set_default!(line, Regex(col_name*"\$"), row[col_name])
            end
        end

        push!(branches, line)
    end

    if add_killswitch
        busses = add_killswitch_to_vm.(busses)
    end

    Network(busses, branches)
end



function set_IEEE39_PF_init(nw)

    

    for bus_idx in 1:39
        if :ZIPLoad₊Vset ∉ nw[VIndex(bus_idx)].psym
            continue
        end

        vi = VIndex(bus_idx) 
        
        # To access parts of network objects we need special indices that distinguish edges/lines and vertices/buses

        set_default!(nw[vi], :ZIPLoad₊Vset, 1.)


        # Bus 31 and 39 contain a load and a generator. 
        # We keep the power and reactive power set points intact, 
        # variations in the powerflow are assigned to the generator during initialization.

        if bus_idx in [31,39]
            continue
        end

        # For all other loads we unset the defaults,
        # P and Q are then determined during initialization
        # from the provided power flow solution

        P = get_default(nw[vi], :ZIPLoad₊Pset)
        Q = get_default(nw[vi], :ZIPLoad₊Qset)
        delete_default!(nw[vi], :ZIPLoad₊Pset)
        delete_default!(nw[vi], :ZIPLoad₊Qset)
        set_guess!(nw[vi], :ZIPLoad₊Pset, P)
        set_guess!(nw[vi], :ZIPLoad₊Qset, Q)
    end

    for ei in 1:ne(nw)
        eidx = EIndex(ei)
        f = copy_pf_parameters(nw[eidx])
        add_pfinitformula!(nw[eidx], f)
    end

    nw
end

"""
    load_ieee39_scenario(row; pert=0.0, rng=Random.default_rng(), BASE_MVA=100.0)
    -> (load_P::Dict{Int,Float64}, load_Q::Dict{Int,Float64})

Return active and reactive power demand dictionaries for one 15-minute time step
from the IEEE39 load scenarios published in:

Tricarico, Gioacchino, et al. "A modified version of the IEEE 39-bus test system for the day-ahead market." 2023 IEEE PES Conference on Innovative Smart Grid Technologies-Middle East (ISGT Middle East). IEEE, 2023.
https://ieeexplore.ieee.org/document/10078548

`row` selects a 15-minute time step (1-indexed; the files cover 96 × 365 = 35040
steps for the year 2014).

The returned dicts map bus index → per-unit power injection (negative = load
consuming power, 100 MVA base).  Active power columns are in MW and reactive
power columns are labelled as *inductive* (absorbed from the network), so both
are negated when converting to per-unit injections:

    load_P[b] = -(P_MW[b] / BASE_MVA)
    load_Q[b] = -(Q_inductive_MVAr[b] / BASE_MVA)

These dicts can be passed directly to [`generate_powerflow_variation`](@ref) or
[`generate_gog_powerflow_variation`](@ref).

If `pert > 0`, each P and Q value is independently multiplied by
`(1 + pert * randn(rng))` to add load uncertainty.
"""
function load_ieee39_scenario(row; pert=0.0, rng=Random.default_rng(), BASE_MVA=100.0)
    p_file = joinpath(DATA_DIR, "Loads_P.csv")
    q_file = joinpath(DATA_DIR, "Loads_Q.csv")

    # P: semicolon-delimited; row 1 = names (first col is timestamp), row 2 = units, rows 3+ = data
    df_P_raw = CSV.read(p_file, DataFrame; delim=';', header=1, skipto=3)
    df_P = df_P_raw[:, 2:end]   # drop timestamp column ("Quasi-Dynamic Simulation AC")

    # Q: comma-delimited; row 1 = names, row 2 = units, rows 3+ = data
    df_Q = CSV.read(q_file, DataFrame; delim=',', header=1, skipto=3)

    n = nrow(df_P)
    1 <= row <= n || throw(ArgumentError("row must be in 1:$n, got $row"))

    bus_of(col) = parse(Int, match(r"\d+", string(col)).match)

    function make_dict(df, row_vals)
        d = Dict{Int, Float64}()
        for col in names(df)
            b = bus_of(col)
            v = -Float64(row_vals[col]) / BASE_MVA   # negate: file stores consumption
            pert > 0 && (v *= 1.0 + pert * randn(rng))
            d[b] = v
        end
        d
    end

    load_P = make_dict(df_P, df_P[row, :])
    load_Q = make_dict(df_Q, df_Q[row, :])

    load_P, load_Q
end

