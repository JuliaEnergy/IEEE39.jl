
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

function get_IEEE39_base(; add_killswitch=false)


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
            pfPQ(P=row.P, Q=row.Q)  ## Load bus: fixed P and Q
        elseif row.bus_type == "PV"
            pfPV(P=row.P, V=row.V)  ## Generator bus: fixed P and V
        elseif row.bus_type == "Slack"
            pfSlack(V=row.V, δ=0)   ## Slack bus: fixed V and angle
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

