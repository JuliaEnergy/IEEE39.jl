
using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using NetworkDynamics
using DataFrames
using CSV
using LinearAlgebra: diag, isdiag

@independent_variables t

DATA_DIR = joinpath(pkgdir(@__MODULE__), "ieee39data")

function get_IEEE39_base_ds(; add_killswitch=false)


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
        #println(i)
        row.has_load && apply_csv_params!(bus, load_df, i)
        row.has_gen && apply_csv_params!(bus, machine_df, i)
        row.has_avr && apply_csv_params!(bus, avr_df, i)
        row.has_gov && apply_csv_params!(bus, gov_df, i)

        ## Set power flow model based on bus type
        pf_model = if row.bus_type == "PQ"
            pfPQ(P=row.P, Q=row.Q)  ## Load bus: fixed P and Q
        elseif row.bus_type == "PV"
            # pfPV(P=row.P, V=row.V)  ## Generator bus: fixed P and V
            if i == 39
                dist_slack_follow_load_current_source(; P=row.P, V=row.V, idx=i)
                # pfPV(P=row.P, V=row.V)
            else
                dist_slack_follow_current_source(; P=row.P, V=row.V, idx=i)
            end
        elseif row.bus_type == "Slack"
            dist_slack_lead_current_source(; P=row.P, V=row.V, δ=0, idx=i)
            # pfSlack(V=row.V, δ=0)   ## Slack bus: fixed V and angle
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


@mtkmodel DistributedSlackLead begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        Pbase, [description="Base Power [pu]"]
        V, [description="Voltage Magnitude [pu]"]
        δ, [description="Voltage Angle [rad]"]
    end
    @variables begin
        P(t), [description="Active Power [pu]"]
        γ(t), [description="Deviation factor from base power"]
    end
    @equations begin
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        PowerDynamics.Library.@no_simplify P ~ Pbase + (1-γ)*abs(Pbase) # equation for gamma
        V^2 ~ terminal.u_r^2 + terminal.u_i^2 + implicit_output(terminal.i_r)
        δ ~ atan(terminal.u_i, terminal.u_r)  + implicit_output(terminal.i_i)
    end
end
function dist_slack_lead_current_source(; P=1, V=1, δ=0, idx)
    lead = DistributedSlackLead(Pbase=P, V=V, δ=δ; name=:slack_lead)
    # load = ZIPLoad(;name=:ZIPLoad)
    # load = Library.PQConstraint(; P, Q)

    mtkbus = MTKBus(lead)
    bus = compile_bus(mtkbus; vidx=idx, current_source=false, assume_io_coupling=false)
    set_default!(bus, r"γ$", 1)

    bus
end

@mtkmodel DistributedSlackFollow begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        Pbase, [description="Base Power [pu]"]
        V, [description="Voltage Magnitude [pu]"]
    end
    @variables begin
        P(t), [description="Active Power [pu]"]
        γ(t), [description="Deviation factor from base power",guess=1]
    end
    @equations begin
        P ~ Pbase + (1-γ)*abs(Pbase)
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        V^2 ~ terminal.u_r^2 + terminal.u_i^2
    end
end
function dist_slack_follow_current_source(; P=1, V=1, idx)
    follow = DistributedSlackFollow(Pbase=P, V=V; name=:slack_follow)
    mtkbus = MTKBus(follow)
    bus = compile_bus(
        mtkbus;
        vidx=idx,
        current_source=false,
        extin=[:slack_follow₊γ => VIndex(:31, :slack_lead₊γ)],
        assume_io_coupling=true,
    )
    bus
end

function dist_slack_follow_load_current_source(; P=1, V=1, idx)
    follow = DistributedSlackFollow(Pbase=P, V=V; name=:slack_follow)

    mtkbus = MTKBus(follow)
    bus = compile_bus(
        mtkbus;
        vidx=idx,
        current_source=false,
        extin=[:slack_follow₊γ => VIndex(:31, :slack_lead₊γ)],
        assume_io_coupling=true,
    )
    bus
end