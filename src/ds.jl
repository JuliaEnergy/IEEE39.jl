using ModelingToolkit

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
        γ(t), [guess=1, description="Deviation factor from base power"]
    end
    @equations begin
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        PowerDynamics.Library.@no_simplify P ~ Pbase + (1-γ)*abs(Pbase) # equation for gamma
        V^2 ~ terminal.u_r^2 + terminal.u_i^2 + implicit_output(terminal.i_r)
        δ ~ atan(terminal.u_i, terminal.u_r)  + implicit_output(terminal.i_i)
    end
end
function dist_slack_lead_current_source(; P, V, δ)
    lead = DistributedSlackLead(Pbase=P, V=V, δ=δ; name=:slack_lead)
    mtkbus = MTKBus(lead)
    bus = compile_bus(mtkbus; current_source=true, assume_io_coupling=false)
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
        γ(t), [description="Deviation factor from base power"]
    end
    @equations begin
        P ~ Pbase + (1-γ)*abs(Pbase)
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        V^2 ~ terminal.u_r^2 + terminal.u_i^2
    end
end
function dist_slack_follow_current_source(; P, V)
    follow = DistributedSlackFollow(Pbase=P, V=V; name=:slack_follow)
    mtkbus = MTKBus(follow)
    bus = compile_bus(
        mtkbus;
        current_source=true,
        extin=[:slack_follow₊γ => VIndex(:undefined, :slack_lead₊γ)],
        assume_io_coupling=true,
    )
    bus
end