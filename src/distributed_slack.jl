
@independent_variables t

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
        PowerDynamics.Library.@no_simplify P ~ Pbase + (1-γ)*abs(Pbase)
        V^2 ~ terminal.u_r^2 + terminal.u_i^2 + implicit_output(terminal.i_r)
        δ ~ atan(terminal.u_i, terminal.u_r)  + implicit_output(terminal.i_i)
    end
end

"""
    dist_slack_lead_current_source(; P, V, δ=0, idx) -> VertexModel

Create the lead bus for a distributed-slack power-flow model.  The lead fixes
voltage magnitude `V` and angle `δ` and exposes the balance variable `γ` as an
output so follower buses can couple to it.
"""
function dist_slack_lead_current_source(; P=1, V=1, δ=0, idx)
    lead = DistributedSlackLead(Pbase=P, V=V, δ=δ; name=:slack_lead)
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
        γ(t), [description="Deviation factor from base power", guess=1]
    end
    @equations begin
        P ~ Pbase + (1-γ)*abs(Pbase)
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        V^2 ~ terminal.u_r^2 + terminal.u_i^2
    end
end

"""
    dist_slack_follow_current_source(; P, V, idx, lead_vidx) -> VertexModel

Create a follower bus for a distributed-slack power-flow model.  The follower
receives the balance variable `γ` from the lead bus (vertex index `lead_vidx`)
and adjusts its net injection proportionally.
"""
function dist_slack_follow_current_source(; P=1, V=1, idx, lead_vidx)
    follow = DistributedSlackFollow(Pbase=P, V=V; name=:slack_follow)
    mtkbus = MTKBus(follow)
    bus = compile_bus(
        mtkbus;
        vidx=idx,
        current_source=false,
        extin=[:slack_follow₊γ => VIndex(lead_vidx, :slack_lead₊γ)],
        assume_io_coupling=true,
    )
    bus
end
