using DiffEqCallbacks
using SciMLBase

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

