using NetworkDynamics

Base.@propagate_inbounds function fhn_electrical_vertex!(dv, v, edges, p, t)
    dv[1] = v[1] - v[1]^3 / 3 - v[2]
    dv[2] = (v[1] - a) * ϵ
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

Base.@propagate_inbounds function electrical_edge!(e, v_s, v_d, p, t)
    e[1] = p * (v_s[1] - v_d[1]) # * σ
    nothing
end

odeelevertex = ODEVertex(; f=fhn_electrical_vertex!, dim=2, sym=[:u, :v]);
electricaledge = StaticEdge(; f=electrical_edge!, dim=1, coupling=:directed)

fhn_network! = network_dynamics(odeelevertex, electricaledge, g_directed)