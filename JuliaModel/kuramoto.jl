using NetworkDynamics, OrdinaryDiffEq, Plots, Graphs


N = 8
g = watts_strogatz(N, 2, 0) # ring network

function kuramoto_edge!(e, θ_s, θ_d, K, t)
    e[1] = K * sin(θ_s[1] - θ_d[1])
end

function kuramoto_vertex!(dθ, θ, edges, ω, t)
    dθ[1] = ω
    sum_coupling!(dθ, edges)
end

vertex! = ODEVertex(; f=kuramoto_vertex!, dim=1, sym=[:θ])
edge!   = StaticEdge(; f=kuramoto_edge!, dim=1)
nd! = network_dynamics(vertex!, edge!, g);