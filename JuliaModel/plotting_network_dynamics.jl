#This will generate a Probability Map For 
using DynamicalSystems
using NetworkDynamics
using Random
using Graphs
using DifferentialEquations
using ChaosTools
using Printf
using GraphPlot
include("./cervera_cells.jl")
include("./plotting.jl")

function make_body(G)
    cell = ODEVertex(; f=SingleCellVertex!, dim=1)
    gap = StaticEdge(;f=GapJunctionEdge!,dim=1,coupling=:directed)
    nd = network_dynamics(cell, gap, G)
    return nd
end

function scatter_attractors!(ax, attractors)
    for k ∈ keys(attractors)
        x, y = columns(attractors[k])
        scatter!(ax, attractors[k].data;
            color = Cycled(k),
            strokewidth = 3, strokecolor = :white
        )
    end
end

function body_dynamics(network, gmax, gpol_head,x,y)
    gref = 0.1
    capacitance = 100
    gpolarization = 4.0 * gref
    epolarization = -70.0
    gdepolarization = 0.5 * gref
    edepolarization = 0.0
    currentpump = 0.0
    gmin = 0.2 * gref
    gmax = gmax * gref #2.0 for cusp region
    dims = x * y


    x0 = zeros(Float64,dims)

    cell_tail=(capacitance, gpolarization, epolarization, gdepolarization, edepolarization, currentpump)
    cell_head=(capacitance, gpol_head * gref, epolarization, gdepolarization, edepolarization, currentpump)
    params = []
    for i in 1:y
        for j in 1:x
            if j <= x/2
                push!(params,cell_head)
            else
                push!(params,cell_tail)
            end
        end
    end

    cell_p = ([cell_tail,cell_head])
    gap_p = (gmin, gmax, 20.0, 2.0)
    ds = ContinuousDynamicalSystem(network, x0,(params,gap_p))
    #ds = ContinuousDynamicalSystem(network, x0,(cell_p,gap_p))
    return ds
end

fs = []
ams = []
Adj = [0 1 ; 1 0]
G = SimpleDiGraph(Adj)
x = 8
y = 1
G = DiGraph(Graphs.grid([x,y]))
gpol_range = range(0.0,0.2,length=2)
gmax_range = range(0.0,0.5,length=2)
body = make_body(G)
gpol_range = [0.1]
gmax_range = [0.1]
plots = []
for gpol in gpol_range
    for gmax in gmax_range
        initial = @SVector zeros(Float64, x*y)

        ds = body_dynamics(body, gmax, gpol, x, y)
        global tr = trajectory(ds, 7000)
        state = (tr[size(tr)[1],:] ./ 70) .+ 1
        global c = get_colors(state,colorant"red",colorant"blue")
        display(gplot(G,nodefillc=c,layout=grid_layout_generic(x,y)))
        p = plot(1:size(tr)[1],tr[:,1])
        for i in 2:x*y
            plot!(1:size(tr)[1],tr[:,i])
        end
        display(p)
    end
end