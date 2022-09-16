#This will generate a Probability Map For 
using DynamicalSystems
using NetworkDynamics
using Random
using Graphs
using DifferentialEquations
using ChaosTools
using Printf
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

function body_dynamics(network, gmax, gpol_head, x0=[0.0,0.0])
    gref = 0.1
    capacitance = 100
    gpolarization = 4.0 * gref
    epolarization = -70.0
    gdepolarization = 0.5 * gref
    edepolarization = 0.0
    currentpump = 0.0
    gmin = 0.2 * gref
    gmax = gmax * gref #2.0 for cusp region


    cell_tail=(capacitance, gpolarization, epolarization, gdepolarization, edepolarization, currentpump)
    cell_head=(capacitance, gpol_head * gref, epolarization, gdepolarization, edepolarization, currentpump)
    cell_p = ([cell_tail,cell_head])
    gap_p = (gmin, gmax, 20.0, 2.0)
    ds = ContinuousDynamicalSystem(network, x0,(cell_p,gap_p))
    return ds
end

function attractor_map()
    fs = []
    ams = []
    Adj = [0 1 ; 1 0]
    G = SimpleDiGraph(Adj)
    gpol_range = range(0.27,0.28,length=20)
    gmax_range = range(1.19,2.0,length=50)
    body = make_body(G)
    gpol_range = [0.11]
    gmax_range = [0.5]
    for gpol in gpol_range
        for gmax in gmax_range
            initial = SVector(0.0,0.0)

            ds = body_dynamics(body, gmax, gpol)
            tr = trajectory(ds, 1000)
            print(tr)
            break
            xg = yg = range(-70.0, 0.0; length = 1000)
            grid = (xg, yg)
            am = AttractorsViaRecurrences(ds, grid)

            rng = Random.MersenneTwister(1234)
            sampler, _ = statespace_sampler(rng; min_bounds=[-70., -70.], max_bounds=[0., 0.])
            push!(fs, basins_fractions(am, sampler; N = 1000, show_progress = false))
            basins, attractors = basins_of_attraction(am)
            if length(attractors) == 2
                tmp = attractors
                f1 = attractors[1][1][1] - attractors[1][1][2]
                f2 = attractors[2][1][1] - attractors[2][1][2]
                if abs(f1) > abs(f2)
                    normal = 1
                    cryptic = 2
                else
                    normal = 2
                    cryptic = 1
                end
                a1 = attractors[normal][1]
                a2 = attractors[cryptic][1]
                ax = [a1[1],a2[1]]
                ay = [a1[2],a2[2]]
            else
                ax = [attractors[1][1][1]]
                ay = [attractors[1][1][2]]
            end

            fig = plot(heatmap(xg,yg,basins,c=[:cyan,:magenta]))
            fig = plot(scatter!(ay,ax,c=[:blue,:red]))

            savefig(fig,"./phase_diagrams/gpol_$(@sprintf("%.4f", gpol)).png")
        end
    end
end


