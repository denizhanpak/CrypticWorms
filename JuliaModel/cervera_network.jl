#This will generate a Probability Map For 
using DynamicalSystems
using NetworkDynamics
using Random
using Graphs
using DifferentialEquations
include("./cervera_cells.jl")
include("./plotting.jl")

function make_body(G)
    cell = ODEVertex(; f=SingleCellVertex!, dim=1)
    gap = StaticEdge(;f=GapJunctionEdge!,dim=1,coupling=:directed)
    nd = network_dynamics(cell, gap, G)
    return nd
end

function body_dynamics(network, gmax, gpol_head)
    gref = 0.1
    capacitance = 100
    gpolarization = 4.0 * gref
    epolarization = -70.0
    gdepolarization = 0.5 * gref
    edepolarization = 0.0
    currentpump = 0.0
    gmin = 0.2 * gref
    gmax = gmax * gref #2.0 for cusp region

    x0 = [0.0,0.0]

    cell_tail=(capacitance, gpolarization, epolarization, gdepolarization, edepolarization, currentpump)
    cell_head=(capacitance, gpol_head * gref, epolarization, gdepolarization, edepolarization, currentpump)
    cell_p = ([cell_tail,cell_head])
    gap_p = (gmin, gmax, 20.0, 2.0)
    ds = ContinuousDynamicalSystem(nd, x0,(cell_p,gap_p))
    return ds
end

fs = []
ams = []
Adj = [0 1 ; 1 0]
G = SimpleDiGraph(Adj)
gpol_range = range(0,0.3,length=5)
gmax_range = range(0.0,2.0,length=5)
body = make_body(G)

for gpol in gpol_range
    for gmax in gmax_range
        initial = SVector(0.0,0.0)

        ds = body_dynamics(body, gmax, gpol)
        xg = yg = range(-70, 70; length = 70)
        grid = (xg, yg)
        am = AttractorsViaRecurrences(ds, grid)

        rng = Random.MersenneTwister(1234)
        sampler, _ = statespace_sampler(rng; min_bounds=[-70., -70.], max_bounds=[10., 70.])
        push!(fs, basins_fractions(am, sampler; N = 1000, show_progress = false))
        basin, attractors = basins_of_attraction(am)
        if length(attractors) == 2
            print(attractors[1][1],"\n")
            print(attractors[2][1],"\n")
        end
    end
end