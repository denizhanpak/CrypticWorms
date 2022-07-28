#This will generate a Probability Map For 
using DynamicalSystems
using NetworkDynamics
using Random
using Graphs
using DifferentialEquations
include("./cervera_cells.jl")
include("./plotting.jl")

Adj = [0 1 ; 1 0]
G = SimpleDiGraph(Adj)
cell = ODEVertex(; f=SingleCellVertex!, dim=1)
gap = StaticEdge(;f=GapJunctionEdge!,dim=1,coupling=:directed)

gref = 0.1
capacitance = 100
gpolarization = 4.0 * gref
epolarization = -70.0
gdepolarization = 0.5 * gref
edepolarization = 0.0
currentpump = 0.0
gmin = 0.2 * gref
gmax = 2.0 * gref #2.0 for cusp region

x0 = [0.0,0.0]

cell_p=(capacitance, gpolarization, epolarization, gdepolarization, edepolarization, currentpump)
gap_p = (gmin, gmax, 20.0, 2.0)
nd = network_dynamics(cell,gap, G)
ode_problem = ODEProblem(nd, x0, (0.0,20.0),(cell_p,gap_p))
sol = solve(ode_problem);