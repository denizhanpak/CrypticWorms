using Pkg
using DifferentialEquations
using Plots
using DynamicalSystems # also exports relevant StaticArrays names
using StaticArrays
theme(:dark)

include("./plotting.jl")
include("./cervera_cells.jl")


#Feeding the inputs to the solver 
#prob = ODEProblem(TwoCellCervera!, u0, tspan, p)
#sol = solve(prob)

#based on graphs setting gpol for Head to 0.011 gets difference to roughly 50mV

#=
for i = 0.01:0.001:0.02
    ph = [capacitance, i, epolarization, gdepolarization, edepolarization, currentpump]
    p = [gmin,gmax,ph,pt]
    pl = MultiPlot(p)
savefig(pl, "Cervera_Plot_$i.png")
end
=#

p = [1.0,2.0]
pl = MultiPlot(p)
savefig(pl, "Cervera_Plot.png")

#=
for i = 0.25:0.25:2.0
    pl = PlotGap!(0.2, i*0.1)
    savefig(pl, "GapJunction_$i.png")
end
=#