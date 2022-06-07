using DynamicalSystems
include("./cervera_cells.jl")

p = [0.011, 2.0]
initial = SVector(0.0,0.0)

ds = ContinuousDynamicalSystem(TwoCellCervera!, initial, p)
xg = yg = range(-70, 70; length = 70)
grid = (xg, yg)
am = AttractorsViaRecurrences(ds, grid)