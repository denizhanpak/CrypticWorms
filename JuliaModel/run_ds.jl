using DynamicalSystems
using JLD2
using Random
include("./cervera_cells.jl")
include("./plotting.jl")

fs = []
gpol_range = range(0,0.5,length=2)
gmax_range = range(0.2,2.0,length=2)
for gpol in gpol_range
    for gmax in gmax_range
        p = [gpol,gmax]
        
        initial = SVector(0.0,0.0)

        ds = ContinuousDynamicalSystem(TwoCellCervera!, initial, p)
        xg = yg = range(-70, 70; length = 70)
        grid = (xg, yg)
        am = AttractorsViaRecurrences(ds, grid)

        rng = Random.MersenneTwister(1234)
        sampler, _ = statespace_sampler(rng; min_bounds=[-70., -70.], max_bounds=[10., 70.])
        push!(fs, basins_fractions(am, sampler; N = 1000, show_progress = false))
    end
end

h = basin_plot(fs, gmax_range, gpol_range)
savefig(h, "Basin Fractions.png")
save "./probs.jdl2" fs