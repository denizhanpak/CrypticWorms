using Plots 

function PlotGap!(gmax::Float64=.20,gmin::Float64=0.02)
    gjs = []
    for i in 60.0:-1.0:-60.0
        gj = GapJunction!(0.0,i,gmin,gmax)
        append!(gjs, gj)
    end
    rv = plot(-60.0:60.0, gjs)
    return rv
end

function MultiPlot(p)
    # main plot settings
    rv = plot(linewidth=4, title="Cellular Voltages",
        xaxis="t",yaxis="Voltage(t)",legend=true)
    step = 10
    total_steps = 10000.0
    plot!([0; total_steps], [50; 50], lw=2, lc=:black, label="",line=(:dot, 4))
    legend_set = false
    for j = -70.0:step:70.0
        for i = -70.0:step:70.0
            u0 = [i;j]
            tspan = (0.00,total_steps)
            prob = ODEProblem(TwoCellCervera!, u0, tspan, p)
            sol = solve(prob)

            # add to existing plot
            if !legend_set
                plot!(sol,vars=(1),color="yellow",linealpha=1,label="Head Cell")
                plot!(sol,vars=(2),color="magenta",linealpha=1,label="Tail Cell")
                #plot!(sol,vars=(3),color="cyan",linealpha=1,label="Voltage Difference")
                legend_set = true
            else
                plot!(sol,vars=(1),color="yellow",linealpha=0.1,label="")
                plot!(sol,vars=(2),color="magenta",linealpha=0.1,label="")
                #plot!(sol,vars=(3),color="cyan",linealpha=0.1,label="")
            end
        end
    end
    return rv
end

function basin_plot(fractions, gmaxs, gpols)
    gr()
    counter = 1
    rv = zeros(Float64, length(gmaxs),length(gpols))
    x = string.(gmaxs)
    y = string.(gpols)
    for (i,gmax) in enumerate(gmaxs)
        for (j,gpol) in enumerate(gpols)
            rv[i,j] = fractions[counter][1]
            counter += 1
        end
    end

    h = heatmap(gmaxs,
    gpols, rv,
    xlabel="Gap Junction Max Capacitance", 
    ylabel="Head Cell Polarization Capacitance",
    title="Probability Map")

    
    return plot(h)
end