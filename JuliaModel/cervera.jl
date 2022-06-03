using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Plots")
using DifferentialEquations
using Plots

function GapJunction!(vh::Float64, vt::Float64, 
    gmin::Float64, gmax::Float64, vth::Float64=20.0, v0::Float64=2.0)
    
    numerator = gmax - gmin
    denominator_term1 = (1 + exp((vh - vt - vth)/v0))
    denominator_term2 = (1 + exp(-(vh - vt + vth)/v0))

    term2 = numerator / (denominator_term1 * denominator_term2)

    return (gmin + term2) * (vh - vt)
end

function SingleCell!(v, p)
    #Parameters for single cell
    capacitance, gpol, epol, gdep, edep, inp, gap = p
    
    #Polarization and depolarization channels
    ipol = gpol * (v - epol)
    idep = gdep * (v - edep)
    
    #Return derivative of voltage
    return (-ipol - idep - gap - inp)/capacitance
end

function TwoCellCervera!(dv,v,p,t)
    
    vh, vt = v   #state vector of voltages     
    gmin, gmax, ph, pt = p   #coefficients are part of vector array p
    #Calculate Gap junction conductance
    #dv[1] = gap = GapJunction!(vh, vt, gmin, gmax,0.0,20.0)
    gap = 0.0
    #gap = GapJunction!(vh, vt, gmin, gmax,20.0,2.0)
    
    #Append gap current to parameters
    push!(ph, gap)
    push!(pt, gap)
    
    #Calculate single cell voltages
    #dv[1] = 0.0``
    #dv[2] = dvh = SingleCell!(vh, ph)
    #dv[3] = dvt = SingleCell!(vt, pt)
    dv[1] = dvh = SingleCell!(vh, ph)
    dv[2] = dvt = SingleCell!(vt, pt)
    dv[3] = dvh - dvt
end

function MultiPlot(p)
    # main plot settings
    rv = plot(linewidth=4, title="Solutions",
        xaxis="t",yaxis="Voltage(t)",legend=false)
    plot!([0; 7000], [30; 30], lw=2, lc=:black, legend=false,line=(:dot, 4))
    for i = -100.00:5:70.00
        for j = -100.0:5:70.0
            u0 = [i;j;i-j]
            tspan = (0.00,7000.0)
            prob = ODEProblem(TwoCellCervera!, u0, tspan, p)
            sol = solve(prob)
            #print(sol[:,0])

            # add to existing plot
            #plot!(sol,vars=(1),color="magenta",linealpha=0.2)
            #plot!(sol,vars=(2),color="cyan",linealpha=0.2)
            plot!(sol,vars=(3),color="yellow",linealpha=0.2)
        end
    end
    return rv
end

#Initial conditions. x=1.0, y=0.0, z=0.0
u0 = [0.0;0.0;0.0]

#Timespan of the simulation. 100 in this case. 
tspan = (0.0, 500.0)

gref = 0.1

capacitance = 100
gpolarization = 4.0 * gref
epolarization = -70.0
gdepolarization = 0.5 * gref
edepolarization = 0.0
currentpump = 0.0
gmin = 0.2 * gref
gmax = 2.0 * gref

#Coefficients of the functions. 
pt = [capacitance, gpolarization, epolarization, gdepolarization, edepolarization, currentpump]

#Feeding the inputs to the solver 
#prob = ODEProblem(TwoCellCervera!, u0, tspan, p)
#sol = solve(prob)
for i = 0.04:0.001:0.05
    ph = [capacitance, i, epolarization, gdepolarization, edepolarization, currentpump]
    p = [gmin,gmax,ph,pt]
    pl = MultiPlot(p)
    savefig(pl, "Cervera_Plot_$i.png")
end