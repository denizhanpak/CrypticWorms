using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Plots")
using DifferentialEquations
using Plots

function GapJunction!(vh::Float64, vt::Float64, 
    gmin::Float64, gmax::Float64, vth::Float64=0.0, v0::Float64=20)
    
    numerator = gmax - gmin
    denominator_term1 = (1 + exp((vh - vt - vth)/v0))
    denominator_term2 = (1 + exp((vh - vt + vth)/v0))

    term2 = numerator / (denominator_term1 * denominator_term2)

    return gmin + term2
end

function SingleCell!(v, p)
    #Parameters for single cell
    cond, gpol, epol, gdep, edep, inp, gap = p
    
    #Polarization and depolarization channels
    ipol = gpol * (v - epol)
    idep = gdep * (v - edep)
    
    #Return derivative of voltage
    return (-ipol - idep - gap - inp)/cond
end

function TwoCellCervera!(dv,v,p,t)
    
    vh, vt = v   #state vector of voltages     
    gmin, gmax, ph, pt = p   #coefficients are part of vector array p
    
    #Calculate Gap junction conductance
    dv[1] = gap = GapJunction!(vh, vt, gmin, gmax)
    
    #Append gap current to parameters
    push!(ph, gap)
    push!(pt, gap)
    
    #Calculate single cell voltages
    dv[2] = dvh = SingleCell!(vh, ph)
    dv[3] = dvt = SingleCell!(vt, pt)
    
end

#Initial conditions. x=1.0, y=0.0, z=0.0
u0 = [1.0, 0.0, 0.0]

#Timespan of the simulation. 100 in this case. 
tspan = (0.0, 100.0)

#Coefficients of the functions. 
p = [10.0, 28.0, 8/3]

#Feeding the inputs to the solver 
prob = ODEProblem(parameterized_lorenz!, u0, tspan, p)
sol = solve(prob)