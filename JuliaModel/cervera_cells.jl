function GapJunction!(vh, vt, 
    gmin, gmax, vth=20.0, v0=2.0)
    
    numerator = gmax - gmin
    denominator_term1 = (1 + exp((vh - vt - vth)/v0))
    denominator_term2 = (1 + exp(-(vh - vt + vth)/v0))

    term2 = numerator / (denominator_term1 * denominator_term2)

    return (gmin + term2)
end

function SingleCell!(v, p,gap)
    #Parameters for single cell
    capacitance, gpol, epol, gdep, edep, inp = p
    
    #Polarization and depolarization channels
    ipol = gpol * (v - epol)
    idep = gdep * (v - edep)
    
    #Return derivative of voltage
    return (-ipol - idep - gap - inp)/capacitance
end

function TwoCellCervera!(v,p::Vector{Float64},t::Float64)

    gref = 0.1
    capacitance = 100
    gpolarization_tail = 4.0 * gref
    gpolarization_head = p[1] * gref #0.011 for cusp region
    epolarization = -70.0
    gdepolarization = 0.5 * gref
    edepolarization = 0.0
    currentpump = 0.0
    gmin = 0.2 * gref
    gmax = p[2] * gref #2.0 for cusp region

    #Coefficients of the functions. 
    pt = [capacitance, gpolarization_tail, epolarization, gdepolarization, edepolarization, currentpump]
    ph = [capacitance, gpolarization_head, epolarization, gdepolarization, edepolarization, currentpump]

    vh, vt = v   #state vector of voltages     
    #Calculate Gap junction conductance
    #dv[1] = gap = GapJunction!(vh, vt, gmin, gmax,0.0,20.0)
    gap = 0.0
    gap = GapJunction!(vh, vt, gmin, gmax, 20.0, 2.0)
    
    #Calculate single cell voltages
    diff = vt - vh
    dvh = SingleCell!(vh, ph, gap * -diff)
    dvt = SingleCell!(vt, pt, gap * diff)
    
    return SVector(dvh, dvt)
end