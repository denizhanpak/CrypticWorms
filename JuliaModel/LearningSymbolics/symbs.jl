using Symbolics
using SymbolicUtils
using Plots

@variables vh vt vd
#vd = vh - vt
gref = 0.1
capacitance = 100
gpolt = 4.0 * gref
gpolh = 1.5 * gref
epolarization = -70.0
gdepolarization = 0.5 * gref
edepolarization = 0.0
currentpump = 0.0
gmin = 0.2 * gref
gmax = 2.0 * gref #2.0 for cusp region
vth = 20.0
v0 = 2.0
gap = gmin + (gmax - gmin)/((1+exp((vd-vth)/v0))*(1 + exp(-(vd - vth)/v0)))
vho = (gpolh * epolarization + gdepolarization * edepolarization)/(gpolh+gdepolarization)
vto = (gpolt * epolarization + gdepolarization * edepolarization)/(gpolt+gdepolarization)
eq = (1 + gap/(gpolh+gdepolarization) + gap/(gpolt+gdepolarization)) * (vd) - (vho-vto)
f = build_function(eq, vd, expression=Val{false})
vds = -70:70
rv = [f(x) for x in vds]
p = plot(vds,rv)
plot!(p,vds,zeros(length(vds)))