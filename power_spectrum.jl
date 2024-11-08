const omegaM0_h2 = 0.1326 
const omegaB0_h2 = 0.0227
const h = 0.719
const ns = 0.963
const sigma8 = 0.798
const c = 300000000 #m/s
const pnl = 0.5
const sig_rz = 0 #from sig_z = 0 assumption
const sig_p = 1 #in km/s

include("EH_transferfunction_v2.jl")

using QuadGK
using Plots
using SpecialFunctions
using .EH_transferfunction

params = EHParamsConstructor(omegaM0_h2/h^2, omegaB0_h2/h^2, h; theta_cmb=2.726/2.7)

function HubbleToHertz(H)
    inHertz = H*3.241*10^-20
    return inHertz
end

function mToMpc(d)
    inMpc = d*3.241*10^-23
    return inMpc
end

function HubbleConst(h = 0.719)
    return h*100.0
end

function HubbleParameter(z, h, omegaM0_h2)
    omegaM0 = omegaM0_h2/h^2
    return HubbleConst(h)*(omegaM0 * (1+z)^3 + (1-omegaM0))^0.5
end

function AngularDiameterDistance(z, h, omegaM0_h2)
    arg(z) = 1/HubbleToHertz(HubbleParameter(z, h, omegaM0_h2))
    return (c/(1+z)) * quadgk(arg, 0, z)[1]
end

function omegaM(z, h, omegaM0_h2)
    omegaM0 = omegaM0_h2/h^2
    omegaL0 = (1-omegaM0)
    return omegaM0 * (1 + z)^3 / ((omegaM0 * (1 + z)^3 + omegaL0))
end

function omegaL(z, h, omegaM0_h2)
    omegaM0 = omegaM0_h2/h^2
    omegaL0 = (1-omegaM0)
    return omegaL0 / ((omegaM0 * (1 + z)^3 + omegaL0))
end