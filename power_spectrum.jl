const omegaM0_h2 = 0.1326 
const omegaB0_h2 = 0.0227
const h = 0.719
const ns = 0.963
const sigma8 = 0.798
const c = 300000000 #m/s
const pnl = 0.5
const sig_rz = 0 #from sig_z = 0 assumption
const sig_p = 1 #in km/s

include("EH_transferfunction.jl")

using QuadGK
using Plots
using SpecialFunctions
using .EH1998

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

function fg(z, h, omegaM0_h2)
    return omegaM(z, h, omegaM0_h2)^0.55
end

function G_z(z, h, omegaM0_h2)
    return (5 * omegaM(z, h, omegaM0_h2)/(2 * (1+z))) * (omegaM(z, h, omegaM0_h2)^(4/7) - omegaL(z, h, omegaM0_h2) + (1 + omegaM(z, h, omegaM0_h2)/2) * (1 + omegaL(z, h, omegaM0_h2)/70))^-1.
end

function sig_rp(z, sig_p, h, omegaM0_h2)
    return sig_rp = (1+z) * sig_p / HubbleParameter(z, h, omegaM0_h2)
end

function I_0(ns, h, omegaM0_h2, omegaB0_h2) #check the kbar
    arg(k) = k ^ (ns + 2) / (2 * pi^2) * (3 * besselj1(8 * k) / (8 * k))^2 * EH_transferfunction(k, h, omegaM0_h2, omegaB0_h2)[3]^2
    return quadgk(arg, 0, Inf)[1]
end

function P_tilde(ns, h, omegaM0_h2, omegaB0_h2, sigma8) 
    return sigma8^2/(I_0(ns, h, omegaM0_h2, omegaB0_h2) * h^ns)
end

function sig_m(z, ns, h, omegaM0_h2, omegaB0_h2, sigma8) 
    return G_z(z, h, omegaM0_h2) * P_tilde(ns, h, omegaM0_h2, omegaB0_h2, sigma8) 
end

function sig_g(z, ns, h, omegaM0_h2, omegaB0_h2, sigma8) #add back in b(z)
    return sig_m(z, ns, h, omegaM0_h2, omegaB0_h2, sigma8)
end

function g_mu(mu, z, h, omegaM0_h2)
    return G_z(z, h, omegaM0_h2) .^2 * ( 1 - mu .^2 + mu .^2 * ( 1 + fg(z, h, omegaM0_h2) ) .^2 )
end

function EH_dewiggled(k, mu, z, h, omegaM0_h2, omegaB0_h2, sigma8, pnl)
    k_star_inv = 8.355 * (sigma8/0.8) * pnl

    T_lin = EH_transferfunction(k, h, omegaM0_h2, omegaB0_h2)[3]
    T_nw = EH_transferfunction(k, h, omegaM0_h2, omegaB0_h2)[1]
    efactor = exp( -g_mu(mu, z, h, omegaM0_h2) * k .^2 * k_star_inv .^2 / 2)
    return T_lin .^2 * efactor + T_nw .^2 * ( 1 - efactor )
end

# function observed_power_spectrum(k, mu, z, ns, h, omegaM0_h2, omegaB0_h2, sigma8, sig_p, sig_z, pnl)
#     term1 = (AngularDiameterDistance(z, 0.719, 0.1326) / AngularDiameterDistance(z, h, omegaM0_h2)) .^2 * (HubbleParameter(z, h, omegaM0_h2) / HubbleParameter(z, 0.719, 0.1326))
#     term2 = (sig_g(z, ns, h, omegaM0_h2, omegaB0_h2, sigma8) + fg(z, h, omegaM0_h2) * sig_m(z, ns, h, omegaM0_h2, omegaB0_h2, sigma8) * mu .^2 ) .^2
#     term3 = EH_dewiggled(k, mu, z, h, omegaM0_h2, omegaB0_h2, sigma8, pnl) .^2
#     term4 = 1 + k .^2 * mu .^2 * sig_rp(z, sig_p, h, omegaM0_h2) .^2/2
#     return k .^ns * term1 * term2 * term3 / term4 
# end 

# k = range(1e-4, stop=1, length=1000)

# observed_power_spectrum.(0.01, 1., 0., ns, h, omegaM0_h2, omegaB0_h2, sigma8, sig_p, 0., pnl)

println(EH_dewiggled(0.01, 0.9, 0, h, omegaM0_h2, omegaB0_h2, sigma8, pnl))


# println(observed_power_spectrum(0.01, 1, 0, ns, h, omegaM0_h2, omegaB0_h2, sigma8, sig_p, 0, pnl))
# savefig("observed_power_spectrum.png")
