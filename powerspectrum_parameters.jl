module ps_parameters
export PS_paramsconstructor, PSz_paramsconstructor

struct PSParams
    # Cosmological parameters
    omm0::Float64
    omb0::Float64
    h::Float64

    H0::Float64
    ns::Float64
    sigma8::Float64
    c::Float64
    pnl::Float64
    sig_rz::Float64
    sig_p::Float64

    kstar::Float64
end

struct PSzParams
    #Density parameters
    omm::Float64
    omL::Float64
    omb::Float64
    omc::Float64

    #Growth parameters
    Gz::Float64
    fg::Float64

    #Hubble
    H::Float64

    #sigmas
    sig_rp::Float64
end

function PSz_paramsconstructor(z, params)

    omm0 = params.omm0
    omL0 = 1-params.omm0
    omm = omm0 * (1+z)^3 / ((omm0 * (1 + z)^3 + omL0))
    omL = omL0  / ((omm0 * (1 + z)^3 + omL0))

    omb0 = params.omb0
    omc0 = omm0 - omb0

    omb = omb0 * (1+z)^3 / ((omm0 * (1 + z)^3 + omL0))
    omc = omc0 * (1+z)^3 / ((omm0 * (1 + z)^3 + omL0))

    Gz = 5 * omm/(2 * (1+z)) * (omm^(4/7) - omL + (1 + omm/2) * (1 + omL/70))^-1.
    fg = omm^0.55

    H = params.H0 * 3.241*10^-20 * (omm0 * (1+z)^3 + (1-omm0))^0.5 #in inverse seconds

    sig_rp = params.sig_p * (1 + z) / H

    PSzParams(omm,omL,omb,omc,Gz,fg,H,sig_rp)

end

function PS_paramsconstructor(ommh2, ombh2, h, ns, sigma8, pnl, sig_rz, sig_p)
    omm0 = ommh2/h^2
    omb0 = ombh2/h^2
    H0 = h*100.
    # c = 300000000
    kstar = (8.355 * pnl * sigma8 / h / 0.8)^(-1) #in Mpc

    PSParams(omm0, omb0, h, H0, ns, sigma8, 300000000, pnl, sig_rz, sig_p, kstar)
end

end