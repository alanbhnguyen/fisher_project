module ps_parameters
export PS_paramsconstructor

struct PSParams
    # Cosmological parameters
    ommh2::Float64
    ombh2::Float64
    h::Float64

    H0::Float64
    ns::Float64
    sigma8::Float64
    c::Float64
    pnl::Float64
    sig_rz::Float64
    sig_p::Float64
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

    #sigmas
    sig_rp::Float64
end

function PS_z_paramsconstructor(z, params)

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

    # sig_rp = sig_p/

    PSzParams(omm,omL,omb,omc,Gz,fg)

end

function PS_paramsconstructor(ommh2, ombh2, h, ns, sigma8, pnl, sig_rz, sig_p)
    omm0 = ommh2/h^2
    omb0 = ombh2/h^2
    H0 = h*100.
    c = 300000000

    PSParams(omm0, omb0, h, H0, ns, sigma8, c, pnl, sig_rz, sig_p)
end

end