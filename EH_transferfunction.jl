module EH1998
export EH_transferfunction

using Plots
using LinearAlgebra

# CONSTANTS
function EH_transferfunction(k, h, omegaM0_h2, omegaB0_h2)

    hsq = h^2
    om_m = omegaM0_h2 / hsq
    om_b = omegaB0_h2 / hsq
    om_c = om_m - om_b
    om_mhsq = omegaM0_h2

    e = exp(1)
    thet = 2.728 / 2.7
    thetsq = thet^2
    thetpf = thetsq^2

    # Equation 4 - redshift of drag epoch
    b1 = 0.313 * om_mhsq^(-0.419) * (1 + 0.607 * om_mhsq^(0.674))
    b2 = 0.238 * om_mhsq^(0.223)
    zd = 1291 * (1 + b1 * (om_b * hsq)^b2) * (om_mhsq^(0.251)) / (1 + 0.659 * om_mhsq^(0.828))

    # Equation 2 - redshift of matter-radiation equality
    ze = 2.5E4 * om_mhsq / thetpf

    # Ratio of baryon-photon momentum density at drag epoch
    rd = 31500 * om_b * hsq / (thetpf * zd)

    # Ratio of baryon-photon momentum density at matter-radiation equality
    re = 31500 * om_b * hsq / (thetpf * ze)

    # Equation 3 - scale of particle horizon at matter-radiation equality
    rke = 7.46E-2 * om_mhsq / thetsq

    # Equation 6 - scale of sound horizon at drag epoch
    s = (2 / (3 * rke)) * sqrt(6 / re) * log((sqrt(1 + rd) + sqrt(rd + re)) / (1 + sqrt(re)))

    # Equation 7 - Silk damping scale
    rks = 1.6 * ((om_b * hsq)^0.52) * ((om_mhsq)^0.73) * (1 + (10.4 * om_mhsq)^-0.95)

    # Equations 11 - CDM transfer function fits
    a1 = ((46.9 * om_mhsq)^0.670) * (1 + (32.1 * om_mhsq)^(-0.532))
    a2 = ((12.0 * om_mhsq)^0.424) * (1 + (45.0 * om_mhsq)^(-0.582))
    ac = (a1^(-om_b / om_m)) * (a2^(-(om_b / om_m)^3))

    # Equations 12 - CDM transfer function fits
    bb1 = 0.944 / (1 + (458 * om_mhsq)^-0.708)
    bb2 = (0.395 * om_mhsq)^-0.0266
    bc = 1 / (1 + bb1 * ((om_c / om_m)^bb2) - 1)

    # Equation 15
    y = (1 + ze) / (1 + zd)
    g = y * (-6 * sqrt(1 + y) + (2 + 3 * y) * log((sqrt(1 + y) + 1) / (sqrt(1 + y) - 1)))

    # Equation 14
    ab = g * 2.07 * rke * s / ((1 + rd)^0.75)

    # Equation 23
    bn = 8.41 * om_mhsq^0.435

    # Equation 24
    bb = 0.5 + (om_b / om_m) + (3 - 2 * om_b / om_m) * sqrt((17.2 * om_mhsq)^2 + 1)


    rk = k * h

    # Equation 10 - define q
    q = rk ./ (13.41 * rke)

    # Equation 20
    c1 = 14.2 .+ 386 ./ (1 .+ 69.9 * q.^1.08)
    c2 = 14.2 ./ ac .+ 386 ./ (1 .+ 69.9 * q.^1.08)

    # Equation 18
    f = 1 ./ (1 .+ (rk * s / 5.4).^4)

    # Equation 17 - CDM transfer function
    tc = f .* log.(e .+ 1.8 * bc .* q) ./ (log.(e .+ 1.8 * bc .* q) .+ c1 .* q.^2) .+ 
        (1 .- f) .* log.(e .+ 1.8 * bc .* q) ./ (log.(e .+ 1.8 * bc .* q) .+ c2 .* q.^2)

    # Equation 22
    ss = s ./ (1 .+ (bn ./ (rk .* s)).^3).^(1/3)

    # Equation 19 & 21
    tb = log.(e .+ 1.8 .* q) ./ (log.(e .+ 1.8 .* q) .+ c1 .* q.^2) ./ (1 .+ (rk .* s / 5.2).^2)
    tb = (tb .+ ab .* exp.(-(rk / rks).^1.4) ./ (1 .+ (bb ./ rk ./ s).^3)) .* sin.(rk .* ss) ./ rk ./ ss

    # Equation 8
    tk_eh = (om_b / om_m) .* tb .+ (1 - om_b / om_m) .* tc

    tk_wb = (om_b / om_m) .* tb #baryon only TF
    tk_nb = (1 - om_b / om_m) .* tc #CDM only TF (no wiggles)

    return [tk_nb, tk_wb, tk_eh] #return cdm, baryon, total
end

end