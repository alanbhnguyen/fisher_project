include("EH_transferfunction_v2.jl")
include("powerspectrum_parameters.jl")

using QuadGK
using Plots
using SpecialFunctions
using .EH_transferfunction
using .ps_parameters

transfer_params = EHParamsConstructor(0.1326/0.719^2, 0.0227/0.719^2, 0.719; theta_cmb=2.726/2.7)
ps_params = PS_paramsconstructor(0.1326, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1)
ps0_params = PSz_paramsconstructor(0, ps_params)

#not a function of z
function I0(t_params, p_params)
    arg(k) = k^(p_params.ns + 2.) * (3 * besselj1(8*k)/8/k)^2 * Tm(k,t_params)^2
    return quadgk(arg, 0, Inf)[1] / 2 / pi^2
end

function P_tilde(I0,ps_params)
    return (ps_params.sigma8)^2/I0/(ps_params.h)^ps_params.ns
end

####
####
####

function AngularDiameterDistance(z, ps_params)
    arg(z) = 1/(PSz_paramsconstructor(z, ps_params).H)
    return (ps_params.c/(1+z)) * quadgk(arg, 0, z)[1]
end

function sigma_m(P_tilde, z, ps_params)
    psz_params = PSz_paramsconstructor(z, ps_params)
    return psz_params.Gz * P_tilde ^ 0.5
end

function sigma_g(sigma_m, b)
    return b * sigma_m
end

function g_μ(μ, z, ps_params)
    psz_params = PSz_paramsconstructor(z, ps_params)
    return psz_params.Gz * (1 - μ^2 + μ^2 * (1+psz_params.fg)^2)
end

####
####
####

function Tdw(μ, k, z, ps_params, transfer_params)
    # psz_params = PSz_paramsconstructor(z, ps_params)
    exp_term = exp(-g_μ(μ,z,ps_params) * k^2 / ps_params.kstar^2 / 2)
    return exp_term * Tm(k,transfer_params)^2 + (1-exp_term) * Tc(k,transfer_params)^2
end

function generate_powerspectrum(μ, k, z, ps_params, transfer_params)

    psz_params = PSz_paramsconstructor(z, ps_params)
    if z == 0
        AP_factor_DA = 1
        AP_factor_H = 1
    else
        AP_factor_DA = (AngularDiameterDistance(z, PS_paramsconstructor(0.1326, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1))/AngularDiameterDistance(z, ps_params))^2
        AP_factor_H = psz_params.H / PSz_paramsconstructor(z, PS_paramsconstructor(0.1326, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1)).H
    end
    sigma_factor = (sigma_g(sigma_m(P_tilde(I0(transfer_params, ps_params), ps_params), z, ps_params), 1)  #what is b here
                   + ps0_params.fg * sigma_m(P_tilde(I0(transfer_params, ps_params), ps_params), z, ps_params) * μ ^2 )^2 #here set b(z) = 1

    T_dw2 = Tdw(μ, k, z, ps_params, transfer_params)^2

    return AP_factor_DA * AP_factor_H * sigma_factor * T_dw2 * (k*ps_params.h)^ps_params.ns / ( 1 + k^2 * μ^2 * ps0_params.sig_rp^2/2)
end

function sigma_m_2(ns, z, ps_params)
    arg(k) = k^(ns + 2.) * (3 * besselj1(8*k)/8/k)^2 * Tm(k,t_params)^2
    I0 = quadgk(arg, 0, Inf)[1] / 2 / pi^2
    P_tilde = (ps_params.sigma8)^2/I0/(ps_params.h)^ns
    psz_params = PSz_paramsconstructor(z, ps_params)
    return psz_params.Gz * P_tilde ^ 0.5
end

function sigma_m_3(z, ps_params)
    arg(k) = k^(ps_params.ns + 2.) * (3 * besselj1(8*k)/8/k)^2 * Tm(k,t_params)^2
    I0 = quadgk(arg, 0, Inf)[1] / 2 / pi^2
    P_tilde = (ps_params.sigma8)^2/I0/(ps_params.h)^ps_params.ns
    psz_params = PSz_paramsconstructor(z, ps_params)
    return psz_params.Gz * P_tilde ^ 0.5
end

function Tdw_2(μ, k, z, ps_params, transfer_params)
    # psz_params = PSz_paramsconstructor(z, ps_params)
    exp_term = exp(-g_μ(μ,z,PS_paramsconstructor(0.1326, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1)) * k^2 / ps_params.kstar^2 / 2)
    return exp_term * Tm(k,transfer_params)^2 + (1-exp_term) * Tc(k,transfer_params)^2
end

function ns_deriv_term(μ, k, z, ns, ps_params)
    psz_params = PSz_paramsconstructor(z, ps_params)
    return 2 * log(sigma_m_2(ns, 0, ps_params)+psz_params.fg*sigma_m_2(ns, 0, ps_params)*μ^2) + ns*log(k*ps_params.h)
end

function ω_m_deriv_term(μ, k, z, ps_params, transfer_params)
    psz_params = PSz_paramsconstructor(z, ps_params)
    return (- 2*log(AngularDiameterDistance(z, ps_params)) + log(psz_params.H)
            + 2*log(sigma_m_3(z,ps_params) + psz_params.fg*sigma_m_3(z,ps_params)*μ^2)
            + 2*log(Tdw_2(μ, k, z, ps_params, transfer_params)))
end

function generate_logpowerspectrum_derivs(μ, k, z, ps_params, index)
    """
    indices
    0: ln H
    1: ln DA
    2: ln (f_g σ_m)
    3: ω_m
    4: ω_b
    5: n_s
    """
    psz_params = PSz_paramsconstructor(z, ps_params)

    if index == 0
        deriv = 1
    elseif index == 1
        deriv = -2
    elseif index == 2
        deriv = 2 * psz_params.fg * sigma_m_2(ps_params.ns, 0, ps_params) * μ^2/(sigma_m_2(ps_params.ns, 0, ps_params) + psz_params.fg * sigma_m_2(ps_params.ns, 0, ps_params) * μ^2)
    elseif index == 3
        deriv = (-ω_m_deriv_term(μ,k,z,PS_paramsconstructor(0.1326*1.02, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1),EHParamsConstructor(0.1326*1.02/0.719^2, 0.0227/0.719^2, 0.719; theta_cmb=2.726/2.7))
                 +8*ω_m_deriv_term(μ,k,z,PS_paramsconstructor(0.1326*1.01, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1),EHParamsConstructor(0.1326*1.01/0.719^2, 0.0227/0.719^2, 0.719; theta_cmb=2.726/2.7))
                 -8*ω_m_deriv_term(μ,k,z,PS_paramsconstructor(0.1326*0.99, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1),EHParamsConstructor(0.1326*0.99/0.719^2, 0.0227/0.719^2, 0.719; theta_cmb=2.726/2.7))
                 +ω_m_deriv_term(μ,k,z,PS_paramsconstructor(0.1326*0.98, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1),EHParamsConstructor(0.1326*0.98/0.719^2, 0.0227/0.719^2, 0.719; theta_cmb=2.726/2.7)))/(0.12*ps_params.omm0)
    elseif index == 4
        deriv = (-ω_m_deriv_term(μ,k,z,PS_paramsconstructor(0.1326, 0.0227*1.02, 0.719, 0.963, 0.798, 0.5, 0, 1),EHParamsConstructor(0.1326/0.719^2, 0.0227*1.02/0.719^2, 0.719; theta_cmb=2.726/2.7))
                 +8*ω_m_deriv_term(μ,k,z,PS_paramsconstructor(0.1326, 0.0227*1.01, 0.719, 0.963, 0.798, 0.5, 0, 1),EHParamsConstructor(0.1326/0.719^2, 0.0227*1.01/0.719^2, 0.719; theta_cmb=2.726/2.7))
                 -8*ω_m_deriv_term(μ,k,z,PS_paramsconstructor(0.1326, 0.0227*0.99, 0.719, 0.963, 0.798, 0.5, 0, 1),EHParamsConstructor(0.1326/0.719^2, 0.0227*0.99/0.719^2, 0.719; theta_cmb=2.726/2.7))
                 +ω_m_deriv_term(μ,k,z,PS_paramsconstructor(0.1326, 0.0227*0.98, 0.719, 0.963, 0.798, 0.5, 0, 1),EHParamsConstructor(0.1326/0.719^2, 0.0227*0.98/0.719^2, 0.719; theta_cmb=2.726/2.7)))/(0.12*ps_params.omb0)
    elseif index == 5
        deriv = (-ns_deriv_term(μ, k, z, 1.02*ps_params.ns, ps_params) + 8*ns_deriv_term(μ, k, z, 1.01*ps_params.ns, ps_params) 
                -8*ns_deriv_term(μ, k, z, 0.99*ps_params.ns, ps_params) + ns_deriv_term(μ, k, z, 0.98*ps_params.ns, ps_params)) / (0.12*ps_params.ns)
    end

    return deriv
end

# transfer_params = EHParamsConstructor(0.1326/0.719^2, 0.0227/0.719^2, 0.719; theta_cmb=2.726/2.7)
# ps_params = PS_paramsconstructor(0.1326, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1)
# function I0(t_params, p_params)
#     arg(k) = k^(p_params.ns + 2.) * (3 * besselj1(8*k)/8/k)^2 * Tm(k,t_params)^2
#     return quadgk(arg, 0, Inf)[1] / 2 / pi^2
# end

# function P_tilde(I0,ps_params)
#     return (ps_params.sigma8)^2/I0/(ps_params.h)^ps_params.ns
# end

function plot_ps()
    ps = []
    k_arr = 10.0.^range(-4,stop=0,length=200)
    for k in k_arr#surely this is not optimal
        push!(ps, generate_powerspectrum(0, k, 0, ps_params, transfer_params))
    end

    plot(k_arr, ps, xaxis=:log, yaxis=:log, label="power_spectrum")
    savefig("./Uni/PhD/fisher_project/ps_test.png")
end

plot_ps()