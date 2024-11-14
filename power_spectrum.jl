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

function AngularDiameterDistance(z, ps_params)
    arg(z) = 1/(PSz_paramsconstructor(z, ps_params).H)
    return (ps_params.c/(1+z)) * quadgk(arg, 0, z)[1]
end

function I0(t_params, p_params)
    arg(k) = k^(p_params.ns + 2.) * (3 * besselj1(8*k)/8/k)^2 * Tm(k,t_params)^2
    return quadgk(arg, 0, Inf)[1] / 2 / pi^2
end

function P_tilde(I0,ps_params)
    return (ps_params.sigma8)^2/I0/(ps_params.h)^ps_params.ns
end

function sigma_m(P_tilde, ps0_params)
    return ps0_params.Gz * P_tilde
end

function sigma_g(sigma_m, b)
    return b * sigma_m
end

function g_μ(μ,ps0_params)
    return ps0_params.Gz * (1 - μ^2 + μ^2 * (1+ps0_params.fg)^2)
end

function Tdw(μ, k, ps0_params, ps_params, transfer_params)
    exp_term = exp(-g_μ(μ,ps0_params) * k^2 / ps_params.kstar^2 / 2)
    return exp_term * Tm(k,transfer_params)^2 + (1-exp_term) * Tc(k,transfer_params)^2
end

function generate_powerspectrum(μ, k, ps0_params, ps_params, transfer_params)
    # AP_factor_DA = (AngularDiameterDistance(0, PS_paramsconstructor(0.1326, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1))/AngularDiameterDistance(0, ps_params))^2
    # AP_factor_H = ps0_params.H / PSz_paramsconstructor(0, PS_paramsconstructor(0.1326, 0.0227, 0.719, 0.963, 0.798, 0.5, 0, 1)).H

    sigma_factor = (sigma_g(sigma_m(P_tilde(I0(transfer_params, ps_params), ps_params), ps0_params), 1) 
                   + ps0_params.fg * sigma_m(P_tilde(I0(transfer_params, ps_params), ps_params), ps0_params) * μ ^2 )^2 #here set b(z) = 1

    T_dw2 = Tdw(μ, k, ps0_params, ps_params, transfer_params)^2

    return sigma_factor * T_dw2 * (k*ps_params.h)^ps_params.ns / ( 1 + k^2 * μ^2 * ps0_params.sig_rp^2/2)
end

function plot_ps()
    ps = []
    k_arr = 10.0.^range(-4,stop=0,length=200)
    for k in k_arr#surely this is not optimal
        push!(ps, generate_powerspectrum(0, k, ps0_params, ps_params, transfer_params))
    end

    plot(k_arr, ps, xaxis=:log, yaxis=:log)
    savefig("./Uni/PhD/fisher_project/ps_test.png")
end

plot_ps()