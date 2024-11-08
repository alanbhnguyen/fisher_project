include("EH_transferfunction_v2.jl")

using .EH_transferfunction
using Plots

omh2 = 0.1326
obh2 = 0.0227
h = 0.719

function plotEH()
    params_struct = EHParamsConstructor(omh2/h^2, obh2/h^2, h; theta_cmb=2.726/2.7)

    T_m = []
    T_c = []
    k_arr = 10.0.^range(-4,stop=0,length=200)
    for k in k_arr#surely this is not optimal
        push!(T_c, Tc(k, params_struct))
        push!(T_m, Tm(k, params_struct))
    end

    plot(k_arr, [T_c, T_m], xaxis=:log, yaxis=:log)

    savefig("EH3.png")
end

plotEH()