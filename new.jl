include("EH_transferfunction_v2.jl")

using .EH_transferfunction
using Plots

obh2 = 0.1326
omh2 = 0.0227
h = 0.719

function plotEH()
    params_struct = EHParamsConstructor(omh2/h^2, obh2/h^2, h; theta_cmb=2.726/2.7)

    T_c = []
    T_b = []

    for k in range(1e-4,stop=1,length=500)#surely this is not optimal
        push!(T_c, Tc(k, params_struct))
        push!(T_b, Tb(k, params_struct))
    end

    plot(range(1e-4,stop=1,length=500), (T_c + T_b), xaxis=:log)
    # plot(range(1e-4,stop=1,length=500), T_b, xaxis=:log)

    savefig("EH3.png")
end

plotEH()