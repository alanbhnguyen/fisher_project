module EH_transferfunction
export EHParamsConstructor, Tc, Tb, Tm

# Define a record to hold all the fitting parameters 
# for the Eisenstein & Hu (1998) transfer function
struct EHParams
    # Cosmological parameters
    omh2::Float64
    ombh2::Float64
    h::Float64
    theta_cmb::Float64

    # Derived parameters
    z_d_b1::Float64
    z_d_b2::Float64
    z_eq::Float64
    k_eq::Float64
    z_d::Float64
    R_eq::Float64
    R_d::Float64
    r_s::Float64
    ksilk::Float64

    # Transfer function parameters
    αc :: Float64
    βc :: Float64
    αb :: Float64
    βb :: Float64
    βnode :: Float64
end

# Helper function to compute the fitting parameters
# Eq. 15 of Eisenstein & Hu (1998)
function Gfunc(y)
    arg = (sqrt(1+y)+1)/(sqrt(1+y)-1)
    return y * (-6 * sqrt(1 + y) + (2 + 3*y) * log(arg))
end

# Write a constructor for the EHParams type
function EHParamsConstructor(om, ob, h; theta_cmb=2.726/2.7)
    omh2 = om * h^2
    ombh2 = ob * h^2

    z_d_b2 = 0.238 * omh2 ^ 0.223
    z_d_b1 = 0.313 * omh2 ^-0.419 * (1 + 0.607 * omh2 ^ 0.674)

    z_eq = 2.50e4 * omh2 / theta_cmb^4
    k_eq = 7.46e-2 * omh2 / theta_cmb^2
    z_d = 1291 * omh2^0.251 * (1 + z_d_b1 * ombh2^z_d_b2) / (1 + 0.659 * omh2^0.828)
    R_eq = 31.5 * ombh2 * theta_cmb^-4 * (z_eq / 1000)^-1
    R_d = 31.5 * ombh2 * theta_cmb^-4 * (z_d / 1000)^-1
    r_s = (2/3/k_eq) * sqrt(6/R_eq) * log(sqrt(1+R_d) + sqrt(R_eq+R_d) / (1 + sqrt(R_eq)))
    ksilk = 1.6 * ombh2^0.52 * omh2^0.73 * (1 + (10.4 * omh2)^(-0.95))

    # Fitting parameters for the transfer function
    a1 = (46.9 * omh2)^0.670 * (1 + (32.1 * omh2)^(-0.532))
    a2 = (12.0 * omh2)^0.424 * (1 + (45.0 * omh2)^(-0.582))
    obom = ob / om
    αc = a1^(-obom) * a2^(-obom^3)
    b2 = (0.395*omh2)^(-0.0266)
    b1 = 0.944/(1 + (458 * omh2)^(-0.708))
    βc = 1 + b1*((obom)^b2 - 1)
    βc = 1/βc

    αb = 2.07 * k_eq * r_s * (1 + R_d)^(-3/4) * Gfunc((1 + z_eq)/(1+z_d))
    βb = 0.5 + obom + (3 - 2*obom) * sqrt((17.2 * omh2)^2 + 1)
    βnode = 8.41 * omh2^0.435

    EHParams(omh2, ombh2, h, theta_cmb, z_d_b1, z_d_b2, z_eq, k_eq, z_d, R_eq, R_d, r_s, ksilk, αc, βc, αb, βb, βnode)
end

function Tc(k, parameters)
    params = parameters

    #recast units of k
    q = k/13.41/params.k_eq

    #Equations 17-20 for calculating T_c
    C_1 = 14.2/1 + 386/(1+69.9*q^1.08)
    T_0_1 = log(exp(1) + 1.8 * params.βc * q)/(log(exp(1) + 1.8 * params.βc * q) + C_1 * q^2)

    C_2 = 14.2/params.αc + 386/(1+69.9*q^1.08)
    T_0_2 = log(exp(1) + 1.8 * params.βc * q)/(log(exp(1) + 1.8 * params.βc * q) + C_2 * q^2)

    f = 1/(1 + (k * params.r_s / 5.4)^4)

    return f * T_0_1 + (1 - f) * T_0_2
end

function Tb(k, parameters)
    params = parameters

    #recast units of k
    q = k/13.41/params.k_eq

    C = 14.2 + 386/(1+69.9*q^1.08)

    term1 = log(exp(1) + 1.8 * q)/(log(exp(1) + 1.8 * q) + C * q^2)/(1 + (k * params.r_s / 5.2)^2)
    term2 = params.αb/(1 + (params.βb/(k * params.r_s))^3) * exp(-(k/params.ksilk)^1.4)

    stilde = params.r_s/(1 + (params.βnode/(k * params.r_s))^3)^(1/3)
    j0 = sin(k * stilde)/k/stilde
    return (term1 + term2) * j0
end

function Tm(k, parameters) #relative to units of k in Mpc^-1
    params = parameters

    Tcdm = Tc(k, params)
    Tbar = Tb(k, params)

    om = params.omh2/(params.h)^2
    ob = params.ombh2/(params.h)^2
    return (ob/om) * Tbar + ((om-ob)/om) * Tcdm
end

end