# Transfer function of Eisenstein & Hu 1998
# (Equation numbers refer to this paper)
import numpy as np
import time
import matplotlib.pyplot as plt

ts = time.time()
"""
k
"""
k = np.linspace(1e-4,1,200)

"""
CONSTANTS
"""

om_c = 0.2685
om_b = 0.049
om_m = om_c + om_b
h = 0.6711
hsq = h*h
om_mhsq = om_m*hsq

e = np.e
thet = 2.728/2.7
thetsq = thet*thet
thetpf = thetsq*thetsq

"""
Equation 4 - redshift of drag epoch
"""
b1 = 0.313*om_mhsq**(-0.419)*(1+0.607*om_mhsq**(0.674))
b2 = 0.238*om_mhsq**(0.223)
zd = 1291*(1+b1*(om_b*hsq)**b2)*(om_mhsq**(0.251))/(1+0.659*om_mhsq**(0.828))

"""
Equation 2 - redshift of matter-radiation equality
"""
ze = 2.5E4*om_mhsq/thetpf

"""
R (ratio of baryon-photon momentum density) at drag epoch
"""
rd = 31500*om_b*hsq/(thetpf*zd)

"""
R (ratio of baryon-photon momentum density) at matter-radiation equality
"""
re = 31500*om_b*hsq/(thetpf*ze)

"""
Equation 3 - scale of particle horizon at matter-radiation equality
"""
rke = 7.46E-2*om_mhsq/thetsq

"""
Equation 6 - scale of sound horizon at drag epoch
"""
s = (2/(3*rke))*np.sqrt(6/re)*np.log((np.sqrt(1+rd)+np.sqrt(rd+re))/(1+np.sqrt(re)))

"""
Equation 7 - Silk damping scale
"""
rks = 1.6*((om_b*hsq)**0.52)*((om_mhsq)**0.73)*(1+(10.4*om_mhsq)**-0.95)

"""
Equations 11 - CDM transfer function fits
"""
a1 = ((46.9*om_mhsq)**0.670)*(1+(32.1*om_mhsq)**(-0.532))
a2 = ((12.0*om_mhsq)**0.424)*(1+(45.0*om_mhsq)**(-0.582))
ac = (a1**(-om_b/om_m))*(a2**(-(om_b/om_m)**3))

"""
Equations 12 - CDM transfer function fits
"""
bb1 = 0.944/(1+(458*om_mhsq)**-0.708)
bb2 = (0.395*om_mhsq)**-0.0266
bc = 1/(1+bb1*((om_c/om_m)**bb2)-1)

"""
Equation 15
"""
y = (1+ze)/(1+zd)
g = y*(-6*np.sqrt(1+y)+(2+3*y)*np.log((np.sqrt(1+y)+1)/(np.sqrt(1+y)-1)))

"""
Equation 14
"""
ab = g*2.07*rke*s/((1+rd)**0.75)

"""
Equation 23
"""
bn = 8.41*om_mhsq**0.435

"""
Equation 24
"""
bb = 0.5 + (om_b/om_m) + (3-2*om_b/om_m)*np.sqrt((17.2*om_mhsq)**2 + 1)

"""
Convert k to Mpc^-1
"""
rk = k*h

"""
Equation 10 - define q 
"""
q = rk/(13.41*rke)

"""
Equation 20
"""
c1=14.2 + 386/(1+69.9*q**1.08)
c2=14.2/ac + 386/(1+69.9*q**1.08)

"""
Equation 18
"""
f = 1/(1+(rk*s/5.4)**4)

"""
Equation 17 - CDM transfer function
"""
tc = f*np.log(e+1.8*bc*q)/(np.log(e+1.8*bc*q)+c1*q*q) + (1-f)*np.log(e+1.8*bc*q)/(np.log(e+1.8*bc*q)+c2*q*q)

"""
Equation 22
"""
ss = s/(1+(bn/(rk*s))**3)**(1/3)

"""
Equation 19 & 21
"""
tb = np.log(e+1.8*q)/(np.log(e+1.8*q)+c1*q*q)/(1+(rk*s/5.2)**2)
tb = (tb+ab*np.exp(-(rk/rks)**1.4)/(1+(bb/rk/s)**3))*np.sin(rk*ss)/rk/ss

"""
Equation 8
"""
tk_eh = (om_b/om_m)*tb+(1-om_b/om_m)*tc

plt.plot(rk, tk_eh)

# # sigg = 

# tk_nb = (1-om_b/om_m)*tc
# tk_wb = tk_nb + (om_b/om_m)*tb

# # g_e06 = np.exp(-0.5*k*k*sigg*sigg)
# g_e06 = 1

# bao = g_e06*tk_wb*tk_wb/(tk_nb*tk_nb)+(1-g_e06)

# kh = k/h

# plt.plot(kh, bao) 
# # plt.yscale('log')
# # plt.xscale('log')
# plt.axhline(1)

# powspec = np.loadtxt("Interlopers/codes/fiducial/CAMB_matterpow_0.dat", unpack=True)
# kk=powspec[0]
# pk=powspec[1]

# def ps_nw(lk, lamb, Plin, qlin):
#     # print(len(log_k))
#     log_q = np.log10(qlin)
#     dlog_q = []
#     for i in range(len(qlin)):
#         if i == 883:
#             break
#         dlog_q.append(log_q[i+1] - log_q[i])
#     dlog_q.append(dlog_q[-1]+0.000003)
#     dlog_q = np.array(dlog_q)
#     p_nw = (1/(np.sqrt(2*np.pi)*lamb))*np.sum(dlog_q*Plin*np.exp((-1/(2*lamb**2))*(lk-log_q)**2))
#     return p_nw


# log10_k = np.log10(kk)


# nw_list = [] 
# for lk in log10_k:
#     nw = ps_nw(lk, 0.25, pk, kk)
#     nw_list.append(nw)

# nw_list = np.array(nw_list)
# nw_interp = np.interp(kh, kk, nw_list)

# plin_interp = np.interp(kh, kk, pk)

# plt.plot(np.log10(kh), np.log10(plin_interp))



# cdm_ps = (tc*tc*kh)
# dif = np.max(np.log10(plin_interp)) - np.max(np.log10(cdm_ps))
# plt.plot(np.log10(kh), plin_interp/(cdm_ps*10**dif))


# print(np.max(np.log10(nw_list)) - np.max(np.log10(pk)))
# plt.plot(log10_k, np.log10(pk))
# plt.plot(log10_k, np.log10(nw_list))
# plt.plot(np.log10(kk), np.log10(pk)+2)
# powspec_interp = np.interp(kh, kk, pk)
# bao_interp = np.interp(kk, kh, bao)

# plt.plot(kh,powspec_interp)
# plt.plot(kk,bao_interp, '-')
# plt.plot(kk, pk/bao_interp)
# plt.plot(kk, pk)
# plt.xlim(1e-4,1)
# plt.ylim(100,10000)
# plt.yscale('log')
plt.xscale('log')


plt.savefig('eh_bao.png')
print(time.time() - ts)