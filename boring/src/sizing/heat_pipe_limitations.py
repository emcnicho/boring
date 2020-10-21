"""
This code is for estimating the performance limitations of a heat pipe.

Equations are from:
    "Heat Pipe Science and Technology", 2nd Edition, 2016
    Amir Faghri

Results are validated against:
    "Mathematical Model for Heat Transfer Limitations of Heat Pipe"
    Patrik Nemec, et al.
    Mathematical and Computer Modelling, 2013
    
Author: Ezra McNichols
        NASA Glenn Research Center
        Turbomachinery and Turboelectric Systems
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', family='serif')
ax=plt.figure(1,figsize=(9,6)).add_subplot(1,1,1)
#################################################################################

# Validation Parameters (Nemec et al, 2013, "Mathemtatical model for heat transfer limitations of heat pipe")
r_i=0.0065
r_hv=0.005
L_e=0.15
L_a=0.2
L_c=0.15
L_eff=0.5*L_e+L_a+0.5*L_c
A_w=0.000054
k_s=393
r_p=0.0001/2
epsilon=0.65
t=0.0015
L_t=0.5
T_w=298
g=9.81
phi=180
r_ce=r_p
r_n=25e-6
K=(r_p*2)**2*epsilon**3/(150*(1-epsilon)**2)
A_v=np.pi*r_hv**2
#################################################################################
def f(T,a_0,a_1,a_2,a_3,a_4,a_5):
    poly=np.exp(a_0+a_1*T+a_2*T**2+a_3*T**3+a_4*T**4+a_5*T**5)
    return poly
def h(T,a_0,a_1,a_2,a_3,a_4,a_5):
    poly=a_0+a_1*T+a_2*T**2+a_3*T**3+a_4*T**4+a_5*T**5
    return poly
T_hpfp=np.arange(0.1,150,0.1)
T_hp=T_hpfp+273.15
######################################## Fluid Properties for Ethanol. From Faghri's book ########################################################
P_v = f(T_hpfp,-4.4114,8.7650e-2,-6.3182e-4,3.9958e-6,-1.4340e-8,2.0359e-11)*1e5
h_fg = h(T_hpfp,1048.6,-1.0921,1.0651e-2,-2.0693e-4,1.1231e-6,-2.4928e-9)*1e3
rho_l = f(T_hpfp,-1.0791e-1,-7.7201e-3,1.5906e-4,-1.6139e-6,7.1873e-9,-1.2075e-11)*1e3    
rho_v = f(T_hpfp,-3.3681,5.2492e-2,5.1630e-5,-1.9542e-6,8.6893e-9,-1.1451e-11)
mu_l = f(T_hpfp,5.8942e-1,-2.2540e-2,1.0283e-4,-8.8574e-7,4.7884e-9,-9.7493e-12)/1e3
mu_v = f(T_hpfp,-2.5759e-1,4.5249e-3,-3.1212e-5,3.9144e-7,-2.3733e-9,5.1450e-12)/1e5
k_l = f(T_hpfp,-1.6976,-1.2505e-3,7.5291e-7,5.2361e-8,-3.4986e-10,6.4599e-13)
k_v = f(T_hpfp,-4.4346,3.3797e-3,2.1001e-4,-3.4778e-6,2.0462e-8,-4.0325e-11)
sigma_l = h(T_hpfp,24.419,-8.1477e-2,-1.1450e-4,8.6540e-7,-7.6432e-9,1.9148e-11)/1e3
cp_l = f(T_hpfp,0.81763,2.6793e-3,1.3888e-5,-4.3856e-11,-4.4424e-10,1.5104e-12)*1e3
cp_v = f(T_hpfp,2.9255e-1,1.2271e-3,8.0938e-5,-1.8513e-6,1.6850e-8,-5.3880e-11)*1e3
v_fg = 1/rho_v-1/rho_l
R_g=P_v/(T_hp*rho_v)
cv_v=cp_v-R_g
gamma=cp_v/cv_v
# print("T= ",T_hpfp)
# print(P_v*1e-5,h_fg,rho_l*1e-3,rho_v,mu_l*1e3,mu_v*1e5,k_l,k_v,sigma_l*1e3,cp_l/1e3,cp_v/1e3)   
##################################################################################

################################ Calculations ################################ 
k_eff=k_s*(2+k_l/k_s-2*epsilon*(1-k_l/k_s))/(2+k_l/k_s+epsilon*(1-k_l/k_s))

q_boiling=4*np.pi*L_eff*k_eff*T_hp*sigma_l/(h_fg*rho_v*np.log(r_i/r_hv))*(1/r_n-1/r_ce)
q_sonic=A_v*rho_v*h_fg*np.sqrt(gamma*R_g*T_hp)*np.sqrt(1+gamma)/(2+gamma)
q_ent=A_v*h_fg*np.sqrt(sigma_l*rho_v/(2*r_p))
q_vis=(r_hv*2)**2*h_fg*rho_v*P_v/(64*mu_v*L_eff)*A_v
q_cap=sigma_l*rho_l*h_fg*K*A_w/(mu_l*L_eff)*(2/r_ce-rho_l*g*L_t*np.cos(phi*np.pi/180)/sigma_l)
##############################################################################
plt.semilogy(T_hpfp,q_cap,color='red')
plt.semilogy(T_hpfp,q_boiling,color='orange')
plt.semilogy(T_hpfp,q_sonic,color='pink')
plt.semilogy(T_hpfp,q_ent,color='green')
plt.semilogy(T_hpfp,q_vis,color='blue')
plt.grid(True,which="both")
plt.xlabel('Temperature [C]')
plt.ylabel('Heat [W]')

ax.legend(['Capillary Limit','Boiling Limit','Sonic Limit','Entrainment Limit','Viscous Inertia Limit'])
plt.show()