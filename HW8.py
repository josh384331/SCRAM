import numpy as np
from matplotlib import pyplot as plt
import Rayleigh_Pitot_Solver as pitot
import Supersonic_external_flow as external
import critical_mach as cm


print("\n\nProblem 1)")
p1 = 26.436
p4 = 30.454
P04 = 143.062
gamma = 1.4
theta = np.deg2rad(10)
alpha = np.deg2rad(12)

# get M_4 from pitot probe
pressure_ratio4 = P04/p4
print("Pressure Ratio:",pressure_ratio4)
M4 = pitot.get_mach_from_pitot(pressure_ratio4,gamma,verbose=False)
print("Mach at 4: ",M4)

# get M_2 from 2theta and v(M_4) slide 28 of 6.2
params = [gamma,2*theta,M4]
M2 = external.get_mach_before_expansion(params)
print("Mach at 2: ",M2)

# newton solver for M_1
eps = 0.001
tol = 1e-6
M1guess = 3
realaxation = 0.6
maxIter = 100
delta = 10
iter = 0
turningAngle = theta+alpha

while abs(delta)>tol:
    # forward difference to get derivative
    params = [gamma,turningAngle,M1guess]
    M2guess, beta = external.get_mach_after_oblique_shock(params)
    params = [gamma,turningAngle,M1guess+eps]
    M2guessL, beta = external.get_mach_after_oblique_shock(params)
    deriv = ((M2guessL-M2) - (M2guess-M2))/eps
    # newton's method
    delta = (M2guess-M2)/deriv
    M1guess -= realaxation*delta
    if iter > maxIter:
        print("!!!! M1 solver did not converge")
        break
    else:
        iter+=1

print("Last guess at M2",M2guess)

print("Freestream Mach Number: ",M1guess)

print("\n\nProblem 2)")
gamma = 1.4
C_P0 = -0.94
alpha = np.deg2rad(1.5)
Gamma = np.deg2rad(30)
M_inf = 0.6

# a) M_crit unswept
params = [gamma,C_P0,"PG",alpha,Gamma]
M_crit_sweep,M_crit = cm.get_critical_mach(params,True)
print("Critical Mach for unswept wing",M_crit)
    
# b) M_crit swept
print("Critical Mach for swept wing",M_crit_sweep)
#c) compare t/c for unswept and swept
ToC = 0.12
print("unswept t/c",ToC)
ToC_sweep = ToC*np.cos(Gamma)
print("swept t/c",ToC_sweep)
#d) conclusion
print("effective t/c drops for swept wings and the critical mach increases")

print("\n\nProblem 3)")
gamma = 1.4
Gamma = 0

alphas = [0,1,2,3,4,5]
alphas = np.deg2rad(alphas)
CP0s = [0.05,-0.01,-0.15,-0.55,-1.1,-1.9]
CP0s = [-0.21,-0.48,-1.05,-2.0,-3.5,-5.0]
M_crit_pg = np.zeros(6)
M_crit_kt = np.zeros(6)
M_crit_lr = np.zeros(6)
for i in range(6):
    alpha = alphas[i]
    C_P0 = CP0s[i]
    params = [gamma,C_P0,"PG",alpha,Gamma]
    M_crit_sweep,M_crit = cm.get_critical_mach(params,True)
    M_crit_pg[i] = M_crit_sweep
    params = [gamma,C_P0,"KT",alpha,Gamma]
    M_crit_sweep,M_crit = cm.get_critical_mach(params,True)
    M_crit_kt[i] = M_crit_sweep
    params = [gamma,C_P0,"LR",alpha,Gamma]
    M_crit_sweep,M_crit = cm.get_critical_mach(params,True)
    M_crit_lr[i] = M_crit_sweep

plt.plot(np.rad2deg(alphas),M_crit_pg,label="PG")
plt.plot(np.rad2deg(alphas),M_crit_kt,label="KT")
plt.plot(np.rad2deg(alphas),M_crit_lr,label="LR")
plt.xlabel("Angle of attack (deg)")
plt.ylabel("Critical Mach Number")
plt.legend()
plt.show()
print("Using data given in problem statement, the zero angle of attack case did non intersect.  C_Pmin values were recalculated using XFOIL.  Attached plot ")



