import numpy as np
import sympy as sp
from matplotlib import pyplot as plt

def newton_solver(f,f_prime,guess,params,tol=1e-8,max_iteration=100,verbose=True,debug=False):
    # set initial guess
    x = guess

    # run newton solver loop
    for i in range(max_iteration):
        # update law
        delta_x = f(x,params) / f_prime(x,params)
        x -= delta_x
        if debug: print(x,delta_x)
        # end conditions
        if abs(delta_x) < tol:
            if verbose:
                print("Converged to {0} after {1} iterations".format(x,i+1))
            return x

    # if max iterations then print this and retrun x
    print("Did not converge after {0} iterations".format(max_iteration))
    return x


def f_crit(M,params):
    '''this function is used to find the critical mach number using different pressure correction methods'''
    gamma = params[0]
    C_P0 = params[1]
    correction = params[2]
    base = 2/(gamma*M**2)*(((1+(gamma-1)/2*M**2)/((gamma+1)/2))**(gamma/(gamma-1)) - 1)
    pg_correction = - C_P0/((1-M**2)**(0.5))
    kt_correction = - C_P0/((1-M**2)**0.5+M**2/(1+(1-M**2)**0.5)*C_P0/2)
    lr_correction = - C_P0/((1-M**2)**0.5+(1+(gamma-1)/2*M**2)*M**2/(1+(1-M**2)**0.5)*C_P0/2)
    if correction == "PG":
        return base + pg_correction
    elif correction == "KT":
        return base + kt_correction
    elif correction == "LR":
        return base + lr_correction 


 # isentropic pitot derivative   
def f_prime_crit(M,params):
    '''this function is used to find the critical mach number using different pressure correction methods'''
    gamma = params[0]
    C_P0 = params[1]
    correction = params[2]
    base = 4*((M**2*(gamma/2 - 1/2) + 1)/(gamma/2 + 1/2))**(gamma/(gamma - 1))*(gamma/2 - 1/2)/(M*(gamma - 1)*(M**2*(gamma/2 - 1/2) + 1)) - 4*(((M**2*(gamma/2 - 1/2) + 1)/(gamma/2 + 1/2))**(gamma/(gamma - 1)) - 1)/(M**3*gamma)
    pg_correction = -1.0*C_P0*M/(1 - M**2)**1.5
    kt_correction = -C_P0*(-0.5*C_P0*M**3/((1 - M**2)**0.5*((1 - M**2)**0.5 + 1)**2) - C_P0*M/((1 - M**2)**0.5 + 1) + 1.0*M/(1 - M**2)**0.5)/(C_P0*M**2/(2*((1 - M**2)**0.5 + 1)) + (1 - M**2)**0.5)**2
    lr_correction = -C_P0*(-0.5*C_P0*M**3*(M**2*(gamma/2 - 1/2) + 1)/((1 - M**2)**0.5*((1 - M**2)**0.5 + 1)**2) - C_P0*M**3*(gamma/2 - 1/2)/((1 - M**2)**0.5 + 1) - C_P0*M*(M**2*(gamma/2 - 1/2) + 1)/((1 - M**2)**0.5 + 1) + 1.0*M/(1 - M**2)**0.5)/(C_P0*M**2*(M**2*(gamma/2 - 1/2) + 1)/(2*((1 - M**2)**0.5 + 1)) + (1 - M**2)**0.5)**2
    if correction == "PG":
        return base + pg_correction
    elif correction == "KT":
        return base + kt_correction
    elif correction == "LR":
        return base + lr_correction
    # return -1.0*C_P0*M/(1 - M**2)**1.5 + 4*((M**2*(gamma/2 - 1/2) + 1)/(gamma/2 + 1/2))**(gamma/(gamma - 1))*(gamma/2 - 1/2)/(M*(gamma - 1)*(M**2*(gamma/2 - 1/2) + 1)) - 4*(((M**2*(gamma/2 - 1/2) + 1)/(gamma/2 + 1/2))**(gamma/(gamma - 1)) - 1)/(M**3*gamma)

def get_critical_mach(params,get_unswept_mach=False):
    '''This function solves for the critical mach number using the input parameters gamma, minimum cp(from xfoil), correction method [PG,KT,LR] reffering to prandtl glauert, Karman-Tsien and Laitone’s rules.  Angle of attack and sweep values in the form params[0-4] '''
    gamma = params[0]
    C_P0 = params[1]
    correction = params[2]
    alpha = params[3]
    Gamma = params[4]
    M_crit = newton_solver(f_crit,f_prime_crit,0.5,params,verbose=False)
    M_crit_sweep = M_crit/((1-(np.sin(Gamma))**2*(np.cos(alpha))**2)**0.5)
    if get_unswept_mach:
        return M_crit_sweep,M_crit
    else:
        return M_crit_sweep


if __name__ == "__main__":
    # M = sp.symbols("M")
    # gamma = sp.symbols("gamma")
    # C_P0 = sp.symbols("C_P0")
    # expression = 2/(gamma*M**2)*(((1 + ((gamma-1)/2)*M**2)/((gamma+1)/2))**(gamma/(gamma-1)) - 1)
    # pg = - C_P0/((1-M**2)**(0.5))
    # kt = - C_P0/((1-M**2)**0.5+M**2/(1+(1-M**2)**0.5)*C_P0/2)
    # lr = - C_P0/((1-M**2)**0.5+(1+(gamma-1)/2*M**2)*M**2/(1+(1-M**2)**0.5)*C_P0/2)

    # # func     = ((2/(gamma*M**2))*((1 + ((gamma-1)/2)*M**2)/((gamma+1)/2))**(gamma/(gamma-1)) - 1) - CPo/((1-M**2)**0.5)
    # derivative = sp.diff(expression,M)
    # pg_dir = sp.diff(pg,M)
    # kt_dir = sp.diff(kt,M)
    # lr_dir = sp.diff(lr,M)
    # # derivative = sp.simplify(derivative)
    # print("Expression",expression)
    # print("Derivative", derivative)
    # print("pg der",pg_dir)
    # print("kt der",kt_dir)
    # print("lr der",lr_dir)

    alpha = np.deg2rad(1.5)
    Gamma = np.deg2rad(30)
    gamma = 1.4
    C_P0 = -0.94
    # params = [gamma,C_P0]
    # M_crit = newton_solver(f_crit,f_prime_crit,0.5,params)


    # # M_infs = np.linspace(0,1,100)
    # # iter = 0
    # # cp_crit = np.zeros(100)
    # # cp = np.zeros(100)
    # # for i in M_infs:
    # #     M = i
    # #     cp_crit[iter] =2/(gamma*M**2)*(((1 + ((gamma-1)/2)*M**2)/((gamma+1)/2))**(gamma/(gamma-1)) - 1)
    # #     cp[iter] = C_P0/((1-M**2)**0.5)
    # #     iter +=1
    
    # # plt.plot(M_infs,cp_crit,label="cp_crit")
    # # plt.plot(M_infs,cp,label="cp")
    # # plt.legend()
    # # plt.show()

   







