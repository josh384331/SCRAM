import numpy as np
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
            return x,True

    # if max iterations then print this and retrun x
    if verbose:
        print("Newton Solver did not converge after {0} iterations".format(max_iteration))
    return x,False


def f_pm(M,params):
    theta = params[1]
    M_1 = params[2]
    return -theta - prandtl_meyer(M_1,params) + prandtl_meyer(M,params)

def prandtl_meyer(M,params):
    gamma = params[0]
    return ((gamma+1)/(gamma-1))**0.5*np.arctan(((gamma-1)/(gamma+1)*(M**2-1))**0.5)-np.arctan((M**2-1)**0.5)

 # isentropic pitot derivative   
def f_prime_pm(M,params):
    gamma = params[0]
    return 1/ M * ((M**2-1)**0.5)/(1+(gamma-1)/2*M**2)

def get_mach_after_expansion(params):
    '''This function returns the mach number after a prandtl-meyer expansion fan when given the gamma, turning angle, and mach before the fan (params[1-3])'''
    M_1 = params[2]
    if M_1 > 1:
        guess = M_1
    else:
        print("There are no expansion fans for subsonic flows")
        return M_1
    M_2,converged = newton_solver(f_pm,f_prime_pm,guess,params,verbose=False)
    if not converged:
        print("get_mach_after_expansion did not converge")
    return M_2

def get_params_across_expansion(params):
    '''This function computes the mach, pressure, and temperature after an oblique shock given gamma, theta, mach before shock, pressure before shock and temperature before shock

        Input
        ---------
        params = [gamma,theta_1,M_inf,p_inf,T_inf]
        
        Output
        ---------
        M_2, p_2, T_2'''
    gamma = params[0]
    theta = params[1]
    M_1 = params[2]
    p_1 = params[3]
    T_1 = params[4]
    M_2 = get_mach_after_expansion(params)

    # find temperatures after the expansion
    T_0 = T_1 * (1+(gamma-1)/2 * M_1**2)
    T_2 = T_0 / (1+(gamma-1)/2 * M_2**2)
    # find pressure after expansion
    P_0 = p_1 * (1+(gamma-1)/2 * M_1**2)**(gamma/(gamma-1))
    p_2 = P_0 / (1+(gamma-1)/2 * M_2**2)**(gamma/(gamma-1))
    # find density after expansion
    

    return M_2, p_2, T_2


def f_pm_backwards(M,params):
    '''This equation is used for the top of the newton solver of get mach before expansion'''
    theta = params[1]
    M_2 = params[2]
    return theta + prandtl_meyer(M,params) - prandtl_meyer(M_2,params)


def get_mach_before_expansion(params):
    '''This function returns the mach number before a prandtl-meyer expansion fan when given the gamma, turning angle, and mach after the fan (params[1-3])'''
    M_2 = params[2]
    if M_2 > 1:
        guess = M_2
    else:
        print("There are no expansion fans for subsonic flows")
        exit(-1)
    M_1 = newton_solver(f_pm_backwards,f_prime_pm,guess,params)
    return M_1

def get_mach_after_oblique_shock(params,strong=False):
    '''Computes the mach number behind and oblique shock as well as the angle of the shock given gamma, theta, and M_1'''
    gamma = params[0]
    theta = params[1]
    M_1 = params[2]
    beta = get_oblique_shock_angle(params,strong)
    M_n2 = ((1+(gamma-1)/2*(M_1*np.sin(beta))**2)/(gamma*(M_1*np.sin(beta))**2-(gamma-1)/2))**0.5
    # M_2 = ((M_1*np.cos(beta))**2*T_1/T_2+(1+(gamma-1)/2*(M_1*np.sin(beta)))/(gamma*(M_1*np.sin(beta))**2-(gamma-1)/2))**0.5
    M_2 = M_n2/(np.sin(beta-theta))
    
    
    return M_2,beta


def get_oblique_shock_angle_newton(params):
    '''This function outputs the oblique shock angle given gamma, theta, and M_1 in the params input'''
    # extract inputs
    gamma = params[0]
    theta = params[1]
    M_1 = params[2]

    # calc constants
    a = (1+(gamma-1)/2*M_1**2)*np.tan(theta)
    b = (M_1**2-1)
    c = (1+(gamma+1)/2*M_1**2)*np.tan(theta)
    params.append(a)
    params.append(b)
    params.append(c)
    guess = 0.8
    t_b = newton_solver(f_oblique_wave_angle,f_prime_oblique_wave_angle,guess,params,verbose=False)
    beta = np.arctan(t_b)

    return beta


def f_oblique_wave_angle(x,params):
    '''Cost function for newton's method to find beta'''
    gamma = params[0]
    theta = params[1]
    M_1 = params[2]
    a = params[3]
    b = params[4]
    c = params[5]
    return a*x**3-b*x**2+c*x+1

def f_prime_oblique_wave_angle(x,params):
    '''derivative of cost function for newton's method to find beta'''
    gamma = params[0]
    theta = params[1]
    M_1 = params[2]
    a = params[3]
    b = params[4]
    c = params[5]
    return 3*a*x**2-2*b*x**2+c

def get_oblique_shock_angle(params,strong=False):
    '''Computes the oblique shock angle given gamma, theta (turning angle), and M_1.  The strong tag is used to specify if the output will be the strong shock angle or not.'''
    gamma = params[0]
    theta = params[1]
    M_1 = params[2]
    if theta == 0:
        beta = np.arcsin(1/M_1)
    else:
        a = (1+(gamma-1)/2*M_1**2)*np.tan(theta)
        b = (M_1**2-1)
        c = (1+(gamma+1)/2*M_1**2)*np.tan(theta)
        lam = (b**2 - 3*a*c)**0.5
        chi = (b**3-9*a*(a+(gamma+1)/4*M_1**4*np.tan(theta)))/lam**3
        if strong:
            delta = 0
        else:
            delta = 1
        t_b = (b + 2 * lam*np.cos((4*np.pi*delta + np.arccos(chi))/3))/(3*a)
        beta = np.arctan(t_b)
    return beta

    
def get_temperature_across_oblique_shock(gamma,beta,T_1,M_1):
    """This function finds """
    left = (1 + 2*gamma/(gamma+1)*((M_1*np.sin(beta))**2-1))
    right = ((2+(gamma-1)*(M_1*np.sin(beta))**2)/((gamma+1)*(M_1*np.sin(beta))**2))
    return T_1 * left * right

def get_params_across_oblique_shock(params):
    '''This function computes the mach, pressure, and temperature after an oblique shock given gamma, theta, mach before shock, pressure before shock and temperature before shock

        Input
        ---------
        params = [gamma,theta_1,M_inf,p_inf,T_inf]
        
        Output
        ---------
        M_2, p_2, T_2'''
    gamma = params[0]
    theta = params[1]
    M_1 = params[2]
    p_1 = params[3]
    T_1 = params[4]
    M_2, beta = get_mach_after_oblique_shock(params)

    # find normal mach before shock
    M_n1 = M_1 * np.sin(beta)

    # get pressure across shock
    p_2 = p_1 * (1+2*gamma/(gamma+1)*(M_n1**2-1))
    P_02 = p_2 * (1+(gamma-1)/2 * M_2**2)**(gamma/(gamma-1)) 

    #get temperature across shock
    T_2 = T_1 * (1+2*gamma/(gamma+1)*(M_n1**2-1))*(2+(gamma-1)*M_n1**2)/((gamma+1)*M_n1**2)
    T_02 = T_2 * (1+(gamma-1)/2 * M_2**2)
    
    return M_2, p_2, T_2


if __name__ == "__main__":
    M_1 = 4
    p_1 = 0.01 * 101325
    T_1 = 217

    gamma = 1.25
    theta_1 = np.deg2rad(15)
    theta_2 = np.deg2rad(15)
    
    params = [gamma,theta_1,M_1]

    # calculate entry mach angle
    mu_1 = np.arcsin(1/M_1)
    # get mach after expansion
    M_2 = get_mach_after_expansion(params)
    # calcualte exit mach angle
    mu_2 = np.arcsin(1/M_2) - theta_1


    # find temperatures after the expansion
    T_0 = T_1 * (1+(gamma-1)/2 * M_1**2)
    T_2 = T_0 / (1+(gamma-1)/2 * M_2**2)
    # find pressure after expansion
    P_0 = p_1 * (1+(gamma-1)/2 * M_1**2)**(gamma/(gamma-1))
    p_2 = P_0 / (1+(gamma-1)/2 * M_2**2)**(gamma/(gamma-1))

    # print out results
    print("Problem 1")
    print("Results\n1-2")
    print(f"Entry Mach Angle[deg]: {round(np.rad2deg(mu_1),4)}°")
    print(f"Exit Mach Angle[deg]: {round(np.rad2deg(mu_2),4)}°")
    print("Mach_2:",round(M_2,4))
    print(f"Stagnation Pressure: {round(P_0,4)}Pa")
    print(f"Static Pressure: {round(p_2,4)}Pa")
    print(f"Stagnation Temperature: {round(T_0,4)}K")
    print(f"Static Temperature: {round(T_2,4)}K")

    
    params = [gamma,theta_1+theta_2,M_2]
    # params = [1.4,np.deg2rad(20),3]
    # p_2 = 101325
    # T_2 = 288
    # M_2 = 3
    # gamma = 1.4
    M_3, beta = get_mach_after_oblique_shock(params)

    # find normal mach before shock
    M_n2 = M_2 * np.sin(beta)

    # get pressure across shock
    p_3 = p_2 * (1+2*gamma/(gamma+1)*(M_n2**2-1))
    P_03 = p_3 * (1+(gamma-1)/2 * M_3**2)**(gamma/(gamma-1)) 

    #get temperature across shock
    T_3 = T_2 * (1+2*gamma/(gamma+1)*(M_n2**2-1))*(2+(gamma-1)*M_n2**2)/((gamma+1)*M_n2**2)
    T_03 = T_3 * (1+(gamma-1)/2 * M_3**2)
    print("2-3")
    print(f"Beta {round(np.rad2deg(beta),4)}°")
    print("Mach_3:",round(M_3,5))
    print(f"Stagnation Pressure: {round(P_03,4)}Pa")
    print(f"Static Pressure: {round(p_3,4)}Pa")
    print(f"Stagnation Temperature: {round(T_03,4)}K")
    print(f"Static Temperature: {round(T_3,4)}K")


    print("\nProblem 2\n")
    M_1 = 1.7
    p_1 = 0.26 * 101325 # Pa
    T_1 = 223
    gamma = 1.4
    theta_1 = np.deg2rad(15)
    theta_2 = 0
    r = 0.4 # m

    # set params
    params = [gamma,theta_1,M_1]
    # get mach after first shock and beta of first shock
    M_2, beta_1 = get_mach_after_oblique_shock(params)

    # find normal mach before shock
    M_n1 = M_1 * np.sin(beta_1)

    # get pressure across shock
    P_01 = p_1 * (1+(gamma-1)/2 * M_1**2)**(gamma/(gamma-1))
    p_2 = p_1 * (1+2*gamma/(gamma+1)*(M_n1**2-1))
    P_02 = p_2 * (1+(gamma-1)/2 * M_2**2)**(gamma/(gamma-1)) 

    #get temperature across shock
    T_01 = T_1 * (1+(gamma-1)/2 * M_1**2)
    T_2 = T_1 * (1+2*gamma/(gamma+1)*(M_n1**2-1))*(2+(gamma-1)*M_n1**2)/((gamma+1)*M_n1**2)
    T_02 = T_2 * (1+(gamma-1)/2 * M_2**2)
    
    print("Results\n")
    print(f"Stagnation Pressure before shock: {round(P_01,4)}Pa")
    print(f"Static Pressure before shock: {round(p_1,4)}Pa")
    print(f"Stagnation Temperature before shock: {round(T_01,4)}K")
    print(f"Static Temperature before shock: {round(T_1,4)}K")
    
    
    print("First Shock")
    print("Mach after shock:",round(M_2,4))
    print(f"Angle of first shock: {round(np.rad2deg(beta_1),4)}°")
    print(f"Stagnation Pressure after shock: {round(P_02,4)}Pa")
    print(f"Static Pressure after shock: {round(p_2,4)}Pa")
    print(f"Stagnation Temperature after shock: {round(T_02,4)}K")
    print(f"Static Temperature after shock: {round(T_2,4)}K")

    # get parameters after normal shock
    params = [gamma,theta_2,M_2]

    M_3, beta_2 = get_mach_after_oblique_shock(params,True)
    
    # get pressure across shock
    p_3 = p_2 * (1+2*gamma/(gamma+1)*(M_2**2-1))
    P_03 = p_3 * (1+(gamma-1)/2 * M_3**2)**(gamma/(gamma-1)) 

    #get temperature across shock
    T_3 = T_2 * (1+2*gamma/(gamma+1)*(M_2**2-1))*(2+(gamma-1)*M_2**2)/((gamma+1)*M_2**2)
    T_03 = T_3 * (1+(gamma-1)/2 * M_3**2)

    print("Results\nSecond Shock")
    print("Mach after shock:",round(M_3,4))
    print(f"Stagnation Pressure after shock: {round(P_03,4)}Pa")
    print(f"Static Pressure after shock: {round(p_3,4)}Pa")
    print(f"Stagnation Temperature after shock: {round(T_03,4)}K")
    print(f"Static Temperature after shock: {round(T_3,4)}K")

    # get pressure recovery
    print("\na) Stagnation Pressure Recovery",P_03/P_01)

    # get top to inlet distance for shock impingement 
    x = r/np.tan(beta_1)
    print("\nc)  X for oblique shock to align with inlet radius",x)

    num = 100
    thetas = np.linspace(1,45,num)
    pressure_ratios = np.zeros([num,1])

    for i in range(num):
        theta = np.deg2rad(thetas[i])
        params = [gamma,theta,M_1]
        # get mach after first shock and beta of first shock
        M_2, beta_1 = get_mach_after_oblique_shock(params)
        M_n1 = M_1 * np.sin(beta_1)
        p_2 = p_1 * (1+2*gamma/(gamma+1)*(M_n1**2-1))
        params = [gamma,theta_2,M_2]
        M_3, beta_2 = get_mach_after_oblique_shock(params,True)
        p_3 = p_2 * (1+2*gamma/(gamma+1)*(M_2**2-1))    
        P_03 = p_3 * (1+(gamma-1)/2 * M_3**2)**(gamma/(gamma-1)) 
        pressure_ratios[i] = P_03/P_01



    thetas_plot = []
    pressure_ratios_plot = []
    print(pressure_ratios)
    for i in range(len(pressure_ratios)):
        if not np.isnan(pressure_ratios[i]):
            pressure_ratios_plot.append(pressure_ratios[i])
            thetas_plot.append(thetas[i])

    pr = max(pressure_ratios_plot)
    I = np.argmax(pressure_ratios_plot) 
    # I = pressure_ratios_plot.argmax()

    print("Best Theta:",thetas_plot[I])
    plt.plot(thetas_plot,pressure_ratios_plot)
    plt.show()



    


