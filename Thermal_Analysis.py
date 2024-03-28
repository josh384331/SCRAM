# 1D solution of the heat equation using the finite difference method

import numpy as np
import matplotlib.pyplot as plt
import Supersonic_external_flow
import wedgeWing




def get_edge_params(M_inf, gamma, wedgeAngles, p_inf, T_inf, alpha, x, length_list, topEdge = False, R_air = 287.05):
    '''
    Calculate the edge parameters for a given set of inputs.

    where points are defined as

        1 2
    inf < >
        3 4

    Args:
        M_inf (float): Mach number at infinity.
        gamma (float): Specific heat ratio.
        wedgeAngles (float,list): Wedge angles in radians.
        p_inf (float): Freestream Pressure in Pa.
        T_inf (float): Freestream Temperature in K.
        alpha (float): Angle of attack in radians.
        x (float): Current x-coordinate in m.
        x_pt_list (float,list): length of each panel in m.
        topEdge (bool, optional): Flag indicating if it's the top edge. Defaults to False.
        R_air (float, optional): Specific gas constant for air. Defaults to 287.05.

    Returns:
        tuple: Tuple containing the edge parameters:
            - M_1 (float): Mach number at point 1.
            - p_1 (float): Pressure at point 1.
            - T_1 (float): Temperature at point 1.
            - M_2 (float): Mach number at point 2.
            - p_2 (float): Pressure at point 2.
            - T_2 (float): Temperature at point 2.
            - M_3 (float): Mach number at point 3.
            - p_3 (float): Pressure at point 3.
            - T_3 (float): Temperature at point 3.
            - M_4 (float): Mach number at point 4.
            - p_4 (float): Pressure at point 4.
            - T_4 (float): Temperature at point 4.
    '''
    # set up input
    params = [0,[],0,0,0,[]]
    params[0] = gamma
    params[1] = wedgeAngles
    params[2] = M_inf
    params[3] = p_inf
    params[4] = T_inf
    params[5] = alpha
    # run the wedge wing code
    M_1,p_1,T_1, M_2,p_2,T_2, M_3,p_3,T_3, M_4,p_4,T_4 = wedgeWing.getArbWingParams(params)
    
    # get the edge density assuming ideal gas
    rho_1 = p_1/(R_air*T_1)
    rho_2 = p_2/(R_air*T_2)
    rho_3 = p_3/(R_air*T_3)
    rho_4 = p_4/(R_air*T_4)

    # get the local velocity at each of the panels assuming calorically perfect
    v_1 = M_1 * np.sqrt(gamma*p_1/rho_1)
    v_2 = M_2 * np.sqrt(gamma*p_2/rho_2)
    v_3 = M_3 * np.sqrt(gamma*p_3/rho_3)
    v_4 = M_4 * np.sqrt(gamma*p_4/rho_4)

    if (x < length_list[0]) or (x < length_list[3]):
        if topEdge:
            return M_1, p_1, T_1, rho_1, v_1
        else:
            return M_3, p_3, T_3, rho_3, v_3
    else:
        if topEdge:
            return M_2, p_2, T_2, rho_2, v_2
        else:
            return M_4, p_4, T_4, rho_4, v_4


    
    
def getConvectiveHeatTransferCoefficient(T_w, M_e, rho_e, v_e, T_e, gamma, T_inf, M_inf,x, R_air = 287.05):
    """
    Calculates the convective heat transfer coefficient for a flat plate in supersonic flow.

    Args:
        T_w (float): Wall temperature in K.
        M_e (float): Mach number at the location.
        rho_e (float): Density at the location in kg/m^3.
        v_e (float): Velocity at the location in m/s.
        T_e (float): Temperature at the location in K.
        gamma (float): Specific heat ratio.
        T_inf (float): Free-stream temperature in K.
        M_inf (float): Free-stream Mach number.
        x (float): Distance along the flat plate in m.

    Returns:
        float: Convective heat transfer coefficient.
    """
    # get properties at the edge
    cpe = gamma * R_air / (gamma - 1) # J/kgK
    mu_e = 1.716*10**(-5) * (T_e/273.1)**(3/2) * 383.1/(T_e+110) # kg/(m*s)
    # get reynolds and prandtl numbers
    Re_x = rho_e * v_e * x / mu_e
    #Pr_e = mu_e * cp / k_e
    Pr_e = 0.72 # as directed for air at 300K

    # get reynolds analogy r_factor
    if Re_x < 5e5:
        r_factor = Pr_e**(1/2)
    else:
        r_factor = Pr_e**(1/3)

    # get adiabatic wall temperature
    T_aw = T_e*(1+r_factor*(gamma-1)/2*M_e**2)

    # get laminar and turbulent boundary layer parameters
    if Re_x < 5e5:
        C_fi = 0.664/Re_x**0.5
        delta_i = 4.64*x/Re_x**0.5
        T_ref = T_e*(0.5+0.037*M_e**2+0.5*T_w/T_e)
        mu_ref = mu_e*(T_ref/T_e)**1.5*(T_e+110)/(T_ref+110)
        theta_ratio = (mu_ref/mu_e)**(0.5)/(T_ref/T_e)**(0.5)
    else:
        C_fi = 0.0592/Re_x**0.2
        delta_i = 0.37*x/Re_x**0.2
        T_ref = T_e*(0.5+0.037*M_e**2+0.5*T_w/T_e)
        mu_ref = mu_e*(T_ref/T_e)**1.5*(T_e+110)/(T_ref+110)
        theta_ratio = (mu_ref/mu_e)**(0.2)/(T_ref/T_e)**(0.8)

    # get skin friction coefficient
    C_f = C_fi * theta_ratio

    # get total temperature and zeta factor
    T_0inf = T_inf * (1+ (gamma-1)/2 * M_inf**2)
    zeta_w = T_w/T_0inf

    # get Reynolds analogy factor
    if zeta_w < 0.2:
        F_RA = 1
    elif zeta_w <= 0.65:
        F_RA = 0.8311+0.9675*zeta_w-0.6142 * zeta_w**2
    else:
        F_RA = 1.2
    
    # get stanton number
    C_H = F_RA / 2 * C_f
    

    # get heat transfer coefficient
    h = rho_e * cpe * v_e * C_H
    return h, T_aw

   
    
def heatEquation1D(dx, T0, Nt, endTime, eps, location, topEdge,rho,cp, length_list,wedgeAngles,k,thicknessList):
    """
    Solves the 1D heat conduction equation numerically using finite difference method.

    Parameters:
    dx (float): Spatial step size in m.
    T0 (float): Initial temperature in K.
    Nt (int): Number of time steps.
    endTime (float): End time of simulation in s.
    eps (float): Emissivity of the material.
    location (float): Location of the point of interest in m from leading edge .
    topEdge (bool): whether the point is on the top edge or not.
    rho (float): Density of the material in kg/m^3.
    cp (float): Specific heat capacity of the material in J/kgK.
    length_list (list): List of lengths of different sections of the geometry in m.
    wedgeAngles (list): List of wedge angles of different sections of the geometry in radians.
    k (float): Thermal conductivity of the material.
    thicknessList (list): List of thicknesses of different sections of the geometry in m.

    Returns:
    numpy.ndarray: Temperature distribution matrix.
    """
    # constants
    sig = 5.67e-8 # W/m^2K^4
    beta = 4 # most conservative value
    x = location
    thickness = thicknessList
    alpha = k/(rho*cp)

    # Number of points
    Nx = int(thickness/dx) + 1

    # Initialize temperature matrix
    T = np.zeros((Nt, Nx))
    # initialize time list
    timeList = np.zeros(Nt)

    # Initial condition
    T[0, :] = T0

    
    #T[p,n]
    # Time loop
    for p in range(0, Nt-1):
        # check if max time is reached
        if timeList[p-1] > endTime:
            break
        # get trajectory parameters
        M_inf = 5
        gamma = 1.4
        p_inf = 1440.396
        T_inf = 225.61674
        AOA = np.deg2rad(4)

        # get edge parameters
        M_e, p_e, T_e, rho_e, v_e = get_edge_params(M_inf, gamma, wedgeAngles, p_inf, T_inf, AOA, x, length_list, topEdge)
        
        # get h
        h,T_aw = getConvectiveHeatTransferCoefficient(T[p,0], M_e, rho_e, v_e, T_e, gamma, T_inf, M_inf,x)
        
        # get Biot number
        Bi_L = h * thickness / k
        # get time step
        dt = dx**2/(2*alpha*(1+dx*sig*eps*beta*(T[p,0])**3/k+Bi_L))
        # save current time
        timeList[p] = timeList[p-1] + dt
        # get forier number
        Fo = alpha * dt / dx**2
        # get boundary conditions
        ## outer boundary
        T[p+1, 0] = dt * (h * (T_aw-T[p,0]) - eps * sig * T[p,0]**4 )/ (rho * cp * dx) - Fo *(T[p,1]-T[p,0]) + T[p,0]

        ## inner boundary
        T[p+1,-1] = Fo*2*T[p,-2]+(1-2*Fo)*T[p,-1]
        

        # Space loop
        for n in range(1, Nx-1):
            T[p+1, n] = Fo*(T[p,n+1]+T[p,n-1]) + (1-2*Fo)*T[p,n]

        # Create new array with timeList and T
        result = np.column_stack((timeList, T))
    return result










if __name__ == "__main__":
    # Parameters
    ## TPS
    kList = [0.0312, 0.5, 0.5] # W/mK
    thicknessList = [0.1, 0.01, 0.01] # m
    rhoList = [3352, 1000, 1000] # kg/m^3
    cpList = [938, 1000, 1000] # J/kgK
    eps = 0.75
    ## Geometry
    wedgeAngles = np.deg2rad([6, -6, 0, 0]) # degrees
    length_list = [2, 2, 2, 2] # m
    x = 1
    topEdge = False
    ## Flow Parameters
    M_inf = 5
    gamma = 1.4
    p_inf = 1440.396
    T_inf = 225.61674
    alpha = 4
    ## settings
    T0 = 300 # K
    dx = thicknessList[0]/100 # m
    Nt = 10000
    endTime = 4000 # s
    
    # heatEquation1D(dx, T0, Nt, endTime, eps, location, topEdge, rho, cp, length_list,wedgeAngles,k,thicknessList):
    
    result = heatEquation1D(dx, T0, Nt, endTime, eps, x, topEdge, rhoList[0],cpList[0], length_list,wedgeAngles,kList[0],thicknessList[0])

    # save text file of results
    np.savetxt('temperature_distribution.txt', result, delimiter='\t')
    