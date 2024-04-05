# 1D solution of the heat equation using the finite difference method

import numpy as np
import matplotlib.pyplot as plt
import Supersonic_external_flow
import wedgeWing
import time
from scipy.optimize import minimize
import standardAtmosphere


def logo():
    print("\n")
    print("┌────────────────────────────────────────────────────────────────────────────────────┐")
    print("│                                                                                    │")
    print("│    O~~ ~~           O~~        O~~~~~~~               O~            O~~       O~~  │")
    print("│  O~~    O~~      O~~   O~~     O~~    O~~            O~ ~~          O~ O~~   O~~~  │")
    print("│   O~~           O~~            O~~    O~~           O~  O~~         O~~ O~~ O O~~  │")
    print("│     O~~         O~~            O~ O~~              O~~   O~~        O~~  O~~  O~~  │")
    print("│        O~~      O~~            O~~  O~~           O~~~~~~ O~~       O~~   O~  O~~  │")
    print("│  O~~    O~~      O~~   O~~     O~~    O~~        O~~       O~~      O~~       O~~  │")
    print("│    O~~ ~~          O~~~~       O~~      O~~     O~~         O~~     O~~       O~~  │")
    print("│                                                                                    │")
    print("└────────────────────────────────────────────────────────────────────────────────────┘")
    print("                                     Supersonic Conduction & Radiative Analysis Model ")
    print("                                                                    by Joshua Hurwitz ")

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
    
    return h, T_aw, C_H
def getQTotal(T_W, eps, sig, M_e, rho_e, v_e, T_e, gamma, T_inf, M_inf,location):
    # get h
    h,T_aw, C_H = getConvectiveHeatTransferCoefficient(T_W, M_e, rho_e, v_e, T_e, gamma, T_inf, M_inf,location)
    # get radiation heating
    q_rad = eps * sig * T_W**4
    # get convective heating
    q_conv = h * (T_aw-T_W)
    # get total heating
    q_total = (q_conv-q_rad)**2
    return q_total

def getRadiationEquilibrium(eps, location, topEdge,length_list,wedgeAngles,freestream):
    """
    Calculates the radiation equilibrium temperature for a given set of parameters.

    Parameters:
    eps (float): Emissivity of the surface.
    location (str): Location of the surface.
    topEdge (float): Top edge of the surface.
    length_list (list): List of lengths.
    wedgeAngles (list): List of wedge angles.
    freestream (list): List of freestream parameters:
        - M_inf (float): Freestream Mach number.
        - AOA (float): Angle of attack in radians.
        - gamma (float): Specific heat ratio.
        - altitude (float): Altitude in m.

    Returns:
    float: The radiation equilibrium temperature.
    """
    sig = 5.67e-8 # W/m^2K^4
    # get trajectory parameters
    M_inf = freestream[0]
    AOA = freestream[1]
    gamma = freestream[2]
    altitude = freestream[3]
    # get freestream parameters
    [Z,T_inf,p_inf,rho_inf,a_inf] = standardAtmosphere.get_atmospheric_properties_si(altitude)

    # get edge parameters
    M_e, p_e, T_e, rho_e, v_e = get_edge_params(M_inf, gamma, wedgeAngles, p_inf, T_inf, AOA, location, length_list, topEdge)
    
    T_W = 1
    # Minimize the objective function subject to the constraints
    result = minimize(getQTotal, T_W, method='SLSQP', args=(eps, sig, M_e, rho_e, v_e, T_e, gamma, T_inf, M_inf, location), constraints={'type': 'ineq', 'fun': lambda T_W: T_W})

    # return the result
    return result.x
    
   
        

    
def heatEquation1D(T0,Nx, Nt, endTime, eps, location, topEdge,rhoList,cpList, length_list,wedgeAngles,kList,thicknessList,freestream):
    """
    Solves the 1D heat conduction equation numerically using finite difference method.

    Parameters:
    T0 (float): Initial temperature in K.
    Nx (int): Number of spatial steps.
    Nt (int): Number of time steps.
    endTime (float): End time of simulation in s.
    eps (float): Emissivity of the material.
    location (float): Location of the point of interest in m from leading edge .
    topEdge (bool): whether the point is on the top edge or not.
    rhoList (float, list): List of densities of different layers of the thickness in kg/m^3.
    cpList (float, list): List of specific heat capacities of different layers of the thickness in J/kgK.
    length_list (float, list): List of lengths of different sections of the geometry in m.
    wedgeAngles (float, list): List of wedge angles of different sections of the geometry in radians.
    kList (float, list): List of thermal conductivities of different layers of the thickness in W/mK.
    thicknessList (float, list): List of thicknesses of different sections of the layers of the skin in m.
    freestream (float, list): List of freestream parameters:
        - M_inf (float): Freestream Mach number.
        - AOA (float): Angle of attack in radians.
        - gamma (float): Specific heat ratio.
        - altitude (float): Altitude in m.
    

    Returns:
    numpy.ndarray: Temperature distribution matrix.
    """
    # constants
    sig = 5.67e-8 # W/m^2K^4
    beta = 4 # most conservative value
    thickness = sum(thicknessList)
    trajectory_Completed = False
    x = location

    # Initialize space and time steps
    dx = thickness/Nx
    # check if dx is too large
    if dx > min(thicknessList):
        dx = min(thicknessList)/2
        print('ERROR: dx too large, setting to half minimum thickness')

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
            print(f"Completed\n Trajectory used {p} time steps")
            trajectory_Completed = True
            break
        # get trajectory parameters
        M_inf = freestream[0]
        AOA = freestream[1]
        gamma = freestream[2]
        altitude = freestream[3]
        # get freestream parameters
        [Z,T_inf,p_inf,rho_inf,a_inf] = standardAtmosphere.get_atmospheric_properties_si(altitude)

        # get edge parameters
        M_e, p_e, T_e, rho_e, v_e = get_edge_params(M_inf, gamma, wedgeAngles, p_inf, T_inf, AOA, x, length_list, topEdge)
    
        # get h
        h,T_aw, C_H = getConvectiveHeatTransferCoefficient(T[p,0], M_e, rho_e, v_e, T_e, gamma, T_inf, M_inf,x)

        # get material properties for surface
        k = kList[0]
        rho = rhoList[0]
        cp = cpList[0]

    
        
        # get time step
        ## check each layer
        timeSteps = np.zeros(len(thicknessList))
        for i in range(0, len(thicknessList)):
            Bi_L = h * thickness / kList[i]
            timeSteps[i] = dx**2/(2*(kList[i]/(rhoList[i]*cpList[i]))*(1+dx*sig*eps*beta*(T[p,0])**3/kList[i]+Bi_L))
        ## pick the smallest time step
        
        #dt = min(timeSteps) #this makes a really small time step
        dt = (timeSteps[0])

        # save current time
        timeList[p+1] = timeList[p] + dt
        # get forier number for surface
        Fo = k/(rho*cp) * dt / dx**2
        # get boundary conditions
        ## outer boundary
        T[p+1, 0] = dt * (h * (T_aw-T[p,0]) - eps * sig * T[p,0]**4 )/ (rho * cp * dx) + Fo * (T[p,1]-T[p,0]) + T[p,0]
        # print((eps * sig * T[p,0]**4),"    ",(h * (T_aw-T[p,0])))
        if T[p+1, 0] < 0:
            print("CRITICAL ERROR: Outer Boundary Temperature is negative")
            return timeList, T
        ## inner boundary
        Foi = kList[-1]/(rhoList[-1]*cpList[-1]) * dt / dx**2
        FoAir = 0.025/(1.23*1005) * dt / dx**2
        T_inf_in = 300
        T[p+1,-1] = Foi*2*T[p,-2]+(1-2*Foi)*T[p,-1]
        #T[p+1,-1] = 0.025*dt*(T_inf_in-T[p,-1])/(1.23*1005*0.0254**2) + 2*Foi*(T[p,-2]-T[p,-1]) + T[p,-1]
        if T[p+1, -1] < 0:
            print("CRITICAL ERROR: Inner Boundary Temperature is negative")
            return timeList, T
        

        # Space loop
        layerCounter = 0
        for n in range(1, Nx-1):
            
            if n*dx > sum(thicknessList[:layerCounter+1]):
                
                Fo_old = Fo
                layerCounter += 1
                k = kList[layerCounter]
                rho = rhoList[layerCounter]
                cp = cpList[layerCounter]
                Fo = k/(rho*cp) * dt / dx**2
                T[p+1, n] = Fo*T[p,n+1]+ Fo*T[p,n-1] + (1-2*Fo)*T[p,n]
            else:
                T[p+1, n] = Fo*(T[p,n+1]+T[p,n-1]) + (1-2*Fo)*T[p,n]
            
            if T[p+1, n] < 0:
                print("CRITICAL ERROR: Inner Temperature is negative")
                return timeList, T

        # Create new array with timeList and T
    if not trajectory_Completed:
        print("Trajectory not completed")
        return timeList, T
    else:
        return timeList[:p], T[:p,:]










if __name__ == "__main__":
    # start time
    start_time = time.time()
    logo()
    # Parameters
    ## TPS
    # kList = [9.28,0.0402] # W/mK
    # thicknessList = [0.001127, 0.05080] # m
    # rhoList = [4540, 30] # kg/m^3
    # cpList = [585.76, 900] # J/kgK
    kList = [0.041153] # W/mK
    thicknessList = [0.05207] # m
    rhoList = [140] # kg/m^3
    cpList = [888.37] # J/kgK
    eps = 0.75
    ## Geometry
    wedgeAngles = np.deg2rad([6, -6, 0, 0]) # degrees
    length_list = [13, 13, 13, 13] # m
    x = 12 # m
    topEdge = False
    ## Flow Parameters
    M_inf = 3
    gamma = 1.4
    alpha = np.deg2rad(4)
    altitude = 29000 # m
    freestream = [M_inf, alpha, gamma, altitude]
    ## settings
    T0 = 300 # K
    Nx = 100 # number of spatial steps
    Nt = 1000000 # number of time steps
    endTime = 1200 # s
    # get radiation equilibrium
    T_0inf = T_inf * (1+ (gamma-1)/2 * M_inf**2)
    RadEquilTemp = getRadiationEquilibrium(eps, x, topEdge,length_list,wedgeAngles,freestream)
    print("Radiation Equilibrium Temperature: ", RadEquilTemp, "K")
    print("Total Temperature: ", T_0inf, "K")
    
    # heatEquation1D(T0,Nx, Nt, endTime, eps, location, topEdge,rhoList,cpList, length_list,wedgeAngles,kList,thicknessList)
    print("Running Thermal Analysis Model for given trajectory...",end="",flush=True)
    timeList, T = heatEquation1D(T0, Nx, Nt, endTime, eps, x, topEdge, rhoList,cpList, length_list,wedgeAngles,kList,thicknessList,freestream)
    
    # output completion time
    end_time = time.time()
    execution_time = end_time - start_time
    hours = int(execution_time // 3600)
    minutes = int((execution_time % 3600) // 60)
    seconds = int(execution_time % 60)
    print(f"Execution time: {hours:02d}:{minutes:02d}:{seconds:02d}")
    # output results
    print("Max Exterior Temperature: ", np.max(T[:,0]), "K")
    print("Max Interior Temperature: ", np.max(T[:,-1]), "K")
    # save text file of results
    
    np.savetxt('time_list.txt', timeList[::10], delimiter='\t')
    np.savetxt('temperature_distribution.txt', T[::10,:], delimiter='\t')
    # end time
    end_save_time = time.time()
    save_time = end_save_time - end_time
    hours = int(save_time // 3600)
    minutes = int((save_time % 3600) // 60)
    seconds = int(save_time % 60)
    print(f"Saving time: {hours:02d}:{minutes:02d}:{seconds:02d}")
    