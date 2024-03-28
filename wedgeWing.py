import numpy as np
from matplotlib import pyplot as plt
import Rayleigh_Pitot_Solver as pitot
import Supersonic_external_flow as external


def getWedgeLD(params,visc=False,TurbulentOnly=False,compressible=True):
    '''This function computes the lift to drag ratio of a double wedge airfoil
    
    Inputs
    __________
    -params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alpha,c*,rho_inf*,mu_inf*,a_inf*]
    -visc = bool

    *{only required if visc==true}
    Outputs
    _________
    -LOD
    '''
    # extract params
    gamma = params[0]
    wedgeAngle = params[1]
    M_inf = params[2]
    p_inf = params[3]
    T_inf = params[4]
    alpha = params[5]
    

    M_1,p_1,T_1, M_2,p_2,T_2, M_3,p_3,T_3, M_4,p_4,T_4 = getWingParams(params)
    CD = 1/(gamma * M_inf**2 * np.cos(wedgeAngle)) * ((p_3 - p_2)/p_inf * np.sin(wedgeAngle+alpha) + (p_1 - p_4)/p_inf * np.sin(wedgeAngle-alpha))
    CL = 1/(gamma * M_inf**2 * np.cos(wedgeAngle)) * ((p_3 - p_2)/p_inf * np.cos(wedgeAngle+alpha) + (p_4 - p_1)/p_inf * np.cos(wedgeAngle-alpha))

    

    if visc and TurbulentOnly:
        Cs = 120
        c = params[6]
        rho_inf = params[7]
        mu_inf = params[8]
        a_inf = params[9]
        Re = rho_inf * a_inf * M_inf * c * np.cos(alpha)/np.cos(wedgeAngle)/ mu_inf
        CDF_inc_t = 14/(225*np.cos(wedgeAngle)*Re**(1/7))
        T_avg = T_inf * (1+(2/9)*((gamma-1)/2*M_inf**2))
        if compressible:
            CompCorrectionT = 1/(((T_inf/T_avg)**(5/2)*((T_avg+Cs)/(T_inf+Cs)))**(1/7))
        else:
            CompCorrectionT = 1
        return CL/(CD+CDF_inc_t*CompCorrectionT), Re
    elif visc:
        Cs = 120
        c = params[6]
        rho_inf = params[7]
        mu_inf = params[8]
        a_inf = params[9]
        Re = rho_inf * a_inf * M_inf * c * (np.cos(alpha)/np.cos(wedgeAngle))/ mu_inf
        Rf = 1

        if Re <= 500000:
            # incompressible laminar 
            CDf_inc = 1.328/((Re)**0.5)  * 2/np.cos(wedgeAngle) # removed this to compare with spencer
            T_avg = T_inf * (1 + (Rf-8/15) * (gamma-1)/2 * M_inf**2)
            if compressible:
                CompCorrectionL = ((T_avg/T_inf)**(5/2)*((T_inf+Cs)/(T_avg+Cs)))**0.5
            else:
                CompCorrectionL = 1
            return CL/(CD+CDf_inc*CompCorrectionL),Re
        else:
            CDF_inc_t = 14/(225*np.cos(wedgeAngle)*(Re**(1/7)))
            T_avg = T_inf * (1+(Rf-7/9)*((gamma-1)/2*M_inf**2))
            if compressible:
                CompCorrectionT = 1/(((T_inf/T_avg)**(5/2)*((T_avg+Cs)/(T_inf+Cs)))**(1/7))
            else:
                CompCorrectionT = 1
            return CL/(CD+CDF_inc_t*CompCorrectionT),Re
    else:
        return CL/CD
    
def getWingParams(params):
    '''This function computes the pressures temperature and mach numbers around a wedge airfoil 
    
    Inputs
    __________
    -params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alpha]

    Outputs
    __________
    -M_1,p_1,T_1, M_2,p_2,T_2, M_3,p_3,T_3, M_4,p_4,T_4

    where points are defined as

        1 2
    inf < >
        3 4
    '''
    # extract params
    gamma = params[0]
    wedgeAngle = params[1]
    M_inf = params[2]
    p_inf = params[3]
    T_inf = params[4]
    alpha = params[5]

    ### find top pressures 
    ## get pressures at point 1
    # get theta 
    theta_inf_1 = wedgeAngle - alpha
    
    # check if expansion fan or oblique shock and call appropriate function
    if theta_inf_1 < 0:
        params = [gamma,-theta_inf_1,M_inf,p_inf,T_inf]
        M_1, p_1, T_1 =  external.get_params_across_expansion(params)
    else:
        params = [gamma,theta_inf_1,M_inf,p_inf,T_inf]
        M_1, p_1, T_1 =  external.get_params_across_oblique_shock(params)

    ## get pressures at point 2
    # get theta 
    theta_1_2 = 2*wedgeAngle
    # set parameters and call function
    params = [gamma,theta_1_2, M_1, p_1, T_1]
    M_2, p_2, T_2 = external.get_params_across_expansion(params)

    ### find bottom pressures
    ## get pressures at point 3
    # get theta 
    theta_inf_3 = wedgeAngle + alpha
    
    # check if expansion fan or oblique shock and call appropriate function
    if theta_inf_3 < 0:
        params = [gamma,-theta_inf_3,M_inf,p_inf,T_inf]
        M_3, p_3, T_3 =  external.get_params_across_expansion(params)
    else:
        params = [gamma,theta_inf_3,M_inf,p_inf,T_inf]
        M_3, p_3, T_3 =  external.get_params_across_oblique_shock(params)
    ## get pressures at point 4
    # get theta 
    theta_3_4 = 2*wedgeAngle
    # set parameters and call function
    params = [gamma,theta_3_4, M_3, p_3, T_3]
    M_4, p_4, T_4 = external.get_params_across_expansion(params)

    return M_1,p_1,T_1, M_2,p_2,T_2, M_3,p_3,T_3, M_4,p_4,T_4


def getArbWingParams(params):
    '''This function computes the pressures temperature and mach numbers around a arbitrary wedge airfoil 
    
    Inputs
    __________
    -params = [gamma,wedgeAngles,M_inf,p_inf,T_inf,alpha]
    wedgeAngles = [wedgeAngle1,wedgeAngle2,wedgeAngle3,wedgeAngle4] 
    where the angles are in radians and are defined as the angle between the wedge and the fuselage reference line at each location.  Negative angles will create expansion fans
    Outputs
    __________
    -M_1,p_1,T_1, M_2,p_2,T_2, M_3,p_3,T_3, M_4,p_4,T_4

    where points are defined as

        1 2
    inf < >
        3 4
    '''
    # extract params
    gamma = params[0]
    wedgeAngle = params[1]
    M_inf = params[2]
    p_inf = params[3]
    T_inf = params[4]
    alpha = params[5]

    ### find top pressures 
    ## get pressures at point 1
    # get theta 
    theta_inf_1 = wedgeAngle[0] - alpha
    
    # check if expansion fan or oblique shock and call appropriate function
    if theta_inf_1 < 0:
        params = [gamma,-theta_inf_1,M_inf,p_inf,T_inf]
        M_1, p_1, T_1 =  external.get_params_across_expansion(params)
    else:
        params = [gamma,theta_inf_1,M_inf,p_inf,T_inf]
        M_1, p_1, T_1 =  external.get_params_across_oblique_shock(params)

    ## get pressures at point 2
    # get theta 
    theta_1_2 = wedgeAngle[1]
    # set parameters and call function
    if theta_1_2 < 0:
        params = [gamma,-theta_1_2, M_1, p_1, T_1]
        M_2, p_2, T_2 =  external.get_params_across_expansion(params)
    else:
        params = [gamma,theta_1_2, M_1, p_1, T_1]
        M_2, p_2, T_2 =  external.get_params_across_oblique_shock(params)


    ### find bottom pressures
    ## get pressures at point 3
    # get theta 
    theta_inf_3 = wedgeAngle[2] + alpha
    
    # check if expansion fan or oblique shock and call appropriate function
    if theta_inf_3 < 0:
        params = [gamma,-theta_inf_3,M_inf,p_inf,T_inf]
        M_3, p_3, T_3 =  external.get_params_across_expansion(params)
    else:
        params = [gamma,theta_inf_3,M_inf,p_inf,T_inf]
        M_3, p_3, T_3 =  external.get_params_across_oblique_shock(params)

    ## get pressures at point 4
    # get theta 
    theta_3_4 = wedgeAngle[3]
    # set parameters and call function
    if theta_inf_3 < 0:
        params = [gamma,-theta_3_4, M_3, p_3, T_3]
        M_4, p_4, T_4 =  external.get_params_across_expansion(params)
    else:
        params = [gamma,theta_3_4, M_3, p_3, T_3]
        M_4, p_4, T_4 =  external.get_params_across_oblique_shock(params)

    return M_1,p_1,T_1, M_2,p_2,T_2, M_3,p_3,T_3, M_4,p_4,T_4


if __name__ == "__main__":
     # get data from standard atmosphere
    wedgeAngle = np.deg2rad(2)
    atm = np.genfromtxt("US_Standard_Atmosphere.txt",skip_header=2)
    atm_alts = atm[:,0]
    atm_press = atm[:,1]
    atm_T = atm[:,2]
    atm_rho = atm[:,3]
    atm_mu = atm[:,4] * 1e-3
    atm_a = atm[:,5]

    print("Problem 2")

    # get data for L/D vs alpha for mach numbers
    numCases = 40
    index = 0
    gamma = atm_a[index]**2*atm_rho[index]/atm_press[index]
    p_inf = atm_press[index]
    T_inf = atm_T[index]
    machs = [1.5,2,4,8]
    alphas = np.deg2rad(np.linspace(0,15,numCases))
    LODS = np.zeros(numCases)
    for M_inf in machs:
        for i in range(numCases):
            params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alphas[i]]
            LODS[i] = getWedgeLD(params,visc=False,TurbulentOnly=False,compressible=True)
        plt.plot(np.rad2deg(alphas),LODS,label=f"Mach {M_inf}")
    plt.xlabel("Alpha(°)")
    plt.ylabel("$C_L/C_D$")
    plt.legend()
    plt.savefig("LODvsAlpha",dpi=300)
    plt.show()

    # get data for L/D max vs mach for altitudes
    numCases = 20
    numCases_alpha = 45
    altitudes = [10,20,30,40,50]#km
    machs = np.linspace(1.75,10,numCases)
    alphas = np.deg2rad(np.linspace(0,9,numCases_alpha))
    LODmaxs = np.zeros(numCases)
    LODS = np.zeros(numCases_alpha)



    for alt in altitudes:
        index = np.where(atm_alts == alt)[0][0]
        gamma = atm_a[index]**2*atm_rho[index]/atm_press[index]
        p_inf = atm_press[index]
        T_inf = atm_T[index]
        for i in range(numCases):
            M_inf = machs[i]
            for j in range(numCases_alpha):
                params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alphas[j]]
                LODS[j] = getWedgeLD(params,visc=False,TurbulentOnly=False,compressible=True)
            LODmaxs[i] = max(LODS)
        plt.plot(machs,LODmaxs,label=f"{alt}km")
    plt.ylim((0,15))
    plt.xlabel("Mach")
    plt.ylabel("$C_L/C_D max$")
    plt.legend()
    plt.savefig("LODvsMach",dpi=300)
    plt.show()


    print("Problem 3")
    # get data for L/D vs alpha for mach numbers
    c = 2
    numCases_alpha = 30
    altitudes = [10,20,30,40,50]#km
    machs = [2,5,10]
    alphas = np.radians(np.linspace(0,15,numCases_alpha))
    LODmaxs = np.zeros(numCases)
    LODS = np.zeros(numCases_alpha)
    for M_inf in machs:
        for alt in altitudes:
            index = np.where(atm_alts == alt)[0][0]
            gamma = atm_a[index]**2*atm_rho[index]/atm_press[index]
            p_inf = atm_press[index]
            T_inf = atm_T[index]
            rho_inf = atm_rho[index]
            a_inf = atm_a[index]
            mu_inf = atm_mu[index]
            LODS = np.zeros(numCases_alpha)
            LODSinv = np.zeros(numCases_alpha)
            for j in range(numCases_alpha):
                params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alphas[j],c,rho_inf,mu_inf,a_inf]
                LODS[j],Re = getWedgeLD(params,visc=True,TurbulentOnly=True,compressible=True)
                LODSinv[j] = getWedgeLD(params,visc=False,TurbulentOnly=False,compressible=True)
            plt.plot(np.degrees(alphas),LODS,label=f"Turbulent {alt}km")
            plt.plot(np.degrees(alphas),LODSinv,label=f"Inviscid {alt}km")
        plt.xlabel("Alpha(°)")
        plt.ylabel("$C_L/C_D$")
        plt.legend()
        saveName = f"LODvsMach_visc_M{M_inf}"
        plt.savefig(saveName.replace(".","_"),dpi=300)
        plt.show()
        

    # L/D max vs mach number
    c = 2
    numCases = 30
    numCases_alpha = 30
    altitudes = [10,20,30,40,50]#km
    machs = np.linspace(1.5,10,numCases)
    alphas = np.radians(np.linspace(0,10,numCases_alpha))
    LODmaxs = np.zeros(numCases)
    LODS = np.zeros(numCases_alpha)

    for alt in altitudes:
        index = np.where(atm_alts == alt)[0][0]
        gamma = atm_a[index]**2*atm_rho[index]/atm_press[index]
        p_inf = atm_press[index]
        T_inf = atm_T[index]
        rho_inf = atm_rho[index]
        a_inf = atm_a[index]
        mu_inf = atm_mu[index]
        for i in range(numCases):
            M_inf = machs[i]
            for j in range(numCases_alpha):
                params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alphas[j],c,rho_inf,mu_inf,a_inf]
                LODS[j],Re = getWedgeLD(params,visc=True,TurbulentOnly=True,compressible=True)
            LODmaxs[i] = max(LODS)
        plt.plot(machs,LODmaxs,label=f"{alt}km")

    plt.xlabel("Mach")
    plt.ylabel("$C_L/C_D$")
    plt.legend()
    plt.savefig("LODvsMach_viscT",dpi=300)
    plt.show()


    #____________________________________________________________________________________________________________________
    print("Problem 4")
    # get data for L/D vs alpha for mach numbers
    # L/D max vs mach number
    c = 2
    numCases = 30
    numCases_alpha = 30
    altitudes = [10,20,30,40,50]#km
    machs = np.linspace(1.5,10,numCases)
    alphas = np.deg2rad(np.linspace(0,10,numCases_alpha))
    LODmaxs = np.zeros(numCases)
    LODS = np.zeros(numCases_alpha)
    Res = np.zeros(numCases)

    

    for alt in altitudes:
        index = np.where(atm_alts == alt)[0][0]
        gamma = atm_a[index]**2*atm_rho[index]/atm_press[index]
        p_inf = atm_press[index]
        T_inf = atm_T[index]
        rho_inf = atm_rho[index]
        a_inf = atm_a[index]
        mu_inf = atm_mu[index]
        for i in range(numCases):
            M_inf = machs[i]
            for j in range(numCases_alpha):
                params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alphas[j],c,rho_inf,mu_inf,a_inf]
                LODS[j], Re = getWedgeLD(params,visc=True,TurbulentOnly=False,compressible=True)
            LODmaxs[i] = max(LODS)
        plt.plot(machs,LODmaxs,label=f"{alt}km")

    plt.xlabel("Mach")
    plt.ylabel("$C_L/C_D$")
    plt.legend()
    plt.savefig("LODvsMach_visc",dpi=300)
    plt.show()



    # part b
    c = 2
    numCases = 40
    numCases_alpha = 40
    altitudes = [10,20,30,40,50]#km
    machs = np.linspace(1.75,10,numCases)
    alphas = np.deg2rad(np.linspace(0,9,numCases_alpha))
    LODmaxs = np.zeros(numCases)
    LODS = np.zeros(numCases_alpha)
    Remaxs = np.zeros(numCases)

    

    for alt in altitudes:
        index = np.where(atm_alts == alt)[0][0]
        gamma = atm_a[index]**2*atm_rho[index]/atm_press[index]
        p_inf = atm_press[index]
        T_inf = atm_T[index]
        print(T_inf)
        rho_inf = atm_rho[index]
        a_inf = atm_a[index]
        mu_inf = atm_mu[index]
        
        for i in range(numCases):
            Res = []
            M_inf = machs[i]
            for j in range(numCases_alpha):
                params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alphas[j],c,rho_inf,mu_inf,a_inf]
                LODS[j], Re = getWedgeLD(params,visc=True,TurbulentOnly=False,compressible=True)
                Res.append(Re)
            LODmaxs[i] = max(LODS)
            index = np.argmax(LODS)
            Remaxs[i] = Res[index]
        plt.plot(Remaxs,LODmaxs,label=f"{alt}km")
    plt.xlabel("Re")
    plt.ylabel("$C_L/C_D$")
    plt.legend()
    plt.xscale("log")
    plt.savefig("LODvsRe_visc",dpi=300)
    plt.show()
























    # c = 2
    # wedgeAngle = np.deg2rad(2)
    # M_inf = 4
    # p_inf = atm_press[1]
    # T_inf = atm_T[1]
    # alpha = np.deg2rad(2)
    # rho_inf = atm_rho[1]
    # a_inf = atm_a[1]
    # mu_inf = atm_mu[1]
    # gamma = a_inf**2*rho_inf/p_inf
    
    

    # params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alpha,c,rho_inf,mu_inf,a_inf]
    # LOD = getWedgeLD(params,visc=True)
    # print(LOD)






    
  # numCases = 45
    # c = 2
    # p_inf = atm_press[0]
    # T_inf = atm_T[0]
    # alpha = np.deg2rad(2)
    # rho_inf = atm_rho[0]
    # a_inf = atm_a[0]
    # mu_inf = atm_mu[0]
    # gamma = a_inf**2*rho_inf/p_inf
    # machs = [2,5,10]
    # alphas = np.deg2rad(np.linspace(0,15,numCases))
    # LODS = np.zeros(numCases)
    # for M_inf in machs:
    #     for i in range(numCases):
    #         params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alphas[i]]
    #         LODS[i] = getWedgeLD(params,visc=False,TurbulentOnly=False,compressible=True)
    #     plt.plot(np.rad2deg(alphas),LODS,label=f"Mach {M_inf} Inviscid")

    # LODS = np.zeros(numCases)
    # for M_inf in machs:
    #     for i in range(numCases):
    #         params = [gamma,wedgeAngle,M_inf,p_inf,T_inf,alphas[i],c,rho_inf,mu_inf,a_inf]
    #         LODS[i],Re = getWedgeLD(params,visc=True,TurbulentOnly=True,compressible=True)
    #     plt.plot(np.rad2deg(alphas),LODS,label=f"Mach {M_inf} Turbulent")


    # plt.xlabel("Alpha(°)")
    # plt.ylabel("$C_L/C_D$")
    # plt.legend()
    # plt.savefig("LODvsAlpha_viscT",dpi=300)
    # plt.show()