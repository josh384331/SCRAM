def newton_solver(f,f_prime,guess,params,tol=1e-8,max_iteration=100,verbose=True):
    # set initial guess
    x = guess

    # run newton solver loop
    for i in range(max_iteration):
        # update law
        delta_x = f(x,params) / f_prime(x,params)
        x -= delta_x
        # end conditions
        if abs(delta_x) < tol:
            if verbose:
                print("Converged to {0} after {1} iterations".format(x,i+1))
            return x

    # if max iterations then print this and retrun x
    print("Did not converge after {0} iterations".format(max_iteration))
    return x

# isentropic pitot equation
def f_iso(M,params):
    gamma = params[0]
    pressure_ratio = params[1]
    return (1 + (gamma-1)/2 * M**2)**((gamma)/(gamma-1))-pressure_ratio

 # isentropic pitot derivative   
def f_prime_iso(M,params):
    gamma = params[0]
    return gamma*M*(1+(gamma-1)/2*M**2)**(1/(gamma-1))

# rayleigh pitot tube equation
def f_rayleigh(M,params):
    gamma = params[0]
    pressure_ratio = params[1]
    return (((gamma+1)/2 * M**2)**(gamma/(gamma-1)))/((2*gamma/(gamma+1)*M**2-(gamma-1)/(gamma+1))**(1/(gamma-1))) - pressure_ratio

# rayleigh pitot tube derivative
def f_prime_rayleigh(M,params):
    gamma = params[0]
    pressure_ratio = params[1]
    return (gamma*M*(2*M**2-1)*(M**2*(gamma+1)/2)**(1/(gamma-1)))/((2*gamma/(gamma+1)*M**2-(gamma-1)/(gamma+1))**(gamma/(gamma-1)))


def get_mach_from_pitot(pressure_ratio,gamma,verbose=True):
    # set parameters
    params = [0,0]
    params[0] = gamma
    params[1] = pressure_ratio
    # check if flow is choked in tube
    test = ((gamma+1)/2)**(gamma/(gamma-1))
    if pressure_ratio >= test:
        # supersonic case
        M = newton_solver(f_rayleigh,f_prime_rayleigh,6,params,verbose=verbose)
    else:
        # subsonic case
        M = newton_solver(f_iso,f_prime_iso,0.1,params,verbose=verbose)
    return M


if __name__ == "__main__":
    # define parameters
    pressure_ratios = [122/101,7222/2116,13107/1020]
    gamma = 1.4
    # find mach number for each case
    for i in range(len(pressure_ratios)):
        print("Pressure Ratio:",pressure_ratios[i])
        M = get_mach_from_pitot(pressure_ratios[i],gamma,verbose=False)
        print("Mach Number: ",M)
