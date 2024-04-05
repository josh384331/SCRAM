import numpy as np
#  This function was written by Joshua Hurwitz and is used to calculate the atmospheric properties at a given altitude using the 1976 Standard atmosphere model
def get_atmospheric_properties_si(altitude):
    """
    Calculates the atmospheric properties at a given altitude in SI units from the 1976 US Standard Atmosphere model.

    Parameters:
    altitude (float): The altitude in meters.

    Returns:
    list: A list containing the following atmospheric properties:
        - Geopotential altitude (Z) in meters
        - Temperature (T) in Kelvin
        - Pressure (p) in Pascal
        - Density (rho) in kg/m^3
        - Speed of sound (a) in m/s
    """
    # variables
    i = 0
    check = False
    ans = [0.0] * 5
    Z = 0.0
    T = 0.0
    temp = 0.0
    rho = 0.0
    a = 0.0
    p = 101325.0
    R_earth = 6356766.0
    g0 = 9.806645
    R = 287.0528
    gamma = 1.4
    Zi = [0, 11000, 20000, 32000, 47000, 52000, 61000, 79000, 90000]
    Ti = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 252.65, 180.65]
    dTi = [-0.0065, 0, 0.001, 0.0028, 0, -0.002, -0.004, 0]
    pi = [0.0] * 8

    # find geopotential altitude
    Z = R_earth * altitude / (R_earth + altitude)

    # find temperature
    for i in range(9):
        if Z < Zi[i + 1]:
            T = Ti[i] + dTi[i] * (Z - Zi[i])
            break

    # find pressure
    j = 0
    while not check:
        if Z < Zi[j + 1]:
            temp = Z
            check = True
        else:
            temp = Zi[j + 1]
        if dTi[j] == 0:
            p *= np.exp(-g0 * (temp - Zi[j]) / (R * Ti[j]))
        else:
            p *= pow(((Ti[j] + dTi[j] * (temp - Zi[j])) / Ti[j]), (-g0 / (R * dTi[j])))
        j += 1

    # find density
    rho = p / R / T

    # find speed of sound
    a = (gamma * R * T)**0.5

    # save answer
    ans[0] = Z
    ans[1] = T
    ans[2] = p
    ans[3] = rho
    ans[4] = a
    return ans
