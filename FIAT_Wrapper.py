from Thermal_Analysis import *
from standardAtmosphere import *
import numpy as np

# Main input file for FIAT
# set 1 Run control parameters
title = "test_title\n"
IDCAL = 0 # mode, 0 neither pyrolysis nor surface recession, 1 pyrolysis without surface recession,2 pyrolysis with surface recession(ablative)
REST = 400. #R back face temperature
RESH = 0. #BTU/ft-s-R back face heat transfer coefficient (0 for adiabatic)
RADTMP = 400. #R T inf for radiation
IUNIT = 0 # units 0 metric, 1 English
IMATL = 0 # material file, 0 use input file matdatabase.inp, !=0 use material file matprp.inp
ICONV = 0 # convergence criteria
L2 = [IDCAL, REST, RESH, RADTMP, IUNIT, IMATL, ICONV]

# set 2 thermal output parameters
NTCS = 3 # number of temperature output nodes
NISO = 2 # number of isotherm outputs
if NTCS + NISO > 11:
    print("Error: NTCS + NISO > 11")
    exit()

DTC = [0,0.5,1]
TISO = [500, 1000]
# thickness optimization terms
IDPLY = 0 # 0 no ply optimization, >0 material ply to be optimized
NOITER = 0 # number of iterations, 0 no iteration
THKMIN = 0
THKMAX = 0
# if IDPLY and NOITER are not 0, then add some other parameters

# set 3 run and output control parameters (if different print timesteps are needed adjust this)
DTPRNT = 60 # print interval in s

# set 4 geometry parameters 
RAFIX = 0 # these are for planar geometry.  read the manual for other geometries
IAEXP = 0
JAFIX = 0


# set 5 material properties
MATID = ["TI-6al-4V", "Tihc1.5_605b-nom"]
TINIT = [540., 540.] # initial temperature in R
THKPLY = [0.032, 3.] # thickness in inches
CONRST = [0.0, 0.0] # contact resistance in ft-s-R/BTU

# Surface Enviroments
# set 1
NUMP = -1
KEYWORD = "Earth   " #leave the extra spaces (8 characters)

# set 2
NTAB2 = -1
HFACT = 1
RFACT = 1
SFACT = 1
FFACT = 1
OFACT = 1
PFACT = 1


# trajectory Definitions
timeList = [0, 240,480,720,960, 1200, 1440, 1680,5680] # s
MachList = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5,5] # []
AltitudeList = [4114.80, 6400, 8686,10972,13258,15544,17830,20116,29000] # m
AOAList = [4, 4, 4, 4, 4, 4, 4, 4, 4] # degrees 
gamma = 1.4 # ratio of specific heats
wedgeAngles = np.deg2rad([6, -6, 0, 0]) # degrees
length_list = [13, 13, 13, 13] # m
location = 1
topEdge = False
R_air = 287.05
cpe = gamma * R_air / (gamma - 1) # J/kgK





# collect all the inputs and save them in the main input file
inputFileMain = title
# set 1
inputFileMain += "{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n".format(float(IDCAL), float(REST), float(RESH), float(RADTMP), float(IUNIT), float(IMATL), float(ICONV))

# set 2
inputFileMain += str(NTCS) + "\t" + str(NISO) + "\n"
for i in DTC:
    inputFileMain += "{:.5f}\t".format(float(i))
for i in TISO:
    inputFileMain += "{:.5f}\t".format(float(i))
inputFileMain += "\n"
inputFileMain += str(IDPLY) + "\t" + str(NOITER) + "\t"
if NOITER != 0:
    inputFileMain += str(THKMIN) + "\t" + str(THKMAX) + "\n"
else:
    inputFileMain += "\n"
if IDPLY != 0 and NOITER != 0:
    # add some other parameters
    pass
else:
    pass

# set 3
inputFileMain += str(DTPRNT) + "\n" # if different print timesteps are needed adjust this

# set 4
inputFileMain += str(RAFIX) + "\t" + str(IAEXP) + "\t" + str(JAFIX) + "\n"

# set 5
for i in range(len(MATID)):
    if i == len(MATID) - 1:
        k = -(i+1)
    else:
        k= i+1
    inputFileMain += str(k) + "\t" + MATID[i] + "\t" + str(TINIT[i]) + "\t" + str(THKPLY[i]) + "\t" + str(CONRST[i]) + "\n"


# write the input file
with open("main.inp", "w") as f:
    f.write(inputFileMain)


# Surface Enviroments
surfaceEnvironments = ""
# set 1
surfaceEnvironments += str(NUMP) + "\t" + KEYWORD + "\n"
# set 2
surfaceEnvironments += title + " Surface Environments\n"
surfaceEnvironments += str(NTAB2) + "\t" + str(HFACT) + "\t" + str(RFACT) + "\t" + str(SFACT) + "\t" + str(FFACT) + "\t" + str(OFACT) + "\t" + str(PFACT) + "\n"

for i in range(len(timeList)):
    # get properties from wedge wing code
    altitude = AltitudeList[i]
    M_inf = MachList[i]
    AOA = np.deg2rad(AOAList[i])
    [Z,T_inf,p_inf,rho_inf,a_inf] = get_atmospheric_properties_si(altitude)

    # get edge parameters
    M_e, p_e, T_e, rho_e, v_e = get_edge_params(M_inf, gamma, wedgeAngles, p_inf, T_inf, AOA, location, length_list, topEdge)

    # get h
    RadEquilTemp = getRadiationEquilibrium(0.7, location, topEdge,length_list,wedgeAngles,[M_inf, AOA, gamma, altitude])
    h, T_aw, C_H = getConvectiveHeatTransferCoefficient(RadEquilTemp, M_e, rho_e, v_e, T_e, gamma, T_inf, M_inf,location)
    C_H = h / (rho_e * cpe * v_e)
    if i == len(timeList) - 1:
        ISEBT = -3
    else:
        ISEBT = 3
    TAB2T = timeList[i]
    TAB2TR = 0
    TAB2RD = 0
    TAB2CT = C_H[0] # convective heat transfer coefficient, BTU/ft^2-s-R
    TAB2PW = p_e/101325.0 # pressure, atm
    TAB2BW = T_aw # adiabatic wall fluid temperature, R

    surfaceEnvironments += "{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n".format(float(ISEBT), float(TAB2T), float(TAB2TR), float(TAB2RD), float(TAB2CT), float(TAB2PW), float(TAB2BW))

# write the input file
with open("envir.inp", "w") as f:
    f.write(surfaceEnvironments)





    