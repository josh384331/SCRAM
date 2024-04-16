# Skin Conduction & Radiative Analysis Model (SCRAM)
This code solves the 1D heat conduction equation numerically using the finite difference method. It calculates the temperature distribution along a surface subjected to supersonic flow and convective heat transfer.

## Usage
To use this code, follow these steps:

1. Install the required dependencies: `numpy`, `time`, `json`, `csv`, `matplotlib`, and `scipy`.
2. Create a JSON input file specifying the desired parameters for the simulation. See the example `SCRAM_Input.json` file for reference.
3. Run the code by executing the `runTrajectory(inputFileName)` function, passing the name of the input file as an argument.

## Input Parameters
The input parameters for the simulation are specified in a JSON file. Here are the key parameters:

- `options`: Contains general options for the simulation, such as the initial temperature, number of spatial and time steps, and verbosity.
- `TPS`: Contains the properties of the thermal protection system (TPS) layers, including thermal conductivity, thickness, density, and specific heat capacity.
- `geometry`: Contains the geometry properties of the surface, such as wedge angles, lengths, and the location of the point of interest.
- `trajectory`: Contains the trajectory data, including freestream velocity, time, angle of attack, and altitude.

## Output
The code outputs two text files:

- `time_list.txt`: Contains the time values at which the temperature distribution was calculated.
- `temperature_distribution.txt`: Contains the temperature distribution along the surface at different time steps.


## Limitations
- This code assumes 1D heat conduction and does not account for radiation heat transfer within the material. (but does account for it leaving the material at the boundaries)
- The code is specifically designed for analyzing thermal conduction in supersonic flow scenarios. Subsonic flow will cause an error.
- The accuracy of the results depends on the chosen spatial and time steps, as well as the chosen properties of the TPS layers.  Time steps are dynamically calculated to ensure accurate results, however this results in VERY small timesteps. The `output_skip` option allows the user to skip rows to reduce output file size.
- This code utilizes a resistance network to speed up and increase stability of the method.  This means that the code creates an equivilent isotropic material and determines the equivalent temperature distribution through that material.  This yields accurate inner and outer surface temperature but would likely not yield good material temperatures.  

## Example Input File
An example input file (`SCRAM_Input.json`) is provided in the repository. You can use this file as a template to create your own input file.

```json
{
    "TPS": [
        {
            "thickness [m]": 0.001127,
            "cp [J/kgK]": 585.76,
            "rho [kg/m^3]": 4540,
            "k [W/mK]": 9.28
        },
        {
            "thickness [m]": 0.0458938,
            "cp [J/kgK]": 900,
            "rho [kg/m^3]": 38,
            "k [W/mK]": 0.05328
        }
    ],
    "geometry": {
        "lengths [m]": [13,13,13,13],
        "angles [deg]" : [6,-6,0,0],
        "location [m]": 3,
        "top_edge" : true
    },
    "options": {
        "initial_temperature [K]": 300,
        "epsilon": 0.75,
        "spatial_steps": 100,
        "time_steps": 10000000,
        "verbose": true,
        "T_infCal [K]": -10,
        "output_file_prefix": "output_",
        "output_skip": 250,
        "epsilonBack": 0.15,
        "hBack [W/m^2K]": 0,
        "backTemp [K]" : 300
    }  
}
```
### TPS
This section is a list of dictionaries where each dictionary contains the following required values 
-	`thickness [m]` : the thickness of this layer of TPS
-	`cp [J/kgK]` : the specific heat capacity of the material
-	`rho [kg/m^3]`: the density of the material 
-	`k [W/mK]`: the thermal conductivity of the material
### geometry
This section contains geometric parameters to determine the heat transfer into the vehicle
-	`lengths [m]`: A list of 4 lengths of each panel of the vehicle
-	 `angles [deg]` : A list of 4 angles between panels
-	`location [m]`: the length down the panels to take temperature data at 0 corresponds to the leading edge (not computed)
-	`top_edge` : a bool value stating whether the point is on the top of the vehicle or the bottom
### options
This section contains the solver options and additional parameters
-	`initial_temperature [K]`: The initial temperature of the material
-	`epsilon`: The emissivity of the material surface
-	`spatial_steps`: The number of spatial steps through the thickness
-	`time_steps`: The maximum number of time steps to take
-	`verbose`: a bool stating whether the program should print information
-	`T_infCal [K]`: an air temperature calibration term.  This is added to the 1976 US standard atmosphere temperature
-	`output_file_prefix`: The prefix for the output files
-	`output_skip`: The number of timesteps to skip when printing as timesteps are VERY small
-	`epsilonBack`: OPTIONAL.  The emissivity of the back of the TPS stackup.  Defaults to adiabatic
-	`hBack [W/m^2K]`: The convective heat transfer coefficient of the back of the TPS stackup. Defaults to adiabatic
-	`backTemp [K]` : The internal temperature of the vehicle used if 	`epsilonBack` and/or `hBack` are present

### trajectory
The trajectory used in this code is specified in the `trajectory.csv` also provided in the repository.  The columns of this .csv file are
- `Velocity[m/s]` : The velocity of the vehicle at each control point
- `Time[s]` : The time at which the vehicle gets to each control point
- `AOA[deg]`: The angle of attack the vehicle is flying at at each control point
- `Altitude[m]`: The altitude of the vehicle at each control point

## Validation (read: "Gut Check")
The example trajectory and input file creates a case that matches NASA TN D-5449 Summary of XB-70 Airplane Cockput Enviromental Data by Kirk S. Irwin and William H. Andrews.  The material properties are assumed properties from current materials matching the description of XB-70 materials.  The calibration temperature was adjusted until the freestream temperature matched the experimental temperatures shown in the report and an assumed interior BC was chosen.  A screenshot of the results overlaied onto the NASA chart is included in the repository. 
