# SCRAM
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
- This code assumes 1D heat conduction and does not account for radiation heat transfer within the material.
- The code is specifically designed for analyzing thermal conduction in supersonic flow scenarios. Subsonic flow will cause an error
- The accuracy of the results depends on the chosen spatial and time steps, as well as the chosen properties of the TPS layers.  Time steps are dynamically calculated to ensure accurate results, however this results in VERY small timesteps. The `output_skip` option allows the user to skip rows to reduce output file size.

## Example
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