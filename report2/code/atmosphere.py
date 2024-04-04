import numpy as np
from dataclasses import dataclass
import os
import sys

sys.setrecursionlimit(10000)

# Constants
R = 8.3144598
STEFANBOLTZMAN = 5.67 * 10**(-8)


@dataclass
class Atmosphere:
    """
        A dataclass representing the atmosphere of a planet

        params:
            incomingFlux: float - the incoming flux in W/m^2
            height: float - the height of the atmosphere in meters
            surfacePressure: float - the surface pressure in pascals
            gravity: float - the gravity at the surface in m/s^2
            groundTemperature: float - the temperature at the surface in Kelvin
            deltaHeight: float - the height of each cell in the atmosphere
            numberOfCells: int - the number of cells in the atmosphere
            surfaceDensity: float - the density at the surface in kg/m^3
            mAir: float - the molar mass of air in kg/mol
            planetAlbedo: float - the albedo of the planet

        returns:
            Atmosphere - the atmosphere of the planet
    """
    incomingFlux: float
    height: float
    surfacePressure: float
    gravity: float
    groundTemperature: float
    deltaHeight: float
    numberOfCells: int
    surfaceDensity: float
    mAir: float = 0.0289644
    planetAlbedo: float = 0.04


def computePressure(h: float) -> float:
    """
        Computes the pressure given the height
        
        params:
            h: float - height in meters
        returns:
            float - pressure in pascals
    """
    return atm.surfacePressure * np.exp(-atm.gravity * h * atm.mAir/ (R * atm.groundTemperature))

def computeDensity(atm: Atmosphere, h: float) -> float:
    """
        Computes the density given the height
        
        params:
            atm: Atmosphere - the atmosphere of the planet
            h: float - height in meters
        returns:
            float - Density in kg/m^3
    """
    return atm.surfaceDensity * np.exp( -atm.gravity * atm.mAir * (h + atm.deltaHeight) / (R * atm.groundTemperature))




def iterate(iteration: int, absorbed, transmittedDownInf, transmittedDownVis, transmittedUpVis, transmittedUpInf, emittedUp, emittedDown, cellFlux: np.ndarray) -> np.ndarray:
    """
        Iterates through the cells to calculate the energy balance, via recursion

        params:
            iteration: int - the current iteration
            absorbed: np.ndarray - the absorbed energy in each cell at the current iteration
            transmittedDownInf: np.ndarray - the transmitted IR radiation going down in each cell at the current iteration
            transmittedDownVis: np.ndarray - the transmitted visible radiation going down in each cell at the current iteration
            transmittedUpVis: np.ndarray - the transmitted visible radiation going up in each cell at the current iteration
            transmittedUpInf: np.ndarray - the transmitted IR radiation going up in each cell at the current iteration
            emittedUp: np.ndarray - the emitted energy in each cell going up at the current iteration
            emittedDown: np.ndarray - the emitted energy in each cell going down at the current iteration
            cellFlux: np.ndarray - the cell flux at the current iteration, used for radiation balance

        returns:
            np.ndarray - the cell flux
        
        Error:
            Exception - if the solution does not converge within 5000 iterations
    """

    # Top Down
    for i in range(len(absorbed) - 2, 0, -1):
        #Visible contribution
        transmittedDownVis[ i ] += transmittedDownVis[i + 1] * epsilonV[ i ]

        #IR contribution
        transmittedDownInf[ i ] += (transmittedDownInf[i + 1] + emittedDown[i + 1]) * epsilonI[ i ] 

        # Attenuated in cell
        absorbed[ i ] += transmittedDownVis[i + 1] * (1 - epsilonV[ i ]) + (transmittedDownInf[i + 1] + emittedDown[i + 1]) * ( 1 - epsilonI[ i ] )

        # Emitting IR energy 50% goes up and 50% goes down
        emittedDown[ i ] += absorbed[ i ] / 2
        emittedUp[ i ] += absorbed[ i ] / 2

        # Keeping track of the cells
        cellFlux[ i ] += emittedUp[ i ] + transmittedUpInf[i] + transmittedUpVis[i] - (emittedDown[i + 1] + transmittedDownInf[ i + 1] + transmittedDownVis[i + 1])  # new

        # Reset, we have 'used' the energies in this iteration
        transmittedDownVis[i + 1] = 0.0
        transmittedDownInf[i + 1] = 0.0
        emittedDown[i + 1] = 0.0
        absorbed[ i ] = 0.0
    
    # Ground
    
    # Visible contribution from reflection
    transmittedUpVis[ 0 ] = transmittedDownVis[ 1 ] * atm.planetAlbedo

    # Absorbed energy
    absorbed[ 0 ] = transmittedDownVis[ 1 ] - transmittedUpVis[ 0 ] + emittedDown[ 1 ] + transmittedDownInf[ 1 ]
    
    # Kirchhoff's law
    emittedUp[ 0 ] = absorbed[ 0 ]

    # fSurface += cells[0].absorbed
    cellFlux[ 0 ] += absorbed[ 0 ]
    
    # Reset
    absorbed[ 0 ] = 0.0
    transmittedDownVis[ 1 ] = 0.0
    transmittedDownInf[ 1 ] = 0.0
    emittedDown[ 1 ] = 0.0

    # Bottom up
    for i in range(1, len(absorbed) - 1):
        # Visible contribution
        transmittedUpVis[ i ] += transmittedUpVis[i - 1] * epsilonV[ i ]

        # IR contribution
        transmittedUpInf[ i ] += (transmittedUpInf[i - 1] + emittedUp[i - 1]) * epsilonI[ i ]

        # Absorbed energy
        absorbed[ i ] += transmittedUpVis[i - 1] * (1 - epsilonV[ i ]) + (transmittedUpInf[i - 1] + emittedUp[i - 1]) * (1 - epsilonI[ i ])

        # Emitting IR energy 50% goes up and 50% goes down
        emittedDown[ i ] += absorbed[ i ] / 2
        emittedUp[ i ] += absorbed[ i ] / 2

        # Keeping track of the cells
        cellFlux[ i ] += emittedUp[ i ] + transmittedUpVis[ i ] + transmittedUpInf[ i ] - (emittedDown[i + 1] + transmittedDownInf[ i ] + transmittedDownVis[ i ])  # new
        
        # Reset, we have 'used' the energies in this iteration
        transmittedUpVis[i - 1] = 0.0
        transmittedUpInf[i - 1] = 0.0
        emittedUp[i - 1] = 0.0
        absorbed[ i ] = 0.0


    # Radiation into space
    absorbed[ -1 ] = transmittedUpVis[ -2 ] + transmittedUpInf[ -2 ] + emittedUp[ -2 ]

    # Reset, we have 'used' the energies in this iteration
    transmittedUpVis[ -2 ] = 0.0
    transmittedUpInf[ -2 ] = 0.0
    emittedUp[ -2 ] = 0.0

    cellFlux[ -1 ] += absorbed[ -1 ] # Emitted energy into space

    if iteration > 5000:
        raise Exception("Did not converge within 5000 iterations")

    if abs(atm.incomingFlux * 0.7 - cellFlux[-1])> 0.001: # Radiation balance condition
        return iterate(iteration + 1, absorbed, transmittedDownInf, transmittedDownVis, transmittedUpVis, transmittedUpInf, emittedUp, emittedDown, cellFlux)
    else:
        return cellFlux
 

if __name__ == '__main__':
    # Global variables
    global atm
    global epsilonI
    global epsilonV

    # Initialize the atmosphere
    numberOfCells: int = 50

    planet = input("Enter the planet name: ")

    if planet == "earth":
        height: int = 100_000 # m

        atm = Atmosphere(
            incomingFlux = 340, #  W/m^2
            height = height, # m
            surfacePressure = 101_325, # Pa
            gravity = 9.81, # m/s^2
            groundTemperature = 288, # K
            deltaHeight = height / (numberOfCells), # m
            numberOfCells = numberOfCells,
            surfaceDensity = 1.225   # kg/m^3
        )

        # Initialize the data that we need to keep track of
        absorbed = np.zeros(numberOfCells + 1, dtype = float)     # Used to keep track of the absorbed energy in each cell, for radiation balance and temperature calculation
        transmittedDownVis = np.zeros(numberOfCells + 1, dtype = float)   # Used to keep track of the transmitted visible radition going down
        transmittedUpVis = np.zeros(numberOfCells + 1, dtype = float)     # Used to keep track of the transmitted visible radiation going up
        transmittedDownInf = np.zeros(numberOfCells + 1, dtype = float)   # Used to keep track of the transmitted ir radiation going down
        transmittedUpInf = np.zeros(numberOfCells + 1, dtype = float)     # Used to keep track of the transmitted ir radiation going up
        emittedUp = np.zeros(numberOfCells + 1, dtype = float)    # Used to keep track of the emitted energy in each cell going up
        emittedDown = np.zeros(numberOfCells + 1, dtype = float)  # Used to keep track of the emitted energy in each cell going down
        
        # Compute the attenuation coefficients for each layer 
        densities = np.zeros(numberOfCells + 1, dtype = float)
        densities[0] = atm.surfaceDensity
        pressures = np.zeros(numberOfCells + 1, dtype = float)
        pressures[0] = atm.surfacePressure

        #fileLayer = open('layerInfo.dat', 'w') # Uncomment the file writing for densities and pressures plot
        #fileLayer.write('z\trho\tp\n')
        #fileLayer.write('0\t{}\t{}\n'.format(densities[0]/densities[0], pressures[0]/pressures[0]))
        Z: np.ndarray = np.zeros(numberOfCells + 1, dtype = float)
        for i in range(0, numberOfCells):
            z = i * atm.deltaHeight
            Z[i + 1] = z
            densities[i + 1] = computeDensity(atm, z)
            pressures[i + 1] = computePressure(z + atm.deltaHeight / 2)
            #fileLayer.write('{}\t{}\t{}\n'.format(z, densities[i+1]/densities[0], pressures[i+1]/pressures[0]))

        #fileLayer.close()

        # Attenuation for visible and IR radiation
        attV: float = 1e-4 / densities[0]
        attInf: float = 1.07e-3 / densities[0]

        # Attenuation for visible and IR radiation
        epsilonV: np.ndarray = np.exp( -attV  * densities * atm.deltaHeight)
        epsilonI: np.ndarray = np.exp( -attInf  * densities * atm.deltaHeight) 

        # Initizalize the incoming radiation
        inc: float = atm.incomingFlux * 0.7     # 30% of the incoming radiation is reflected
        transmittedDownVis[ -1 ] = inc * 0.25   # 25 % of the incoming is visible
        transmittedDownInf[ -1 ] = inc * 0.75   # 75 % of the incoming is IR

        cellFlux = np.zeros(numberOfCells + 1, dtype=float) # Used to keep track of the emitted energy in each cell, for radiation balance and temperature calculation

    elif planet == "venus":
        height: int = 250_000 # m

        atm = Atmosphere(
            incomingFlux = 2622/4, #  W/m^2
            height = height, # m
            surfacePressure = 9.3 * 10*6, # Pa
            gravity = 8.87, # m/s^2
            groundTemperature = 737, # K
            deltaHeight = height / (numberOfCells), # m
            numberOfCells = numberOfCells,
            surfaceDensity = 67 # kg/m^3
        )

        # Initialize the data that we need to keep track of
        absorbed = np.zeros(numberOfCells + 1, dtype = float)     # Used to keep track of the absorbed energy in each cell, for radiation balance and temperature calculation
        transmittedDownVis = np.zeros(numberOfCells + 1, dtype = float)   # Used to keep track of the transmitted visible radition going down
        transmittedUpVis = np.zeros(numberOfCells + 1, dtype = float)     # Used to keep track of the transmitted visible radiation going up
        transmittedDownInf = np.zeros(numberOfCells + 1, dtype = float)   # Used to keep track of the transmitted ir radiation going down
        transmittedUpInf = np.zeros(numberOfCells + 1, dtype = float)     # Used to keep track of the transmitted ir radiation going up
        emittedUp = np.zeros(numberOfCells + 1, dtype = float)    # Used to keep track of the emitted energy in each cell going up
        emittedDown = np.zeros(numberOfCells + 1, dtype = float)  # Used to keep track of the emitted energy in each cell going down
        
        # Compute the attenuation coefficients for each layer 
        densities = np.zeros(numberOfCells + 1, dtype = float)
        densities[0] = atm.surfaceDensity
        pressures = np.zeros(numberOfCells + 1, dtype = float)
        pressures[0] = atm.surfacePressure

        #fileLayer = open('layerInfo.dat', 'w') # Uncomment the file writing for densities and pressures plot
        #fileLayer.write('z\trho\tp\n')
        #fileLayer.write('0\t{}\t{}\n'.format(densities[0]/densities[0], pressures[0]/pressures[0]))
        Z: np.ndarray = np.zeros(numberOfCells + 1, dtype = float)
        for i in range(0, numberOfCells):
            z = i * atm.deltaHeight
            Z[i + 1] = z
            densities[i + 1] = computeDensity(atm, z)
            pressures[i + 1] = computePressure(z + atm.deltaHeight / 2)
            #fileLayer.write('{}\t{}\t{}\n'.format(z, densities[i+1]/densities[0], pressures[i+1]/pressures[0]))

        #fileLayer.close()

        # Attenuation for visible and IR radiation
        attV: float = 1e-6 / densities[0]
        attInf: float = 5e-1 / densities[0]

        # Attenuation for visible and IR radiation
        epsilonV: np.ndarray = np.exp( -attV  * densities * atm.deltaHeight)
        epsilonI: np.ndarray = np.exp( -attInf  * densities * atm.deltaHeight) 
        print(densities)

        # Initizalize the incoming radiation
        inc: float = atm.incomingFlux * 0.7     # 30% of the incoming radiation is reflected
        transmittedDownVis[ -1 ] = inc * 0.25   # 25 % of the incoming is visible
        transmittedDownInf[ -1 ] = inc * 0.75   # 75 % of the incoming is IR

        cellFlux = np.zeros(numberOfCells + 1, dtype=float) # Used to keep track of the emitted energy in each cell, for radiation balance and temperature calculation
    else:
        print("Invalid planet name")
        sys.exit(1)

    
    # Iterate until the solution converges or the maximum number of iterations is reached
    flux = iterate(0, absorbed, transmittedDownInf, transmittedDownVis, transmittedUpVis, transmittedUpInf, emittedUp, emittedDown, cellFlux)

    # Convert the flux to temperature via the Stefan-Boltzman law
    temperature = ( flux / STEFANBOLTZMAN ) ** 0.25

    file = open("testoutput{}.dat".format(planet), "w") # For the temperature at different heights
    for i in range( len( flux ) - 1 ):
        current_height = i * atm.deltaHeight / 1000 # Convert to km
        file.write(f"{current_height} {temperature[i]}\n")

    file.close()

    #file = open("ground_temperature2.dat", "a+")
    #file.write('{}\t{}\n'.format(numberOfCells, temperature[0]))
    #file.close()

    print("-"*50)
    print("Temperatures [K]")
    print(temperature)

    