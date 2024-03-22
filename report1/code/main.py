import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
from dataclasses import dataclass
import scipy

HBAR_C = 197.3269631 # MeV fm
ALPHA = 1.0/137.035399 # Unit-less and dimension-less
V0 = 134 # MeV
HELIUM = 4 # Mass number
DELTA_R = 1.4 # fm
atomicMassUnit = 931.5 # MeV/c^2
C_SQUARED = (299792458 * 1e-15) ** 2 # fm^2/s^2

mc_squred = 0.511 # MeV2

@dataclass
class Preference:
    numberOfBins: int
    protonNumber: int
    neutronNumber: int

def potentialEnergy(r: float, Z: int, radii) -> float:
    if r <= radii:
        return -V0
    return +ALPHA / (r) * Z * HELIUM * HBAR_C # Some scaling factor for now

def potentialEnergySolver(r: float, Z: int, radii, kineticEnergy: float) -> float:
    iterations = 0
    while abs(ALPHA / (r) * Z * HELIUM * HBAR_C - kineticEnergy) > 0.01 and iterations < 2000:
        r += 0.001
        iterations += 1
    if iterations == 1000:
        raise ValueError("Could not find a solution")
    else:   return r
    


def computeK(energy: float, potential: float) -> float:
    # Return the wave number
    if energy < potential:
        return np.sqrt( 2.0 * (potential - energy) * mc_squred ) / HBAR_C

    return 1j * np.sqrt( 2.0 * (energy - potential) * mc_squred ) / HBAR_C
    

def createMatrix(pref: Preference) -> Tuple[np.ndarray, np.ndarray]:
    """
        We have 2 * numberOfBins + 2 equations
        Therefore
    """
    size = 2 * numberOfBins + 1
    matrix: np.ndarray = np.zeros((size, size), dtype="complex128")

    matrix[0,0] = 1 # We force the first constant to be 1.

    R = 1.4 * (pref.protonNumber + pref.neutronNumber) ** (1/3) # fm

    maxiumDistance = potentialEnergySolver(R, pref.protonNumber, R, 47.3089)

    boundaries = np.linspace(R, maxiumDistance, numberOfBins)

    index = 0

    KList = []

    for row in range(1, size - 1, 2):
        boundary = boundaries[row // 2]

        # The energies will differer, as well as the potential
        K_left = computeK(0.0, potentialEnergy(boundary, pref.protonNumber, R))
        K_right = computeK(0.0, potentialEnergy(boundary, pref.protonNumber, R))

        if row == 1:
            KList.append(K_left)

        if row > size - 3:
            matrix[row][index] = np.exp(K_left * boundary)
            matrix[row][index + 1] = np.exp(-K_left * boundary)
            matrix[row][index + 2] = np.exp(K_right * boundary)

            matrix[row + 1][index] = K_left * np.exp(K_left * boundary)
            matrix[row + 1][index + 1] = -K_left * np.exp(-K_left * boundary)
            matrix[row + 1][index + 2] = K_right * np.exp(K_right * boundary)
            KList.append(K_right)

            break
 
        matrix[row][index] = np.exp(K_left * boundary)
        matrix[row][index + 1] = np.exp(-K_left * boundary)
        matrix[row][index + 2] = np.exp(K_right * boundary)
        matrix[row][index + 3] = np.exp(-K_right * boundary)

        matrix[row + 1][index] = K_left * np.exp(K_left * boundary)
        matrix[row + 1][index + 1] = -K_left * np.exp(-K_left * boundary)
        matrix[row + 1][index + 2] = K_right * np.exp(K_right * boundary)
        matrix[row + 1][index + 3] = -K_right * np.exp(-K_right * boundary)

        index += 2
        
    solution = np.zeros((size,  1))    
    solution[0] = 1.0
    
    for i in range(len(matrix)):
        print(matrix[i])

    return matrix, solution





if __name__ == "__main__":
    numberOfBins = 2
    pref: Preference = Preference(numberOfBins, 2, 2)

    matrix, sol = createMatrix(pref)

    #print(matrix)

    #constants = np.real(matrix) @ np.real(sol) + 1.j * np.imag(matrix) @ np.imag(sol)
    #print(constants)

    #constants = matrix @ sol

    #constants = np.linalg.solve(np.real(matrix), np.real(sol)) + 1.j*np.linalg.solve(np.imag(matrix), np.imag(sol))

    #print(constants)
