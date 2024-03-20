import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple

RADII = 1.0
HBAR_C = 197.3269631 # MeV fm
ALPHA = 1.0/137.035399 # Unit-less and dimension-less
V0 = -14.1 # MeV

mc_squred = 0.511 # MeV

def potentialEnergy(r: np.ndarray) -> float:
    if r <= RADII:
        return -V0
    return 1e-2 / (r) # Some scaling factor for now


def computeK(energy: float, potential: float) -> float:
    # Return the wave number
    if energy < potential:
        return np.sqrt( 2.0 * (potential - energy) * mc_squred ) / HBAR_C

    return 1j * np.sqrt( 2.0 * (energy - potential) * mc_squred ) / HBAR_C
    

def createMatrix(numberOfBins: int) -> Tuple[np.ndarray, np.ndarray]:
    """
        We have 2 * numberOfBins + 2 equations
        Therefore
    """
    size = 2 * numberOfBins + 1
    matrix: np.ndarray = np.zeros((size, size), dtype="complex128")

    matrix[0,0] = 1 # We force the first constant to be 1.

    boundaries: np.ndarray = np.linspace(RADII, 3 * RADII, numberOfBins)

    index = 0

    for row in range(1, size - 1, 2):
        boundary = boundaries[row // 2]

        # The energies will differer, as well as the potential
        K_left = computeK(0.0, potentialEnergy(boundary))
        K_right = computeK(0.0, potentialEnergy(boundary))

        if row > size - 3:
            matrix[row][index] = np.exp(K_left * boundary)
            matrix[row][index + 1] = np.exp(-K_left * boundary)
            matrix[row][index + 2] = np.exp(K_right * boundary)

            matrix[row + 1][index] = K_left * np.exp(K_left * boundary)
            matrix[row + 1][index + 1] = -K_left * np.exp(-K_left * boundary)
            matrix[row + 1][index + 2] = K_right * np.exp(K_right * boundary)

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
    
    return matrix, solution



if __name__ == "__main__":
    numberOfBins = 4

    matrix, sol = createMatrix(numberOfBins)

    #print(matrix)


    constants = np.real(matrix) @ np.real(sol) + 1.j * np.imag(matrix) @ np.imag(sol)
    print(constants)

    #constants = matrix @ sol

    #constants = np.linalg.solve(np.real(matrix), np.real(sol)) + 1.j*np.linalg.solve(np.imag(matrix), np.imag(sol))

    #print(constants)
