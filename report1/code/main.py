import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
from dataclasses import dataclass
import scipy
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

HBAR_C = 197.3269631 # MeV fm
ALPHA = 1.0/137.035399 # Unit-less and dimension-less
V0 = 134 # MeV
HELIUM = 2 # Mass number
DELTA_R = 1.4 # fm
atomicMassUnit = 931.5 # MeV/c^2
C_SQUARED = (299792458 * 10 ** (-15)) ** 2 # fm^2/s^2
MASSOFHELIUM = 3727.379 # MeV
amuToMev = 931.5 # MeV/c^2


@dataclass
class Preference:
    numberOfBins: int
    protonNumber: int
    neutronNumber: int
    massParent: float
    massDaughter: float

def findBoundary(kineticEnergy: float, Z) -> float:
    """
        Finds the boundary where the potential energy is equal to the kinetic energy,
        i.e. where the particle exists the forbidden region.
        This is then the the region from the atom radii to this point is where 
        we discretize the potential.
    """
    return ALPHA / kineticEnergy * HBAR_C * HELIUM * Z
    

def potentialEnergy(r: float, Z: int, radii: float) -> float:
    """
        Computes the potential energy at a given distance r from the nuclei
    """
    if r < radii:
        return -V0
    return ALPHA / (r) * Z * HELIUM * HBAR_C # Some scaling factor for now


def computeK(energy: float, potential: float) -> float:
    """
        Returns the wave-number for a given energy and potential
        the imaginary part is baked into the wave-number
    """
    # Return the wave number
    if energy < potential:
        return np.sqrt( 2.0 * (potential - energy) * MASSOFHELIUM) / HBAR_C

    return 1j * np.sqrt( 2.0 * (energy - potential) * MASSOFHELIUM) / HBAR_C

def computeAlphaEnergy(pref: Preference) -> float:
    """
        Computes the alpha energy from the masses of the parent and daughter nuclei
    """
    return (pref.massParent * amuToMev - (pref.massDaughter * amuToMev + MASSOFHELIUM))
    

def createMatrix(pref: Preference) -> Tuple[np.ndarray, np.ndarray]:
    """
        Computes the matrix equation for the given preference,
        i.e. the number of bins and the number of protons and neutrons
    """
    size = 2*( pref.numberOfBins - 1 ) + 1
    print(size)
    m: np.ndarray = np.zeros((size, size), dtype="complex128")

    m[0,0] = 1 # We force the first constant to be 1.

    kineticEnergy: float = computeAlphaEnergy(pref)

    minimumDistance = 1.4 * (pref.protonNumber + pref.neutronNumber) ** (1/3) # fm

    maxiumDistance = findBoundary(kineticEnergy, pref.protonNumber) # fm

    print('Maxium distance: ', maxiumDistance)

    deltaR = (maxiumDistance - minimumDistance) / (pref.numberOfBins - 2)

    offset = 0.1
    boundaries = [0] * (2 * pref.numberOfBins - 2)
    for i in range(pref.numberOfBins - 1):
        b = minimumDistance + i * deltaR
        boundaries[2 * i] = b - offset
        boundaries[2 * i + 1] = b + offset

    kList = [0, 0]

    for idx in range(0, len(boundaries), 2):
        boundary = boundaries[idx]
        print(idx)
        print('Boundary: ', boundary)
        
        kLeft = computeK(kineticEnergy, potentialEnergy(boundaries[idx], pref.protonNumber, minimumDistance))
        kRight = computeK(kineticEnergy, potentialEnergy(boundaries[idx + 1], pref.protonNumber, minimumDistance))

        print('kLeft: ', kLeft)
        print('kLeft: ', kLeft)
        
        # \phi_i A_i
        m[idx + 1, idx] = np.exp(kLeft * boundary + offset)
        m[idx + 2, idx] = kLeft * np.exp(kLeft * boundary + offset)

        # \phi_i B_i
        m[idx + 1, idx + 1] = np.exp(-kLeft * boundary + offset)
        m[idx + 2, idx + 1] = - kLeft * np.exp(-kLeft * boundary + offset)

        # \phi_{i+1} A_{i+1}
        m[idx + 1, idx + 2] = - np.exp(kRight * boundary - offset)
        m[idx + 2, idx + 2] = - kRight * np.exp(kRight * boundary - offset)

        if idx < len(boundaries) - 2:
            # \phi_{i+1} B_{i+1}
            m[idx + 1, idx + 3] = np.exp(-kRight * boundary - offset)
            m[idx + 2, idx + 3] = -kRight * np.exp(-kRight * boundary - offset)

    solution = np.zeros((size,  1))    
    solution[0] = 1.0

    for i in range(size):
        string = '[' 
        for j in range(size):
            string += '{:.3f} + {:.3f}i'.format(m[i,j].real, m[i,j].imag)
            if j == size - 1:
                string += ']\n'
            else:
                string += ', '
        print(string)

    

    return m, solution





if __name__ == "__main__":
    numberOfBins = 3

    pref: Preference = Preference(numberOfBins, 92, 146, 238.050788, 234.043601)

    matrix, sol = createMatrix(pref)

    coefficients = matrix @ sol
    coefficients2 = np.linalg.solve(matrix, sol)


    m = csc_matrix(matrix)

    coeff = spsolve(m, sol)


    print(coefficients2)
    #print(coeff)

    print(coeff[-1] / coeff[0])


    #constants = np.linalg.solve(np.real(matrix), np.real(sol)) + 1.j*np.linalg.solve(np.imag(matrix), np.imag(sol))

    #print(constants)

