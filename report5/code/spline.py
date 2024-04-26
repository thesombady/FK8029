import numpy as np
import matplotlib.pyplot as plt
import math
from typing import Callable

import bspline as bs
import bspline.splinelab as bsl

a_0 = 1

@np.vectorize
def uniformSphere(r: float) -> float:
    if r > bSpline.outerRadius:
        return 0
    else:
        return 3 / (bSpline.outerRadius**3) 


@np.vectorize
def exactSphere(r: float) -> float:
    rprime = r + 1e-5
    r0 = bSpline.outerRadius
    return 1 /(2 * r0) * (3 - rprime**2 / r0**2)

@np.vectorize
def sphericalShell(r: float) -> float:
    if r > bSpline.outerRadius:
        return 0
    elif r < bSpline.innerRadius:
        return 0
    else:
        return 3 / ((np.power(bSpline.outerRadius, 3) - np.power(bSpline.innerRadius, 3)))

@np.vectorize
def exactShell(r: float) -> float:
    r0 = bSpline.innerRadius
    r1 = bSpline.outerRadius
    rprime = r + 1e-5
    return 3 / (r1 ** 3 - r0**2) * (r1 ** 2 / 2 - 1 / 3 * (r0**3/rprime + rprime ** 2 / 2))

@np.vectorize
def hydrogenGround(r: float) -> float:
    return np.exp(-2 * r / a_0 ) / (np.pi * a_0)

@np.vectorize
def exactHydrogen(r: float) -> float:
    rprime = r + 1e-5
    return (1/rprime - np.exp(-2*rprime) * (1 / rprime + 1)) / (2 *np.pi**2)

class bSpline:
    # Class that contains solution methods to the Poisson equation
    
    outerRadius = 10
    innerRadius = 5

    def __init__(self, order = 4):
        # __init__(order = 4) -> None; Constructor
        # order: int -> Order of the spline (k = 4 for cubic splines, k = 3 for quadratic splines, etc.)
        self.knots: np.ndarray = np.array([])
        self.physicalKnots: np.ndarray = np.array([])
        self.matrix: np.ndarray = np.array([[]])
        self.rhs: np.ndarray = np.array([])
        self.coeff: np.ndarray = np.array([])
        self.size: int = 0
        self.splineOrder = order


    def setKnots(self, knots: np.ndarray):
        # setKnot(knots: np.ndarray) -> None
        # Set the knots of the spline,
        # Not implemented this in python since I didnt need it
        # Look in the C++ code for implementation
        pass

    def initialize(self, numberOfPhysicalKnots: int):
        # initialize(numberOfPhysicalKnots: int) -> None
        # Set up the spline with the given number of physical knots
        k: int = self.splineOrder
        self.physicalKnots = [i / 10 for i in range(0, numberOfPhysicalKnots * 10)] + [numberOfPhysicalKnots - 1 - 0.00001]
        self.physicalKnots = np.array(sorted(list(self.physicalKnots)))
        self.knots = bsl.augknt(self.physicalKnots, self.splineOrder - 1)
        self.size = len(self.knots)
        self.spl = bs.Bspline(self.knots, self.splineOrder - 1)

    def save_matrix(self, matrix: np.ndarray, filename: str):
        # save_matrix(matrix: np.ndarray, filename: str) -> None
        # Save the matrix to the file (filename)
        # The matrix is saved in a tab separated format
        np.savetxt(
            filename,
            matrix,
            delimiter = '\t'
        )        

    def constructMatrix(self):
        # constructMatrix() -> None
        # Given initialized spline, construct the matrix
        # used to obtain the matrix to solve the Poisson equation
        temp = self.spl.collmat(self.physicalKnots, deriv_order = 2)
        temp[-1,:] = self.spl.collmat(self.physicalKnots[-2])
        print(self.spl.collmat(self.physicalKnots[-2], deriv_order = 0))
        temp = temp[:, :-1]
        temp = temp[:,1:]
        self.matrix = temp
        

    def solve(self, f: Callable, exact: Callable, filename: str):
        # solve(f: Callable, exact: Callable, filename: str) -> None
        # f: Callable -> Function that returns the value of the function at a given point
        # exact: Callable -> Function that returns the exact value of the function at a given point
        # filename: str -> Name of the file to save the matrix
        # exact is only used for plotting the exact solution
        # Solve the Poisson equation using the given function f, and save the matrix to the file
        rhs: np.ndarray = np.zeros(len(self.matrix))
        for i in range(len(self.matrix) - 1):
            rhs[i] = -f(self.physicalKnots[i]) * self.physicalKnots[i]
        rhs[-1] = 1
        res = np.linalg.solve(self.matrix, rhs)
        x = np.linspace(self.physicalKnots[0], self.physicalKnots[-2], 100)
        matrix = self.spl.collmat(x)
        matrix = matrix[:, :-1]
        matrix = matrix[:,1:]

        y = matrix @ res / x
        mat: np.ndarray = np.zeros((len(x), 3))
        mat[:,0] = x
        mat[:,1] = y
        ex = exact(x)
        mat[:,2] = ex 
        self.save_matrix(mat, filename)
        


if __name__ == '__main__':

    spline = bSpline(4)
    spline.initialize(15)
    spline.constructMatrix()
    spline.solve(hydrogenGround, "hydrogen.dat", exactHydrogen)


