#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>

#include "Eigen/Dense" // Add to path by PATH=$PATH:/PATH_TO_EIGEN?
#include "Eigen/src/Core/CommaInitializer.h"
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Eigenvalues/EigenSolver.h"
#include "Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h"

const float xmin = -6.0;
const float xmax = 6.0;

/**
 * Potential energy function
 * @param[in] x Position
 * @return Potential energy at position x, V(x)
*/
float potential(float x) {
    return x * x;
}

struct EigenResult {

    Eigen::VectorXd vector;
    float value;
    int iterations;
    std::vector<float> x;
};

enum Method {
    none,
    three_point,
    five_point
};

/**
 * Creates a matrix for the 1D Schrödinger equation
 * @param[in] n Number of points
 * @param[in] ve Pointer to the potential energy function
 * 
 * @return Matrix for the 1D Schrödinger equation, LHS
*/
Eigen::MatrixXd createMatrix(int n, float (*ve)(float), Method m) {
    Eigen::MatrixXd hamiltonian = Eigen::MatrixXd::Zero(n, n);

    float dx = (xmax - xmin) / (n - 1);

    float x = xmin;

    if (m == three_point) {
        for (int i = 0; i < n; i++) {
            if (i == 0) {
                hamiltonian(i, 1) = -1.0 / (dx * dx);
            } else if (i == n - 1) {
                hamiltonian(i, i - 1) = -1.0 / (dx * dx);
            } else {
                hamiltonian(i, i - 1) = -1.0 / (dx * dx);
                hamiltonian(i, i + 1) = -1.0 / (dx * dx);
            }
            hamiltonian(i, i) = 2.0 / (dx * dx) + ve(x);

            x += dx;
        }
    } else if (m == five_point) {
        float c_2 = 1; 
        float c_1 = -16; 
        float c1 = -16;
        float c2 = 1;
        for (int i = 0; i < n; i++) {
            if (i == 0) {
                hamiltonian(i, 1) = c1; // c_1
                hamiltonian(i, 2) = c2; // c_2
            } else if (i == 1) {
                hamiltonian(i, i - 1) = c_1; // c_(-1)
                hamiltonian(i, i + 1) = c1; // c_1
                hamiltonian(i, i + 2) = c2; // c_2
            } else if (i == n -2) {
                hamiltonian(i, i - 2) = c_2; // c_(-2)
                hamiltonian(i, i - 1) = c_1; // c_(-1)
                hamiltonian(i, i + 1) = c_1; // c_1
                
            } else if (i == n - 1) {
                hamiltonian(i, i - 2) = c_2; // c_(-2)
                hamiltonian(i, i - 1) = c_1; // c_(-1)
                
            } else {
                hamiltonian(i, i - 2) = c_2; // c_(-2)
                hamiltonian(i, i - 1) = c_1; // c_(-1)
                hamiltonian(i, i + 1) = c1; // c_(1)
                hamiltonian(i, i + 2) = c_2; // c_(2)
                
            }
            hamiltonian(i, i) = 30 + (12 * dx *dx) * ve(x); // c_0
            x += dx;
        }
        hamiltonian *= 1/(12 * dx * dx);
    } else {
        throw std::invalid_argument("Invalid method for computing the matrix");
    }
    //std::cout << hamiltonian << std::endl;

    return hamiltonian;
}


/**
    Computes the eigen-values to the schrödinergs hamiltonian
    using a three point stencil
    @param[in] int n - The discretization number
    @param[in] float (*ve)(float) - The potential function
    @return Eigen::VectorXd - The vector contaning all the eigen-values
*/
void ThreePointTest(int n, float (*ve)(float)) {
    
    Method method = three_point;

    int numberOfPoints = n;

    Eigen::MatrixXd hamiltonian = createMatrix(numberOfPoints, ve, method) / 2;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonian);

    Eigen::VectorXd values = solver.eigenvalues();
    
    Eigen::MatrixXd temp = solver.eigenvectors().reshaped(n,n);
    
    Eigen::VectorXd xlist(n);

    float dx = (xmax - xmin) / (n - 1);
    float x = xmin;
    for (int i = 0; i < n; i++) {
        xlist(i) = x;
        x += dx;
    }

    std::string filename = "3pt.dat";

    std::ofstream file(filename);

    file << "x\tv1\tv2\tv3" << std::endl;
    for (int i = 0; i < n; i++) {
        file << xlist(i) << "\t" << temp(0, i) << "\t" << temp(1, i) << "\t" << temp(2, i) << std::endl;
    }
    for (int i = 0; i < 4; i++) {
        std::cout << "Eigen-Value: " << i << "=" << values(i) << std::endl;
    }
    file.close();
}

/**
    Computes the eigen-values to the schrödinergs hamiltonian
    using a five point stencil
    @param[in] int n - The discretization number
    @param[in] float (*ve)(float) - The potential function
*/
void FivePointTest(int n, float (*ve)(float)) {
    
    Method method = five_point;

    int numberOfPoints = n;

    Eigen::MatrixXd hamiltonian = createMatrix(numberOfPoints, ve, method) / 2;
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonian);

    Eigen::VectorXd values = solver.eigenvalues();

    Eigen::MatrixXd temp = solver.eigenvectors().reshaped(n,n);
    
    Eigen::VectorXd xlist(n);

    float dx = (xmax - xmin) / (n - 1);
    float x = xmin;
    for (int i = 0; i < n; i++) {
        xlist(i) = x;
        x += dx;
    }

    std::string filename = "5pt.dat";

    std::ofstream file(filename);

    file << "x\tv1\tv2\tv3" << std::endl;
    for (int i = 0; i < n; i++) {
        file << xlist(i) << "\t" << temp(0, i) << "\t" << temp(1, i) << "\t" << temp(2, i) << std::endl;
    }
    for (int i = 0; i < 4; i++) {
        std::cout << "Eigen-Value: " << i << "=" << values(i) << std::endl;
    }
    file.close();
}

/**
    InversePowerIteration; solve the schrödinger equation
    @param[in] int n - The discreziation number
    @param[in] Method m - Which stencil to use (three_point or five_point)
    @param[in] float (*ve)(float) - Potential function to use
    @param[in] flaot shift_ - Shift value to find the closest eigen-value
    @return flaot - The eigenvalue given a shift
*/
EigenResult InversePowerIteration(int n, Method m, float (*ve)(float), float shift_) {
    double tol = 1e-6;
    
    Eigen::MatrixXd hamiltonian = createMatrix(n, ve, m); 
    Eigen::MatrixXd shift = Eigen::MatrixXd::Identity(n, n) * shift_;

    Eigen::MatrixXd matrix = hamiltonian/2 - shift; // Corrected hamiltonian to get correct eigenvalues

    Eigen::VectorXd guess = Eigen::VectorXd::Ones(n).normalized();

    Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd> > Lu(matrix); // Change?

    int counter = 0;

    Eigen::VectorXd newGuess(n);

    double norm = 0;

    while (true) {

        newGuess = Lu.solve(guess).normalized(); // FIX: Implement lu outside loop

        if ((guess - newGuess).norm() < tol) {
            EigenResult res;
            res.value = newGuess.transpose() * hamiltonian / 2 * newGuess ;
            res.vector = newGuess;
            res.iterations = counter;
            return res;
        }

        if (counter > 1000) {
            throw std::invalid_argument("Solution did not converge within limit");
            break;
        }

        guess = newGuess.normalized();

        counter ++;        
    }
}

int main() {
    ThreePointTest(100, &potential);   
    FivePointTest(100, &potential);   
    
    //std::cout << test1 << std::endl;
    // std::cout << test2 << std::endl;

    /* EigenResult res;

    // We know the eigenvalues approximately:

    float shift[5] = {0.1, 1.4, 2.4, 3.6, 4.3};
    // TODO: Save information in dat file
    
    for (int i = 0; i < 5; i++ ) {
        res = InversePowerIteration(1500, ::five_point, &potential, shift[i]);
        std::cout << "Eigen value: " << res.value << std::endl;
        //std::cout << res.vector << std::endl;
        std::cout << "No. iterations: " << res.iterations << std::endl;
    } */
    return 0;
}
