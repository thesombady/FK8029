#include <arm_neon.h>
#include <codecvt>
#include <iostream>
#include <mutex>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "Eigen/Dense" // Add to path by PATH=$PATH:/PATH_TO_EIGEN?
#include "Eigen/src/Core/CommaInitializer.h"
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Eigenvalues/EigenSolver.h"
#include "Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h"

const float xmin = -6.0;
const float xmax = 6.0;

/**
    Potential energy function
    @param[in] x Position
    @return Potential energy at position x, V(x)
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

enum Save {
    probability,
    wave
};

/**
    Creates a matrix for the 1D Schrödinger equation
    @param[in] n Number of points
    @param[in] ve Pointer to the potential energy function
    @return Matrix for the 1D Schrödinger equation, LHS
*/
Eigen::MatrixXd createMatrix(int n, float (*ve)(float), Method m) {
    Eigen::MatrixXd hamiltonian = Eigen::MatrixXd::Zero(n, n);

    float dx = (xmax - xmin) / (n - 1);

    float x = xmin;

    if (m == three_point) {
        for (int i = 0; i < n; i++) {
            float c1 = -1;
            if (i == 0) {
                hamiltonian(i, 1) = c1;
            } else if (i == n - 1) {
                hamiltonian(i, i - 1) = c1;
            } else {
                hamiltonian(i, i - 1) = c1;
                hamiltonian(i, i + 1) = c1;
            }
            hamiltonian(i, i) = 2.0 + dx * dx * ve(x);

            x += dx;
        }
        return hamiltonian / (dx *dx);
    } else if (m == five_point) {
        float c1 = -16;
        float c2 = 1;
        for (int i = 0; i < n; i++) {
            if (i == 0) {
                hamiltonian(i, 1) = c1; // c_1
                hamiltonian(i, 2) = c2; // c_2
            } else if (i == 1) {
                hamiltonian(i, i - 1) = c1; // c_(-1)
                hamiltonian(i, i + 1) = c1; // c_1
                hamiltonian(i, i + 2) = c2; // c_2
            } else if (i == n - 2) {
                hamiltonian(i, i - 2) = c2; // c_(-2)
                hamiltonian(i, i - 1) = c1; // c_(-1)
                hamiltonian(i, i + 1) = c1; // c_1
                
            } else if (i == n - 1) {
                hamiltonian(i, i - 2) = c2; // c_(-2)
                hamiltonian(i, i - 1) = c1; // c_(-1)
                
            } else {
                hamiltonian(i, i - 2) = c2; // c_(-2)
                hamiltonian(i, i - 1) = c1; // c_(-1)
                hamiltonian(i, i + 1) = c1; // c_(1)
                hamiltonian(i, i + 2) = c2; // c_(2)
                
            }
            hamiltonian(i, i) = 30 + (12 * dx *dx) * ve(x); // c_0
            x += dx;
        }
        return hamiltonian / (12 * dx * dx);
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
    @param[in] std::string filename - The filename to save the data
    @param[in] Save s - Which type of information to save
*/
void ThreePointTest(int n, float (*ve)(float), std::string filename, Save s) {
    
    Method method = three_point;

    int numberOfPoints = n;

    Eigen::MatrixXd hamiltonian = createMatrix(numberOfPoints, ve, method) / 2;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonian);

    Eigen::VectorXd values = solver.eigenvalues();
    
    Eigen::MatrixXd temp = solver.eigenvectors();
    
    Eigen::VectorXd xlist(n);

    float dx = (xmax - xmin) / (n - 1);
    float x = xmin;
    for (int i = 0; i < n; i++) {
        xlist(i) = x;
        x += dx;
    }

    std::ofstream file(filename);

    file << "x\tv1\tv2\tv3" << std::endl;
    for (int i = 0; i < n; i++) {
        file << xlist(i);
        for (int j = 0; j < 4; j++) {
            file << "\t" << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
            if (s == probability) {
                file << std::pow(temp(i, j), 2);
            } else {
                file << temp(i, j);
            }
        }
        file << std::endl;
    }
    for (int i = 0; i < 4; i++) {
        std::cout << "Eigen-Value: " << i << "=" <<  std::setprecision(std::numeric_limits<long double>::digits10 + 1) << values(i) << std::endl;
    }
    file.close();
}

/**
    Computes the eigen-values to the schrödinergs hamiltonian
    using a five point stencil
    @param[in] int n - The discretization number
    @param[in] float (*ve)(float) - The potential function
    @param[in] std::string filename - The filename to save the data
    @param[in] Save s - Which type of information to save
*/
void FivePointTest(int n, float (*ve)(float), std::string filename, Save s) {
    
    Method method = five_point;

    int numberOfPoints = n;

    Eigen::MatrixXd hamiltonian = createMatrix(numberOfPoints, ve, method) / 2;
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonian);

    Eigen::VectorXd values = solver.eigenvalues();

    Eigen::MatrixXd temp = solver.eigenvectors();
    
    Eigen::VectorXd xlist(n);

    float dx = (xmax - xmin) / (n - 1);
    float x = xmin;
    for (int i = 0; i < n; i++) {
        xlist(i) = x;
        x += dx;
    }

    std::ofstream file(filename);
    int sign = 1;

    file << "x\tv1\tv2\tv3" << std::endl;
    for (int i = 0; i < n; i++) {
        file << xlist(i);
        for (int j = 0; j < 4; j++) {
            if (j == 1) {// FIXED: Sign of wave-funciton when plotting
                sign = 1;
            } else if (j == 3) {
                sign = -1;
            }

            file << "\t" << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
            if (s == probability) {
                file << std::pow(temp(i, j), 2);
            } else {
                file << sign * temp(i, j);
            }
        }
        file << std::endl;
    }

    for (int i = 0; i < 4; i++) {
        std::cout << "Eigen-Value: " << i << "=" << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << values(i) << std::endl;
    }

    file.close();
}

/**
    InversePowerIteration; solve the schrödinger equation
    @param[in] int n - The discreziation number
    @param[in] Method m - Which stencil to use (three_point or five_point)
    @param[in] float (*ve)(float) - Potential function to use
    @param[in] flaot shift_ - Shift value to find the closest eigen-value
    @return EigenResult - The eigen-value, eigen-vector and the number of iterations
*/
EigenResult InversePowerIteration(int n, Method m, float (*ve)(float), float shift_) {
    double tol = 1e-6;
    
    Eigen::MatrixXd hamiltonian = createMatrix(n, ve, m) / 2.0; 
    Eigen::MatrixXd shift = Eigen::MatrixXd::Identity(n, n) * shift_;

    Eigen::MatrixXd matrix = hamiltonian - shift; // Corrected hamiltonian to get correct eigenvalues

    Eigen::VectorXd guess = Eigen::VectorXd::Random(n).normalized();

    Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd> > Lu(matrix); // Change?

    int counter = 0;

    Eigen::VectorXd newGuess(n);

    double norm = 0;

    while (true) {
        newGuess = Lu.solve(guess).normalized(); // FIX: Implement lu outside loop
        // std::cout << (guess - newGuess).norm() << std::endl;

        if (counter > 1000) {
            throw std::invalid_argument("Solution did not converge within limit");
            break;
        }

        if ((guess - newGuess).norm() < tol ) {
            break;
        }
        
        guess = newGuess.normalized();

        counter ++;        
        
    }

    EigenResult res;
    res.value = newGuess.transpose() * hamiltonian * newGuess ;
    res.vector = newGuess;
    res.iterations = counter;
    return res;

}


/**
    InversePowerTest; solve the schrödinger equation using the inverse power method
    Using the five-point stencil in order to parametrize the Schrödinger equation.
    @param[in] int n - The discreziation number
    @param[in] float (*ve)(float) - Potential function to use
    @param[in] std::string filename - The filename to save the data
    @param[in] std::vector<float> shift - Shift values to find the closest eigen-value
*/
void inversePowerTest(int n, float (*ve)(float), std::string filename, std::vector<float> shift, Save s) {
    
    std::vector<EigenResult> res(5);

    std::ofstream file(filename);
    
    for (int i = 0; i < shift.size(); i++ ) {
        res[i] = InversePowerIteration(n, ::five_point, ve, shift[i]);
        std::cout << "Eigen value: " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << res[i].value << std::endl;
        //std::cout << res.vector << std::endl;
        std::cout << "No. iterations: " << res[i].iterations << std::endl;
    } 


    std::vector<std::string> shift_;
    //std::string shift_[6] = {"x", "v0", "v1", "v2", "v3", "v4"};
    for (int i = 0; i < shift.size() + 1; i++) {
        if (i == 0) {
            file << "x";
        } else {
            file << "v" + std::to_string(i - 1);
        }
        if (i < shift.size() + 1) {
            file << "\t";
        } else {
            file << std::endl;
        }
    }

    double dx = (xmax - xmin) / ( n - 1 );

    file << std::endl;

    for (int i = 0; i < n; i++){
        for (int j = 0; j < shift.size() + 1; j++) {
            if ( j == 0 ) {
                // We are at x and thus do that a bit differently
                file << xmin + i * dx;
            } else {
                // psi^2 = |psi|^2
                file << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
                if (s == probability) {
                    file << std::pow(res[j - 1].vector(i), 2);
                } else {
                    file  << res[j - 1].vector(i);
                }
            }
            if (j < shift.size() ) {
                file << "\t";
            } else {
                file << std::endl;
            }
        }
    }

    file.close();

    std::ofstream file2("compare_" + filename);

    file2 << "i\tN" << std::endl;

    for (int i = 0; i < shift.size(); i++) {
        file2 << i << "\t" <<  std::setprecision(std::numeric_limits<long double>::digits10 + 1) << std::abs(res[i].value - (i + 0.5)) << std::endl;
    }

    file2.close();

}


/**
    compareDirect; compares the eigenvalues of the direct method and the exact values
    Compares both the three-point and five-point stencil
    @param[in] int n - The number of points
    @param[in] float (*ve)(float) - The potential function
    @param[in] std::string filename - The filename to save the data
*/
void compareDirect(int n, float (*ve)(float), std::string filename) {
    
    Eigen::VectorXd exactValues(n);
    for (int i = 0; i < n; i++) {
        exactValues(i) = (1.0 / 2.0 + i); // * hbar omega but we have it in units of hbar omega
    }

    Eigen::MatrixXd threePointMatrix = createMatrix(n, ve, ::three_point) / 2.0;
    Eigen::MatrixXd fivePointMatrix = createMatrix(n, ve, ::five_point) / 2.0;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> threePointSolver(threePointMatrix);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> fivePointSolver(fivePointMatrix);

    Eigen::VectorXd threePointValues = threePointSolver.eigenvalues();
    Eigen::VectorXd fivePointValues = fivePointSolver.eigenvalues();


    std::ofstream file(filename);

    file << "state\tdE_3pt\tdE_5pt" << std::endl;

    for (int i = 0; i < n; i++) {
        file << i << "\t" <<std::setprecision(std::numeric_limits<long double>::digits10 + 1) << std::abs( exactValues(i) - threePointValues(i) ) << "\t" << std::abs( exactValues(i) - fivePointValues(i) ) << std::endl;
    }

    file.close();
}


/**
    compareExact; compares the exact eigenvalues with the computed eigenvalues
    using the five-point stencil
    @param[in] int n - The number of points
    @param[in] float (*ve)(float) - The potential function
    @param[in] std::string filename - The filename to save the data
*/
void compareExact(int n, float (*ve)(float), std::string filename) {
    
    Eigen::VectorXd exactValues(n);
    for (int i = 0; i < n; i++) {
        exactValues(i) = (1.0 / 2.0 + i); // * hbar omega but we have it in units of hbar omega
    }

    Eigen::MatrixXd fivePointMatrix = createMatrix(n, ve, ::five_point) / 2.0;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> fivePointSolver(fivePointMatrix);

    Eigen::VectorXd fivePointValues = fivePointSolver.eigenvalues();


    std::ofstream file(filename);

    file << "state\texact\tdm" << std::endl;

    for (int i = 0; i < n; i++) {
        file << i << "\t" << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << exactValues(i) << "\t";
        file << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << fivePointValues(i) << std::endl;
    }

    file.close();
}

/**
    potential2; modified potential function
    @param[in] z - Position
    @return float - the potential at the point `z`
*/
float potential2(float z) {
    return z * z + 10 * std::exp(-std::abs(z));
}

/* // TODO: Implement ladder operators in spair time.
std::vector<Eigen::MatrixXd> createOperators(int n) {
    Eigen::MatrixXd creationOp = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd destructOp = Eigen::MatrixXd::Zero(n, n);

    float dx = (xmax - xmin) / (n - 1);

    float pre = 1 / (2 * dx);
    
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            creationOp(i,i + 1) = pre;
            destructOp(i, i + 1) = -pre;
        } else if (i == n - 1) {
            creationOp(i, i - 1) = -pre;
            destructOp(i, i - 1) = pre;
        } else {
            creationOp(i, i) = xmin + i * dx;
            destructOp(i, i) = xmin + 1 * dx;
        }
    }

    std::vector<Eigen::MatrixXd> operators(2);
    operators[0] = creationOp / std::sqrt(2);
    operators[1] = destructOp / std::sqrt(2);

    return operators;
}
*/

int main() {

    // ThreePointTest(500, &potential, "3pt.dat", ::probability);   
    // FivePointTest(500, &potential, "5pt.dat", ::probability);   
    // compareDirect(500, &potential, "compare.dat");
    // compareExact(1000, &potential, "exact.dat");
    
    std::vector<float> shift(5);
    shift[0] = 0.4;
    shift[1] = 1.2;
    shift[2] = 2.3;
    shift[3] = 3.3;
    shift[4] = 4.4;
    // inversePowerTest(500, &potential, "inverse.dat", shift, ::probability);
    std::vector<float> newShift(2);
    newShift[0] = 2.8;
    newShift[1] = 4.1;
    // FivePointTest(500, &potential2, "5pt_modified.dat", ::probability);
    //FivePointTest(500, &potential2, "dm_wave.dat", ::wave);
    //inversePowerTest(500, &potential2, "ipm_wave.dat", newShift, ::wave);

    FivePointTest(500, &potential2, "prob_pert.dat", ::probability);
    FivePointTest(500, &potential, "prob_act.dat", ::probability);
    
    return 0;
}
