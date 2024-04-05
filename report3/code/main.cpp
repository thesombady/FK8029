#include <iostream>
#include <string>

#include "Eigen/Dense"

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

/**
 * Creates a matrix for the 1D Schrödinger equation
 * @param[in] n Number of points
 * @param[in] ve Pointer to the potential energy function
 * 
 * @return Matrix for the 1D Schrödinger equation, LHS
*/
Eigen::MatrixXd createMatrix(int n, float (*ve)(float) ) {
    Eigen::MatrixXd hamiltonian = Eigen::MatrixXd::Zero(n, n);

    float dx = (xmax - xmin) / (n - 1);
    std::cout << "dx: " << dx << std::endl;

    float x = xmin;

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

    std::cout << hamiltonian << std::endl;

    return hamiltonian;
}


int main() {

    std::cout << "Hello, World!" << std::endl;

    Eigen::MatrixXd hamiltonian = createMatrix(10, &potential);


    return 0;
}