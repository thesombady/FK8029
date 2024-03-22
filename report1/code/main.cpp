// Your First C++ Program

#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include "Eigen/Dense"

const double hbar_c = 197.326; // MeV fm
const double alpha = 1.0 / 137.035399;
const double V0 = 134; // MeV
const int helium = 4;
const double deltaR = 1.4;
const double atomicMass = 931.5; // MeV
const double c_squared = std::pow(299792458 * std::pow(10, 15), 2);

struct Preferences {
    int numberOfBins;
    int protonNumber;
    int neutronNumber;

    Preferences(int numberOfBins, int protonNumber, int neutronNumber) {
        this->numberOfBins = numberOfBins;
        this->protonNumber = protonNumber;
        this->neutronNumber = neutronNumber;
    };
};

struct Buffer {
    std::vector<std::vector<std::complex<double> > > matrix;

    std::vector<std::complex<double> >  constants;
    
    // Constructor
    Buffer(std::vector<std::vector<std::complex<double> > > matrix) {
        this->matrix = matrix;
    };
};

double findBoundary(double R, int Z, float radii, double Kinetic) {
    int iterations = 0;
    double r = R;
    
    while (std::abs(alpha / r * hbar_c * helium * Z - Kinetic) > 0.01 && iterations < 3000)  {
        r += 0.001;
        iterations += 1;
    }
    if (iterations == 3000) {
        throw std::invalid_argument("Could not find boundary");
    }
    return r - 0.001;
}

double potentialEnergy(double r, int Z, float radii) {
    if (r < radii) {
        return -V0;
    } else {
        return alpha / r * hbar_c * helium * Z;
    }
}

std::complex<double> computeK(double energy, double potential) {
    if (energy < potential) {
        return std::complex<double> (std::sqrt(2 * (potential - energy) * energy) / hbar_c, 0.0);
    }
    return std::complex<double> (0.0, std::sqrt(2 * (energy - potential) * energy)/hbar_c);
}

void computeMatrix(Preferences pref) {

    int size = 2 * pref.numberOfBins + 1;

    std::vector<std::vector<std::complex<double> > > matrix (size, std::vector<std::complex<double> >(size, std::complex<double> (0.0, 0.0)));

    //
    double energy = 47.3089;

    double boundaries [pref.numberOfBins];

    int massNumber = pref.protonNumber + pref.neutronNumber;

    double R = 1.4 * std::pow(massNumber, 1.0 / 3.0);

    double maximumDistance = findBoundary(R + 0.1, pref.protonNumber, R, energy);

    for (int i = 0; i < pref.numberOfBins; i++) {
        boundaries[i] = R + i * deltaR;
        //std::cout << boundaries[i] << std::endl;
    }

    matrix[0][0] = std::complex<double> (1.0, 0.0);

    int index = 0;

    std::complex<double> kList [2];

    for (int row = 1; row < size - 1; row += 2) {
        double boundary = boundaries[row / 2];
        std::complex<double> k_left = computeK(energy, potentialEnergy(boundary, pref.protonNumber, R));
        std::complex<double> k_right = computeK(energy, potentialEnergy(boundary, pref.protonNumber, R));


        if (row == 1) {
            kList[0] = k_left;
        }

        if (row > size - 3) {
            matrix[row][index] = std::exp(k_left * boundary);
            matrix[row][index + 1] = std::exp(- k_left * boundary);
            matrix[row][index + 2] = std::exp(k_right * boundary);

            matrix[row + 1][index] = k_left * std::exp(k_left * boundary);
            matrix[row + 1][index + 1] = -k_left * std::exp(-k_left * boundary);
            matrix[row + 1][index + 2] = k_right * std::exp(k_right * boundary);

            kList[1] = k_right;

            break;
        }

        matrix[row][index] = std::exp(k_left * boundary);
        matrix[row][index + 1] = std::exp(- k_left * boundary);
        matrix[row][index + 2] = std::exp(k_right * boundary);
        matrix[row][index + 3] = std::exp(-k_right * boundary);

        matrix[row + 1][index] = k_left * std::exp(k_left * boundary);
        matrix[row + 1][index + 1] = -k_left * std::exp(-k_left * boundary);
        matrix[row + 1][index + 2] = k_right * std::exp(k_right * boundary);
        matrix[row + 1][index + 3] = -k_right * std::exp(-k_right * boundary);

        index += 2;
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (j != size - 1) {
                std::cout << matrix[i][j] << ", ";
            } else {
                std::cout << matrix[i][j] << std::endl;;
            }
        }
        std::cout << std::endl;
    }

    //std::complex<double> k_0 = computeK(0.0, potentialEnergy(R, pref.protonNumber, R)); // First boundary

    //std::complex<double> k_n = computeK(0.0, potentialEnergy(r, pref.protonNumber, R)); // Last Boundary


    //return Buffer(matrix);
}


int main() {
    int numberOfBins = 2;

    //double disc [numberOfBins][numberOfBins];

    Preferences pref = Preferences(numberOfBins, 92, 146);
    
    computeMatrix(pref);    

    return 0;
}