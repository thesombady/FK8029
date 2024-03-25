#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Sparse"

const double hbar_c = 197.326; // MeV fm
const double alpha = 1.0 / 137.035399;
const float V0 = 134; // MeV
const int protonHelium = 2; // Two protons
const double speedOfLight = 299792458; // m/s
const long double c_squared = std::pow(speedOfLight * std::pow(10, 15), 2); // fm^2 / s^2
const double massOfHelium = 3727.379; // MeV 
const double amuToMeV = 931.5; // MeV

/**
 * Prints a string to the console
 * @param[in] s the string to print
*/
void println(std::string s) {
    std::cout << s << std::endl;
}

/**
 * Preferences for the decay
 * @param[in] numberOfBins the number of bins used for the discretisation
 * @param[in] protonNumber the number of protons in the nuclei
 * @param[in] neutronNumber the number of neutrons in the nuclei
 * @param[in] massParent the mass of the parent nuclei in amu
 * @param[in] massDaughter the mass of the daughter nuclei in amu
 * @param[in] name the name of the decay
 * @param[in] r0 the scaling factor for the radii
*/
struct Preferences {
    int numberOfBins;
    int protonNumber;
    int neutronNumber;
    long double massParent; // MeV
    long double massDaughter; // MeV
    double r0;
    std::string name;

    Preferences(int numberOfBins, int protonNumber, int neutronNumber, double massParent, double massDaughter, std::string name, double r0) {
        this->numberOfBins = numberOfBins;
        this->protonNumber = protonNumber;
        this->neutronNumber = neutronNumber;
        this->massParent = massParent * amuToMeV;
        this->massDaughter = massDaughter * amuToMeV;
        this->name = name;
        this->r0 = r0;
    };
};

struct Buffer {
    double long kineticEnergy; // MeV
    double long halfLife; // s
    std::string name;
};

/**
 * Computes the distance of at which the potential energy is equal to the kinetic energy
 * @param[in] Z the number of protons in the nuclei
 * @param[in] kineticEnergy the kinetic energy of the alpha particle
 * @returns the distance of at which the potential energy is equal to the kinetic energy
*/
long double findBoundary(int Z, long double kineticEnergy) {
    // V = K
    // V = alpha / r * hbar_c * helium * Z = K
    // alpha / energy * hbar_c * helium * Z = r

    return alpha / kineticEnergy * hbar_c * protonHelium* (Z);
}

/**
 * Computes the potential energy at a given distance r from the nuclei
 * @param[in] r the distance from the nuclei
 * @param[in] Z the number of protons in the nuclei
 * @param[in] radii the radii of the nuclei
 * @param[in] maximumDistance where the potential energy is equal to the kinetic energy
 * @returns the potential energy at a given distance r from the nuclei
*/
long double potentialEnergy(long double r, int Z, float radii, long double maximumDistance) {
    if (r < radii) {
        return -V0 + alpha / r * hbar_c * protonHelium * Z;
    } else if (r > maximumDistance) {
        return 0;
    } else {
        return alpha / r * hbar_c * protonHelium * Z;
    }
}

/**
 * Computes the wave number k for a given energy and potential energy
 * @param[in] energy the kinetic energy of the alpha particle
 * @param[in] potentialEnergy the potential energy at a given distance r from the nuclei
 * @returns the wave number k which might be complex
*/
std::complex<long double> computeK(long double energy, long double potentialEnergy) {
    if (energy < potentialEnergy) {
        return std::complex<long double>(std::sqrt(2 * massOfHelium * (potentialEnergy - energy)) / hbar_c, 0.0);
    }
    return std::complex<long double> (0.0, std::sqrt(2 * massOfHelium * (energy - potentialEnergy)) / hbar_c);
}

/**
 * Computes the energy of the alpha particle
 * @param[in] pref the preferences for the decay, e.g. number of bins, proton number, neutron number, mass of parent, mass of daughter, name
 * @returns the energy of the alpha particle
*/
long double computeAlphaEnergy(Preferences pref) {
    return pref.massParent - (pref.massDaughter + massOfHelium);
}

/**
 * Computes the matrix for the given potential and energy
 *
 * @param[in] pref the preferences for the decay, e.g. number of bins, proton number, neutron number, mass of parent, mass of daughter, name
 * @returns Buffer the buffer containing the half-life, kinetic energy and name of the decay
*/
Buffer solve(Preferences pref) {

    int size = 2 * (pref.numberOfBins - 1) + 1;
    
    // We compute the energy of the alpha particle
    long double energy = computeAlphaEnergy(pref); // MeV

    int massNumber = pref.protonNumber + pref.neutronNumber;

    long double minimumDistance = pref.r0 * std::pow(massNumber - 4, 1.0 / 3.0); // massNumber - 4 is the massNumber of the daughter nuclei

    long double maximumDistance =  findBoundary(pref.protonNumber - 2, energy);

    long double stepSize = ( maximumDistance - minimumDistance ) / ( pref.numberOfBins - 2 );

    long double boundaries [2 * pref.numberOfBins - 2];
    double offset = stepSize * 0.01;

    long double b;
    for (int i = 0; i < pref.numberOfBins - 1; i++) {
        b = minimumDistance + i * stepSize;
        boundaries[2 * i] = b - offset;
        boundaries[2 * i + 1] = b + offset;
    }
    boundaries[pref.numberOfBins - 1] = maximumDistance;

    Eigen::MatrixXcf m = Eigen::MatrixXcf::Zero(size, size);
    m(0,0) = 1;

    // We have the size 2N - 1, where N is the number of bins
    // We can fill from m(1,0) -> m(2*N - 1, 2*N - 2)

    for (int idx = 0; idx < 2 * pref.numberOfBins - 2; idx += 2) {
        long double boundary = boundaries[idx] + offset;
        
        std::complex<long double> kLeft, kRight;

        kLeft = computeK(energy, potentialEnergy(boundaries[idx], pref.protonNumber - 2, minimumDistance, maximumDistance));
        kRight = computeK(energy, potentialEnergy(boundaries[idx + 1], pref.protonNumber - 2, minimumDistance, maximumDistance));

        // Now we fill the matrix elements

        // \phi_i A_i
        m(idx + 1, idx) = std::exp(kLeft * (boundary));
        m(idx + 2, idx) = kLeft * std::exp(kLeft * (boundary));

        // \phi_i B_i
        m(idx + 1, idx + 1) = std::exp(-kLeft * (boundary));
        m(idx + 2, idx + 1) = - kLeft * std::exp(-kLeft * (boundary));

        // \phi_{i+1} A_{i+1}
        m(idx + 1, idx + 2) = - std::exp(kRight * (boundary));
        m(idx + 2, idx + 2) = - kRight * std::exp(kRight * (boundary));

        if (idx < 2 * pref.numberOfBins - 4){
            // Only the last row is different, because the matrix is constructed such that;
            // phi_n = a_n e^{ik_nx*)} + 0 * b_n e^{-ik_nx}
            // \phi_{i+1} B_{i+1}
            m(idx + 1, idx + 3) = std::exp(-kRight * (boundary));
            m(idx + 2, idx + 3) = -kRight * std::exp(-kRight * (boundary));
        }
    }

    // Solving the system of equations and producing the half life
    Eigen::VectorXcf RHS = Eigen::VectorXcf::Zero(size);

    RHS(0) = 1.0; // The "initial" condition

    Eigen::VectorXcf coefficients;
    coefficients = m.partialPivLu().solve(RHS);

    std::complex<long double> rate;

    rate = computeK(energy, potentialEnergy(boundaries[2 * pref.numberOfBins - 3], pref.protonNumber - 2, minimumDistance, maximumDistance));
    rate /= computeK(energy, potentialEnergy(boundaries[0], pref.protonNumber - 2, minimumDistance, maximumDistance));

    long double transmission;
   
    transmission = std::pow(std::abs(coefficients[size - 1] / coefficients[0]), 2) * (rate.real());// + rate.imag());

    /*
        The half-life is given by the following formula
        T = ln(2) * 2 * R / (v_eff * Transmission)
        where v is the velocity of the alpha particle
    */
    double effectiveMass = massOfHelium * pref.massDaughter / (pref.massDaughter + massOfHelium);

    long double v_eff = std::sqrt(2 * energy / effectiveMass); // c

    v_eff = v_eff * speedOfLight; // m / s

    long double tau = 2 * (minimumDistance * std::pow(10, -15)) / ( transmission * v_eff); // s

    long double halfLife = std::log(2) * tau;

    Buffer buffer;
    buffer.halfLife = halfLife;
    buffer.kineticEnergy = energy;
    buffer.name = pref.name;

    return buffer;
}

int main() {
    println("Input whether you want to perform the calculations for the isotope chain (y/n):");
    std::string input;

    std::cin >> input;
    if (input == "y") {

        int numberOfBins = 100;

        double deltaR = 1.4;

        Preferences pref1 = Preferences(numberOfBins, 82, 132, 213.9952014, 209.9841885, "Po-214", deltaR); // Po-214 -> Pb-210
        // Preferences pref2 = Preferences(numberOfBins, 94, 148, 242.059, 238.05078826, "Pu-242", deltaR); // Pu-242 -> U-238
        //Preferences pref2 = Preferences(numberOfBins, 88, 135, 223.0185007,219.0094802, "Ra-223", deltaR); // Ra-223 -> Rn-219
        Preferences pref2 = Preferences(numberOfBins, 86, 133, 219.00948, 214.9994200, "Rn-219", deltaR); // Radon-219 -> Polonium-215

        Preferences pref3 = Preferences(numberOfBins, 92, 146, 238.050788, 234.043601, "U-238", deltaR); // Uranium-238 -> Thorium-234
        Preferences pref4 = Preferences(numberOfBins, 90, 142, 230.03313, 226.025408, "Th-232", deltaR); // Thorium-230 -> Radium-226
        Preferences pref5 = Preferences(numberOfBins, 88, 136, 226.025408, 222.0175763, "Ra-226", deltaR); // Radium-226 -> Radon-222
        Preferences pref6 = Preferences(numberOfBins, 86, 136, 222.0175763, 218.0089730, "Rn-222", deltaR); // Radon-218 -> Polonium-214

        Buffer buffers [6];

        buffers[0] = solve(pref1);
        buffers[1] = solve(pref2);
        buffers[2] = solve(pref3);
        buffers[3] = solve(pref4);
        buffers[4] = solve(pref5);
        buffers[5] = solve(pref6);
        
        std::ofstream file("data.dat");

        file << "K\t t\t name" << std::endl;

        for (int i = 0; i < 6; i++) {
            file << buffers[i].kineticEnergy <<  "\t" << buffers[i].halfLife << "\t" << buffers[i].name << std::endl;
        }

        file.close();

    } else {
        println("Input the decay: ");
        std::vector<Preferences> pref;
        std::cin >> input;

        println("Input the radii scaling: ");
        double r0;
        std::cin >> r0;

        int maximumBins = 21;


        if (input == "cf-250") {
            for (int i = 1; i < maximumBins + 1; i++) {
                pref.push_back(Preferences(i * 10, 98, 152, 250.0764045, 246.0672237, "Cf-250", r0));
            }
        } else if (input == "pu-242") {
            for (int i = 1; i < maximumBins + 1; i++) {
                pref.push_back(Preferences(i * 10, 94, 148, 242.059, 238.05078826, "Pu-242", r0));
            }
        } else if (input == "u-238") {
            for (int i = 1; i < maximumBins + 1; i++) {
                pref.push_back(Preferences(i * 10, 92, 146, 238.050788, 234.043601, "U-238", r0));
            }
        } else if (input == "th-232") {
            for (int i = 1; i < maximumBins + 1; i++) {
                pref.push_back(Preferences(i * 10, 90, 142, 230.03313, 226.025408, "Th-232", r0));
            }
        } else if (input == "ra-226") {
            for (int i = 1; i < maximumBins + 1; i++) {
                pref.push_back(Preferences(i * 10, 88, 136, 226.025408, 222.0175763, "Ra-226", r0));
            }
        } else if (input == "rn-222") {
            for (int i = 1; i < maximumBins + 1; i++) {
                pref.push_back(Preferences(i * 10, 86, 136, 222.0175763, 218.0089730, "Rn-222", r0));
            }
        } else {
            println("Invalid input");
            throw std::invalid_argument("Invalid input");
            return 0;
        }

        Buffer b [pref.size()];

        for (int i = 0; i < pref.size(); i++) {
            Buffer buffer = solve(pref[i]);
            //println("Half-life for " + pref[i].name + " with " + std::to_string(pref[i].numberOfBins) + " bins: " + std::to_string(buffer.halfLife));
            b[i] = buffer;
        }

        std::string filename;
        
        filename = "data";

        filename.append("_");
        filename.append(input);
        filename.append("_");
        filename.append(std::to_string(float(r0)));
        filename.append(".dat");

        std::ofstream file(filename);

        file << "binSize\t t\tr0" << std::endl;

        for (int i = 0; i < pref.size(); i++) {
            file << ((i + 1) * 10) <<  "\t" << b[i].halfLife << "\t" << r0 << std::endl;
        }

        file.close();
    }

    return 0;
}