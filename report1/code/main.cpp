// Your First C++ Program

#include <iostream>
#include <string>
#include <vector>

const float radii = 1.0;

float potentialEnergy(float r) {
    if (r <= radii) {
        return 0.0;
    }
    return 1.0 / r;
}
/*
    We want to discretize our potential energy to small widts
*/

int main() {
    int numberOfBins = 1;

    float disc [numberOfBins][numberOfBins];
    
    for (int i = 0; i < numberOfBins; i++) {
        for (int j = 0; j < numberOfBins; j++) {
            disc[i][j] = potentialEnergy(3.0);
        }
    }

    std::cout << *disc[0];


    return 0;
}

