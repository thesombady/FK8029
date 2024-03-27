#include <iostream>
#include <cmath>
#include <vector>
#include <string>


const float height = 10000; // m
const float sigmaV = 0.1; // placeholder
const float sigmaIf = 0.1; // placeholder

struct Cell {
    double energy;

    double tV[2];
    double tIf[2];

    double density;
    double height;

    Cell() { // Constructor
        this->energy = 0;
        this->density = 0;
        this->height = 0;

        this->tV[0] = 0.0; // In
        this->tV[1] = 0.0; // Out
        
        this->tIf[0] = 0.0; // In
        this->tIf[1] = 0.0; // Out
    };

    ~Cell() {}; // Destructor
};

/**
 * This function will create the cells and calculate the height and density of the cells
 * @param[in] numberofCells The number of cells to create
 * @returns A vector of cells
*/
std::vector<Cell> makeCells(int numberofCells) {
    double deltaHeight = height / ( numberofCells - 1 );
    std::vector<Cell> cells;
    for (int i = 0; i < numberofCells; i++) {
        // Here we want to create a new cell and actually calculate the height and density before we return it
        Cell cell = Cell();
        cell.height = i * deltaHeight;
        cells.push_back(cell);
    }
    // Calculate the density if the cells here.
    for (int i = 0; i < numberofCells; i++) {
        std::cout << cells[i].height << std::endl;
    }
    return cells;
}

/**
 * This function will compute the energy in the cells from above
 * @param[in] cells The cells to compute the energy in
 * @returns A vector of cells
*/
std::vector<Cell> computeDown(std::vector<Cell> cells) {
    // First we compute T
    double tV = 0; // We will overide this value in each of the loops but we one want initialize once.
    double tIf = 0; // We will overide this value in each of the loops but we one want initialize once.
    double energy = 0; // We will overide this value in each of the loops but we one want initialize once.

    for (int i = 0; i < cells.size(); i++) {
        if (i == 0) {
            // We are at the top cell

        } else if (i == cells.size() - 1) {
            // We are at the bottom cell
            tV = cells[i - 1].tV[1] * std::exp(-sigmaV * cells[i].height * cells[i].density);
            cells[i].tV[0] = tV;
            tIf = (cells[i - 1].tIf[1] + cells[i - 1].energy / 2 ) * std::exp(-sigmaIf * cells[i].height * cells[i].density);
            cells[i].tIf[0] = tIf;
            energy = ((cells[i - 1].energy + cells[i + 1].energy) / 2 + cells[i + 1].tV[1]); // Note done
            // Do we need in and out for both tV and tIf or or just T, tV and tIf?
        } else {
            // We are in the middle cells
        }
    }

    return cells;
}

/**
 * This function will compute the energy in the cells from below
 * @param[in] cells The cells to compute the energy in
 * @returns A vector of cells
*/
std::vector<Cell> computeUp(std::vector<Cell> cells) {

    return cells;
}

void calculate() {
    int numberOfCells = 10;
    float deltaHeight = height / numberOfCells;    

    std::vector<Cell> cells = makeCells(numberOfCells);
    // We want to go from above to below and then back and so on...   
    // We determine index 0 to be at the top of the atomosphere and index numberOfCells - 1 to be at the bottom

    int maxIterations, index;
    maxIterations = 1000;
    index = 1;

    while (index < maxIterations) {
        // Change to pointer and overide the content instead of returning a new vector?
        // Compute the energy in the cells from above;
        cells = computeDown(cells);
        
        // Compute the energy in the cells from below;
        cells = computeUp(cells);
        
        index++;
    }
    std::cout << "Done" << std::endl;

}



int main() {
    calculate();
    std::cout << "Hello, World!" << std::endl;
    return 0;
}