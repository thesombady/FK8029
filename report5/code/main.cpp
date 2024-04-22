#include <iostream>
#include <vector>
#include <string>
#include "Eigen/Dense"
#include "Eigen/src/Core/Matrix.h"
#include <fstream>


struct cSpline { // TODO: Make such that i can have inheritance?
  static const int splineOrder = 4;
  static const double charge = 10;


  int knotsNumber; // Argument for number of knots in constructor

  int knotsPhyiscal; // Calculated from the number of knots and the spline order

  Eigen::VectorXd knots; // The knots indicies
  Eigen::VectorXd splines; // The spline values
  Eigen::MatrixXd matrix; // Matrix constructed form
  

  /**
    Constructor: Creates the splines used
    @param[in] knotsNumber: The number of knots to use
    @param[in] filename: The file of where to save the spline data
  */
  cSpline(int knotsNumber, std::string filename) {
    this -> knotsNumber = knotsNumber;

    double dx = 0.01;
    // dx: The increment we will use to find all the elements in the space

    int knotsPhyiscal = knotsNumber - 2 * (this -> splineOrder - 1); // Number of physical points,
    // the remaining points are ghost points
    int factor = 1 / dx;

    this -> knotsPhyiscal = knotsPhyiscal;
    // Points of intrests

    this -> knots = Eigen::VectorXd::Zero(knotsNumber, 1);
    // The knots, including ghost points

    this -> splines = Eigen::VectorXd::Zero(factor * knotsPhyiscal + 1, 1);
    // The splines we will fill

    this -> matrix = Eigen::MatrixXd::Zero(knotsPhyiscal, knotsPhyiscal); 
    // The matrix of coefficients, square matrix. We will fill this later on

    // Initialize the ghost points
    for (int i = 0; i < this -> splineOrder; i++) {
      // Left side
      this -> knots(i) = 0;//- this->splineOrder + (i + 1);
    }

    // Right side
    for (int i = 0; i < this -> splineOrder - 1; i++) {
        this -> knots(knotsPhyiscal + this -> splineOrder -1 + i) = knotsPhyiscal - 1;
    }

    std::ofstream file(filename);

    file << "x\t";


    for (int i = 0; i < this -> knotsPhyiscal; i++) {
      this -> knots[i + 3] = i;
    }

    std::cout << this -> knots << std::endl;

    for (int i = 0; i < this -> knotsPhyiscal + 2; i++) {
      file << "y" + std::to_string(i) << "\t";
    }

    file << "tot" << std::endl;

    double x = 0;
    double data = 0;
    for (int i = 0; i < factor * knotsPhyiscal + 1; i++) {
      file << x << "\t";
      for (int j = 0; j < knotsPhyiscal + 2; j++){
        data = build(x, j, this -> splineOrder, this -> knots);
        this -> splines(i) += data;
        file << data << "\t";
      }
      file << this -> splines(i) << std::endl;
      x += dx;
    }

    file.close();

  };

  /**
    Build: Builds the values of the corresponding spline knots
    @param[in] x_pos: The current position of the i:th spline
    @param[in] knot: The knot investigated
    @param[in] k: The k:th spline, starts at 4, then recursivly goes towards 1
    @param[in] knots: All the knots that we have, including the ghost points
    @return the value of the knot at knot 'knot'
  */
  double build(double x_pos, int knot, int k, Eigen::VectorXd knots) { // Recursive call for k

    if (k == 1) { // Fírst order spline
      if (knots(knot) <= x_pos && knots(knot + 1) > x_pos) { // Initial condition if k = 1, i.e. linear spline
        return 1.0;                         
      } else {
        return 0.0;
      }
    }

    double sol = 0.0;

    if (knots(knot + k - 1) != knots(knot)) { // Avoid division by zero
      sol += (x_pos - knots(knot)) / (knots(knot + k - 1) - knots(knot)) * build(x_pos, knot, k - 1, knots);
    }

    if (knots(knot + k) != knots[knot + 1]) { // Avoid division by zero
      sol += (knots(knot + k) - x_pos) / (knots(knot + k) - knots(knot + 1)) * build(x_pos, knot + 1, k - 1, knots);
    }

    return sol;
  
  };

  /**
    uniformSphere: Return the charge density from a uniformly charged sphere
    @param[in] r: The position r
    @param[in] R: The radii of the sphere
    @return The charge density at a point r
  */
  double uniformSphere(double r, double R) {
    if (r <= R) {
      return this -> charge / (4.0 / 3.0 * 3.1415 * R * R * R);
    }
    return 0.0;
  };

  /**
    sphericalShell: Return the charge density from a spherically charged shell
    @param[in] r: The position r
    @param[in] r1: The inner radii
    @param[in] r2: The outer radii
    @return The charge density at a point r
  */
  double sphericalShell(double r, double r1, double r2) {
    if (r>= r1 && r <= r2) {
      return this -> charge / (4.0 / 3.0 * 3.1415 * r2 * r2 * r2);
    }
    return 0;
  };

};



int main() {

  cSpline spline(11, "test.dat");

  return 0;
}
