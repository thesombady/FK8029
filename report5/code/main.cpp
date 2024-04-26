#include <iostream>
#include <vector>
#include <string>
#include "Eigen/Dense"
#include "Eigen/src/Core/Matrix.h"
#include <fstream>


struct cSpline {

  static const int splineOrder = 4; // Spline order, we can change to 5, 6 and so on
  static const double charge = 1.0; // Unity from exercise
  static const double dx = 0.01;    // ONLY USED ONCE, in the testCase method

  int knotsNumber;

  int knotsPhyiscal;

  Eigen::VectorXd knots;  // The knots points

  Eigen::VectorXd splines;// The different splines

  Eigen::MatrixXd matrix; // Left hand side matrix.

  Eigen::VectorXd rhs;    // Right hand side of the matrix eq.

  Eigen::VectorXd coeff;  // Coeffecients in f(t_i) = sum_(n = i - k + 1)^i c_n B_i_k(t_i)

  /**
    Constructor: Creates the splines used
  */
  cSpline() {};

  // Deconstructor
  ~cSpline(){};

  /**
    setKnots: Manually set the knot-points
    @param[in] knots: an vector of the knot points
  */
  void setKnots(std::vector<double> knots) {
    int knotsPhysical = knots.size();

    // If we provide 4 physical knots, we will have 5 + 2 * 3 = 11 knots in total
    if (knotsPhysical < 5) {
      std::cout << "To few knot points, please include more" << std::endl;
      exit(1);
    }
    this -> knotsPhyiscal = knotsPhysical;
    this -> knotsNumber = knotsPhysical + 2 * (this -> splineOrder - 1);
    // If we have N knots, then knots is N, 1 vector / array
    this -> knots = Eigen::VectorXd::Zero(knotsNumber, 1);

    // The splines that we want to save
    this -> splines = Eigen::VectorXd::Zero(this -> knotsNumber - this -> splineOrder, 1);

    // Initialize the ghost points
    for (int i = 0; i < this -> splineOrder - 1; i++) {
      // Left side
      this -> knots(i) = knots[0];
      // Right side
      this -> knots(knotsPhyiscal + this -> splineOrder -1 + i) = knots[knotsPhysical - 1];
    }

    // Physical points
    for (int i = 0; i < this -> knotsPhyiscal; i++) {
      this -> knots(i + this -> splineOrder - 1) = knots[i];
    } 
   
    // We now need to compute the laplacian, which solves the problem
    getMatrix();

  };
  
  /**
    testCase: Visual representation of the splines
    @param[in] knotsNumber: The number of knots to use
    @param[in] filename: The file of where to save the spline data
  */
  void testCase(int knotsNumber, std::string filename) {
    this -> knotsNumber = knotsNumber;

    // dx: The increment we will use to find all the elements in the space

    int knotsPhyiscal = knotsNumber - 2 * (this -> splineOrder - 1); // Number of physical points,
    // the remaining points are ghost points
    int factor = 1 / this ->  dx;

    this -> knotsPhyiscal = knotsPhyiscal;
    // Points of intrests

    this -> knots = Eigen::VectorXd::Zero(knotsNumber, 1);
    // The knots, including ghost points

    this -> splines = Eigen::VectorXd::Zero(factor * knotsPhyiscal + 1, 1);
    // The splines we will fill

    // Initialize the ghost points
    for (int i = 0; i < this -> splineOrder - 1; i++) {
      // Left side
      this -> knots(i) = 0;//- this->splineOrder + (i + 1);
    }

    // Right side
    for (int i = 0; i < this -> splineOrder - 1; i++) {
        this -> knots(knotsPhyiscal + this -> splineOrder -1 + i) = knotsPhyiscal - 1;
    }

    // Physical points
    for (int i = 0; i < this -> knotsPhyiscal; i++) {
      this -> knots[i + 3] = i;
    }

    std::ofstream file(filename);

    file << "x\t";

    for (int i = 0; i < this -> knotsPhyiscal + 2; i++) {
      file << "y" + std::to_string(i) << "\t";
    }

    file << "tot" << std::endl;

    double x = 0;
    double data = 0;
    for (int i = 0; i < factor * knotsPhyiscal + 1; i++) {
      file << x << "\t";
      for (int j = 0; j < knotsPhyiscal + 2; j++){
        data = buildSpline(x, j, this -> splineOrder);
        this -> splines(i) += data;
        file << data << "\t";
      }
      file << this -> splines(i) << std::endl;
      x += this -> dx;
    }

    file.close();

  };

  /**
    secondDerivative: Computes the second derivative of Bspline_i^k,
    @param[in] i - The index
    @param[in] k - The spline order k
    @param[in] x - The position to evalute the derivative at
    @return Value of the second derivative at x
  */
  double secondDerivative(int i, int k, double x) {
    // Checked, formula is okey
    
    double sol = 0.0;

    Eigen::VectorXd t = this -> knots;

    double factor = (k - 1) * (k - 2);

    double cond1 = (t(i + k - 1) - t(i)) * (t(i + k - 2) - t(i));

    if (cond1 != 0.0){
      sol += factor / cond1 * buildSpline(x, i, k - 2);
    }

    double cond2 = (t(i + k - 1) - t(i)) * (t(i + k - 1) - t(i + 1));

    if (cond2 != 0.0){
      sol -= factor / cond2 * buildSpline(x, i + 1, k - 2);
    }

    double cond3 = (t(i + k) - t(i + 1)) * (t(i + k - 1) - t(i + 1));

    if (cond3 != 0.0){
      sol -=  factor / cond3 * buildSpline(x, i + 1, k - 2);
    }

    double cond4 = (t(i + k) - t(i + 1)) * (t(i + k) - t(i + 2));

    if (cond4 != 0.0){
      sol += factor / cond4 * buildSpline(x, i + 2, k - 2);
    }

    return sol;
  };

  /**
    firstDerivative: Computes the first derivative of Bspline_i^k,
    @param[in] i - The index
    @param[in] k - The spline order k
    @param[in] x - The position to evalute the derivative at
    @return Value of the first derivative at x
  */
  double firstDerivative(int i, int k, double x) {
    double sol = 0.0;
    double cond1 = this -> knots(i + k - 1) - this -> knots(i); 
    if (cond1 != 0) {
      sol += (k - 1) * buildSpline(x, i, k - 1) / (cond1);
    }

    double cond2 = this -> knots(i + k) - this -> knots(i + 1);
    if (cond2 != 0) {
      sol -= (k - 1) * buildSpline(x, i + 1, k) / (cond2);
    }
    return sol;
  };
  
  /**
    getMatrix: Computes the left hand side matrix of the equation
  */
  void getMatrix() {
  
    int k = this -> splineOrder;
    int N = this -> knotsNumber - k; // N - k equations
  
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(N, N);
    Eigen::MatrixXd rhs = Eigen::VectorXd::Zero(N, 1);

    matrix(0,0) = 1;
    matrix(N - 1, N - 1) = 1;

    double x;
    for (int i = 1; i < N - 1; i++) {
      x = this -> knots(i + this -> splineOrder);
      for (int j = i - k + 1; j <= i; j ++ ) {
        if (j > N || j < 0) {
          continue;
        }
        matrix(i, j + 1) = secondDerivative(j, k , x);
        rhs = uniformSphere(x);
      }
    }
    rhs(N - 1) = 1; // Boundary condition
    this -> matrix = matrix;
    this -> rhs = rhs;
  };

  /**
    solve: Solves the equation
    Solves the Poisson equation using LU factorization
    on the matrix obtained from getMatrix
  */
  void solve() {
    Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd> > Lu(this -> matrix);
    Eigen::VectorXd coeff;
    coeff = Lu.solve(this -> rhs);
    std::cout << "Coefficients" << std::endl;
    std::cout << coeff << std::endl;
    this -> coeff = coeff;
  };

  /**
    save: Saves the spline values to a file
    @param[in] filename: The file to save the spline values to
  */
  void save(std::string filename) {

    int factor = 100;

    this -> splines = Eigen::VectorXd(this -> knotsPhyiscal * factor + 1, 1);
    
    std::ofstream file(filename);

    Eigen::VectorXd t = this -> knots;
    double dx;
    int counter = 0;
    double y = 0;
    for (int i = 0; i < this -> knotsPhyiscal; i++) {
      dx = (t(i + this -> splineOrder -1) - t(i - this -> splineOrder + 1)) / factor;
      
      for (double x = t(i - this -> splineOrder + 1); x < t(i + this -> splineOrder - 1); x += dx) {
        buildSpline(x, counter, this -> splineOrder);
        
        counter++;
      }
    }

    file.close();
    
  }

  /**
    Build: Builds the values of the corresponding spline knots
    @param[in] x_pos: The current position of the i:th spline
    @param[in] i: The knot investigated
    @param[in] k: The k:th spline, starts form k then recursivly goes towards 1
    @return the value of the knot at knot 'knot'
  */
  double buildSpline(double x_pos, int i, int k) { // Recursive call for k

    Eigen::VectorXd t = this -> knots;
    if (k == 1) { // Fírst order spline
      if (t(i) <= x_pos && t(i + 1) > x_pos) { // Initial condition if k = 1, i.e. linear spline
        return 1.0;                         
      } else {
        return 0.0;
      }
    }

    double sol = 0.0;

    if (t(i + k - 1) != t(i)) { // Avoid division by zero
      sol += (x_pos - t(i)) / (t(i + k - 1) - t(i)) * buildSpline(x_pos, i, k - 1);
    }

    if (t(i + k) != t(i + 1)) { // Avoid division by zero
      sol += (t(i + k) - x_pos) / (t(i + k) - t(i + 1)) * buildSpline(x_pos, i + 1, k - 1);
    }

    return sol;  

  };

  /**
    uniformSphere: Return the charge density from a uniformly charged sphere
    @param[in] r: The position r
    @return The charge density at a point r
  */
  double uniformSphere(double r) {
    double R = 10.0;
    if (r <= R) {
      return this -> charge / (4.0 / 3.0 * 3.1415 * R * R * R);
    }
    return 0.0;
  };

  /**
    sphericalShell: Return the charge density from a spherically charged shell
    @param[in] r: The position r
    @return The charge density at a point r
  */
  double sphericalShell(double r) {
    double r1 = 5.0;
    double r2 = 10.0;
    if (r>= r1 && r <= r2) {
      return this -> charge / (4.0 / 3.0 * 3.1415 * (std::pow(r2, 3) - std::pow(r1, 3)));
    }
    return 0.0;
  };

};



int main() {
  //cSpline *test = new cSpline;

  //test -> testCase(11, "test.dat");
  
  std::vector<double> knots;

  knots.push_back(0.0);
  knots.push_back(1.0);
  knots.push_back(2.0);
  knots.push_back(3.0);
  knots.push_back(4.0);
  knots.push_back(5.0);
  knots.push_back(6.0);
  // knots.push_back(7.0);
  // knots.push_back(8.0);
  // knots.push_back(9.0);
  // knots.push_back(10.0);
  // knots.push_back(11.0);
  // knots.push_back(12.0);
  // knots.push_back(13.0);
  // knots.push_back(14.0);

  cSpline *knot = new cSpline;

  knot -> setKnots(knots);
  knot -> solve();

  //knot -> save("Value.dat");
  

  return 0;
}
