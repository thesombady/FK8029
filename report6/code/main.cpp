#include <cinttypes>
#include <ios>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include <vector>
#include <string>
#include <functional>
#include <tuple>
#include <fstream>

#include "../../report5/code/Eigen/Dense"
#include "../../report5/code/Eigen/Eigen"
#include "../../report5/code/Eigen/Eigenvalues"

#include "legendre.h" // https://github.com/haranjackson/LegendreGauss/tree/master
// The above include is not my own work, but the work of the author Hari Jackson, MIT Licence

const double charge = 1.0;
const double hbar = 1.0;
const double m = 1.0;
const double epsilon_0 = 1.0;

struct cSpline {


  int knotsNumber;

  int knotsPhyiscal;

  Eigen::VectorXd knots;  // The knots points

  std::vector<std::function<double(double)> > splines;// Splines[i](x)
  std::vector<std::function<double(double)> > dsplines;// d/dx Splines[i](x)
  std::vector<std::function<double(double)> > ddsplines;// d^2/dx^2 Splines[i](x)

  int splineOrder;

  /**
    Constructor: Creates the splines used
  */
  cSpline(int k = 4) {
    this -> splineOrder = k; // Spline order, we can change to 5, 6 and so on
  };

  // Deconstructor
  ~cSpline(){};

  void initialize(int physicalKnotsNumber) {
    this -> knotsPhyiscal = physicalKnotsNumber;
    this -> knotsNumber = physicalKnotsNumber + 2 * (this -> splineOrder - 1);

    Eigen::VectorXd knots = Eigen::VectorXd::Zero(knotsNumber, 1);

    for (int i = 0; i < this -> splineOrder - 1; i++) {
      // Left side
      knots(i) = 0;
      // Right side
      knots(knotsPhyiscal + this -> splineOrder -1 + i) = physicalKnotsNumber - 1;
    }

    // Physical points
    for (int i = 0; i < this -> knotsPhyiscal; i++) {
      knots(i + this -> splineOrder - 1) = i;
      if (i == this -> knotsNumber - 1) {

      }
    }

    // Prepare the splines and the first derivatives
    this -> createSplineFunctions();

    this -> knots = knots;

  };

  /**
    collmat(func, order): Computes the collocation matrix of the ghost physical points
    @param[in] func: The value on the right hand side,
    @return {collmat, rhs}
  */
  std::tuple<Eigen::MatrixXd, Eigen::VectorXd> collMat(std::function<double(double)> func) {
    // N - k splines, but we set the spline to be zero, since  we set the first value to be zero 
    int n = this -> knotsNumber - this -> splineOrder - 1;
    int k = this -> splineOrder;

    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n,n);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);

    
    std::vector<std::function<double(double)> > spl = this -> ddsplines;
    
    Eigen::VectorXd t = this -> knots;
    double shift = 1e-12;
    
    for (int i = k - 1; i < n + k - 1; i++) {
      if (i == n + 2) {
        for (int j = 1; j < spl.size(); j++) {
          mat(i - (k + 1), j - 1) = spl[j](t(i) - shift);
        }
      } else {
        for (int j = 0; j < n; j++) {
          mat(i - (k - 1), j) = spl[j + 1](t(i));
        }
      }
      rhs(i - (k - 1)) = func(t(i + 1));
    }
    // Right boundary
    mat(n - 1, n - 1) = this -> dsplines[spl.size() - 1](t(t.size() - 1) - shift);
    mat(n - 1, n - 2) = this -> dsplines[spl.size() - 2](t(t.size() - 1) - shift);
    rhs(rhs.size() - 1) = 0; // last element

    std::cout << mat << std::endl;


    return {mat, rhs};
    
  };

  /**
    saveSplines: Save the splines in the domain
    @param[in] order: The derivative order
  */
  void saveSplines(int order = 0) {
    Eigen::VectorXd t = this -> knots;
    double xmin, xmax;
    xmin = t(0);
    xmax = t(t.size() - 1);

    std::ofstream file("spline" + std::to_string(order) + ".dat");
    file << "x\t";

    std::vector<std::function<double(double)> > spl;
    if (order == 0){
      spl = this -> splines; // Splines[i](x)
    } else if (order == 1) {
      spl = this -> dsplines; // dSplines[i](x)
    } else if (order == 2) {
      spl = this -> ddsplines;
    } else {
      return; // We should throw an error
    }
    for (int i = 0; i < spl.size(); i++) {
      file << "B" + std::to_string(i);
      if (i < spl.size() - 1) {
        file << "\t";
      } else {
        file << "\n";
      }
    }

    double dx = 0.01;
    for (double x = xmin; x < xmax; x += dx){
         file << x << "\t";
      for (int i = 0; i < spl.size(); i++) {
        file << spl[i](x);
         
        if (i < spl.size() - 1) {
          file << "\t";
        } else {
          file << "\n";
        }
      }
    }

    file.close();
  };
  /**
    createSplineFunctions: Creates the spline functions
    Creates the splines in an array of n-k length
    Creates the first derivative splines in an array of n-k length
    Creates the second derivative splines in an array of n-k length
  */
  void createSplineFunctions() {
    int k = this -> splineOrder;
    std::vector<std::function<double(double)> > splines;
    std::vector<std::function<double(double)> > dspline;
    std::vector<std::function<double(double)> > ddspline;
    for (int i = 0; i < this -> knotsNumber - this -> splineOrder; i++) { //n - k splines
      splines.push_back([i, k, this](double x) -> double {return this -> buildSpline(x, i, k);});
      dspline.push_back([i, k, this](double x) -> double {return this -> firstDerivative(i, k, x);});
      ddspline.push_back([i, k, this](double x) -> double {return this -> secondDerivative(i, k, x);});
    }
    this -> splines = splines;
    this -> dsplines = dspline;
    this -> ddsplines = ddspline;

  }

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

    int k = this -> splineOrder;
    this -> createSplineFunctions(); 
   
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
      sol -= (k - 1) * buildSpline(x, i + 1, k - 1) / (cond2);
    }
    return sol;
  };

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
    gaussianQuad: Computes the integral of a function using Gaussian quadrature
    @param[in] a: The lower limit
    @param[in] b: The upper limit
    @param[in] f: The function to integrate
    @returns The integral of the function

    Note: Implemented from numerical recipes
  */
  double gaussianQuad(float a, float b, std::function<double(double)> f) {
    
    int n = this -> splineOrder;

    Eigen::MatrixXd gauss = leggauss(n);

    Eigen::VectorXd roots = gauss.row(0);

    Eigen::VectorXd weights = gauss.row(1);

    double xm = 1.0 / 2.0 * (b + a);

    double xr = 1.0 / 2.0 * (b - a);

    double s = 0.0;

    for (int i = 0; i < n; i++) {
      s += weights(i) * f(xr * roots(i) + xm);
    }

    s *= xr;

    return s;
  };

  /**
    getB: Computes the matrix B for the linear system of equations
    @param[in] func: additional function to compute
    @returns The matrix B
  */
  Eigen::MatrixXd getB(std::function<double(double)> func = [](double r) -> double {return 1.0;}) {

    int n = this -> knotsNumber - this -> splineOrder - 2; // -2 from bc

    // We will Have  n - k splines,
    // we disregard the first and last spline due to bc

    int k = this -> splineOrder;

    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, n);

    std::vector<std::function<double(double)>> splines = this -> splines; // All splines as function of x
    // thus splines[index](x) = B_i^k(x)

    double ti;
    double ti1;
    for (int i = 1; i < n + 1; i++) {
      for (int j = 1; j < n + 1; j++) {
        // The function to integrate, which is B_i^k(x) * f(x) * B_j^k(x), where f(x) is either a constant
        // or the potential term on the lhs of the equation
        // We get spline [i + 1] due to the fact that we disregard the first spline, and the last spline
        int min = std::max(i, j);
        int max = std::min(i, j) + k - 1;
        double element = 0.0;
        for (int m = min; m <= max; m++) {

          std::function<double(double)> f = [i, j, splines, func](double x) -> double {return splines[i](x) * func(x) * splines[j](x);}; // To skip the first and last spline
          ti = this -> knots(m);
          ti1 = this -> knots(m + 1);

          element += this -> gaussianQuad(ti, ti1, f);

        }
        mat(i - 1, j - 1) = element;
      }
    }
    return mat;
  };


  /**
    V: Computes the electron repulsion and angular part of the hamiltonian
    @param[in] l: The angular momentum
    @param[in] z: proton number
    @param[in] r: The position
    @returns the potential term
  */
  static double V(int l, int z, double r) {
    if (r == 0) {
      r += 1e-12;
    }
    return + l * (l + 1) / (2 * m * r * r) - z * charge * charge / (hbar * hbar * r);
  };


  /**
    V: Computes the electron repulsion of a uniform sphere and angular mometum part of the hamiltonian
    @param[in] l: The angular momentum
    @param[in] z: proton number
    @param[in] r: The position
    @returns the potential term
  */
  static double VcV(int l, int z, double r) {
  if (r == 0) {
    r+= 1e-12;
  }
  double sphereRadii = 1.2 * std::pow(z, 1.0/3.0);
  
  double v = l * (l + 1) / (2 * m * r * r);
  if (r < sphereRadii) {
    v += -z * charge * charge / (hbar * hbar * r); // * 4 * pi * epsilon_0 
  } else {
    v+= - z * charge * charge / (2 * sphereRadii * hbar * hbar) * (3 - std::pow(r / sphereRadii, 2.0));
  }
  return v;
   };
  /**
    getH1: Computes the matrix H for db_i * db_j for the linear system of equations
    @returns The matrix H1
  */
  Eigen::MatrixXd getH1() {

    
    int n = this -> knotsNumber - this -> splineOrder - 2; // -2 from bc

    // We will Have  n - k splines,
    // we disregard the first and last spline due to bc

    int k = this -> splineOrder;

    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, n);

    std::vector<std::function<double(double)>> dsplines = this -> dsplines; // All splines as function of x 
    // thus dsplines[index](x) = B_i^k(x)

    double ti;
    double ti1;
    for (int i = 1; i < n + 1; i++) {
      for (int j = 1; j < n + 1; j++) {
        // The function to integrate, which is B_i^k(x) * f(x) * B_j^k(x), where f(x) is either a constant
        // or the potential term on the lhs of the equation
        // We get spline [i + 1] due to the fact that we disregard the first spline, and the last spline
        int min = std::max(i, j);
        int max = std::min(i, j) + k - 1;
        double element = 0.0;
        for (int m = min; m <= max; m++) {

          std::function<double(double)> f = [i, j, dsplines](double x) -> double {return dsplines[i](x) * dsplines[j](x);}; // To skip the first and last spline
          ti = this -> knots(m);
          ti1 = this -> knots(m + 1);

          element += this -> gaussianQuad(ti, ti1, f);

        }
        mat(i - 1, j - 1) = element;
      }
    }
    return mat * 1/2;
  };

  /**
    getH: Computes the matrix H for the linear system of equations
    @param[in] l: The angular momentum
    @param[in] z: The charge
    @returns The matrix H
  */
  Eigen::MatrixXd getH(int l, int z, bool cV) {
    Eigen::MatrixXd H = this -> getH1();
    
    if (cV == false) {
      H += this -> getB([l, z](double x) -> double {return cSpline::V(l, z, x);});
    } else {
      H += this -> getB([l, z](double x) -> double {return cSpline::VcV(l, z, x);});
    }

    return H;
  };

  /**
    save: Saves the function P_{nl}(r) and R_{nl}(r)
    @param[in] filename: The filename to save to
    @param[in] ges: The general eigen value solver, with the computed eigen-pair
  */
  void save(std::string filename, Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges) {
    filename += ".dat";

    std::ofstream file(filename);

    Eigen::VectorXd t = this -> knots;
    std::vector<std::function<double(double)>> splines = this -> splines;
    Eigen::MatrixXd mat = ges.eigenvectors();

    double r = 0.00001;
    double dr = 0.1;
    file << "r\t1s\t2s\t2p\t3s\t3p\t3d\n";

    double p;
    double p1;
    double p2;
    double R;
    double R1;
    double R2;


    while (r < this -> knots(this -> knots.size() - 1)) {
      p = 0;
      p1 = 0;
      p2 = 0;
      R = 0;
      R1 = 0;
      R2 = 0;

      for (int i = 1; i < splines.size() - 1; i++) {
        p += splines[i](r) * mat(i - 1, 0);
        p1 += splines[i](r) * mat(i - 1, 1);
        p2 += splines[i](r) * mat(i - 1, 2);

        R += splines[i](r) * mat(i - 1, 0) / r;
        R1 += splines[i](r) * mat(i - 1, 1) / r;
        R2 += splines[i](r) * mat(i - 1, 2) / r;
      }
      file << r << "\t" << p << "\t" << p1 << "\t" << p2 << "\t" << R << "\t" << R1 << "\t" << R2 << "\n";
      r += dr;
    }

    std::cout << ges.eigenvalues() << std::endl;
    

    file.close();

  }

  /**
    solve: Solves the linear system of equation
    @param[in] l: angular momentum 
    @param[in] z: Particle 
  */
  void solve(int l = 0, int z = 1.0, bool cV = false) {
    Eigen::MatrixXd B = this -> getB();
    Eigen::MatrixXd H = this -> getH(l, z, cV);
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges;
    
    // std::cout << "B: \n";
    // std::cout << B << std::endl;
    // std::cout << "H: \n";
    // std::cout << H << std::endl; 

    ges.compute(H, B);
    
    std::string filename = "l" + std::to_string(l) + "z" + std::to_string(z);
    if (cV) {
      filename += "VcV";
    }
    this -> save(filename, ges);
    
  };

};


int main() {
  cSpline *test = new cSpline;
  std::vector<double> knots;
  ///*
  knots.push_back(0.0);
  knots.push_back(.5);
  knots.push_back(1.0);
  knots.push_back(1.5);
  knots.push_back(2.0);
  knots.push_back(2.5);
  knots.push_back(3.5);
  knots.push_back(4.0);
  knots.push_back(4.5);
  knots.push_back(5.0);
  knots.push_back(5.5);
  knots.push_back(6.0);
  knots.push_back(6.5);
  knots.push_back(7.0);
  knots.push_back(7.5);
  knots.push_back(8.0);
  knots.push_back(8.5);
  knots.push_back(9.0);
  knots.push_back(9.5);
  knots.push_back(10.0);
  knots.push_back(10.5);
  knots.push_back(11.0);
  knots.push_back(11.5);
  knots.push_back(12.0);
  knots.push_back(12.5);
  knots.push_back(13.0);
  knots.push_back(13.5);
  knots.push_back(14.0);
  knots.push_back(14.5);
  knots.push_back(15.0);
  knots.push_back(15.5);
  knots.push_back(16.0);
  knots.push_back(16.5);
  knots.push_back(17.0);
  knots.push_back(17.5);
  knots.push_back(18.0);
  knots.push_back(18.5);
  knots.push_back(19.0);
  knots.push_back(19.5);
  knots.push_back(20.0);
  knots.push_back(20.5);
  knots.push_back(21.0);
  knots.push_back(21.5);
  knots.push_back(22.0);
  knots.push_back(22.5);
  knots.push_back(23.0);
  knots.push_back(23.5);
  knots.push_back(24.0);
  knots.push_back(24.5);
  knots.push_back(25.0);
  //*/
  int numberOfPoints = 50;
  for (int i = 0; i <= numberOfPoints; i++) {
    //knots.push_back(double(i) / (double(numberOfPoints)/2.0));
  }

  //test -> initialize(11);
  test -> setKnots(knots);
  int l, z;
  z = 1;
  std::cout << "Input l: ";
  std::cin >> l;
  std::cout << std::endl;
  //test -> saveSplines(0); // They look fine
  //test -> saveSplines(1); // They look fine
  //test -> saveSplines(2); // They look fine
  //test -> collMat([](double x) -> double {return 1.0;}); // collMat now works fine
  //std::cout << test -> gaussianQuad(0, 1, [](double x) -> double {return x;}) << std::endl;
 
  test -> solve(l, z, false);
  

  return 0;
}
