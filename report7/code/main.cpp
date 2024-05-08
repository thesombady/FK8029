#include <cinttypes> // int a; cin >> a works
#include <cmath> // map
#include <ios> // os
#include <iostream> // io
#include <stdexcept> // error
#include <system_error> // error
#include <vector> // on call heap array
#include <string> // string
#include <functional> // function type
#include <tuple> // tuple type
#include <fstream> // write and read to file

#include "../../report5/code/Eigen/Dense"
#include "../../report5/code/Eigen/Eigen"
#include "../../report5/code/Eigen/Eigenvalues"

#include "legendre.h" // https://github.com/haranjackson/LegendreGauss/tree/master
// The above include is not my own work, but the work of the author Hari Jackson, MIT Licence


typedef Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> Solver;
typedef Eigen::EigenSolver<Eigen::MatrixXd> Solver2;
typedef Eigen::MatrixXd Mat;
typedef std::function<double(double)> lambda;

const double charge = 1.0;
const double hbar = 1.0;
const double m = 1.0;
const double PI = 3.14159265358979323846;
const double a0 = 1.0;
const double eps0 = 1/ (4.0 * PI);


double uniformSphere(double r) {
  if (r > 10) {
    return 0.0; 
  } else {
    double z = 1;
    return charge * z / (4 * PI / 3 * (10.0*10.0*10.0));
  }
}

double shell(double r) {
  if (r < 5) {
    return 0.0;
  } else if (r > 10) {
    return 0.0;
  } else {
    double z = 1;
    return charge * z / (4*PI/3*(10.0 * 10.0 * 10.0 - 5.0 * 5.0 * 5.0));
  }
}


double exactShell(double r) {
  // if (r < 5 || r > 10) {
  //   return 0.0;
  // }
  if (r == 0.0) {
    r += 1e-12;
  }
  double R3 = std::pow(10.0, 3.0);
  double r3 = std::pow(5.0, 3.0);
  double volume = 4.0 * PI * (R3 - r3);
  return charge / ( eps0 * volume ) * (std::pow(10.0, 2.0) / 2.0 - 1.0 / 3.0*(std::pow(5.0, 3.0) / r + std::pow(r, 2.0)/2.0));
}

double wf(double r) {
  return 1.0 / (PI * a0) * std::exp(-2.0*r / a0);
}

double wfe(double r) {
  if (r == 0.0) {
    r += 1e-12;
  }
  return charge / (4*PI*eps0) * (1.0 / r - std::exp(-2.0*r)*(1.0/r + 1.0)); 
}

struct cSpline {

  int knotsNumber;

  int knotsPhyiscal;

  Eigen::VectorXd knots;  // The knots points

  std::vector<std::function<double(double)> > splines;// Splines[i](x)
  std::vector<std::function<double(double)> > dsplines;// d/dx Splines[i](x)
  std::vector<std::function<double(double)> > ddsplines;// d^2/dx^2 Splines[i](x)

  Eigen::VectorXd coeff;

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
    @param[in] z: The electron number
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
    
    double element;
    for (int i = k - 1; i < n + k - 1; i++) {
      if (i == n + 2) { // "Last" row is special
        for (int j = 1; j < spl.size(); j++) {
          element = spl[j](t(i) - shift);
          if (std::abs(element) >= 1e-10) {// Due to the shift we actually have 4 points,
            mat(i - (k), j - 1) = element; // We simply remove it
          }
        }
      }

      for (int j = 0; j < n; j++) {
        mat(i - (k - 1), j) = spl[j + 1](t(i));
      }     
    }

    for (int j = 1; j < n; j++) {
      rhs(j) = -func(t(j + k - 1)) * t(j + k - 1) * 4.0 * PI / (4.0 * PI * eps0);
    }
    // Right boundary
    mat(n - 1, n - 1) = this -> dsplines[spl.size() - 1](t(t.size() - 1) - shift);
    mat(n - 1, n - 2) = this -> dsplines[spl.size() - 2](t(t.size() - 1) - shift);
    rhs(rhs.size() - 1) = 0; // last element

    //std::cout << mat << std::endl;
    return std::make_tuple(mat, rhs);
    
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
    @param[in] n: Order of legendre polynomial
    @returns The integral of the function
    Note: Implemented from numerical recipes
  */
  static double gaussianQuad(float a, float b, std::function<double(double)> f, int n) {

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
    @param[in] func: additional function to compute, used for part of the lhs.
    @returns The matrix B
  */
  Eigen::MatrixXd getB(std::function<double(double)> func = [](double r) -> double {return 1.0;}) {
    int n = this -> knotsNumber - this -> splineOrder - 2; // -2 from bc

    // We will Have  n - k splines,
    // we disregard the first and last spline due to bc
    int k = this -> splineOrder;

    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, n);
    std::vector<std::function<double(double)>> splines = this -> splines; // B_i^k(x) 

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

          element += cSpline::gaussianQuad(ti, ti1, f, this -> splineOrder);
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
    @param[in] vEE: The additional term
    @returns the potential term
  */
  static double V(int l, int z, double r, std::function<double(double)> vEE) {
    if (r == 0) {
      r += 1e-12;
    }
    return hbar * hbar * l * (l + 1) / (2 * m * r * r) - z * charge * charge / (4.0 * PI * eps0 * r) + charge * vEE(r);
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
    std::vector<std::function<double(double)>> dsplines = this -> dsplines; // d B_i^k(x) / dx

    double ti;
    double ti1;
    for (int i = 1; i < n + 1; i++) {
      for (int j = 1; j < n + 1; j++) {
        // The function to integrate, which is dB_i^k(x) * dB_j^k(x)
        // We get spline [i + 1] due to the fact that we disregard the first spline, and the last spline
        int min = std::max(i, j);
        int max = std::min(i, j) + k - 1;
        double element = 0.0;
        for (int m = min; m <= max; m++) {

          std::function<double(double)> f = [i, j, dsplines](double x) -> double {return dsplines[i](x) * dsplines[j](x);}; // To skip the first and last spline
          ti = this -> knots(m);
          ti1 = this -> knots(m + 1);

          element += cSpline::gaussianQuad(ti, ti1, f, this -> splineOrder);

        }
        mat(i - 1, j - 1) = element;
      }
    }
    return mat * hbar * hbar / (2 * m);
  };

  /**
    getH: Computes the matrix H for the linear system of equations
    @param[in] l: The angular momentum
    @param[in] z: The charge
    @param[in] vEE: The additional term
    @returns The matrix H
  */
  Eigen::MatrixXd getH(int l, int z, std::function<double(double)> vEE) {
    Eigen::MatrixXd H = ( this -> getH1() );
    //H += this -> getB([l, z, vEE](double x) -> double {return cSpline::V(l, z, x, vEE);});
    return H + this -> getB([l, z, vEE](double x) -> double {return cSpline::V(l, z, x, vEE);}) ;
  };

  /**
    saveWave: Saves the function P_{nl}(r) and R_{nl}(r)
    @param[in] filename: The filename to save to
    @param[in] ges: The general eigen value solver, with the computed eigen-pair
  */
  void saveWave(std::string filename, Solver ges) {
    std::string path = "wave" + filename + ".dat";

    std::ofstream file(path);

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
      file << r << "\t" << p << "\t" << p1 << "\t" << p2 << "\t";
      file << R << "\t" << R1 << "\t" << R2 << "\n";
      r += dr;
    }

    std::cout << ges.eigenvalues() << std::endl;

    file.close();
  };

  /**
    solveColl: Solves the collocation problem
    @param[in] f - The charge density function
    @param[in] z - The number of electrons
    @returns func - The solution V(r) to the collocation problem
  */
  std::function<double(double)> solveColl(std::function<double(double)> f){
    Eigen::MatrixXd mat;
    Eigen::VectorXd res;
    Eigen::VectorXd coeff;
    std::tie(mat, res) = this -> collMat(f);
    coeff = mat.partialPivLu().solve(res);
    //this -> coeff = coeff;
    std::vector<std::function<double(double)>> splines = this -> splines;

    // Test to plot the collocation problem, to view that it works.
    // double r = 0.0;
    // double dr = 0.1;
    // std::cout << "Spline len: " << this -> splines.size() << std::endl;
    // std::ofstream file("test.dat");
    // while (r < 15) {
    //   file << r << "\t" << getV(r) << std::endl;

    //   r += dr;
    // }
    // file.close();

    // Returns V = P_{nl}(r) / r
    return [coeff, splines](double r) -> double {
      // Optimize by finding the nearest spline and just
      // evaluate close to it instead of evaluating all of them
      double y = 0;
      if (r == 0) {
        r +=1e-12;
      }
      for (int i = 0; i < coeff.size(); i++) {
        y+= coeff[i] * splines[i + 1](r) / r;
      }
      return y;
    };
    
  };

  // Remove, moved into the file, the return of solveColl
  double getV(Eigen::VectorXd coeff, std::vector<std::function<double(double)>> splines, double r) {
    double y = 0;
    if (r == 0) {
      r += 1e-12;
    }
    for (int i = 0; i < coeff.size(); i++) {
      //  V = phi(r) / r = sum_i b_i^k(r) / r
      y +=  coeff[i] *  splines[i + 1](r) / r;
    }

    return y;
  };  

  /**
    solve: Solves the linear system of equation
    @param[in] l: angular momentum 
    @param[in] z: Particle 
    @param[in] vEE: The additiuonal function, vEE(x) = 0 by default
    @returns func: Returns the function R(r, j) where r is the position, and j denotes the level n
  */
  std::function<double(double, int)> solveAtom(int l = 0, int z = 1.0, lambda vEE = [](double x) -> double {return 0;}) {
    Eigen::MatrixXd B = this -> getB();
    Eigen::MatrixXd H = this -> getH(l, z, vEE);
    Solver ges;
    
    // std::cout << "B: \n";
    // std::cout << B << std::endl;
    // std::cout << "H: \n";
    // std::cout << H << std::endl; 

    ges.compute(H, B);
    // this -> saveWave("test3", ges);

    std::vector<std::function<double(double)>> splines = this -> splines;
    Eigen::MatrixXd mat = ges.eigenvectors();
    // We want to solve the thing in the end
    // Returns R(r)
    return [mat, splines](double r, int j) -> double {
      double y = 0;
      if (r == 0) {
        r += 1e-12;
      }
      for (int i = 1; i < splines.size() - 1; i++) {
        y += mat(i - 1, j) * splines[i](r) / r; // phi(r) / ( r ) = c_i * B_i^k(r) / r
      }
      return y;
    };
  };
  
  /**
    solveAtomPref: Solves the atom with the additional potential
    @param[in] ges: The general eigen value solver
    @param[in] B: The matrix B
    @param[in] H1: The matrix H1, the constant part
    @param[in] l: The angular momentum
    @param[in] z: The number of protons
    @param[in] vEE: The additional potential
    @returns The eigen vectors
  */
  Mat solveAtomPref(Solver ges, Mat B, Mat H1, int l = 0, double z = 1.0, lambda vEE = [](double x) -> double {return 0;}) {
    // make tuple and return ges.eigenvalues()
    Mat H = H1 + this -> getB([l, z, vEE](double x) -> double {return cSpline::V(l, z, x, vEE);}); // The term that changes
    ges.compute(H, B);
    
    std::vector<std::function<double(double)>> splines = this -> splines;
    Mat basis = ges.eigenvectors();

    return basis; 
  };
  
  Mat solveAtomPref2(Mat Binv, Mat H1, int l = 0, double z = 1.0, lambda vEE = [](double x) -> double {return 0;}) {
    // make tuple and return ges.eigenvalues()
    Mat H = H1 + this -> getB([l, z, vEE](double x) -> double {return cSpline::V(l, z, x, vEE);}); // The term that changes
    Solver2 ges(Binv * H)
    
    std::vector<std::function<double(double)>> splines = this -> splines;
    Mat basis = ges.eigenvectors();

    return basis; 
  };

  std::function<double(double)> getP(int k, Mat basis, Mat B) {
    std::vector<std::function<double(double)>> splines = this -> splines;

    double norm = basis.col(k).transpose() * B * basis.col(k); // Normalize
    std::cout << "Inner norm: " << norm << std::endl;

    // Can we do like this?
    // return basis.col(k) * splines(r)
    return [basis, splines, k, norm](double r) -> double {
      double y = 0;
      for (int j = 1; j < splines.size() - 1; j++) {
        y += basis(j - 1, k) * splines[j](r);
      }
      return y / std::sqrt(norm);
    };
  };
};

/**
  Simple Riemann integration
*/
double Riemann(double a, double b, std::function<double(double)> func){
  double sum = 0;
  double dx = 0.05;
  while (a < b) {
    sum += func(a) * dx;

    a+= dx;
  }

  return sum;
}

/**
  Simple Trapezoidal method for integration
*/
double trapz(double a, double b, std::function<double(double)> func) {
  double dx = 0.05;
  double sum = 2 * func(a);
  while (a < b) {
    sum += 2 * func(a);
    a += dx;
  }
  sum += func(b);
  return sum * dx / 2.0;
  
}


std::function<double(double)> normalize(double a, double b, std::function<double(double)> func) {
  double norm = cSpline::gaussianQuad(a, b, [func](double x) -> double {return func(x) * func(x);}, 8);
  std::cout << "Norm: " << std::sqrt(norm) << std::endl;
  return [func, norm](double x) -> double {return func(x) / std::sqrt(norm);};
}

void test() {
  double Z = 1.0;

  std::vector<double> t;

  int xmax = 15;
  double dx = 0.2;
  double x = 0.0;
  int i = 0;
  while (x <= xmax) {
    t.push_back(x);
    //std::cout << knots[i] << std::endl;
    x += dx;
    i++;
  }

  cSpline *solver = new cSpline;
  solver -> setKnots(t);
  Mat B = solver -> getB();

  Mat H1 = solver -> getH1();


  // We first solve the atom without any additional potential, i.e. vEE(r) = 0;
  std::function<double(double)> p_nl; // int -> the state 1s, 2s, 3s and so on

  Mat basis;

  basis = solver -> solveAtomPref2(B, H1, 0, Z); // l = 0, Z = 2 for helium
  p_nl = solver -> getP(0, basis, B); // 1s state

  std::ofstream file("refactor.dat");
  double r = 0;
  while (r < xmax) {
    file << r << "\t" << p_nl(r) << std::endl;

    r += 0.1;
  }
  file.close();

}


// Solve the helium
void Helium(){

  double Z = 2.0;

  std::vector<double> t;

  int xmax = 15;
  double dx = 0.2;
  double x = 0.0;
  int i = 0;
  while (x <= xmax) {
    t.push_back(x);
    //std::cout << knots[i] << std::endl;
    x += dx;
    i++;
  }

  cSpline *solver = new cSpline;
  solver -> setKnots(t);
  Mat B = solver -> getB();

  Mat H1 = solver -> getH1();
  Solver ges;

  // We first solve the atom without any additional potential, i.e. vEE(r) = 0;
  std::function<double(double)> p_nl; // int -> the state 1s, 2s, 3s and so on
  Mat basis;

  basis = solver -> solveAtomPref(ges, B, H1, 0, Z); // l = 0, Z = 2 for helium
  p_nl = solver -> getP(0, basis, B); // 1s state

  //normalize(0, xmax, p_nl);

  std::function<double(double)> rho;

  rho = [p_nl](double r) -> double {
    if (r == 0) {
      r+= 1e-12;
    }
    return charge / (4 * PI) * 2 * (std::pow(p_nl(r) / r, 2.0)); // 2 comes form N_j, only the 1s state occupied.
  };

  std::cout << "Trapz: " << 4.0 * PI * trapz(0, xmax, [rho](double r) -> double {return rho(r) * r * r;}) << std::endl;
  std::cout << "Riemann: " << 4.0 * PI * Riemann(0, xmax, [rho](double r) -> double {return rho(r) * r * r;}) << std::endl;
  std::cout << "Quad: " << 4.0 * PI * cSpline::gaussianQuad(0, xmax, [rho](double r) -> double {return rho(r) * r * r;}, 15) << std::endl;

  // All yields 2.000x, which is the correct value

  lambda vEE_dir, vEE_exc, vEE, vEE_old; // All potential functions that we will use

  double eta = 0.4;

  vEE_old = [](double r) -> double {return 0.0;}; // Placeholder which we will replace with the last iterations value

  ///*
  for (int i = 0; i < 15; i++) {
    // Replace with tolerance condition abs(E_nl - E_nl_old) < tol
    std::cout << "Iteration: " << i << std::endl;

    vEE_dir = solver -> solveColl(rho);

    vEE_exc = [rho](double r) -> double {
      return -3.0 * charge / (4.0 * PI * eps0) * std::pow(3.0 * rho(r) / (charge * 8.0 * PI), 1.0 / 3.0);
    };

    vEE = [vEE_dir, vEE_exc, eta, vEE_old](double r) -> double {
      return (vEE_dir(r) + vEE_exc(r))*(1-eta) + eta * vEE_old(r);
    };
    vEE_old = vEE;

    basis = solver -> solveAtomPref(ges, B, H1, 0, Z, vEE);
    p_nl = solver -> getP(0, basis, B);

    rho = [p_nl](double r) -> double {
      if (r == 0) {
        r+= 1e-12;
      }
      return charge / (4 * PI) * 2 * (std::pow(p_nl(r) / r, 2.0)); // 2 comes form N_j, only the 1s state occupied.
    };
    
  }
  //*/
  /*
  vEE_dir = solver -> solveColl(rho);
  vEE_exc = [rho](double r) -> double {
      return -3.0 * charge / (4.0 * PI * eps0) * std::pow(3.0 * rho(r) / (charge * 8.0 * PI), 1.0 / 3.0);
  };
  std::cout << vEE_dir(1.0) << ", " << vEE_exc(1.0) << std::endl;
  vEE = [vEE_dir, vEE_exc, eta, vEE_old](double r) -> double {
      return (vEE_dir(r) + vEE_exc(r))*(1-eta) + eta * vEE_old(r);
  };
  std::cout << vEE(1.0) << std::endl;
  basis = solver -> solveAtomPref(ges, B, H1, 0, Z, vEE);
  p_nl = solver -> getP(0, basis, B);
  std::cout << p_nl(1.0) << std::endl;
  rho = [p_nl](double r) -> double {
      if (r == 0) {
        r+= 1e-12;
      }
      return charge / (4 * PI) * 2 * (std::pow(p_nl(r) / r, 2.0)); // 2 comes form N_j, only the 1s state occupied.
  };
  std::cout << rho(1.0) << std::endl;
  */
}


int main() {
  /*
  cSpline *test = new cSpline;
  std::vector<double> knots;
  // Creating the phyiscal knots
  int xmax = 15;
  double dx = 0.5;
  double x = 0.0;
  int i = 0;
  while (x <= xmax) {
    knots.push_back(x);
    //std::cout << knots[i] << std::endl;
    x += dx;
    i++;
  }
  */
  //Helium();
  test();

  // test -> initialize(11);
  /*
  test -> setKnots(knots);

  test -> solveColl(wf);
  std::ofstream file("test2.dat");
  double r = 0.0;
  double dr = 0.1;
  while (r < 15) {
    file << r << "\t" << wfe(r) << std::endl;
    r += dr;
  }
  file.close();
  // int l, z;
  // z = 1;
  // std::cout << "Input l: ";
  // std::cin >> l;
  // std::cout << std::endl;
  // test -> collMat([](double x, int z) -> double {return 1.0;}); // collMat now works fine
  // std::cout << test -> gaussianQuad(0, 1, [](double x) -> double {return x;}) << std::endl;
 
  //test -> solveAtom(l, z, false);
  */

  return 0;
}
