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
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef std::function<double(double)> lambda;

const double charge = 1.0;
const double hbar = 1.0;
const double m = 1.0;
const double PI = 3.14159265358979323846;
const double a0 = 1.0;
const double eps0 = 1/ (4.0 * PI);

inline double zeroFunc(double x) { return 0.0; }

struct cSpline {

  int knotsNumber;
  int knotsPhyiscal;

  Eigen::VectorXd knots;  // The knots points
  std::vector<lambda> splines;// Splines[i](x)
  std::vector<lambda> dsplines;// d/dx Splines[i](x)
  std::vector<lambda> ddsplines;// d^2/dx^2 Splines[i](x)

  Eigen::VectorXd coeff;
  Mat H2;

  int splineOrder;

  /**
    Constructor: Creates the splines used
  */
  cSpline(int k = 4) {
    this -> splineOrder = k; // Spline order, we can change to 5, 6 and so on
  };

  // Deconstructor
  ~cSpline(){};

  /**
    collmat(func, order): Computes the collocation matrix of the ghost physical points
    @param[in] func: The value on the right hand side,
    @param[in] z: The electron number
    @return {collmat, rhs}
  */
  Mat collMat() {
    // N - k splines, but we set the spline to be zero, since  we set the first value to be zero 
    int n = this -> knotsNumber - this -> splineOrder - 1;
    int k = this -> splineOrder;
    Mat mat = Eigen::MatrixXd::Zero(n,n);
    std::vector<std::function<double(double)> > spl = this -> ddsplines;
    Vec t = this -> knots;
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
    // Right boundary
    mat(n - 1, n - 1) = this -> dsplines[spl.size() - 1](t(t.size() - 1) - shift);
    mat(n - 1, n - 2) = this -> dsplines[spl.size() - 2](t(t.size() - 1) - shift);
    
    return mat;
  };

  Vec collFunc(std::function<double(double)> func) {
    int n = this -> knotsNumber - this -> splineOrder - 1;
    int k = this -> splineOrder;
    Vec rhs = Eigen::VectorXd::Zero(n);
    Vec t = this -> knots;
    double shift = 1e-12;

    for (int j = 1; j < n; j++) {
      rhs(j) = -func(t(j + k - 1)) * t(j + k - 1) * 4.0 * PI / (4.0 * PI * eps0);
    }
    rhs(rhs.size() - 1) = 0;
    return rhs;
  };

  /**
    saveSplines: Save the splines in the domain
    @param[in] order: The derivative order
  */
  void saveSplines(int order = 0) {
    Vec t = this -> knots;
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
    0'th, 1'th and 2'th derivative
  */
  void createSplineFunctions() {
    int k = this -> splineOrder;
    std::vector<lambda> splines;
    std::vector<lambda> dspline;
    std::vector<lambda> ddspline;
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

    if (knotsPhysical < 5) {
      std::cout << "To few knot points, please include more" << std::endl;
      exit(1);
    }

    this -> knotsPhyiscal = knotsPhysical;
    this -> knotsNumber = knotsPhysical + 2 * (this -> splineOrder - 1);
    // If we have N knots, then knots is N, 1 vector / array
    this -> knots = Eigen::VectorXd::Zero(knotsNumber, 1);

    for (int i = 0; i < this -> splineOrder - 1; i++) {// Initialize the ghost points
      this -> knots(i) = knots[0]; // left side
      this -> knots(knotsPhyiscal + this -> splineOrder -1 + i) = knots[knotsPhysical - 1]; // right side
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
    Vec t = this -> knots;
    double sol = 0.0;
    if (k == 1) { // Fírst order spline
      if (t(i) <= x_pos && t(i + 1) > x_pos) { // Initial condition if k = 1, i.e. linear spline
        return 1.0;                         
      } else {
        return 0.0;
      }
    }
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

    Mat gauss = leggauss(n);
    Vec roots = gauss.row(0);
    Vec weights = gauss.row(1);

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
  Mat getB(std::function<double(double)> func = [](double r) -> double {return 1.0;}) {
    int n = this -> knotsNumber - this -> splineOrder - 2;
    // We will Have  n - k splines,
    // we disregard the first and last spline due to bc
    int k = this -> splineOrder;

    Mat mat = Eigen::MatrixXd::Zero(n, n);
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
    @returns the potential term
  */
  static double V(int l, int z, double r) {
    if (r == 0) {
      r += 1e-12;
    }
    return hbar * hbar * l * (l + 1) / (2 * m * r * r) - z * charge * charge / (4.0 * PI * eps0 * r);
  };
  
  /**
    getH1: Computes the matrix H for db_i * db_j for the linear system of equations
    @returns The matrix H1
  */
  Mat getH1() {
    int n = this -> knotsNumber - this -> splineOrder - 2; // -2 from bc
    // We will Have  n - k splines,
    // we disregard the first and last spline due to bc
    int k = this -> splineOrder;

    Mat mat = Eigen::MatrixXd::Zero(n, n);
    std::vector<std::function<double(double)>> dsplines = this -> dsplines; // d B_i^k(x) / dx

    double ti, ti1, element;

    for (int i = 1; i < n + 1; i++) {
      for (int j = 1; j < n + 1; j++) {
        // The function to integrate, which is dB_i^k(x) * dB_j^k(x)
        // We get spline [i + 1] due to the fact that we disregard the first spline, and the last spline
        int min = std::max(i, j);
        int max = std::min(i, j) + k - 1;
        element = 0.0;
        for (int m = min; m <= max; m++) {

          lambda f = [i, j, dsplines](double x) -> double {return dsplines[i](x) * dsplines[j](x);}; // To skip the first and last spline
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
    @returns The matrix H
  */
  Mat getH(int l, int z) {
    Mat H = ( this -> getH1() );

    return H + this -> getB([l, z](double x) -> double {return cSpline::V(l, z, x);});
  };

  /**
    getPertubedHamiltonian: Computes the pertubed hamiltonian
    @param[in] vEE: The additional potential
    @returns The pertubed hamiltonian
  */
  Mat getPertubedHamiltonian(lambda vEE) {
    return this -> getB(vEE);
  };

  /**
    solveColl: Solves the collocation problem
    @param[in] collmat - The collocation matrix
    @returns func - The solution V(r) to the collocation problem
  */
  lambda solveColl(Mat collmat, lambda f){
    Vec res = this -> collFunc(f);
    Vec coeff;
    coeff = collmat.partialPivLu().solve(res);
    std::vector<lambda> splines = this -> splines;

    // Returns V = phi_{nl}(r) / r
    return [coeff, splines](double r) -> double {
      // Optimize by finding the nearest spline and just
      // evaluate close to it instead of evaluating all of them
      double y = 0;
      if (r == 0) {
        r += 1e-12;
      }
      for (int i = 0; i < coeff.size(); i++) {
        y+= coeff[i] * splines[i + 1](r) / r;
      }
      return y;
    };
  };

  /**
    solveAtomPref: Solves the atom with the additional potential
    @param[in] ges: The general eigen solver
    @param[in] B: The matrix B
    @param[in] H1: The matrix H1, the constant part
    @param[in] vEE: The additional potential
    @returns The eigen-vectors & eigen-values
  */
  std::tuple<Mat, Vec> solveAtomPref(Solver ges, Mat B, Mat H1, lambda vEE) {
    Mat H = H1 + charge * (this -> getPertubedHamiltonian(vEE)) ; // The term that changes
    ges.compute(H, B);
    
    return std::make_tuple(ges.eigenvectors(), ges.eigenvalues());
  };
  
  /**
    getP: Computes the P_{nl}(r) function, normalized
    @param[in] k: The k:th order, i.e. the k -s/p/d/f state
    @param[in] basis: The coefficients matrix (spline coefficients for each state)
    @param[in] B: The rhs matrix
    @returns The function P_{nl}(r)/sqrt(norm)
  */
  lambda getP(int k, Mat basis, Mat B) {
    // Can we do like this?
    // return basis.col(k) * splines(r)
    std::vector<lambda> splines = this -> splines;
    double norm = basis.col(k).transpose() * B * basis.col(k); // Normalize
    //std::cout << "Inner norm: " << norm << std::endl;

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
double Riemann(double a, double b, lambda func){
  double sum = 0;
  double dx = 0.01;
  while (a < b) {
    sum += func(a) * dx;

    a+= dx;
  }
  return sum;
}

/**
  Simple Trapezoidal method for integration
*/
double trapz(double a, double b, lambda func) {
  double dx = 0.01;
  double sum = 2 * func(a);
  while (a < b) {
    sum += 2 * func(a);
    a += dx;
  }
  sum += func(b);
  return sum * dx / 2.0;
  
}

lambda normalize(double a, double b, lambda func) {
  double norm = cSpline::gaussianQuad(a, b, [func](double x) -> double {return func(x) * func(x);}, 20);
  //double norm = Riemann(a, b, func);
  std::cout << "Norm: " << std::sqrt(norm) << std::endl;
  return [func, norm](double x) -> double {return func(x) / std::sqrt(norm);};
}

std::vector<double> linspace(double a, double b, int n) {
  std::vector<double> vec;
  double dx = (b - a) / (n - 1);
  for (int i = 0; i < n; i++) {
    vec.push_back(a + i * dx);
  }
  return vec;
}

std::vector<double> denselinspace(double a, double b) {
  std::vector<double> vec;
  while (a < b) {
    vec.push_back(a);
    if (a < b * 0.05) {
      a += 0.05;
    } else if (a < b * 0.3) {
      a += 0.1;
    } else {
      a += 0.5;
    }
  }
  if (vec[vec.size() - 1] > b) {
    vec.pop_back();
    vec.push_back(b);
  } else if (vec[vec.size() - 1] < b) {
    vec.push_back(b);
  }
  return vec;
}

struct Atom {

  int Z, N, maxL; // N_p, N_e l_max;
  double rmax;
  std::map<double, std::pair<int, int>> states;

  cSpline *atomSolver;
  cSpline *collSolver;

  std::vector<lambda> p_n0;
  std::vector<lambda> p_n1;
  std::vector<lambda> p_n2;
  //Vec hydrogenLike;

  Mat Energies;  
  
  // Constructor for the atom struct.
  Atom(int Z, int N, double rmax = 30, int Lmax = 2) {
    this -> Z = Z;
    this -> N = N;
    this -> rmax = rmax;
    this -> maxL = Lmax;
  };

  // makeIon: Makes the atom an ion by reducing the number of electron by one
  void makeIon() {
    this -> N--; // Decrement the number of electrons
  };

  // bindSpline: Binds the collocation and atom splines to the atom.
  void bindSpline(cSpline *a_spl, cSpline *c_spl) {
    this -> atomSolver = a_spl;
    this -> collSolver = c_spl;
  };

  /**
    solveWavef: Solves the Schrödinger equation for the atom.
    @param[in] ges: The generalized eigen value solver
    @param[in] B: The rhs matrix of the problem
    @param[in] vEE: The electron interaction potential
    @param[in] iteration: Internal for self-consistancy usage
  */
  void solveWavef(Solver ges, Mat B, lambda vEE, int iteration = 0) {
    lambda p_nl;
    Mat H1;
    Mat basis;
    Vec eigenvals;

    std::vector<lambda> p_n0;
    std::vector<lambda> p_n1;
    std::vector<lambda> p_n2;

    int nMax = this -> atomSolver -> knotsNumber - this -> atomSolver -> splineOrder - 2;

    Mat Energies = Eigen::MatrixXd::Zero(nMax, this -> maxL + 2);
    
    for (int l = 0; l < this -> maxL + 1; l++) {
      H1 = this -> atomSolver -> getH(l, this -> Z);
      std::tie(basis, eigenvals) = this -> atomSolver -> solveAtomPref(ges, B, H1, vEE);
      Energies.col(l) = eigenvals;
      /*
      if (l == 0 && iteration == 0) {
        this -> hydrogenLike = eigenvals;
      }
      */
      for (int n = 0; n < eigenvals.size() - 1; n ++) {
        p_nl = this -> atomSolver -> getP(n, basis, B);
        if ( eigenvals(n) > 0 ) {break;} // Not a bound state

        if (l == 0) {
          p_n0.push_back(p_nl); // 1s 2s 3s 4s ...
        } else if (l == 1) {
          p_n1.push_back(p_nl); // 2p 3p 4p ...
        } else if (l == 2) {
          p_n2.push_back(p_nl); // 3d 4d ...
        } else {
          break;
        }
      }
    }

    this -> p_n0 = p_n0;
    this -> p_n1 = p_n1;
    this -> p_n2 = p_n2;
    this -> Energies = Energies;
  }

  /**
    occupancy: Computes the lowest energy states and fill the respecitly with electrons
    @returns A matrix containing all the states that have electrons in them
  */
  Mat occupancy() {

    Vec energies;
    int hState = 0;
    for (int l = 0; l < this -> maxL + 1; l++) {
      energies = this -> Energies.col(l);
      for (int i = 0; i < energies.size() - 1; i++) {
        if (energies(i) > 0) {if (i > hState){hState=i;};break;}
        this -> states[energies(i)] = {i + l, l};
      }
    }
    
    // Disregard since we fill with atlest 2
    Mat occ = Eigen::MatrixXd::Zero(hState / 3, this -> maxL + 1); // HighestState x Lmax + 1
    int nOcc = 0;
    int n, l, e;

    for (const auto& [e, nl]: this -> states) {
      int n = nl.first;
      int l = nl.second;
      double mul = 2 * (2 * l + 1);
      for (int i = 0; i < mul; i++) {
        if (nOcc >= this -> N) {break;}
        occ(n, l) += 1;
        nOcc++;
      }
    }    
    std::cout << occ << std::endl;
    // Now it's exescilly large, can we reduice it show how? 
    return occ;
  };

  /**
    computeRho: Computes the charge density rho
    @param[in] occ: The occupancy matrix
    @returns the charge density function rho(r)
  */
  lambda computeRho(Mat occ) {
    return [occ, this](double r) -> double {
     double sum = 0;
     if (r == 0) {
       r += 1e-12;
     }
     for (int n = 0; n < occ.rows(); n++) {
       if (occ(n, 0) != 0) {
        sum += occ(n, 0) * std::pow(this -> p_n0[n - l](r) / r, 2.0); // -l is the shift due to occupaction
       }
       if (occ(n, 1) != 0 ) {
         sum += occ(n, 1) * std::pow(this -> p_n1[n - l](r) / r, 2.0);
        }
        if (occ(n, 2) != 0) {
        sum += occ(n, 2) * std::pow(this -> p_n2[n - l](r) / r, 2.0);
        }
     } 
     return sum * charge / (4.0 * PI);
    };
  };

  /**
    selfConsistancy: Solve the self-consistancy equation
    @param[in] ges: Generalized eigen value solver
    @param[in] occ: The occupation matrix
    @param[in] B: The rhs matrix of problem
  */
  void selfConsistancy(Solver ges, Mat occ, Mat B) {
    lambda vEE_dir, vEE_exc, vEE, vEE_old;
    //vEE_old = [](double r) -> double {return 0.0;};
    vEE_old = zeroFunc;
    double eta = 0.4;
    double e_00 = 20.0; // we use this to check the convergence, a random number so we dont converge on the first run
    double tolerance = 1e-5;
    lambda rho;
    Mat collmat = this -> collSolver -> collMat();

    for (int i = 0; i < 20; i++) {
      std::cout << "Iteration: " << i << std::endl;;
      rho = this -> computeRho(occ);

      vEE_dir = this -> collSolver -> solveColl(collmat, rho);

      vEE_exc = [rho](double r) -> double {
        return -3.0 * charge / (4.0 * PI * eps0) * std::pow(3.0 * rho(r) / (charge * 8.0 * PI), 1.0 / 3.0);
      };

      vEE = [vEE_dir, vEE_exc, eta, vEE_old](double r) -> double {
        return (vEE_dir(r) + vEE_exc(r))*(1-eta) + eta * vEE_old(r);
      };

      vEE_old = vEE;

      this -> solveWavef(ges, B, vEE, i);

      if (std::abs(this -> Energies(0,0) - e_00) < tolerance) {break;}

      e_00 = this -> Energies(0,0);
    }

    // We want to save: VEE, all pnls and the total energy
    std::string filename = "Atom" + std::to_string(this -> Z); 
    if (this -> Z != this -> N) {
      filename += "Ion";
    }
    filename += ".dat";
    std::ofstream vFile("Potential_" + filename);
    std::ofstream p_squared("Squred" + filename);
    std::ofstream eFile("Energy" + filename);
    double r = 0;
    double dr = 0.1;

    while (r < this -> rmax) {

      vFile << r << "\t" << vEE(r) << std::endl;
      p_squared << r << "\t" << 4.0 * PI * r * r * rho(r) << std::endl;
      
      r += dr;      
    }

    vFile.close();
    p_squared.close();
    eFile << "Total energy: " << this -> totalEnergy(occ, vEE) << std::endl;;
    eFile.close();

  };

  double totalEnergy(Mat occ, lambda vEE) {

    double sum = 0;
    lambda f;
    // std::cout << this -> Energies << std::endl;
    // std::cout << occ << std::endl;
    for (int i = 0; i < occ.rows(); i++) {
      for (int l = 0; l < occ.cols(); l++) {
        if (occ(i, l) == 0) {continue;}
        if (l == 0) {
           f = [this, i, vEE](double r) -> double {return std::pow(this -> p_n0[i](r), 2.0) * vEE(r);};
        }else if (l == 1) {
           f = [this, i, vEE](double r) -> double {return std::pow(this -> p_n1[i](r), 2.0) * vEE(r);};
        }else if (l == 2) {
           f = [this, i, vEE](double r) -> double {return std::pow(this -> p_n2[i](r), 2.0) * vEE(r);};
        }
        sum += occ(i,l) * (this -> Energies(i, l) - 0.5 * trapz(0, this -> rmax, f));
      }
    }
    return sum;
  };

  /**
    solve: Solves the atom by computing the self-consistancy equation
  */
  void solve() {

    std::vector<double> t1 = denselinspace(0, this -> rmax); // for the collocation solver
    std::vector<double> t2 = linspace(0, this -> rmax, int(this -> rmax * 2)); // for the atom solver

    cSpline *collSolver = new cSpline(4); // cubic bspline
    cSpline *atomSolver = new cSpline(4); // cubic bspline
    collSolver -> setKnots(t1);
    atomSolver -> setKnots(t2);
    this -> bindSpline(atomSolver, collSolver);

    Solver ges;

    lambda vEE = zeroFunc;

    Mat B = this -> atomSolver -> getB();

    this -> solveWavef(ges, B, vEE);
    Mat occ = this -> occupancy();
    lambda rho = this -> computeRho(occ);

    this -> checkRho(rho);
    this -> selfConsistancy(ges, occ, B);

    std::cout << "Computing Ion now\n" << std::endl;
    
    if (this -> Z != this -> N) { return; }

    this -> makeIon(); // this -> N--;
    this -> solve(); // Recursive call but making it an ion instead
    
  };

  /**
    checkRho: Computes the total number of electrons in the system
    @param[in] f: The charge density function
  */
  void checkRho(lambda f) {
    double rmin = 0;
    double rmax = this -> rmax;
    std::cout << "Total number of electrons: " << 4.0 * PI * Riemann(rmin, rmax, [f](double r) -> double {return f(r) * r* r;}) << std::endl;
  };
};

void HeliumProblem() {

  Atom *Helium = new Atom(2, 2, 30);

  Helium -> solve();

  delete Helium;
}

void NeonProblem() {

  Atom *Neon = new Atom(10, 10, 75);

  Neon -> solve();

  delete Neon;
}

void ArgonProblem() {

  Atom *Argon = new Atom(18, 18, 150);

  Argon-> solve();

  delete Argon;
}

void GeneralProblem(int z, double rmax) {
  Atom *General = new Atom(z, z, rmax);

  General -> solve();

  delete General;
}


int main() {

  //HeliumProblem();

  NeonProblem();

  //ArgonProblem();

  /*
  int z;
  int rmax;
  std::cout << "Input the number of protons: ";
  std::cin >> z;
  std::cout << "Input the maximum radii: ";
  std::cin >> rmax;
  GeneralProblem(z, rmax);
  */

  return 0;
}
