#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <fstream>

#include "../../report5/code/Eigen/Dense"
#include "../../report5/code/Eigen/Eigen"
#include "../../report5/code/Eigen/Eigenvalues"



typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::ArrayXd Arr;

const double charge = 1.0;
const double hbar = 1.0;
const double m = 1.0;
const double epsilon_0 = 1.0;

// https://github.com/haranjackson/LegendreGauss/tree/master
// NOT MY WORK, CREDITS TO THE AUTHOR
Arr legval(Vec x, Vec c)
{   /*
    Evaluate a Legendre series at points x.
    If `c` is of length `n + 1`, this function returns the value:

    p(x) = c_0 * L_0(x) + c_1 * L_1(x) + ... + c_n * L_n(x)

    `p(x)` has the same shape as `x`.
    Trailing zeros in the coefficients will be used in the evaluation, so
    they should be avoided if efficiency is a concern.

    Input
    ----------
    x : array

    c : array
        Array of coefficients ordered so that the coefficients for terms of
        degree n are contained in c[n].

    Output
    -------
    values : array

    Notes
    -----
    The evaluation uses Clenshaw recursion, aka synthetic division.
    */

    int nc = c.size();
    int n = x.size();
    Arr ret(n);
    ret.setZero(n);

    if (nc == 1)
    {
        ret += c(0);
    }
    else if (nc == 2)
    {
        ret = c(0) + c(1) * x.array();
    }
    else
    {
        int nd = nc-1;
        double c0 = c(nc-3) - (c(nc-1) * (nd-1)) / nd;
        Arr c10 = c(nc-2) + (c(nc-1) * x.array() * (2*nd-1)) / nd;

        if (nc == 3)
        {
            ret = c0 + c10 * x.array();
        }
        else
        {
            nd -= 1;
            Arr c00 = c(nc-4) - (c10 * (nd-1)) / nd;
            Arr c11 = c0 + (c10 * x.array() * (2*nd-1)) / nd;

            for (int i=5; i<nc+1; i++)
            {
                Arr tmp = c00;
                nd -= 1;
                c00 = c(nc-i) - (c11 * (nd-1)) / nd;
                c11 = tmp + (c11 * x.array() * (2*nd-1)) / nd;
            }
            ret = c00 + c11 * x.array();
        }
    }
    return ret;
}

// NOT MY WORK, CREDITS TO THE AUTHOR
Mat legcompanion(Vec c)
{   /*
    Return the scaled companion matrix of c.
    The basis polynomials are scaled so that the companion matrix is
    symmetric when `c` is an Legendre basis polynomial. This provides
    better eigenvalue estimates than the unscaled case and for basis
    polynomials the eigenvalues are guaranteed to be real if
    `numpy.linalg.eigvalsh` is used to obtain them.

    Input
    ----------
    c : array
        1-D array of Legendre series coefficients ordered from low to high
        degree.

    Output
    -------
    mat : array
        Scaled companion matrix of dimensions (deg, deg).
    */

    int n = c.size()-1;
    Mat mat(n,n);
    mat.setZero(n,n);
    Vec scl(n);
    for (int i=0; i<n; i++)
    {
        scl(i) = 1./sqrt(2*i+1);
    }
    for (int i=0; i<n-1; i++)
    {
        double tmp = (i+1) * scl(i) * scl(i+1);
        mat(1+i*(n+1)) = tmp;
        mat(n+i*(n+1)) = tmp;
    }
    return mat;
}

Vec legder(Vec c)
{   /*
    Differentiate a Legendre series.
    Returns the Legendre series coefficients `c` differentiated once.
    he argument `c` is an array of coefficients from low to high degree,
    e.g. [1,2,3] represents the series ``1*L_0 + 2*L_1 + 3*L_2``.

    Input
    ----------
    c : array
        Array of Legendre series coefficients.

    Output
    -------
    der : array
        Legendre series of the derivative.

    Notes
    -----
    In general, the result of differentiating a Legendre series does not
    resemble the same operation on a power series. Thus the result of this
    function may be "unintuitive," albeit correct.
    */

    int n = c.size()-1;
    Vec der(n);
    der.setZero(n);
    for (int j=n; j>2; j--)
    {
        der(j-1) = (2*j - 1) * c(j);
        c(j-2) += c(j);
    }
    if (n > 1)
    {
        der(1) = 3*c(2);
    }
    der(0) = c(1);
    return der;
}

// NOT MY WORK, CREDITS TO THE AUTHOR
Mat leggauss(int deg)
{   /*
    Computes the nodes and weights for Gauss-Legendre quadrature.
    These nodes and weights will correctly integrate polynomials of
    degree < 2*deg over the interval [-1, 1] with the weight function w(x) = 1.

    Input
    ----------
    deg : int
        Number of sample points and weights (must be >= 1)

    Output
    -------
    x : array
        1D array containing the nodes
    w : array
        1D array containing the weights

    Notes
    -----
    The results have only been tested up to degree 100, higher degrees may
    be problematic. The weights are determined by using the fact that
    w_k = c / (L'_n(x_k) * L_{n-1}(x_k))
    where c is a constant independent of k and x_k is the kth root of L_n,
    and then scaling the results to get the right value when integrating 1.
    */

    // First approximation of roots. We use the fact that the companion
    // matrix is symmetric in this case in order to obtain better zeros.
    Vec c(deg+1);
    c.setZero(deg+1);
    c(deg) = 1;

    Mat m = legcompanion(c);
    Eigen::SelfAdjointEigenSolver<Mat> eigs(m);
    Vec x = eigs.eigenvalues();

    // Improve roots by one application of Newton.
    Arr dy = legval(x, c);
    Arr df = legval(x, legder(c));
    x -= (dy/df).matrix();

    // Compute the weights. Factor is scaled to avoid numerical overflow.
    Arr fm = legval(x, c.tail(deg));
    fm /= fm.matrix().lpNorm<Eigen::Infinity>();
    df /= df.matrix().lpNorm<Eigen::Infinity>();
    Vec w(deg);
    w = (1. / (fm * df)).matrix();

    // Symmetrize
    w = (w + w.reverse())/2;
    x = (x - x.reverse())/2;

    // Scale w to get the right value
    w *= 2. / w.sum();

    Mat ret(2, deg);
    for (int i=0; i<deg; i++)
    {
        ret(0,i) = x(i);
        ret(1,i) = w(i);
    }
    return ret;
}

struct cSpline {


  int knotsNumber;

  int knotsPhyiscal;

  Eigen::VectorXd knots;  // The knots points

  std::vector<std::function<double(double)> > splines;// Splines[i](x)
  std::vector<std::function<double(double)> > dsplines;// d/dx Splines[i](x)

  Eigen::MatrixXd matrix; // Left hand side matrix.

  Eigen::VectorXd rhs;    // Right hand side of the matrix eq.

  Eigen::VectorXd coeff;  // Coeffecients in f(t_i) = sum_(n = i - k + 1)^i c_n B_i_k(t_i)
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
    createSplineFunctions: Creates the spline functions
  */
  void createSplineFunctions() {
    int k = this -> splineOrder;
    std::vector<std::function<double(double)> > splines;
    std::vector<std::function<double(double)> > dspline;
    for (int i = 0; i < this -> knotsNumber - this -> splineOrder; i++) {
      splines.push_back([i, k, this](double x) -> double {return this -> buildSpline(x, i, k);});
      dspline.push_back([i, k, this](double x) -> double {return this -> firstDerivative(i, k, x);});
    }
    this -> splines = splines;
    this -> dsplines = dspline;

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
    // Prepare the splines and the first derivatives
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
      sol -= (k - 1) * buildSpline(x, i + 1, k) / (cond2);
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


  static double V(int l, double z, double r) {
    return - l * (l + 1) / (2 * r * r) + z * charge * charge / (hbar * hbar * 4 * 3.1415 * epsilon_0 * r);
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

          std::function<double(double)> f = [i, j, dsplines](double x) -> double {return dsplines[i](x) * dsplines[j](x);}; // To skip the first and last spline
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
    getH: Computes the matrix H for the linear system of equations
    @param[in] l: The angular momentum
    @param[in] z: The charge
    @returns The matrix H
  */
  // private
  Eigen::MatrixXd getH(int l, double z) {
    Eigen::MatrixXd H = this -> getH1();

    if (l != 0 && z != 0) {
      H += this -> getB([l, z](double x) -> double {return cSpline::V(l, z, x);});
    }
    return H;
  };

  void save(std::string filename, Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges) {
    filename += ".dat";

    std::ofstream file(filename);

    Eigen::VectorXd t = this -> knots;
    std::vector<std::function<double(double)>> splines = this -> splines;
    Eigen::MatrixXd mat = ges.eigenvectors().real();

    double r = 0.0;
    double dr = 0.1;
    file << "r\t1s\t2s\t3s\n";

    double f;
    double f1;
    double f2;

    while (r < 10) {
      f = 0;
      for (int i = 1; i < splines.size() - 1; i++) {
        //f += splines[i](r) * mat(0, i - 1) * r;
        //f1 += splines[i](r) * mat(1, i - 1) * r;
        // f2+= splines[i](r) * mat(2, i - 1) * r;
        f += splines[i](r) * mat(i - 1, 0) / r;
        f1 += splines[i](r) * mat(i - 1, 1) / r;
        f2+= splines[i](r) * mat(i - 1, 2) / r;
      }
      file << r << "\t" << f << "\t" << f1 << "\t" << f2 << "\n";
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
  void solve(int l = 0, double z = 1.0) {
    Eigen::MatrixXd B = this -> getB();
    Eigen::MatrixXd H = this -> getH(l, z);
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
    /*
    std::cout << "B: \n";
    std::cout << B << std::endl;
    std::cout << "H: \n";
    std::cout << H << std::endl;
    */

    ges.compute(H, B);
    //std::cout << ges.eigenvalues() << "\n" << std::endl;
    //std::cout << ges.alphas() << std::endl;
    //std::cout << ges.eigenvectors() << std::endl;
    
    std::string filename = "l" + std::to_string(l) + "z" + std::to_string(z);
    this -> save(filename, ges);
    
  };

};


double linear(double x) {
  return x;
}


int main() {
  cSpline *test = new cSpline;
  test -> initialize(11);
  //test -> testCase(11, "test.dat");
  test -> solve();
  //std::cout << test -> gaussianQuad(0, 1, &linear) << std::endl;

  //test -> testCase(11, "test.dat");
  
  /*
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
  */

  //knot -> save("Value.dat");
  

  return 0;
}
