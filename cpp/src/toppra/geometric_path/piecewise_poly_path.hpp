#ifndef TOPPRA_PIECEWISE_POLY_PATH_HPP
#define TOPPRA_PIECEWISE_POLY_PATH_HPP

#include <toppra/geometric_path.hpp>
#include <toppra/toppra.hpp>
#include <variant>

namespace toppra {

/// Numerical boundary condition
typedef struct {
    int order;
    Vector values;
} NumericalBoundaryCond;

/// Boundary condition
typedef struct {
    enum {string, numeric} type;
    std::variant<NumericalBoundaryCond, std::string> value;
} BoundaryCond;

/**
 * \brief Piecewise polynomial geometric path.
 *
 * An implementation of a piecewise polynomial geometric path.
 *
 * The coefficient vector has shape (N, P, D), where N is the number
 * segments. For each segment, the i-th row (P index) denotes the
 * power, while the j-th column is the degree of freedom. In
 * particular,
 *
 *  coeff(0) * dt ^ 3 + coeff(1) * dt ^ 2 + coeff(2) * dt + coeff(3)
 *
 *
 */
class PiecewisePolyPath : public GeometricPath {
 public:
  PiecewisePolyPath() = default;

  /**
   * \brief Construct new piecewise polynomial.
   *
   * See class docstring for details.
   *
   * @param coefficients Polynomial coefficients.
   * @param breakpoints Vector of breakpoints.
   */
  PiecewisePolyPath(const Matrices &coefficients, std::vector<value_type> breakpoints);

  /**
   * \brief Construct a new piecewise 3rd degree polynomial.
   * @param positions Vectors of path positions at given times.
   * @param times Vector of times in a strictly increasing order.
   * @param bc_type Boundary conditions at the curve start and end.
   */
  PiecewisePolyPath(const Vectors &positions, const Vector &times, const std::array<NumericalBoundaryCond, 2> &bc_type);

  /**
   * \brief Construct a new piecewise 3rd degree polynomial.
   * @param positions Vectors of path positions at given times.
   * @param times Vector of times in a strictly increasing order.
   * @param bc_type the boundary condition type to be applied at both curve endpoints. Available conditions are:
   *    - 'not-a-knot' (default): the first and second segment at a curve end are the same polynomial.
   *    - 'clamped': the first derivatives at curves ends are zero.
   *    - 'natural': the second derivatives at curves end are zero.
   */
  PiecewisePolyPath(const Vectors &positions, const Vector &times, const std::string bc_type = "not-a-knot");

  /**
   * /brief Evaluate the path at given position.
   */
  Vector eval_single(value_type, int order = 0) const override;

  /**
   * /brief Evaluate the path at given positions (vector).
   */
  Vectors eval(const Vector &, int order = 0) const override;

  /**
   * Return the starting and ending path positions.
   */
  Bound pathInterval() const override;
  void serialize(std::ostream &O) const override;
  void deserialize(std::istream &I) override;

  /**
   * \brief Construct a new Hermite polynomial.
   */
  static PiecewisePolyPath constructHermite(const Vectors &positions,
                                            const Vectors &velocities,
                                            const std::vector<value_type> times);

 protected:
  void initAsHermite(const Vectors &positions, const Vectors &velocities,
                     const std::vector<value_type> times);
  Matrices computeCubicSplineCoefficients(const Vectors &positions, const Vector &times,
          const std::array<BoundaryCond, 2> &bc_type);
  std::vector<value_type> computeBreakpoints(const Vector &times,
          const std::array<BoundaryCond, 2> &bc);
  void reset();
  size_t findSegmentIndex(value_type pos) const;
  void checkInputArgs();
  void checkInputArgs(const Vectors &positions, const Vector &times,
                      const std::array<BoundaryCond, 2> &bc_type);
  void computeDerivativesCoefficients();
  const Matrix &getCoefficient(int seg_index, int order) const;
  Matrices m_coefficients, m_coefficients_1, m_coefficients_2, m_coefficients_3;
  std::vector<value_type> m_breakpoints;
  int m_degree;
};

}  // namespace toppra

#endif
