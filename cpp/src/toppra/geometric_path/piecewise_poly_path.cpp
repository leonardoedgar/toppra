#include <iostream>
#include <ostream>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/toppra.hpp>
#include <Eigen/Dense>

#ifdef TOPPRA_OPT_MSGPACK
#include <msgpack.hpp>
#endif

namespace toppra {

Matrix differentiateCoefficients(const Matrix &coefficients) {
  Matrix deriv(coefficients.rows(), coefficients.cols());
  deriv.setZero();
  for (size_t i = 1; i < coefficients.rows(); i++) {
    deriv.row(i) = coefficients.row(i - 1) * (coefficients.rows() - i);
  }
  return deriv;
}

PiecewisePolyPath::PiecewisePolyPath(const Matrices & coefficients,
                                     std::vector<value_type> breakpoints)
    : GeometricPath (coefficients[0].cols()),
      m_coefficients(coefficients), m_breakpoints(std::move(breakpoints)),
      m_degree(coefficients[0].rows() - 1) {

  checkInputArgs();
  computeDerivativesCoefficients();
}

PiecewisePolyPath::PiecewisePolyPath(const Vectors &positions, const Vector &times,
        const std::array<NumericalBoundaryCond, 2> &bc_type) {
    // Prepare input
    BoundaryCond bc_start = {BoundaryCond::numeric, bc_type[0]}, bc_end = {BoundaryCond::numeric, bc_type[1]};
    std::array<BoundaryCond, 2> bcs {bc_start, bc_end};
    checkInputArgs(positions, times, bcs);
    Matrices coefficients = computeCubicSplineCoefficients(positions, times, bcs);
    std::vector<value_type> breakpoints = computeBreakpoints(times, bcs);
    *this = PiecewisePolyPath(coefficients, breakpoints);
}

PiecewisePolyPath::PiecewisePolyPath(const Vectors &positions, const Vector &times, const std::string bc_type) {
    // Prepare input
    BoundaryCond bc = {BoundaryCond::string, bc_type};
    std::transform(std::get<std::string>(bc.value).begin(), std::get<std::string>(bc.value).end(),
            std::get<std::string>(bc.value).begin(), ::tolower);
    std::array<BoundaryCond, 2> bcs = {bc, bc};
    checkInputArgs(positions, times, bcs);
    if (std::get<std::string>(bc.value).compare("clamped") == 0) {
        bc.type = BoundaryCond::numeric;
        bc.value = NumericalBoundaryCond{1, Vector::Zero(positions[0].size())};
        bcs = {bc, bc};
    }
    else if (std::get<std::string>(bc.value).compare("natural") == 0) {
        bc.type = BoundaryCond::numeric;
        bc.value = NumericalBoundaryCond{2, Vector::Zero(positions[0].size())};
        bcs = {bc, bc};
    }
    Matrices coefficients = computeCubicSplineCoefficients(positions, times, bcs);
    std::vector<value_type> breakpoints = computeBreakpoints(times, bcs);
    *this = PiecewisePolyPath(coefficients, breakpoints);
}

Matrices PiecewisePolyPath::computeCubicSplineCoefficients(const Vectors &positions, const Vector &times,
        const std::array<BoundaryCond, 2> &bc_type) {
    // h(i) = t(i+1) - t(i)
    Vector h (times.rows() - 1);
    for (size_t i = 0; i < h.rows(); i++){
        h(i) = times(i + 1) - times(i);
    }

    // Special case for not-a-knot with times.rows() == 2 (1st order derivatives at curve ends = slope)
    std::vector<BoundaryCond> bc(std::begin(bc_type), std::end(bc_type));
    if (times.rows() == 2) {
        for (size_t i = 0; i < bc.size(); i++) {
            if (bc[i].type == BoundaryCond::string && std::get<std::string>(bc[i].value).compare("not-a-knot") == 0) {
                bc[i].type = BoundaryCond::numeric;
                Vector slope (positions[0].rows());
                for (size_t j = 0; j < positions[0].rows(); j++) {
                    slope(j) = (positions[positions.size() - 1](j) - positions[0](j)) / h(0);
                }
                bc[i].value = NumericalBoundaryCond{1, slope};
            }
        }
    }

    // Special case for not-a-knot with times.rows() == 3 (parabola)
    Matrices coefficients;
    if (times.rows() == 3 && bc[0].type == BoundaryCond::string && bc[1].type == BoundaryCond::string &&
        std::get<std::string>(bc[0].value).compare("not-a-knot") == 0 &&
        std::get<std::string>(bc[0].value).compare("not-a-knot") == 0) {

        // Construct the A matrix
        Matrix A (3, 3);
        A << std::pow(times(0), 2), times(0), 1,
            std::pow(times(1), 2), times(1), 1,
            std::pow(times(2), 2), times(2), 1;

        // Construct the B vectors
        Vectors B (positions[0].rows());
        for (size_t i = 0; i < B.size(); i++) {
            B[i].resize(3);
            for (size_t j = 0; j < B[i].rows(); j++) {
                B[i](j) = positions[j](i);
            }
        }

        // Solve AX = B
        Vectors X (positions[0].rows());
        for (size_t i = 0; i < X.size(); i++) {
            X[i].resize(3);
            X[i] = A.colPivHouseholderQr().solve(B[i]);
        }

        // Insert spline coefficients
        coefficients.resize(1);
        coefficients[0].resize(4, positions[0].rows());
        for (size_t i = 0; i < coefficients[0].cols(); i++) {
            coefficients[0](0, i) = 0;
            coefficients[0](1, i) = X[i](0);
            coefficients[0](2, i) = X[i](1);
            coefficients[0](3, i) = X[i](2);
        }

    }
    else {
        // Construct the tri-diagonal matrix A based on spline continuity criteria
        Matrix A = Matrix::Zero(times.rows(), times.rows());
        for (size_t i = 1; i < A.rows() - 1; i++) {
            A.row(i).segment(i - 1, 3) << h(i - 1), 2 * (h(i - 1) + h(i)), h(i);
        }

        // Construct B based on spline continuity criteria
        Vectors B (positions.at(0).rows());
        for (size_t i = 0; i < B.size(); i++) {
            B[i].resize(times.rows());
            for (size_t j = 1; j < B[i].rows() - 1; j++) {
                B[i](j) = 3 * (positions[j + 1](i) - positions[j](i)) / h(j) -
                          3 * (positions[j](i) - positions[j - 1](i)) / h(j - 1);
            }
        }

        // Insert boundary conditions to A and B
        if (bc[0].type == BoundaryCond::string) {
            if (std::get<std::string>(bc[0].value).compare("not-a-knot") == 0) {
                A.row(0).segment(0, 3) << h(1), -(h(0) + h(1)), h(0);
                for (size_t i = 0; i < B.size(); i++) {
                    B[i](0) = 0;
                }
            }
        }
        else {
            NumericalBoundaryCond numeric_bc = std::get<NumericalBoundaryCond>(bc[0].value);
            if (numeric_bc.order == 1) {
                A.row(0).segment(0, 2) << 2 * h(0), h(0);
                for (size_t i = 0; i < B.size(); i++) {
                    B[i](0) = 3 * (positions[1](i) - positions[0](i)) / h(0) - 3 * numeric_bc.values(i);
                }
            }
            else if (numeric_bc.order == 2) {
                A(0, 0) = 2;
                for (size_t i = 0; i < B.size(); i++) {
                    B[i](0) = numeric_bc.values(i);
                }
            }
        }
        if (bc[1].type == BoundaryCond::string) {
            if (std::get<std::string>(bc[1].value).compare("not-a-knot") == 0) {
                A.row(A.rows() - 1).segment(A.cols() - 3, 3) << h(h.rows() - 1), -(h(h.rows() - 2) + h(h.rows() - 1)),
                        h(h.rows() - 2);
                for (size_t i = 0; i < B.size(); i++) {
                    B[i](B[i].rows() - 1) = 0;
                }
            }
        }
        else {
            NumericalBoundaryCond numeric_bc = std::get<NumericalBoundaryCond>(bc[1].value);
            if (numeric_bc.order == 1) {
                A.row(A.rows() - 1).segment(A.cols() - 2, 2) << h(h.rows() - 1), 2 * h(h.rows() - 1);
                for (size_t i = 0; i < B.size(); i++) {
                    B[i](B[i].rows() - 1) = 3 * numeric_bc.values(i) -
                            3 * (positions[positions.size() - 1](i) - positions[positions.size() - 2](i)) / h(h.rows() - 1);
                }
            }
            else if (numeric_bc.order == 2) {
                A(A.rows() - 1, A.cols() - 1) = 2;
                for (size_t i = 0; i < B.size(); i++) {
                    B[i](B[i].rows() - 1) = numeric_bc.values(i);
                }
            }
        }

        // Solve AX = B
        Vectors X (positions[0].rows());
        for (size_t i = 0; i < X.size(); i++) {
            X[i].resize(times.rows());
            X[i] = A.colPivHouseholderQr().solve(B[i]);
        }

        // Insert spline coefficients
        coefficients.resize(times.rows() - 1);
        for (size_t i = 0; i < coefficients.size() ; i++) {
            coefficients[i].resize(4, positions[0].rows());
            for (size_t j = 0; j < coefficients[i].cols(); j++) {
                coefficients[i](0, j) = (X[j](i + 1) - X[j](i)) / (3 * h(i));
                coefficients[i](1, j) = X[j](i);
                coefficients[i](2, j) = (positions[i + 1](j) - positions[i](j)) / h(i) -
                                        h(i) / 3 * (2 * X[j](i) + X[j](i + 1));
                coefficients[i](3, j) = positions[i](j);
            }
        }
    }
    return coefficients;
}

std::vector<value_type> PiecewisePolyPath::computeBreakpoints(const Vector &times,
        const std::array<BoundaryCond, 2> &bc) {
    std::vector<value_type> breakpoints;
    if (times.rows() == 3 && bc[0].type == BoundaryCond::string && bc[1].type == BoundaryCond::string &&
        std::get<std::string>(bc[0].value).compare("not-a-knot") == 0 &&
        std::get<std::string>(bc[0].value).compare("not-a-knot") == 0) {
        breakpoints = {times(0), times(times.rows() - 1)};
    }
    else {
        breakpoints = std::vector<value_type>(times.data(), times.data() + times.rows());
    }
    return breakpoints;
}

Bound PiecewisePolyPath::pathInterval() const {
  Bound v;
  v << m_breakpoints.front(), m_breakpoints.back();
  return v;
};

Vector PiecewisePolyPath::eval_single(value_type pos, int order) const {
  assert(order <= 3 && order >= 0);
  Vector v(m_dof);
  v.setZero();
  size_t seg_index = findSegmentIndex(pos);
  auto coeff = getCoefficient(seg_index, order);
  for (int power = 0; power < m_degree + 1; power++) {
    v += coeff.row(power) *
         pow(pos - m_breakpoints[seg_index], m_degree - power);
  }
  return v;
}

// Not the most efficient implementation. Coefficients are
// recomputed. Should be refactored.
Vectors PiecewisePolyPath::eval(const Vector &positions, int order) const {
  assert(order <= 3 && order >= 0);
  Vectors outputs;
  outputs.resize(positions.size());
  for (size_t i = 0; i < positions.size(); i++) {
    outputs[i] = eval_single(positions(i), order);
  }
  return outputs;
}

size_t PiecewisePolyPath::findSegmentIndex(value_type pos) const {
  size_t seg_index = -1;
  for (size_t i = 0; i < m_coefficients.size(); i++) {
    if (m_breakpoints[i] <= pos && pos <= m_breakpoints[i + 1]) {
      seg_index = i;
      break;
    }
  }
  if (seg_index == -1) {
    std::ostringstream oss;
    oss << "Position " << pos << " is outside of range [ " << m_breakpoints[0]
      << ", " << m_breakpoints[m_breakpoints.size()-1] << ']';
    throw std::runtime_error(oss.str());
  }
  return seg_index;
}

void PiecewisePolyPath::checkInputArgs() {
  assert(m_coefficients[0].cols() == m_dof);
  assert(m_coefficients[0].rows() == (m_degree + 1));
  if ((1 + m_coefficients.size()) != m_breakpoints.size()) {
    throw std::runtime_error(
        "Number of breakpoints must equals number of segments plus 1.");
  }
  for (size_t seg_index = 0; seg_index < m_coefficients.size(); seg_index++) {
    if (m_breakpoints[seg_index] >= m_breakpoints[seg_index + 1]) {
      throw std::runtime_error("Require strictly increasing breakpoints");
    }
  }
}

void PiecewisePolyPath::checkInputArgs(const Vectors &positions, const Vector &times,
        const std::array<BoundaryCond, 2> &bc_type) {
    if (positions.size() != times.rows()) {
        throw std::runtime_error("The length of 'positions' doesn't match the length of 'times'.");
    }
    if (times.rows() < 2) {
        throw std::runtime_error("'times' must contain at least 2 elements.");
    }
    for (size_t i = 1; i < positions.size(); i++) {
        if (positions[i].rows() != positions[i - 1].rows()) {
            throw std::runtime_error("The number of elements in each position has to be equal.");
        }
    }
    Vector dtimes (times.rows() - 1);
    for (size_t i = 1; i < times.rows(); i++) {
        dtimes(i - 1) = times(i) - times(i - 1);
        if (dtimes(i - 1) <= 0) {
            throw std::runtime_error("'times' must be a strictly increasing sequence.");
        }
    }

    // Validate boundary conditions
    int expected_deriv_size = positions[0].rows();
    for (const BoundaryCond &bc: bc_type) {
        if (bc.type == BoundaryCond::string) {
            std::string parsed_bc;
            try {
                parsed_bc = std::get<std::string>(bc.value);
            }
            catch (const std::bad_variant_access&) {
                throw std::runtime_error("`bc.type` string is not consistent with the value stored.");
            }
            if (parsed_bc.compare("clamped") != 0 && parsed_bc.compare("natural") != 0 &&
                parsed_bc.compare("not-a-knot") != 0) {
                throw std::runtime_error("`bc` " + parsed_bc + " is not a valid boundary condition.");
            }
        }
        else {
            NumericalBoundaryCond numeric_bc;
            try {
                numeric_bc = std::get<NumericalBoundaryCond>(bc.value);
            }
            catch (const std::bad_variant_access&) {
                throw std::runtime_error("`bc.type` numeric is not consistent with the value stored.");
            }
            if (numeric_bc.order != 1 && numeric_bc.order != 2) {
                throw std::runtime_error("The specified derivative order must be 1 or 2.");
            }
            if (numeric_bc.values.rows() != expected_deriv_size) {
                throw std::runtime_error(
                        "`deriv_value` size " + std::to_string(numeric_bc.values.rows()) + " is not the expected one " +
                        std::to_string(expected_deriv_size) + ".");
            }
        }
    }
}

void PiecewisePolyPath::computeDerivativesCoefficients() {
  m_coefficients_1.reserve(m_coefficients.size());
  m_coefficients_2.reserve(m_coefficients.size());
  m_coefficients_3.reserve(m_coefficients.size());
  for (size_t seg_index = 0; seg_index < m_coefficients.size(); seg_index++) {
    m_coefficients_1.push_back(
        differentiateCoefficients(m_coefficients[seg_index]));
    m_coefficients_2.push_back(
        differentiateCoefficients(m_coefficients_1[seg_index]));
    m_coefficients_3.push_back(
        differentiateCoefficients(m_coefficients_2[seg_index]));
  }
}

const Matrix &PiecewisePolyPath::getCoefficient(int seg_index, int order) const {
  if (order == 0) {
    return m_coefficients.at(seg_index);
  } else if (order == 1) {
    return m_coefficients_1.at(seg_index);
  } else if (order == 2) {
    return m_coefficients_2.at(seg_index);
  } else {
    return m_coefficients_3.at(seg_index);
  }
}

void PiecewisePolyPath::serialize(std::ostream &O) const {
#ifdef TOPPRA_OPT_MSGPACK
  MatricesData allraw;
  allraw.reserve(m_coefficients.size());
  for (const auto &c : m_coefficients) {
    MatrixData raw{c.rows(), c.cols(), {c.data(), c.data() + c.size()}};
    allraw.push_back(raw);
  }
  msgpack::pack(O, allraw);
  msgpack::pack(O, m_breakpoints);
#endif
}

void PiecewisePolyPath::deserialize(std::istream &I) {
#ifdef TOPPRA_OPT_MSGPACK
  std::stringstream buffer;
  buffer << I.rdbuf();
  std::size_t offset = 0;

  auto oh = msgpack::unpack(buffer.str().data(), buffer.str().size(), offset);
  auto obj = oh.get();
  TOPPRA_LOG_DEBUG(obj << "at offset:=" << offset << "/" << buffer.str().size());
  MatricesData x;
  toppra::Matrices new_coefficients;
  obj.convert(x);
  for (auto const &y : x) {
    int nrow, ncol;
    nrow = std::get<0>(y);
    ncol = std::get<1>(y);
    std::vector<value_type> mdata = std::get<2>(y);
    toppra::Matrix m(nrow, ncol);
    for (size_t i = 0; i < mdata.size(); i++) m(i) = mdata[i];
    TOPPRA_LOG_DEBUG(nrow << ncol << mdata.size() << m);
    new_coefficients.push_back(m);
  }

  reset();
  m_coefficients = new_coefficients;
  oh = msgpack::unpack(buffer.str().data(), buffer.str().size(), offset);
  obj = oh.get();
  TOPPRA_LOG_DEBUG(obj << "at offset:=" << offset << "/" << buffer.str().size());
  assert(offset == buffer.str().size());
  obj.convert(m_breakpoints);

  TOPPRA_LOG_DEBUG("degree: " << m_degree);
  m_dof = new_coefficients[0].cols();
  m_degree = new_coefficients[0].rows() - 1;
  checkInputArgs();
  computeDerivativesCoefficients();
#endif
}

void PiecewisePolyPath::reset() {
  m_breakpoints.clear();
  m_coefficients.clear();
  m_coefficients_1.clear();
  m_coefficients_2.clear();
  m_coefficients_3.clear();
}

void PiecewisePolyPath::initAsHermite(const Vectors &positions,
                                      const Vectors &velocities,
                                      const std::vector<value_type> times) {
  reset();
  assert(positions.size() == times.size());
  assert(velocities.size() == times.size());
  TOPPRA_LOG_DEBUG("Constructing new Hermite polynomial");
  m_configSize = m_dof = positions[0].size();
  m_degree = 3;  // cubic spline
  m_breakpoints = times;
  for (std::size_t i = 0; i < times.size() - 1; i++) {
    TOPPRA_LOG_DEBUG("Processing segment index: " << i << "/" << times.size() - 1);
    Matrix c(4, m_dof);
    auto dt = times[i + 1] - times[i];
    assert(dt > 0);
    // ... after some derivations
    c.row(3) = positions.at(i);
    c.row(2) = velocities.at(i);
    c.row(0) = (velocities.at(i + 1).transpose() * dt -
                2 * positions.at(i + 1).transpose() + c.row(2) * dt + 2 * c.row(3)) /
               pow(dt, 3);
    c.row(1) = (velocities.at(i + 1).transpose() - c.row(2) - 3 * c.row(0) * dt * dt) /
               (2 * dt);
    m_coefficients.push_back(c);
  }
  checkInputArgs();
  computeDerivativesCoefficients();
}

PiecewisePolyPath PiecewisePolyPath::constructHermite(
    const Vectors &positions, const Vectors &velocities,
    const std::vector<value_type> times) {
  PiecewisePolyPath path;
  path.initAsHermite(positions, velocities, times);
  return path;
}

} // namespace toppra
