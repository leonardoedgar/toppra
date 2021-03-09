#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/toppra.hpp>
#include <iostream>
#include <algorithm>
#include "gtest/gtest.h"
#include "utils.hpp"

namespace toppra {

NumericalBoundaryCond makeNumBoundaryCond(const int order, const std::vector<value_type> &values) {
    NumericalBoundaryCond cond;
    cond.order = order;
    cond.values.resize(values.size());
    for (std::size_t i = 0; i < values.size(); i++) cond.values(i) = values[i];
    return cond;
}

class CubicSpline : public testing::Test {
protected:
    CubicSpline() {}
    virtual ~CubicSpline() {}
    static void SetUpTestSuite() {
        positions =
                makeVectors({{1.3, 2.1, 4.35, 2.14, -7.31, 4.31},
                             {1.5, -4.3, 1.23, -4.3, 2.13, 6.24},
                             {-3.78, 1.53, 8.12, 12.75, 9.11, 5.42},
                             {6.25, 8.12, 9.52, 20.42, 5.21, 8.31},
                             {7.31, 3.53, 8.41, 9.56, -3.15, 4.83}});
        times.resize (5);
        times << 0.9, 1.3, 2.2, 2.4, 2.6;
    }

    void ConstructCubicSpline(const Vectors &positions, const Vector &times,
                              const std::array<NumericalBoundaryCond, 2> &bc_type) {
        path = toppra::PiecewisePolyPath(positions, times, bc_type);
    }

    void ConstructCubicSpline(const Vectors &positions, const Vector &times, const std::string bc_type) {
        path = toppra::PiecewisePolyPath(positions, times, bc_type);
    }

    void AssertSplineKnots(const Vectors &positions, const Vector &times) {
        Vectors actual_positions = path.eval(times, 0);
        for (size_t i = 0; i < actual_positions.size(); i++) {
            ASSERT_TRUE(positions[i].isApprox(actual_positions[i]));
        }
    }

    void AssertSplineBoundaryConditions(const std::array<NumericalBoundaryCond, 2> &bc_type) {
        ASSERT_TRUE((path.eval_single(times(0), bc_type[0].order) - bc_type[0].values).norm() < TOPPRA_NEARLY_ZERO);
        ASSERT_TRUE((path.eval_single(times(times.rows() - 1), bc_type[1].order) - bc_type[1].values).norm()
            < TOPPRA_NEARLY_ZERO);
    }

    void AssertSplineBoundaryConditions(const Vector &times, const std::string bc_type) {
        std::string bc_str;
        bc_str.resize(bc_type.size());
        std::transform(bc_type.begin(), bc_type.end(), bc_str.begin(), ::tolower);
        if (bc_str.compare("clamped") == 0) {
            NumericalBoundaryCond bc = NumericalBoundaryCond{1, Vector::Zero(positions[0].size())};
            AssertSplineBoundaryConditions(std::array<NumericalBoundaryCond, 2>{bc, bc});
        }
        else if (bc_str.compare("natural") == 0) {
            NumericalBoundaryCond bc = NumericalBoundaryCond{2, Vector::Zero(positions[0].size())};
            AssertSplineBoundaryConditions(std::array<NumericalBoundaryCond, 2>{bc, bc});
        }
        else if (bc_str.compare("not-a-knot") == 0) {
            if (times.rows() == 2) {
                ASSERT_TRUE((path.eval_single(times(0), 3) - path.eval_single(times(1), 3)).norm() < TOPPRA_NEARLY_ZERO);
            }
            else if (times.rows() == 3) {
                ASSERT_TRUE((path.eval_single(times(0), 3) - path.eval_single(times(1), 3)).norm() < TOPPRA_NEARLY_ZERO);
                ASSERT_TRUE((path.eval_single(times(1), 3) - path.eval_single(times(2), 3)).norm() < TOPPRA_NEARLY_ZERO);
            }
            else {
                ASSERT_TRUE((path.eval_single(times(0), 3) - path.eval_single(times(1), 3)).norm() < TOPPRA_NEARLY_ZERO);
                ASSERT_TRUE((path.eval_single(times(times.rows() - 2), 3) -
                    path.eval_single(times(times.rows() - 1), 3)).norm() < TOPPRA_NEARLY_ZERO);
            }
        }
    }

    toppra::PiecewisePolyPath path;
    static Vectors positions;
    static Vector times;
};

Vectors CubicSpline::positions;
Vector CubicSpline::times;

TEST_F(CubicSpline, ClampedSpline) {
    NumericalBoundaryCond bc = makeNumBoundaryCond(1, {0, 0, 0, 0, 0, 0});
    std::array<NumericalBoundaryCond, 2> bc_type {bc, bc};
    ConstructCubicSpline(positions, times, bc_type);
    AssertSplineKnots(positions, times);
    AssertSplineBoundaryConditions(bc_type);

    ConstructCubicSpline(positions, times, "clamped");
    AssertSplineKnots(positions, times);
    AssertSplineBoundaryConditions(times, "clamped");
}

TEST_F(CubicSpline, NaturalSpline) {
    NumericalBoundaryCond bc = makeNumBoundaryCond(2, {0, 0, 0, 0, 0, 0});
    std::array<NumericalBoundaryCond, 2> bc_type {bc, bc};
    ConstructCubicSpline(positions, times, bc_type);
    AssertSplineKnots(positions, times);
    AssertSplineBoundaryConditions(bc_type);

    ConstructCubicSpline(positions, times, "natural");
    AssertSplineKnots(positions, times);
    AssertSplineBoundaryConditions(times, "natural");
}

TEST_F(CubicSpline, NotAKnotSpline) {
    ConstructCubicSpline(positions, times, "not-a-knot");
    AssertSplineKnots(positions, times);
    AssertSplineBoundaryConditions(times, "not-a-knot");

    Vectors positions = makeVectors({{1.3, 2.1, 4.35, 2.14, -7.31, 4.31},
                                     {1.5, -4.3, 1.23, -4.3, 2.13, 6.24}});
    Vector times (2);
    times << 0.0, 1.5;
    ConstructCubicSpline(positions, times, "not-a-knot");
    AssertSplineKnots(positions, times);
    AssertSplineBoundaryConditions(times, "not-a-knot");

    positions = makeVectors({{1.3, 2.1, 4.35, 2.14, -7.31, 4.31},
                             {1.5, -4.3, 1.23, -4.3, 2.13, 6.24},
                             {7.31, 3.53, 8.41, 9.56, -3.15, 4.83}});
    times.resize(3);
    times << 0.0, 0.5, 1.2;
    ConstructCubicSpline(positions, times, "not-a-knot");
    AssertSplineKnots(positions, times);
    AssertSplineBoundaryConditions(times, "not-a-knot");
}

TEST_F(CubicSpline, FirstOrderBoundaryConditions) {
    std::array<NumericalBoundaryCond, 2> bc_type {makeNumBoundaryCond(1, {1.25, 0, 4.12, 1.75, 7.43, 5.31}),
                                         makeNumBoundaryCond(1, {3.51, 5.32, 4.63, 0, -3.12, 3.53})};
    ConstructCubicSpline(positions, times, bc_type);
    AssertSplineKnots(positions, times);
    AssertSplineBoundaryConditions(bc_type);
}

TEST_F(CubicSpline, SecondOrderBoundaryConditions) {
    std::array<NumericalBoundaryCond, 2> bc_type {makeNumBoundaryCond(2, {1.52, -4.21, 7.21, 9.31, -1.53, 7.54}),
                                         makeNumBoundaryCond(2, {-5.12, 8.21, 9.12, 5.12, 24.12, 9.42})};
    ConstructCubicSpline(positions, times, bc_type);
    AssertSplineKnots(positions, times);
    AssertSplineBoundaryConditions(bc_type);
}

TEST_F(CubicSpline, FirstOrderAndSecondOrderBoundaryConditions) {
    std::array<NumericalBoundaryCond, 2> bc_type {makeNumBoundaryCond(1, {1.52, -4.21, 7.21, 9.31, -1.53, 7.54}),
                                         makeNumBoundaryCond(2, {-5.12, 8.21, 9.12, 5.12, 24.12, 9.42})};
    ConstructCubicSpline(positions, times, bc_type);
    AssertSplineKnots(positions, times);
    AssertSplineBoundaryConditions(bc_type);
}

TEST_F(CubicSpline, BadPositions) {
    Vectors bad_positions = makeVectors({{1.3, 2.1, 4.35, 2.14, -7.31, 4.31},
                                         {1.5, -4.3, 1.23, -4.3, 2.13, 6.24},
                                         {-3.78, 1.53, 8.12, 12.75, 9.11, 5.42},
                                         {6.25, 8.12, 9.52, 5.21, 8.31},
                                         {7.31, 3.53, 8.41, 9.56, -3.15, 4.83}});
    std::array<NumericalBoundaryCond, 2> bc_type {makeNumBoundaryCond(2, {1.52, -4.21, 7.21, 9.31, -1.53, 7.54}),
                                         makeNumBoundaryCond(2, {-5.12, 8.21, 9.12, 5.12, 24.12, 9.42})};
    ASSERT_THROW(ConstructCubicSpline(bad_positions, times, bc_type), std::runtime_error);

    bad_positions = makeVectors({{1.3, 2.1, 4.35, 2.14, -7.31, 4.31},
                                 {1.5, -4.3, 1.23, -4.3, 2.13, 6.24},
                                 {-3.78, 1.53, 8.12, 12.75, 9.11, 5.42},
                                 {6.25, 8.12, 11.23, 9.52, 5.21, 8.31}});
    ASSERT_THROW(ConstructCubicSpline(bad_positions, times, bc_type), std::runtime_error);
}

TEST_F(CubicSpline, BadTimes) {
    std::array<NumericalBoundaryCond, 2> bc_type {makeNumBoundaryCond(2, {1.52, -4.21, 7.21, 9.31, -1.53, 7.54}),
                                         makeNumBoundaryCond(2, {-5.12, 8.21, 9.12, 5.12, 24.12, 9.42})};
    Vector bad_times (1);
    bad_times << 1;
    ASSERT_THROW(ConstructCubicSpline(positions, bad_times, bc_type), std::runtime_error);

    bad_times.resize(4);
    bad_times << 0.11, 1.23, 4.35, 6.75;
    ASSERT_THROW(ConstructCubicSpline(positions, bad_times, bc_type), std::runtime_error);

    bad_times.resize(5);
    bad_times << 0.11, 1.23, 4.35, 6.75, 1.23;
    ASSERT_THROW(ConstructCubicSpline(positions, bad_times, bc_type), std::runtime_error);
}

TEST_F(CubicSpline, BadBoundaryConditions) {
    std::array<NumericalBoundaryCond, 2> bc_type {makeNumBoundaryCond(2, {1.52, 7.21, 9.31, -1.53, 7.54}),
                                         makeNumBoundaryCond(2, {-5.12, 8.21, 9.12, 5.12, 24.12, 9.42})};
    ASSERT_THROW(ConstructCubicSpline(positions, times, bc_type), std::runtime_error);

    bc_type = {makeNumBoundaryCond(2, {1.52, 7.21, 9.31}),
               makeNumBoundaryCond(2, {-5.12, 8.21, 9.12})};
    ASSERT_THROW(ConstructCubicSpline(positions, times, bc_type), std::runtime_error);

    bc_type = {makeNumBoundaryCond(3, {1.52, 7.21, 9.31, 2.52, 4.41, 5.54}),
               makeNumBoundaryCond(3, {-5.12, 8.21, 9.12, -1.32, 3.53, 9.21})};
    ASSERT_THROW(ConstructCubicSpline(positions, times, bc_type), std::runtime_error);

    ASSERT_THROW(ConstructCubicSpline(positions, times, "periodic"), std::runtime_error);
}

}  // namespace toppra
