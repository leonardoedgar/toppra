#include <toppra/parametrizer/spline.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>

namespace toppra {

namespace parametrizer {

Spline::Spline(GeometricPathPtr path, const Vector& gridpoints, const Vector& vsquared)
    : Parametrizer(path, gridpoints, vsquared) {
        std::vector<value_type> t_grid (gridpoints.rows(), 0);
        std::vector<value_type> copied_gridpoints (gridpoints.data(), gridpoints.data() + gridpoints.rows());
        std::vector<int> skip_ent;
        value_type sd_average, delta_t;
        for (int i = 1; i < t_grid.size(); i++) {
            sd_average = (m_vs[i - 1] + m_vs[i]) / 2;
            if (sd_average > TOPPRA_NEARLY_ZERO) {
                delta_t = (gridpoints[i] - gridpoints[i - 1]) / sd_average;
            }
            else {
                delta_t = 5;
            }
            t_grid[i] = t_grid[i - 1] + delta_t;
            if (delta_t < TOPPRA_NEARLY_ZERO) {
                skip_ent.push_back(i);
            }
        }

        for (int i = skip_ent.size() - 1; i > -1; i--) {
            t_grid.erase(t_grid.begin() + skip_ent[i]);
            copied_gridpoints.erase(copied_gridpoints.begin() + skip_ent[i]);
        }

        Vectors q_grid = path->eval(Eigen::Map<Vector>(copied_gridpoints.data(), copied_gridpoints.size()));
        m_ts = Eigen::Map<Vector>(t_grid.data(), t_grid.size());

        Bound path_interval = m_path->pathInterval();
        NumericalBoundaryCond init_bc, final_bc;
        init_bc.order = final_bc.order = 1;
        init_bc.values = m_path->eval_single(path_interval[0], 1) * m_vs[0];
        final_bc.values = m_path->eval_single(path_interval[1], 1) * m_vs[m_vs.rows() - 1];
        m_path = std::make_shared<PiecewisePolyPath>(q_grid, m_ts, std::array<NumericalBoundaryCond, 2>{init_bc, final_bc});
}

Vectors Spline::eval_impl(const Vector &times, int order) const {
    return m_path->eval(times, order);
}

bool Spline::validate_impl() const {
    return true;
}

Bound Spline::pathInterval_impl() const {
    Bound b;
    b << m_ts[0], m_ts[m_ts.size() - 1];
    return b;
}

#ifdef OPENRAVE_FOUND
void Spline::computeRaveTrajectory(const OpenRAVE::RobotBasePtr probot, OpenRAVE::TrajectoryBasePtr ptraj) {
    OpenRAVE::ConfigurationSpecification spec = probot->GetActiveConfigurationSpecification("cubic");
    spec.AddDerivativeGroups(1, false);
    spec.AddDerivativeGroups(2, true);
    ptraj->Init(spec);
    std::vector<OpenRAVE::dReal> deltas (m_ts.rows());
    deltas[0] = 0;

    for (size_t i = 0; i < m_ts.rows() - 1; i++) {
        deltas[i + 1] = m_ts[i + 1] - m_ts[i];
    }

    Vectors qs = eval(m_ts), qds = eval(m_ts, 1),
        qdds = eval(m_ts, 2);
    std::vector<OpenRAVE::dReal> traj_data (qs[0].rows() * 3 + 1), q, qd, qdd;
    for (size_t i = 0; i < m_ts.rows(); i++) {
        q = std::vector<OpenRAVE::dReal>(qs[i].data(), qs[i].data() + qs[i].rows());
        qd = std::vector<OpenRAVE::dReal>(qds[i].data(), qds[i].data() + qds[i].rows());
        qdd = std::vector<OpenRAVE::dReal>(qdds[i].data(), qdds[i].data() + qdds[i].rows());
        std::move(q.begin(), q.end(), traj_data.begin());
        std::move(qd.begin(), qd.end(), traj_data.begin() + q.size());
        std::move(qdd.begin(), qdd.end(), traj_data.begin() + qd.size());
        traj_data[traj_data.size() - 1] = deltas[i];
        ptraj->Insert(ptraj->GetNumWaypoints(), traj_data);
    }
}
#endif

}  // namespace parametrizer

}  // namespace toppra
