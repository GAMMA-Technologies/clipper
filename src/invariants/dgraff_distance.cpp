/**
 * @file pointnormal_distance.cpp
 * @brief Pairwise geometric invariant between point-normals (e.g., planes)
 * @author Parker Lusk <plusk@mit.edu>
 * @date 15 May 2021
 */

#include "clipper/invariants/dgraff_distance.h"
#include <fstream>

namespace clipper {
namespace invariants {

inline void log_to_file(const std::string &message) {
  std::ofstream log("/tmp/gamma_slam_cpp_debug_log.txt", std::ios::app);
  if (log.is_open()) {
    log << message << std::endl;
  }
}

inline Eigen::Matrix<double, 4, 2>
NormalDistance::getGraffMatrixFromDatum(const Eigen::Matrix<double, 6, 1> &d) {
  Eigen::Vector3d n = d.tail<3>().normalized();

  // Find two orthonormal basis vectors spanning the plane
  Eigen::Vector3d t1;
  if (std::abs(n.x()) > std::abs(n.z()))

    t1 = Eigen::Vector3d(-n.y(), n.x(), 0.0);
  else
    t1 = Eigen::Vector3d(0.0, -n.z(), n.y());
  t1.normalize();
  Eigen::Vector3d t2 = n.cross(t1).normalized();

  // Point on plane (first  3 components)
  Eigen::Vector3d p0 = d.head<3>();

  // Build 3×2 basis
  Eigen::Matrix<double, 3, 2> U;
  U.col(0) = t1;
  U.col(1) = t2;

  // Build 4×2 matrix
  Eigen::Matrix<double, 4, 2> Y;
  Y.block<3, 2>(0, 0) = U;
  Eigen::RowVector2d lastRow = (-p0.transpose() * U); // as in paper
  Y.row(3) = lastRow;

  return Y;
}

inline double NormalDistance::dGraff(const Eigen::Matrix<double, 4, 2> &Y1,
                                     const Eigen::Matrix<double, 4, 2> &Y2) {

  Eigen::HouseholderQR<Eigen::Matrix<double, 4, 2>> qr1(Y1);
  Eigen::Matrix<double, 4, 2> Q1 =
      qr1.householderQ() * Eigen::Matrix<double, 4, 2>::Identity();
  Eigen::HouseholderQR<Eigen::Matrix<double, 4, 2>> qr2(Y2);
  Eigen::Matrix<double, 4, 2> Q2 =
      qr2.householderQ() * Eigen::Matrix<double, 4, 2>::Identity();

  Eigen::Matrix2d M = Q1.transpose() * Q2;
  Eigen::JacobiSVD<Eigen::Matrix2d> svd(M);
  Eigen::Vector2d s = svd.singularValues();

  for (int i = 0; i < 2; ++i) {
    if (s[i] > 1.0)
      s[i] = 1.0;
    if (s[i] < -1.0)
      s[i] = -1.0;
  }

  double theta1 = std::acos(s[0]);
  double theta2 = std::acos(s[1]);
  return std::sqrt(theta1 * theta1 + theta2 * theta2);
}

double NormalDistance::operator()(const Datum &ai, const Datum &aj,
                                  const Datum &bi, const Datum &bj) {
  std::ostringstream log;

  log << "NormalDistance::operator() called with:\n"
      << "ai: " << ai.transpose() << "\naj: " << aj.transpose()
      << "\nbi: " << bi.transpose() << "\nbj: " << bj.transpose() << std::endl;

  // Build Graff matrices for the two planes in first frame
  Eigen::Matrix<double, 4, 2> Y_ai = getGraffMatrixFromDatum(ai);
  Eigen::Matrix<double, 4, 2> Y_aj = getGraffMatrixFromDatum(aj);
  log << "Y_ai:\n" << Y_ai << "\nY_aj:\n" << Y_aj << std::endl;

  // Build Graff matrices for the two planes in second frame
  Eigen::Matrix<double, 4, 2> Y_bi = getGraffMatrixFromDatum(bi);
  Eigen::Matrix<double, 4, 2> Y_bj = getGraffMatrixFromDatum(bj);
  log << "Y_bi:\n" << Y_bi << "\nY_bj:\n" << Y_bj << std::endl;

  double d1 = dGraff(Y_ai, Y_aj);
  double d2 = dGraff(Y_bi, Y_bj);

  log << "d1: " << d1 << ", d2: " << d2 << std::endl;

  // check consistency
  double dp = std::abs(d1 - d2);
  log << "dp: " << dp << std::endl;
  // log_to_file(log.str());

  if (dp < params_.epsp) {
    const double sp = std::exp(-0.5 * dp * dp / (params_.sigp * params_.sigp));
    return sp * sp;
  } else {
    return 0.0;
  }
}

} // namespace invariants
} // namespace clipper
