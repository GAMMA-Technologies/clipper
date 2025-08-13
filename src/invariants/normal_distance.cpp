/**
 * @file pointnormal_distance.cpp
 * @brief Pairwise geometric invariant between point-normals (e.g., planes)
 * @author Parker Lusk <plusk@mit.edu>
 * @date 15 May 2021
 */

#include "clipper/invariants/normal_distance.h"
#include <fstream>

namespace clipper {
namespace invariants {

inline void NormalDistance::log_to_file(const std::string &message) {
  std::ofstream log(filepath_ + "/gamma_slam_cpp_debug_log.txt", std::ios::app);
  if (log.is_open()) {
    log << message << std::endl;
  }
}

double NormalDistance::operator()(const Datum &ai, const Datum &aj,
                                  const Datum &bi, const Datum &bj) {
  std::ostringstream log;
  log << "ai: " << ai.transpose() << "\n"
      << "aj: " << aj.transpose() << "\n"
      << "bi: " << bi.transpose() << "\n"
      << "bj: " << bj.transpose() << "\n";

  // point distance
  const double d_ai = -1 * ai.head<3>().dot(ai.tail<3>());
  const double d_aj = -1 * aj.head<3>().dot(aj.tail<3>());

  log << "d_ai: " << d_ai << "\n"
      << "d_aj: " << d_aj << "\n";

  const double d_bi = -1 * bi.head<3>().dot(bi.tail<3>());
  const double d_bj = -1 * bj.head<3>().dot(bj.tail<3>());

  log << "d_bi: " << d_bi << "\n"
      << "d_bj: " << d_bj << "\n";

  const double l1 = std::abs(d_ai - d_aj);
  const double l2 = std::abs(d_bi - d_bj);

  log << "l1: " << l1 << "\n"
      << "l2: " << l2 << "\n";

  // normal distance
  const double alpha1 = std::acos(ai.tail<3>().transpose() * aj.tail<3>());
  const double alpha2 = std::acos(bi.tail<3>().transpose() * bj.tail<3>());

  log << "alpha1: " << alpha1 << "\n"
      << "alpha2: " << alpha2 << "\n";

  // check consistency
  const double dp = std::abs(l1 - l2);
  const double dn = std::abs(alpha1 - alpha2);

  log << "dp: " << dp << "\n"
      << "dn: " << dn << "\n";

  log_to_file(log.str());

  if (dp < params_.epsp && dn < params_.epsn) {
    const double sp = std::exp(-0.5 * dp * dp / (params_.sigp * params_.sigp));
    const double sn = std::exp(-0.5 * dn * dn / (params_.sign * params_.sign));
    return sp * sn;
  } else {
    return 0.0;
  }
}

} // namespace invariants
} // namespace clipper
