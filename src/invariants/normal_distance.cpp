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


double NormalDistance::operator()(const Datum& ai, const Datum& aj,
                                       const Datum& bi, const Datum& bj)
{

  // point distance
  const double d_ai = -1 * ai.head<3>().dot(ai.tail<3>());
  const double d_aj = -1 * aj.head<3>().dot(aj.tail<3>());
 
  const double d_bi = -1 * bi.head<3>().dot(bi.tail<3>());
  const double d_bj = -1 * bj.head<3>().dot(bj.tail<3>());

  const double l1 = std::abs(d_ai - d_aj);
  const double l2 = std::abs(d_bi - d_bj);

  // normal distance
  const double alpha1 = std::acos(ai.tail<3>().transpose() * aj.tail<3>());
  const double alpha2 = std::acos(bi.tail<3>().transpose() * bj.tail<3>());

  // check consistency
  const double dp = std::abs(l1 - l2);
  const double dn = std::abs(alpha1 - alpha2);
 
  if (dp < params_.epsp && dn < params_.epsn) {
    const double sp = std::exp(-0.5*dp*dp/(params_.sigp*params_.sigp));
    const double sn = std::exp(-0.5*dn*dn/(params_.sign*params_.sign));
    return sp * sn;
  } else {
    return 0.0;
  }
}

} // ns invariants
} // ns clipper
