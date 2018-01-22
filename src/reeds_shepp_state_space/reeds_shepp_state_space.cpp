/*********************************************************************
*  Copyright (c) 2017 - for information on the respective copyright
*  owner see the NOTICE file and/or the repository
*
*      https://github.com/hbanzhaf/steering_functions.git
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*      http://www.apache.org/licenses/LICENSE-2.0
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
*  implied. See the License for the specific language governing
*  permissions and limitations under the License.

*  This source code is derived from the Open Motion Planning Library
*  (OMPL) V 1.3.1 (https://github.com/ompl/ompl).
*  Copyright (c) 2010, Rice University, licensed under the BSD license,
*  cf. 3rd-party-licenses.txt file in the root directory of this source
*  tree.
**********************************************************************/

/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2010, Rice University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Rice University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

#include "steering_functions/reeds_shepp_state_space/reeds_shepp_state_space.hpp"

namespace
{
// The comments, variable names, etc. use the nomenclature from the Reeds & Shepp paper.
const double pi = M_PI;
const double twopi = 2. * pi;
const double RS_EPS = 1e-6;
const double ZERO = 10 * numeric_limits<double>::epsilon();

inline int sgn(double x)
{
  return (((x) < 0) ? -1 : 1);
}

inline double mod2pi(double x)
{
  double v = fmod(x, twopi);
  if (v < -pi)
    v += twopi;
  else if (v > pi)
    v -= twopi;
  return v;
}
inline void polar(double x, double y, double &r, double &theta)
{
  r = sqrt(x * x + y * y);
  theta = atan2(y, x);
}
inline void tauOmega(double u, double v, double xi, double eta, double phi, double &tau, double &omega)
{
  double delta = mod2pi(u - v), A = sin(u) - sin(delta), B = cos(u) - cos(delta) - 1.;
  double t1 = atan2(eta * A - xi * B, xi * A + eta * B), t2 = 2. * (cos(delta) - cos(v) - cos(u)) + 3;
  tau = (t2 < 0) ? mod2pi(t1 + pi) : mod2pi(t1);
  omega = mod2pi(tau - u + v - phi);
}

// formula 8.1 in Reeds-Shepp paper
inline bool LpSpLp(double x, double y, double phi, double &t, double &u, double &v)
{
  polar(x - sin(phi), y - 1. + cos(phi), u, t);
  if (t >= -ZERO)
  {
    v = mod2pi(phi - t);
    if (v >= -ZERO)
    {
      assert(fabs(u * cos(t) + sin(phi) - x) < RS_EPS);
      assert(fabs(u * sin(t) - cos(phi) + 1 - y) < RS_EPS);
      assert(fabs(mod2pi(t + v - phi)) < RS_EPS);
      return true;
    }
  }
  return false;
}
// formula 8.2
inline bool LpSpRp(double x, double y, double phi, double &t, double &u, double &v)
{
  double t1, u1;
  polar(x + sin(phi), y - 1. - cos(phi), u1, t1);
  u1 = u1 * u1;
  if (u1 >= 4.)
  {
    double theta;
    u = sqrt(u1 - 4.);
    theta = atan2(2., u);
    t = mod2pi(t1 + theta);
    v = mod2pi(t - phi);
    assert(fabs(2 * sin(t) + u * cos(t) - sin(phi) - x) < RS_EPS);
    assert(fabs(-2 * cos(t) + u * sin(t) + cos(phi) + 1 - y) < RS_EPS);
    assert(fabs(mod2pi(t - v - phi)) < RS_EPS);
    return t >= -ZERO && v >= -ZERO;
  }
  return false;
}
void CSC(double x, double y, double phi, Reeds_Shepp_State_Space::Reeds_Shepp_Path &path)
{
  double t, u, v, Lmin = path.length(), L;
  if (LpSpLp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[14], t, u, v);
    Lmin = L;
  }
  if (LpSpLp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[14], -t, -u, -v);
    Lmin = L;
  }
  if (LpSpLp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[15], t, u, v);
    Lmin = L;
  }
  if (LpSpLp(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[15], -t, -u, -v);
    Lmin = L;
  }
  if (LpSpRp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[12], t, u, v);
    Lmin = L;
  }
  if (LpSpRp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[12], -t, -u, -v);
    Lmin = L;
  }
  if (LpSpRp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[13], t, u, v);
    Lmin = L;
  }
  if (LpSpRp(-x, -y, phi, t, u, v) && Lmin > (fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[13], -t, -u, -v);
}
// formula 8.3 / 8.4  *** TYPO IN PAPER ***
inline bool LpRmL(double x, double y, double phi, double &t, double &u, double &v)
{
  double xi = x - sin(phi), eta = y - 1. + cos(phi), u1, theta;
  polar(xi, eta, u1, theta);
  if (u1 <= 4.)
  {
    u = -2. * asin(.25 * u1);
    t = mod2pi(theta + .5 * u + pi);
    v = mod2pi(phi - t + u);
    assert(fabs(2 * (sin(t) - sin(t - u)) + sin(phi) - x) < RS_EPS);
    assert(fabs(2 * (-cos(t) + cos(t - u)) - cos(phi) + 1 - y) < RS_EPS);
    assert(fabs(mod2pi(t - u + v - phi)) < RS_EPS);
    return t >= -ZERO && u <= ZERO;
  }
  return false;
}
void CCC(double x, double y, double phi, Reeds_Shepp_State_Space::Reeds_Shepp_Path &path)
{
  double t, u, v, Lmin = path.length(), L;
  if (LpRmL(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[0], t, u, v);
    Lmin = L;
  }
  if (LpRmL(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[0], -t, -u, -v);
    Lmin = L;
  }
  if (LpRmL(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[1], t, u, v);
    Lmin = L;
  }
  if (LpRmL(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[1], -t, -u, -v);
    Lmin = L;
  }

  // backwards
  double xb = x * cos(phi) + y * sin(phi), yb = x * sin(phi) - y * cos(phi);
  if (LpRmL(xb, yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[0], v, u, t);
    Lmin = L;
  }
  if (LpRmL(-xb, yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[0], -v, -u, -t);
    Lmin = L;
  }
  if (LpRmL(xb, -yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[1], v, u, t);
    Lmin = L;
  }
  if (LpRmL(-xb, -yb, phi, t, u, v) && Lmin > (fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[1], -v, -u, -t);
}
// formula 8.7
inline bool LpRupLumRm(double x, double y, double phi, double &t, double &u, double &v)
{
  double xi = x + sin(phi), eta = y - 1. - cos(phi), rho = .25 * (2. + sqrt(xi * xi + eta * eta));
  if (rho <= 1.)
  {
    u = acos(rho);
    tauOmega(u, -u, xi, eta, phi, t, v);
    assert(fabs(2 * (sin(t) - sin(t - u) + sin(t - 2 * u)) - sin(phi) - x) < RS_EPS);
    assert(fabs(2 * (-cos(t) + cos(t - u) - cos(t - 2 * u)) + cos(phi) + 1 - y) < RS_EPS);
    assert(fabs(mod2pi(t - 2 * u - v - phi)) < RS_EPS);
    return t >= -ZERO && v <= ZERO;
  }
  return false;
}
// formula 8.8
inline bool LpRumLumRp(double x, double y, double phi, double &t, double &u, double &v)
{
  double xi = x + sin(phi), eta = y - 1. - cos(phi), rho = (20. - xi * xi - eta * eta) / 16.;
  if (rho >= 0 && rho <= 1)
  {
    u = -acos(rho);
    if (u >= -.5 * pi)
    {
      tauOmega(u, u, xi, eta, phi, t, v);
      assert(fabs(4 * sin(t) - 2 * sin(t - u) - sin(phi) - x) < RS_EPS);
      assert(fabs(-4 * cos(t) + 2 * cos(t - u) + cos(phi) + 1 - y) < RS_EPS);
      assert(fabs(mod2pi(t - v - phi)) < RS_EPS);
      return t >= -ZERO && v >= -ZERO;
    }
  }
  return false;
}
void CCCC(double x, double y, double phi, Reeds_Shepp_State_Space::Reeds_Shepp_Path &path)
{
  double t, u, v, Lmin = path.length(), L;
  if (LpRupLumRm(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[2], t, u, -u, v);
    Lmin = L;
  }
  if (LpRupLumRm(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[2], -t, -u, u, -v);
    Lmin = L;
  }
  if (LpRupLumRm(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[3], t, u, -u, v);
    Lmin = L;
  }
  if (LpRupLumRm(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // timeflip + reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[3], -t, -u, u, -v);
    Lmin = L;
  }

  if (LpRumLumRp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[2], t, u, u, v);
    Lmin = L;
  }
  if (LpRumLumRp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[2], -t, -u, -u, -v);
    Lmin = L;
  }
  if (LpRumLumRp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[3], t, u, u, v);
    Lmin = L;
  }
  if (LpRumLumRp(-x, -y, phi, t, u, v) && Lmin > (fabs(t) + 2. * fabs(u) + fabs(v)))  // timeflip + reflect
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[3], -t, -u, -u, -v);
}
// formula 8.9
inline bool LpRmSmLm(double x, double y, double phi, double &t, double &u, double &v)
{
  double xi = x - sin(phi), eta = y - 1. + cos(phi), rho, theta;
  polar(xi, eta, rho, theta);
  if (rho >= 2.)
  {
    double r = sqrt(rho * rho - 4.);
    u = 2. - r;
    t = mod2pi(theta + atan2(r, -2.));
    v = mod2pi(phi - .5 * pi - t);
    assert(fabs(2 * (sin(t) - cos(t)) - u * sin(t) + sin(phi) - x) < RS_EPS);
    assert(fabs(-2 * (sin(t) + cos(t)) + u * cos(t) - cos(phi) + 1 - y) < RS_EPS);
    assert(fabs(mod2pi(t + pi / 2 + v - phi)) < RS_EPS);
    return t >= -ZERO && u <= ZERO && v <= ZERO;
  }
  return false;
}
// formula 8.10
inline bool LpRmSmRm(double x, double y, double phi, double &t, double &u, double &v)
{
  double xi = x + sin(phi), eta = y - 1. - cos(phi), rho, theta;
  polar(-eta, xi, rho, theta);
  if (rho >= 2.)
  {
    t = theta;
    u = 2. - rho;
    v = mod2pi(t + .5 * pi - phi);
    assert(fabs(2 * sin(t) - cos(t - v) - u * sin(t) - x) < RS_EPS);
    assert(fabs(-2 * cos(t) - sin(t - v) + u * cos(t) + 1 - y) < RS_EPS);
    assert(fabs(mod2pi(t + pi / 2 - v - phi)) < RS_EPS);
    return t >= -ZERO && u <= ZERO && v <= ZERO;
  }
  return false;
}
void CCSC(double x, double y, double phi, Reeds_Shepp_State_Space::Reeds_Shepp_Path &path)
{
  double t, u, v, Lmin = path.length() - .5 * pi, L;
  if (LpRmSmLm(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
  {
    path =
        Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[4], t, -.5 * pi, u, v);
    Lmin = L;
  }
  if (LpRmSmLm(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[4], -t, .5 * pi, -u,
                                                     -v);
    Lmin = L;
  }
  if (LpRmSmLm(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
  {
    path =
        Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[5], t, -.5 * pi, u, v);
    Lmin = L;
  }
  if (LpRmSmLm(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[5], -t, .5 * pi, -u,
                                                     -v);
    Lmin = L;
  }

  if (LpRmSmRm(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
  {
    path =
        Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[8], t, -.5 * pi, u, v);
    Lmin = L;
  }
  if (LpRmSmRm(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[8], -t, .5 * pi, -u,
                                                     -v);
    Lmin = L;
  }
  if (LpRmSmRm(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
  {
    path =
        Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[9], t, -.5 * pi, u, v);
    Lmin = L;
  }
  if (LpRmSmRm(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[9], -t, .5 * pi, -u,
                                                     -v);
    Lmin = L;
  }

  // backwards
  double xb = x * cos(phi) + y * sin(phi), yb = x * sin(phi) - y * cos(phi);
  if (LpRmSmLm(xb, yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
  {
    path =
        Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[6], v, u, -.5 * pi, t);
    Lmin = L;
  }
  if (LpRmSmLm(-xb, yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[6], -v, -u, .5 * pi,
                                                     -t);
    Lmin = L;
  }
  if (LpRmSmLm(xb, -yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
  {
    path =
        Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[7], v, u, -.5 * pi, t);
    Lmin = L;
  }
  if (LpRmSmLm(-xb, -yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[7], -v, -u, .5 * pi,
                                                     -t);
    Lmin = L;
  }

  if (LpRmSmRm(xb, yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[10], v, u, -.5 * pi,
                                                     t);
    Lmin = L;
  }
  if (LpRmSmRm(-xb, yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[10], -v, -u,
                                                     .5 * pi, -t);
    Lmin = L;
  }
  if (LpRmSmRm(xb, -yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[11], v, u, -.5 * pi,
                                                     t);
    Lmin = L;
  }
  if (LpRmSmRm(-xb, -yb, phi, t, u, v) && Lmin > (fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[11], -v, -u,
                                                     .5 * pi, -t);
}
// formula 8.11 *** TYPO IN PAPER ***
inline bool LpRmSLmRp(double x, double y, double phi, double &t, double &u, double &v)
{
  double xi = x + sin(phi), eta = y - 1. - cos(phi), rho, theta;
  polar(xi, eta, rho, theta);
  if (rho >= 2.)
  {
    u = 4. - sqrt(rho * rho - 4.);
    if (u <= ZERO)
    {
      t = mod2pi(atan2((4 - u) * xi - 2 * eta, -2 * xi + (u - 4) * eta));
      v = mod2pi(t - phi);
      assert(fabs(4 * sin(t) - 2 * cos(t) - u * sin(t) - sin(phi) - x) < RS_EPS);
      assert(fabs(-4 * cos(t) - 2 * sin(t) + u * cos(t) + cos(phi) + 1 - y) < RS_EPS);
      assert(fabs(mod2pi(t - v - phi)) < RS_EPS);
      return t >= -ZERO && v >= -ZERO;
    }
  }
  return false;
}
void CCSCC(double x, double y, double phi, Reeds_Shepp_State_Space::Reeds_Shepp_Path &path)
{
  double t, u, v, Lmin = path.length() - pi, L;
  if (LpRmSLmRp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[16], t, -.5 * pi, u,
                                                     -.5 * pi, v);
    Lmin = L;
  }
  if (LpRmSLmRp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[16], -t, .5 * pi,
                                                     -u, .5 * pi, -v);
    Lmin = L;
  }
  if (LpRmSLmRp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
  {
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[17], t, -.5 * pi, u,
                                                     -.5 * pi, v);
    Lmin = L;
  }
  if (LpRmSLmRp(-x, -y, phi, t, u, v) && Lmin > (fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
    path = Reeds_Shepp_State_Space::Reeds_Shepp_Path(Reeds_Shepp_State_Space::reeds_shepp_path_type[17], -t, .5 * pi,
                                                     -u, .5 * pi, -v);
}

Reeds_Shepp_State_Space::Reeds_Shepp_Path reeds_shepp(double x, double y, double phi)
{
  Reeds_Shepp_State_Space::Reeds_Shepp_Path path;
  CSC(x, y, phi, path);
  CCC(x, y, phi, path);
  CCCC(x, y, phi, path);
  CCSC(x, y, phi, path);
  CCSCC(x, y, phi, path);
  return path;
}
}

const Reeds_Shepp_State_Space::Reeds_Shepp_Path_Segment_Type Reeds_Shepp_State_Space::reeds_shepp_path_type[18][5] = {
  { RS_LEFT, RS_RIGHT, RS_LEFT, RS_NOP, RS_NOP },         // 0
  { RS_RIGHT, RS_LEFT, RS_RIGHT, RS_NOP, RS_NOP },        // 1
  { RS_LEFT, RS_RIGHT, RS_LEFT, RS_RIGHT, RS_NOP },       // 2
  { RS_RIGHT, RS_LEFT, RS_RIGHT, RS_LEFT, RS_NOP },       // 3
  { RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_NOP },    // 4
  { RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_NOP },   // 5
  { RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_LEFT, RS_NOP },    // 6
  { RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_RIGHT, RS_NOP },   // 7
  { RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP },   // 8
  { RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_NOP },    // 9
  { RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_LEFT, RS_NOP },   // 10
  { RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_RIGHT, RS_NOP },    // 11
  { RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_NOP, RS_NOP },     // 12
  { RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_NOP, RS_NOP },     // 13
  { RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_NOP, RS_NOP },      // 14
  { RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP, RS_NOP },    // 15
  { RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_RIGHT },  // 16
  { RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_LEFT }   // 17
};

Reeds_Shepp_State_Space::Reeds_Shepp_Path::Reeds_Shepp_Path(const Reeds_Shepp_Path_Segment_Type *type, double t,
                                                            double u, double v, double w, double x)
  : type_(type)
{
  length_[0] = t;
  length_[1] = u;
  length_[2] = v;
  length_[3] = w;
  length_[4] = x;
  total_length_ = fabs(t) + fabs(u) + fabs(v) + fabs(w) + fabs(x);
}

Reeds_Shepp_State_Space::Reeds_Shepp_Path Reeds_Shepp_State_Space::reeds_shepp(const State &state1,
                                                                               const State &state2) const
{
  double dx = state2.x - state1.x, dy = state2.y - state1.y, dth = state2.theta - state1.theta;
  double c = cos(state1.theta), s = sin(state1.theta);
  double x = c * dx + s * dy, y = -s * dx + c * dy;
  return ::reeds_shepp(x * kappa_, y * kappa_, dth);
}

double Reeds_Shepp_State_Space::get_distance(const State &state1, const State &state2) const
{
  return kappa_inv_ * this->reeds_shepp(state1, state2).length();
}

vector<Control> Reeds_Shepp_State_Space::get_controls(const State &state1, const State &state2) const
{
  vector<Control> controls;
  controls.reserve(5);
  Reeds_Shepp_State_Space::Reeds_Shepp_Path path = this->reeds_shepp(state1, state2);
  for (unsigned int i = 0; i < 5; ++i)
  {
    Control control;
    switch (path.type_[i])
    {
      case RS_NOP:
        return controls;
      case RS_LEFT:
        control.delta_s = kappa_inv_ * path.length_[i];
        control.kappa = kappa_;
        control.sigma = 0.0;
        break;
      case RS_RIGHT:
        control.delta_s = kappa_inv_ * path.length_[i];
        control.kappa = -kappa_;
        control.sigma = 0.0;
        break;
      case RS_STRAIGHT:
        control.delta_s = kappa_inv_ * path.length_[i];
        control.kappa = 0.0;
        control.sigma = 0.0;
        break;
    }
    controls.push_back(control);
  }
  return controls;
}

vector<State> Reeds_Shepp_State_Space::get_path(const State &state1, const State &state2) const
{
  vector<Control> controls = get_controls(state1, state2);
  return integrate(state1, controls);
}

vector<State> Reeds_Shepp_State_Space::integrate(const State &state, const vector<Control> &controls) const
{
  vector<State> path;
  State state_curr, state_next;
  // reserve capacity of path
  int n_states(0);
  for (const auto &control : controls)
  {
    double abs_delta_s(fabs(control.delta_s));
    n_states += ceil(abs_delta_s / discretization_);
  }
  path.reserve(n_states + 1);
  // get first state
  state_curr.x = state.x;
  state_curr.y = state.y;
  state_curr.theta = state.theta;

  for (const auto &control : controls)
  {
    double delta_s(control.delta_s);
    double abs_delta_s(fabs(delta_s));
    double kappa(control.kappa);
    double s_seg(0.0);
    double integration_step(0.0);
    // push_back current state
    state_curr.kappa = kappa;
    state_curr.d = sgn(delta_s);
    path.push_back(state_curr);

    while (s_seg < abs_delta_s)
    {
      // get integration step
      s_seg += discretization_;
      if (s_seg > abs_delta_s)
      {
        integration_step = discretization_ - (s_seg - abs_delta_s);
        s_seg = abs_delta_s;
      }
      else
      {
        integration_step = discretization_;
      }
      state_next = integrate_ODE(state_curr, control, integration_step);
      path.push_back(state_next);
      state_curr = state_next;
    }
  }
  return path;
}

State Reeds_Shepp_State_Space::interpolate(const State &state, const vector<Control> &controls, double t) const
{
  State state_curr, state_next;
  // get first state
  state_curr.x = state.x;
  state_curr.y = state.y;
  state_curr.theta = state.theta;
  state_curr.kappa = controls.front().kappa;
  state_curr.d = sgn(controls.front().delta_s);
  // get arc length at t
  double s_path(0.0);
  double s_inter(0.0);
  for (const auto &control : controls)
  {
    s_path += fabs(control.delta_s);
  }
  if (t <= 0.0)
    return state_curr;
  else if (t > 1.0)
    s_inter = s_path;
  else
    s_inter = t * s_path;

  double s(0.0);
  bool interpolated = false;
  for (const auto &control : controls)
  {
    if (interpolated)
      break;

    double delta_s(control.delta_s);
    double abs_delta_s(fabs(delta_s));
    double s_seg(0.0);
    double integration_step(0.0);

    s += abs_delta_s;
    if (s > s_inter)
    {
      abs_delta_s = abs_delta_s - (s - s_inter);
      interpolated = true;
    }

    while (s_seg < abs_delta_s)
    {
      // get integration step
      s_seg += discretization_;
      if (s_seg > abs_delta_s)
      {
        integration_step = discretization_ - (s_seg - abs_delta_s);
        s_seg = abs_delta_s;
      }
      else
      {
        integration_step = discretization_;
      }
      state_next = integrate_ODE(state_curr, control, integration_step);
      state_curr = state_next;
    }
  }
  return state_curr;
}

inline State Reeds_Shepp_State_Space::integrate_ODE(const State &state, const Control &control,
                                                    double integration_step) const
{
  State state_next;
  double kappa(control.kappa);
  double d(sgn(control.delta_s));
  if (fabs(kappa) > RS_EPS)
  {
    state_next.x = state.x + (1 / kappa) * (-sin(state.theta) + sin(state.theta + d * integration_step * kappa));
    state_next.y = state.y + (1 / kappa) * (cos(state.theta) - cos(state.theta + d * integration_step * kappa));
    state_next.theta = mod2pi(state.theta + d * integration_step * kappa);
    state_next.kappa = kappa;
    state_next.d = d;
  }
  else
  {
    state_next.x = state.x + d * integration_step * cos(state.theta);
    state_next.y = state.y + d * integration_step * sin(state.theta);
    state_next.theta = state.theta;
    state_next.kappa = kappa;
    state_next.d = d;
  }
  return state_next;
}
