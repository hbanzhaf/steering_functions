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

#include "steering_functions/dubins_state_space/dubins_state_space.hpp"

namespace
{
const double twopi = 2. * M_PI;
const double DUBINS_EPS = 1e-6;
const double DUBINS_ZERO = -1e-9;

inline int sgn(double x)
{
  return (((x) < 0) ? -1 : 1);
}

inline double mod2pi(double x)
{
  if (x < 0 && x > DUBINS_ZERO)
    return 0;
  return x - twopi * floor(x / twopi);
}

Dubins_State_Space::Dubins_Path dubinsLSL(double d, double alpha, double beta)
{
  double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
  double tmp = 2. + d * d - 2. * (ca * cb + sa * sb - d * (sa - sb));
  if (tmp >= DUBINS_ZERO)
  {
    double theta = atan2(cb - ca, d + sa - sb);
    double t = mod2pi(-alpha + theta);
    double p = sqrt(max(tmp, 0.));
    double q = mod2pi(beta - theta);
    assert(fabs(p * cos(alpha + t) - sa + sb - d) < DUBINS_EPS);
    assert(fabs(p * sin(alpha + t) + ca - cb) < DUBINS_EPS);
    assert(mod2pi(alpha + t + q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
    return Dubins_State_Space::Dubins_Path(Dubins_State_Space::dubins_path_type[0], t, p, q);
  }
  return Dubins_State_Space::Dubins_Path();
}

Dubins_State_Space::Dubins_Path dubinsRSR(double d, double alpha, double beta)
{
  double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
  double tmp = 2. + d * d - 2. * (ca * cb + sa * sb - d * (sb - sa));
  if (tmp >= DUBINS_ZERO)
  {
    double theta = atan2(ca - cb, d - sa + sb);
    double t = mod2pi(alpha - theta);
    double p = sqrt(max(tmp, 0.));
    double q = mod2pi(-beta + theta);
    assert(fabs(p * cos(alpha - t) + sa - sb - d) < DUBINS_EPS);
    assert(fabs(p * sin(alpha - t) - ca + cb) < DUBINS_EPS);
    assert(mod2pi(alpha - t - q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
    return Dubins_State_Space::Dubins_Path(Dubins_State_Space::dubins_path_type[1], t, p, q);
  }
  return Dubins_State_Space::Dubins_Path();
}

Dubins_State_Space::Dubins_Path dubinsRSL(double d, double alpha, double beta)
{
  double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
  double tmp = d * d - 2. + 2. * (ca * cb + sa * sb - d * (sa + sb));
  if (tmp >= DUBINS_ZERO)
  {
    double p = sqrt(max(tmp, 0.));
    double theta = atan2(ca + cb, d - sa - sb) - atan2(2., p);
    double t = mod2pi(alpha - theta);
    double q = mod2pi(beta - theta);
    assert(fabs(p * cos(alpha - t) - 2. * sin(alpha - t) + sa + sb - d) < DUBINS_EPS);
    assert(fabs(p * sin(alpha - t) + 2. * cos(alpha - t) - ca - cb) < DUBINS_EPS);
    assert(mod2pi(alpha - t + q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
    return Dubins_State_Space::Dubins_Path(Dubins_State_Space::dubins_path_type[2], t, p, q);
  }
  return Dubins_State_Space::Dubins_Path();
}

Dubins_State_Space::Dubins_Path dubinsLSR(double d, double alpha, double beta)
{
  double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
  double tmp = -2. + d * d + 2. * (ca * cb + sa * sb + d * (sa + sb));
  if (tmp >= DUBINS_ZERO)
  {
    double p = sqrt(max(tmp, 0.));
    double theta = atan2(-ca - cb, d + sa + sb) - atan2(-2., p);
    double t = mod2pi(-alpha + theta);
    double q = mod2pi(-beta + theta);
    assert(fabs(p * cos(alpha + t) + 2. * sin(alpha + t) - sa - sb - d) < DUBINS_EPS);
    assert(fabs(p * sin(alpha + t) - 2. * cos(alpha + t) + ca + cb) < DUBINS_EPS);
    assert(mod2pi(alpha + t - q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
    return Dubins_State_Space::Dubins_Path(Dubins_State_Space::dubins_path_type[3], t, p, q);
  }
  return Dubins_State_Space::Dubins_Path();
}

Dubins_State_Space::Dubins_Path dubinsRLR(double d, double alpha, double beta)
{
  double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
  double tmp = .125 * (6. - d * d + 2. * (ca * cb + sa * sb + d * (sa - sb)));
  if (fabs(tmp) < 1.)
  {
    double p = twopi - acos(tmp);
    double theta = atan2(ca - cb, d - sa + sb);
    double t = mod2pi(alpha - theta + .5 * p);
    double q = mod2pi(alpha - beta - t + p);
    assert(fabs(2. * sin(alpha - t + p) - 2. * sin(alpha - t) - d + sa - sb) < DUBINS_EPS);
    assert(fabs(-2. * cos(alpha - t + p) + 2. * cos(alpha - t) - ca + cb) < DUBINS_EPS);
    assert(mod2pi(alpha - t + p - q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
    return Dubins_State_Space::Dubins_Path(Dubins_State_Space::dubins_path_type[4], t, p, q);
  }
  return Dubins_State_Space::Dubins_Path();
}

Dubins_State_Space::Dubins_Path dubinsLRL(double d, double alpha, double beta)
{
  double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
  double tmp = .125 * (6. - d * d + 2. * (ca * cb + sa * sb - d * (sa - sb)));
  if (fabs(tmp) < 1.)
  {
    double p = twopi - acos(tmp);
    double theta = atan2(-ca + cb, d + sa - sb);
    double t = mod2pi(-alpha + theta + .5 * p);
    double q = mod2pi(beta - alpha - t + p);
    assert(fabs(-2. * sin(alpha + t - p) + 2. * sin(alpha + t) - d - sa + sb) < DUBINS_EPS);
    assert(fabs(2. * cos(alpha + t - p) - 2. * cos(alpha + t) + ca - cb) < DUBINS_EPS);
    assert(mod2pi(alpha + t - p + q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
    return Dubins_State_Space::Dubins_Path(Dubins_State_Space::dubins_path_type[5], t, p, q);
  }
  return Dubins_State_Space::Dubins_Path();
}

Dubins_State_Space::Dubins_Path dubins(double d, double alpha, double beta)
{
  if (d < DUBINS_EPS && fabs(alpha - beta) < DUBINS_EPS)
    return Dubins_State_Space::Dubins_Path(Dubins_State_Space::dubins_path_type[0], 0, d, 0);

  Dubins_State_Space::Dubins_Path path(dubinsLSL(d, alpha, beta)), tmp(dubinsRSR(d, alpha, beta));
  double len, minLength = path.length();

  if ((len = tmp.length()) < minLength)
  {
    minLength = len;
    path = tmp;
  }
  tmp = dubinsRSL(d, alpha, beta);
  if ((len = tmp.length()) < minLength)
  {
    minLength = len;
    path = tmp;
  }
  tmp = dubinsLSR(d, alpha, beta);
  if ((len = tmp.length()) < minLength)
  {
    minLength = len;
    path = tmp;
  }
  tmp = dubinsRLR(d, alpha, beta);
  if ((len = tmp.length()) < minLength)
  {
    minLength = len;
    path = tmp;
  }
  tmp = dubinsLRL(d, alpha, beta);
  if ((tmp.length()) < minLength)
    path = tmp;
  return path;
}
}

const Dubins_State_Space::Dubins_Path_Segment_Type Dubins_State_Space::dubins_path_type[6][3] = {
  { DUBINS_LEFT, DUBINS_STRAIGHT, DUBINS_LEFT },  { DUBINS_RIGHT, DUBINS_STRAIGHT, DUBINS_RIGHT },
  { DUBINS_RIGHT, DUBINS_STRAIGHT, DUBINS_LEFT }, { DUBINS_LEFT, DUBINS_STRAIGHT, DUBINS_RIGHT },
  { DUBINS_RIGHT, DUBINS_LEFT, DUBINS_RIGHT },    { DUBINS_LEFT, DUBINS_RIGHT, DUBINS_LEFT }
};

Dubins_State_Space::Dubins_Path Dubins_State_Space::dubins(const State &state1, const State &state2) const
{
  double dx = state2.x - state1.x, dy = state2.y - state1.y, th = atan2(dy, dx), d = sqrt(dx * dx + dy * dy) * kappa_;
  double alpha = mod2pi(state1.theta - th), beta = mod2pi(state2.theta - th);
  return ::dubins(d, alpha, beta);
}

double Dubins_State_Space::get_distance(const State &state1, const State &state2) const
{
  if (forwards_)
    return kappa_inv_ * this->dubins(state1, state2).length();
  else
    return kappa_inv_ * this->dubins(state2, state1).length();
}

vector<Control> Dubins_State_Space::get_controls(const State &state1, const State &state2) const
{
  vector<Control> dubins_controls;
  dubins_controls.reserve(3);
  Dubins_State_Space::Dubins_Path path;
  if (forwards_)
    path = this->dubins(state1, state2);
  else
    path = this->dubins(state2, state1);
  for (unsigned int i = 0; i < 3; ++i)
  {
    Control control;
    switch (path.type_[i])
    {
      case DUBINS_LEFT:
        control.delta_s = kappa_inv_ * path.length_[i];
        control.kappa = kappa_;
        control.sigma = 0.0;
        break;
      case DUBINS_STRAIGHT:
        control.delta_s = kappa_inv_ * path.length_[i];
        control.kappa = 0.0;
        control.sigma = 0.0;
        break;
      case DUBINS_RIGHT:
        control.delta_s = kappa_inv_ * path.length_[i];
        control.kappa = -kappa_;
        control.sigma = 0.0;
        break;
    }
    dubins_controls.push_back(control);
  }
  // reverse controls
  if (!forwards_)
  {
    reverse(dubins_controls.begin(), dubins_controls.end());
    for (auto &control : dubins_controls)
      control.delta_s = -control.delta_s;
  }
  return dubins_controls;
}

vector<State> Dubins_State_Space::get_path(const State &state1, const State &state2) const
{
  vector<Control> dubins_controls = get_controls(state1, state2);
  return integrate(state1, dubins_controls);
}

vector<State> Dubins_State_Space::integrate(const State &state, const vector<Control> &controls) const
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

State Dubins_State_Space::interpolate(const State &state, const vector<Control> &controls, double t) const
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

inline State Dubins_State_Space::integrate_ODE(const State &state, const Control &control,
                                               double integration_step) const
{
  State state_next;
  double kappa(control.kappa);
  double d(sgn(control.delta_s));
  if (kappa != 0.0)
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
