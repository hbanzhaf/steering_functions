/*********************************************************************
*  Copyright (c) 2017 Robert Bosch GmbH.
*  All rights reserved.
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*      http://www.apache.org/licenses/LICENSE-2.0
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
***********************************************************************/

#include "steering_functions/hc_cc_state_space/hc_cc_state_space.hpp"

HC_CC_State_Space::HC_CC_State_Space(double kappa, double sigma, double discretization)
  : kappa_(kappa), sigma_(sigma), discretization_(discretization)
{
  // intermediate configuration after first clothoid
  double length = kappa / sigma;
  double x_i, y_i, theta_i;
  if (length > get_epsilon())
  {
    double kappa_i;
    end_of_clothoid(0, 0, 0, 0, sigma, 1, length, &x_i, &y_i, &theta_i, &kappa_i);
  }
  else
  {
    x_i = 0;
    y_i = 0;
    theta_i = 0;
  }
  // radius
  double xc, yc;
  xc = x_i - sin(theta_i) / kappa;
  yc = y_i + cos(theta_i) / kappa;
  double radius = point_distance(xc, yc, 0.0, 0.0);
  // mu
  double mu = atan(fabs(xc / yc));
  double sin_mu = sin(mu);
  double cos_mu = cos(mu);
  // delta_min
  double delta_min = twopify(pow(kappa, 2) / sigma);
  // assign
  hc_cc_circle_param_.set_param(kappa, sigma, radius, mu, sin_mu, cos_mu, delta_min);
}

vector<State> HC_CC_State_Space::get_path(const State &state1, const State &state2) const
{
  vector<Control> controls = get_controls(state1, state2);
  return integrate(state1, controls);
}

vector<State> HC_CC_State_Space::integrate(const State &state, const vector<Control> &controls) const
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
  // push back first state
  state_curr.x = state.x;
  state_curr.y = state.y;
  state_curr.theta = state.theta;
  state_curr.kappa = controls.front().kappa;
  state_curr.d = sgn(controls.front().delta_s);
  path.push_back(state_curr);

  for (const auto &control : controls)
  {
    double delta_s(control.delta_s);
    double abs_delta_s(fabs(delta_s));
    double kappa(control.kappa);
    double d(sgn(delta_s));
    double s_seg(0.0);
    double integration_step(0.0);
    // push_back current state if curvature discontinuity
    if (fabs(kappa - state_curr.kappa) > get_epsilon())
    {
      state_curr.kappa = kappa;
      state_curr.d = d;
      path.push_back(state_curr);
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
      path.push_back(state_next);
      state_curr = state_next;
    }
  }
  return path;
}

State HC_CC_State_Space::interpolate(const State &state, const vector<Control> &controls, double t) const
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

inline State HC_CC_State_Space::integrate_ODE(const State &state, const Control &control, double integration_step) const
{
  State state_next;
  double kappa(control.kappa);
  double sigma(control.sigma);
  double d(sgn(control.delta_s));
  if (fabs(sigma) > get_epsilon())
  {
    end_of_clothoid(state.x, state.y, state.theta, state.kappa, sigma, d, integration_step, &state_next.x,
                    &state_next.y, &state_next.theta, &state_next.kappa);
    state_next.d = d;
  }
  else
  {
    if (fabs(kappa) > get_epsilon())
    {
      end_of_circular_arc(state.x, state.y, state.theta, kappa, d, integration_step, &state_next.x, &state_next.y,
                          &state_next.theta);
      state_next.kappa = kappa;
      state_next.d = d;
    }
    else
    {
      end_of_straight_line(state.x, state.y, state.theta, d, integration_step, &state_next.x, &state_next.y,
                           &state_next.theta);
      state_next.kappa = kappa;
      state_next.d = d;
    }
  }
  return state_next;
}
