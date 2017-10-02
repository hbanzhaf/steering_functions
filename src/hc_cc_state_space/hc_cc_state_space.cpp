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
  double x_i, y_i, theta_i, kappa_i;
  if (length > get_epsilon())
  {
    end_of_clothoid(0, 0, 0, 0, sigma, true, length, &x_i, &y_i, &theta_i, &kappa_i);
  }
  else
  {
    x_i = 0;
    y_i = 0;
    theta_i = 0;
    kappa_i = kappa;
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
  vector<Control> controls = this->get_controls(state1, state2);
  return this->forward_euler(state1, controls);
}

vector<State> HC_CC_State_Space::forward_euler(const State &state, const vector<Control> &controls) const
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
    double sigma(control.sigma);
    double d(sgn(delta_s));
    double s(0.0);
    double integration_step(0.0);
    // push_back current state if curvature discontinuity
    if (fabs(kappa - state_curr.kappa) > get_epsilon())
    {
      state_curr.kappa = kappa;
      state_curr.d = d;
      path.push_back(state_curr);
    }

    while (s < abs_delta_s)
    {
      // get integration step
      s += discretization_;
      if (s > abs_delta_s)
      {
        integration_step = discretization_ - (s - abs_delta_s);
        s = abs_delta_s;
      }
      else
      {
        integration_step = discretization_;
      }
      // forward euler
      if (sigma != 0.0)
      {
        end_of_clothoid(state_curr.x, state_curr.y, state_curr.theta, state_curr.kappa, sigma, d > 0, integration_step,
                        &state_next.x, &state_next.y, &state_next.theta, &state_next.kappa);
        state_next.d = d;
      }
      else
      {
        if (kappa != 0.0)
        {
          state_next.x = state_curr.x +
                         (1 / kappa) * (-sin(state_curr.theta) + sin(state_curr.theta + d * integration_step * kappa));
          state_next.y = state_curr.y +
                         (1 / kappa) * (cos(state_curr.theta) - cos(state_curr.theta + d * integration_step * kappa));
          state_next.theta = pify(state_curr.theta + d * integration_step * kappa);
          state_next.kappa = kappa;
          state_next.d = d;
        }
        else
        {
          state_next.x = state_curr.x + d * integration_step * cos(state_curr.theta);
          state_next.y = state_curr.y + d * integration_step * sin(state_curr.theta);
          state_next.theta = state_curr.theta;
          state_next.kappa = kappa;
          state_next.d = d;
        }
      }
      path.push_back(state_next);
      state_curr = state_next;
    }
  }
  return path;
}
