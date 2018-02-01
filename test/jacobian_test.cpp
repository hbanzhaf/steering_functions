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
* *********************************************************************/

#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "steering_functions/filter/ekf.hpp"
#include "steering_functions/steering_functions.hpp"
#include "steering_functions/utilities/utilities.hpp"

using namespace std;
using namespace steer;

#define EPS_JACOBI 1e-4                  // [-]
#define EPS_PERTURB 1e-7                 // [-]
#define SAMPLES 1e6                      // [-]
#define OPERATING_REGION_X 20.0          // [m]
#define OPERATING_REGION_Y 20.0          // [m]
#define OPERATING_REGION_THETA 2 * M_PI  // [rad]
#define OPERATING_REGION_KAPPA 2.0       // [1/m]
#define OPERATING_REGION_DELTA_S 0.4     // [m]
#define OPERATING_REGION_SIGMA 2.0       // [1/m^2]
#define random(lower, upper) (rand() * (upper - lower) / RAND_MAX + lower)
#define random_boolean() rand() % 2

typedef Eigen::Matrix<double, 2, 2> Matrix2d;
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 4, 2> Matrix42d;

State get_random_state()
{
  State state;
  state.x = random(-OPERATING_REGION_X / 2.0, OPERATING_REGION_X / 2.0);
  state.y = random(-OPERATING_REGION_Y / 2.0, OPERATING_REGION_Y / 2.0);
  state.theta = random(-OPERATING_REGION_THETA / 2.0, OPERATING_REGION_THETA / 2.0);
  state.kappa = random(-OPERATING_REGION_KAPPA / 2.0, OPERATING_REGION_KAPPA / 2.0);
  state.d = 0.0;
  return state;
}

Control get_random_control()
{
  Control control;
  control.delta_s = random(-OPERATING_REGION_DELTA_S / 2.0, OPERATING_REGION_DELTA_S / 2.0);
  control.sigma = random(-OPERATING_REGION_SIGMA / 2.0, OPERATING_REGION_SIGMA / 2.0);
  return control;
}

State integrate_ODE(const State &state, const Control &control, double integration_step)
{
  State state_next;
  double sigma(control.sigma);
  double d(sgn(control.delta_s));
  if (fabs(sigma) > get_epsilon())
  {
    end_of_clothoid(state.x, state.y, state.theta, state.kappa, sigma, d, integration_step, &state_next.x,
                    &state_next.y, &state_next.theta, &state_next.kappa);
  }
  else
  {
    if (fabs(state.kappa) > get_epsilon())
    {
      end_of_circular_arc(state.x, state.y, state.theta, state.kappa, d, integration_step, &state_next.x, &state_next.y,
                          &state_next.theta);
    }
    else
    {
      end_of_straight_line(state.x, state.y, state.theta, d, integration_step, &state_next.x, &state_next.y);
    }
  }
  state_next.theta =
      pify(state.theta + state.kappa * d * integration_step + 0.5 * sigma * d * integration_step * integration_step);
  state_next.kappa = state.kappa + sigma * integration_step;
  state_next.d = d;

  return state_next;
}

Matrix4d get_num_motion_jacobi_x(const State &state, const Control &control, double integration_step)
{
  Matrix4d F_x;
  Matrix4d perturb(EPS_PERTURB * Matrix4d::Identity());

  State f_x = integrate_ODE(state, control, integration_step);
  for (int i = 0; i < F_x.cols(); ++i)
  {
    // perturb state
    State state_perturb;
    state_perturb.x = state.x + perturb(0, i);
    state_perturb.y = state.y + perturb(1, i);
    state_perturb.theta = state.theta + perturb(2, i);
    state_perturb.kappa = state.kappa + perturb(3, i);

    State f_x_perturb = integrate_ODE(state_perturb, control, integration_step);
    F_x(0, i) = (f_x_perturb.x - f_x.x) / perturb(i, i);
    F_x(1, i) = (f_x_perturb.y - f_x.y) / perturb(i, i);
    F_x(2, i) = (f_x_perturb.theta - f_x.theta) / perturb(i, i);
    F_x(3, i) = (f_x_perturb.kappa - f_x.kappa) / perturb(i, i);
  }
  return F_x;
}

Matrix42d get_num_motion_jacobi_u(const State &state, const Control &control, double integration_step)
{
  Matrix42d F_u;
  Matrix2d perturb(EPS_PERTURB * Matrix2d::Identity());

  State f_u = integrate_ODE(state, control, integration_step);
  for (int i = 0; i < F_u.cols(); ++i)
  {
    // perturb control
    Control control_perturb;
    control_perturb.delta_s = control.delta_s + perturb(0, i);
    control_perturb.sigma = control.sigma + perturb(1, i);
    double integration_step_perturb = fabs(control_perturb.delta_s);

    State f_u_perturb = integrate_ODE(state, control_perturb, integration_step_perturb);
    F_u(0, i) = (f_u_perturb.x - f_u.x) / perturb(i, i);
    F_u(1, i) = (f_u_perturb.y - f_u.y) / perturb(i, i);
    F_u(2, i) = (f_u_perturb.theta - f_u.theta) / perturb(i, i);
    F_u(3, i) = (f_u_perturb.kappa - f_u.kappa) / perturb(i, i);
  }
  return F_u;
}

TEST(Jacobian, F_x)
{
  EKF ekf;
  for (int i = 0; i < SAMPLES; i++)
  {
    State state = get_random_state();
    Control control = get_random_control();

    state.kappa = random_boolean() * state.kappa;
    control.sigma = random_boolean() * control.sigma;
    double integration_step = fabs(control.delta_s);

    Matrix4d F_x_ana(Matrix4d::Zero());
    Matrix42d F_u_ana(Matrix42d::Zero());
    ekf.get_motion_jacobi(state, control, integration_step, F_x_ana, F_u_ana);
    Matrix4d F_x_num = get_num_motion_jacobi_x(state, control, integration_step);
    for (int i = 0; i < F_x_ana.rows(); ++i)
    {
      for (int j = 0; j < F_x_ana.cols(); ++j)
      {
        EXPECT_LE(fabs(F_x_ana(i, j) - F_x_num(i, j)), EPS_JACOBI);
      }
    }
  }
}

TEST(Jacobian, F_u)
{
  EKF ekf;
  for (int i = 0; i < SAMPLES; i++)
  {
    State state = get_random_state();
    Control control = get_random_control();

    state.kappa = random_boolean() * state.kappa;
    control.sigma = random_boolean() * control.sigma;
    double integration_step = fabs(control.delta_s);

    Matrix4d F_x_ana(Matrix4d::Zero());
    Matrix42d F_u_ana(Matrix42d::Zero());
    ekf.get_motion_jacobi(state, control, integration_step, F_x_ana, F_u_ana);
    Matrix42d F_u_num = get_num_motion_jacobi_u(state, control, integration_step);
    for (int i = 0; i < F_u_ana.rows(); ++i)
    {
      for (int j = 0; j < F_u_ana.cols(); ++j)
      {
        EXPECT_LE(fabs(F_u_ana(i, j) - F_u_num(i, j)), EPS_JACOBI);
      }
    }
  }
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
