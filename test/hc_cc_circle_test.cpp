/*********************************************************************
*  Copyright (c) 2019 Robert Bosch GmbH.
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

#include "steering_functions/steering_functions.hpp"
#include "steering_functions/hc_cc_state_space/configuration.hpp"
#include "steering_functions/hc_cc_state_space/hc_cc_circle.hpp"
#include "steering_functions/hc_cc_state_space/hc_cc_state_space.hpp"
#include "steering_functions/hc_cc_state_space/paths.hpp"
#include "steering_functions/utilities/utilities.hpp"

#define EPS_DISTANCE 0.01    // [m]
#define EPS_KAPPA 1e-6       // [1/m]
#define EPS_SIGMA 1e-6       // [1/m^2]
#define KAPPA 1.0            // [1/m]
#define DISCRETIZATION 0.05  // [m]

using namespace std;
using namespace steering;

class Test_HC_CC_State_Space : public HC_CC_State_Space
{
public:
  Test_HC_CC_State_Space(double kappa, double sigma, double discretization)
    : HC_CC_State_Space(kappa, sigma, discretization)
  {
  }

  vector<Control> get_controls(const State &state1, const State &state2) const
  {
    throw runtime_error("Test_HC_CC_State_Space.get_controls() not implemented");
  }

  HC_CC_Circle_Param get_hc_cc_circle_param()
  {
    return hc_cc_circle_param_;
  }
};

Configuration get_cc_goal_config(double delta, const HC_CC_Circle_Param &hc_cc_circle_param)
{
  double radius = hc_cc_circle_param.radius;
  double mu = hc_cc_circle_param.mu;
  double xc = radius * sin(mu);
  double yc = radius * cos(mu);

  // goal config of cc-turn given in local frame (aligned with start config) assuming left = forward = true
  Configuration goal;
  goal.x = xc + radius * sin(delta + mu);
  goal.y = yc - radius * cos(delta + mu);
  goal.theta = delta;
  goal.kappa = 0.0;
  return goal;
}

Configuration get_hc_goal_config(double delta, const HC_CC_Circle_Param &hc_cc_circle_param)
{
  double radius = hc_cc_circle_param.radius;
  double mu = hc_cc_circle_param.mu;
  double xc = radius * sin(mu);
  double yc = radius * cos(mu);
  double kappa = hc_cc_circle_param.kappa;

  // goal config of hc-turn given in local frame (aligned with start config) assuming left = forward = true
  Configuration goal;
  goal.x = xc + (1 / kappa) * sin(delta);
  goal.y = yc - (1 / kappa) * cos(delta);
  goal.theta = delta;
  goal.kappa = kappa;
  return goal;
}

vector<double> get_linear_samples(double start, double stop, size_t num)
{
  double step = (stop - start) / static_cast<double>(num - 1);
  vector<double> samples(num);
  double sample = start;
  for (vector<double>::iterator iter = samples.begin(); iter != samples.end(); ++iter, sample += step)
    *iter = sample;
  return samples;
}

double get_distance(const State &state1, const State &state2)
{
  return sqrt(pow(state2.x - state1.x, 2) + pow(state2.y - state1.y, 2));
}

double get_path_length(const vector<State> &path)
{
  double path_length = 0;
  State state1 = path.front();
  for (const auto &state2 : path)
  {
    path_length += get_distance(state1, state2);
    state1 = state2;
  }
  return path_length;
}

Configuration start_config(0.0, 0.0, 0.0, 0.0);                  // do not change
vector<double> sigmas = get_linear_samples(0.02, 200.0, 1000);   // [1/m^2]
vector<double> deltas = get_linear_samples(0.0, 2 * M_PI, 1e3);  // [rad]

TEST(HC_CC_Circle, pathLength)
{
  State start_state;
  start_state.x = start_config.x;
  start_state.y = start_config.y;
  start_state.theta = start_config.theta;
  start_state.kappa = start_config.kappa;
  start_state.d = 0.0;

  for (const auto &sigma : sigmas)
  {
    Test_HC_CC_State_Space hc_cc_ss(KAPPA, sigma, DISCRETIZATION);
    HC_CC_Circle_Param hc_cc_circle_param = hc_cc_ss.get_hc_cc_circle_param();

    HC_CC_Circle reg_hc_cc_circle = HC_CC_Circle(start_config, true, true, true, hc_cc_circle_param);
    HC_CC_Circle irreg_hc_cc_circle = HC_CC_Circle(start_config, true, true, false, hc_cc_circle_param);

    for (const auto &delta : deltas)
    {
      // cc-turn
      Configuration cc_goal_config = get_cc_goal_config(delta, hc_cc_circle_param);

      // regular cc-turn
      vector<Control> reg_cc_controls;
      cc_turn_controls(reg_hc_cc_circle, cc_goal_config, true, reg_cc_controls);
      vector<State> reg_cc_path = hc_cc_ss.integrate(start_state, reg_cc_controls);
      double reg_cc_path_length =
          accumulate(reg_cc_controls.begin(), reg_cc_controls.end(), 0.0,
                     [](double sum, const Control &control) { return sum + fabs(control.delta_s); });
      EXPECT_LT(fabs(reg_cc_path_length - get_path_length(reg_cc_path)), EPS_DISTANCE);

      // irregular cc-turn
      vector<Control> irreg_cc_controls;
      cc_turn_controls(irreg_hc_cc_circle, cc_goal_config, true, irreg_cc_controls);
      vector<State> irreg_cc_path = hc_cc_ss.integrate(start_state, irreg_cc_controls);
      double irreg_cc_path_length =
          accumulate(irreg_cc_controls.begin(), irreg_cc_controls.end(), 0.0,
                     [](double sum, const Control &control) { return sum + fabs(control.delta_s); });
      EXPECT_LT(fabs(irreg_cc_path_length - get_path_length(irreg_cc_path)), EPS_DISTANCE);

      // hc-turn
      Configuration hc_goal_config = get_hc_goal_config(delta, hc_cc_circle_param);

      // regular hc-turn
      vector<Control> reg_hc_controls;
      hc_turn_controls(reg_hc_cc_circle, hc_goal_config, true, reg_hc_controls);
      vector<State> reg_hc_path = hc_cc_ss.integrate(start_state, reg_hc_controls);
      double reg_hc_path_length =
          accumulate(reg_hc_controls.begin(), reg_hc_controls.end(), 0.0,
                     [](double sum, const Control &control) { return sum + fabs(control.delta_s); });
      EXPECT_LT(fabs(reg_hc_path_length - get_path_length(reg_hc_path)), EPS_DISTANCE);

      // irregular hc-turn
      vector<Control> irreg_hc_controls;
      hc_turn_controls(irreg_hc_cc_circle, hc_goal_config, true, irreg_hc_controls);
      vector<State> irreg_hc_path = hc_cc_ss.integrate(start_state, irreg_hc_controls);
      double irreg_hc_path_length =
          accumulate(irreg_hc_controls.begin(), irreg_hc_controls.end(), 0.0,
                     [](double sum, const Control &control) { return sum + fabs(control.delta_s); });
      EXPECT_LT(fabs(irreg_hc_path_length - get_path_length(irreg_hc_path)), EPS_DISTANCE);
    }
  }
}

TEST(HC_CC_Circle, reachingGoal)
{
  State start_state;
  start_state.x = start_config.x;
  start_state.y = start_config.y;
  start_state.theta = start_config.theta;
  start_state.kappa = start_config.kappa;
  start_state.d = 0.0;

  for (const auto &sigma : sigmas)
  {
    Test_HC_CC_State_Space hc_cc_ss(KAPPA, sigma, DISCRETIZATION);
    HC_CC_Circle_Param hc_cc_circle_param = hc_cc_ss.get_hc_cc_circle_param();

    HC_CC_Circle reg_hc_cc_circle = HC_CC_Circle(start_config, true, true, true, hc_cc_circle_param);
    HC_CC_Circle irreg_hc_cc_circle = HC_CC_Circle(start_config, true, true, false, hc_cc_circle_param);

    for (const auto &delta : deltas)
    {
      // cc-turn
      Configuration cc_goal_config = get_cc_goal_config(delta, hc_cc_circle_param);
      State cc_goal_state;
      cc_goal_state.x = cc_goal_config.x;
      cc_goal_state.y = cc_goal_config.y;
      cc_goal_state.theta = cc_goal_config.theta;
      cc_goal_state.kappa = cc_goal_config.kappa;
      cc_goal_state.d = 0.0;

      // regular cc-turn
      vector<Control> reg_cc_controls;
      cc_turn_controls(reg_hc_cc_circle, cc_goal_config, true, reg_cc_controls);
      vector<State> reg_cc_path = hc_cc_ss.integrate(start_state, reg_cc_controls);
      EXPECT_LT(get_distance(cc_goal_state, reg_cc_path.back()), EPS_DISTANCE);

      // irregular cc-turn
      vector<Control> irreg_cc_controls;
      cc_turn_controls(irreg_hc_cc_circle, cc_goal_config, true, irreg_cc_controls);
      vector<State> irreg_cc_path = hc_cc_ss.integrate(start_state, irreg_cc_controls);
      EXPECT_LT(get_distance(cc_goal_state, irreg_cc_path.back()), EPS_DISTANCE);

      // hc-turn
      Configuration hc_goal_config = get_hc_goal_config(delta, hc_cc_circle_param);
      State hc_goal_state;
      hc_goal_state.x = hc_goal_config.x;
      hc_goal_state.y = hc_goal_config.y;
      hc_goal_state.theta = hc_goal_config.theta;
      hc_goal_state.kappa = hc_goal_config.kappa;
      hc_goal_state.d = 0.0;

      // regular hc-turn
      vector<Control> reg_hc_controls;
      hc_turn_controls(reg_hc_cc_circle, hc_goal_config, true, reg_hc_controls);
      vector<State> reg_hc_path = hc_cc_ss.integrate(start_state, reg_hc_controls);
      EXPECT_LT(get_distance(hc_goal_state, reg_hc_path.back()), EPS_DISTANCE);

      // irregular hc-turn
      vector<Control> irreg_hc_controls;
      hc_turn_controls(irreg_hc_cc_circle, hc_goal_config, true, irreg_hc_controls);
      vector<State> irreg_hc_path = hc_cc_ss.integrate(start_state, irreg_hc_controls);
      EXPECT_LT(get_distance(hc_goal_state, irreg_hc_path.back()), EPS_DISTANCE);
    }
  }
}

TEST(HC_CC_Circle, maxSharpness)
{
  for (const auto &sigma : sigmas)
  {
    Test_HC_CC_State_Space hc_cc_ss(KAPPA, sigma, DISCRETIZATION);
    HC_CC_Circle_Param hc_cc_circle_param = hc_cc_ss.get_hc_cc_circle_param();
    HC_CC_Circle hc_cc_circle = HC_CC_Circle(start_config, true, true, true, hc_cc_circle_param);

    for (const auto &delta : deltas)
    {
      // cc-turn
      Configuration goal_config = get_cc_goal_config(delta, hc_cc_circle_param);
      State goal_state;
      goal_state.x = goal_config.x;
      goal_state.y = goal_config.y;
      goal_state.theta = goal_config.theta;
      goal_state.kappa = goal_config.kappa;
      goal_state.d = 0.0;

      vector<Control> cc_controls;
      cc_turn_controls(hc_cc_circle, goal_config, true, cc_controls);
      for (const auto &cc_control : cc_controls)
        EXPECT_LT(fabs(cc_control.sigma), sigma + EPS_SIGMA);
    }
  }
}

TEST(HC_CC_Circle, maxCurvature)
{
  State start_state;
  start_state.x = start_config.x;
  start_state.y = start_config.y;
  start_state.theta = start_config.theta;
  start_state.kappa = start_config.kappa;
  start_state.d = 0.0;

  for (const auto &sigma : sigmas)
  {
    Test_HC_CC_State_Space hc_cc_ss(KAPPA, sigma, DISCRETIZATION);
    HC_CC_Circle_Param hc_cc_circle_param = hc_cc_ss.get_hc_cc_circle_param();
    HC_CC_Circle hc_cc_circle = HC_CC_Circle(start_config, true, true, true, hc_cc_circle_param);

    for (const auto &delta : deltas)
    {
      // cc-turn
      Configuration goal_config = get_cc_goal_config(delta, hc_cc_circle_param);
      State goal_state;
      goal_state.x = goal_config.x;
      goal_state.y = goal_config.y;
      goal_state.theta = goal_config.theta;
      goal_state.kappa = goal_config.kappa;
      goal_state.d = 0.0;

      vector<Control> cc_controls;
      cc_turn_controls(hc_cc_circle, goal_config, true, cc_controls);
      vector<State> cc_path = hc_cc_ss.integrate(start_state, cc_controls);
      for (const auto &state : cc_path)
        EXPECT_LT(fabs(state.kappa), KAPPA + EPS_KAPPA);
    }
  }
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
