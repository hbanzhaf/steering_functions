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
#include <time.h>
#include <fstream>
#include <iostream>

#include "steering_functions/dubins_state_space/dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/cc0pm_dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/cc_dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/cc_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/ccpm0_dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/ccpmpm_dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hc00_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hc0pm_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hc_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hcpm0_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hcpmpm_reeds_shepp_state_space.hpp"
#include "steering_functions/reeds_shepp_state_space/reeds_shepp_state_space.hpp"
#include "steering_functions/steering_functions.hpp"

using namespace std;
using namespace steer;

#define EPS_DISTANCE 0.01                // [m]
#define EPS_YAW 0.01                     // [rad]
#define EPS_KAPPA 1e-6                   // [1/m]
#define KAPPA 1.0                        // [1/m]
#define SIGMA 1.0                        // [1/m^2]
#define DISCRETIZATION 0.1               // [m]
#define SAMPLES 1e5                      // [-]
#define OPERATING_REGION_X 20.0          // [m]
#define OPERATING_REGION_Y 20.0          // [m]
#define OPERATING_REGION_THETA 2 * M_PI  // [rad]
#define random(lower, upper) (rand() * (upper - lower) / RAND_MAX + lower)
#define random_boolean() rand() % 2

State get_random_state()
{
  State state;
  state.x = random(-OPERATING_REGION_X / 2.0, OPERATING_REGION_X / 2.0);
  state.y = random(-OPERATING_REGION_Y / 2.0, OPERATING_REGION_Y / 2.0);
  state.theta = random(-OPERATING_REGION_THETA / 2.0, OPERATING_REGION_THETA / 2.0);
  state.kappa = random_boolean() * random(-KAPPA, KAPPA);
  state.d = 0.0;
  return state;
}

double get_distance(const State& state1, const State& state2)
{
  return sqrt(pow(state2.x - state1.x, 2) + pow(state2.y - state1.y, 2));
}

double get_path_length(const vector<State>& path)
{
  double path_length = 0;
  State state1 = path.front();
  for (const auto& state2 : path)
  {
    path_length += get_distance(state1, state2);
    state1 = state2;
  }
  return path_length;
}

double get_path_length(const vector<Control>& controls)
{
  double path_length = 0;
  for (const auto& control : controls)
    path_length += fabs(control.delta_s);
  return path_length;
}

bool is_equal(const State& state1, const State& state2)
{
  if (state1.x == state2.x && state1.y == state2.y && state1.theta == state2.theta && state1.kappa == state2.kappa &&
      state1.d == state2.d)
  {
    return true;
  }
  return false;
}

CC_Dubins_State_Space cc_dubins_forwards_ss(KAPPA, SIGMA, DISCRETIZATION, true);
CC_Dubins_State_Space cc_dubins_backwards_ss(KAPPA, SIGMA, DISCRETIZATION, false);
CC0pm_Dubins_State_Space cc0pm_dubins_forwards_ss(KAPPA, SIGMA, DISCRETIZATION, true);
CC0pm_Dubins_State_Space cc0pm_dubins_backwards_ss(KAPPA, SIGMA, DISCRETIZATION, false);
CCpm0_Dubins_State_Space ccpm0_dubins_forwards_ss(KAPPA, SIGMA, DISCRETIZATION, true);
CCpm0_Dubins_State_Space ccpm0_dubins_backwards_ss(KAPPA, SIGMA, DISCRETIZATION, false);
CCpmpm_Dubins_State_Space ccpmpm_dubins_forwards_ss(KAPPA, SIGMA, DISCRETIZATION, true);
CCpmpm_Dubins_State_Space ccpmpm_dubins_backwards_ss(KAPPA, SIGMA, DISCRETIZATION, false);
Dubins_State_Space dubins_forwards_ss(KAPPA, DISCRETIZATION, true);
Dubins_State_Space dubins_backwards_ss(KAPPA, DISCRETIZATION, false);
CC_Reeds_Shepp_State_Space cc_rs_ss(KAPPA, SIGMA, DISCRETIZATION);
HC_Reeds_Shepp_State_Space hc_ss(KAPPA, SIGMA, DISCRETIZATION);
HC00_Reeds_Shepp_State_Space hc00_ss(KAPPA, SIGMA, DISCRETIZATION);
HC0pm_Reeds_Shepp_State_Space hc0pm_ss(KAPPA, SIGMA, DISCRETIZATION);
HCpm0_Reeds_Shepp_State_Space hcpm0_ss(KAPPA, SIGMA, DISCRETIZATION);
HCpmpm_Reeds_Shepp_State_Space hcpmpm_ss(KAPPA, SIGMA, DISCRETIZATION);
Reeds_Shepp_State_Space rs_ss(KAPPA, DISCRETIZATION);
int seed(time(nullptr));

TEST(SteeringFunctions, pathLength)
{
  srand(seed);
  for (int i = 0; i < SAMPLES; i++)
  {
    State start = get_random_state();
    State goal;
    if (i == 0)
      goal = start;
    else if (i == 1)
    {
      goal.x = start.x + EPS_DISTANCE;
      goal.y = start.y + EPS_DISTANCE;
      goal.theta = start.theta + EPS_YAW;
      goal.kappa = start.kappa;
      goal.d = start.d;
    }
    else
      goal = get_random_state();

    vector<State> cc_dubins_path_forwards = cc_dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(cc_dubins_forwards_ss.get_distance(start, goal) - get_path_length(cc_dubins_path_forwards)),
              EPS_DISTANCE);

    vector<State> cc_dubins_path_backwards = cc_dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(cc_dubins_backwards_ss.get_distance(start, goal) - get_path_length(cc_dubins_path_backwards)),
              EPS_DISTANCE);

    vector<State> cc0pm_dubins_path_forwards = cc0pm_dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(cc0pm_dubins_forwards_ss.get_distance(start, goal) - get_path_length(cc0pm_dubins_path_forwards)),
              EPS_DISTANCE);

    vector<State> cc0pm_dubins_path_backwards = cc0pm_dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(cc0pm_dubins_backwards_ss.get_distance(start, goal) - get_path_length(cc0pm_dubins_path_backwards)),
              EPS_DISTANCE);

    vector<State> ccpm0_dubins_path_forwards = ccpm0_dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(ccpm0_dubins_forwards_ss.get_distance(start, goal) - get_path_length(ccpm0_dubins_path_forwards)),
              EPS_DISTANCE);

    vector<State> ccpm0_dubins_path_backwards = ccpm0_dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(ccpm0_dubins_backwards_ss.get_distance(start, goal) - get_path_length(ccpm0_dubins_path_backwards)),
              EPS_DISTANCE);

    vector<State> ccpmpm_dubins_path_forwards = ccpmpm_dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(ccpmpm_dubins_forwards_ss.get_distance(start, goal) - get_path_length(ccpmpm_dubins_path_forwards)),
              EPS_DISTANCE);

    vector<State> ccpmpm_dubins_path_backwards = ccpmpm_dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(
        fabs(ccpmpm_dubins_backwards_ss.get_distance(start, goal) - get_path_length(ccpmpm_dubins_path_backwards)),
        EPS_DISTANCE);

    vector<State> dubins_path_forwards = dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(dubins_forwards_ss.get_distance(start, goal) - get_path_length(dubins_path_forwards)), EPS_DISTANCE);

    vector<State> dubins_path_backwards = dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(dubins_backwards_ss.get_distance(start, goal) - get_path_length(dubins_path_backwards)),
              EPS_DISTANCE);

    vector<State> cc_rs_path = cc_rs_ss.get_path(start, goal);
    EXPECT_LT(fabs(cc_rs_ss.get_distance(start, goal) - get_path_length(cc_rs_path)), EPS_DISTANCE);

    vector<State> hc_path = hc_ss.get_path(start, goal);
    EXPECT_LT(fabs(hc_ss.get_distance(start, goal) - get_path_length(hc_path)), EPS_DISTANCE);

    vector<State> hc00_path = hc00_ss.get_path(start, goal);
    EXPECT_LT(fabs(hc00_ss.get_distance(start, goal) - get_path_length(hc00_path)), EPS_DISTANCE);

    vector<State> hc0pm_path = hc0pm_ss.get_path(start, goal);
    EXPECT_LT(fabs(hc0pm_ss.get_distance(start, goal) - get_path_length(hc0pm_path)), EPS_DISTANCE);

    vector<State> hcpm0_path = hcpm0_ss.get_path(start, goal);
    EXPECT_LT(fabs(hcpm0_ss.get_distance(start, goal) - get_path_length(hcpm0_path)), EPS_DISTANCE);

    vector<State> hcpmpm_path = hcpmpm_ss.get_path(start, goal);
    EXPECT_LT(fabs(hcpmpm_ss.get_distance(start, goal) - get_path_length(hcpmpm_path)), EPS_DISTANCE);

    vector<State> rs_path = rs_ss.get_path(start, goal);
    EXPECT_LT(fabs(rs_ss.get_distance(start, goal) - get_path_length(rs_path)), EPS_DISTANCE);
  }
}

TEST(SteeringFunctions, reachingGoal)
{
  srand(seed);
  for (int i = 0; i < SAMPLES; i++)
  {
    State start = get_random_state();
    State goal;
    if (i == 0)
      goal = start;
    else if (i == 1)
    {
      goal.x = start.x + EPS_DISTANCE;
      goal.y = start.y + EPS_DISTANCE;
      goal.theta = start.theta + EPS_YAW;
      goal.kappa = start.kappa;
      goal.d = start.d;
    }
    else
      goal = get_random_state();

    vector<State> cc_dubins_path_forwards = cc_dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, cc_dubins_path_forwards.back()), EPS_DISTANCE);

    vector<State> cc_dubins_path_backwards = cc_dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, cc_dubins_path_backwards.back()), EPS_DISTANCE);

    vector<State> cc0pm_dubins_path_forwards = cc0pm_dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, cc0pm_dubins_path_forwards.back()), EPS_DISTANCE);

    vector<State> cc0pm_dubins_path_backwards = cc0pm_dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, cc0pm_dubins_path_backwards.back()), EPS_DISTANCE);

    vector<State> ccpm0_dubins_path_forwards = ccpm0_dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, ccpm0_dubins_path_forwards.back()), EPS_DISTANCE);

    vector<State> ccpm0_dubins_path_backwards = ccpm0_dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, ccpm0_dubins_path_backwards.back()), EPS_DISTANCE);

    vector<State> ccpmpm_dubins_path_forwards = ccpmpm_dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, ccpmpm_dubins_path_forwards.back()), EPS_DISTANCE);

    vector<State> ccpmpm_dubins_path_backwards = ccpmpm_dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, ccpmpm_dubins_path_backwards.back()), EPS_DISTANCE);

    vector<State> dubins_path_forwards = dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, dubins_path_forwards.back()), EPS_DISTANCE);

    vector<State> dubins_path_backwards = dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, dubins_path_backwards.back()), EPS_DISTANCE);

    vector<State> cc_rs_path = cc_rs_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, cc_rs_path.back()), EPS_DISTANCE);

    vector<State> hc_path = hc_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, hc_path.back()), EPS_DISTANCE);

    vector<State> hc00_path = hc00_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, hc00_path.back()), EPS_DISTANCE);

    vector<State> hc0pm_path = hc0pm_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, hc0pm_path.back()), EPS_DISTANCE);

    vector<State> hcpm0_path = hcpm0_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, hcpm0_path.back()), EPS_DISTANCE);

    vector<State> hcpmpm_path = hcpmpm_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, hcpmpm_path.back()), EPS_DISTANCE);

    vector<State> rs_path = rs_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, rs_path.back()), EPS_DISTANCE);
  }
}

TEST(SteeringFunctions, curvatureContinuity)
{
  srand(seed);
  for (int i = 0; i < SAMPLES; i++)
  {
    State start = get_random_state();
    State goal;
    if (i == 0)
      goal = start;
    else if (i == 1)
    {
      goal.x = start.x + EPS_DISTANCE;
      goal.y = start.y + EPS_DISTANCE;
      goal.theta = start.theta + EPS_YAW;
      goal.kappa = start.kappa;
      goal.d = start.d;
    }
    else
      goal = get_random_state();

    State state1;

    vector<State> cc_dubins_forwards_path = cc_dubins_forwards_ss.get_path(start, goal);
    state1 = cc_dubins_forwards_path.front();
    for (const auto& state2 : cc_dubins_forwards_path)
    {
      EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      state1 = state2;
    }

    vector<State> cc_dubins_backwards_path = cc_dubins_backwards_ss.get_path(start, goal);
    state1 = cc_dubins_backwards_path.front();
    for (const auto& state2 : cc_dubins_backwards_path)
    {
      EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      state1 = state2;
    }

    vector<State> cc0pm_dubins_forwards_path = cc0pm_dubins_forwards_ss.get_path(start, goal);
    state1 = cc0pm_dubins_forwards_path.front();
    for (const auto& state2 : cc0pm_dubins_forwards_path)
    {
      EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      state1 = state2;
    }

    vector<State> cc0pm_dubins_backwards_path = cc0pm_dubins_backwards_ss.get_path(start, goal);
    state1 = cc0pm_dubins_backwards_path.front();
    for (const auto& state2 : cc0pm_dubins_backwards_path)
    {
      EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      state1 = state2;
    }

    vector<State> ccpm0_dubins_forwards_path = ccpm0_dubins_forwards_ss.get_path(start, goal);
    state1 = ccpm0_dubins_forwards_path.front();
    for (const auto& state2 : ccpm0_dubins_forwards_path)
    {
      EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      state1 = state2;
    }

    vector<State> ccpm0_dubins_backwards_path = ccpm0_dubins_backwards_ss.get_path(start, goal);
    state1 = ccpm0_dubins_backwards_path.front();
    for (const auto& state2 : ccpm0_dubins_backwards_path)
    {
      EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      state1 = state2;
    }

    vector<State> ccpmpm_dubins_forwards_path = ccpmpm_dubins_forwards_ss.get_path(start, goal);
    state1 = ccpmpm_dubins_forwards_path.front();
    for (const auto& state2 : ccpmpm_dubins_forwards_path)
    {
      EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      state1 = state2;
    }

    vector<State> ccpmpm_dubins_backwards_path = ccpmpm_dubins_backwards_ss.get_path(start, goal);
    state1 = ccpmpm_dubins_backwards_path.front();
    for (const auto& state2 : ccpmpm_dubins_backwards_path)
    {
      EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      state1 = state2;
    }

    vector<State> cc_rs_path = cc_rs_ss.get_path(start, goal);
    state1 = cc_rs_path.front();
    for (const auto& state2 : cc_rs_path)
    {
      EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      state1 = state2;
    }

    vector<State> hc_path = hc_ss.get_path(start, goal);
    state1 = hc_path.front();
    for (const auto& state2 : hc_path)
    {
      if (state1.d * state2.d >= 0)
      {
        EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      }
      state1 = state2;
    }

    vector<State> hc00_path = hc00_ss.get_path(start, goal);
    state1 = hc00_path.front();
    for (const auto& state2 : hc00_path)
    {
      if (state1.d * state2.d >= 0)
      {
        EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      }
      state1 = state2;
    }

    vector<State> hc0pm_path = hc0pm_ss.get_path(start, goal);
    state1 = hc0pm_path.front();
    for (const auto& state2 : hc0pm_path)
    {
      if (state1.d * state2.d >= 0)
      {
        EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      }
      state1 = state2;
    }

    vector<State> hcpm0_path = hcpm0_ss.get_path(start, goal);
    state1 = hcpm0_path.front();
    for (const auto& state2 : hcpm0_path)
    {
      if (state1.d * state2.d >= 0)
      {
        EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      }
      state1 = state2;
    }

    vector<State> hcpmpm_path = hcpmpm_ss.get_path(start, goal);
    state1 = hcpmpm_path.front();
    for (const auto& state2 : hcpmpm_path)
    {
      if (state1.d * state2.d >= 0)
      {
        EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
      }
      state1 = state2;
    }
  }
}

TEST(SteeringFunctions, interpolation)
{
  srand(seed);
  double s(0.0);
  for (int i = 0; i < SAMPLES; i++)
  {
    double t = ((double)rand() / (RAND_MAX));
    State start = get_random_state();
    State goal;
    if (i == 0)
      goal = start;
    else if (i == 1)
    {
      goal.x = start.x + EPS_DISTANCE;
      goal.y = start.y + EPS_DISTANCE;
      goal.theta = start.theta + EPS_YAW;
      goal.kappa = start.kappa;
      goal.d = start.d;
    }
    else
      goal = get_random_state();

    vector<Control> cc_dubins_forwards_controls = cc_dubins_forwards_ss.get_controls(start, goal);
    double cc_dubins_forwards_s_path = get_path_length(cc_dubins_forwards_controls);
    double cc_dubins_forwards_s_inter = t * cc_dubins_forwards_s_path;
    s = 0.0;
    vector<Control> cc_dubins_forwards_controls_inter;
    cc_dubins_forwards_controls_inter.reserve(cc_dubins_forwards_controls.size());
    for (const auto& control : cc_dubins_forwards_controls)
    {
      double abs_delta_s = fabs(control.delta_s);
      s += abs_delta_s;
      if (s < cc_dubins_forwards_s_inter)
        cc_dubins_forwards_controls_inter.push_back(control);
      else
      {
        Control control_inter;
        control_inter.delta_s = sgn(control.delta_s) * (abs_delta_s - (s - cc_dubins_forwards_s_inter));
        control_inter.kappa = control.kappa;
        control_inter.sigma = control.sigma;
        cc_dubins_forwards_controls_inter.push_back(control_inter);
        break;
      }
    }
    vector<State> cc_dubins_forwards_path = cc_dubins_forwards_ss.integrate(start, cc_dubins_forwards_controls_inter);
    State cc_dubins_forwards_state_inter = cc_dubins_forwards_ss.interpolate(start, cc_dubins_forwards_controls, t);
    EXPECT_EQ(is_equal(cc_dubins_forwards_path.back(), cc_dubins_forwards_state_inter), true);

    vector<Control> cc0pm_dubins_forwards_controls = cc0pm_dubins_forwards_ss.get_controls(start, goal);
    double cc0pm_dubins_forwards_s_path = get_path_length(cc0pm_dubins_forwards_controls);
    double cc0pm_dubins_forwards_s_inter = t * cc0pm_dubins_forwards_s_path;
    s = 0.0;
    vector<Control> cc0pm_dubins_forwards_controls_inter;
    cc0pm_dubins_forwards_controls_inter.reserve(cc0pm_dubins_forwards_controls.size());
    for (const auto& control : cc0pm_dubins_forwards_controls)
    {
      double abs_delta_s = fabs(control.delta_s);
      s += abs_delta_s;
      if (s < cc0pm_dubins_forwards_s_inter)
        cc0pm_dubins_forwards_controls_inter.push_back(control);
      else
      {
        Control control_inter;
        control_inter.delta_s = sgn(control.delta_s) * (abs_delta_s - (s - cc0pm_dubins_forwards_s_inter));
        control_inter.kappa = control.kappa;
        control_inter.sigma = control.sigma;
        cc0pm_dubins_forwards_controls_inter.push_back(control_inter);
        break;
      }
    }
    vector<State> cc0pm_dubins_forwards_path =
        cc0pm_dubins_forwards_ss.integrate(start, cc0pm_dubins_forwards_controls_inter);
    State cc0pm_dubins_forwards_state_inter =
        cc0pm_dubins_forwards_ss.interpolate(start, cc0pm_dubins_forwards_controls, t);
    EXPECT_EQ(is_equal(cc0pm_dubins_forwards_path.back(), cc0pm_dubins_forwards_state_inter), true);

    vector<Control> ccpm0_dubins_forwards_controls = ccpm0_dubins_forwards_ss.get_controls(start, goal);
    double ccpm0_dubins_forwards_s_path = get_path_length(ccpm0_dubins_forwards_controls);
    double ccpm0_dubins_forwards_s_inter = t * ccpm0_dubins_forwards_s_path;
    s = 0.0;
    vector<Control> ccpm0_dubins_forwards_controls_inter;
    ccpm0_dubins_forwards_controls_inter.reserve(ccpm0_dubins_forwards_controls.size());
    for (const auto& control : ccpm0_dubins_forwards_controls)
    {
      double abs_delta_s = fabs(control.delta_s);
      s += abs_delta_s;
      if (s < ccpm0_dubins_forwards_s_inter)
        ccpm0_dubins_forwards_controls_inter.push_back(control);
      else
      {
        Control control_inter;
        control_inter.delta_s = sgn(control.delta_s) * (abs_delta_s - (s - ccpm0_dubins_forwards_s_inter));
        control_inter.kappa = control.kappa;
        control_inter.sigma = control.sigma;
        ccpm0_dubins_forwards_controls_inter.push_back(control_inter);
        break;
      }
    }
    vector<State> ccpm0_dubins_forwards_path =
        ccpm0_dubins_forwards_ss.integrate(start, ccpm0_dubins_forwards_controls_inter);
    State ccpm0_dubins_forwards_state_inter =
        ccpm0_dubins_forwards_ss.interpolate(start, ccpm0_dubins_forwards_controls, t);
    EXPECT_EQ(is_equal(ccpm0_dubins_forwards_path.back(), ccpm0_dubins_forwards_state_inter), true);

    vector<Control> ccpmpm_dubins_forwards_controls = ccpmpm_dubins_forwards_ss.get_controls(start, goal);
    double ccpmpm_dubins_forwards_s_path = get_path_length(ccpmpm_dubins_forwards_controls);
    double ccpmpm_dubins_forwards_s_inter = t * ccpmpm_dubins_forwards_s_path;
    s = 0.0;
    vector<Control> ccpmpm_dubins_forwards_controls_inter;
    ccpmpm_dubins_forwards_controls_inter.reserve(ccpmpm_dubins_forwards_controls.size());
    for (const auto& control : ccpmpm_dubins_forwards_controls)
    {
      double abs_delta_s = fabs(control.delta_s);
      s += abs_delta_s;
      if (s < ccpmpm_dubins_forwards_s_inter)
        ccpmpm_dubins_forwards_controls_inter.push_back(control);
      else
      {
        Control control_inter;
        control_inter.delta_s = sgn(control.delta_s) * (abs_delta_s - (s - ccpmpm_dubins_forwards_s_inter));
        control_inter.kappa = control.kappa;
        control_inter.sigma = control.sigma;
        ccpmpm_dubins_forwards_controls_inter.push_back(control_inter);
        break;
      }
    }
    vector<State> ccpmpm_dubins_forwards_path =
        ccpmpm_dubins_forwards_ss.integrate(start, ccpmpm_dubins_forwards_controls_inter);
    State ccpmpm_dubins_forwards_state_inter =
        ccpmpm_dubins_forwards_ss.interpolate(start, ccpmpm_dubins_forwards_controls, t);
    EXPECT_EQ(is_equal(ccpmpm_dubins_forwards_path.back(), ccpmpm_dubins_forwards_state_inter), true);

    vector<Control> dubins_forwards_controls = dubins_forwards_ss.get_controls(start, goal);
    double dubins_forwards_s_path = get_path_length(dubins_forwards_controls);
    double dubins_forwards_s_inter = t * dubins_forwards_s_path;
    s = 0.0;
    vector<Control> dubins_forwards_controls_inter;
    dubins_forwards_controls_inter.reserve(dubins_forwards_controls.size());
    for (const auto& control : dubins_forwards_controls)
    {
      double abs_delta_s = fabs(control.delta_s);
      s += abs_delta_s;
      if (s < dubins_forwards_s_inter)
        dubins_forwards_controls_inter.push_back(control);
      else
      {
        Control control_inter;
        control_inter.delta_s = sgn(control.delta_s) * (abs_delta_s - (s - dubins_forwards_s_inter));
        control_inter.kappa = control.kappa;
        control_inter.sigma = control.sigma;
        dubins_forwards_controls_inter.push_back(control_inter);
        break;
      }
    }
    vector<State> dubins_forwards_path = dubins_forwards_ss.integrate(start, dubins_forwards_controls_inter);
    State dubins_forwards_state_inter = dubins_forwards_ss.interpolate(start, dubins_forwards_controls, t);
    EXPECT_EQ(is_equal(dubins_forwards_path.back(), dubins_forwards_state_inter), true);

    vector<Control> cc_rs_controls = cc_rs_ss.get_controls(start, goal);
    double cc_rs_s_path = get_path_length(cc_rs_controls);
    double cc_rs_s_inter = t * cc_rs_s_path;
    s = 0.0;
    vector<Control> cc_rs_controls_inter;
    cc_rs_controls_inter.reserve(cc_rs_controls.size());
    for (const auto& control : cc_rs_controls)
    {
      double abs_delta_s = fabs(control.delta_s);
      s += abs_delta_s;
      if (s < cc_rs_s_inter)
        cc_rs_controls_inter.push_back(control);
      else
      {
        Control control_inter;
        control_inter.delta_s = sgn(control.delta_s) * (abs_delta_s - (s - cc_rs_s_inter));
        control_inter.kappa = control.kappa;
        control_inter.sigma = control.sigma;
        cc_rs_controls_inter.push_back(control_inter);
        break;
      }
    }
    vector<State> cc_rs_path = cc_rs_ss.integrate(start, cc_rs_controls_inter);
    State cc_rs_state_inter = cc_rs_ss.interpolate(start, cc_rs_controls, t);
    EXPECT_EQ(is_equal(cc_rs_path.back(), cc_rs_state_inter), true);

    vector<Control> hc_controls = hc_ss.get_controls(start, goal);
    double hc_s_path = get_path_length(hc_controls);
    double hc_s_inter = t * hc_s_path;
    s = 0.0;
    vector<Control> hc_controls_inter;
    hc_controls_inter.reserve(hc_controls.size());
    for (const auto& control : hc_controls)
    {
      double abs_delta_s = fabs(control.delta_s);
      s += abs_delta_s;
      if (s < hc_s_inter)
        hc_controls_inter.push_back(control);
      else
      {
        Control control_inter;
        control_inter.delta_s = sgn(control.delta_s) * (abs_delta_s - (s - hc_s_inter));
        control_inter.kappa = control.kappa;
        control_inter.sigma = control.sigma;
        hc_controls_inter.push_back(control_inter);
        break;
      }
    }
    vector<State> hc_path = hc_ss.integrate(start, hc_controls_inter);
    State hc_state_inter = hc_ss.interpolate(start, hc_controls, t);
    EXPECT_EQ(is_equal(hc_path.back(), hc_state_inter), true);

    vector<Control> hc00_controls = hc00_ss.get_controls(start, goal);
    double hc00_s_path = get_path_length(hc00_controls);
    double hc00_s_inter = t * hc00_s_path;
    s = 0.0;
    vector<Control> hc00_controls_inter;
    hc00_controls_inter.reserve(hc00_controls.size());
    for (const auto& control : hc00_controls)
    {
      double abs_delta_s = fabs(control.delta_s);
      s += abs_delta_s;
      if (s < hc00_s_inter)
        hc00_controls_inter.push_back(control);
      else
      {
        Control control_inter;
        control_inter.delta_s = sgn(control.delta_s) * (abs_delta_s - (s - hc00_s_inter));
        control_inter.kappa = control.kappa;
        control_inter.sigma = control.sigma;
        hc00_controls_inter.push_back(control_inter);
        break;
      }
    }
    vector<State> hc00_path = hc00_ss.integrate(start, hc00_controls_inter);
    State hc00_state_inter = hc00_ss.interpolate(start, hc00_controls, t);
    EXPECT_EQ(is_equal(hc00_path.back(), hc00_state_inter), true);

    vector<Control> hcpmpm_controls = hcpmpm_ss.get_controls(start, goal);
    double hcpmpm_s_path = get_path_length(hcpmpm_controls);
    double hcpmpm_s_inter = t * hcpmpm_s_path;
    s = 0.0;
    vector<Control> hcpmpm_controls_inter;
    hcpmpm_controls_inter.reserve(hcpmpm_controls.size());
    for (const auto& control : hcpmpm_controls)
    {
      double abs_delta_s = fabs(control.delta_s);
      s += abs_delta_s;
      if (s < hcpmpm_s_inter)
        hcpmpm_controls_inter.push_back(control);
      else
      {
        Control control_inter;
        control_inter.delta_s = sgn(control.delta_s) * (abs_delta_s - (s - hcpmpm_s_inter));
        control_inter.kappa = control.kappa;
        control_inter.sigma = control.sigma;
        hcpmpm_controls_inter.push_back(control_inter);
        break;
      }
    }
    vector<State> hcpmpm_path = hcpmpm_ss.integrate(start, hcpmpm_controls_inter);
    State hcpmpm_state_inter = hcpmpm_ss.interpolate(start, hcpmpm_controls, t);
    EXPECT_EQ(is_equal(hcpmpm_path.back(), hcpmpm_state_inter), true);

    vector<Control> rs_controls = rs_ss.get_controls(start, goal);
    double rs_s_path = get_path_length(rs_controls);
    double rs_s_inter = t * rs_s_path;
    s = 0.0;
    vector<Control> rs_controls_inter;
    rs_controls_inter.reserve(rs_controls.size());
    for (const auto& control : rs_controls)
    {
      double abs_delta_s = fabs(control.delta_s);
      s += abs_delta_s;
      if (s < rs_s_inter)
        rs_controls_inter.push_back(control);
      else
      {
        Control control_inter;
        control_inter.delta_s = sgn(control.delta_s) * (abs_delta_s - (s - rs_s_inter));
        control_inter.kappa = control.kappa;
        control_inter.sigma = control.sigma;
        rs_controls_inter.push_back(control_inter);
        break;
      }
    }
    vector<State> rs_path = rs_ss.integrate(start, rs_controls_inter);
    State rs_state_inter = rs_ss.interpolate(start, rs_controls, t);
    EXPECT_EQ(is_equal(rs_path.back(), rs_state_inter), true);
  }
}

TEST(SteeringFunctions, symmetry)
{
  srand(seed);
  for (int i = 0; i < SAMPLES; ++i)
  {
    State start = get_random_state();
    State goal = get_random_state();

    double cc_rs_distance_forwards = cc_rs_ss.get_distance(start, goal);
    double cc_rs_distance_backwards = cc_rs_ss.get_distance(goal, start);
    EXPECT_LT(fabs(cc_rs_distance_forwards - cc_rs_distance_backwards), EPS_DISTANCE);

    double hc00_distance_forwards = hc00_ss.get_distance(start, goal);
    double hc00_distance_backwards = hc00_ss.get_distance(goal, start);
    EXPECT_LT(fabs(hc00_distance_forwards - hc00_distance_backwards), EPS_DISTANCE);

    double hcpmpm_distance_forwards = hcpmpm_ss.get_distance(start, goal);
    double hcpmpm_distance_backwards = hcpmpm_ss.get_distance(goal, start);
    EXPECT_LT(fabs(hcpmpm_distance_forwards - hcpmpm_distance_backwards), EPS_DISTANCE);

    double rs_distance_forwards = rs_ss.get_distance(start, goal);
    double rs_distance_backwards = rs_ss.get_distance(goal, start);
    EXPECT_LT(fabs(rs_distance_forwards - rs_distance_backwards), EPS_DISTANCE);
  }
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
