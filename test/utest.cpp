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
#include "steering_functions/hc_cc_state_space/cc_dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/cc_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hc00_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hc0pm_reeds_shepp_state_space.hpp"
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
#define DISCRETIZATION 0.01              // [m]
#define SAMPLES 1e5                      // [-]
#define OPERATING_REGION_X 20.0          // [m]
#define OPERATING_REGION_Y 20.0          // [m]
#define OPERATING_REGION_THETA 2 * M_PI  // [rad]
#define random(lower, upper) (rand() * (upper - lower) / RAND_MAX + lower)

State get_random_state()
{
  State state;
  state.x = random(-OPERATING_REGION_X / 2.0, OPERATING_REGION_X / 2.0);
  state.y = random(-OPERATING_REGION_Y / 2.0, OPERATING_REGION_Y / 2.0);
  state.theta = random(-OPERATING_REGION_THETA / 2.0, OPERATING_REGION_THETA / 2.0);
  state.kappa = 0.0;
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

double get_mean(const vector<double>& v)
{
  double sum = accumulate(v.begin(), v.end(), 0.0);
  return sum / v.size();
}

double get_std(const vector<double>& v)
{
  double mean = get_mean(v);
  double diff_sq = 0;
  for (const auto& x : v)
  {
    diff_sq += pow(x - mean, 2);
  }
  return sqrt(diff_sq / v.size());
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
Dubins_State_Space dubins_forwards_ss(KAPPA, DISCRETIZATION, true);
Dubins_State_Space dubins_backwards_ss(KAPPA, DISCRETIZATION, false);
CC_Reeds_Shepp_State_Space cc_rs_ss(KAPPA, SIGMA, DISCRETIZATION);
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

    vector<State> dubins_path_forwards = dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(dubins_forwards_ss.get_distance(start, goal) - get_path_length(dubins_path_forwards)), EPS_DISTANCE);

    vector<State> dubins_path_backwards = dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(fabs(dubins_backwards_ss.get_distance(start, goal) - get_path_length(dubins_path_backwards)),
              EPS_DISTANCE);

    vector<State> cc_rs_path = cc_rs_ss.get_path(start, goal);
    EXPECT_LT(fabs(cc_rs_ss.get_distance(start, goal) - get_path_length(cc_rs_path)), EPS_DISTANCE);

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

    vector<State> dubins_path_forwards = dubins_forwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, dubins_path_forwards.back()), EPS_DISTANCE);

    vector<State> dubins_path_backwards = dubins_backwards_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, dubins_path_backwards.back()), EPS_DISTANCE);

    vector<State> cc_rs_path = cc_rs_ss.get_path(start, goal);
    EXPECT_LT(get_distance(goal, cc_rs_path.back()), EPS_DISTANCE);

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

    vector<State> cc_rs_path = cc_rs_ss.get_path(start, goal);
    state1 = cc_rs_path.front();
    for (const auto& state2 : cc_rs_path)
    {
      EXPECT_LE(fabs(state1.kappa - state2.kappa), DISCRETIZATION * SIGMA + EPS_KAPPA);
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

#include <ros/package.h>

struct Statistic
{
  State start;
  State goal;
  double computation_time;
  double path_length;
};

void write_to_file(const string& id, const vector<Statistic>& stats)
{
  string path_to_output = ros::package::getPath("steering_functions") + "/test/" + id + "_stats.csv";
  remove(path_to_output.c_str());
  string header = "start,goal,computation_time,path_length";
  fstream f;
  f.open(path_to_output, ios::app);
  f << header << endl;
  for (const auto& stat : stats)
  {
    State start = stat.start;
    State goal = stat.goal;
    f << start.x << " " << start.y << " " << start.theta << " " << start.kappa << " " << start.d << ",";
    f << goal.x << " " << goal.y << " " << goal.theta << " " << goal.kappa << " " << goal.d << ",";
    f << stat.computation_time << "," << stat.path_length << endl;
  }
}

vector<Statistic> get_stats(const string& id, const vector<State>& starts, const vector<State>& goals)
{
  clock_t clock_start;
  clock_t clock_finish;
  Statistic stat;
  vector<Statistic> stats;
  stats.reserve(SAMPLES);
  double path_length;
  for (auto start = starts.begin(), goal = goals.begin(); start != starts.end(); ++start, ++goal)
  {
    if (id == "CC_Dubins")
    {
      clock_start = clock();
      cc_dubins_forwards_ss.get_controls(*start, *goal);
      clock_finish = clock();
      path_length = cc_dubins_forwards_ss.get_distance(*start, *goal);
    }
    else if (id == "Dubins")
    {
      clock_start = clock();
      dubins_forwards_ss.get_controls(*start, *goal);
      clock_finish = clock();
      path_length = dubins_forwards_ss.get_distance(*start, *goal);
    }
    else if (id == "CC_RS")
    {
      clock_start = clock();
      cc_rs_ss.get_controls(*start, *goal);
      clock_finish = clock();
      path_length = cc_rs_ss.get_distance(*start, *goal);
    }
    else if (id == "HC00")
    {
      clock_start = clock();
      hc00_ss.get_controls(*start, *goal);
      clock_finish = clock();
      path_length = hc00_ss.get_distance(*start, *goal);
    }
    else if (id == "HC0pm")
    {
      clock_start = clock();
      hc0pm_ss.get_controls(*start, *goal);
      clock_finish = clock();
      path_length = hc0pm_ss.get_distance(*start, *goal);
    }
    else if (id == "HCpm0")
    {
      clock_start = clock();
      hcpm0_ss.get_controls(*start, *goal);
      clock_finish = clock();
      path_length = hcpm0_ss.get_distance(*start, *goal);
    }
    else if (id == "HCpmpm")
    {
      clock_start = clock();
      hcpmpm_ss.get_controls(*start, *goal);
      clock_finish = clock();
      path_length = hcpmpm_ss.get_distance(*start, *goal);
    }
    else if (id == "RS")
    {
      clock_start = clock();
      rs_ss.get_controls(*start, *goal);
      clock_finish = clock();
      path_length = rs_ss.get_distance(*start, *goal);
    }
    stat.start = *start;
    stat.goal = *goal;
    stat.computation_time = double(clock_finish - clock_start) / CLOCKS_PER_SEC;
    stat.path_length = path_length;
    stats.push_back(stat);
  }
  return stats;
}

TEST(SteeringFunctions, stats)
{
  srand(0);
  vector<double> computation_times;
  computation_times.reserve(SAMPLES);
  vector<State> starts;
  starts.reserve(SAMPLES);
  vector<State> goals;
  goals.reserve(SAMPLES);
  for (int i = 0; i < SAMPLES; i++)
  {
    State start = get_random_state();
    State goal = get_random_state();
    starts.push_back(start);
    goals.push_back(goal);
  }

  string cc_dubins_id = "CC_Dubins";
  vector<Statistic> cc_dubins_stats = get_stats(cc_dubins_id, starts, goals);
  computation_times.clear();
  for (const auto& stat : cc_dubins_stats)
  {
    computation_times.push_back(stat.computation_time);
  }
  cout << "[----------] " + cc_dubins_id + " mean [s] +/- std [s]: " << get_mean(computation_times) << " +/- "
       << get_std(computation_times) << endl;

  string dubins_id = "Dubins";
  vector<Statistic> dubins_stats = get_stats(dubins_id, starts, goals);
  computation_times.clear();
  for (const auto& stat : dubins_stats)
  {
    computation_times.push_back(stat.computation_time);
  }
  cout << "[----------] " + dubins_id + " mean [s] +/- std [s]: " << get_mean(computation_times) << " +/- "
       << get_std(computation_times) << endl;

  string cc_rs_id = "CC_RS";
  vector<Statistic> cc_rs_stats = get_stats(cc_rs_id, starts, goals);
  computation_times.clear();
  for (const auto& stat : cc_rs_stats)
  {
    computation_times.push_back(stat.computation_time);
  }
  cout << "[----------] " + cc_rs_id + " mean [s] +/- std [s]: " << get_mean(computation_times) << " +/- "
       << get_std(computation_times) << endl;

  string hc00_id = "HC00";
  vector<Statistic> hc00_stats = get_stats(hc00_id, starts, goals);
  computation_times.clear();
  for (const auto& stat : hc00_stats)
  {
    computation_times.push_back(stat.computation_time);
  }
  cout << "[----------] " + hc00_id + " mean [s] +/- std [s]: " << get_mean(computation_times) << " +/- "
       << get_std(computation_times) << endl;

  string hc0pm_id = "HC0pm";
  vector<Statistic> hc0pm_stats = get_stats(hc0pm_id, starts, goals);
  computation_times.clear();
  for (const auto& stat : hc0pm_stats)
  {
    computation_times.push_back(stat.computation_time);
  }
  cout << "[----------] " + hc0pm_id + " mean [s] +/- std [s]: " << get_mean(computation_times) << " +/- "
       << get_std(computation_times) << endl;

  string hcpm0_id = "HCpm0";
  vector<Statistic> hcpm0_stats = get_stats(hcpm0_id, starts, goals);
  computation_times.clear();
  for (const auto& stat : hcpm0_stats)
  {
    computation_times.push_back(stat.computation_time);
  }
  cout << "[----------] " + hcpm0_id + " mean [s] +/- std [s]: " << get_mean(computation_times) << " +/- "
       << get_std(computation_times) << endl;

  string hcpmpm_id = "HCpmpm";
  vector<Statistic> hcpmpm_stats = get_stats(hcpmpm_id, starts, goals);
  computation_times.clear();
  for (const auto& stat : hcpmpm_stats)
  {
    computation_times.push_back(stat.computation_time);
  }
  cout << "[----------] " + hcpmpm_id + " mean [s] +/- std [s]: " << get_mean(computation_times) << " +/- "
       << get_std(computation_times) << endl;

  string rs_id = "RS";
  vector<Statistic> rs_stats = get_stats(rs_id, starts, goals);
  computation_times.clear();
  for (const auto& stat : rs_stats)
  {
    computation_times.push_back(stat.computation_time);
  }
  cout << "[----------] " + rs_id + " mean [s] +/- std [s]: " << get_mean(computation_times) << " +/- "
       << get_std(computation_times) << endl;

  //  write_to_file(cc_dubins_id, cc_dubins_stats);
  //  write_to_file(dubins_id, dubins_stats);
  //  write_to_file(cc_rs_id, cc_rs_stats);
  //  write_to_file(hc00_id, hc00_stats);
  //  write_to_file(hc0pm_id, hc0pm_stats);
  //  write_to_file(hcpm0_id, hcpm0_stats);
  //  write_to_file(hcpmpm_id, hcpmpm_stats);
  //  write_to_file(rs_id, rs_stats);
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
