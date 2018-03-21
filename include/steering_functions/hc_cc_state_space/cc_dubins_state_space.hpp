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

#ifndef CC_DUBINS_STATE_SPACE_HPP
#define CC_DUBINS_STATE_SPACE_HPP

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "configuration.hpp"
#include "hc_cc_circle.hpp"
#include "hc_cc_state_space.hpp"
#include "paths.hpp"
#include "steering_functions/hc_cc_state_space/cc00_dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/cc0pm_dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/ccpm0_dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/ccpmpm_dubins_state_space.hpp"
#include "steering_functions/steering_functions.hpp"
#include "steering_functions/utilities/utilities.hpp"

using namespace std;
using namespace steer;

/** \brief
    An implementation of continuous curvature (CC) steer for a Dubins car
    with arbitrary curvature at the start and goal configuration.
    */
class CC_Dubins_State_Space : public HC_CC_State_Space
{
public:
  /** \brief Constructor */
  CC_Dubins_State_Space(double kappa, double sigma, double discretization = 0.1, bool forwards = true);

  /** \brief Predicts a state forwards or backwards to zero and max. curvature */
  vector<pair<State, Control>> predict_state(const State& state, bool forwards) const;

  /** \brief Returns shortest path length from state1 to state2 */
  double get_distance(const State& state1, const State& state2) const;

  /** \brief Returns controls of the shortest path from state1 to state2 */
  vector<Control> get_controls(const State& state1, const State& state2) const;

private:
  /** \brief Driving direction */
  bool forwards_;

  /** \brief Required state spaces */
  CC00_Dubins_State_Space cc00_dubins_state_space_;
  CC0pm_Dubins_State_Space cc0pm_dubins_state_space_;
  CCpm0_Dubins_State_Space ccpm0_dubins_state_space_;
  CCpmpm_Dubins_State_Space ccpmpm_dubins_state_space_;
};

#endif
