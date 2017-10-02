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

#ifndef HC_CC_STATE_SPACE_HPP
#define HC_CC_STATE_SPACE_HPP

#include <cmath>
#include <vector>

#include "steering_functions/hc_cc_state_space/hc_cc_circle.hpp"
#include "steering_functions/steering_functions.hpp"
#include "utilities.hpp"

using namespace std;
using namespace steer;

class HC_CC_State_Space
{
public:
  /** \brief Constructor */
  HC_CC_State_Space(double kappa, double sigma, double discretization);

  /** \brief Virtual function that returns controls of the shortest path from state1 to state2 */
  virtual vector<Control> get_controls(const State& state1, const State& state2) const = 0;

  /** \brief Returns path from state1 to state2 */
  vector<State> get_path(const State& state1, const State& state2) const;

  /** \brief Numeric integration using the forward euler method */
  vector<State> forward_euler(const State& state, const vector<Control>& controls) const;

protected:
  /** \brief Curvature, sharpness of clothoid */
  double kappa_, sigma_;

  /** \brief Discretization of path */
  double discretization_;

  /** \brief Parameters of a hc-/cc-circle */
  HC_CC_Circle_Param hc_cc_circle_param_;
};

#endif
