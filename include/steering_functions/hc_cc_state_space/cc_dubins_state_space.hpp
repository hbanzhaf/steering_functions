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

*  This source code is derived from Continuous Curvature (CC) Steer.
*  Copyright (c) 2016, Thierry Fraichard and Institut national de
*  recherche en informatique et en automatique (Inria), licensed under
*  the BSD license, cf. 3rd-party-licenses.txt file in the root
*  directory of this source tree.
**********************************************************************/

#ifndef CC_DUBINS_STATE_SPACE_HPP
#define CC_DUBINS_STATE_SPACE_HPP

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include "configuration.hpp"
#include "hc_cc_circle.hpp"
#include "hc_cc_state_space.hpp"
#include "paths.hpp"
#include "steering_functions/steering_functions.hpp"
#include "steering_functions/utilities/utilities.hpp"

using namespace std;
using namespace steer;

/** \brief
    An implementation of a Dubins car with continuous curvature (CC) steer
    as described in: T. Fraichard and A. Scheuer, "From Reeds and Shepp's
    to continuous-curvature paths," IEEE Transactions on Robotics (Volume
    20, Issue: 6, Dec. 2004).
    It evaluates all Dubins families and returns the shortest path.
    */
class CC_Dubins_State_Space : public HC_CC_State_Space
{
public:
  /** \brief Constructor */
  CC_Dubins_State_Space(double kappa, double sigma, double discretization = 0.1, bool forwards = true)
    : HC_CC_State_Space(kappa, sigma, discretization), forwards_(forwards)
  {
  }

  /** \brief Returns a sequence of turns and straight lines connecting a start and an end configuration */
  CC_Dubins_Path* cc_dubins(const State& state1, const State& state2) const;

  /** \brief Returns shortest path length from state1 to state2 */
  double get_distance(const State& state1, const State& state2) const;

  /** \brief Returns controls of the shortest path from state1 to state2 */
  vector<Control> get_controls(const State& state1, const State& state2) const;

private:
  /** \brief Driving direction */
  bool forwards_;
};

#endif
