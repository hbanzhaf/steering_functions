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

#ifndef CC00_DUBINS_STATE_SPACE_HPP
#define CC00_DUBINS_STATE_SPACE_HPP

#include <memory>
#include <vector>

#include "steering_functions/hc_cc_state_space/hc_cc_circle.hpp"
#include "steering_functions/hc_cc_state_space/hc_cc_state_space.hpp"
#include "steering_functions/hc_cc_state_space/paths.hpp"
#include "steering_functions/steering_functions.hpp"

namespace steering
{

/** \brief
    An implementation of continuous curvature (CC) steer for a Dubins car
    with zero curvature at the start and goal configuration as described in:
    T. Fraichard and A. Scheuer, "From Reeds and Shepp's to continuous-curvature
    paths," IEEE Transactions on Robotics (Volume 20, Issue: 6, Dec. 2004).
    It evaluates all Dubins families and returns the shortest path.
    */
class CC00_Dubins_State_Space : public HC_CC_State_Space
{
public:
  /** \brief Constructor */
  CC00_Dubins_State_Space(double kappa, double sigma, double discretization = 0.1, bool forwards = true);

  /** \brief Destructor */
  ~CC00_Dubins_State_Space();

  /** \brief Returns a sequence of turns and straight lines connecting the two circles c1 and c2 */
  CC_Dubins_Path* cc00_circles_dubins_path(const HC_CC_Circle& c1, const HC_CC_Circle& c2) const;

  /** \brief Returns a sequence of turns and straight lines connecting a start and an end configuration */
  CC_Dubins_Path* cc00_dubins(const State& state1, const State& state2) const;

  /** \brief Returns shortest path length from state1 to state2 */
  double get_distance(const State& state1, const State& state2) const;

  /** \brief Returns controls of the shortest path from state1 to state2 */
  std::vector<Control> get_controls(const State& state1, const State& state2) const;

private:
  /** \brief Driving direction */
  bool forwards_;

  /** \brief Pimpl Idiom: class that contains functions to compute the families  */
  class CC00_Dubins;

  /** \brief Pimpl Idiom: unique pointer on class with families  */
  std::unique_ptr<CC00_Dubins> cc00_dubins_;
};

} // namespace steering

#endif
