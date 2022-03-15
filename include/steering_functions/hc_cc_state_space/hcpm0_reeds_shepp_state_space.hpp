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

#ifndef HCPM0_REEDS_SHEPP_STATE_SPACE_HPP
#define HCPM0_REEDS_SHEPP_STATE_SPACE_HPP

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "configuration.hpp"
#include "hc_cc_circle.hpp"
#include "hc_cc_state_space.hpp"
#include "paths.hpp"
#include "steering_functions/steering_functions.hpp"
#include "steering_functions/utilities/utilities.hpp"

namespace steering
{

/** \brief
    An implementation of hybrid curvature (HC) steer with either positive
    (p) or negative (n) max. curvature at the start configuration and zero
    curvature at the goal configuration, also see: H. Banzhaf et al.,
    "Hybrid Curvature Steer: A Novel Extend Function for Sampling-Based Non-
    holonomic Motion Planning in Tight Environments," IEEE International
    Conference on Intelligent Transportation Systems (Oct. 2017).
    It evaluates all Reeds-Shepp families plus the four families TTT, TcST,
    TScT, TcScT, where "T" stands for a turn, "S" for a straight line and
    "c" for a cusp, and returns the shortest path.
    */
class HCpm0_Reeds_Shepp_State_Space : public HC_CC_State_Space
{
public:
  /** \brief Constructor */
  HCpm0_Reeds_Shepp_State_Space(double kappa, double sigma, double discretization = 0.1);

  /** \brief Destructor */
  ~HCpm0_Reeds_Shepp_State_Space();

  /** \brief Returns a sequence of turns and straight lines connecting the two circles c1 and c2 */
  HC_CC_RS_Path* hcpm0_circles_rs_path(const HC_CC_Circle& c1, const HC_CC_Circle& c2) const;

  /** \brief Returns a sequence of turns and straight lines connecting a start and an end configuration */
  HC_CC_RS_Path* hcpm0_reeds_shepp(const State& state1, const State& state2) const;

  /** \brief Returns shortest path length from state1 to state2 */
  double get_distance(const State& state1, const State& state2) const;

  /** \brief Returns controls of the shortest path from state1 to state2 */
  std::vector<Control> get_controls(const State& state1, const State& state2) const;

private:
  /** \brief Pimpl Idiom: class that contains functions to compute the families  */
  class HCpm0_Reeds_Shepp;

  /** \brief Pimpl Idiom: unique pointer on class with families  */
  std::unique_ptr<HCpm0_Reeds_Shepp> hcpm0_reeds_shepp_;

  /** \brief Parameter of a rs-circle */
  HC_CC_Circle_Param rs_circle_param_;

  /** \brief Outer radius of a hc-/cc-circle */
  double radius_;

  /** \brief Angle between a configuration on the hc-/cc-circle and the tangent to the circle at that position */
  double mu_;
};

} // namespace steering

#endif
