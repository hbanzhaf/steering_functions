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

*  This source code is derived from the Open Motion Planning Library
*  (OMPL) V 1.3.1 (https://github.com/ompl/ompl).
*  Copyright (c) 2010, Rice University, licensed under the BSD license,
*  cf. 3rd-party-licenses.txt file in the root directory of this source
*  tree.
**********************************************************************/

/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2010, Rice University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Rice University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

#ifndef __REEDS_SHEPP_STATE_SPACE_HPP_
#define __REEDS_SHEPP_STATE_SPACE_HPP_

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "steering_functions/steering_functions.hpp"

using namespace std;
using namespace steer;

/** \brief An SE(2) state space where distance is measured by the
    length of Reeds-Shepp curves.
    The notation and solutions are taken from:
    J.A. Reeds and L.A. Shepp, “Optimal paths for a car that goes both
    forwards and backwards,” Pacific Journal of Mathematics,
    145(2):367–393, 1990.
    This implementation explicitly computes all 48 Reeds-Shepp curves
    and returns the shortest valid solution. This can be improved by
    using the configuration space partition described in:
    P. Souères and J.-P. Laumond, “Shortest paths synthesis for a
    car-like robot,” IEEE Trans. on Automatic Control, 41(5):672–688,
    May 1996.
    */
class Reeds_Shepp_State_Space
{
public:
  /** \brief The Reeds-Shepp path segment types */
  enum Reeds_Shepp_Path_Segment_Type
  {
    RS_NOP = 0,
    RS_LEFT = 1,
    RS_STRAIGHT = 2,
    RS_RIGHT = 3
  };

  /** \brief Reeds-Shepp path types */
  static const Reeds_Shepp_Path_Segment_Type reeds_shepp_path_type[18][5];

  /** \brief Complete description of a ReedsShepp path */
  class Reeds_Shepp_Path
  {
  public:
    /** \brief Constructor of the Reeds_Shepp_Path */
    Reeds_Shepp_Path(const Reeds_Shepp_Path_Segment_Type *type = reeds_shepp_path_type[0],
                     double t = numeric_limits<double>::max(), double u = 0., double v = 0., double w = 0.,
                     double x = 0.);

    double length() const
    {
      return total_length_;
    }

  private:
    /** Path segment types */
    const Reeds_Shepp_Path_Segment_Type *type_;

    /** Path segment lengths */
    double length_[5];

    /** Total length */
    double total_length_;

    friend class Reeds_Shepp_State_Space;
  };

  /** \brief Constructor of the Reeds_Shepp_State_Space */
  Reeds_Shepp_State_Space(double kappa, double discretization = 0.1) : kappa_(kappa), discretization_(discretization)
  {
    kappa_inv_ = 1 / kappa;
  }

  /** \brief Returns type and length of segments of path from state1 to state2 with curvature = 1.0 */
  Reeds_Shepp_Path reeds_shepp(const State &state1, const State &state2) const;

  /** \brief Returns shortest path length from state1 to state2 with curvature = kappa_ */
  double get_distance(const State &state1, const State &state2) const;

  /** \brief Returns controls of the shortest path from state1 to state2 with curvature = kappa_ */
  vector<Control> get_controls(const State &state1, const State &state2) const;

  /** \brief Returns shortest path from state1 to state2 with curvature = kappa_ */
  vector<State> get_path(const State &state1, const State &state2) const;

  /** \brief Returns integrated states given a start state and controls with curvature = kappa_ */
  vector<State> integrate(const State &state, const vector<Control> &controls) const;

  /** \brief Returns interpolated state at distance t in [0,1] (percent of total path length with curvature = kappa_) */
  State interpolate(const State &state, const vector<Control> &controls, double t) const;

  /** \brief Numeric integration using the forward euler method */
  inline State forward_euler(const State &state, const Control &control, double integration_step) const;

private:
  /** \brief Curvature */
  double kappa_;

  /** \brief Inverse of curvature */
  double kappa_inv_;

  /** \brief Discretization of path */
  double discretization_;
};

#endif
