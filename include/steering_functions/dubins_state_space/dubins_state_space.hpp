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
*  Software License Agreement (BSD License)
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

#ifndef __DUBINS_STATE_SPACE_HPP_
#define __DUBINS_STATE_SPACE_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "steering_functions/filter/ekf.hpp"
#include "steering_functions/steering_functions.hpp"
#include "steering_functions/utilities/utilities.hpp"

namespace steering
{

/** \brief An SE(2) state space where distance is measured by the
    length of Dubins curves.
    Note that this Dubins distance is \b not a proper distance metric,
    so nearest neighbor methods that rely on distance() being a metric
    (such as ompl::NearestNeighborsGNAT) will not always return the
    true nearest neighbors or get stuck in an infinite loop.
    The notation and solutions in the code are taken from:<br>
    A.M. Shkel and V. Lumelsky, “Classification of the Dubins set,”
    Robotics and Autonomous Systems, 34(4):179-202, 2001.
    DOI: <a href="http://dx.doi.org/10.1016/S0921-8890(00)00127-5">10.1016/S0921-8890(00)00127-5</a>
    The classification scheme described there is not actually used,
    since it only applies to “long” paths.
    */
class Dubins_State_Space
{
public:
  /** \brief The Dubins path segment type */
  enum Dubins_Path_Segment_Type
  {
    DUBINS_LEFT = 0,
    DUBINS_STRAIGHT = 1,
    DUBINS_RIGHT = 2
  };

  /** \brief Dubins path types */
  static const Dubins_Path_Segment_Type dubins_path_type[6][3];

  /** \brief Complete description of a Dubins path */
  class Dubins_Path
  {
  public:
    /** \brief Constructor of the Dubins_Path */
    Dubins_Path(const Dubins_Path_Segment_Type *type = dubins_path_type[0], double t = 0.,
                double p = std::numeric_limits<double>::max(), double q = 0.)
      : type_(type)
    {
      length_[0] = t;
      length_[1] = p;
      length_[2] = q;
      assert(t >= 0.);
      assert(p >= 0.);
      assert(q >= 0.);
    }

    double length() const
    {
      return length_[0] + length_[1] + length_[2];
    }

  private:
    /** Path segment types */
    const Dubins_Path_Segment_Type *type_;

    /** Path segment lengths */
    double length_[3];

    friend class Dubins_State_Space;
  };

  /** \brief Constructor of the Dubins_State_Space */
  Dubins_State_Space(double kappa, double discretization = 0.1, bool forwards = true)
    : kappa_(kappa), discretization_(discretization), forwards_(forwards)
  {
    assert(kappa > 0.0 && discretization > 0.0);
    kappa_inv_ = 1 / kappa;
  }

  /** \brief Sets the parameters required by the filter */
  void set_filter_parameters(const Motion_Noise &motion_noise, const Measurement_Noise &measurement_noise,
                             const Controller &controller);

  /** \brief Returns type and length of segments of path from state1 to state2 with curvature = 1.0 */
  Dubins_Path dubins(const State &state1, const State &state2) const;

  /** \brief Returns shortest path length from state1 to state2 with curvature = kappa_ */
  double get_distance(const State &state1, const State &state2) const;

  /** \brief Returns controls of the shortest path from state1 to state2 with curvature = kappa_ */
  std::vector<Control> get_controls(const State &state1, const State &state2) const;

  /** \brief Returns shortest path from state1 to state2 with curvature = kappa_ */
  std::vector<State> get_path(const State &state1, const State &state2) const;

  /** \brief Returns shortest path including covariances from state1 to state2 with curvature = kappa_ */
  std::vector<State_With_Covariance> get_path_with_covariance(const State_With_Covariance &state1,
                                                              const State &state2) const;

  /** \brief Returns integrated states given a start state and controls with curvature = kappa_ */
  std::vector<State> integrate(const State &state, const std::vector<Control> &controls) const;

  /** \brief Returns integrated states including covariance given a start state and controls with curvature = kappa_ */
  std::vector<State_With_Covariance> integrate_with_covariance(const State_With_Covariance &state,
                                                               const std::vector<Control> &controls) const;

  /** \brief Returns interpolated state at distance t in [0,1] (percent of total path length with curvature = kappa_) */
  State interpolate(const State &state, const std::vector<Control> &controls, double t) const;

  /** \brief Returns integrated state given a start state, a control, and an integration step */
  inline State integrate_ODE(const State &state, const Control &control, double integration_step) const;

private:
  /** \brief Curvature */
  double kappa_;

  /** \brief Inverse of curvature */
  double kappa_inv_;

  /** \brief Discretization of path */
  double discretization_;

  /** \brief Driving direction */
  bool forwards_;

  /** \brief Extended Kalman Filter for uncertainty propagation */
  EKF ekf_;
};

} // namespace steering

#endif
