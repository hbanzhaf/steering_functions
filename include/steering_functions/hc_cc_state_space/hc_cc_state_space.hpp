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

#include "steering_functions/filter/ekf.hpp"
#include "steering_functions/hc_cc_state_space/hc_cc_circle.hpp"
#include "steering_functions/steering_functions.hpp"
#include "steering_functions/utilities/utilities.hpp"

namespace steering
{

class HC_CC_State_Space
{
public:
  /** \brief Constructor */
  HC_CC_State_Space(double kappa, double sigma, double discretization);

  /** \brief Sets the parameters required by the filter */
  void set_filter_parameters(const Motion_Noise& motion_noise, const Measurement_Noise& measurement_noise,
                             const Controller& controller);

  /** \brief Virtual function that returns controls of the shortest path from state1 to state2 */
  virtual std::vector<Control> get_controls(const State& state1, const State& state2) const = 0;

  /** \brief Returns path from state1 to state2 */
  std::vector<State> get_path(const State& state1, const State& state2) const;

  /** \brief Returns path including covariances from state1 to state2 */
  std::vector<State_With_Covariance> get_path_with_covariance(const State_With_Covariance& state1,
                                                              const State& state2) const;

  /** \brief Returns integrated states given a start state and controls */
  std::vector<State> integrate(const State& state, const std::vector<Control>& controls) const;

  /** \brief Returns integrated states including covariance given a start state and controls */
  std::vector<State_With_Covariance> integrate_with_covariance(const State_With_Covariance& state,
                                                               const std::vector<Control>& controls) const;

  /** \brief Returns interpolated state at distance t in [0,1] (percentage of total path length) */
  State interpolate(const State& state, const std::vector<Control>& controls, double t) const;

  /** \brief Returns integrated state given a start state, a control, and an integration step */
  inline State integrate_ODE(const State& state, const Control& control, double integration_step) const;

protected:
  /** \brief Curvature, sharpness of clothoid */
  double kappa_, sigma_;

  /** \brief Discretization of path */
  double discretization_;

  /** \brief Parameters of a hc-/cc-circle */
  HC_CC_Circle_Param hc_cc_circle_param_;

  /** \brief Extended Kalman Filter for uncertainty propagation */
  EKF ekf_;
};

} // namespace steering

#endif
