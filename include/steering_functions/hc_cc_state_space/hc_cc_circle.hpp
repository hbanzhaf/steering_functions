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

#ifndef HC_CC_CIRCLE_HPP
#define HC_CC_CIRCLE_HPP

#include "steering_functions/hc_cc_state_space/configuration.hpp"

namespace steering
{

class HC_CC_Circle_Param
{
public:
  /** \brief Set parameters */
  void set_param(double _kappa, double _sigma, double _radius, double _mu, double _sin_mu, double _cos_mu,
                 double _delta_min);

  /** \brief Max. curvature, inverse of max. curvature, max. sharpness */
  double kappa, kappa_inv, sigma;

  /** \brief Radius of the outer circle */
  double radius;

  /** \brief Angle between the initial orientation and the tangent to the circle at the initial position */
  double mu;

  /** \brief Sine and cosine of mu */
  double sin_mu, cos_mu;

  /** \brief Minimal deflection */
  double delta_min;
};

class HC_CC_Circle : public HC_CC_Circle_Param
{
public:
  /** \brief Constructor */
  HC_CC_Circle(const Configuration &_start, bool _left, bool _forward, bool _regular, const HC_CC_Circle_Param &_param);

  /** \brief Constructor */
  HC_CC_Circle(double _xc, double _yc, bool _left, bool _forward, bool _regular, const HC_CC_Circle_Param &_param);

  /** \brief Computation of deflection (angle between start configuration of circle and configuration q) */
  double deflection(const Configuration &q) const;

  /** \brief Calculation of D1 for the evaluation of an elementary path */
  double D1(double alpha) const;

  /** \brief Computation of a rs-turn's circular deflection */
  double rs_circular_deflection(double delta) const;

  /** \brief Length of a rs-turn */
  double rs_turn_length(const Configuration &q) const;

  /** \brief Computation of a hc-turn's circular deflection */
  double hc_circular_deflection(double delta) const;

  /** \brief Length of a hc-turn */
  double hc_turn_length(const Configuration &q) const;

  /** \brief Computation of an elementary path's sharpness */
  bool cc_elementary_sharpness(const Configuration &q, double delta, double &sigma0) const;

  /** \brief Computation of a cc-turn's circular deflection */
  double cc_circular_deflection(double delta) const;

  /** \brief Length of a cc-turn */
  double cc_turn_length(const Configuration &q) const;

  /** \brief Alphanumeric display */
  void print(bool eol) const;

  /** \brief Start configuration */
  Configuration start;

  /** \brief Turning direction: left/right */
  bool left;

  /** \brief Driving direction: forwards/backwards */
  bool forward;

  /** \brief Type of the circle: regular/irregular */
  bool regular;

  /** \brief Center of the circle */
  double xc, yc;
};

/** \brief Cartesian distance between the centers of two circles */
double center_distance(const HC_CC_Circle &c1, const HC_CC_Circle &c2);

/** \brief Configuration on the circle? */
bool configuration_on_hc_cc_circle(const HC_CC_Circle &c, const Configuration &q);

} // namespace steering

#endif
