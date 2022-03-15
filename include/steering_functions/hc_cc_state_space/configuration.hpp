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

#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

namespace steering
{

class Configuration
{
public:
  /** \brief Constructor */
  Configuration(double _x = 0.0, double _y = 0.0, double _theta = 0.0, double _kappa = 0.0);

  /** \brief Alphanumeric display */
  void print(bool eol) const;

  /** \brief Position */
  double x, y;

  /** \brief Orientation in rad between [0, 2*pi[ */
  double theta;

  /** \brief Curvature */
  double kappa;
};

/** \brief Cartesian distance between two configurations */
double configuration_distance(const Configuration &q1, const Configuration &q2);

/** \brief Are two configurations aligned? */
bool configuration_aligned(const Configuration &q1, const Configuration &q2);

/** \brief Are two configurations equal? */
bool configuration_equal(const Configuration &q1, const Configuration &q2);

} // namespace steering

#endif
