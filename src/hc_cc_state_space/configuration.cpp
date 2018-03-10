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

#include "steering_functions/hc_cc_state_space/configuration.hpp"

Configuration::Configuration(double _x, double _y, double _theta, double _kappa)
{
  x = _x;
  y = _y;
  theta = twopify(_theta);
  kappa = _kappa;
}

void Configuration::print(bool eol) const
{
  cout << "(" << x << ", " << y << ", " << theta << ", " << kappa << ")";
  if (eol)
  {
    cout << endl;
  }
}

double configuration_distance(const Configuration &q1, const Configuration &q2)
{
  return point_distance(q1.x, q1.y, q2.x, q2.y);
}

bool configuration_aligned(const Configuration &q1, const Configuration &q2)
{
  if (fabs(q2.theta - q1.theta) > get_epsilon())
  {
    return false;
  }
  double angle = twopify(atan2(q2.y - q1.y, q2.x - q1.x));
  return fabs(angle - q1.theta) <= get_epsilon();
}

bool configuration_equal(const Configuration &q1, const Configuration &q2)
{
  if (fabs(q2.theta - q1.theta) > get_epsilon())
    return false;
  if (configuration_distance(q1, q2) > get_epsilon())
    return false;
  return true;
}
