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

#include "steering_functions/hc_cc_state_space/hc_cc_circle.hpp"

void HC_CC_Circle_Param::set_param(double _kappa, double _sigma, double _radius, double _mu, double _sin_mu,
                                   double _cos_mu, double _delta_min)
{
  kappa = _kappa;
  sigma = _sigma;
  radius = _radius;
  mu = _mu;
  sin_mu = _sin_mu;
  cos_mu = _cos_mu;
  delta_min = _delta_min;
}

HC_CC_Circle::HC_CC_Circle(const Configuration &_start, bool _left, bool _forward, bool _regular,
                           const HC_CC_Circle_Param &_param)
{
  start = _start;
  start.kappa = 0;
  left = _left;
  forward = _forward;
  regular = _regular;
  double delta_x = _param.radius * _param.sin_mu;
  double delta_y = _param.radius * _param.cos_mu;
  if (left)
  {
    kappa = fabs(_param.kappa);
    sigma = fabs(_param.sigma);
    if (forward)
      global_frame_change(_start.x, _start.y, _start.theta, delta_x, delta_y, &xc, &yc);
    else
      global_frame_change(_start.x, _start.y, _start.theta, -delta_x, delta_y, &xc, &yc);
  }
  else
  {
    kappa = -fabs(_param.kappa);
    sigma = -fabs(_param.sigma);
    if (forward)
      global_frame_change(_start.x, _start.y, _start.theta, delta_x, -delta_y, &xc, &yc);
    else
      global_frame_change(_start.x, _start.y, _start.theta, -delta_x, -delta_y, &xc, &yc);
  }
  radius = _param.radius;
  mu = _param.mu;
  sin_mu = _param.sin_mu;
  cos_mu = _param.cos_mu;
  delta_min = _param.delta_min;
}

HC_CC_Circle::HC_CC_Circle(double _xc, double _yc, bool _left, bool _forward, bool _regular,
                           const HC_CC_Circle_Param &_param)
{
  start = Configuration(0, 0, 0, 0);
  left = _left;
  forward = _forward;
  regular = _regular;
  if (left)
  {
    kappa = fabs(_param.kappa);
    sigma = fabs(_param.sigma);
  }
  else
  {
    kappa = -fabs(_param.kappa);
    sigma = -fabs(_param.sigma);
  }
  xc = _xc;
  yc = _yc;
  radius = _param.radius;
  mu = _param.mu;
  sin_mu = _param.sin_mu;
  cos_mu = _param.cos_mu;
  delta_min = _param.delta_min;
}

void HC_CC_Circle::deflection(const Configuration &q, double *delta) const
{
  double alpha_c = this->start.theta;
  double alpha_q = q.theta;
  if (this->left && this->forward)
  {
    *delta = twopify(alpha_q - alpha_c);
  }
  if (this->left && !this->forward)
  {
    *delta = twopify(alpha_c - alpha_q);
  }
  if (!this->left && this->forward)
  {
    *delta = twopify(alpha_c - alpha_q);
  }
  if (!this->left && !this->forward)
  {
    *delta = twopify(alpha_q - alpha_c);
  }
}

double HC_CC_Circle::rs_turn_length(const Configuration &q) const
{
  assert(fabs(this->kappa) - fabs(q.kappa) < get_epsilon() &&
         fabs(this->sigma) - numeric_limits<double>::max() < get_epsilon());
  double delta;
  this->deflection(q, &delta);
  // irregular rs-turn
  if (!this->regular && (delta > PI))
  {
    return fabs((TWO_PI - delta) / this->kappa);
  }
  // regular rs-turn
  else
  {
    return fabs(delta / this->kappa);
  }
}

double HC_CC_Circle::hc_turn_length(const Configuration &q) const
{
  assert(fabs(this->kappa) - fabs(q.kappa) < get_epsilon());
  double delta;
  this->deflection(q, &delta);
  double length_min = fabs(this->kappa / this->sigma);
  double length_arc;
  // regular hc-turn
  if (this->regular && (delta < delta_min / 2.0))
  {
    length_arc = fabs((TWO_PI + delta - delta_min / 2.0) / this->kappa);
  }
  // irregular hc-turn
  else if (!this->regular && (delta < delta_min / 2.0))
  {
    length_arc = fabs((-delta + delta_min / 2.0) / this->kappa);
  }
  // irregular hc-turn
  else if (!this->regular && (delta > delta_min / 2.0 + PI))
  {
    length_arc = fabs((TWO_PI - delta + delta_min / 2.0) / this->kappa);
  }
  // regular hc-turn
  else
  {
    length_arc = fabs((delta - delta_min / 2.0) / this->kappa);
  }
  return length_min + length_arc;
}

double HC_CC_Circle::cc_turn_length(const Configuration &q) const
{
  assert(fabs(q.kappa) < get_epsilon());
  double delta;
  this->deflection(q, &delta);
  // straight line
  if (fabs(delta) < get_epsilon())
  {
    return fabs(2 * this->radius * this->sin_mu);
  }
  // elementary path
  if (delta < delta_min)
  {
    double d1 = D1(delta / 2);
    double d2 = point_distance(this->start.x, this->start.y, q.x, q.y);
    double sharpness = 4 * PI * pow(d1, 2) / pow(d2, 2);
    return 2 * sqrt(delta / sharpness);
  }
  double length_min = fabs(2 * this->kappa / this->sigma);
  // irregular cc-turn
  if (!this->regular && (delta > delta_min + PI))
  {
    return length_min + fabs((TWO_PI - delta + delta_min) / this->kappa);
  }
  // regular cc-turn
  else
  {
    return length_min + fabs((delta - delta_min) / this->kappa);
  }
}

void HC_CC_Circle::print(bool eol) const
{
  cout << "HC_CC_Circle: ";
  cout << "start: ";
  start.print(false);
  if (left)
  {
    cout << ", left";
  }
  else
  {
    cout << ", right";
  }
  if (forward)
  {
    cout << ", forward";
  }
  else
  {
    cout << ", backward";
  }
  if (regular)
  {
    cout << ", regular";
  }
  else
  {
    cout << ", irregular";
  }
  cout << ", kappa: " << kappa << ", sigma: " << sigma;
  cout << ", centre: (" << xc << ", " << yc << "), radius " << radius << ", mu: " << mu;
  if (eol)
  {
    cout << endl;
  }
}

double center_distance(const HC_CC_Circle &c1, const HC_CC_Circle &c2)
{
  return sqrt(pow(c2.xc - c1.xc, 2) + pow(c2.yc - c1.yc, 2));
}

bool configuration_on_hc_cc_circle(const HC_CC_Circle &c, const Configuration &q)
{
  double distance = point_distance(c.xc, c.yc, q.x, q.y);
  if (fabs(distance - c.radius) > get_epsilon())
  {
    return false;
  }
  double angle = atan2(q.y - c.yc, q.x - c.xc);
  if (c.left && c.forward)
  {
    angle = angle + HALF_PI - c.mu;
  }
  if (c.left && !c.forward)
  {
    angle = angle + HALF_PI + c.mu;
  }
  if (!c.left && c.forward)
  {
    angle = angle - HALF_PI + c.mu;
  }
  if (!c.left && !c.forward)
  {
    angle = angle - HALF_PI - c.mu;
  }
  angle = twopify(angle);
  return fabs(q.theta - angle) < get_epsilon();
}
