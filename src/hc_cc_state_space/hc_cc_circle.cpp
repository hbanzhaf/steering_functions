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

using namespace std;

namespace steering
{

void HC_CC_Circle_Param::set_param(double _kappa, double _sigma, double _radius, double _mu, double _sin_mu,
                                   double _cos_mu, double _delta_min)
{
  kappa = _kappa;
  kappa_inv = 1 / _kappa;
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
  left = _left;
  forward = _forward;
  regular = _regular;
  double delta_x = _param.radius * _param.sin_mu;
  double delta_y = _param.radius * _param.cos_mu;
  if (left)
  {
    kappa = _param.kappa;
    kappa_inv = _param.kappa_inv;
    sigma = _param.sigma;
    if (forward)
      global_frame_change(_start.x, _start.y, _start.theta, delta_x, delta_y, &xc, &yc);
    else
      global_frame_change(_start.x, _start.y, _start.theta, -delta_x, delta_y, &xc, &yc);
  }
  else
  {
    kappa = -_param.kappa;
    kappa_inv = -_param.kappa_inv;
    sigma = -_param.sigma;
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
    kappa = _param.kappa;
    kappa_inv = _param.kappa_inv;
    sigma = _param.sigma;
  }
  else
  {
    kappa = -_param.kappa;
    kappa_inv = -_param.kappa_inv;
    sigma = -_param.sigma;
  }
  xc = _xc;
  yc = _yc;
  radius = _param.radius;
  mu = _param.mu;
  sin_mu = _param.sin_mu;
  cos_mu = _param.cos_mu;
  delta_min = _param.delta_min;
}

double HC_CC_Circle::deflection(const Configuration &q) const
{
  double alpha_c = this->start.theta;
  double alpha_q = q.theta;
  if ((this->left && this->forward) || (!this->left && !this->forward))
  {
    return twopify(alpha_q - alpha_c);
  }
  else
  {
    return twopify(alpha_c - alpha_q);
  }
}

double HC_CC_Circle::D1(double alpha) const
{
  double fresnel_s, fresnel_c;
  double s = sqrt(2 * alpha / PI);
  fresnel(s, fresnel_s, fresnel_c);
  return cos(alpha) * fresnel_c + sin(alpha) * fresnel_s;
}

double HC_CC_Circle::rs_circular_deflection(double delta) const
{
  // regular rs-turn
  if (this->regular)
    return delta;
  // irregular rs-turn
  else
    return (delta <= PI) ? delta : delta - TWO_PI;
}

double HC_CC_Circle::rs_turn_length(const Configuration &q) const
{
  assert(fabs(fabs(this->kappa) - fabs(q.kappa)) < get_epsilon() &&
         fabs(fabs(this->sigma) - numeric_limits<double>::max()) < get_epsilon());
  double delta = this->deflection(q);
  return fabs(this->kappa_inv * this->rs_circular_deflection(delta));
}

double HC_CC_Circle::hc_circular_deflection(double delta) const
{
  double delta_min_twopified = twopify(this->delta_min);
  // regular hc-turn
  if (this->regular)
  {
    if (delta < delta_min_twopified)
      return TWO_PI + delta - delta_min_twopified;
    else
      return delta - delta_min_twopified;
  }
  // irregular hc-turn
  else
  {
    double delta_arc1, delta_arc2;
    if (delta < delta_min_twopified)
    {
      delta_arc1 = delta - delta_min_twopified;  // negative
      delta_arc2 = delta_arc1 + TWO_PI;          // positive
    }
    else
    {
      delta_arc1 = delta - delta_min_twopified;  // positive
      delta_arc2 = delta_arc1 - TWO_PI;          // negative
    }
    return (fabs(delta_arc1) < fabs(delta_arc2)) ? delta_arc1 : delta_arc2;
  }
}

double HC_CC_Circle::hc_turn_length(const Configuration &q) const
{
  assert(fabs(fabs(this->kappa) - fabs(q.kappa)) < get_epsilon());
  double delta = this->deflection(q);
  return fabs(this->kappa / this->sigma) + fabs(this->kappa_inv * this->hc_circular_deflection(delta));
}

bool HC_CC_Circle::cc_elementary_sharpness(const Configuration &q, double delta, double &sigma0) const
{
  double distance = point_distance(this->start.x, this->start.y, q.x, q.y);
  //   existence conditions for an elementary path (also see: A. Scheuer and T. Fraichard, "Continuous-Curvature Path
  //   Planning for Car-Like Vehicles," in IEEE/RSJ International Conference on Intelligent Robots and Systems, 1997.
  if (delta < 4.5948 && distance > get_epsilon())
  {
    sigma0 = 4 * PI * pow(this->D1(0.5 * delta), 2) / pow(distance, 2);
    if (!this->left)
    {
      sigma0 = -sigma0;
    }
    return true;
  }
  return false;
}

double HC_CC_Circle::cc_circular_deflection(double delta) const
{
  double two_delta_min_twopified = twopify(2 * this->delta_min);
  // regular cc-turn
  if (this->regular)
  {
    if (delta < two_delta_min_twopified)
      return TWO_PI + delta - two_delta_min_twopified;
    else
      return delta - two_delta_min_twopified;
  }
  // irregular cc-turn
  else
  {
    double delta_arc1, delta_arc2;
    if (delta < two_delta_min_twopified)
    {
      delta_arc1 = delta - two_delta_min_twopified;  // negative
      delta_arc2 = delta_arc1 + TWO_PI;              // positive
    }
    else
    {
      delta_arc1 = delta - two_delta_min_twopified;  // positive
      delta_arc2 = delta_arc1 - TWO_PI;              // negative
    }
    return (fabs(delta_arc1) < fabs(delta_arc2)) ? delta_arc1 : delta_arc2;
  }
}

double HC_CC_Circle::cc_turn_length(const Configuration &q) const
{
  assert(fabs(q.kappa) < get_epsilon());
  double delta = this->deflection(q);
  // delta = 0
  if (delta < get_epsilon())
  {
    return 2 * this->radius * this->sin_mu;
  }
  // 0 < delta < 2 * delta_min
  double length_min = fabs(this->kappa / this->sigma);
  double length_default = 2 * length_min + fabs(this->kappa_inv * this->cc_circular_deflection(delta));
  if (delta < 2 * this->delta_min)
  {
    double sigma0;
    if (this->cc_elementary_sharpness(q, delta, sigma0))
    {
      double length_elementary = 2 * sqrt(delta / fabs(sigma0));
      return (length_elementary < length_default) ? length_elementary : length_default;
    }
    else
    {
      return length_default;
    }
  }
  // delta >= 2 * delta_min
  else
  {
    return length_default;
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

} // namespace steering
