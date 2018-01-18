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

#include "steering_functions/hc_cc_state_space/utilities.hpp"

double get_epsilon()
{
  return epsilon;
}

double sgn(double x)
{
  return ((x < 0) ? -1.0 : 1.0);
}

double point_distance(double x1, double y1, double x2, double y2)
{
  return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

double twopify(double alpha)
{
  while (alpha >= TWO_PI)
  {
    alpha = alpha - (2 * PI);
  }
  while (alpha < 0)
  {
    alpha = alpha + (2 * PI);
  }
  return alpha;
}

double pify(double alpha)
{
  double v = fmod(alpha, TWO_PI);
  if (v < -PI)
    v += TWO_PI;
  else if (v > PI)
    v -= TWO_PI;
  return v;
}

double fresnel(double x, bool fresnelc)
{
  if (fabs(x) > Fresnel_Length)
  {
    cerr << "Fresnel integral out of range" << endl;
    return 0;
  }
  int sign = sgn(x);
  double diter = fabs(x) / Fresnel_Length * Fresnel_Samples;
  int iter = (int)diter;
  double rest = diter - iter;
  double inf, sup;
  if (fresnelc)
  {
    inf = FresnelC[iter];
    sup = FresnelC[iter + 1];
  }
  else
  {
    inf = FresnelS[iter];
    sup = FresnelS[iter + 1];
  }
  // interpolation
  return (sign * ((1 - rest) * inf + rest * sup));
}

double fresnelc(double s)
{
  return fresnel(s, true);
}

double fresnels(double s)
{
  return fresnel(s, false);
}

void end_of_clothoid(double x_i, double y_i, double theta_i, double kappa_i, double sigma, bool forward, double length,
                     double *x_f, double *y_f, double *theta_f, double *kappa_f)
{
  double x, y, theta, kappa, d;
  if (forward)
    d = 1;
  else
    d = -1;

  // assume initial configuration is at origin (x, y, theta = 0)
  if (fabs(sigma) < get_epsilon())
  {
    if (fabs(kappa_i) < get_epsilon())
    {
      x = d * length;
      y = 0;
      theta = 0;
      kappa = 0;
    }
    else
    {
      x = d * 1 / kappa_i * sin(kappa_i * length);
      y = 1 / kappa_i * (1 - cos(kappa_i * length));
      theta = d * kappa_i * length;
      kappa = kappa_i;
    }
  }
  else
  {
    // kappa(s) = kappa(0) + sigma * s
    kappa = kappa_i + sigma * length;
    // theta(s) = 1/2 * sigma * s^2 + kappa(0) * s
    theta = 0.5 * sigma * pow(length, 2) + kappa_i * length;
    // x(s) = int(cos(1/2 * sigma * u^2 + kappa * u), u = 0...s)
    // y(s) = int(sin(1/2 * sigma * u^2 + kappa * u), u= 0...s)
    int ssigma = sgn(sigma);
    double usigma = fabs(sigma);
    int skappa = sgn(kappa_i);
    double ukappa = fabs(kappa_i);
    double k1 = 0.5 * pow(ukappa, 2) / usigma;
    double k2 = (usigma * length + ssigma * skappa * ukappa) / sqrt(PI * usigma);
    double k3 = ukappa / sqrt(PI * usigma);
    x = sqrt(PI / usigma) * (cos(k1) * fresnelc(k2) + sin(k1) * fresnels(k2) -
                             ssigma * skappa * cos(k1) * fresnelc(k3) - ssigma * skappa * sin(k1) * fresnels(k3));
    y = sqrt(PI / usigma) * (ssigma * cos(k1) * fresnels(k2) - ssigma * sin(k1) * fresnelc(k2) -
                             skappa * cos(k1) * fresnels(k3) + skappa * sin(k1) * fresnelc(k3));
    x = d * x;
    theta = d * theta;
  }

  // translation and rotation to account for initial configuration
  global_frame_change(x_i, y_i, theta_i, x, y, x_f, y_f);
  *theta_f = pify(theta_i + theta);
  *kappa_f = kappa;
}

void end_of_circular_arc(double x_i, double y_i, double theta_i, double kappa, bool forward, double length, double *x_f,
                         double *y_f, double *theta_f)
{
  double x, y, theta, d;
  if (forward)
    d = 1;
  else
    d = -1;

  // assume initial configuration is at origin (x, y, theta = 0)
  if (fabs(kappa) < get_epsilon())
  {
    x = d * length;
    y = 0;
    theta = 0;
  }
  else
  {
    x = d * 1 / kappa * sin(kappa * length);
    y = 1 / kappa * (1 - cos(kappa * length));
    theta = d * kappa * length;
  }

  // translation and rotation to account for initial configuration
  *x_f = x * cos(theta_i) - y * sin(theta_i) + x_i;
  *y_f = x * sin(theta_i) + y * cos(theta_i) + y_i;
  *theta_f = pify(theta_i + theta);
}

void global_frame_change(double x, double y, double theta, double local_x, double local_y, double *global_x,
                         double *global_y)
{
  double sin_th = sin(theta);
  double cos_th = cos(theta);
  *global_x = local_x * cos_th - local_y * sin_th + x;
  *global_y = local_x * sin_th + local_y * cos_th + y;
}

void local_frame_change(double x, double y, double theta, double global_x, double global_y, double *local_x,
                        double *local_y)
{
  double sin_th = sin(theta);
  double cos_th = cos(theta);
  *local_x = (global_x - x) * cos_th + (global_y - y) * sin_th;
  *local_y = -(global_x - x) * sin_th + (global_y - y) * cos_th;
}

double D1(double alpha)
{
  return cos(alpha) * fresnelc(sqrt(2 * alpha / PI)) + sin(alpha) * fresnels(sqrt(2 * alpha / PI));
}

int array_index_min(double array[], int size)
{
  double min = array[0];
  int index_min = 0;
  for (int i = 1; i < size; i++)
  {
    if (array[i] < min)
    {
      index_min = i;
      min = array[i];
    }
  }
  return index_min;
}

void double_array_init(double array[], int size, double value)
{
  for (int i = 0; i < size; i++)
  {
    array[i] = value;
  }
}

void pointer_array_init(void *array[], int size)
{
  for (int i = 0; i < size; i++)
  {
    array[i] = nullptr;
  }
}
