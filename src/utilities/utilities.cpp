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

#include <cassert>
#include <cmath>

#include "steering_functions/utilities/utilities.hpp"

namespace steering
{

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

void polar(double x, double y, double &r, double &theta)
{
  r = sqrt(x * x + y * y);
  theta = atan2(y, x);
}

double twopify(double alpha)
{
  return alpha - TWO_PI * floor(alpha / TWO_PI);
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

void fresnel_0_8(double x, double &S_f, double &C_f)
{
  // T2n(x/8) = Tn(2*(x/8)*(x/8)-1.0)
  // T2n_p1(x/8) = 2*(x/8)*T2n(x/8)-T2n_m1(x/8)
  double quarter_x = 0.25 * x;
  double arg = 0.03125 * x * x - 1.0;
  double T0 = 1;
  double T1 = 0.125 * x;
  double T2 = arg;
  double T3 = quarter_x * T2 - T1;
  double A = chebev_a[0] * T0 + chebev_a[1] * T2;
  double B = chebev_b[0] * T1 + chebev_b[1] * T3;
  double T2n_m4 = T0;
  double T2n_m2 = T2;
  double T2n_m1 = T3;
  for (int n = 2; n < 17; n++)
  {
    double T2n = 2.0 * arg * T2n_m2 - T2n_m4;
    double T2n_p1 = quarter_x * T2n - T2n_m1;
    A += chebev_a[n] * T2n;
    B += chebev_b[n] * T2n_p1;
    T2n_m4 = T2n_m2;
    T2n_m2 = T2n;
    T2n_m1 = T2n_p1;
  }
  double T34 = 2.0 * arg * T2n_m2 - T2n_m4;
  A += chebev_a[17] * T34;

  double sqrt_x = sqrt(x);
  C_f = SQRT_TWO_PI_INV * sqrt_x * A;
  S_f = SQRT_TWO_PI_INV * sqrt_x * B;
  return;
}

void fresnel_8_inf(double x, double &S_f, double &C_f)
{
  // T2n(8/x) = Tn(2*(8/x)*(8/x)-1.0)
  double arg = 128.0 / (x * x) - 1.0;
  double T0 = 1;
  double T2 = arg;
  double E = chebev_e[0] * T0 + chebev_e[1] * T2;
  double F = chebev_f[0] * T0 + chebev_f[1] * T2;
  double T2n_m4 = T0;
  double T2n_m2 = T2;
  for (int n = 2; n < 35; n++)
  {
    double T2n = 2.0 * arg * T2n_m2 - T2n_m4;
    E += chebev_e[n] * T2n;
    F += chebev_f[n] * T2n;
    T2n_m4 = T2n_m2;
    T2n_m2 = T2n;
  }
  for (int n = 35; n < 41; n++)
  {
    double T2n = 2.0 * arg * T2n_m2 - T2n_m4;
    E += chebev_e[n] * T2n;
    T2n_m4 = T2n_m2;
    T2n_m2 = T2n;
  }

  double sin_x = sin(x);
  double cos_x = cos(x);
  double sqrt_x = sqrt(x);
  C_f = 0.5 - SQRT_TWO_PI_INV * (E * cos_x / (2 * x) - F * sin_x) / sqrt_x;
  S_f = 0.5 - SQRT_TWO_PI_INV * (E * sin_x / (2 * x) + F * cos_x) / sqrt_x;
  return;
}

void fresnel(double s, double &S_f, double &C_f)
{
  double x = HALF_PI * s * s;
  if (x <= 8.0)
    fresnel_0_8(x, S_f, C_f);
  else
    fresnel_8_inf(x, S_f, C_f);
  if (s < 0)
  {
    S_f = -S_f;
    C_f = -C_f;
  }
  return;
}

void end_of_clothoid(double x_i, double y_i, double theta_i, double kappa_i, double sigma, double direction,
                     double length, double *x_f, double *y_f, double *theta_f, double *kappa_f)
{
  // x_f = x_i + int_0_length(cos(theta_i + kappa_i*s + 0.5*sigma_i*s^2)ds)
  // y_f = y_i + int_0_length(sin(theta_i + kappa_i*s + 0.5*sigma_i*s^2)ds)
  double sgn_sigma = sgn(sigma);
  double abs_sigma = fabs(sigma);
  double sqrt_sigma_inv = 1 / sqrt(abs_sigma);
  double k1 = theta_i - 0.5 * direction * kappa_i * kappa_i / sigma;
  double k2 = SQRT_PI_INV * sqrt_sigma_inv * (abs_sigma * length + sgn_sigma * kappa_i);
  double k3 = SQRT_PI_INV * sqrt_sigma_inv * sgn_sigma * kappa_i;
  double cos_k1 = cos(k1);
  double sin_k1 = sin(k1);
  double fresnel_s_k2;
  double fresnel_c_k2;
  double fresnel_s_k3;
  double fresnel_c_k3;
  fresnel(k2, fresnel_s_k2, fresnel_c_k2);
  fresnel(k3, fresnel_s_k3, fresnel_c_k3);
  *x_f = x_i +
         SQRT_PI * sqrt_sigma_inv *
             (direction * cos_k1 * (fresnel_c_k2 - fresnel_c_k3) - sgn_sigma * sin_k1 * (fresnel_s_k2 - fresnel_s_k3));
  *y_f = y_i +
         SQRT_PI * sqrt_sigma_inv *
             (direction * sin_k1 * (fresnel_c_k2 - fresnel_c_k3) + sgn_sigma * cos_k1 * (fresnel_s_k2 - fresnel_s_k3));

  // theta_f = theta_i + kappa_i*length + 0.5*sigma*length^2
  *theta_f = pify(theta_i + kappa_i * direction * length + 0.5 * sigma * direction * length * length);
  // kappa_f = kappa_i + sigma * d * length
  *kappa_f = kappa_i + sigma * length;
}

void end_of_circular_arc(double x_i, double y_i, double theta_i, double kappa, double direction, double length,
                         double *x_f, double *y_f, double *theta_f)
{
  *x_f = x_i + (1 / kappa) * (-sin(theta_i) + sin(theta_i + direction * length * kappa));
  *y_f = y_i + (1 / kappa) * (cos(theta_i) - cos(theta_i + direction * length * kappa));
  *theta_f = pify(theta_i + kappa * direction * length);
}

void end_of_straight_line(double x_i, double y_i, double theta, double direction, double length, double *x_f,
                          double *y_f)
{
  *x_f = x_i + direction * length * cos(theta);
  *y_f = y_i + direction * length * sin(theta);
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

} // namespace steering
