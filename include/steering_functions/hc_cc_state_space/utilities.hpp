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

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cassert>
#include <cmath>
#include <iostream>

#include "steering_functions/hc_cc_state_space/fresnel.data"

using namespace std;

#define PI 3.1415926535897932384
#define HALF_PI 1.5707963267948966192
#define TWO_PI 6.2831853071795864770

const double epsilon = 1e-4;

/** \brief Return value of epsilon */
double get_epsilon();

/** \brief Return sign of a number */
double sgn(double x);

/** \brief Cartesian distance between two points */
double point_distance(double x1, double y1, double x2, double y2);

/** \brief Conversion of arbitrary angle given in [rad] to [0, 2*pi[ */
double twopify(double alpha);

/** \brief Conversion of arbitrary angle given in [rad] to [-pi, pi[ */
double pify(double alpha);

/** \brief Computation of the factorial */
unsigned int factorial(unsigned int n);

/** \brief Approximation of sine with a 7th order polynomial (absolute error is bounded by 1.7e-4 in [-pi/2, pi/2]) */
float approxSin(const float x);

/** \brief Approximation of cosine with a 8th order polynomial (absolute error is bounded by 2.6e-5 in [-pi/2, pi/2]) */
float approxCos(const float x);

/** \brief Fresnel integrals */
double fresnelc(double s);
double fresnels(double s);

/** \brief Computation of the end point on a clothoid
    x_i, y_i, theta_i, kappa_i: initial configuration
    sigma: sharpness of clothoid
    forward: driving direction
    length: length of clothoid (positive)
    x_f, y_f, theta_f, kappa_f: final configuration on clothoid
    */
void end_of_clothoid(double x_i, double y_i, double theta_i, double kappa_i, double sigma, bool forward, double length,
                     double *x_f, double *y_f, double *theta_f, double *kappa_f);

/** \brief Computation of the end point on a circular arc
    x_i, y_i, theta_i, kappa_i: initial configuration
    kappa: curvature of circular arc
    forward: driving direction
    length: arc length (positive)
    x_f, y_f, theta_f, kappa_f: final configuration on circular arc
    */
void end_of_circular_arc(double x_i, double y_i, double theta_i, double kappa, bool forward, double length, double *x_f,
                         double *y_f, double *theta_f);

/** \brief Transformation of (local_x, local_y) from local coordinate system to global one */
void global_frame_change(double x, double y, double theta, double local_x, double local_y, double *global_x,
                         double *global_y);

/** \brief Transformation of (global_x, global_y) from global coordinate system to local one */
void local_frame_change(double x, double y, double theta, double global_x, double global_y, double *local_x,
                        double *local_y);

/** \brief Calculation of D1 for an angle alpha between [0, pi] */
double D1(double alpha);

/** \brief Find index with minimal value in double array */
int array_index_min(double array[], int size);

/** \brief Initialize an array with a given value */
void double_array_init(double array[], int size, double value);

/** \brief Initialize an array with nullptr */
void pointer_array_init(void *array[], int size);

#endif
