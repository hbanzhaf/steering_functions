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
***********************************************************************/

#ifndef STEERING_FUNCTIONS_HPP
#define STEERING_FUNCTIONS_HPP

namespace steer
{
/** \brief Description of a kinematic car's state */
struct State
{
  /** \brief Position in x of the robot */
  double x;

  /** \brief Position in y of the robot */
  double y;

  /** \brief Orientation of the robot */
  double theta;

  /** \brief Curvature at position (x,y) */
  double kappa;

  /** \Driving direction {-1,0,1} */
  double d;
};

/** \brief Description of a kinematic car's state with covariance */
struct State_With_Covariance
{
  /** \brief Expected state of the robot */
  State state;

  /** \brief Covariance of the state: (x_x      x_y      x_theta      x_kappa
                                       y_x      y_y      y_theta      y_kappa
                                       theta_x  theta_y  theta_theta  theta_kappa
                                       kappa_x  kappa_y  kappa_theta  kappa_kappa) */
  double covariance[16] = { 0.0 };
};

/** \brief Description of a path segment with its corresponding control inputs */
struct Control
{
  /** \brief Signed arc length of a segment */
  double delta_s;

  /** \brief Curvature at the beginning of a segment */
  double kappa;

  /** \brief Sharpness (derivative of curvature with respect to arc length) of a segment */
  double sigma;
};

/** \brief Parameters of the motion noise model according to the book:
    Probabilistic Robotics, S. Thrun and others, MIT Press, 2006, p. 127-128 and p.204-206. */
struct Motion_Noise
{
  /** \brief Variance in longitudinal direction: alpha1*delta_s*delta_s + alpha2*kappa*kappa + alpha3*sigma*sigma  */
  double alpha1;
  double alpha2;
  double alpha3;

  /** \brief Variance in lateral direction: alpha4*delta_s*delta_s + alpha5*kappa*kappa + alpha6*sigma*sigma  */
  double alpha4;
  double alpha5;
  double alpha6;
};

/** \brief Parameters of the measurement noise */
struct Measurement_Noise
{
  /** \brief Standard deviation of localization in x */
  double std_x;

  /** \brief Standard deviation of localization in y */
  double std_y;

  /** \brief Standard deviation of localization in theta */
  double std_theta;

  /** \brief Standard deviation of curvature measurement */
  double std_kappa;
};

/** \brief Parameters of the feedback controller */
struct Controller
{
  /** \brief Weight on longitudinal error */
  double k1;

  /** \brief Weight on lateral error */
  double k2;

  /** \brief Weight on heading error */
  double k3;

  /** \brief Weight on curvature error */
  double k4;
};
}

#endif
