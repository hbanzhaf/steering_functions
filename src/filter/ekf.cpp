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
* *********************************************************************/

#include <cmath>

#include "steering_functions/filter/ekf.hpp"
#include "steering_functions/utilities/utilities.hpp"

namespace steering
{

EKF::EKF()
{
  I_ = Matrix3d::Identity();
}

void EKF::set_parameters(const Motion_Noise &motion_noise, const Measurement_Noise &measurement_noise,
                         const Controller &controller)
{
  motion_noise_ = motion_noise;
  measurement_noise_ = measurement_noise;
  controller_ = controller;
}

Matrix3d EKF::covariance_to_eigen(const double covariance[16]) const
{
  Matrix3d covariance_eigen;
  for (int i = 0, n = covariance_eigen.rows(); i < n; ++i)
  {
    for (int j = 0, m = covariance_eigen.cols(); j < m; ++j)
    {
      covariance_eigen(i, j) = covariance[(i << 2) + j];  // bitwise operation for: covariance[i * 4 + j]
    }
  }
  return covariance_eigen;
}

void EKF::eigen_to_covariance(const Matrix3d &covariance_eigen, double covariance[16]) const
{
  for (int i = 0, n = covariance_eigen.rows(); i < n; ++i)
  {
    for (int j = 0, m = covariance_eigen.cols(); j < m; ++j)
    {
      covariance[(i << 2) + j] = covariance_eigen(i, j);  // bitwise operation for: covariance[i * 4 + j]
    }
  }
}

void EKF::get_motion_jacobi(const State &state, const Control &control, double integration_step, Matrix3d &F_x,
                            Matrix32d &F_u) const
{
  double d(sgn(control.delta_s));

  if (fabs(control.sigma) > get_epsilon())
  {
    double sgn_sigma = sgn(control.sigma);
    double abs_sigma = fabs(control.sigma);
    double sqrt_sigma_inv = 1 / sqrt(abs_sigma);
    double k1 = state.theta - 0.5 * d * state.kappa * state.kappa / control.sigma;
    double k2 = SQRT_PI_INV * sqrt_sigma_inv * (abs_sigma * integration_step + sgn_sigma * state.kappa);
    double k3 = SQRT_PI_INV * sqrt_sigma_inv * sgn_sigma * state.kappa;
    double k4 = k1 + d * sgn_sigma * HALF_PI * k2 * k2;
    double k5 = k1 + d * sgn_sigma * HALF_PI * k3 * k3;
    double cos_k1 = cos(k1);
    double sin_k1 = sin(k1);
    double cos_k4 = cos(k4);
    double sin_k4 = sin(k4);
    double cos_k5 = cos(k5);
    double sin_k5 = sin(k5);
    double fresnel_s_k2;
    double fresnel_c_k2;
    double fresnel_s_k3;
    double fresnel_c_k3;
    fresnel(k2, fresnel_s_k2, fresnel_c_k2);
    fresnel(k3, fresnel_s_k3, fresnel_c_k3);

    // df/dx
    F_x(0, 0) = 1.0;
    F_x(1, 1) = 1.0;

    F_x(0, 2) = SQRT_PI * sqrt_sigma_inv *
                (-d * sin_k1 * (fresnel_c_k2 - fresnel_c_k3) - sgn_sigma * cos_k1 * (fresnel_s_k2 - fresnel_s_k3));
    F_x(1, 2) = SQRT_PI * sqrt_sigma_inv *
                (d * cos_k1 * (fresnel_c_k2 - fresnel_c_k3) - sgn_sigma * sin_k1 * (fresnel_s_k2 - fresnel_s_k3));
    F_x(2, 2) = 1.0;

    // df/du
    F_u(0, 0) = cos_k4;
    F_u(1, 0) = sin_k4;
    F_u(2, 0) = state.kappa + control.sigma * integration_step;

    F_u(0, 1) = SQRT_PI * sqrt_sigma_inv * state.kappa *
                    (sgn_sigma * sin_k1 * (fresnel_c_k2 - fresnel_c_k3) + d * cos_k1 * (fresnel_s_k2 - fresnel_s_k3)) /
                    abs_sigma +
                d * (cos_k4 - cos_k5) / control.sigma;
    F_u(1, 1) = SQRT_PI * sqrt_sigma_inv * state.kappa *
                    (-sgn_sigma * cos_k1 * (fresnel_c_k2 - fresnel_c_k3) + d * sin_k1 * (fresnel_s_k2 - fresnel_s_k3)) /
                    abs_sigma +
                d * (sin_k4 - sin_k5) / control.sigma;
    F_u(2, 1) = d * integration_step;
  }
  else
  {
    if (fabs(state.kappa) > get_epsilon())
    {
      double sin_th = sin(state.theta);
      double cos_th = cos(state.theta);
      double k9 = state.theta + d * integration_step * state.kappa;
      double cos_k9 = cos(k9);
      double sin_k9 = sin(k9);

      // df/dx
      F_x(0, 0) = 1.0;
      F_x(1, 1) = 1.0;

      F_x(0, 2) = (-cos_th + cos_k9) / state.kappa;
      F_x(1, 2) = (-sin_th + sin_k9) / state.kappa;
      F_x(2, 2) = 1.0;

      // df/du
      F_u(0, 0) = cos_k9;
      F_u(1, 0) = sin_k9;
      F_u(2, 0) = state.kappa;

      F_u(0, 1) = (sin_th - sin_k9) / (state.kappa * state.kappa) + d * integration_step * cos_k9 / state.kappa;
      F_u(1, 1) = (-cos_th + cos_k9) / (state.kappa * state.kappa) + d * integration_step * sin_k9 / state.kappa;
      F_u(2, 1) = d * integration_step;
    }
    else
    {
      double cos_th = cos(state.theta);
      double sin_th = sin(state.theta);

      // df/dx
      F_x(0, 0) = 1.0;
      F_x(1, 1) = 1.0;

      F_x(0, 2) = -d * integration_step * sin_th;
      F_x(1, 2) = d * integration_step * cos_th;
      F_x(2, 2) = 1.0;

      // df/du
      F_u(0, 0) = cos_th;
      F_u(1, 0) = sin_th;
      F_u(2, 0) = state.kappa;

      F_u(2, 1) = d * integration_step;
    }
  }
}

Matrix3d EKF::get_observation_jacobi() const
{
  return I_;
}

Matrix2d EKF::get_motion_covariance(const State &state, const Control &control, double integration_step) const
{
  Matrix2d Q(Matrix2d::Zero());
  Q(0, 0) =
      motion_noise_.alpha1 * integration_step * integration_step + motion_noise_.alpha2 * state.kappa * state.kappa;
  Q(1, 1) =
      motion_noise_.alpha3 * integration_step * integration_step + motion_noise_.alpha4 * state.kappa * state.kappa;
  return Q;
}

Matrix3d EKF::get_observation_covariance() const
{
  Matrix3d R(Matrix3d::Zero());
  R(0, 0) = measurement_noise_.std_x * measurement_noise_.std_x;
  R(1, 1) = measurement_noise_.std_y * measurement_noise_.std_y;
  R(2, 2) = measurement_noise_.std_theta * measurement_noise_.std_theta;
  return R;
}

Matrix23d EKF::get_controller_gain(const Control &control) const
{
  Matrix23d K(Matrix23d::Zero());
  K(0, 0) = controller_.k1;
  K(1, 1) = controller_.k2;
  K(1, 2) = sgn(control.delta_s) * controller_.k3;
  return K;
}

Matrix3d EKF::get_rotation_matrix(double angle) const
{
  Matrix3d R(Matrix3d::Zero());
  double cos_angle = cos(angle);
  double sin_angle = sin(angle);
  R(0, 0) = cos_angle;
  R(1, 0) = -sin_angle;
  R(0, 1) = sin_angle;
  R(1, 1) = cos_angle;
  R(2, 2) = 1.0;
  return R;
}

void EKF::predict(const State_With_Covariance &state, const Control &control, double integration_step,
                  State_With_Covariance &state_pred) const
{
  Matrix3d Sigma = covariance_to_eigen(state.Sigma);
  Matrix3d Lambda = covariance_to_eigen(state.Lambda);

  // Sigma_pred
  Matrix3d F_x(Matrix3d::Zero());
  Matrix32d F_u(Matrix32d::Zero());
  get_motion_jacobi(state.state, control, integration_step, F_x, F_u);
  Matrix2d Q = get_motion_covariance(state.state, control, integration_step);
  Matrix3d Sigma_pred = F_x * Sigma * F_x.transpose() + F_u * Q * F_u.transpose();

  // Lambda_pred
  Matrix23d K = get_controller_gain(control);
  Matrix3d R_1I = get_rotation_matrix(state.state.theta);
  Matrix3d F_K = F_x - F_u * K * R_1I;
  Matrix3d Lambda_pred = F_K * Lambda * F_K.transpose();

  eigen_to_covariance(Sigma_pred, state_pred.Sigma);
  eigen_to_covariance(Lambda_pred, state_pred.Lambda);
  eigen_to_covariance(Sigma_pred + Lambda_pred, state_pred.covariance);
}

void EKF::update(const State_With_Covariance &state_pred, State_With_Covariance &state_corr) const
{
  Matrix3d Sigma_pred = covariance_to_eigen(state_pred.Sigma);
  Matrix3d Lambda_pred = covariance_to_eigen(state_pred.Lambda);

  // Kalman gain L
  Matrix3d H_x = get_observation_jacobi();
  Matrix3d R = get_observation_covariance();
  Matrix3d Sigma_xz = Sigma_pred * H_x.transpose();
  Matrix3d Sigma_z = H_x * Sigma_xz + R;
  Matrix3d L = Sigma_xz * Sigma_z.inverse();

  // Sigma_corr (Joseph's form)
  Matrix3d I_LH_x = I_ - L * H_x;
  Matrix3d Sigma_corr = I_LH_x * Sigma_pred * I_LH_x.transpose() + L * R * L.transpose();

  // Lambda_corr
  Matrix3d Lambda_corr = Lambda_pred + L * Sigma_xz.transpose();

  eigen_to_covariance(Sigma_corr, state_corr.Sigma);
  eigen_to_covariance(Lambda_corr, state_corr.Lambda);
  eigen_to_covariance(Sigma_corr + Lambda_corr, state_corr.covariance);
}

} // namespace steering
