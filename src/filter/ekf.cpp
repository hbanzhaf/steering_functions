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

#include "steering_functions/filter/ekf.hpp"

void EKF::set_parameters(const Motion_Noise &motion_noise, const Measurement_Noise &measurement_noise,
                         const Controller &controller)
{
  motion_noise_ = motion_noise;
  measurement_noise_ = measurement_noise;
  controller_ = controller;
}

Matrix4d EKF::covariance_to_eigen(const double covariance[16]) const
{
  Matrix4d covariance_eigen;
  for (int i = 0; i < covariance_eigen.rows(); ++i)
  {
    for (int j = 0; j < covariance_eigen.cols(); ++j)
    {
      covariance_eigen(i, j) = covariance[i + 4 * j];
    }
  }
  return covariance_eigen;
}

void EKF::eigen_to_covariance(const Matrix4d &covariance_eigen, double covariance[16]) const
{
  for (int i = 0; i < covariance_eigen.rows(); ++i)
  {
    for (int j = 0; j < covariance_eigen.cols(); ++j)
    {
      covariance[i + 4 * j] = covariance_eigen(i, j);
    }
  }
}

Matrix4d EKF::get_motion_jacobi_x(const State &state, const Control &control, double integration_step) const
{
  Matrix4d F_x(Matrix4d::Zero());
  double d(sgn(control.delta_s));
  double sin_th = sin(state.theta);
  double cos_th = cos(state.theta);

  if (fabs(state.kappa) > get_epsilon())
  {
    F_x(0, 0) = 1.0;
    F_x(1, 1) = 1.0;

    F_x(0, 2) = (1 / state.kappa) * (-cos_th + cos(state.theta + d * integration_step * state.kappa));
    F_x(1, 2) = (1 / state.kappa) * (-sin_th + sin(state.theta + d * integration_step * state.kappa));
    F_x(2, 2) = 1.0;

    F_x(0, 3) = -(-sin_th + sin(state.theta + d * integration_step * state.kappa)) / (state.kappa * state.kappa) +
                (d * integration_step / state.kappa) * cos(state.theta + d * integration_step * state.kappa);
    F_x(1, 3) = -(cos_th - cos(state.theta + d * integration_step * state.kappa)) / (state.kappa * state.kappa) +
                (d * integration_step / state.kappa) * sin(state.theta + d * integration_step * state.kappa);
    F_x(2, 3) = d * integration_step;
    F_x(3, 3) = 1.0;
  }
  else
  {
    F_x(0, 0) = 1.0;
    F_x(1, 1) = 1.0;

    F_x(0, 2) = -d * integration_step * sin_th;
    F_x(1, 2) = d * integration_step * cos_th;
    F_x(2, 2) = 1.0;

    F_x(2, 3) = d * integration_step;
    F_x(3, 3) = 1.0;
  }
  return F_x;
}

Matrix42d EKF::get_motion_jacobi_u(const State &state, const Control &control, double integration_step) const
{
  Matrix42d F_u(Matrix42d::Zero());
  double d(sgn(control.delta_s));
  if (fabs(state.kappa) > get_epsilon())
  {
    F_u(0, 0) = cos(state.theta + d * integration_step * state.kappa);
    F_u(1, 0) = sin(state.theta + d * integration_step * state.kappa);
    F_u(2, 0) = state.kappa + control.sigma * integration_step;
    F_u(3, 0) = d * control.sigma;

    F_u(2, 1) = 0.5 * d * integration_step * integration_step;
    F_u(3, 1) = integration_step;
  }
  else
  {
    F_u(0, 0) = cos(state.theta);
    F_u(1, 0) = sin(state.theta);
    F_u(2, 0) = state.kappa + control.sigma * integration_step;
    F_u(3, 0) = d * control.sigma;

    F_u(2, 1) = 0.5 * d * integration_step * integration_step;
    F_u(3, 1) = integration_step;
  }
  return F_u;
}

Matrix4d EKF::get_observation_jacobi_x() const
{
  return Matrix4d::Identity();
}

Matrix2d EKF::get_motion_covariance(const State &state, const Control &control, double integration_step) const
{
  Matrix2d Q(Matrix2d::Zero());
  Q(0, 0) = motion_noise_.alpha1 * integration_step * integration_step +
            motion_noise_.alpha2 * state.kappa * state.kappa + motion_noise_.alpha3 * control.sigma * control.sigma;
  Q(1, 1) = motion_noise_.alpha4 * integration_step * integration_step +
            motion_noise_.alpha5 * state.kappa * state.kappa + motion_noise_.alpha6 * control.sigma * control.sigma;
  return Q;
}

Matrix4d EKF::get_observation_covariance() const
{
  Matrix4d R(Matrix4d::Zero());
  R(0, 0) = measurement_noise_.std_x * measurement_noise_.std_x;
  R(1, 1) = measurement_noise_.std_y * measurement_noise_.std_y;
  R(2, 2) = measurement_noise_.std_theta * measurement_noise_.std_theta;
  R(3, 3) = measurement_noise_.std_kappa * measurement_noise_.std_kappa;
  return R;
}

Matrix24d EKF::get_controller_gain(const Control &control) const
{
  Matrix24d K(Matrix24d::Zero());
  K(0, 0) = controller_.k1;
  K(1, 1) = controller_.k2;
  K(1, 2) = sgn(control.delta_s) * controller_.k3;
  K(1, 3) = controller_.k4;
  return K;
}

Matrix4d EKF::get_rotation_matrix(double angle) const
{
  Matrix4d R(Matrix4d::Zero());
  double cos_angle = cos(angle);
  double sin_angle = sin(angle);
  R(0, 0) = cos_angle;
  R(1, 0) = -sin_angle;
  R(0, 1) = sin_angle;
  R(1, 1) = cos_angle;
  R(2, 2) = 1.0;
  R(3, 3) = 1.0;
  return R;
}

void EKF::predict(const State_With_Covariance &state, const Control &control, double integration_step,
                  State_With_Covariance &state_pred) const
{
  Matrix4d Sigma = covariance_to_eigen(state.Sigma);
  Matrix4d Lambda = covariance_to_eigen(state.Lambda);

  // Sigma_pred
  Matrix4d F_x = get_motion_jacobi_x(state.state, control, integration_step);
  Matrix42d F_u = get_motion_jacobi_u(state.state, control, integration_step);
  Matrix2d Q = get_motion_covariance(state.state, control, integration_step);
  Matrix4d Sigma_pred = F_x * Sigma * F_x.transpose() + F_u * Q * F_u.transpose();

  // Lambda_pred
  Matrix24d K = get_controller_gain(control);
  Matrix4d R_1I = get_rotation_matrix(state.state.theta);
  Matrix4d F_K = F_x - F_u * K * R_1I;
  Matrix4d Lambda_pred = F_K * Lambda * F_K.transpose();

  eigen_to_covariance(Sigma_pred, state_pred.Sigma);
  eigen_to_covariance(Lambda_pred, state_pred.Lambda);
  eigen_to_covariance(Sigma_pred + Lambda_pred, state_pred.covariance);
}

void EKF::update(const State_With_Covariance &state_pred, State_With_Covariance &state_corr) const
{
  Matrix4d Sigma_pred = covariance_to_eigen(state_pred.Sigma);
  Matrix4d Lambda_pred = covariance_to_eigen(state_pred.Lambda);

  // Kalman gain L
  Matrix4d H_x = get_observation_jacobi_x();
  Matrix4d R = get_observation_covariance();
  Matrix4d Sigma_xz = Sigma_pred * H_x.transpose();
  Matrix4d Sigma_z = H_x * Sigma_xz + R;
  Matrix4d L = Sigma_xz * Sigma_z.inverse();

  // Sigma_corr (Joseph's form)
  Matrix4d I = Matrix4d::Identity();
  Matrix4d LH_x = L * H_x;
  Matrix4d Sigma_corr = (I - LH_x) * Sigma_pred * (I - LH_x).transpose() + L * R * L.transpose();

  // Lambda_corr
  Matrix4d Lambda_corr = Lambda_pred + L * Sigma_xz.transpose();

  eigen_to_covariance(Sigma_corr, state_corr.Sigma);
  eigen_to_covariance(Lambda_corr, state_corr.Lambda);
  eigen_to_covariance(Sigma_corr + Lambda_corr, state_corr.covariance);
}
