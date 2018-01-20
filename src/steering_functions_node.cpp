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

#include <iostream>

#include <costmap_2d/footprint.h>
#include <geometry_msgs/PoseArray.h>
#include <nav_msgs/Path.h>
#include <ros/ros.h>
#include <tf/transform_datatypes.h>
#include <visualization_msgs/MarkerArray.h>

#include <Eigen/Dense>

#include "steering_functions/dubins_state_space/dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/cc_dubins_state_space.hpp"
#include "steering_functions/hc_cc_state_space/cc_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hc00_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hc0pm_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hcpm0_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_cc_state_space/hcpmpm_reeds_shepp_state_space.hpp"
#include "steering_functions/reeds_shepp_state_space/reeds_shepp_state_space.hpp"
#include "steering_functions/steering_functions.hpp"

#define FRAME_ID "/world"
#define DISCRETIZATION 0.1               // [m]
#define VISUALIZATION_DURATION 2         // [s]
#define ANIMATE false                    // [-]
#define OPERATING_REGION_X 20.0          // [m]
#define OPERATING_REGION_Y 20.0          // [m]
#define OPERATING_REGION_THETA 2 * M_PI  // [rad]
#define random(lower, upper) (rand() * (upper - lower) / RAND_MAX + lower)

using namespace std;

class PathClass
{
public:
  // publisher
  ros::Publisher pub_path_;
  ros::Publisher pub_poses_;
  ros::Publisher pub_text_;
  ros::Publisher pub_covariances_;

  // path properties
  string id_;
  string path_type_;
  double discretization_;
  State_With_Covariance state_start_;
  State state_goal_;
  double kappa_max_;
  double sigma_max_;
  vector<State_With_Covariance> path_;

  // visualization
  string frame_id_;
  nav_msgs::Path nav_path_;
  geometry_msgs::PoseArray pose_array_;
  visualization_msgs::MarkerArray marker_array_text_;
  visualization_msgs::MarkerArray marker_array_covariance_;
  visualization_msgs::Marker marker_text_;
  visualization_msgs::Marker marker_covariance_;

  // constructor
  PathClass(const string& path_type, const State_With_Covariance& state_start, const State& state_goal,
            const double kappa_max, const double sigma_max)
    : path_type_(path_type)
    , discretization_(DISCRETIZATION)
    , state_start_(state_start)
    , state_goal_(state_goal)
    , kappa_max_(kappa_max)
    , sigma_max_(sigma_max)
    , frame_id_(FRAME_ID)
  {
    ros::NodeHandle pnh("~");

    // publisher
    pub_path_ = pnh.advertise<nav_msgs::Path>("visualization_path", 10);
    pub_poses_ = pnh.advertise<geometry_msgs::PoseArray>("visualization_poses", 10);
    pub_text_ = pnh.advertise<visualization_msgs::MarkerArray>("visualization_text", 10);
    pub_covariances_ = pnh.advertise<visualization_msgs::MarkerArray>("visualization_covariances", 10);
    while (pub_path_.getNumSubscribers() == 0 || pub_poses_.getNumSubscribers() == 0 ||
           pub_text_.getNumSubscribers() == 0 || pub_covariances_.getNumSubscribers() == 0)
      ros::Duration(0.001).sleep();

    // path
    //    if (path_type_ == "CC_Dubins")
    //    {
    //      id_ = "1";
    //      CC_Dubins_State_Space state_space(kappa_max_, sigma_max_, discretization_, true);
    //      path_ = state_space.get_path(state_start_, state_goal_);
    //    }
    //    else if (path_type_ == "Dubins")
    //    {
    //      id_ = "2";
    //      Dubins_State_Space state_space(kappa_max_, discretization_, true);
    //      path_ = state_space.get_path(state_start_, state_goal_);
    //    }
    //    else if (path_type_ == "CC_RS")
    //    {
    //      id_ = "3";
    //      CC_Reeds_Shepp_State_Space state_space(kappa_max_, sigma_max_, discretization_);
    //      path_ = state_space.get_path(state_start_, state_goal_);
    //    }
    //    else if (path_type_ == "HC00")
    //    {
    //      id_ = "4";
    //      HC00_Reeds_Shepp_State_Space state_space(kappa_max_, sigma_max_, discretization_);
    //      path_ = state_space.get_path(state_start_, state_goal_);
    //    }
    //    else if (path_type_ == "HC0pm")
    //    {
    //      id_ = "5";
    //      HC0pm_Reeds_Shepp_State_Space state_space(kappa_max_, sigma_max_, discretization_);
    //      path_ = state_space.get_path(state_start_, state_goal_);
    //    }
    //    else if (path_type_ == "HCpm0")
    //    {
    //      id_ = "6";
    //      HCpm0_Reeds_Shepp_State_Space state_space(kappa_max_, sigma_max_, discretization_);
    //      path_ = state_space.get_path(state_start_, state_goal_);
    //    }
    //    else if (path_type_ == "HCpmpm")
    //    {
    //      id_ = "7";life
    //      HCpmpm_Reeds_Shepp_State_Space state_space(kappa_max_, sigma_max_, discretization_);
    //      path_ = state_space.get_path(state_start_, state_goal_);
    //    }
    if (path_type_ == "RS")
    {
      id_ = "8";
      Reeds_Shepp_State_Space state_space(kappa_max_, discretization_);
      path_ = state_space.get_path_with_covariance(state_start_, state_goal_);
    }

    // nav_path
    nav_path_.header.frame_id = frame_id_;

    // pose_array
    pose_array_.header.frame_id = frame_id_;

    // marker_text
    marker_text_.header.frame_id = frame_id_;
    marker_text_.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
    marker_text_.action = visualization_msgs::Marker::ADD;
    marker_text_.color.r = 1.0;
    marker_text_.color.g = 1.0;
    marker_text_.color.b = 1.0;
    marker_text_.color.a = 1.0;

    // marker_covariance
    marker_covariance_.header.frame_id = frame_id_;
    marker_covariance_.type = visualization_msgs::Marker::SPHERE;
    marker_covariance_.action = visualization_msgs::Marker::ADD;
    marker_covariance_.lifetime = ros::Duration(VISUALIZATION_DURATION);
    marker_covariance_.color.r = 0.6;
    marker_covariance_.color.g = 0.0;
    marker_covariance_.color.b = 0.6;
    marker_covariance_.color.a = 0.4;
  }

  // visualization
  void covariance_to_marker(const State_With_Covariance& state, visualization_msgs::Marker& marker)
  {
    marker.pose.position.x = state.state.x;
    marker.pose.position.y = state.state.y;
    marker.pose.position.z = 0.0;

    // eigenvalues/-vectors
    Eigen::Vector3d eigenvalues(Eigen::Vector3d::Identity());
    Eigen::Matrix3d eigenvectors(Eigen::Matrix3d::Zero());
    Eigen::Matrix3d covariance;
    covariance << state.covariance[0 + 4 * 0], state.covariance[1 + 4 * 0], state.covariance[2 + 4 * 0],
        state.covariance[0 + 4 * 1], state.covariance[1 + 4 * 1], state.covariance[2 + 4 * 1],
        state.covariance[0 + 4 * 2], state.covariance[1 + 4 * 2], state.covariance[2 + 4 * 2];
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covariance);
    if (eigensolver.info() == Eigen::Success)
    {
      eigenvalues = eigensolver.eigenvalues();
      eigenvectors = eigensolver.eigenvectors();
    }
    else
    {
      ROS_WARN("Failed to visualize covariance because eigenvalues can not be computed.");
    }

    // make right-handed system
    if (eigenvectors.determinant() < 0)
      eigenvectors.col(0) *= -1.0;

    // orientation
    tf::Matrix3x3 rotation(eigenvectors(0, 0), eigenvectors(0, 1), eigenvectors(0, 2), eigenvectors(1, 0),
                           eigenvectors(1, 1), eigenvectors(1, 2), eigenvectors(2, 0), eigenvectors(2, 1),
                           eigenvectors(2, 2));
    tf::Quaternion quaternion;
    rotation.getRotation(quaternion);
    quaternionTFToMsg(quaternion, marker.pose.orientation);

    // scale with 3 * standard deviation
    marker.scale.x = 2 * 3 * sqrt(eigenvalues[0]);
    marker.scale.y = 2 * 3 * sqrt(eigenvalues[1]);
    marker.scale.z = 2 * 3 * sqrt(eigenvalues[2]);
  }

  void visualize()
  {
    // start
    marker_text_.id = 1;
    marker_text_.pose.position.x = state_start_.state.x;
    marker_text_.pose.position.y = state_start_.state.y;
    marker_text_.scale.z = 0.7;
    marker_text_.text = "start";
    marker_array_text_.markers.push_back(marker_text_);

    geometry_msgs::Pose pose_start;
    pose_start.position.x = state_start_.state.x;
    pose_start.position.y = state_start_.state.y;
    pose_start.orientation.z = sin(state_start_.state.theta / 2.0);
    pose_start.orientation.w = cos(state_start_.state.theta / 2.0);
    pose_array_.poses.push_back(pose_start);

    // goal
    marker_text_.id = 2;
    marker_text_.pose.position.x = state_goal_.x;
    marker_text_.pose.position.y = state_goal_.y;
    marker_text_.scale.z = 0.7;
    marker_text_.text = "goal";
    marker_array_text_.markers.push_back(marker_text_);

    geometry_msgs::Pose pose_goal;
    pose_goal.position.x = state_goal_.x;
    pose_goal.position.y = state_goal_.y;
    pose_goal.orientation.z = sin(state_goal_.theta / 2.0);
    pose_goal.orientation.w = cos(state_goal_.theta / 2.0);
    pose_array_.poses.push_back(pose_goal);

    // path description
    marker_text_.id = 3;
    marker_text_.pose.position.x = 0.0;
    marker_text_.pose.position.y = 12.0;
    marker_text_.scale.z = 1.0;
    marker_text_.text = id_ + ") " + path_type_ + " Steer";
    marker_array_text_.markers.push_back(marker_text_);

    // path
    for (const auto& state : path_)
    {
      geometry_msgs::PoseStamped pose;
      pose.pose.position.x = state.state.x;
      pose.pose.position.y = state.state.y;
      pose.pose.orientation.z = sin(state.state.theta / 2.0);
      pose.pose.orientation.w = cos(state.state.theta / 2.0);
      nav_path_.poses.push_back(pose);
    }

    // covariances
    for (int i = 0; i < path_.size(); i += 10)
    {
      marker_covariance_.id = i;
      covariance_to_marker(path_[i], marker_covariance_);
      marker_array_covariance_.markers.push_back(marker_covariance_);
    }

    // publish
    pub_path_.publish(nav_path_);
    pub_poses_.publish(pose_array_);
    pub_text_.publish(marker_array_text_);
    pub_covariances_.publish(marker_array_covariance_);
    ros::spinOnce();
  }
};

class RobotClass
{
public:
  // publisher
  ros::Publisher pub_swath_;

  // robot config
  double kappa_max_;
  double sigma_max_;
  double wheel_base_;
  double track_width_;
  geometry_msgs::Point wheel_fl_pos_, wheel_fr_pos_;
  vector<geometry_msgs::Point> footprint_;
  vector<geometry_msgs::Point> wheel_;
  double wheel_radius_;
  double wheel_width_;

  // visualization
  string frame_id_;
  bool animate_;
  visualization_msgs::MarkerArray marker_array_swath_;
  visualization_msgs::Marker marker_chassis_, marker_wheels_;

  // constructor
  explicit RobotClass() : frame_id_(FRAME_ID), animate_(ANIMATE)
  {
    ros::NodeHandle nh;
    ros::NodeHandle pnh("~");

    // publisher
    pub_swath_ = pnh.advertise<visualization_msgs::MarkerArray>("visualization_swath", 10);
    while (pub_swath_.getNumSubscribers() == 0)
      ros::Duration(0.001).sleep();

    // robot kinematic
    nh.param<double>("kappa_max", kappa_max_, 1.0);
    nh.param<double>("sigma_max", sigma_max_, 1.0);

    // footprint
    footprint_ = costmap_2d::makeFootprintFromParams(nh);

    // wheel
    nh.getParam("wheel_base", wheel_base_);
    nh.getParam("track_width", track_width_);
    wheel_fl_pos_.x = wheel_base_;
    wheel_fl_pos_.y = track_width_ / 2.0;
    wheel_fr_pos_.x = wheel_base_;
    wheel_fr_pos_.y = -track_width_ / 2.0;

    nh.getParam("wheel_radius", wheel_radius_);
    nh.getParam("wheel_width", wheel_width_);
    geometry_msgs::Point point1, point2, point3, point4;
    point1.x = wheel_radius_;
    point1.y = wheel_width_ / 2.0;
    point2.x = -wheel_radius_;
    point2.y = wheel_width_ / 2.0;
    point3.x = -wheel_radius_;
    point3.y = -wheel_width_ / 2.0;
    point4.x = wheel_radius_;
    point4.y = -wheel_width_ / 2.0;
    wheel_.push_back(point1);
    wheel_.push_back(point2);
    wheel_.push_back(point3);
    wheel_.push_back(point4);

    // marker chassis
    marker_chassis_.header.frame_id = frame_id_;
    marker_chassis_.action = visualization_msgs::Marker::ADD;
    marker_chassis_.id = 1;
    marker_chassis_.type = visualization_msgs::Marker::LINE_LIST;
    marker_chassis_.scale.x = 0.03;
    marker_chassis_.color.r = 0.6;
    marker_chassis_.color.g = 0.6;
    marker_chassis_.color.b = 0.6;
    marker_chassis_.color.a = 0.5;

    // marker wheels
    marker_wheels_.header.frame_id = frame_id_;
    marker_wheels_.action = visualization_msgs::Marker::ADD;
    marker_wheels_.id = 2;
    marker_wheels_.type = visualization_msgs::Marker::LINE_LIST;
    marker_wheels_.scale.x = 0.03;
    marker_wheels_.color.r = 0.9;
    marker_wheels_.color.g = 0.9;
    marker_wheels_.color.b = 0.9;
    marker_wheels_.color.a = 1.0;
  }

  // translate and rotate polygon
  void get_oriented_polygon(double x, double y, double theta, const vector<geometry_msgs::Point>& polygon,
                            vector<geometry_msgs::Point>& oriented_polygon) const
  {
    oriented_polygon.clear();
    double cos_th = cos(theta);
    double sin_th = sin(theta);
    geometry_msgs::Point new_point;
    for (const auto& point : polygon)
    {
      new_point.x = x + (point.x * cos_th - point.y * sin_th);
      new_point.y = y + (point.x * sin_th + point.y * cos_th);
      oriented_polygon.push_back(new_point);
    }
  }

  // visualization
  void polygon_to_marker(const vector<geometry_msgs::Point>& polygon, visualization_msgs::Marker& marker)
  {
    auto point1 = polygon.back();
    for (auto const& point2 : polygon)
    {
      marker.points.push_back(point1);
      marker.points.push_back(point2);
      point1 = point2;
    }
  }

  void visualize(const vector<State_With_Covariance>& path)
  {
    marker_array_swath_.markers.clear();
    marker_chassis_.points.clear();
    marker_wheels_.points.clear();
    for (const auto& state : path)
    {
      double steer_angle_fl, steer_angle_fr;
      vector<geometry_msgs::Point> wheel_fl, wheel_fr;
      vector<geometry_msgs::Point> oriented_wheel_fl, oriented_wheel_fr;
      vector<geometry_msgs::Point> oriented_footprint;

      // steering angle
      if (fabs(state.state.kappa) > 1e-4)
      {
        steer_angle_fl = atan(wheel_fl_pos_.x / ((1 / state.state.kappa) - wheel_fl_pos_.y));
        steer_angle_fr = atan(wheel_fr_pos_.x / ((1 / state.state.kappa) - wheel_fr_pos_.y));
      }
      else
      {
        steer_angle_fl = 0.0;
        steer_angle_fr = 0.0;
      }

      // transform wheels and footprint
      get_oriented_polygon(wheel_fl_pos_.x, wheel_fl_pos_.y, steer_angle_fl, wheel_, wheel_fl);
      get_oriented_polygon(wheel_fr_pos_.x, wheel_fr_pos_.y, steer_angle_fr, wheel_, wheel_fr);
      get_oriented_polygon(state.state.x, state.state.y, state.state.theta, wheel_fl, oriented_wheel_fl);
      get_oriented_polygon(state.state.x, state.state.y, state.state.theta, wheel_fr, oriented_wheel_fr);
      get_oriented_polygon(state.state.x, state.state.y, state.state.theta, footprint_, oriented_footprint);

      polygon_to_marker(oriented_footprint, marker_chassis_);
      polygon_to_marker(oriented_wheel_fl, marker_wheels_);
      polygon_to_marker(oriented_wheel_fr, marker_wheels_);

      // animate
      if (animate_)
      {
        marker_array_swath_.markers.clear();
        marker_array_swath_.markers.push_back(marker_chassis_);
        marker_array_swath_.markers.push_back(marker_wheels_);
        pub_swath_.publish(marker_array_swath_);
        ros::spinOnce();
        ros::Duration(0.08).sleep();
      }
    }

    // publish
    if (!animate_)
    {
      marker_array_swath_.markers.push_back(marker_chassis_);
      marker_array_swath_.markers.push_back(marker_wheels_);
      pub_swath_.publish(marker_array_swath_);
      ros::spinOnce();
    }
  }
};

int main(int argc, char** argv)
{
  ros::init(argc, argv, "steering_functions");

  RobotClass robot;

  // sample
  int seed(ros::Time::now().toSec());
  srand(seed);
  while (ros::ok())
  {
    State_With_Covariance start;
    start.state.x = random(-OPERATING_REGION_X / 2.0, OPERATING_REGION_X / 2.0);
    start.state.y = random(-OPERATING_REGION_Y / 2.0, OPERATING_REGION_Y / 2.0);
    start.state.theta = random(-OPERATING_REGION_THETA / 2.0, OPERATING_REGION_THETA / 2.0);
    start.state.kappa = 0.0;
    start.state.d = 0.0;
    start.covariance[0 + 4 * 0] = 0.1 * 0.1;
    start.covariance[1 + 4 * 1] = 0.1 * 0.1;
    start.covariance[2 + 4 * 2] = 0.01 * 0.01;

    State goal;
    goal.x = random(-OPERATING_REGION_X / 2.0, OPERATING_REGION_X / 2.0);
    goal.y = random(-OPERATING_REGION_Y / 2.0, OPERATING_REGION_Y / 2.0);
    goal.theta = random(-OPERATING_REGION_THETA / 2.0, OPERATING_REGION_THETA / 2.0);
    goal.kappa = 0.0;
    goal.d = 0.0;

    //    PathClass cc_dubins_path("CC_Dubins", start, goal, robot.kappa_max_, robot.sigma_max_);
    //    PathClass dubins_path("Dubins", start, goal, robot.kappa_max_, robot.sigma_max_);
    //    PathClass cc_rs_path("CC_RS", start, goal, robot.kappa_max_, robot.sigma_max_);
    //    PathClass hc00_path("HC00", start, goal, robot.kappa_max_, robot.sigma_max_);
    //    PathClass hc0pm_path("HC0pm", start, goal, robot.kappa_max_, robot.sigma_max_);
    //    PathClass hcpm0_path("HCpm0", start, goal, robot.kappa_max_, robot.sigma_max_);
    //    PathClass hcpmpm_path("HCpmpm", start, goal, robot.kappa_max_, robot.sigma_max_);
    PathClass rs_path("RS", start, goal, robot.kappa_max_, robot.sigma_max_);

    // visualize
    //    cc_dubins_path.visualize();
    //    robot.visualize(cc_dubins_path.path_);
    //    ros::Duration(VISUALIZATION_DURATION).sleep();

    //    dubins_path.visualize();
    //    robot.visualize(dubins_path.path_);
    //    ros::Duration(VISUALIZATION_DURATION).sleep();

    //    cc_rs_path.visualize();
    //    robot.visualize(cc_rs_path.path_);
    //    ros::Duration(VISUALIZATION_DURATION).sleep();

    //    hc00_path.visualize();
    //    robot.visualize(hc00_path.path_);
    //    ros::Duration(VISUALIZATION_DURATION).sleep();

    //    hc0pm_path.visualize();
    //    robot.visualize(hc0pm_path.path_);
    //    ros::Duration(VISUALIZATION_DURATION).sleep();

    //    hcpm0_path.visualize();
    //    robot.visualize(hcpm0_path.path_);
    //    ros::Duration(VISUALIZATION_DURATION).sleep();

    //    hcpmpm_path.visualize();
    //    robot.visualize(hcpmpm_path.path_);
    //    ros::Duration(VISUALIZATION_DURATION).sleep();

    rs_path.visualize();
    robot.visualize(rs_path.path_);
    ros::Duration(VISUALIZATION_DURATION).sleep();
  }
  return 0;
}
