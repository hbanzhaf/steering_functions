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

#include "steering_functions/hc_cc_state_space/hcpm0_reeds_shepp_state_space.hpp"

#define HC_REGULAR false
#define CC_REGULAR false

class HCpm0_Reeds_Shepp_State_Space::HCpm0_Reeds_Shepp
{
private:
  HCpm0_Reeds_Shepp_State_Space *parent_;

public:
  explicit HCpm0_Reeds_Shepp(HCpm0_Reeds_Shepp_State_Space *parent)
  {
    parent_ = parent;
  }

  double distance = 0.0;
  double angle = 0.0;

  // ##### TT ###################################################################
  bool TT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return fabs(distance - 2 * c2.radius) < get_epsilon();
  }

  void TT_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q) const
  {
    double x = (c1.xc + c2.xc) / 2;
    double y = (c1.yc + c2.yc) / 2;
    double angle = atan2(c2.yc - c1.yc, c2.xc - c1.xc);
    double theta;
    if (c1.left)
    {
      if (c1.forward)
      {
        theta = angle + HALF_PI - c2.mu;
      }
      else
      {
        theta = angle + HALF_PI + c2.mu;
      }
    }
    else
    {
      if (c1.forward)
      {
        theta = angle - HALF_PI + c2.mu;
      }
      else
      {
        theta = angle - HALF_PI - c2.mu;
      }
    }
    *q = new Configuration(x, y, theta, 0);
  }

  double TT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                 Configuration **q1, Configuration **q2) const
  {
    TT_tangent_circles(c1, c2, q2);
    *cstart = new HC_CC_Circle(**q2, !c2.left, c2.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    return (*cstart)->hc_turn_length(**q1) + (*cend)->cc_turn_length(**q2);
  }

  // ##### TcT ##################################################################
  bool TcT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return fabs(distance - fabs(2 / c1.kappa)) < get_epsilon();
  }

  void TcT_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q) const
  {
    double distance = center_distance(c1, c2);
    double delta_x = 0.5 * distance;
    double delta_y = 0.0;
    double angle = atan2(c2.yc - c1.yc, c2.xc - c1.xc);
    double x, y, theta;
    if (c1.left)
    {
      if (c1.forward)
      {
        theta = angle + HALF_PI;
        global_frame_change(c1.xc, c1.yc, angle, delta_x, delta_y, &x, &y);
      }
      else
      {
        theta = angle + HALF_PI;
        global_frame_change(c1.xc, c1.yc, angle, delta_x, -delta_y, &x, &y);
      }
    }
    else
    {
      if (c1.forward)
      {
        theta = angle - HALF_PI;
        global_frame_change(c1.xc, c1.yc, angle, delta_x, -delta_y, &x, &y);
      }
      else
      {
        theta = angle - HALF_PI;
        global_frame_change(c1.xc, c1.yc, angle, delta_x, delta_y, &x, &y);
      }
    }
    *q = new Configuration(x, y, theta, c1.kappa);
  }

  double TcT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                  Configuration **q) const
  {
    TcT_tangent_circles(c1, c2, q);
    *cstart = new HC_CC_Circle(c1);
    *cend = new HC_CC_Circle(c2);
    return (*cstart)->rs_turn_length(**q) + (*cend)->hc_turn_length(**q);
  }

  // ##### Reeds-Shepp families: ################################################

  // ##### TcTcT ################################################################
  bool TcTcT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return distance <= fabs(4 / c1.kappa);
  }

  void TcTcT_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2,
                             Configuration **q3, Configuration **q4) const
  {
    double theta = angle;
    double r = fabs(2 / c1.kappa);
    double delta_x = 0.5 * distance;
    double delta_y = sqrt(fabs(pow(r, 2) - pow(delta_x, 2)));
    double x, y;

    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c2.left, c2.forward, c2.regular, parent_->rs_circle_param_);
    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    HC_CC_Circle tgt2(x, y, !c2.left, c2.forward, c2.regular, parent_->rs_circle_param_);

    TcT_tangent_circles(c1, tgt1, q1);
    TcT_tangent_circles(tgt1, c2, q2);
    TcT_tangent_circles(c1, tgt2, q3);
    TcT_tangent_circles(tgt2, c2, q4);
  }

  double TcTcT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                    Configuration **q1, Configuration **q2, HC_CC_Circle **ci) const
  {
    Configuration *qa, *qb, *qc, *qd;
    TcTcT_tangent_circles(c1, c2, &qa, &qb, &qc, &qd);
    HC_CC_Circle *middle1, *middle2;
    middle1 = new HC_CC_Circle(*qa, !c2.left, c2.forward, true, parent_->rs_circle_param_);
    middle2 = new HC_CC_Circle(*qc, !c2.left, c2.forward, true, parent_->rs_circle_param_);

    *cstart = new HC_CC_Circle(c1);
    *cend = new HC_CC_Circle(c2);

    // select shortest connection
    double length1 = (*cstart)->rs_turn_length(*qa) + middle1->rs_turn_length(*qb) + (*cend)->hc_turn_length(*qb);
    double length2 = (*cstart)->rs_turn_length(*qc) + middle2->rs_turn_length(*qd) + (*cend)->hc_turn_length(*qd);
    if (length1 < length2)
    {
      *q1 = qa;
      *q2 = qb;
      *ci = middle1;
      delete qc;
      delete qd;
      delete middle2;
      return length1;
    }
    else
    {
      *q1 = qc;
      *q2 = qd;
      *ci = middle2;
      delete qa;
      delete qb;
      delete middle1;
      return length2;
    }
    return numeric_limits<double>::max();
  }

  // ##### TcTT #################################################################
  bool TcTT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return (distance <= 2 * c2.radius + 2 / fabs(c2.kappa)) && (distance >= 2 * c2.radius - 2 / fabs(c2.kappa));
  }

  void TcTT_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2,
                            Configuration **q3, Configuration **q4) const
  {
    double theta = angle;
    double r1 = 2 / fabs(c2.kappa);
    double r2 = 2 * c2.radius;
    double delta_x = (pow(r1, 2) + pow(distance, 2) - pow(r2, 2)) / (2 * distance);
    double delta_y = sqrt(fabs(pow(r1, 2) - pow(delta_x, 2)));
    double x, y;

    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c1.left, !c1.forward, c1.regular, parent_->hc_cc_circle_param_);
    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    HC_CC_Circle tgt2(x, y, !c1.left, !c1.forward, c1.regular, parent_->hc_cc_circle_param_);

    TcT_tangent_circles(c1, tgt1, q1);
    TT_tangent_circles(tgt1, c2, q2);
    TcT_tangent_circles(c1, tgt2, q3);
    TT_tangent_circles(tgt2, c2, q4);
  }

  double TcTT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                   Configuration **q1, Configuration **q2, HC_CC_Circle **ci) const
  {
    Configuration *qa, *qb, *qc, *qd;
    TcTT_tangent_circles(c1, c2, &qa, &qb, &qc, &qd);
    HC_CC_Circle *middle1, *middle2;
    middle1 = new HC_CC_Circle(*qb, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    middle2 = new HC_CC_Circle(*qd, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);

    *cstart = new HC_CC_Circle(c1);
    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);

    // select shortest connection
    double length1 = (*cstart)->rs_turn_length(*qa) + middle1->hc_turn_length(*qa) + (*cend)->cc_turn_length(*qb);
    double length2 = (*cstart)->rs_turn_length(*qc) + middle2->hc_turn_length(*qc) + (*cend)->cc_turn_length(*qd);
    if (length1 < length2)
    {
      *q1 = qa;
      *q2 = qb;
      *ci = middle1;
      delete qc;
      delete qd;
      delete middle2;
      return length1;
    }
    else
    {
      *q1 = qc;
      *q2 = qd;
      *ci = middle2;
      delete qa;
      delete qb;
      delete middle1;
      return length2;
    }
    return numeric_limits<double>::max();
  }

  // ##### TTcT #################################################################
  bool TTcT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return (distance <= 2 * c2.radius + 2 / fabs(c2.kappa)) && (distance >= 2 * c2.radius - 2 / fabs(c2.kappa));
  }

  void TTcT_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2,
                            Configuration **q3, Configuration **q4) const
  {
    double theta = angle;
    double r1 = 2 * c2.radius;
    double r2 = 2 / fabs(c2.kappa);
    double delta_x = (pow(r1, 2) + pow(distance, 2) - pow(r2, 2)) / (2 * distance);
    double delta_y = sqrt(fabs(pow(r1, 2) - pow(delta_x, 2)));
    double x, y;

    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c1.left, c1.forward, c1.regular, parent_->hc_cc_circle_param_);
    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    HC_CC_Circle tgt2(x, y, !c1.left, c1.forward, c1.regular, parent_->hc_cc_circle_param_);

    TT_tangent_circles(c1, tgt1, q1);
    TcT_tangent_circles(tgt1, c2, q2);
    TT_tangent_circles(c1, tgt2, q3);
    TcT_tangent_circles(tgt2, c2, q4);
  }

  double TTcT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                   Configuration **q1, Configuration **q2, HC_CC_Circle **ci) const
  {
    Configuration *qa, *qb, *qc, *qd;
    TTcT_tangent_circles(c1, c2, &qa, &qb, &qc, &qd);
    HC_CC_Circle *start1, *start2, *middle1, *middle2;
    start1 = new HC_CC_Circle(*qa, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    middle1 = new HC_CC_Circle(*qa, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    start2 = new HC_CC_Circle(*qc, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    middle2 = new HC_CC_Circle(*qc, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);

    *cend = new HC_CC_Circle(c2);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);

    // select shortest connection
    double length1 = start1->hc_turn_length(**q1) + middle1->hc_turn_length(*qb) + (*cend)->hc_turn_length(*qb);
    double length2 = start2->hc_turn_length(**q1) + middle2->hc_turn_length(*qd) + (*cend)->hc_turn_length(*qd);
    if (length1 < length2)
    {
      *cstart = start1;
      *q2 = qb;
      *ci = middle1;
      delete qa;
      delete qc;
      delete qd;
      delete start2;
      delete middle2;
      return length1;
    }
    else
    {
      *cstart = start2;
      *q2 = qd;
      *ci = middle2;
      delete qa;
      delete qb;
      delete qc;
      delete start1;
      delete middle1;
      return length2;
    }
    return numeric_limits<double>::max();
  }

  // ##### TST ##################################################################
  bool TiST_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return (distance >= 2 * c2.radius);
  }

  bool TeST_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return (distance >= 2 * c2.radius * c2.sin_mu);
  }

  bool TST_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    return TiST_exists(c1, c2) || TeST_exists(c1, c2);
  }

  void TiST_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1,
                            Configuration **q2) const
  {
    double distance = center_distance(c1, c2);
    double angle = atan2(c2.yc - c1.yc, c2.xc - c1.xc);
    double alpha = fabs(asin(2 * c2.radius * c2.cos_mu / distance));
    double delta_x = fabs(c2.radius * c2.sin_mu);
    double delta_y = fabs(c2.radius * c2.cos_mu);
    double x, y, theta;
    if (c1.left && c1.forward)
    {
      theta = angle + alpha;
      global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
    }
    if (c1.left && !c1.forward)
    {
      theta = angle - alpha;
      global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
    }
    if (!c1.left && c1.forward)
    {
      theta = angle - alpha;
      global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
    }
    if (!c1.left && !c1.forward)
    {
      theta = angle + alpha;
      global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
    }
  }

  void TeST_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1,
                            Configuration **q2) const
  {
    double delta_x = fabs(c2.radius * c2.sin_mu);
    double delta_y = fabs(c2.radius * c2.cos_mu);
    double theta = atan2(c2.yc - c1.yc, c2.xc - c1.xc);
    double x, y;
    if (c1.left && c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
    }
    if (c1.left && !c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
    }
    if (!c1.left && c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
    }
    if (!c1.left && !c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
    }
  }

  double TiST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                   Configuration **q1, Configuration **q2, Configuration **q3) const
  {
    TiST_tangent_circles(c1, c2, q2, q3);
    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    return (*cstart)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) + (*cend)->cc_turn_length(**q3);
  }

  double TeST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                   Configuration **q1, Configuration **q2, Configuration **q3) const
  {
    TeST_tangent_circles(c1, c2, q2, q3);
    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    return (*cstart)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) + (*cend)->cc_turn_length(**q3);
  }

  double TST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                  Configuration **q1, Configuration **q2, Configuration **q3) const
  {
    if (TiST_exists(c1, c2))
    {
      return TiST_path(c1, c2, cstart, cend, q1, q2, q3);
    }
    if (TeST_exists(c1, c2))
    {
      return TeST_path(c1, c2, cstart, cend, q1, q2, q3);
    }
    return numeric_limits<double>::max();
  }

  // ##### TSTcT ################################################################
  bool TiSTcT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return (distance >=
            sqrt(pow(2 * c2.radius * c2.sin_mu + 2 / fabs(c2.kappa), 2) + pow(2 * c2.radius * c2.cos_mu, 2)));
  }

  bool TeSTcT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return (distance >= 2 * (1 / fabs(c2.kappa) + c2.radius * c2.sin_mu));
  }

  bool TSTcT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    return TiSTcT_exists(c1, c2) || TeSTcT_exists(c1, c2);
  }

  double TiSTcT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                     Configuration **q1, Configuration **q2, Configuration **q3, Configuration **q4,
                     HC_CC_Circle **ci) const
  {
    double theta = angle;
    double delta_y = (4 * c2.radius * c2.cos_mu) / (fabs(c2.kappa) * distance);
    double delta_x = sqrt(pow(2 / c2.kappa, 2) - pow(delta_y, 2));
    double x, y;

    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c2.left, c2.forward, c2.regular, parent_->hc_cc_circle_param_);

    TiST_tangent_circles(c1, tgt1, q2, q3);
    TcT_tangent_circles(tgt1, c2, q4);

    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(c2);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    *ci = new HC_CC_Circle(**q3, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);

    return (*cstart)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) + (*ci)->hc_turn_length(**q4) +
           (*cend)->hc_turn_length(**q4);
  }

  double TeSTcT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                     Configuration **q1, Configuration **q2, Configuration **q3, Configuration **q4,
                     HC_CC_Circle **ci) const
  {
    double theta = angle;
    double delta_x = 2 / fabs(c2.kappa);
    double delta_y = 0;
    double x, y;

    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c2.left, c2.forward, c2.regular, parent_->hc_cc_circle_param_);

    TeST_tangent_circles(c1, tgt1, q2, q3);
    TcT_tangent_circles(tgt1, c2, q4);

    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(c2);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    *ci = new HC_CC_Circle(**q3, !c2.left, c2.forward, true, parent_->hc_cc_circle_param_);

    return (*cstart)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) + (*ci)->hc_turn_length(**q4) +
           (*cend)->hc_turn_length(**q4);
  }

  double TSTcT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                    Configuration **q1, Configuration **q2, Configuration **q3, Configuration **q4,
                    HC_CC_Circle **ci) const
  {
    if (TiSTcT_exists(c1, c2))
    {
      return TiSTcT_path(c1, c2, cstart, cend, q1, q2, q3, q4, ci);
    }
    if (TeSTcT_exists(c1, c2))
    {
      return TeSTcT_path(c1, c2, cstart, cend, q1, q2, q3, q4, ci);
    }
    return numeric_limits<double>::max();
  }

  // ##### TcTST ################################################################
  bool TcTiST_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return (distance >=
            sqrt(pow(2 * c2.radius * c2.sin_mu + 2 / fabs(c2.kappa), 2) + pow(2 * c2.radius * c2.cos_mu, 2)));
  }

  bool TcTeST_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return (distance >= 2 * (1 / fabs(c2.kappa) + c2.radius * c2.sin_mu));
  }

  bool TcTST_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    return TcTiST_exists(c1, c2) || TcTeST_exists(c1, c2);
  }

  double TcTiST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                     Configuration **q1, Configuration **q2, Configuration **q3, HC_CC_Circle **ci) const
  {
    double theta = angle;
    double delta_y = (4 * c2.radius * c2.cos_mu) / (fabs(c2.kappa) * distance);
    double delta_x = sqrt(pow(2 / c2.kappa, 2) - pow(delta_y, 2));
    double x, y;

    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c1.left, !c1.forward, c1.regular, parent_->hc_cc_circle_param_);

    TcT_tangent_circles(c1, tgt1, q1);
    TiST_tangent_circles(tgt1, c2, q2, q3);

    *cstart = new HC_CC_Circle(c1);
    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *ci = new HC_CC_Circle(**q2, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);

    return (*cstart)->rs_turn_length(**q1) + (*ci)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) +
           (*cend)->cc_turn_length(**q3);
  }

  double TcTeST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                     Configuration **q1, Configuration **q2, Configuration **q3, HC_CC_Circle **ci) const
  {
    double theta = angle;
    double delta_x = 2 / fabs(c2.kappa);
    double delta_y = 0;
    double x, y;

    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c1.left, !c1.forward, c1.regular, parent_->hc_cc_circle_param_);

    TcT_tangent_circles(c1, tgt1, q1);
    TeST_tangent_circles(tgt1, c2, q2, q3);

    *cstart = new HC_CC_Circle(c1);
    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *ci = new HC_CC_Circle(**q2, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);

    return (*cstart)->rs_turn_length(**q1) + (*ci)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) +
           (*cend)->cc_turn_length(**q3);
  }

  double TcTST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                    Configuration **q1, Configuration **q2, Configuration **q3, HC_CC_Circle **ci) const
  {
    if (TcTiST_exists(c1, c2))
    {
      return TcTiST_path(c1, c2, cstart, cend, q1, q2, q3, ci);
    }
    if (TcTeST_exists(c1, c2))
    {
      return TcTeST_path(c1, c2, cstart, cend, q1, q2, q3, ci);
    }
    return numeric_limits<double>::max();
  }

  // ##### TcTSTcT ##############################################################
  bool TcTiSTcT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return (distance >=
            sqrt(pow(2 * c2.radius, 2) + 16 * c2.radius * c2.sin_mu / fabs(c2.kappa) + pow(4 / c2.kappa, 2)));
  }

  bool TcTeSTcT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return (distance >= 4 / fabs(c2.kappa) + 2 * c2.radius * c2.sin_mu);
  }

  bool TcTSTcT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    return TcTiSTcT_exists(c1, c2) || TcTeSTcT_exists(c1, c2);
  }

  double TcTiSTcT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                       Configuration **q1, Configuration **q2, Configuration **q3, Configuration **q4,
                       HC_CC_Circle **ci1, HC_CC_Circle **ci2) const
  {
    double theta = angle;
    double delta_y = (4 * c2.radius * c2.cos_mu) / (distance * fabs(c2.kappa));
    double delta_x = sqrt(pow(2 / c2.kappa, 2) - pow(delta_y, 2));
    double x, y;

    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c1.left, !c1.forward, c1.regular, parent_->hc_cc_circle_param_);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
    HC_CC_Circle tgt2(x, y, !c2.left, c2.forward, c2.regular, parent_->hc_cc_circle_param_);

    TcT_tangent_circles(c1, tgt1, q1);
    TiST_tangent_circles(tgt1, tgt2, q2, q3);
    TcT_tangent_circles(tgt2, c2, q4);

    *cstart = new HC_CC_Circle(c1);
    *cend = new HC_CC_Circle(c2);
    *ci1 = new HC_CC_Circle(**q2, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    *ci2 = new HC_CC_Circle(**q3, !c2.left, c2.forward, true, parent_->hc_cc_circle_param_);

    return (*cstart)->rs_turn_length(**q1) + (*ci1)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) +
           (*ci2)->hc_turn_length(**q4) + (*cend)->hc_turn_length(**q4);
  }

  double TcTeSTcT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                       Configuration **q1, Configuration **q2, Configuration **q3, Configuration **q4,
                       HC_CC_Circle **ci1, HC_CC_Circle **ci2) const
  {
    double theta = angle;
    double delta_x = 2 / fabs(c1.kappa);
    double delta_y = 0;
    double x, y;

    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c1.left, !c1.forward, c1.regular, parent_->hc_cc_circle_param_);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt2(x, y, !c2.left, c2.forward, c2.regular, parent_->hc_cc_circle_param_);

    TcT_tangent_circles(c1, tgt1, q1);
    TeST_tangent_circles(tgt1, tgt2, q2, q3);
    TcT_tangent_circles(tgt2, c2, q4);

    *cstart = new HC_CC_Circle(c1);
    *cend = new HC_CC_Circle(c2);
    *ci1 = new HC_CC_Circle(**q2, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    *ci2 = new HC_CC_Circle(**q3, !c2.left, c2.forward, true, parent_->hc_cc_circle_param_);

    return (*cstart)->rs_turn_length(**q1) + (*ci1)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) +
           (*ci2)->hc_turn_length(**q4) + (*cend)->hc_turn_length(**q4);
  }

  double TcTSTcT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                      Configuration **q1, Configuration **q2, Configuration **q3, Configuration **q4,
                      HC_CC_Circle **ci1, HC_CC_Circle **ci2) const
  {
    if (TcTiSTcT_exists(c1, c2))
    {
      return TcTiSTcT_path(c1, c2, cstart, cend, q1, q2, q3, q4, ci1, ci2);
    }
    if (TcTeSTcT_exists(c1, c2))
    {
      return TcTeSTcT_path(c1, c2, cstart, cend, q1, q2, q3, q4, ci1, ci2);
    }
    return numeric_limits<double>::max();
  }

  // ##### TTcTT ###############################################################
  bool TTcTT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return (distance <= 4 * c2.radius + 2 / fabs(c2.kappa));
  }

  void TTcTT_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2,
                             Configuration **q3, Configuration **q4, Configuration **q5, Configuration **q6) const
  {
    double theta = angle;
    double r1, r2, delta_x, delta_y, x, y;
    r1 = 2 / fabs(c2.kappa);
    r2 = 2 * c2.radius;
    if (distance < 4 * c2.radius - 2 / fabs(c2.kappa))
    {
      delta_x = (distance + r1) / 2;
      delta_y = sqrt(fabs((pow(r2, 2) - pow(delta_x, 2))));
    }
    else
    {
      delta_x = (distance - r1) / 2;
      delta_y = sqrt(fabs((pow(r2, 2) - pow(delta_x, 2))));
    }

    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c1.left, c1.forward, c1.regular, parent_->hc_cc_circle_param_);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt2(x, y, !c2.left, !c2.forward, c2.regular, parent_->hc_cc_circle_param_);

    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    HC_CC_Circle tgt3(x, y, !c1.left, c1.forward, c1.regular, parent_->hc_cc_circle_param_);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
    HC_CC_Circle tgt4(x, y, !c2.left, !c2.forward, c2.regular, parent_->hc_cc_circle_param_);

    TT_tangent_circles(c1, tgt1, q1);
    TcT_tangent_circles(tgt1, tgt2, q2);
    TT_tangent_circles(tgt2, c2, q3);

    TT_tangent_circles(c1, tgt3, q4);
    TcT_tangent_circles(tgt3, tgt4, q5);
    TT_tangent_circles(tgt4, c2, q6);
  }

  double TTcTT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                    Configuration **q1, Configuration **q2, Configuration **q3, HC_CC_Circle **ci1,
                    HC_CC_Circle **ci2) const
  {
    Configuration *qa, *qb, *qc, *qd, *qe, *qf;
    TTcTT_tangent_circles(c1, c2, &qa, &qb, &qc, &qd, &qe, &qf);
    HC_CC_Circle *start1, *start2, *middle1, *middle2, *middle3, *middle4;
    start1 = new HC_CC_Circle(*qa, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    middle1 = new HC_CC_Circle(*qa, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    middle2 = new HC_CC_Circle(*qc, !c2.left, c2.forward, true, parent_->hc_cc_circle_param_);
    start2 = new HC_CC_Circle(*qd, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    middle3 = new HC_CC_Circle(*qd, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    middle4 = new HC_CC_Circle(*qf, !c2.left, c2.forward, true, parent_->hc_cc_circle_param_);

    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);

    // select shortest connection
    double length1 = start1->hc_turn_length(**q1) + middle1->hc_turn_length(*qb) + middle2->hc_turn_length(*qb) +
                     (*cend)->cc_turn_length(*qc);
    double length2 = start2->hc_turn_length(**q1) + middle3->hc_turn_length(*qe) + middle4->hc_turn_length(*qe) +
                     (*cend)->cc_turn_length(*qf);
    if (length1 < length2)
    {
      *cstart = start1;
      *q2 = qb;
      *q3 = qc;
      *ci1 = middle1;
      *ci2 = middle2;
      delete qa;
      delete qd;
      delete qe;
      delete qf;
      delete start2;
      delete middle3;
      delete middle4;
      return length1;
    }
    else
    {
      *cstart = start2;
      *q2 = qe;
      *q3 = qf;
      *ci1 = middle3;
      *ci2 = middle4;
      delete qa;
      delete qb;
      delete qc;
      delete qd;
      delete start1;
      delete middle1;
      delete middle2;
      return length2;
    }

    return numeric_limits<double>::max();
  }

  // ##### TcTTcT ###############################################################
  bool TcTTcT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return (distance <= 4 / fabs(c2.kappa) + 2 * c2.radius) && (distance >= 4 / fabs(c2.kappa) - 2 * c2.radius);
  }

  void TcTTcT_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2,
                              Configuration **q3, Configuration **q4, Configuration **q5, Configuration **q6) const
  {
    double theta = angle;
    double r1 = 2 / fabs(c2.kappa);
    double r2 = c2.radius;
    double delta_x = (pow(r1, 2) + pow(distance / 2, 2) - pow(r2, 2)) / distance;
    double delta_y = sqrt(fabs(pow(r1, 2) - pow(delta_x, 2)));
    double x, y;

    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c1.left, !c1.forward, c1.regular, parent_->hc_cc_circle_param_);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
    HC_CC_Circle tgt2(x, y, !c2.left, c2.forward, c2.regular, parent_->hc_cc_circle_param_);

    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    HC_CC_Circle tgt3(x, y, !c1.left, !c1.forward, c1.regular, parent_->hc_cc_circle_param_);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt4(x, y, !c2.left, c2.forward, c2.regular, parent_->hc_cc_circle_param_);

    TcT_tangent_circles(c1, tgt1, q1);
    TT_tangent_circles(tgt1, tgt2, q2);
    TcT_tangent_circles(tgt2, c2, q3);

    TcT_tangent_circles(c1, tgt3, q4);
    TT_tangent_circles(tgt3, tgt4, q5);
    TcT_tangent_circles(tgt4, c2, q6);
  }

  double TcTTcT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                     Configuration **q1, Configuration **q2, HC_CC_Circle **ci1, HC_CC_Circle **ci2) const
  {
    Configuration *qa, *qb, *qc, *qd, *qe, *qf;
    TcTTcT_tangent_circles(c1, c2, &qa, &qb, &qc, &qd, &qe, &qf);
    HC_CC_Circle *middle1, *middle2, *middle3, *middle4;
    middle1 = new HC_CC_Circle(*qb, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    middle2 = new HC_CC_Circle(*qb, c1.left, !c1.forward, true, parent_->hc_cc_circle_param_);
    middle3 = new HC_CC_Circle(*qe, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    middle4 = new HC_CC_Circle(*qe, c1.left, !c1.forward, true, parent_->hc_cc_circle_param_);

    *cstart = new HC_CC_Circle(c1);
    *cend = new HC_CC_Circle(c2);

    // select shortest connection
    double length1 = (*cstart)->rs_turn_length(*qa) + middle1->hc_turn_length(*qa) + middle2->hc_turn_length(*qc) +
                     (*cend)->hc_turn_length(*qc);
    double length2 = (*cstart)->rs_turn_length(*qd) + middle3->hc_turn_length(*qd) + middle4->hc_turn_length(*qf) +
                     (*cend)->hc_turn_length(*qf);
    if (length1 < length2)
    {
      *q1 = qa;
      *q2 = qc;
      *ci1 = middle1;
      *ci2 = middle2;
      delete qb;
      delete qd;
      delete qe;
      delete qf;
      delete middle3;
      delete middle4;
      return length1;
    }
    else
    {
      *q1 = qd;
      *q2 = qf;
      *ci1 = middle3;
      *ci2 = middle4;
      delete qa;
      delete qb;
      delete qc;
      delete qe;
      delete middle1;
      delete middle2;
      return length2;
    }
    return numeric_limits<double>::max();
  }

  // ############################################################################

  // ##### TTT ##################################################################
  bool TTT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return distance <= 4 * c2.radius;
  }

  void TTT_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2,
                           Configuration **q3, Configuration **q4) const
  {
    double theta = angle;
    double r = 2 * c2.radius;
    double delta_x = 0.5 * distance;
    double delta_y = sqrt(fabs(pow(delta_x, 2) - pow(r, 2)));
    double x, y;

    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    HC_CC_Circle tgt1(x, y, !c1.left, c1.forward, c1.regular, parent_->hc_cc_circle_param_);
    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    HC_CC_Circle tgt2(x, y, !c1.left, c1.forward, c1.regular, parent_->hc_cc_circle_param_);

    TT_tangent_circles(c1, tgt1, q1);
    TT_tangent_circles(tgt1, c2, q2);
    TT_tangent_circles(c1, tgt2, q3);
    TT_tangent_circles(tgt2, c2, q4);
  }

  double TTT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                  Configuration **q1, Configuration **q2, HC_CC_Circle **ci) const
  {
    Configuration *qa, *qb, *qc, *qd;
    TTT_tangent_circles(c1, c2, &qa, &qb, &qc, &qd);
    HC_CC_Circle *start1, *start2, *middle1, *middle2;
    start1 = new HC_CC_Circle(*qa, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    middle1 = new HC_CC_Circle(*qa, !c1.left, c1.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    start2 = new HC_CC_Circle(*qc, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    middle2 = new HC_CC_Circle(*qc, !c1.left, c1.forward, CC_REGULAR, parent_->hc_cc_circle_param_);

    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);

    // select shortest connection
    double length1 = start1->hc_turn_length(**q1) + middle1->cc_turn_length(*qb) + (*cend)->cc_turn_length(*qb);
    double length2 = start2->hc_turn_length(**q1) + middle2->cc_turn_length(*qd) + (*cend)->cc_turn_length(*qd);
    if (length1 < length2)
    {
      *cstart = start1;
      *q2 = qb;
      *ci = middle1;
      delete qa;
      delete qc;
      delete qd;
      delete start2;
      delete middle2;
      return length1;
    }
    else
    {
      *cstart = start2;
      *q2 = qd;
      *ci = middle2;
      delete qa;
      delete qb;
      delete qc;
      delete start1;
      delete middle1;
      return length2;
    }
    return numeric_limits<double>::max();
  }

  // ##### TcST ################################################################
  bool TciST_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return distance >= 2 * c2.radius * c2.cos_mu;
  }

  bool TceST_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return distance >= get_epsilon();
  }

  bool TcST_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    return TciST_exists(c1, c2) || TceST_exists(c1, c2);
  }

  double TciST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                    Configuration **q1, Configuration **q2, Configuration **q3) const
  {
    double alpha = fabs(asin(2 * c2.radius * c2.cos_mu / distance));
    double delta_x = fabs(c2.radius * c2.sin_mu);
    double delta_y = fabs(c2.radius * c2.cos_mu);
    double x, y, theta;
    if (c1.left && c1.forward)
    {
      theta = angle - alpha;
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta + PI, 0);
    }
    if (c1.left && !c1.forward)
    {
      theta = angle + alpha;
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta, 0);
    }
    if (!c1.left && c1.forward)
    {
      theta = angle + alpha;
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta + PI, 0);
    }
    if (!c1.left && !c1.forward)
    {
      theta = angle - alpha;
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta, 0);
    }
    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    return (*cstart)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) + (*cend)->cc_turn_length(**q3);
  }

  double TceST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                    Configuration **q1, Configuration **q2, Configuration **q3) const
  {
    double delta_x = fabs(c2.radius * c2.sin_mu);
    double delta_y = fabs(c2.radius * c2.cos_mu);
    double theta = angle;
    double x, y;
    if (c1.left && c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta + PI, 0);
    }
    if (c1.left && !c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta, 0);
    }
    if (!c1.left && c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta + PI, 0);
    }
    if (!c1.left && !c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta, 0);
    }
    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    return (*cstart)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) + (*cend)->cc_turn_length(**q3);
  }

  double TcST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                   Configuration **q1, Configuration **q2, Configuration **q3) const
  {
    if (TciST_exists(c1, c2))
    {
      return TciST_path(c1, c2, cstart, cend, q1, q2, q3);
    }
    if (TceST_exists(c1, c2))
    {
      return TceST_path(c1, c2, cstart, cend, q1, q2, q3);
    }
    return numeric_limits<double>::max();
  }

  // ##### TScT #################################################################
  bool TiScT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return distance >= 2 * c2.radius * c2.cos_mu;
  }

  bool TeScT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward != c2.forward)
    {
      return false;
    }
    return distance >= get_epsilon();
  }

  bool TScT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    return TiScT_exists(c1, c2) || TeScT_exists(c1, c2);
  }

  double TiScT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                    Configuration **q1, Configuration **q2, Configuration **q3) const
  {
    double alpha = fabs(asin(2 * c2.radius * c2.cos_mu / distance));

    double delta_x = fabs(c2.radius * c2.sin_mu);
    double delta_y = fabs(c2.radius * c2.cos_mu);

    double x, y, theta;
    if (c1.left && c1.forward)
    {
      theta = angle + alpha;
      global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta, 0);
    }
    if (c1.left && !c1.forward)
    {
      theta = angle - alpha;
      global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, -delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta + PI, 0);
    }
    if (!c1.left && c1.forward)
    {
      theta = angle - alpha;
      global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, -delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta, 0);
    }
    if (!c1.left && !c1.forward)
    {
      theta = angle + alpha;
      global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta + PI, 0);
    }
    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    return (*cstart)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) + (*cend)->cc_turn_length(**q3);
  }

  double TeScT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                    Configuration **q1, Configuration **q2, Configuration **q3) const
  {
    double delta_x = fabs(c2.radius * c2.sin_mu);
    double delta_y = fabs(c2.radius * c2.cos_mu);
    double theta = angle;
    double x, y;
    if (c1.left && c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, -delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta, 0);
    }
    if (c1.left && !c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta + PI, 0);
    }
    if (!c1.left && c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, 0);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta, 0);
    }
    if (!c1.left && !c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, 0);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, -delta_y, &x, &y);
      *q3 = new Configuration(x, y, theta + PI, 0);
    }
    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, HC_REGULAR, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(c2.start, c2.left, c2.forward, CC_REGULAR, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    return (*cstart)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) + (*cend)->cc_turn_length(**q3);
  }

  double TScT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                   Configuration **q1, Configuration **q2, Configuration **q3) const
  {
    if (TiScT_exists(c1, c2))
    {
      return TiScT_path(c1, c2, cstart, cend, q1, q2, q3);
    }
    if (TeScT_exists(c1, c2))
    {
      return TeScT_path(c1, c2, cstart, cend, q1, q2, q3);
    }
    return numeric_limits<double>::max();
  }

  // ##### TcScT ################################################################
  bool TciScT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return distance >= fabs(2 / c2.kappa);
  }

  bool TceScT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left != c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return distance >= get_epsilon();
  }

  bool TcScT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    return TciScT_exists(c1, c2) || TceScT_exists(c1, c2);
  }

  double TciScT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                     Configuration **q1, Configuration **q2) const
  {
    double alpha = fabs(asin(2 / (c2.kappa * distance)));
    double delta_x = 0.0;
    double delta_y = fabs(1 / c2.kappa);
    double x, y, theta;
    if (c1.left && c1.forward)
    {
      theta = angle - alpha;
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta + PI, c1.kappa);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, c2.kappa);
    }
    if (c1.left && !c1.forward)
    {
      theta = angle + alpha;
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, -delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta, c1.kappa);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, c2.kappa);
    }
    if (!c1.left && c1.forward)
    {
      theta = angle + alpha;
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, -delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta + PI, c1.kappa);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, c2.kappa);
    }
    if (!c1.left && !c1.forward)
    {
      theta = angle - alpha;
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta, c1.kappa);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, c2.kappa);
    }
    *cstart = new HC_CC_Circle(c1);
    *cend = new HC_CC_Circle(c2);
    return (*cstart)->rs_turn_length(**q1) + configuration_distance(**q1, **q2) + (*cend)->hc_turn_length(**q2);
  }

  double TceScT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                     Configuration **q1, Configuration **q2) const
  {
    double delta_x = 0.0;
    double delta_y = fabs(1 / c2.kappa);
    double theta = angle;
    double x, y;
    if (c1.left && c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta + PI, c1.kappa);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, c2.kappa);
    }
    if (c1.left && !c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, -delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta, c1.kappa);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, c2.kappa);
    }
    if (!c1.left && c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, -delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta + PI, c1.kappa);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, -delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta + PI, c2.kappa);
    }
    if (!c1.left && !c1.forward)
    {
      global_frame_change(c1.xc, c1.yc, theta, -delta_x, delta_y, &x, &y);
      *q1 = new Configuration(x, y, theta, c1.kappa);
      global_frame_change(c2.xc, c2.yc, theta, delta_x, delta_y, &x, &y);
      *q2 = new Configuration(x, y, theta, c2.kappa);
    }
    *cstart = new HC_CC_Circle(c1);
    *cend = new HC_CC_Circle(c2);
    return (*cstart)->rs_turn_length(**q1) + configuration_distance(**q1, **q2) + (*cend)->hc_turn_length(**q2);
  }

  double TcScT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                    Configuration **q1, Configuration **q2) const
  {
    if (TciScT_exists(c1, c2))
    {
      return TciScT_path(c1, c2, cstart, cend, q1, q2);
    }
    if (TceScT_exists(c1, c2))
    {
      return TceScT_path(c1, c2, cstart, cend, q1, q2);
    }
    return numeric_limits<double>::max();
  }
};

// ############################################################################

HCpm0_Reeds_Shepp_State_Space::HCpm0_Reeds_Shepp_State_Space(double kappa, double sigma, double discretization)
  : HC_CC_State_Space(kappa, sigma, discretization)
  , hcpm0_reeds_shepp_{ unique_ptr<HCpm0_Reeds_Shepp>(new HCpm0_Reeds_Shepp(this)) }
{
  rs_circle_param_.set_param(kappa_, numeric_limits<double>::max(), 1 / kappa_, 0.0, 0.0, 1.0, 0.0);
  radius_ = hc_cc_circle_param_.radius;
  mu_ = hc_cc_circle_param_.mu;
}

HCpm0_Reeds_Shepp_State_Space::~HCpm0_Reeds_Shepp_State_Space() = default;

HC_CC_RS_Path *HCpm0_Reeds_Shepp_State_Space::hcpm0_circles_rs_path(const HC_CC_Circle &c1,
                                                                    const HC_CC_Circle &c2) const
{
  // table containing the lengths of the paths, the intermediate configurations and circles
  double length[nb_hc_cc_rs_paths];
  double_array_init(length, nb_hc_cc_rs_paths, numeric_limits<double>::max());
  Configuration *qi1[nb_hc_cc_rs_paths];
  pointer_array_init((void **)qi1, nb_hc_cc_rs_paths);
  Configuration *qi2[nb_hc_cc_rs_paths];
  pointer_array_init((void **)qi2, nb_hc_cc_rs_paths);
  Configuration *qi3[nb_hc_cc_rs_paths];
  pointer_array_init((void **)qi3, nb_hc_cc_rs_paths);
  Configuration *qi4[nb_hc_cc_rs_paths];
  pointer_array_init((void **)qi4, nb_hc_cc_rs_paths);
  HC_CC_Circle *cstart[nb_hc_cc_rs_paths];
  pointer_array_init((void **)cstart, nb_hc_cc_rs_paths);
  HC_CC_Circle *ci1[nb_hc_cc_rs_paths];
  pointer_array_init((void **)ci1, nb_hc_cc_rs_paths);
  HC_CC_Circle *ci2[nb_hc_cc_rs_paths];
  pointer_array_init((void **)ci2, nb_hc_cc_rs_paths);
  HC_CC_Circle *cend[nb_hc_cc_rs_paths];
  pointer_array_init((void **)cend, nb_hc_cc_rs_paths);

  // precomputations
  hcpm0_reeds_shepp_->distance = center_distance(c1, c2);
  hcpm0_reeds_shepp_->angle = atan2(c2.yc - c1.yc, c2.xc - c1.xc);

  // case E
  if (configuration_equal(c1.start, c2.start))
  {
    length[hc_cc_rs::E] = 0;
    goto label_end;
  }
  // case T
  if (hcpm0_reeds_shepp_->distance < get_epsilon())
  {
    cend[hc_cc_rs::T] = new HC_CC_Circle(c2.start, c2.left, c2.forward, HC_REGULAR, hc_cc_circle_param_);
    length[hc_cc_rs::T] = cend[hc_cc_rs::T]->hc_turn_length(c1.start);
    goto label_end;
  }
  // case TT
  if (hcpm0_reeds_shepp_->TT_exists(c1, c2))
  {
    length[hc_cc_rs::TT] = hcpm0_reeds_shepp_->TT_path(c1, c2, &cstart[hc_cc_rs::TT], &cend[hc_cc_rs::TT],
                                                       &qi1[hc_cc_rs::TT], &qi2[hc_cc_rs::TT]);
  }
  // case TcT
  if (hcpm0_reeds_shepp_->TcT_exists(c1, c2))
  {
    length[hc_cc_rs::TcT] =
        hcpm0_reeds_shepp_->TcT_path(c1, c2, &cstart[hc_cc_rs::TcT], &cend[hc_cc_rs::TcT], &qi1[hc_cc_rs::TcT]);
  }
  // ##### Reeds-Shepp families: ############################################
  // case TcTcT
  if (hcpm0_reeds_shepp_->TcTcT_exists(c1, c2))
  {
    length[hc_cc_rs::TcTcT] =
        hcpm0_reeds_shepp_->TcTcT_path(c1, c2, &cstart[hc_cc_rs::TcTcT], &cend[hc_cc_rs::TcTcT], &qi1[hc_cc_rs::TcTcT],
                                       &qi2[hc_cc_rs::TcTcT], &ci1[hc_cc_rs::TcTcT]);
  }
  // case TcTT
  if (hcpm0_reeds_shepp_->TcTT_exists(c1, c2))
  {
    length[hc_cc_rs::TcTT] =
        hcpm0_reeds_shepp_->TcTT_path(c1, c2, &cstart[hc_cc_rs::TcTT], &cend[hc_cc_rs::TcTT], &qi1[hc_cc_rs::TcTT],
                                      &qi2[hc_cc_rs::TcTT], &ci1[hc_cc_rs::TcTT]);
  }
  // case TTcT
  if (hcpm0_reeds_shepp_->TTcT_exists(c1, c2))
  {
    length[hc_cc_rs::TTcT] =
        hcpm0_reeds_shepp_->TTcT_path(c1, c2, &cstart[hc_cc_rs::TTcT], &cend[hc_cc_rs::TTcT], &qi1[hc_cc_rs::TTcT],
                                      &qi2[hc_cc_rs::TTcT], &ci1[hc_cc_rs::TTcT]);
  }
  // case TST
  if (hcpm0_reeds_shepp_->TST_exists(c1, c2))
  {
    length[hc_cc_rs::TST] = hcpm0_reeds_shepp_->TST_path(c1, c2, &cstart[hc_cc_rs::TST], &cend[hc_cc_rs::TST],
                                                         &qi1[hc_cc_rs::TST], &qi2[hc_cc_rs::TST], &qi3[hc_cc_rs::TST]);
  }
  // case TSTcT
  if (hcpm0_reeds_shepp_->TSTcT_exists(c1, c2))
  {
    length[hc_cc_rs::TSTcT] = hcpm0_reeds_shepp_->TSTcT_path(
        c1, c2, &cstart[hc_cc_rs::TSTcT], &cend[hc_cc_rs::TSTcT], &qi1[hc_cc_rs::TSTcT], &qi2[hc_cc_rs::TSTcT],
        &qi3[hc_cc_rs::TSTcT], &qi4[hc_cc_rs::TSTcT], &ci1[hc_cc_rs::TSTcT]);
  }
  // case TcTST
  if (hcpm0_reeds_shepp_->TcTST_exists(c1, c2))
  {
    length[hc_cc_rs::TcTST] =
        hcpm0_reeds_shepp_->TcTST_path(c1, c2, &cstart[hc_cc_rs::TcTST], &cend[hc_cc_rs::TcTST], &qi1[hc_cc_rs::TcTST],
                                       &qi2[hc_cc_rs::TcTST], &qi3[hc_cc_rs::TcTST], &ci1[hc_cc_rs::TcTST]);
  }
  // case TcTSTcT
  if (hcpm0_reeds_shepp_->TcTSTcT_exists(c1, c2))
  {
    length[hc_cc_rs::TcTSTcT] = hcpm0_reeds_shepp_->TcTSTcT_path(
        c1, c2, &cstart[hc_cc_rs::TcTSTcT], &cend[hc_cc_rs::TcTSTcT], &qi1[hc_cc_rs::TcTSTcT], &qi2[hc_cc_rs::TcTSTcT],
        &qi3[hc_cc_rs::TcTSTcT], &qi4[hc_cc_rs::TcTSTcT], &ci1[hc_cc_rs::TcTSTcT], &ci2[hc_cc_rs::TcTSTcT]);
  }
  // case TTcTT
  if (hcpm0_reeds_shepp_->TTcTT_exists(c1, c2))
  {
    length[hc_cc_rs::TTcTT] = hcpm0_reeds_shepp_->TTcTT_path(
        c1, c2, &cstart[hc_cc_rs::TTcTT], &cend[hc_cc_rs::TTcTT], &qi1[hc_cc_rs::TTcTT], &qi2[hc_cc_rs::TTcTT],
        &qi3[hc_cc_rs::TTcTT], &ci1[hc_cc_rs::TTcTT], &ci2[hc_cc_rs::TTcTT]);
  }
  // case TcTTcT
  if (hcpm0_reeds_shepp_->TcTTcT_exists(c1, c2))
  {
    length[hc_cc_rs::TcTTcT] = hcpm0_reeds_shepp_->TcTTcT_path(
        c1, c2, &cstart[hc_cc_rs::TcTTcT], &cend[hc_cc_rs::TcTTcT], &qi1[hc_cc_rs::TcTTcT], &qi2[hc_cc_rs::TcTTcT],
        &ci1[hc_cc_rs::TcTTcT], &ci2[hc_cc_rs::TcTTcT]);
  }
  // ############################################################################
  // case TTT
  if (hcpm0_reeds_shepp_->TTT_exists(c1, c2))
  {
    length[hc_cc_rs::TTT] = hcpm0_reeds_shepp_->TTT_path(c1, c2, &cstart[hc_cc_rs::TTT], &cend[hc_cc_rs::TTT],
                                                         &qi1[hc_cc_rs::TTT], &qi2[hc_cc_rs::TTT], &ci1[hc_cc_rs::TTT]);
  }
  // case TcST
  if (hcpm0_reeds_shepp_->TcST_exists(c1, c2))
  {
    length[hc_cc_rs::TcST] =
        hcpm0_reeds_shepp_->TcST_path(c1, c2, &cstart[hc_cc_rs::TcST], &cend[hc_cc_rs::TcST], &qi1[hc_cc_rs::TcST],
                                      &qi2[hc_cc_rs::TcST], &qi3[hc_cc_rs::TcST]);
  }
  // case TScT
  if (hcpm0_reeds_shepp_->TScT_exists(c1, c2))
  {
    length[hc_cc_rs::TScT] =
        hcpm0_reeds_shepp_->TScT_path(c1, c2, &cstart[hc_cc_rs::TScT], &cend[hc_cc_rs::TScT], &qi1[hc_cc_rs::TScT],
                                      &qi2[hc_cc_rs::TScT], &qi3[hc_cc_rs::TScT]);
  }
  // case TcScT
  if (hcpm0_reeds_shepp_->TcScT_exists(c1, c2))
  {
    length[hc_cc_rs::TcScT] = hcpm0_reeds_shepp_->TcScT_path(c1, c2, &cstart[hc_cc_rs::TcScT], &cend[hc_cc_rs::TcScT],
                                                             &qi1[hc_cc_rs::TcScT], &qi2[hc_cc_rs::TcScT]);
  }
label_end:
  // select shortest path
  hc_cc_rs::path_type best_path = (hc_cc_rs::path_type)array_index_min(length, nb_hc_cc_rs_paths);
  HC_CC_RS_Path *path;
  path = new HC_CC_RS_Path(c1.start, c2.start, best_path, kappa_, sigma_, qi1[best_path], qi2[best_path],
                           qi3[best_path], qi4[best_path], cstart[best_path], cend[best_path], ci1[best_path],
                           ci2[best_path], length[best_path]);

  // clean up
  for (int i = 0; i < nb_hc_cc_rs_paths; i++)
  {
    if (i != best_path)
    {
      delete qi1[i];
      delete qi2[i];
      delete qi3[i];
      delete qi4[i];
      delete cstart[i];
      delete ci1[i];
      delete ci2[i];
      delete cend[i];
    }
  }
  return path;
}

HC_CC_RS_Path *HCpm0_Reeds_Shepp_State_Space::hcpm0_reeds_shepp(const State &state1, const State &state2) const
{
  // compute the 4 circles at the intial and final configuration
  Configuration start1(state1.x, state1.y, state1.theta, kappa_);
  Configuration start2(state1.x, state1.y, state1.theta, -kappa_);
  Configuration end(state2.x, state2.y, state2.theta, 0.0);

  HC_CC_Circle *start_circle[4];
  HC_CC_Circle *end_circle[4];
  start_circle[0] = new HC_CC_Circle(start1, true, true, true, rs_circle_param_);
  start_circle[1] = new HC_CC_Circle(start2, false, true, true, rs_circle_param_);
  start_circle[2] = new HC_CC_Circle(start1, true, false, true, rs_circle_param_);
  start_circle[3] = new HC_CC_Circle(start2, false, false, true, rs_circle_param_);
  end_circle[0] = new HC_CC_Circle(end, true, true, true, hc_cc_circle_param_);
  end_circle[1] = new HC_CC_Circle(end, false, true, true, hc_cc_circle_param_);
  end_circle[2] = new HC_CC_Circle(end, true, false, true, hc_cc_circle_param_);
  end_circle[3] = new HC_CC_Circle(end, false, false, true, hc_cc_circle_param_);

  // compute the shortest path for the 16 combinations (4 circles at the beginning and 4 at the end)
  HC_CC_RS_Path *path[] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };

  double lg[] = { numeric_limits<double>::max(), numeric_limits<double>::max(), numeric_limits<double>::max(),
                  numeric_limits<double>::max(), numeric_limits<double>::max(), numeric_limits<double>::max(),
                  numeric_limits<double>::max(), numeric_limits<double>::max(), numeric_limits<double>::max(),
                  numeric_limits<double>::max(), numeric_limits<double>::max(), numeric_limits<double>::max(),
                  numeric_limits<double>::max(), numeric_limits<double>::max(), numeric_limits<double>::max(),
                  numeric_limits<double>::max() };

  for (int i = 0; i < 4; i++)
  {
    // skip circle at the beginning for curvature continuity
    if (i == 0 && state1.kappa < 0)
      continue;
    else if (i == 1 && state1.kappa > 0)
      continue;
    else if (i == 2 && state1.kappa < 0)
      continue;
    else if (i == 3 && state1.kappa > 0)
      continue;
    for (int j = 0; j < 4; j++)
    {
      path[4 * i + j] = hcpm0_circles_rs_path(*start_circle[i], *end_circle[j]);
      if (path[4 * i + j])
      {
        lg[4 * i + j] = path[4 * i + j]->length;
      }
    }
  }

  // determine shortest path
  int best_path = array_index_min(lg, 16);

  //  // display calculations
  //  cout << "HCpm0_Reeds_Shepp_State_Space" << endl;
  //  for (int i = 0; i < 16; i++)
  //  {
  //    cout << i << ": ";
  //    if (path[i])
  //    {
  //      path[i]->print(true);
  //    }
  //  }
  //  cout << "shortest path: " << (int)best_path << endl;
  //  path[best_path]->print(true);

  // clean up
  for (int i = 0; i < 4; i++)
  {
    delete start_circle[i];
    delete end_circle[i];
  }
  for (int i = 0; i < 16; i++)
  {
    if (i != best_path)
    {
      delete path[i];
    }
  }
  return path[best_path];
}

double HCpm0_Reeds_Shepp_State_Space::get_distance(const State &state1, const State &state2) const
{
  HC_CC_RS_Path *p = this->hcpm0_reeds_shepp(state1, state2);
  double length = p->length;
  delete p;
  return length;
}

vector<Control> HCpm0_Reeds_Shepp_State_Space::get_controls(const State &state1, const State &state2) const
{
  vector<Control> hc_rs_controls;
  hc_rs_controls.reserve(9);
  HC_CC_RS_Path *p = this->hcpm0_reeds_shepp(state1, state2);
  switch (p->type)
  {
    case hc_cc_rs::E:
      empty_controls(hc_rs_controls);
      break;
    case hc_cc_rs::T:
      hc_turn_controls(*(p->cend), p->start, false, hc_rs_controls);
      break;
    case hc_cc_rs::TT:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, hc_rs_controls);
      cc_turn_controls(*(p->cend), *(p->qi2), false, hc_rs_controls);
      break;
    case hc_cc_rs::TcT:
      rs_turn_controls(*(p->cstart), *(p->qi1), true, hc_rs_controls);
      hc_turn_controls(*(p->cend), *(p->qi1), false, hc_rs_controls);
      break;
    case hc_cc_rs::TcTcT:
      rs_turn_controls(*(p->cstart), *(p->qi1), true, hc_rs_controls);
      rs_turn_controls(*(p->ci1), *(p->qi2), true, hc_rs_controls);
      hc_turn_controls(*(p->cend), *(p->qi2), false, hc_rs_controls);
      break;
    case hc_cc_rs::TcTT:
      rs_turn_controls(*(p->cstart), *(p->qi1), true, hc_rs_controls);
      hc_turn_controls(*(p->ci1), *(p->qi1), false, hc_rs_controls);
      cc_turn_controls(*(p->cend), *(p->qi2), false, hc_rs_controls);
      break;
    case hc_cc_rs::TTcT:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, hc_rs_controls);
      hc_turn_controls(*(p->ci1), *(p->qi2), true, hc_rs_controls);
      hc_turn_controls(*(p->cend), *(p->qi2), false, hc_rs_controls);
      break;
    case hc_cc_rs::TST:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, hc_rs_controls);
      straight_controls(*(p->qi2), *(p->qi3), hc_rs_controls);
      cc_turn_controls(*(p->cend), *(p->qi3), false, hc_rs_controls);
      break;
    case hc_cc_rs::TSTcT:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, hc_rs_controls);
      straight_controls(*(p->qi2), *(p->qi3), hc_rs_controls);
      hc_turn_controls(*(p->ci1), *(p->qi4), true, hc_rs_controls);
      hc_turn_controls(*(p->cend), *(p->qi4), false, hc_rs_controls);
      break;
    case hc_cc_rs::TcTST:
      rs_turn_controls(*(p->cstart), *(p->qi1), true, hc_rs_controls);
      hc_turn_controls(*(p->ci1), *(p->qi1), false, hc_rs_controls);
      straight_controls(*(p->qi2), *(p->qi3), hc_rs_controls);
      cc_turn_controls(*(p->cend), *(p->qi3), false, hc_rs_controls);
      break;
    case hc_cc_rs::TcTSTcT:
      rs_turn_controls(*(p->cstart), *(p->qi1), true, hc_rs_controls);
      hc_turn_controls(*(p->ci1), *(p->qi1), false, hc_rs_controls);
      straight_controls(*(p->qi2), *(p->qi3), hc_rs_controls);
      hc_turn_controls(*(p->ci2), *(p->qi4), true, hc_rs_controls);
      hc_turn_controls(*(p->cend), *(p->qi4), false, hc_rs_controls);
      break;
    case hc_cc_rs::TTcTT:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, hc_rs_controls);
      hc_turn_controls(*(p->ci1), *(p->qi2), true, hc_rs_controls);
      hc_turn_controls(*(p->ci2), *(p->qi2), false, hc_rs_controls);
      cc_turn_controls(*(p->cend), *(p->qi3), false, hc_rs_controls);
      break;
    case hc_cc_rs::TcTTcT:
      rs_turn_controls(*(p->cstart), *(p->qi1), true, hc_rs_controls);
      hc_turn_controls(*(p->ci1), *(p->qi1), false, hc_rs_controls);
      hc_turn_controls(*(p->ci2), *(p->qi2), true, hc_rs_controls);
      hc_turn_controls(*(p->cend), *(p->qi2), false, hc_rs_controls);
      break;
    // ########################################################################
    case hc_cc_rs::TTT:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, hc_rs_controls);
      cc_turn_controls(*(p->ci1), *(p->qi2), true, hc_rs_controls);
      cc_turn_controls(*(p->cend), *(p->qi2), false, hc_rs_controls);
      break;
    case hc_cc_rs::TcST:
    case hc_cc_rs::TScT:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, hc_rs_controls);
      straight_controls(*(p->qi2), *(p->qi3), hc_rs_controls);
      cc_turn_controls(*(p->cend), *(p->qi3), false, hc_rs_controls);
      break;
    case hc_cc_rs::TcScT:
      rs_turn_controls(*(p->cstart), *(p->qi1), true, hc_rs_controls);
      straight_controls(*(p->qi1), *(p->qi2), hc_rs_controls);
      hc_turn_controls(*(p->cend), *(p->qi2), false, hc_rs_controls);
      break;
    default:
      break;
  }
  delete p;
  return hc_rs_controls;
}
