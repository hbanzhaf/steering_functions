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

#include "steering_functions/hc_cc_state_space/ccpmpm_dubins_state_space.hpp"

class CCpmpm_Dubins_State_Space::CCpmpm_Dubins
{
private:
  CCpmpm_Dubins_State_Space *parent_;

public:
  explicit CCpmpm_Dubins(CCpmpm_Dubins_State_Space *parent)
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
    return fabs(distance - 2 * parent_->radius_) < get_epsilon();
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
        theta = angle + HALF_PI - parent_->mu_;
      }
      else
      {
        theta = angle + HALF_PI + parent_->mu_;
      }
    }
    else
    {
      if (c1.forward)
      {
        theta = angle - HALF_PI + parent_->mu_;
      }
      else
      {
        theta = angle - HALF_PI - parent_->mu_;
      }
    }
    *q = new Configuration(x, y, theta, 0);
  }

  double TT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                 Configuration **q1, Configuration **q2, Configuration **q3) const
  {
    TT_tangent_circles(c1, c2, q2);
    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, true, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(**q2, c2.left, !c2.forward, true, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    *q3 = new Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
    return (*cstart)->hc_turn_length(**q1) + (*cend)->hc_turn_length(**q3);
  }

  // ##### Dubins families: #####################################################

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
    return distance >= 2 * parent_->radius_;
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
    return distance >= 2 * parent_->radius_ * parent_->sin_mu_;
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
    double alpha = asin(2 * parent_->radius_ * parent_->cos_mu_ / distance);
    double delta_x = parent_->radius_ * parent_->sin_mu_;
    double delta_y = parent_->radius_ * parent_->cos_mu_;
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
    double delta_x = parent_->radius_ * parent_->sin_mu_;
    double delta_y = parent_->radius_ * parent_->cos_mu_;
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
                   Configuration **q1, Configuration **q2, Configuration **q3, Configuration **q4) const
  {
    TiST_tangent_circles(c1, c2, q2, q3);
    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, true, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(**q3, c2.left, !c2.forward, true, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    *q4 = new Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
    return (*cstart)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) + (*cend)->hc_turn_length(**q4);
  }

  double TeST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                   Configuration **q1, Configuration **q2, Configuration **q3, Configuration **q4) const
  {
    TeST_tangent_circles(c1, c2, q2, q3);
    *cstart = new HC_CC_Circle(**q2, c1.left, !c1.forward, true, parent_->hc_cc_circle_param_);
    *cend = new HC_CC_Circle(**q3, c2.left, !c2.forward, true, parent_->hc_cc_circle_param_);
    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    *q4 = new Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
    return (*cstart)->hc_turn_length(**q1) + configuration_distance(**q2, **q3) + (*cend)->hc_turn_length(**q4);
  }

  double TST_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                  Configuration **q1, Configuration **q2, Configuration **q3, Configuration **q4) const
  {
    if (TiST_exists(c1, c2))
    {
      return TiST_path(c1, c2, cstart, cend, q1, q2, q3, q4);
    }
    if (TeST_exists(c1, c2))
    {
      return TeST_path(c1, c2, cstart, cend, q1, q2, q3, q4);
    }
    return numeric_limits<double>::max();
  }

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
    return distance <= 4 * parent_->radius_;
  }

  void TTT_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2,
                           Configuration **q3, Configuration **q4) const
  {
    double theta = angle;
    double r = 2 * parent_->radius_;
    double delta_x = 0.5 * distance;
    double delta_y = sqrt(pow(r, 2) - pow(delta_x, 2));
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
                  Configuration **q1, Configuration **q2, Configuration **q3, HC_CC_Circle **ci) const
  {
    Configuration *qa, *qb, *qc, *qd;
    TTT_tangent_circles(c1, c2, &qa, &qb, &qc, &qd);
    HC_CC_Circle *start1, *start2, *end1, *end2, *middle1, *middle2;
    start1 = new HC_CC_Circle(*qa, c1.left, !c1.forward, true, parent_->hc_cc_circle_param_);
    middle1 = new HC_CC_Circle(*qa, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    end1 = new HC_CC_Circle(*qb, c2.left, !c2.forward, true, parent_->hc_cc_circle_param_);
    start2 = new HC_CC_Circle(*qc, c1.left, !c1.forward, true, parent_->hc_cc_circle_param_);
    middle2 = new HC_CC_Circle(*qc, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    end2 = new HC_CC_Circle(*qd, c2.left, !c2.forward, true, parent_->hc_cc_circle_param_);

    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    *q3 = new Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);

    // select shortest connection
    double length1 = start1->hc_turn_length(**q1) + middle1->cc_turn_length(*qb) + end1->hc_turn_length(**q3);
    double length2 = start2->hc_turn_length(**q1) + middle2->cc_turn_length(*qd) + end2->hc_turn_length(**q3);
    if (length1 < length2)
    {
      *cstart = start1;
      *ci = middle1;
      *cend = end1;
      *q2 = qb;
      delete qa;
      delete qc;
      delete qd;
      delete start2;
      delete middle2;
      delete end2;
      return length1;
    }
    else
    {
      *cstart = start2;
      *ci = middle2;
      *cend = end2;
      *q2 = qd;
      delete qa;
      delete qb;
      delete qc;
      delete start1;
      delete middle1;
      delete end1;
      return length2;
    }
    return numeric_limits<double>::max();
  }

  // ############################################################################

  // ##### TTTT #################################################################
  bool TTTT_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2) const
  {
    if (c1.left == c2.left)
    {
      return false;
    }
    if (c1.forward == c2.forward)
    {
      return false;
    }
    return distance <= 6 * parent_->radius_;
  }

  void TTTT_tangent_circles(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2,
                            Configuration **q3, Configuration **q4, Configuration **q5, Configuration **q6) const
  {
    double theta = angle;
    double r1, delta_x, delta_y, x, y;
    r1 = 2 * parent_->radius_;
    if (distance < r1)
    {
      delta_x = (distance + r1) / 2;
      delta_y = sqrt(pow(r1, 2) - pow(delta_x, 2));
    }
    else
    {
      delta_x = (distance - r1) / 2;
      delta_y = sqrt(pow(r1, 2) - pow(delta_x, 2));
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
    TT_tangent_circles(tgt1, tgt2, q2);
    TT_tangent_circles(tgt2, c2, q3);

    TT_tangent_circles(c1, tgt3, q4);
    TT_tangent_circles(tgt3, tgt4, q5);
    TT_tangent_circles(tgt4, c2, q6);
  }

  double TTTT_path(const HC_CC_Circle &c1, const HC_CC_Circle &c2, HC_CC_Circle **cstart, HC_CC_Circle **cend,
                   Configuration **q1, Configuration **q2, Configuration **q3, HC_CC_Circle **ci1,
                   HC_CC_Circle **ci2) const
  {
    Configuration *qa, *qb, *qc, *qd, *qe, *qf;
    TTTT_tangent_circles(c1, c2, &qa, &qb, &qc, &qd, &qe, &qf);
    HC_CC_Circle *start1, *start2, *end1, *end2, *middle1, *middle2, *middle3, *middle4;
    start1 = new HC_CC_Circle(*qa, c1.left, !c1.forward, true, parent_->hc_cc_circle_param_);
    middle1 = new HC_CC_Circle(*qa, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    middle2 = new HC_CC_Circle(*qc, !c2.left, c2.forward, true, parent_->hc_cc_circle_param_);
    end1 = new HC_CC_Circle(*qc, c2.left, !c2.forward, true, parent_->hc_cc_circle_param_);
    start2 = new HC_CC_Circle(*qd, c1.left, !c1.forward, true, parent_->hc_cc_circle_param_);
    middle3 = new HC_CC_Circle(*qd, !c1.left, c1.forward, true, parent_->hc_cc_circle_param_);
    middle4 = new HC_CC_Circle(*qf, !c2.left, c2.forward, true, parent_->hc_cc_circle_param_);
    end2 = new HC_CC_Circle(*qf, c2.left, !c2.forward, true, parent_->hc_cc_circle_param_);

    *q1 = new Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
    *q3 = new Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);

    // select shortest connection
    double length1 = start1->hc_turn_length(**q1) + middle1->cc_turn_length(*qb) + middle2->cc_turn_length(*qb) +
                     end1->hc_turn_length(**q3);
    double length2 = start2->hc_turn_length(**q1) + middle3->cc_turn_length(*qe) + middle4->cc_turn_length(*qe) +
                     end2->hc_turn_length(**q3);
    if (length1 < length2)
    {
      *cstart = start1;
      *cend = end1;
      *ci1 = middle1;
      *ci2 = middle2;
      *q2 = qb;
      delete qa;
      delete qc;
      delete qd;
      delete qe;
      delete qf;
      delete start2, delete end2, delete middle3;
      delete middle4;
      return length1;
    }
    else
    {
      *cstart = start2;
      *cend = end2;
      *ci1 = middle3;
      *ci2 = middle4;
      *q2 = qe;
      delete qa;
      delete qb;
      delete qc;
      delete qd;
      delete qf;
      delete start1;
      delete end1;
      delete middle1;
      delete middle2;
      return length2;
    }
    return numeric_limits<double>::max();
  }
};

// ############################################################################

CCpmpm_Dubins_State_Space::CCpmpm_Dubins_State_Space(double kappa, double sigma, double discretization, bool forwards)
  : HC_CC_State_Space(kappa, sigma, discretization)
  , forwards_(forwards)
  , ccpmpm_dubins_{ unique_ptr<CCpmpm_Dubins>(new CCpmpm_Dubins(this)) }
{
  rs_circle_param_.set_param(kappa_, numeric_limits<double>::max(), 1 / kappa_, 0.0, 0.0, 1.0, 0.0);
  radius_ = hc_cc_circle_param_.radius;
  mu_ = hc_cc_circle_param_.mu;
  sin_mu_ = hc_cc_circle_param_.sin_mu;
  cos_mu_ = hc_cc_circle_param_.cos_mu;
}

CCpmpm_Dubins_State_Space::~CCpmpm_Dubins_State_Space() = default;

CC_Dubins_Path *CCpmpm_Dubins_State_Space::ccpmpm_circles_dubins_path(const HC_CC_Circle &c1,
                                                                      const HC_CC_Circle &c2) const
{
  // table containing the lengths of the paths, the intermediate configurations and circles
  double length[nb_cc_dubins_paths];
  double_array_init(length, nb_cc_dubins_paths, numeric_limits<double>::max());
  Configuration *qi1[nb_cc_dubins_paths];
  pointer_array_init((void **)qi1, nb_cc_dubins_paths);
  Configuration *qi2[nb_cc_dubins_paths];
  pointer_array_init((void **)qi2, nb_cc_dubins_paths);
  Configuration *qi3[nb_cc_dubins_paths];
  pointer_array_init((void **)qi3, nb_cc_dubins_paths);
  Configuration *qi4[nb_cc_dubins_paths];
  pointer_array_init((void **)qi4, nb_cc_dubins_paths);
  HC_CC_Circle *cstart[nb_cc_dubins_paths];
  pointer_array_init((void **)cstart, nb_cc_dubins_paths);
  HC_CC_Circle *ci1[nb_cc_dubins_paths];
  pointer_array_init((void **)ci1, nb_cc_dubins_paths);
  HC_CC_Circle *ci2[nb_cc_dubins_paths];
  pointer_array_init((void **)ci2, nb_cc_dubins_paths);
  HC_CC_Circle *cend[nb_cc_dubins_paths];
  pointer_array_init((void **)cend, nb_cc_dubins_paths);

  // precomputations
  ccpmpm_dubins_->distance = center_distance(c1, c2);
  ccpmpm_dubins_->angle = atan2(c2.yc - c1.yc, c2.xc - c1.xc);

  // case E
  if (configuration_equal(c1.start, c2.start))
  {
    length[cc_dubins::E] = 0;
    goto label_end;
  }
  // case T
  if (configuration_on_hc_cc_circle(c1, c2.start))
  {
    cstart[cc_dubins::T] = new HC_CC_Circle(c1.start, c1.left, c1.forward, false, rs_circle_param_);
    length[cc_dubins::T] = cstart[cc_dubins::T]->rs_turn_length(c2.start);
    goto label_end;
  }
  // case TT
  if (ccpmpm_dubins_->TT_exists(c1, c2))
  {
    length[cc_dubins::TT] = ccpmpm_dubins_->TT_path(c1, c2, &cstart[cc_dubins::TT], &cend[cc_dubins::TT],
                                                    &qi1[cc_dubins::TT], &qi2[cc_dubins::TT], &qi3[cc_dubins::TT]);
  }
  // ##### Dubins families: ######################################################
  // case TST
  if (ccpmpm_dubins_->TST_exists(c1, c2))
  {
    length[cc_dubins::TST] =
        ccpmpm_dubins_->TST_path(c1, c2, &cstart[cc_dubins::TST], &cend[cc_dubins::TST], &qi1[cc_dubins::TST],
                                 &qi2[cc_dubins::TST], &qi3[cc_dubins::TST], &qi4[cc_dubins::TST]);
  }
  // case TTT
  if (ccpmpm_dubins_->TTT_exists(c1, c2))
  {
    length[cc_dubins::TTT] =
        ccpmpm_dubins_->TTT_path(c1, c2, &cstart[cc_dubins::TTT], &cend[cc_dubins::TTT], &qi1[cc_dubins::TTT],
                                 &qi2[cc_dubins::TTT], &qi3[cc_dubins::TTT], &ci1[cc_dubins::TTT]);
  }
  //  ############################################################################
  // case TTTT
  if (ccpmpm_dubins_->TTTT_exists(c1, c2))
  {
    length[cc_dubins::TTTT] = ccpmpm_dubins_->TTTT_path(
        c1, c2, &cstart[cc_dubins::TTTT], &cend[cc_dubins::TTTT], &qi1[cc_dubins::TTTT], &qi2[cc_dubins::TTTT],
        &qi3[cc_dubins::TTTT], &ci1[cc_dubins::TTTT], &ci2[cc_dubins::TTTT]);
  }
label_end:
  // select shortest path
  cc_dubins::path_type best_path = (cc_dubins::path_type)array_index_min(length, nb_cc_dubins_paths);
  CC_Dubins_Path *path;
  path = new CC_Dubins_Path(c1.start, c2.start, best_path, kappa_, sigma_, qi1[best_path], qi2[best_path],
                            qi3[best_path], qi4[best_path], cstart[best_path], cend[best_path], ci1[best_path],
                            ci2[best_path], length[best_path]);

  // clean up
  for (int i = 0; i < nb_cc_dubins_paths; i++)
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

CC_Dubins_Path *CCpmpm_Dubins_State_Space::ccpmpm_dubins(const State &state1, const State &state2) const
{
  // compute the 2 circles at the intial and final configuration
  Configuration start1(state1.x, state1.y, state1.theta, kappa_);
  Configuration start2(state1.x, state1.y, state1.theta, -kappa_);
  Configuration end1(state2.x, state2.y, state2.theta, kappa_);
  Configuration end2(state2.x, state2.y, state2.theta, -kappa_);

  HC_CC_Circle *start_circle[2];
  HC_CC_Circle *end_circle[2];
  if (forwards_)
  {
    start_circle[0] = new HC_CC_Circle(start1, true, true, true, rs_circle_param_);
    start_circle[1] = new HC_CC_Circle(start2, false, true, true, rs_circle_param_);
    end_circle[0] = new HC_CC_Circle(end1, true, false, true, rs_circle_param_);
    end_circle[1] = new HC_CC_Circle(end2, false, false, true, rs_circle_param_);
  }
  else
  {
    start_circle[0] = new HC_CC_Circle(start1, true, false, true, rs_circle_param_);
    start_circle[1] = new HC_CC_Circle(start2, false, false, true, rs_circle_param_);
    end_circle[0] = new HC_CC_Circle(end1, true, true, true, rs_circle_param_);
    end_circle[1] = new HC_CC_Circle(end2, false, true, true, rs_circle_param_);
  }

  // compute the shortest path for the 4 combinations (2 circles at the beginning and 2 at the end)
  CC_Dubins_Path *path[] = { nullptr, nullptr, nullptr, nullptr };

  double lg[] = { numeric_limits<double>::max(), numeric_limits<double>::max(), numeric_limits<double>::max(),
                  numeric_limits<double>::max() };

  for (int i = 0; i < 2; i++)
  {
    // skip circle at the beginning for curvature continuity
    if (i == 0 && state1.kappa < 0)
      continue;
    else if (i == 1 && state1.kappa > 0)
      continue;
    for (int j = 0; j < 2; j++)
    {
      // skip circle at the end for curvature continuity
      if (j == 0 && state2.kappa < 0)
        continue;
      else if (j == 1 && state2.kappa > 0)
        continue;
      path[2 * i + j] = ccpmpm_circles_dubins_path(*start_circle[i], *end_circle[j]);
      if (path[2 * i + j])
      {
        lg[2 * i + j] = path[2 * i + j]->length;
      }
    }
  }

  // select shortest path
  int best_path = array_index_min(lg, 4);

  //  // display calculations
  //  cout << endl << "CCpmpm_Dubins_State_Space" << endl;
  //  for (int i = 0; i < 4; i++)
  //  {
  //    cout << i << ": ";
  //    if (path[i])
  //    {
  //      path[i]->print(true);
  //      cout << endl;
  //    }
  //  }
  //  cout << "shortest path: " << (int)best_path << endl;
  //  path[best_path]->print(true);

  // clean up
  for (int i = 0; i < 2; i++)
  {
    delete start_circle[i];
    delete end_circle[i];
  }
  for (int i = 0; i < 4; i++)
  {
    if (i != best_path)
    {
      delete path[i];
    }
  }
  return path[best_path];
}

double CCpmpm_Dubins_State_Space::get_distance(const State &state1, const State &state2) const
{
  CC_Dubins_Path *p = this->ccpmpm_dubins(state1, state2);
  double length = p->length;
  delete p;
  return length;
}

vector<Control> CCpmpm_Dubins_State_Space::get_controls(const State &state1, const State &state2) const
{
  vector<Control> cc_dubins_controls;
  cc_dubins_controls.reserve(10);
  CC_Dubins_Path *p = this->ccpmpm_dubins(state1, state2);
  switch (p->type)
  {
    case cc_dubins::E:
      empty_controls(cc_dubins_controls);
      break;
    case cc_dubins::T:
      rs_turn_controls(*(p->cstart), p->end, true, cc_dubins_controls);
      break;
    case cc_dubins::TT:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, cc_dubins_controls);
      hc_turn_controls(*(p->cend), *(p->qi3), true, cc_dubins_controls);
      break;
    // ##### Dubins families: #####################################################
    case cc_dubins::TST:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, cc_dubins_controls);
      straight_controls(*(p->qi2), *(p->qi3), cc_dubins_controls);
      hc_turn_controls(*(p->cend), *(p->qi4), true, cc_dubins_controls);
      break;
    case cc_dubins::TTT:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, cc_dubins_controls);
      cc_turn_controls(*(p->ci1), *(p->qi2), true, cc_dubins_controls);
      hc_turn_controls(*(p->cend), *(p->qi3), true, cc_dubins_controls);
      break;
    // ############################################################################
    case cc_dubins::TTTT:
      hc_turn_controls(*(p->cstart), *(p->qi1), false, cc_dubins_controls);
      cc_turn_controls(*(p->ci1), *(p->qi2), true, cc_dubins_controls);
      cc_turn_controls(*(p->ci2), *(p->qi2), false, cc_dubins_controls);
      hc_turn_controls(*(p->cend), *(p->qi3), true, cc_dubins_controls);
      break;
    default:
      break;
  }
  delete p;
  return cc_dubins_controls;
}
