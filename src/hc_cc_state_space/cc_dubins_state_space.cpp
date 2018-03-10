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

#include "steering_functions/hc_cc_state_space/cc_dubins_state_space.hpp"

namespace cc_dubins
{
bool external_mu_tangent_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2)
{
  if (fabs(c1.radius - c2.radius) > get_epsilon())
  {
    return false;
  }
  if (fabs(c1.mu - c2.mu) > get_epsilon())
  {
    return false;
  }
  if (c1.left != c2.left)
  {
    return false;
  }
  if (c1.forward == c2.forward)
  {
    return false;
  }
  double distance = point_distance(c1.xc, c1.yc, c2.xc, c2.yc);
  return (distance >= 2 * c1.radius * c1.sin_mu);
}

void external_mu_tangent(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2)
{
  double theta = atan2(c2.yc - c1.yc, c2.xc - c1.xc);
  double delta_x = fabs(c1.radius * c1.sin_mu);
  double delta_y = fabs(c1.radius * c1.cos_mu);
  double x, y;
  if (c1.left && c1.forward)
  {
    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new Configuration(x, y, theta, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new Configuration(x, y, theta, 0);
  }
  if (!c1.left && c1.forward)
  {
    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new Configuration(x, y, theta, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new Configuration(x, y, theta, 0);
  }
}

bool internal_mu_tangent_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2)
{
  if (c1.left == c2.left)
  {
    return false;
  }
  if (c1.forward == c2.forward)
  {
    return false;
  }
  double distance = point_distance(c1.xc, c1.yc, c2.xc, c2.yc);
  return (distance >= 2 * c1.radius);
}

void internal_mu_tangent(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2)
{
  double distance = point_distance(c1.xc, c1.yc, c2.xc, c2.yc);
  double theta = atan2(c2.yc - c1.yc, c2.xc - c1.xc);
  double alpha = fabs(asin(2 * c1.radius * c1.cos_mu / distance));
  double delta_x = fabs(c1.radius * c1.sin_mu);
  double delta_y = fabs(c1.radius * c1.cos_mu);
  double x, y;
  if (c1.left && c1.forward)
  {
    theta = theta + alpha;
    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new Configuration(x, y, theta, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new Configuration(x, y, theta, 0);
  }
  if (!c1.left && c1.forward)
  {
    theta = theta - alpha;
    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new Configuration(x, y, theta, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new Configuration(x, y, theta, 0);
  }
}

bool tangent_circle_exists(const HC_CC_Circle &c1, const HC_CC_Circle &c2)
{
  if (c1.left != c2.left)
  {
    return false;
  }
  if (c1.forward == c2.forward)
  {
    return false;
  }
  double distance = point_distance(c1.xc, c1.yc, c2.xc, c2.yc);
  return (distance <= 4 * c1.radius);
}

void tangent_circle(const HC_CC_Circle &c1, const HC_CC_Circle &c2, Configuration **q1, Configuration **q2,
                    Configuration **q3, Configuration **q4)
{
  double distance = point_distance(c1.xc, c1.yc, c2.xc, c2.yc);
  double theta = atan2(c2.yc - c1.yc, c2.xc - c1.xc);
  double h = sqrt(fabs(pow(2 * c1.radius, 2) - pow(0.5 * distance, 2)));
  double alpha = fabs(atan(2 * h / distance));
  double delta_x = fabs(c1.radius * cos(alpha));
  double delta_y = fabs(c1.radius * sin(alpha));
  double x, y;
  if (c1.left && c1.forward)
  {
    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new Configuration(x, y, theta + alpha + HALF_PI - c1.mu, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new Configuration(x, y, theta - alpha + 1.5 * PI + c1.mu, 0);
    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    *q3 = new Configuration(x, y, theta - alpha + HALF_PI - c1.mu, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
    *q4 = new Configuration(x, y, theta + alpha + 1.5 * PI + c1.mu, 0);
  }
  if (c1.left && !c1.forward)
  {
    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new Configuration(x, y, theta - alpha + HALF_PI + c1.mu, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new Configuration(x, y, theta + alpha + 1.5 * PI - c1.mu, 0);
    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    *q3 = new Configuration(x, y, theta + alpha + HALF_PI + c1.mu, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    *q4 = new Configuration(x, y, theta - alpha + 1.5 * PI - c1.mu, 0);
  }
  if (!c1.left && c1.forward)
  {
    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new Configuration(x, y, theta - alpha - HALF_PI + c1.mu, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new Configuration(x, y, theta + alpha + HALF_PI - c1.mu, 0);
    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    *q3 = new Configuration(x, y, theta + alpha - HALF_PI + c1.mu, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    *q4 = new Configuration(x, y, theta - alpha + HALF_PI - c1.mu, 0);
  }
  if (!c1.left && !c1.forward)
  {
    global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new Configuration(x, y, theta + alpha - HALF_PI - c1.mu, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new Configuration(x, y, theta - alpha + HALF_PI + c1.mu, 0);
    global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
    *q3 = new Configuration(x, y, theta - alpha - HALF_PI - c1.mu, 0);
    global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
    *q4 = new Configuration(x, y, theta + alpha - HALF_PI + c1.mu, 0);
  }
}
}

CC_Dubins_Path *CC_Dubins_State_Space::cc_dubins(const State &state1, const State &state2) const
{
  // table containing the lengths of the paths, the intermediate configurations and circles
  double length[nb_cc_dubins_paths];
  double_array_init(length, nb_cc_dubins_paths, numeric_limits<double>::max());
  Configuration *qi1[nb_cc_dubins_paths];
  pointer_array_init((void **)qi1, nb_cc_dubins_paths);
  Configuration *qi2[nb_cc_dubins_paths];
  pointer_array_init((void **)qi2, nb_cc_dubins_paths);
  HC_CC_Circle *cstart[nb_cc_dubins_paths];
  pointer_array_init((void **)cstart, nb_cc_dubins_paths);
  HC_CC_Circle *ci1[nb_cc_dubins_paths];
  pointer_array_init((void **)ci1, nb_cc_dubins_paths);
  HC_CC_Circle *cend[nb_cc_dubins_paths];
  pointer_array_init((void **)cend, nb_cc_dubins_paths);

  // compute the 2 circles at the intial and final configuration
  Configuration start(state1.x, state1.y, state1.theta, 0.0);
  Configuration end(state2.x, state2.y, state2.theta, 0.0);

  HC_CC_Circle *start_left_forward, *start_right_forward;
  HC_CC_Circle *end_left_backward, *end_right_backward;
  start_left_forward = new HC_CC_Circle(start, true, true, true, hc_cc_circle_param_);
  start_right_forward = new HC_CC_Circle(start, false, true, true, hc_cc_circle_param_);
  end_left_backward = new HC_CC_Circle(end, true, false, true, hc_cc_circle_param_);
  end_right_backward = new HC_CC_Circle(end, false, false, true, hc_cc_circle_param_);

  // case Empty
  if (configuration_equal(start, end))
  {
    length[E] = 0;
    goto label_end;
  }
  // case Straight
  if (configuration_aligned(start, end))
  {
    length[S] = configuration_distance(start, end);
    goto label_end;
  }
  // case Left
  if (configuration_on_hc_cc_circle(*start_left_forward, end))
  {
    cstart[L] = new HC_CC_Circle(*start_left_forward);
    length[L] = start_left_forward->cc_turn_length(end);
    goto label_end;
  }
  // case Right
  if (configuration_on_hc_cc_circle(*start_right_forward, end))
  {
    cstart[R] = new HC_CC_Circle(*start_right_forward);
    length[R] = start_right_forward->cc_turn_length(end);
    goto label_end;
  }
  // case LSL and subcases LeS, eSL
  if (cc_dubins::external_mu_tangent_exists(*start_left_forward, *end_left_backward))
  {
    Configuration *qa, *qb;
    cc_dubins::external_mu_tangent(*start_left_forward, *end_left_backward, &qa, &qb);
    // subcase LeS
    if (configuration_aligned(*qb, end))
    {
      cstart[LeS] = new HC_CC_Circle(*start_left_forward);
      qi1[LeS] = new Configuration(*qa);
      length[LeS] = start_left_forward->cc_turn_length(*qa) + configuration_distance(*qa, end);
    }
    // subcase eSL
    else if (configuration_aligned(start, *qa))
    {
      cend[eSL] = new HC_CC_Circle(*end_left_backward);
      qi1[eSL] = new Configuration(*qb);
      length[eSL] = configuration_distance(start, *qb) + end_left_backward->cc_turn_length(*qb);
    }
    // case LSL
    else
    {
      cstart[LSL] = new HC_CC_Circle(*start_left_forward);
      cend[LSL] = new HC_CC_Circle(*end_left_backward);
      qi1[LSL] = new Configuration(*qa);
      qi2[LSL] = new Configuration(*qb);
      length[LSL] = start_left_forward->cc_turn_length(*qa) + configuration_distance(*qa, *qb) +
                    end_left_backward->cc_turn_length(*qb);
    }
    delete qa;
    delete qb;
  }
  // case LSR and subcases LiS, iSR
  if (cc_dubins::internal_mu_tangent_exists(*start_left_forward, *end_right_backward))
  {
    Configuration *qa, *qb;
    cc_dubins::internal_mu_tangent(*start_left_forward, *end_right_backward, &qa, &qb);
    // subcase LiS
    if (configuration_aligned(*qb, end))
    {
      cstart[LiS] = new HC_CC_Circle(*start_left_forward);
      qi1[LiS] = new Configuration(*qa);
      length[LiS] = start_left_forward->cc_turn_length(*qa) + configuration_distance(*qa, end);
    }
    // subcase iSR
    else if (configuration_aligned(start, *qa))
    {
      cend[iSR] = new HC_CC_Circle(*end_right_backward);
      qi1[iSR] = new Configuration(*qb);
      length[iSR] = configuration_distance(start, *qb) + end_right_backward->cc_turn_length(*qb);
    }
    // case LSR
    else
    {
      cstart[LSR] = new HC_CC_Circle(*start_left_forward);
      cend[LSR] = new HC_CC_Circle(*end_right_backward);
      qi1[LSR] = new Configuration(*qa);
      qi2[LSR] = new Configuration(*qb);
      length[LSR] = start_left_forward->cc_turn_length(*qa) + configuration_distance(*qa, *qb) +
                    end_right_backward->cc_turn_length(*qb);
    }
    delete qa;
    delete qb;
  }
  // case RSL and subcases RiS, iSL
  if (cc_dubins::internal_mu_tangent_exists(*start_right_forward, *end_left_backward))
  {
    Configuration *qa, *qb;
    cc_dubins::internal_mu_tangent(*start_right_forward, *end_left_backward, &qa, &qb);
    // subcase RiS
    if (configuration_aligned(*qb, end))
    {
      cstart[RiS] = new HC_CC_Circle(*start_right_forward);
      qi1[RiS] = new Configuration(*qa);
      length[RiS] = start_right_forward->cc_turn_length(*qa) + configuration_distance(*qa, end);
    }
    // subcase iSL
    else if (configuration_aligned(start, *qa))
    {
      cend[iSL] = new HC_CC_Circle(*end_left_backward);
      qi1[iSL] = new Configuration(*qb);
      length[iSL] = configuration_distance(start, *qb) + end_left_backward->cc_turn_length(*qb);
    }
    // case RSL
    else
    {
      cstart[RSL] = new HC_CC_Circle(*start_right_forward);
      cend[RSL] = new HC_CC_Circle(*end_left_backward);
      qi1[RSL] = new Configuration(*qa);
      qi2[RSL] = new Configuration(*qb);
      length[RSL] = start_right_forward->cc_turn_length(*qa) + configuration_distance(*qa, *qb) +
                    end_left_backward->cc_turn_length(*qb);
    }
    delete qa;
    delete qb;
  }
  // case RSR and subcases ReS, eSR
  if (cc_dubins::external_mu_tangent_exists(*start_right_forward, *end_right_backward))
  {
    Configuration *qa, *qb;
    cc_dubins::external_mu_tangent(*start_right_forward, *end_right_backward, &qa, &qb);
    // subcase ReS
    if (configuration_aligned(*qb, end))
    {
      cstart[ReS] = new HC_CC_Circle(*start_right_forward);
      qi1[ReS] = new Configuration(*qa);
      length[ReS] = start_right_forward->cc_turn_length(*qa) + configuration_distance(*qa, end);
    }
    // subcase eSR
    else if (configuration_aligned(start, *qa))
    {
      cend[eSR] = new HC_CC_Circle(*end_right_backward);
      qi1[eSR] = new Configuration(*qb);
      length[eSR] = configuration_distance(start, *qb) + end_right_backward->cc_turn_length(*qb);
    }
    // case RSR
    else
    {
      cstart[RSR] = new HC_CC_Circle(*start_right_forward);
      cend[RSR] = new HC_CC_Circle(*end_right_backward);
      qi1[RSR] = new Configuration(*qa);
      qi2[RSR] = new Configuration(*qb);
      length[RSR] = start_right_forward->cc_turn_length(*qa) + configuration_distance(*qa, *qb) +
                    end_right_backward->cc_turn_length(*qb);
    }
    delete qa;
    delete qb;
  }
  // case LRL
  if (cc_dubins::tangent_circle_exists(*start_left_forward, *end_left_backward))
  {
    HC_CC_Circle *middle_right_forward = nullptr;
    cc_dubins::tangent_circle(*start_left_forward, *end_left_backward, &qi1[LR1L], &qi2[LR1L], &qi1[LR2L], &qi2[LR2L]);
    cstart[LR1L] = new HC_CC_Circle(*start_left_forward);
    cend[LR1L] = new HC_CC_Circle(*end_left_backward);
    middle_right_forward = new HC_CC_Circle(*qi1[LR1L], false, true, true, hc_cc_circle_param_);
    ci1[LR1L] = middle_right_forward;
    length[LR1L] = start_left_forward->cc_turn_length(*qi1[LR1L]) + middle_right_forward->cc_turn_length(*qi2[LR1L]) +
                   end_left_backward->cc_turn_length(*qi2[LR1L]);

    cstart[LR2L] = new HC_CC_Circle(*start_left_forward);
    cend[LR2L] = new HC_CC_Circle(*end_left_backward);
    middle_right_forward = new HC_CC_Circle(*qi1[LR2L], false, true, true, hc_cc_circle_param_);
    ci1[LR2L] = middle_right_forward;
    length[LR2L] = start_left_forward->cc_turn_length(*qi1[LR2L]) + middle_right_forward->cc_turn_length(*qi2[LR2L]) +
                   end_left_backward->cc_turn_length(*qi2[LR2L]);
  }
  // case RLR
  if (cc_dubins::tangent_circle_exists(*start_right_forward, *end_right_backward))
  {
    HC_CC_Circle *middle_left_forward = nullptr;
    cc_dubins::tangent_circle(*start_right_forward, *end_right_backward, &qi1[RL1R], &qi2[RL1R], &qi1[RL2R],
                              &qi2[RL2R]);
    cstart[RL1R] = new HC_CC_Circle(*start_right_forward);
    cend[RL1R] = new HC_CC_Circle(*end_right_backward);
    middle_left_forward = new HC_CC_Circle(*qi1[RL1R], true, true, true, hc_cc_circle_param_);
    ci1[RL1R] = middle_left_forward;
    length[RL1R] = start_right_forward->cc_turn_length(*qi1[RL1R]) + middle_left_forward->cc_turn_length(*qi2[RL1R]) +
                   end_right_backward->cc_turn_length(*qi2[RL1R]);

    cstart[RL2R] = new HC_CC_Circle(*start_right_forward);
    cend[RL2R] = new HC_CC_Circle(*end_right_backward);
    middle_left_forward = new HC_CC_Circle(*qi1[RL2R], true, true, true, hc_cc_circle_param_);
    ci1[RL2R] = middle_left_forward;
    length[RL2R] = start_right_forward->cc_turn_length(*qi1[RL2R]) + middle_left_forward->cc_turn_length(*qi2[RL2R]) +
                   end_right_backward->cc_turn_length(*qi2[RL2R]);
  }
label_end:
  // select shortest path
  cc_dubins_path_type best_path = (cc_dubins_path_type)array_index_min(length, nb_cc_dubins_paths);
  CC_Dubins_Path *path;
  path = new CC_Dubins_Path(start, end, best_path, kappa_, sigma_, qi1[best_path], qi2[best_path], cstart[best_path],
                            cend[best_path], ci1[best_path], length[best_path]);

  //  // display calculations
  //  cout << endl << "CC_Dubins_State_Space" << endl;
  //  for (int i = 0; i < nb_cc_dubins_paths; i++)
  //  {
  //    cout << i << ": ";
  //    if (qi1[i])
  //    {
  //      qi1[i]->print(false);
  //    }
  //    cout << ", ";
  //    if (qi2[i])
  //    {
  //      qi2[i]->print(false);
  //    }
  //    cout << ", " << length[i] << endl;
  //  }
  //  cout << "shortest path: " << (int)best_path << endl;
  //  path->print(true);

  // clean up
  delete start_left_forward;
  delete start_right_forward;
  delete end_left_backward;
  delete end_right_backward;
  for (int i = 0; i < nb_cc_dubins_paths; i++)
  {
    if (i != best_path)
    {
      delete qi1[i];
      delete qi2[i];
      delete cstart[i];
      delete ci1[i];
      delete cend[i];
    }
  }
  return path;
}

double CC_Dubins_State_Space::get_distance(const State &state1, const State &state2) const
{
  CC_Dubins_Path *p;
  if (forwards_)
    p = this->cc_dubins(state1, state2);
  else
    p = this->cc_dubins(state2, state1);
  double length = p->length;
  delete p;
  return length;
}

vector<Control> CC_Dubins_State_Space::get_controls(const State &state1, const State &state2) const
{
  vector<Control> cc_dubins_controls;
  cc_dubins_controls.reserve(9);
  CC_Dubins_Path *p;
  if (forwards_)
    p = this->cc_dubins(state1, state2);
  else
    p = this->cc_dubins(state2, state1);
  switch (p->type)
  {
    case E:
      empty_controls(cc_dubins_controls);
      break;
    case S:
      straight_controls(p->start, p->end, cc_dubins_controls);
      break;
    case R:
    case L:
      cc_turn_controls(*(p->cstart), p->end, true, cc_dubins_controls);
      break;
    case eSL:
    case iSL:
    case eSR:
    case iSR:
      straight_controls(p->start, *(p->qi1), cc_dubins_controls);
      cc_turn_controls(*(p->cend), *(p->qi1), false, cc_dubins_controls);
      break;
    case ReS:
    case RiS:
    case LeS:
    case LiS:
      cc_turn_controls(*(p->cstart), *(p->qi1), true, cc_dubins_controls);
      straight_controls(*(p->qi1), p->end, cc_dubins_controls);
      break;
    case LSL:
    case LSR:
    case RSL:
    case RSR:
      cc_turn_controls(*(p->cstart), *(p->qi1), true, cc_dubins_controls);
      straight_controls(*(p->qi1), *(p->qi2), cc_dubins_controls);
      cc_turn_controls(*(p->cend), *(p->qi2), false, cc_dubins_controls);
      break;
    case LR1L:
    case LR2L:
    case RL1R:
    case RL2R:
      cc_turn_controls(*(p->cstart), *(p->qi1), true, cc_dubins_controls);
      cc_turn_controls(*(p->ci1), *(p->qi2), true, cc_dubins_controls);
      cc_turn_controls(*(p->cend), *(p->qi2), false, cc_dubins_controls);
      break;
    default:
      break;
  }
  // reverse controls
  if (!forwards_)
  {
    reverse(cc_dubins_controls.begin(), cc_dubins_controls.end());
    for (auto &control : cc_dubins_controls)
    {
      reverse_control(control);
    }
  }
  delete p;
  return cc_dubins_controls;
}
