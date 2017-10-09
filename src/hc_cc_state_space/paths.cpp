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

#include "steering_functions/hc_cc_state_space/paths.hpp"

Path::Path(const Configuration &_start, const Configuration &_end, double _kappa, double _sigma, double _length)
{
  start = _start;
  end = _end;
  kappa = _kappa;
  sigma = _sigma;
  length = _length;
}

CC_Dubins_Path::CC_Dubins_Path(const Configuration &_start, const Configuration &_end, cc_dubins_path_type _type,
                               double _kappa, double _sigma, Configuration *_qi1, Configuration *_qi2,
                               HC_CC_Circle *_cstart, HC_CC_Circle *_cend, HC_CC_Circle *_ci1, double _length)
  : Path(_start, _end, _kappa, _sigma, _length)
{
  type = _type;
  qi1 = _qi1;
  qi2 = _qi2;
  cstart = _cstart;
  cend = _cend;
  ci1 = _ci1;
  ci2 = nullptr;
}

CC_Dubins_Path::~CC_Dubins_Path()
{
  delete qi1;
  delete qi2;
  delete cstart;
  delete cend;
  delete ci1;
  delete ci2;
}

void CC_Dubins_Path::print(bool eol) const
{
  cout << "CC_Dubins_Path: type ";
  switch (type)
  {
    case E:
      cout << "EMPTY";
      break;
    case LSL:
      cout << "LSL";
      break;
    case LSR:
      cout << "LSR";
      break;
    case RSL:
      cout << "RSL";
      break;
    case RSR:
      cout << "RSR";
      break;
    case LR1L:
      cout << "LR1L";
      break;
    case LR2L:
      cout << "LR2L";
      break;
    case RL1R:
      cout << "RL1R";
      break;
    case RL2R:
      cout << "RL2R";
      break;
    case S:
      cout << "S";
      break;
    case L:
      cout << "L";
      break;
    case R:
      cout << "R";
      break;
    case LeS:
      cout << "LeS";
      break;
    case LiS:
      cout << "LiS";
      break;
    case eSL:
      cout << "eSL";
      break;
    case iSL:
      cout << "iSL";
      break;
    case eSR:
      cout << "eSR";
      break;
    case iSR:
      cout << "iSR";
      break;
    case ReS:
      cout << "ReS";
      break;
    case RiS:
      cout << "RiS";
      break;
    default:
      cout << "?";
      break;
  }
  cout << ", length " << length << ", configurations ";
  start.print(false);
  cout << " -> ";
  if (qi1)
  {
    qi1->print(false);
    cout << " -> ";
  }
  if (qi2)
  {
    qi2->print(false);
    cout << " -> ";
  }
  end.print(false);
  if (eol)
  {
    cout << endl;
  }
}

HC_CC_RS_Path::HC_CC_RS_Path(const Configuration &_start, const Configuration &_end, hc_cc_rs_path_type _type,
                             double _kappa, double _sigma, Configuration *_qi1, Configuration *_qi2,
                             Configuration *_qi3, Configuration *_qi4, HC_CC_Circle *_cstart, HC_CC_Circle *_cend,
                             HC_CC_Circle *_ci1, HC_CC_Circle *_ci2, double _length)
  : Path(_start, _end, _kappa, _sigma, _length)
{
  type = _type;
  qi1 = _qi1;
  qi2 = _qi2;
  qi3 = _qi3;
  qi4 = _qi4;
  cstart = _cstart;
  cend = _cend;
  ci1 = _ci1;
  ci2 = _ci2;
  ci3 = nullptr;
}

HC_CC_RS_Path::~HC_CC_RS_Path()
{
  delete qi1;
  delete qi2;
  delete qi3;
  delete qi4;
  delete cstart;
  delete ci1;
  delete ci2;
  delete ci3;
  delete cend;
}

void HC_CC_RS_Path::print(bool eol) const
{
  cout << "HC_CC_RS_Path: type ";
  switch (type)
  {
    case EMPTY:
      cout << "EMPTY";
      break;
    case STRAIGHT:
      cout << "STRAIGHT";
      break;
    case T:
      cout << "T";
      break;
    case TT:
      cout << "TT";
      break;
    case TcT:
      cout << "TcT";
      break;
    // Reeds-Shepp families:
    case TcTcT:
      cout << "TcTcT";
      break;
    case TcTT:
      cout << "TcTT";
      break;
    case TTcT:
      cout << "TTcT";
      break;
    case TST:
      cout << "TST";
      break;
    case TSTcT:
      cout << "TSTcT";
      break;
    case TcTST:
      cout << "TcTST";
      break;
    case TcTSTcT:
      cout << "TcTSTcT";
      break;
    case TTcTT:
      cout << "TTcTT";
      break;
    case TcTTcT:
      cout << "TcTTcT";
      break;
    // #####################
    case TTT:
      cout << "TTT";
      break;
    case TcST:
      cout << "TcST";
      break;
    case TScT:
      cout << "TScT";
      break;
    case TcScT:
      cout << "TcScT";
      break;
    default:
      cout << "?";
      break;
  }
  cout << ", length " << length << ", configurations ";
  start.print(false);
  cout << " -> ";
  if (qi1)
  {
    qi1->print(false);
    cout << " -> ";
  }
  if (qi2)
  {
    qi2->print(false);
    cout << " -> ";
  }
  if (qi3)
  {
    qi3->print(false);
    cout << " -> ";
  }
  if (qi4)
  {
    qi4->print(false);
    cout << " -> ";
  }
  end.print(false);
  if (eol)
  {
    cout << endl;
  }
}

void empty_controls(vector<Control> &controls)
{
  Control control;
  control.delta_s = 0.0;
  control.kappa = 0.0;
  control.sigma = 0.0;
  controls.push_back(control);
}

void straight_controls(const Configuration &q1, const Configuration &q2, vector<Control> &controls)
{
  double length = point_distance(q1.x, q1.y, q2.x, q2.y);
  double dot_product = cos(q1.theta) * (q2.x - q1.x) + sin(q1.theta) * (q2.y - q1.y);
  int d = sgn(dot_product);
  Control control;
  control.delta_s = d * length;
  control.kappa = 0.0;
  control.sigma = 0.0;
  controls.push_back(control);
}

int direction(bool forward, bool order)
{
  if (forward && order)
    return 1;
  else if (forward && !order)
    return -1;
  else if (!forward && order)
    return -1;
  else if (!forward && !order)
    return 1;
}

void rs_turn_controls(const HC_CC_Circle &c, const Configuration &q, bool order, vector<Control> &controls)
{
  assert(fabs(c.kappa) - fabs(q.kappa) < get_epsilon() &&
         fabs(c.sigma) - numeric_limits<double>::max() < get_epsilon());
  Control arc;
  double delta;
  c.deflection(q, &delta);
  int d = direction(c.forward, order);
  double length_arc;
  int shift;
  // irregular rs-turn
  if (!c.regular && (delta > PI))
  {
    shift = -d;
    length_arc = fabs((TWO_PI - delta) / c.kappa);
  }
  // regular rs-turn
  else
  {
    shift = d;
    length_arc = fabs(delta / c.kappa);
  }
  arc.delta_s = shift * length_arc;
  arc.kappa = c.kappa;
  arc.sigma = 0.0;
  controls.push_back(arc);
  return;
}

void hc_turn_controls(const HC_CC_Circle &c, const Configuration &q, bool order, vector<Control> &controls)
{
  assert(fabs(c.kappa) - fabs(q.kappa) < get_epsilon());
  Control clothoid, arc;
  double delta;
  c.deflection(q, &delta);
  int d = direction(c.forward, order);
  double length_min = fabs(c.kappa / c.sigma);
  double length_arc;
  int shift;
  // regular hc-turn
  if (c.regular && (delta < c.delta_min / 2.0))
  {
    shift = d;
    length_arc = (TWO_PI + delta - c.delta_min / 2.0) / fabs(c.kappa);
  }
  // irregular hc-turn
  else if (!c.regular && (delta > c.delta_min / 2.0 + PI))
  {
    shift = -d;
    length_arc = (TWO_PI - delta + c.delta_min / 2.0) / fabs(c.kappa);
  }
  // regular hc-turn
  else
  {
    shift = d;
    length_arc = (delta - c.delta_min / 2.0) / fabs(c.kappa);
  }
  if (order)
  {
    clothoid.delta_s = d * length_min;
    clothoid.kappa = 0.0;
    clothoid.sigma = c.sigma;
    controls.push_back(clothoid);
  }
  arc.delta_s = shift * length_arc;
  arc.kappa = c.kappa;
  arc.sigma = 0.0;
  controls.push_back(arc);
  if (!order)
  {
    clothoid.delta_s = d * length_min;
    clothoid.kappa = c.kappa;
    clothoid.sigma = -c.sigma;
    controls.push_back(clothoid);
  }
  return;
}

void cc_turn_controls(const HC_CC_Circle &c, const Configuration &q, bool order, vector<Control> &controls)
{
  assert(fabs(q.kappa) < get_epsilon());
  Control clothoid1, arc, clothoid2;
  double delta;
  c.deflection(q, &delta);
  int d = direction(c.forward, order);
  // straight line
  if (fabs(delta) < get_epsilon())
  {
    if (order)
      straight_controls(c.start, q, controls);
    else
      straight_controls(q, c.start, controls);
    return;
  }
  // elementary path
  if (delta < c.delta_min)
  {
    double d1 = D1(delta / 2);
    double d2 = point_distance(c.start.x, c.start.y, q.x, q.y);
    double sigma = 4 * PI * pow(d1, 2) / pow(d2, 2);
    double length = sqrt(delta / sigma);
    if (!c.left)
    {
      sigma = -sigma;
    }
    clothoid1.delta_s = d * length;
    clothoid1.kappa = 0.0;
    clothoid1.sigma = sigma;
    controls.push_back(clothoid1);

    clothoid2.delta_s = d * length;
    clothoid2.kappa = sigma*length;
    clothoid2.sigma = -sigma;
    controls.push_back(clothoid2);
    return;
  }
  // regular and irregular turn
  double length_min = fabs(c.kappa / c.sigma);
  double length_arc;
  int shift;
  // irregular
  if (!c.regular && (delta > c.delta_min + PI))
  {
    shift = -d;
    length_arc = (TWO_PI - delta + c.delta_min) / fabs(c.kappa);
  }
  // regular
  else
  {
    shift = d;
    length_arc = (delta - c.delta_min) / fabs(c.kappa);
  }
  clothoid1.delta_s = d * length_min;
  clothoid1.kappa = 0.0;
  clothoid1.sigma = c.sigma;
  controls.push_back(clothoid1);

  arc.delta_s = shift * length_arc;
  arc.kappa = c.kappa;
  arc.sigma = 0.0;
  controls.push_back(arc);

  clothoid2.delta_s = d * length_min;
  clothoid2.kappa = c.kappa;
  clothoid2.sigma = -c.sigma;
  controls.push_back(clothoid2);
  return;
}
