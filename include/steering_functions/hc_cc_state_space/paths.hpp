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

#ifndef PATHS_HPP
#define PATHS_HPP

#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

#include "configuration.hpp"
#include "hc_cc_circle.hpp"
#include "steering_functions/steering_functions.hpp"
#include "steering_functions/utilities/utilities.hpp"

using namespace std;
using namespace steer;

class Path
{
public:
  /** \brief Constructor */
  Path(const Configuration &_start, const Configuration &_end, double _kappa, double _sigma, double _length);

  /** \brief Start and end configuration */
  Configuration start, end;

  /** \brief Max. curvature (unsigned), max. sharpness (unsigned) */
  double kappa, sigma;

  /** \brief Path length */
  double length;
};

/** \brief cc-dubins path types: E (Empty), L (Left), S (Straight), R (Right) */
enum cc_dubins_path_type
{
  E,
  LSL,
  LSR,
  RSL,
  RSR,
  LR1L,
  LR2L,
  RL1R,
  RL2R,
  S,
  L,
  R,
  LeS,
  LiS,
  eSL,
  iSL,
  eSR,
  iSR,
  ReS,
  RiS
};
const int nb_cc_dubins_paths = 20;

class CC_Dubins_Path : public Path
{
public:
  /** \brief Constructor */
  CC_Dubins_Path(const Configuration &_start, const Configuration &_end, cc_dubins_path_type _type, double _kappa,
                 double _sigma, Configuration *_qi1, Configuration *_qi2, HC_CC_Circle *_cstart, HC_CC_Circle *_cend,
                 HC_CC_Circle *_ci1, double _length);

  /** \brief Destructor */
  ~CC_Dubins_Path();

  /** \brief Alphanumeric display */
  void print(bool eol) const;

  /** \brief Path type */
  cc_dubins_path_type type;

  /** \brief Intermediate configurations */
  Configuration *qi1, *qi2;

  /** \brief Start, end and intermediate circles */
  HC_CC_Circle *cstart, *cend, *ci1, *ci2;
};

/** \brief hc-/cc-reeds-shepp path types: T (Turn), S (Straight), c (Cusp) */
enum hc_cc_rs_path_type
{
  EMPTY,
  STRAIGHT,
  T,
  TT,
  TcT,
  // Reeds-Shepp families:
  TcTcT,
  TcTT,
  TTcT,
  TST,
  TSTcT,
  TcTST,
  TcTSTcT,
  TTcTT,
  TcTTcT,
  // #####################
  TTT,
  TcST,
  TScT,
  TcScT
};
const int nb_hc_cc_rs_paths = 18;

class HC_CC_RS_Path : public Path
{
public:
  /** \brief Constructor */
  HC_CC_RS_Path(const Configuration &_start, const Configuration &_end, hc_cc_rs_path_type _type, double _kappa,
                double _sigma, Configuration *_qi1, Configuration *_qi2, Configuration *_qi3, Configuration *_qi4,
                HC_CC_Circle *_cstart, HC_CC_Circle *_cend, HC_CC_Circle *_ci1, HC_CC_Circle *_ci2, double _length);

  /** \brief Destructor */
  ~HC_CC_RS_Path();

  /** \brief Alphanumeric display */
  void print(bool eol) const;

  /** \brief Path type */
  hc_cc_rs_path_type type;

  /** \brief Intermediate configurations */
  Configuration *qi1, *qi2, *qi3, *qi4;

  /** \brief Start, end and intermediate circles */
  HC_CC_Circle *cstart, *cend, *ci1, *ci2, *ci3;
};

/** \brief Reverses a control */
void reverse_control(Control &control);

/** \brief Appends controls with 0 input */
void empty_controls(vector<Control> &controls);

/** \brief Appends controls with a straight line */
void straight_controls(const Configuration &q1, const Configuration &q2, vector<Control> &controls);

/** \brief Appends controls with a rs-turn */
void rs_turn_controls(const HC_CC_Circle &c, const Configuration &q, bool order, vector<Control> &controls);

/** \brief Appends controls with a hc-turn */
void hc_turn_controls(const HC_CC_Circle &c, const Configuration &q, bool order, vector<Control> &controls);

/** \brief Appends controls with a cc-turn */
void cc_turn_controls(const HC_CC_Circle &c, const Configuration &q, bool order, vector<Control> &controls);

#endif
