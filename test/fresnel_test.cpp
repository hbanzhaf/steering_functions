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
#include <gtest/gtest.h>

#include "steering_functions/utilities/utilities.hpp"

using namespace std;
using namespace steering;

#define EPS_FRESNEL 1e-14

// Comparison with values from scipy.special.fresnel
TEST(Fresnel, fresnel)
{
  double fresnel_s, fresnel_c;

  fresnel(-40.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4920422537902731), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4999984168574456), EPS_FRESNEL);

  fresnel(-39.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49255430877943795), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5030823268038184), EPS_FRESNEL);

  fresnel(-39.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49999829192809686), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5081617908810524), EPS_FRESNEL);

  fresnel(-38.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5076391196617077), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.496837694806001), EPS_FRESNEL);

  fresnel(-38.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49162342526889974), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.499998153500694), EPS_FRESNEL);

  fresnel(-37.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5078428670985601), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49675345773238233), EPS_FRESNEL);

  fresnel(-37.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49999799970186465), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5086029685015766), EPS_FRESNEL);

  fresnel(-36.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49194221801222643), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5033353872735516), EPS_FRESNEL);

  fresnel(-36.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4911580603172578), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4999978283373622), EPS_FRESNEL);

  fresnel(-35.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.491715191566429), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5034292286980901), EPS_FRESNEL);

  fresnel(-35.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49999763682609794), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5090945663345068), EPS_FRESNEL);

  fresnel(-34.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5085250000581694), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4964715000692896), EPS_FRESNEL);

  fresnel(-34.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4906379466534976), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49999742211814463), EPS_FRESNEL);

  fresnel(-33.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5087795363583066), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4963663140599485), EPS_FRESNEL);

  fresnel(-33.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4999971805923199), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5096457516544846), EPS_FRESNEL);

  fresnel(-32.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4909502579162105), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5037453311765392), EPS_FRESNEL);

  fresnel(-32.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49005281894025865), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4999969079273435), EPS_FRESNEL);

  fresnel(-31.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49066288968973953), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5038640489726605), EPS_FRESNEL);

  fresnel(-31.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49999659893870957), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.510268057465068), EPS_FRESNEL);

  fresnel(-30.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5096433300548732), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4960094672174812), EPS_FRESNEL);

  fresnel(-30.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4893896744421938), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4999962473706091), EPS_FRESNEL);

  fresnel(-29.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5099703195139875), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49587443030573236), EPS_FRESNEL);

  fresnel(-29.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4999958456285237), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5109761982547072), EPS_FRESNEL);

  fresnel(-28.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4896797336947276), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5042700567697976), EPS_FRESNEL);

  fresnel(-28.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4886317954009967), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4999953844327101), EPS_FRESNEL);

  fresnel(-27.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4893043235007156), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5044250209259784), EPS_FRESNEL);

  fresnel(-27.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4999948523652942), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5117892483008945), EPS_FRESNEL);

  fresnel(-26.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5110994347652468), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49540835671456757), EPS_FRESNEL);

  fresnel(-26.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4877573202131747), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49999423527272), EPS_FRESNEL);

  fresnel(-25.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5115348786861477), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49522871078678943), EPS_FRESNEL);

  fresnel(-25.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49999351546947596), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5127323855397702), EPS_FRESNEL);

  fresnel(-24.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4879941087024924), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5049655456487026), EPS_FRESNEL);

  fresnel(-24.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.48673710022664085), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.499992670665544), EPS_FRESNEL);

  fresnel(-23.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4874829827132583), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5051762678838345), EPS_FRESNEL);

  fresnel(-23.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4999916725048582), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5138395452365079), EPS_FRESNEL);

  fresnel(-22.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5130736102268015), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4945943613060027), EPS_FRESNEL);

  fresnel(-22.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4855313875836048), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4999904845486192), EPS_FRESNEL);

  fresnel(-21.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5136820209670363), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49434375620493276), EPS_FRESNEL);

  fresnel(-21.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4999890594544999), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5151575899376805), EPS_FRESNEL);

  fresnel(-20.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.48565015873667594), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5059311691275604), EPS_FRESNEL);

  fresnel(-20.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4840845359259539), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49998733497234443), EPS_FRESNEL);

  fresnel(-19.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4849137774794187), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5062341276891209), EPS_FRESNEL);

  fresnel(-19.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4999852281670667), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5167531128300553), EPS_FRESNEL);

  fresnel(-18.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5159022981370349), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4934303730848653), EPS_FRESNEL);

  fresnel(-18.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.48231616863711996), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49998262693469964), EPS_FRESNEL);

  fresnel(-17.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5168117510006868), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4930568074654182), EPS_FRESNEL);

  fresnel(-17.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4999793772969549), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5187240428109737), EPS_FRESNEL);

  fresnel(-16.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4821684121023745), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5073616723760539), EPS_FRESNEL);

  fresnel(-16.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.48010572438090104), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4999752639565001), EPS_FRESNEL);

  fresnel(-15.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4810167854210039), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5078336554052517), EPS_FRESNEL);

  fresnel(-15.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4999699798097023), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5212205316743734), EPS_FRESNEL);

  fresnel(-14.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5202939571560238), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49162993982669234), EPS_FRESNEL);

  fresnel(-14.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4772637594418203), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49996307683096536), EPS_FRESNEL);

  fresnel(-13.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.521799262187777), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4910150221176157), EPS_FRESNEL);

  fresnel(-13.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4999538844819125), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5244851153043584), EPS_FRESNEL);

  fresnel(-12.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4764540427410991), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5096969076697003), EPS_FRESNEL);

  fresnel(-12.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.47347456491993545), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49994136935201117), EPS_FRESNEL);

  fresnel(-11.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.47440277911222917), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5105306129965809), EPS_FRESNEL);

  fresnel(-11.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4999238837937511), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5289366617557965), EPS_FRESNEL);

  fresnel(-10.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5280404079981299), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4884800073027095), EPS_FRESNEL);

  fresnel(-10.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.46816997858488224), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49989869420551575), EPS_FRESNEL);

  fresnel(-9.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5309998491513985), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4872873310265677), EPS_FRESNEL);

  fresnel(-9.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4998610456296846), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5353661274681198), EPS_FRESNEL);

  fresnel(-8.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.46534124898107443), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5141775985837332), EPS_FRESNEL);

  fresnel(-8.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.46021421439301446), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.499802180377197), EPS_FRESNEL);

  fresnel(-7.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4607012329468305), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5160182501523362), EPS_FRESNEL);

  fresnel(-7.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4997047894534464), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5454670925469698), EPS_FRESNEL);

  fresnel(-6.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5453764552432336), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.48160345989096404), EPS_FRESNEL);

  fresnel(-6.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4469607612369303), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.499531467855501), EPS_FRESNEL);

  fresnel(-5.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.5536840627790217), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4784214149253144), EPS_FRESNEL);

  fresnel(-5.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4991913819171169), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5636311887040122), EPS_FRESNEL);

  fresnel(-4.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.43427297504870355), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5260259150535386), EPS_FRESNEL);

  fresnel(-4.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.42051575424692844), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.49842603303817756), EPS_FRESNEL);

  fresnel(-3.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.41524801197243744), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.5325724350280007), EPS_FRESNEL);

  fresnel(-3.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.49631299896737496), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.6057207892976857), EPS_FRESNEL);

  fresnel(-2.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.6191817558195929), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.45741300964177706), EPS_FRESNEL);

  fresnel(-2.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.34341567836369824), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.48825340607534073), EPS_FRESNEL);

  fresnel(-1.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.697504960082093), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.44526117603982157), EPS_FRESNEL);

  fresnel(-1.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.4382591473903547), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.779893400376823), EPS_FRESNEL);

  fresnel(-0.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - -0.06473243285999929), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - -0.4923442258714464), EPS_FRESNEL);

  fresnel(0.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.0), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.0), EPS_FRESNEL);

  fresnel(0.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.06473243285999929), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4923442258714464), EPS_FRESNEL);

  fresnel(1.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4382591473903547), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.779893400376823), EPS_FRESNEL);

  fresnel(1.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.697504960082093), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.44526117603982157), EPS_FRESNEL);

  fresnel(2.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.34341567836369824), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.48825340607534073), EPS_FRESNEL);

  fresnel(2.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.6191817558195929), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.45741300964177706), EPS_FRESNEL);

  fresnel(3.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49631299896737496), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.6057207892976857), EPS_FRESNEL);

  fresnel(3.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.41524801197243744), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5325724350280007), EPS_FRESNEL);

  fresnel(4.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.42051575424692844), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49842603303817756), EPS_FRESNEL);

  fresnel(4.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.43427297504870355), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5260259150535386), EPS_FRESNEL);

  fresnel(5.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4991913819171169), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5636311887040122), EPS_FRESNEL);

  fresnel(5.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5536840627790217), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4784214149253144), EPS_FRESNEL);

  fresnel(6.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4469607612369303), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.499531467855501), EPS_FRESNEL);

  fresnel(6.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5453764552432336), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.48160345989096404), EPS_FRESNEL);

  fresnel(7.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4997047894534464), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5454670925469698), EPS_FRESNEL);

  fresnel(7.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4607012329468305), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5160182501523362), EPS_FRESNEL);

  fresnel(8.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.46021421439301446), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.499802180377197), EPS_FRESNEL);

  fresnel(8.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.46534124898107443), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5141775985837332), EPS_FRESNEL);

  fresnel(9.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4998610456296846), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5353661274681198), EPS_FRESNEL);

  fresnel(9.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5309998491513985), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4872873310265677), EPS_FRESNEL);

  fresnel(10.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.46816997858488224), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49989869420551575), EPS_FRESNEL);

  fresnel(10.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5280404079981299), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4884800073027095), EPS_FRESNEL);

  fresnel(11.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4999238837937511), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5289366617557965), EPS_FRESNEL);

  fresnel(11.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.47440277911222917), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5105306129965809), EPS_FRESNEL);

  fresnel(12.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.47347456491993545), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49994136935201117), EPS_FRESNEL);

  fresnel(12.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4764540427410991), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5096969076697003), EPS_FRESNEL);

  fresnel(13.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4999538844819125), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5244851153043584), EPS_FRESNEL);

  fresnel(13.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.521799262187777), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4910150221176157), EPS_FRESNEL);

  fresnel(14.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4772637594418203), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49996307683096536), EPS_FRESNEL);

  fresnel(14.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5202939571560238), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49162993982669234), EPS_FRESNEL);

  fresnel(15.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4999699798097023), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5212205316743734), EPS_FRESNEL);

  fresnel(15.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4810167854210039), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5078336554052517), EPS_FRESNEL);

  fresnel(16.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.48010572438090104), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4999752639565001), EPS_FRESNEL);

  fresnel(16.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4821684121023745), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5073616723760539), EPS_FRESNEL);

  fresnel(17.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4999793772969549), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5187240428109737), EPS_FRESNEL);

  fresnel(17.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5168117510006868), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4930568074654182), EPS_FRESNEL);

  fresnel(18.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.48231616863711996), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49998262693469964), EPS_FRESNEL);

  fresnel(18.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5159022981370349), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4934303730848653), EPS_FRESNEL);

  fresnel(19.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4999852281670667), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5167531128300553), EPS_FRESNEL);

  fresnel(19.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4849137774794187), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5062341276891209), EPS_FRESNEL);

  fresnel(20.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4840845359259539), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49998733497234443), EPS_FRESNEL);

  fresnel(20.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.48565015873667594), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5059311691275604), EPS_FRESNEL);

  fresnel(21.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4999890594544999), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5151575899376805), EPS_FRESNEL);

  fresnel(21.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5136820209670363), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49434375620493276), EPS_FRESNEL);

  fresnel(22.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4855313875836048), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4999904845486192), EPS_FRESNEL);

  fresnel(22.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5130736102268015), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4945943613060027), EPS_FRESNEL);

  fresnel(23.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4999916725048582), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5138395452365079), EPS_FRESNEL);

  fresnel(23.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4874829827132583), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5051762678838345), EPS_FRESNEL);

  fresnel(24.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.48673710022664085), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.499992670665544), EPS_FRESNEL);

  fresnel(24.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4879941087024924), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5049655456487026), EPS_FRESNEL);

  fresnel(25.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49999351546947596), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5127323855397702), EPS_FRESNEL);

  fresnel(25.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5115348786861477), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49522871078678943), EPS_FRESNEL);

  fresnel(26.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4877573202131747), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49999423527272), EPS_FRESNEL);

  fresnel(26.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5110994347652468), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49540835671456757), EPS_FRESNEL);

  fresnel(27.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4999948523652942), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5117892483008945), EPS_FRESNEL);

  fresnel(27.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4893043235007156), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5044250209259784), EPS_FRESNEL);

  fresnel(28.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4886317954009967), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4999953844327101), EPS_FRESNEL);

  fresnel(28.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4896797336947276), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5042700567697976), EPS_FRESNEL);

  fresnel(29.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4999958456285237), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5109761982547072), EPS_FRESNEL);

  fresnel(29.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5099703195139875), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49587443030573236), EPS_FRESNEL);

  fresnel(30.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4893896744421938), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4999962473706091), EPS_FRESNEL);

  fresnel(30.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5096433300548732), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4960094672174812), EPS_FRESNEL);

  fresnel(31.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49999659893870957), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.510268057465068), EPS_FRESNEL);

  fresnel(31.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49066288968973953), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5038640489726605), EPS_FRESNEL);

  fresnel(32.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49005281894025865), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4999969079273435), EPS_FRESNEL);

  fresnel(32.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4909502579162105), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5037453311765392), EPS_FRESNEL);

  fresnel(33.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4999971805923199), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5096457516544846), EPS_FRESNEL);

  fresnel(33.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5087795363583066), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4963663140599485), EPS_FRESNEL);

  fresnel(34.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4906379466534976), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49999742211814463), EPS_FRESNEL);

  fresnel(34.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5085250000581694), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4964715000692896), EPS_FRESNEL);

  fresnel(35.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49999763682609794), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5090945663345068), EPS_FRESNEL);

  fresnel(35.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.491715191566429), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5034292286980901), EPS_FRESNEL);

  fresnel(36.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4911580603172578), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4999978283373622), EPS_FRESNEL);

  fresnel(36.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49194221801222643), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5033353872735516), EPS_FRESNEL);

  fresnel(37.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49999799970186465), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5086029685015766), EPS_FRESNEL);

  fresnel(37.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5078428670985601), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.49675345773238233), EPS_FRESNEL);

  fresnel(38.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49162342526889974), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.499998153500694), EPS_FRESNEL);

  fresnel(38.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.5076391196617077), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.496837694806001), EPS_FRESNEL);

  fresnel(39.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49999829192809686), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5081617908810524), EPS_FRESNEL);

  fresnel(39.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.49255430877943795), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5030823268038184), EPS_FRESNEL);

  fresnel(40.0, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4920422537902731), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.4999984168574456), EPS_FRESNEL);

  fresnel(40.5, fresnel_s, fresnel_c);
  EXPECT_LE(fabs(fresnel_s - 0.4927381828482828), EPS_FRESNEL);
  EXPECT_LE(fabs(fresnel_c - 0.5030062922561774), EPS_FRESNEL);
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
