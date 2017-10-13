#*********************************************************************
#  Copyright (c) 2017 Robert Bosch GmbH.
#  All rights reserved.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# *********************************************************************/

import os.path
import numpy as np
import matplotlib.pyplot as plt

CSV_HEADER = "x,y,theta,kappa,d".split(",")

def read_csv(filename):
    with open(filename, 'r') as f:
        next(f)
        lines = [line.strip().split(",") for line in f.readlines()]
        f.close
        return lines


def distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)


class State:
    def __init__(self, x, y, theta, kappa, d):
        self.x = x
        self.y = y
        self.theta = theta
        self.kappa = kappa
        self.d = d


class Path:
    def __init__(self):
        self.states = []
        self.x = None
        self.y = None
        self.theta = None
        self.kappa = []
        self.d = []
        self.s = []

    def load(self, fpath, fname):
        if not os.path.exists(fpath + fname):
            return 0
        else:
            output = read_csv(fpath + fname)
            for line in output:
                x = np.float(line[CSV_HEADER.index("x")])
                y = np.float(line[CSV_HEADER.index("y")])
                theta = np.float(line[CSV_HEADER.index("theta")])
                kappa = np.float(line[CSV_HEADER.index("kappa")])
                d = np.float(line[CSV_HEADER.index("d")])
                self.states.append(State(x, y, theta, kappa, d))
        return 1

    def postprocess(self):
        self.x = [state.x for state in self.states]
        self.y = [state.y for state in self.states]
        self.theta = [state.theta for state in self.states]
        self.kappa = [state.kappa for state in self.states]
        self.d = [state.d for state in self.states]

        s = 0.0
        state1 = self.states[0]
        for state2 in self.states:
            s += distance(state1.x, state1.y, state2.x, state2.y)
            self.s.append(s)
            state1 = state2


if __name__ == "__main__":
    # load data
    filepath = "./"

    path = Path()
    path.load(filepath, "path.csv")
    path.postprocess()

    # plot
    f, axarr = plt.subplots(2, 2, figsize=(16, 8))

    axarr[0, 0].plot(path.x, path.y)
    axarr[0, 0].set_xlabel('$x$ [m]')
    axarr[0, 0].set_ylabel('$y$ [m]')
    axarr[0, 0].set_aspect('equal')

    axarr[0, 1].plot(path.s, path.theta)
    axarr[0, 1].set_xlabel('$s$ [m]')
    axarr[0, 1].set_ylabel('$theta$ [rad]')

    ax1 = axarr[1, 0]
    ax2 = axarr[1, 0].twinx()
    handle1, = ax1.plot(path.s, path.kappa, label='$kappa$')
    handle2, = ax2.plot(path.s, path.d, label='$d$', color='red', alpha=0.4)
    axarr[1, 0].set_xlabel('$s$ [m]')
    ax1.set_ylabel('$kappa$ [1/m]')
    ax2.set_ylabel('$d$ [-]')
    axarr[1, 0].legend(handles=[handle1, handle2])

    plt.show()
