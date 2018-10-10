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
from matplotlib.ticker import FuncFormatter, MaxNLocator

CSV_HEADER = "start,goal,computation_time,path_length,curv_discont".split(",")

def to_percent(x, pos):
    p = "%.1f" % (100 * x)
    return p

def read_csv(filename):
    with open(filename, 'r') as f:
        next(f)
        lines = [line.strip().split(",") for line in f.readlines()]
        f.close
        return lines

class Output:
    def __init__(self):
        self.id = None
        self.n_samples = None
        self.path_length = []
        self.comp_time = []
        self.curv_discont = []

    def load(self, fpath, fname):
        if not os.path.exists(fpath + fname):
            return 0
        if "CC00_Dubins" in fname:
            self.id = "CC$^{00}$-Dubins"
        elif "CC0pm_Dubins" in fname:
            self.id = "CC$^{0\pm}$-Dubins"
        elif "CCpm0_Dubins" in fname:
            self.id = "CC$^{\pm0}$-Dubins"
        elif "CCpmpm_Dubins" in fname:
            self.id = "CC$^{\pm\pm}$-Dubins"
        elif "CC_Dubins" in fname:
            self.id = "CC-Dubins"
        elif "Dubins" in fname:
            self.id = "Dubins"
        elif "CC00_RS" in fname:
            self.id = "CC$^{00}$-RS"
        elif "HC00_RS" in fname:
            self.id = "HC$^{00}$-RS"
        elif "HC0pm_RS" in fname:
            self.id = "HC$^{0\pm}$-RS"
        elif "HCpm0_RS" in fname:
            self.id = "HC$^{\pm0}$-RS"
        elif "HCpmpm_RS" in fname:
            self.id = "HC$^{\pm\pm}$-RS"
        elif "HC_RS" in fname:
            self.id = "HC-RS"
        elif "RS" in fname:
            self.id = "RS"

        stats = read_csv(fpath + fname)
        self.n_samples = len(stats)
        for stat in stats:
            self.path_length.append(np.float(stat[CSV_HEADER.index("path_length")]))
            self.comp_time.append(np.float(stat[CSV_HEADER.index("computation_time")]))
            self.curv_discont.append(np.int(stat[CSV_HEADER.index("curv_discont")]))
        self.path_length = np.array(self.path_length)
        self.comp_time = np.array(self.comp_time)
        self.curv_discont = np.array(self.curv_discont)
        return 1

if __name__ == "__main__":
    # load data
    filepath = "../test/"

    dubins_outputs = []
    Dubins = Output()
    CC_Dubins = Output()
    CC00_Dubins = Output()
    CC0pm_Dubins = Output()
    CCpm0_Dubins = Output()
    CCpmpm_Dubins = Output()
    if Dubins.load(filepath, "Dubins_stats.csv"):
        dubins_outputs.append(Dubins)
    if CCpmpm_Dubins.load(filepath, "CCpmpm_Dubins_stats.csv"):
        dubins_outputs.append(CCpmpm_Dubins)
    if CCpm0_Dubins.load(filepath, "CCpm0_Dubins_stats.csv"):
        dubins_outputs.append(CCpm0_Dubins)
    if CC0pm_Dubins.load(filepath, "CC0pm_Dubins_stats.csv"):
        dubins_outputs.append(CC0pm_Dubins)
    if CC00_Dubins.load(filepath, "CC00_Dubins_stats.csv"):
        dubins_outputs.append(CC00_Dubins)
    if CC_Dubins.load(filepath, "CC_Dubins_stats.csv"):
        dubins_outputs.append(CC_Dubins)

    rs_outputs = []
    CC00_RS = Output()
    HC_RS = Output()
    HC00_RS = Output()
    HC0pm_RS = Output()
    HCpm0_RS = Output()
    HCpmpm_RS = Output()
    RS = Output()
    if RS.load(filepath, "RS_stats.csv"):
        rs_outputs.append(RS)
    if HCpmpm_RS.load(filepath, "HCpmpm_RS_stats.csv"):
        rs_outputs.append(HCpmpm_RS)
    if HCpm0_RS.load(filepath, "HCpm0_RS_stats.csv"):
        rs_outputs.append(HCpm0_RS)
    if HC0pm_RS.load(filepath, "HC0pm_RS_stats.csv"):
        rs_outputs.append(HC0pm_RS)
    if HC00_RS.load(filepath, "HC00_RS_stats.csv"):
        rs_outputs.append(HC00_RS)
    if HC_RS.load(filepath, "HC_RS_stats.csv"):
        rs_outputs.append(HC_RS)
    if CC00_RS.load(filepath, "CC00_RS_stats.csv"):
        rs_outputs.append(CC00_RS)

    # computation time
    print("\nComputation Times [µs]: mean ± std\n")
    n_samples = rs_outputs[0].n_samples
    for output in dubins_outputs + rs_outputs:
        assert output.n_samples == n_samples
        print(output.id + ": %.2f ± %.2f" % (np.average(output.comp_time) * 1e6, np.std(output.comp_time) * 1e6))

    # path length
    dubins_path_length_hist = []
    dubins_path_length_labels = []
    for output in dubins_outputs:
        if output.id != "Dubins":
            rel_path_length = (output.path_length - Dubins.path_length) / Dubins.path_length
            dubins_path_length_hist.append(rel_path_length)
            dubins_path_length_labels.append(output.id)

    rs_path_length_hist = []
    rs_path_length_labels = []
    for output in rs_outputs:
        if output.id != "RS":
            rel_path_length = (output.path_length - RS.path_length) / RS.path_length
            rs_path_length_hist.append(rel_path_length)
            rs_path_length_labels.append(output.id)

    # plot
    f, axarr = plt.subplots(2, 2, figsize=(10, 8), sharey=True)
    f.subplots_adjust(wspace=.4, hspace=.3)
    weights = np.ones(n_samples, dtype='float') / n_samples

    # path length histograms
    xlim = (0.0, 0.4)
    axarr[0, 0].plot([], [])  # advance color cycle once
    axarr[0, 0].hist(dubins_path_length_hist, bins=9, range=xlim, weights=[weights]*len(dubins_path_length_hist),
                     label=dubins_path_length_labels, linewidth=.1)
    axarr[0, 0].legend(loc='best')
    axarr[0, 0].set_xlim(xlim)
    axarr[0, 0].grid('on')
    axarr[0, 0].set_xlabel('Rel. Difference in Path Length to Dubins [%]')
    axarr[0, 0].set_ylabel('Normalized Frequency [%]')
    axarr[0, 0].xaxis.set_major_formatter(FuncFormatter(to_percent))
    axarr[0, 0].yaxis.set_major_formatter(FuncFormatter(to_percent))

    axarr[0, 1].plot([], [])  # advance color cycle once
    axarr[0, 1].hist(rs_path_length_hist, bins=9, range=xlim, weights=[weights]*len(rs_path_length_hist),
                     label=rs_path_length_labels, linewidth=.1)
    axarr[0, 1].legend(loc='best')
    axarr[0, 1].set_xlim(xlim)
    axarr[0, 1].grid('on')
    axarr[0, 1].set_xlabel('Rel. Difference in Path Length to Reeds-Shepp [%]')
    axarr[0, 1].set_ylabel('Normalized Frequency [%]')
    plt.setp(axarr[0, 1].get_yticklabels(), visible=True)
    axarr[0, 1].xaxis.set_major_formatter(FuncFormatter(to_percent))
    axarr[0, 1].yaxis.set_major_formatter(FuncFormatter(to_percent))

    # curvature discont. histrograms
    width = 0.06
    for i, output in enumerate(dubins_outputs):
        hist, bin_edges = np.histogram(output.curv_discont, weights=weights, bins=range(4))
        axarr[1, 0].bar(bin_edges[:-1] + width*(i - len(dubins_outputs)/2.0), hist, width=width, label=output.id)
    axarr[1, 0].legend(loc='best')
    axarr[1, 0].set_xlim(right=4)
    axarr[1, 0].grid('on')
    axarr[1, 0].set_xlabel('Number of Curvature Discontinuities [-]')
    axarr[1, 0].set_ylabel('Normalized Frequency [%]')
    axarr[1, 0].xaxis.set_major_locator(MaxNLocator(integer=True))
    axarr[1, 0].yaxis.set_major_formatter(FuncFormatter(to_percent))

    width = 0.08
    for i, output in enumerate(rs_outputs):
        hist, bin_edges = np.histogram(output.curv_discont, weights=weights, bins=range(6))
        axarr[1, 1].bar(bin_edges[:-1] + width*(i - len(rs_outputs)/2.0), hist, width=width, label=output.id)
    axarr[1, 1].legend(loc='best')
    axarr[1, 1].set_xlim(right=5)
    axarr[1, 1].grid('on')
    axarr[1, 1].set_xlabel('Number of Curvature Discontinuities [-]')
    axarr[1, 1].set_ylabel('Normalized Frequency [%]')
    plt.setp(axarr[1, 1].get_yticklabels(), visible=True)
    axarr[1, 1].xaxis.set_major_locator(MaxNLocator(integer=True))
    axarr[1, 1].yaxis.set_major_formatter(FuncFormatter(to_percent))

    # f.savefig('../doc/images/statistics.png', bbox_inches='tight', pad_inches=0)
    plt.show()
