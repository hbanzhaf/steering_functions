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
from matplotlib.ticker import FuncFormatter

CSV_HEADER = "start,goal,computation_time,path_length".split(",")

def to_percent(x, pos):
    p = str(100 * x)
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
        self.path_length = np.array(self.path_length)
        self.comp_time = np.array(self.comp_time)
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

    print("\nComputation Times [µs]: mean ± std\n")
    n_samples = rs_outputs[0].n_samples
    for output in dubins_outputs + rs_outputs:
        assert output.n_samples == n_samples
        print(output.id + ":", np.average(output.comp_time) * 1e6, "±", np.std(output.comp_time) * 1e6)

    dubins_hist = []
    dubins_labels = []
    for output in dubins_outputs:
        if output.id != "Dubins":
            rel_path_length = (output.path_length - Dubins.path_length) / Dubins.path_length
            dubins_hist.append(rel_path_length)
            dubins_labels.append(output.id)

    rs_hist = []
    rs_labels = []
    for output in rs_outputs:
        if output.id != "RS":
            rel_path_length = (output.path_length - RS.path_length) / RS.path_length
            rs_hist.append(rel_path_length)
            rs_labels.append(output.id)


    # plot histograms
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4), sharey=True)
    f.subplots_adjust(wspace=.5)
    formatter = FuncFormatter(to_percent)
    weight = np.ones(n_samples, dtype='float')/n_samples

    ax1.hist(dubins_hist, bins=370, weights=[weight]*len(dubins_hist), label=dubins_labels, linewidth=.1)
    ax1.legend(loc='best')
    ax1.set_xlim([0, 0.4])
    ax1.grid('on')
    ax1.set_xlabel('Rel. Difference in Path Length to Dubins [%]')
    ax1.set_ylabel('Normalized Frequency [%]')
    ax1.xaxis.set_major_formatter(formatter)
    ax1.yaxis.set_major_formatter(formatter)

    ax2.hist(rs_hist, bins=300, weights=[weight]*len(rs_hist), label=rs_labels, linewidth=.1)
    ax2.legend(loc='best')
    ax2.set_xlim([0, 0.4])
    ax2.grid('on')
    ax2.set_xlabel('Rel. Difference in Path Length to Reeds-Shepp [%]')
    ax2.set_ylabel('Normalized Frequency [%]')
    plt.setp(ax2.get_yticklabels(), visible=True)
    ax2.xaxis.set_major_formatter(formatter)
    ax2.yaxis.set_major_formatter(formatter)

    # f.savefig('../doc/images/path_length_comparison.png', bbox_inches='tight', pad_inches=0)
    plt.show()
