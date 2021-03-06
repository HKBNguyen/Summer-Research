import csv
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from pandas import read_csv
import pandas as pd
import re
import glob



s_List = []
l_List = []
OR_proper = ['nC-0_5', 'nC-4_5', 'nC-8_5']
AO_proper = ['nC-0_5','nC-1_5','nC-2_5']
AO_nonproper = ['nC-0','nC-4','nC-8']
OR_proper_res = []
AO_proper_res = []
AO_nonproper_res = []



#no_ext is a list of filenames without their extensions. spiltext is able to separate them for every file in the directory iff that certain file ends with the .out extension
no_ext = [os.path.splitext(x)[0] for x in os.listdir('.') if x.endswith(".out")]
#list of problems with their respective tree types
pTree = []
pNames = []
for filename in no_ext:
    name_trunc = filename.split('-i-')[0]
    if name_trunc not in pTree:
        pTree.append(name_trunc)

for file in os.listdir('.'):
        if file.endswith('.out'):
            pn = file.split('.uai')[0]
            pNames.append(pn)

def findZ (problem_name):
    Z_Dict ={
        #pedigrees
        "pedigree1": -14.11,
        "pedigree13": -31.18,
        "pedigree18": -78.14,
        "pedigree19": -59.02,
        "pedigree20": -29.63,
        "pedigree23": -38.56,
        "pedigree25": -115.78,
        "pedigree30": -83.73,
        "pedigree31": -69.70,
        "pedigree33": -54.28,
        "pedigree34": -64.23,
        "pedigree37": -116.58,
        "pedigree38": -54.25,
        "pedigree39": -102.20,
        "pedigree40": -87.88,
        "pedigree41": -76.04,
        "pedigree42": -30.76,
        "pedigree44": -63.47,
        "pedigree50": -22.88,
        "pedigree51": -77.27,
        "pedigree7": -64.82,
        "pedigree9": -78.52,
        "Family2Dominant.1.5loci": -17.22,
        "Family2Dominant.20.5loci": -11.46,
        "Family2Recessive.15.5loci": -15.21,
        #grids
        "grid10x10.f10": 303.09,
        "grid10x10.f10.wrap": 333.32,
        "grid10x10.f15.wrap": 497.76,
        "grid10x10.f5.wrap": 169.41,
        "grid20x20.f10": 1311.98,
        "grid20x20.f15": 1962.98,
        "grid20x20.f2": 291.73,
        "grid20x20.f5": 665.12,
        #promedas
        "or_chain_1.fg": -10.76,
        "or_chain_10.fg": -8.39,
        "or_chain_102.fg": -10.04,
        "or_chain_106.fg": -10.06,
        "or_chain_107.fg": -11.57,
        "or_chain_111.fg": -5.86,
        "or_chain_12.fg": -12.22,
        "or_chain_128.fg": -10.91,
        "or_chain_129.fg": -9.42,
        "or_chain_132.fg": -8.43,
        "or_chain_133.fg": -6.71,
        "or_chain_138.fg": -10.05,
        "or_chain_139.fg": -7.00,
        "or_chain_140.fg": -5.54,
        "or_chain_147.fg": -8.14,
        "or_chain_148.fg": -8.13,
        "or_chain_149.fg": -8.60,
        "or_chain_15.fg": -11.85,
        "or_chain_150.fg": -7.25,
        "or_chain_152.fg": -10.15,
        "or_chain_153.fg": -8.69,
        "or_chain_154.fg": -1.19,
        "or_chain_155.fg": -9.95,
        "or_chain_161.fg": -9.24,
        "or_chain_165.fg": -21.80,
        "or_chain_168.fg": -1.80,
        "or_chain_17.fg": -3.16,
        "or_chain_175.fg": -1.27,
        "or_chain_176.fg": -2.20,
        "or_chain_18.fg": -4.05,
        "or_chain_180.fg": -2.81,
        "or_chain_182.fg": -3.08,
        "or_chain_186.fg": -9.47,
        "or_chain_188.fg": -12.31,
        "or_chain_197.fg": -1.19,
        "or_chain_198.fg": -4.30,
        "or_chain_201.fg": -1.74,
        "or_chain_203.fg": -1.25,
        "or_chain_209.fg": -11.46,
        "or_chain_214.fg": -14.97,
        "or_chain_218.fg": -1.49,
        "or_chain_220.fg": -10.33,
        "or_chain_225.fg": -1.74,
        "or_chain_236.fg": -4.17,
        "or_chain_242.fg": -5.35,
        "or_chain_247.fg": -4.98,
        "or_chain_248.fg": -11.55,
        "or_chain_29.fg": -11.93,
        "or_chain_38.fg": -7.51,
        "or_chain_4.fg": -11.72,
        "or_chain_42.fg": -3.34,
        "or_chain_45.fg": -6.98,
        "or_chain_48.fg": -4.79,
        "or_chain_53.fg": -9.18,
        "or_chain_61.fg": -8.00,
        "or_chain_62.fg": -9.46,
        "or_chain_64.fg": -8.84,
        "or_chain_65.fg": -4.67,
        "or_chain_68.fg": -4.34,
        "or_chain_72.fg": -7.06,
        "or_chain_74.fg": -5.58,
        "or_chain_85.fg": -2.49,
        "or_chain_90.fg": -8.55,
        "or_chain_91.fg": -11.32,
        "or_chain_93.fg": -8.72,
        #DBN
        "rbm_20": 58.53,
        "rbm_21": 63.15,
        "rbm_22": 66.55,
        "rbm_ferro_20": 151.16,
        "rbm_ferro_21": 152.62,
        "rbm_ferro_22": 166.11,
        "rus_20_40_0_1": 617.31,
        "rus_20_40_0_2": 791.57,
        "rus_20_40_0_3": 903.04,
        "rus_20_40_1_1": 1004.68,
        "rus_20_40_1_2": 910.22,
        "rus_20_40_1_3": 900.64,
        "rus_20_40_2_1": 797.82,
        "rus_20_40_2_2": 689.28,
        "rus_20_40_2_3": 749.30,
        "rus_20_40_3_1": 825.86,
        "rus_20_40_3_2": 839.01,
        "rus_20_40_3_3": 841.92,
        "rus_20_40_4_1": 931.45,
        "rus_20_40_4_2": 935.72,
        "rus_20_40_4_3": 853.26,
        "rus_20_40_5_1": 875.91,
        "rus_20_40_5_2": 839.66,
        "rus_20_40_5_3": 881.95,
        "rus_20_40_6_1": 723.87,
        "rus_20_40_6_2": 982.42,
        "rus_20_40_6_3": 966.71,
        "rus_20_40_7_1": 752.12,
        "rus_20_40_7_2": 756.19,
        "rus_20_40_7_3": 869.75,
        "rus_20_40_8_1": 726.88,
        "rus_20_40_8_2": 898.68,
        "rus_20_40_8_3": 703.83,
        "rus_20_40_9_1": 866.61,
        "rus_20_40_9_2": 834.53,
        "rus_20_40_9_3": 864.01,
        #break
        "rus2_20_40_0_1": 106.92,
        "rus2_20_40_0_2": 120.62,
        "rus2_20_40_0_3": 88.20,
        "rus2_20_40_1_1": 107.43,
        "rus2_20_40_1_2": 104.48,
        "rus2_20_40_1_3": 101.62,
        "rus2_20_40_2_1": 119.26,
        "rus2_20_40_2_2": 118.57,
        "rus2_20_40_2_3": 112.44,
        "rus2_20_40_3_1": 110.14,
        "rus2_20_40_3_2": 141.20,
        "rus2_20_40_3_3": 123.54,
        "rus2_20_40_4_1": 99.16,
        "rus2_20_40_4_2": 107.05,
        "rus2_20_40_4_3": 95.18,
        "rus2_20_40_5_1": 95.46,
        "rus2_20_40_5_2": 88.32,
        "rus2_20_40_5_3": 169.89,
        "rus2_20_40_6_1": 119.14,
        "rus2_20_40_6_2": 68.35,
        "rus2_20_40_6_3": 134.12,
        "rus2_20_40_7_1": 88.48,
        "rus2_20_40_7_2": 111.55,
        "rus2_20_40_7_3": 81.44,
        "rus2_20_40_8_1": 106.90,
        "rus2_20_40_8_2": 99.86,
        "rus2_20_40_8_3": 100.59,
        "rus2_20_40_9_1": 94.55,
        "rus2_20_40_9_2": 94.58,
        "rus2_20_40_9_3": 63.61
    }
    return Z_Dict.get(problem_name, 0)
def create_csv():
    for filename in no_ext:
        out_file = filename + '.out'
        csv_file = filename + '.csv'
        in_file = csv.reader(open(out_file,'rb'), delimiter ='\t')
        out_csv = csv.writer(open(csv_file, 'wb'))
        #column header
        out_csv.writerow(['Col1', 'Col2', 'Col3' ,'Col4' ,'Col5' ,'Col6', 'Col7', 'Col8', 'Col9', 'Col10', 'Col11', 'Col12', 'Col13', 'Col14', 'Col15'])
        out_csv.writerows(in_file)

def problem_sized():
    for filename in no_ext:
        problem = filename.split('.uai')[0]
        if findZ(problem):
            s_List.append(filename)
        else:
            l_List.append(filename)

def mean_errors_AO():
    for abstraction in AO_proper:
        number_of_prob = 0
        running_sum_1min = 0
        running_sum_20min = 0
        running_sum_60min = 0
        running_sum_nodes = 0
        number_of_probl = 0
        running_sum_1minl = 0
        running_sum_20minl = 0
        running_sum_60minl = 0
        running_sum_nodesl = 0
        for small in s_List:
            problem = small.split('.uai')[0]
            if abstraction in small and 'AO' in small:
                df = pd.read_csv(small + '.csv')
                lastRow = int(df.shape[0])
                running_sum_1min += abs(findZ(problem) - float('%.2f' % float(df.iloc[60,2])))
                running_sum_20min += abs(findZ(problem) - float('%.2f' % float(df.iloc[1201,2])))
                running_sum_60min += abs(findZ(problem) - float('%.2f' % float(df.iloc[lastRow - 2, 2])))
                running_sum_nodes += float('%.2f' % float(df.iloc[lastRow - 2, 11]))
                number_of_prob += 1
        if number_of_prob == 0:
            continue
        else:
            print(abstraction + ' times for AO_proper small problems' + '\n'
                'error at 1 min: ' + str(running_sum_1min/number_of_prob) + '\n'
                'error at 20 min: ' + str(running_sum_20min/number_of_prob) + '\n'
                'error at 60 min: ' + str(running_sum_60min/number_of_prob) + '\n'
                'avg # of nodes: ' + str(running_sum_nodes/number_of_prob) + '\n')
        for large in l_List:
            problem = small.split('.uai')[0]
            if abstraction in large and 'AO' in large:
                df = pd.read_csv(large + '.csv')
                lastRow = int(df.shape[0])
                running_sum_1minl += abs(float(df.iloc[0,9]) - float('%.2f' % float(df.iloc[60,2])))
                running_sum_20minl += abs(float(df.iloc[0,9]) - float('%.2f' % float(df.iloc[1201,2])))
                running_sum_60minl += abs(float(df.iloc[0,9]) - float('%.2f' % float(df.iloc[lastRow - 2, 2])))
                running_sum_nodesl += float('%.2f' % float(df.iloc[lastRow - 2, 11]))
                number_of_probl += 1
        if number_of_probl == 0:
            continue
        else:
            print(abstraction + ' times for AO_proper large problems' + '\n'
                'error at 1 min: ' + str(running_sum_1minl/number_of_probl) + '\n'
                'error at 20 min: ' + str(running_sum_20minl/number_of_probl) + '\n'
                'error at 60 min: ' + str(running_sum_60minl/number_of_probl) + '\n'
                'avg # of nodes: ' + str(running_sum_nodesl/number_of_probl) + '\n')
#def mean_errors_OR():
    for abstraction in OR_proper:
        number_of_prob = 0
        running_sum_1min = 0
        running_sum_20min = 0
        running_sum_60min = 0
        running_sum_nodes = 0
        number_of_probl = 0
        running_sum_1minl = 0
        running_sum_20minl = 0
        running_sum_60minl = 0
        running_sum_nodesl = 0
        for small in s_List:
            problem = small.split('.uai')[0]
            if abstraction in small and 'OR' in small:
                df = pd.read_csv(small + '.csv')
                lastRow = int(df.shape[0])
                running_sum_1min += abs(findZ(problem) - float('%.2f' % float(df.iloc[60,2])))
                running_sum_20min += abs(findZ(problem) - float('%.2f' % float(df.iloc[1201,2])))
                running_sum_60min += abs(findZ(problem) - float('%.2f' % float(df.iloc[lastRow - 2, 2])))
                running_sum_nodes += float('%.2f' % float(df.iloc[lastRow - 2, 11]))
                number_of_prob += 1
        if number_of_prob == 0:
            continue
        else:
            print(abstraction + ' times for OR_proper small problems' + '\n'
                'error at 1 min: ' + str(running_sum_1min/number_of_prob) + '\n'
                'error at 20 min: ' + str(running_sum_20min/number_of_prob) + '\n'
                'error at 60 min: ' + str(running_sum_60min/number_of_prob) + '\n'
                'avg # of nodes: ' + str(running_sum_nodes/number_of_prob) + '\n')
        for large in l_List:
            problem = small.split('.uai')[0]
            if abstraction in large and 'OR' in large:
                df = pd.read_csv(large + '.csv')
                lastRow = int(df.shape[0])
                running_sum_1minl += abs(float(df.iloc[0,9]) - float('%.2f' % float(df.iloc[60,2])))
                running_sum_20minl += abs(float(df.iloc[0,9]) - float('%.2f' % float(df.iloc[1201,2])))
                running_sum_60minl += abs(float(df.iloc[0,9]) - float('%.2f' % float(df.iloc[lastRow - 2, 2])))
                running_sum_nodesl += float('%.2f' % float(df.iloc[lastRow - 2, 11]))
                number_of_probl += 1
        if number_of_probl == 0:
            continue
        else:
            print(abstraction + ' times for OR_proper large problems' + '\n'
                'error at 1 min: ' + str(running_sum_1minl/number_of_probl) + '\n'
                'error at 20 min: ' + str(running_sum_20minl/number_of_probl) + '\n'
                'error at 60 min: ' + str(running_sum_60minl/number_of_probl) + '\n'
                'avg # of nodes: ' + str(running_sum_nodesl/number_of_probl) + '\n')

def main():
    create_csv()
    problem_sized()
    mean_errors_AO()
  #  mean_errors_OR()

main()