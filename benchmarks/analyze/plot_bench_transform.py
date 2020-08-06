# -*- coding: utf-8 -*-
# @Author: Jose Tascon
# @Date:   2020-06-14 19:00:30
# @Last Modified by:   Jose Tascon
# @Last Modified time: 2020-08-04 00:08:12

import os
import numpy as np
import matplotlib.pyplot as plt

file = './logs/transform_interpolation_compare_float.log'
f = open(file, 'r')
count = 0
item_list = []


for x in f:
    if count>2:
        # print(x)
        line_list = str.split(x)
        func_label = line_list[0]
        func_split = func_label.split('/') 
        func_name = func_split[0]
        exp_name = func_split[1]
        time = line_list[1]
        units = line_list[2]
        # print(line_list)
        if(not ('_' in exp_name) and not ('3d' in func_name)):
            item_list.append([func_name, exp_name, time, units])
            # print('{:20}{:18}{:>15}{:^15}'.format(func_name, exp_name, time, units))
    count += 1
f.close()

a = np.array(item_list)
print(a)
print(a.shape)

funcs = {}
idxf = []

exps = {}
idx = []

# Functions Benchmarks
for k in range(a.shape[0]):
    s = a[k,0]
    if(funcs.get(s) == None):
        funcs[s] = 1;
        idxf.append(k)
    else:
        funcs[s] = funcs.get(s)+1

# Experiment argument
for k in range(a.shape[0]):
    s = a[k,1]
    d = a[k,0]
    if(list(funcs.keys())[0] == d):
        if (exps.get(s) == None):
            exps[s] = 1;
            idx.append(k)
        else:
            exps[s] = exps.get(s)+1

print(funcs)
print(idxf)

print(exps)
print(idx)

m = np.array(list(exps.values()))
data = np.zeros((np.max(m),len(exps.keys())*len(funcs.keys())))

ll = len(exps.keys())
for j,f in enumerate(funcs.keys()):
    for i,k in enumerate(exps.keys()):
        print(f,i,k)
        b = (a[idxf[j]+idx[i]:idxf[j]+idx[i]+exps.get(k),2])
        print(b)
        data[0:exps.get(k),ll*j + i] = b.astype(float)*(10e-6)
# print(data)

plt.figure(figsize=(7,4))
# plt.boxplot(data)
# plt.boxplot(data, labels=exps.keys())
# plt.title(os.path.basename(file))


labelsm=['8','80','800','8000','8','80','800','8000','8','80','800','8000']
plt.boxplot(data, labels=labelsm, )
plt.title('Benchmark of transform and interpolation 2D', fontsize=16)
plt.ylabel('Time [ms]', fontsize=15)
plt.yscale('log')
plt.xlabel('Image size $^2$', fontsize=15)
# plt.rc('xtick',labelsize=14)
# plt.rc('ytick',labelsize=18)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
