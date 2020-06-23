# -*- coding: utf-8 -*-
# @Author: Jose Tascon
# @Date:   2020-06-14 19:00:30
# @Last Modified by:   Jose Tascon
# @Last Modified time: 2020-06-14 20:20:52

import os
import numpy as np
import matplotlib.pyplot as plt

file = './logs/mem_cpu2gpu.log'
file = './logs/mem_gpu2cpu.log'
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
        if(not ('_' in exp_name)):
            item_list.append([func_name, exp_name, time, units])
            # print('{:20}{:18}{:>15}{:^15}'.format(func_name, exp_name, time, units))
    count += 1
f.close()

a = np.array(item_list)
print(a)
print(a.shape)

exps = {}
idx = []

for k in range(a.shape[0]):
    s = a[k,1]
    if(exps.get(s) == None):
        exps[s] = 1;
        idx.append(k)
    else:
        exps[s] = exps.get(s)+1

print(exps)
print(idx)

m = np.array(list(exps.values()))
data = np.zeros((np.max(m),len(exps.keys())))

for i,k in enumerate(exps.keys()):
    print(i,k)
    b = (a[idx[i]:idx[i]+exps.get(k),2])
    print(b)
    data[0:exps.get(k),i] = b.astype(float)*(10e-6)
print(data)

plt.figure(figsize=(6,4))
plt.boxplot(data, labels=exps.keys())
plt.title(os.path.basename(file))
plt.ylabel('Time [ms]')
plt.xlabel('Vector size')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
