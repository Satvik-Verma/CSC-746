#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: wes
Created: Thu Sep 30 05:51:28 PDT 2021

Description: this code generates a 2D "heatmap" style plot using sample data that
is hard-coded into the code.

Inputs: none, all problem parameters are hard-coded.

Outputs: a plot showing the heatmap, displayed to the screen

Dependencies: matplotlib, numpy

Assumptions: Developed and Tested with Python 3.8.8 on MacOS 11.6
'''

import numpy as np
import matplotlib.pyplot as plt

threads_per_block = ['32', '64', '128', '256', '512', '1024'] # y axis, 6 of them
thread_blocks = ["1", "4", "16", "64", "256", "1024", "4096"] # x axis, 7 of them

runtime = np.array([
    [1.56, 1.56, 1.56, 1.56, 3.70, 14.77, 32.74],
    [3.12, 3.12, 3.12, 3.12, 7.41, 29.07, 70.32],
    [6.25, 6.25, 6.25, 6.25, 14.79, 55.03, 85.00],
    [12.50, 12.50, 12.50, 12.49, 29.20, 71.16, 91.48],
    [24.98, 24.99, 24.80, 24.79, 55.99, 84.87, 92.66],
    [49.93, 48.83, 48.69, 48.57, 84.27, 91.05, 91.41]
])


fig, ax = plt.subplots()
im = ax.imshow(runtime, cmap="coolwarm")

# We want to show all ticks...
ax.set_xticks(np.arange(len(thread_blocks)))
ax.set_yticks(np.arange(len(threads_per_block)))
# ... and label them with the respective list entries
ax.set_xticklabels(thread_blocks)
ax.set_yticklabels(threads_per_block)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
for i in range(len(threads_per_block)): # y axis
    for j in range(len(thread_blocks)): # x axis
        text = ax.text(j, i, runtime[i, j],
                       ha="center", va="center", color="k")

ax.set_title("Achieved Occupancy on CUDA-GPU at Varying Block Size and Number of Blocks")
ax.set_ylabel('Threads per block')
ax.set_xlabel('Block Sizes')
fig.colorbar(im, ax=ax)
fig.tight_layout()
plt.show()

# EOF
