#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Copyright (C) 2015 João Pedro Jericó.

from numpy import *
from matplotlib.pyplot import *
import os
import re
import sys
import glob
import pandas as pd

own_files = glob.glob('own*')
ps_series_files = glob.glob('ps_time_series*')
ps_files = glob.glob('ps_*')

fig1 = figure(figsize=(10, 6))
fig2 = figure()

n_goods = 1

ps_full_data = zeros(n_goods + 2)
for files in zip(own_files, ps_series_files, ps_files):

    own_file = re.search('own\S*beta=(\d\.\d\d)_e=(\d\.\d\d)\S*.dat', files[0])
    ps_series_file = re.search('ps_time_series_\S*beta=(\d\.\d\d)_e=(\d\.\d\d)\S*.dat', files[1])
    ps_file = re.search('ps_\S*beta=(\d\.\d\d)_e=(\d\.\d\d)\S*.dat', files[2])
    beta = own_file.group(1)
    ee = own_file.group(2)
    print('matched files:')
    print(own_file.group(0))
    print(ps_series_file.group(0))
    print('beta = ', beta, ee)

    own_data = np.loadtxt(own_file.group(0))
    ps_ts_data = np.loadtxt(ps_series_file.group(0))
    ps_data = np.loadtxt(ps_file.group(0))

    ps_full_data = np.vstack((ps_full_data, ps_data[:-4]))
    
    own_fig_name = 'ownership_beta='+beta+'_e='+ee
    ps_ts_fig_name = 'ps_ts_beta='+beta+'_e='+ee
    n_goods = len(own_data[1]) - 1

    ax2 = fig2.add_subplot(111)
    for i in range(n_goods):

        ax1 = fig1.add_subplot(n_goods,1,i+1) 
        ax1.semilogx(own_data[:,0], own_data[:,i+1], 'o', color = 'blue')        
        ax1.set_xlabel('$c$')
        ax1.set_ylabel('<Z> of good %d' % (i+1))

        ax2.plot(ps_ts_data, 'o-', markersize = 2.)
        ax2.set_xlabel('$t$')
        ax2.set_ylabel('$p_s$')
        
    fig1.savefig(own_fig_name+'.png', bbox_inches = 'tight')
    fig2.savefig(ps_ts_fig_name+'.png', bbox_inches = 'tight')
    fig1.clf()
    fig2.clf()
ps_full_data = ps_full_data[1:]

plot_data = pd.DataFrame(ps_full_data)
len_fd = len(ps_full_data)
df = plot_data.groupby([n_goods, n_goods+1]).agg([mean, std, min, max])
df.index.names = ['beta', 'e']
df = df.reset_index()
df.fillna(0, inplace = True)

#plot main figure
fig3 = figure(figsize=(10,6))
ax3 = fig3.add_subplot(111)
for i in range(n_goods):
    ax3.errorbar(df['beta'], df[i]['mean'], yerr = df[i]['std'], fmt = '--o')
ax3.axvline(x=1, color = 'black', linestyle = '--')
ax3.set_xlabel(r'$\beta$', size = 26)
ax3.set_ylabel(r'$p_s$', size = 26)
ax3.set_ylim([0,1])
fig3.savefig('ps_plot.png', bbox_inches = 'tight')
#fig3.savefig('ps_plot.pdf', bbox_inches = 'tight', format = 'pdf')    
    
