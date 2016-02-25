#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Copyright (C) 2015 João Pedro Jericó.

from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt #c est la meme chose que plt
import os
import re
import sys
import glob

plt.ion()
############available in matplotlib 1.1.rc
params = {'backend': 'pdf',
          'axes.labelsize': 25,
          'axes.labelweight': 'roman',
          'axes.titlesize': 16,
          'text.fontsize': 12,
          'legend.fontsize': 16,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
         'font.serif': ['computer modern roman'],
          'text.usetex': True,
          'font.family': 'serif'
            }
          
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True          
plt.rcParams.update(params)





rcParams['text.latex.unicode']=True
G = 1.5
pi_0 = 0.005

file_list = glob.glob('own*.dat')


def get_prices(own_data):
    n_goods = own_data.shape[-1] - 1
    exponents = array(range(n_goods))
    prices = pi_0 * G ** exponents
    return prices

fig = figure(1,[10,6])
ax = fig.add_subplot(111)

beta_list = arange(1.1, 2.6, 0.1)
#beta_list = append(beta_list, array([3.,6.,10.,20.,100., 1000.]))
#cmap = cm.jet  # sort of rainbow
cmap = cm.cool_r # shades of blue(s)
splitter=1	#use larger split to have very different colors 
gradient=cmap(np.linspace(0.0,1.0, 15))

    
for own_file in file_list:
    own_params = re.search('own\S*beta=(\d+\.\d\d)_e=(\d\.\d\d)\S*.dat', own_file)
    beta = float(own_params.group(1))
    ee = own_params.group(2)
    print('matched files:')
    print(own_params.group(0))
    print('beta = ', beta, 'e =', ee)
    i_beta = argwhere(abs(beta_list - beta) < 0.01).flatten()
    i_beta = i_beta[0]
    own_data = np.loadtxt(own_file)
    cash_left = own_data[:,0] - own_data[:,1:].dot(get_prices(own_data))
    n_goods = own_data.shape[-1] - 1
    if(float(beta) < min(beta_list)+0.02):
        ax.plot(own_data[:,0], cash_left, linestyle = '-', color = tuple(gradient[(i_beta+0.)*splitter]), label = r'$\beta = %.1f$' % (beta))            
        ax.plot(own_data[:,0], own_data[:,0], 'r--')
    elif(float(beta) > max(beta_list)-0.02):
        ax.plot(own_data[:,0], cash_left, linestyle = '-', color = tuple(gradient[(i_beta+0.)*splitter]), label = r'$\beta = %.1f$' % (beta))            
    else:
        ax.plot(own_data[:,0], cash_left, linestyle = '-', color = tuple(gradient[(i_beta+0.)*splitter]))
    
        




#ax2.loglog(data[:,0], data[:,0], 'r--')
for i in range(n_goods):
    ax.axhline(y = get_prices(own_data), color = 'black', linestyle = '--')
ax.axhline(y = 1, color = 'green', linestyle = '--')
ax.set_xlabel('$c$')
ax.set_ylabel('$l$')
ax.set_xscale('log')
ax.set_yscale('log')
#leg = ax.legend(loc = 'upper left', frameon = False)
fig.savefig('liquidity_cluster_nolog.png', bbox_inches = 'tight')

show()
