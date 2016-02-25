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


def calc_gini(x): #follow transformed formula
    """
    Return computed Gini coefficient.
    """
    xsort = sorted(x) # increasing order
    y = np.cumsum(xsort)
    B = sum(y) / (y[-1] * len(x))
    return 1 + 1./len(x) - 2*B

def get_prices(own_data):
    n_goods = own_data.shape[-1] - 1
    exponents = array(range(n_goods))
    prices = pi_0 * G ** exponents
    return prices

fig = figure(1,[10,6])
ax = fig.add_subplot(111)

beta_list = arange(1.1, 2.6, 0.1)
beta_list = append(beta_list, array([3.,6.,10.,20.,100., 1000.]))
#cmap = cm.jet  # sort of rainbow
cmap = cm.cool_r # shades of blue(s)
splitter=1	#use larger split to have very different colors 
gradient=cmap(np.linspace(0.0,1.0, 15))

gini_cap = []
gini_liq = []

files_list = glob.glob('own*')
for own_file in files_list:
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
    gini_cap.append(calc_gini(own_data[:,0]))
    gini_liq.append(calc_gini(cash_left))
    if beta == 1.1 or beta == 2.5 or beta == 6.0:
      ax.text(gini_cap[-1]+0.01, gini_liq[-1]-0.05, r'$\beta$ = %.1f' % beta)

gini_cap = np.array(gini_cap)
gini_liq = np.array(gini_liq)
order = np.argsort(gini_cap)
gini_cap = gini_cap[order]
gini_liq = gini_liq[order]
        
x = arange(0,1, 0.01)
ax.spines['top'].set_color('none'); ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom'); ax.yaxis.set_ticks_position('left')
ax.plot(x, x, 'k--')
ax.plot(gini_cap, gini_liq, 'r-o', markersize = 7)
ax.axhline(y = 1, linestyle = '--', alpha = 0.5, color = 'black')
ax.set_xlim([0.,1.])
ax.set_ylim([0.,1.])
ax.set_xlabel('$G_c$')
ax.set_ylabel('$G_l$')
fig.savefig('gini_coef.png', bbox_inches = 'tight')
