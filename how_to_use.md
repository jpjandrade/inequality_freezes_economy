Granular Economies code
=======================

How to compile
--------------

simply type:
make
in the folder containing gran_econ.c, gran_econ.h and Makefile, and the code will be compiled into  "gran_econ"


How to run the simulation
-------------------------

The code should be run with 5 arguments:

```
./gran_econ N beta e flag_thermalized flag_hist {optional: M_user}
```

e.g.

./gran_econ 10000  1.4 1.2  1  0 

or:

./gran_econ 100000 1.2 2.0  0  0 


__N__ is the number of agents in the economy. We call attention to the fact that for cluster considerations running time is fixed inside the code via the T_TRIAL_STEPS argument. Therefore, one must be careful not to set N too large that it will not thermalize (forget the initial condition). As a notion of lenght, N ~ 10^5 is good for codes that run around 2h for beta not too close to 1

__beta__ is the exponent for the Pareto distribution of initial capitals. This will directly impact the time it takes for convergence. For beta > 2 it will converge orders of magnitude faster than beta close to 1. For beta less than 1, the average capital is ill-defined and every run will yirld a priori very different results. 

__e__ is the capital / goods ratio (ratio C/\Pi in the paper), ie, how many goods to generate in proportion to the capital alocated. Anything from 1.2 to 2 is fine, but this will have an impact on the p^suc as a function of beta (ie, choose one and stick with it to compare results). All the figures in the paper were done with _e = 1.2_

__flag_thermalized__ should be either 0 or 1, and is whether the code should consider the results thermalized or not. This will impact a) how much of the simulation it should regarded as thermalized and record statistics and b) whether the code will look for files of a previous run (see M_user entry)

__flag_hist__ should be either 0 or 1, and will dictate whether the code will keep track of good ownerships and generate histograms such as in the insets of figure 4 of the paper (blue dots). This slows down the code considerably, so we recommend setting it 1 only when the figure is desired.

__M_user__ should be equivalent to the number of goods M generated in a previous run, for which there are _capital_M=M_user.dat_ and _goods_M=M_user.dat_ files in the same directory. The code will look for those to generate initial capital distribution and good alocation, allowing one to resume a run from where it left off. Useful for low betas where thermalization requires a long time, but can be disregarded most of the time. Requires flag_thermalized to be 1.

Outputs
-------

The code output is the following files:

- capital_M.dat
- goods_M.dat
- ownership_averages_{args}.dat
- ps_{args}.dat
- ps_time_series_{args}.dat

Where {args} is a list of arguments used, as for example _ps_K=1_N=100_beta=1.80_e=1.20.dat_.

__capital_M.dat__ and __goods_M.dat__ are the capital distribution and final goods allocation after the simulation is done. They serve in case one wants to run the simulation again from where it left. For it's usage, see __M_user__ in the __How to run__ section.

__ownership_averages_{args}.dat__ is a file with N rows, each with K+1 columns corresponding to the capital and average ownership of one agent for each type of good. So the lines look like this:

```
...
23.429355 544.131668 58.232979 8.745685 0.008059 
23.367232 544.402810 58.228270 8.709992 0.008805 
23.270794 544.230535 58.181916 8.685585 0.008111
...
```

This means the consumer with capital 23.270794 had on average 544 goods of type 1, 58 of type 2, etc.

__ps_{args}.dat__ contains the final K values of average p_suc plus some extra auxiliary quantities: beta, e, N, pi_0, G and total_capital. For example, in this sample output:

```
0.5836556011 0.3556625681 0.1523882596 0.0482661422 0.0171140729  1.400000 1.200000 1000 0.00500 2.0 3876.416016 
```

There are 5 types of goods, and the final p_suc probabilities are given by: 0.58, 0.35, ..., 0.02. beta is 1.4, e is 1.2, N is 1000, the initial price of good 0 is 0.005, each subsequent good is 2.0 times the price of the previous one and the total capital of all the agents (samples from a Pareto distribution with beta = 1.4) is 3876.42.

__ps_time_series_{args}.dat__

This is a histogram for the evolution of p_suc along the thermalization window simulation, during which statistics is collected. It should be used as a guarantee of convergence, specially for lower betas. As long as this is not stationary, the thermalization window should be extended. If it only midly fluctuates, you are typivally in the steady-state regime. Each row is simply the p_suc for the current Monte Carlo step (1 MC step = M random trades) and they look like this:

```
 0.931962 0.875309 0.777982 0.604457 0.384932 
 0.936982 0.879936 0.774672 0.601598 0.382448 
 0.934086 0.874522 0.780484 0.612656 0.398350 
```


Plotting scripts
----------------

The figures generated from the simulation results were done by the following Python scripts.

__plot_gini.py__

__plot_liquidity.py__


Both of these scripts will look for all the _ownership_averages_{args}.dat_ files in the current directory and generate the Gini figure (figure 1 of the paper) and the liquidity figure (figure 2, right)

If one wants to put the data in another directory, the `files_list = glob.glob('own*')` line must be changed to something appropriate, for ex `files_list = glob.glob('data/own*')`.

Also, the scripts don't differentiate different parameters, so if you have several files from different runs with e, N and other things changed, they will plot them all together!

__plot_ps.py__

This will look at all the _ownership_averages_{args}.dat_, _ps_{args}.dat_ and _ps_time_series_{args}.dat_ in the current directory and generate the p_suc file (figure 2, left) and the average ownership plots (figure 3, main).

The same caveats as the former script apply: if the data is not on the same directory as the script, the code has to be changed. Also, don't mix different parameters!

The number of goods has to be set by hand, inside the code, in the line `n_goods = 5`. If there are several files in the directory with different number of goods, the code will fail.

