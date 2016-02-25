__Original Git repository for this code [here](https://bitbucket.org/flandes/ineqfreezeseconomy_montecarlo_mastereq/src). This repository is just a GitHub mirror. Please go there for up-to-date code.__

Codes for: When does inequality freeze an economy?
=================================================

This repository contains code originally implemented for the paper
*"When does inequality freezes an economy?"*,
João Pedro Jericó, François P. Landes, Matteo Marsili, Isaac Pérez Castillo and Valerio Volpati. http://arxiv.org/abs/1602.07300
**Please cite when appropriate**.  

*C program: gran_econ, a sequential Monte-Carlo simulation*  

gran_econ.c is a Monte Carlo simulation of the trades of indivisible (granular) goods, between agents whom have a constant, Pareto-distriubted capital. Trading happens randomly (see paper) with the constraint that each agent's cash must be positive at all times.
Copyright (C) 2015 João Pedro Jericó, François P. Landes, Valerio Volpati.


*Mathematic script: iteSol_MM=2.nb, an iterative integration scheme*

iteSol_MM=2.nb is an iterative integration scheme that provides the self-consistent solution to Eqs. (11) and (4). See appendix A.3.2, "Algorithm computing self-consistent solution ...". It could be re-written in any formal mathematics program, including scilab, as it does not make use of any Mathematica-specific tools. 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.  
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  
You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.


Dependencies
------------

The simulation was written in C and requires the GSL library.

The plotting scripts were written in Python and require the Numpy, Matplotlib and Pandas libraries, but you may also make your own.

The Mathematica script requires Mathematica.

How to Use
----------

Please see how_to_use.md for details.

