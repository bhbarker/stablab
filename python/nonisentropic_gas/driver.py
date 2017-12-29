#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 09:00:21 2017

@author: blake
"""

# -----------------------------------------------------------------------------
# imports
# -----------------------------------------------------------------------------

import sys
sys.path.append('../bin/')

from sympy import symbols
from time_evolution import create_code

# -----------------------------------------------------------------------------
# setup
# -----------------------------------------------------------------------------

# parameters
nu, mu, eta, cnu, G = symbols('nu, mu, eta, cnu, G')

# system variables
rho, u, e = symbols('rho, u, e')

# vector containing the system variables
U = [rho,u,e]

# system parameters
parameters = [nu, mu, eta, cnu, G]

#
# f_0(U)_t + F_1(U)_x = (B(U)U_x)_x
#

f0 = [rho,rho*u,rho*(e+u**2/2)]

f1 = [rho*u, rho*u**2+G*rho*e, rho*u*(e+u**2/2)+G*rho*u*e]

BU = [[0,0,0],[0,(2*mu+eta),0],[0,(2*mu+eta)*u, nu]]


# file path to location where code will be written
file_path = 'FILE PATH'

create_code(U,parameters,f0,f1,BU,file_path)

