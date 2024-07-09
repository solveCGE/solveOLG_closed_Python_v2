#------------------------------------------------------------------------------#
# "solveOLG_closed_Python2"
# 
# Solves an AK-OLG-Model, closed economy, with income effects in Python
# Philip Schuster, July, 2023
#
# run.py: main script, shock definition, then computes transition path
#------------------------------------------------------------------------------#

import numpy as np
from numpy import sum, linspace, arange, array, kron, append, concatenate, mean, min, max
from scipy.optimize import fsolve
from math import floor, log, inf, isfinite
from timeit import default_timer
import copy
import matplotlib.pyplot as plt

from numba import jit
from loadfunctions import *

print("\nSIMPLE AUERBACH-KOTLIKOFF CLOSED ECONOMY MODEL IN PYTHON\n")

# load main functions
exec(open("algo.py").read())
exec(open("hh.py").read())
exec(open("firm.py").read())

# control center
tend            = 300   # number of periods
nag             = 100   # number of age groups (nag = x => max age = x-1)
budget_bal      = 3     # budget closing instrument (1.. tauW, 2.. tauF, 3.. tauC, 4.. taul, 5.. tauprof, 6.. cG) to fulfill given debt path
genplots        = False  # generate plots

# initializes all variables globally
exec(open("initdata.py").read())

# run calibration routine
exec(open("calib.py").read())   # <- change parameters in here

##########################
# POLICY SHOCK SECTION
# ========================
# Note: it is typically easiest to introduce shocks to period-view variables 
#       and then convert them to cohort-view variables using per2coh()

## tax shocks
#tauprof = tauprof*0+0.15                    # profit tax rate is set to 15%
#tauprof[9:tend] = tauprof[9:tend]*0+0.15    # profit tax rate is set to 15%, starting in 10 years (announced today)
#tauprof[0:10] = tauprof[0:10]*0+0.15        # profit tax rate is set to 15% temporarily for next 10 years
#tauWv = tauWv*1.02; tauWz = per2coh(tauWv)  # wage tax is increased by 2%

## delayed retirement (uncomment whole block)
# rag[0:10] = linspace(rag0,rag0+2,10); rag[10:tend]=rag0+2
# notretv[:,:] = 0
# for tt in range(0,tend):
#   notretv[0:floor(rag[tt]),tt] = 1
#   notretv[floor(rag[tt]),tt]   = rag[tt]-floor(rag[tt])
# notretz = per2coh(notretv)                    # effective retirement age is increased linearly by 2 years over the next 10 years
 
## pension cut
#pv = pv*0.95; pz = per2coh(pv) # pensions are cut by 5%

## productivity shocks
#thetav = thetav*1.02; thetaz = per2coh(thetav);                      # individual productivity increases permanently by 2%
#thetav[:,0:30] = thetav[:,0:30]*1.02; thetaz = per2coh(thetav)       # individual productivity increases by 2% in the next 30 years
#TFP = TFP*1.02                                                       # total factor productivity increases permanently by 2%

## fertility shocks
#NB = NB*1.02             # 2% more newborns every year
NB[0:30] = NB[0:30]*1.02 # 2% more newborns every year over next 50 years

## mortality shocks
gamv[59:nag,:] = 1-(1-gamv[59:nag,:])*0.9; gamv[nag-1,:]=0; gamz = per2coh(gamv)  # reduction of old-age mortality by 10%

## shock to the initial capital stock
#K[0] = K0*0.99                                                       # 1% of capital stock is lost in first period

## change in targeted debt path
#DG[4:20] = linspace(DG0,DG0*0.9,20-5+1); DG[20:tend] = DG0*0.9 # reduce public debt by 10% over 15 years starting in 5 years

##########################

# Solve transition path to new steady state
solveOLG(starttime = 1, maxiter = 200, tol = 1e-4)

# some transition plots
if genplots:

  plt.figure()
  plt.plot(range(0,tend+1),append(r0,r))
  plt.xlabel("time")
  plt.ylabel("real interest rate")
  plt.show()
  
  plt.figure()
  plt.plot(range(0,tend+1),append(w0,w))
  plt.xlabel("time")
  plt.ylabel("wage rate")
  plt.show()

  y_min = min(concatenate((N/N0,Y/Y0,Inv/Inv0,Cons/Cons0,CG/CG0,A/A0,P/P0)))
  y_max = max(concatenate((N/N0,Y/Y0,Inv/Inv0,Cons/Cons0,CG/CG0,A/A0,P/P0)))
  
  plt.figure()
  plt.plot(range(0,tend+1),append(1,N/N0)*100)
  plt.plot(range(0,tend+1),append(1,Y/Y0)*100)
  plt.plot(range(0,tend+1),append(1,Inv/Inv0)*100)
  plt.plot(range(0,tend+1),append(1,Cons/Cons0)*100)
  plt.plot(range(0,tend+1),append(1,CG/CG0)*100)
  plt.plot(range(0,tend+1),append(1,A/A0)*100)
  plt.plot(range(0,tend+1),append(1,P/P0)*100)
  plt.xlabel("time")
  plt.ylabel("period 0 = 100")
  plt.ylim(y_min*100,y_max*100)
  plt.legend(loc=1,labels=["population","GDP","investment","consumption","public consumption","aggregate assets","pension expend."])
  plt.show()
