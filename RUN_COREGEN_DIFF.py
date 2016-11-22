import numpy as np
import math
import os
import PLOTS
import EVOLVE
import EOS
from CONST import RE, ME, G

####################################################
#This run is designed to generate a core given
#the profile at the interface, assuming that the
#core density is higher at the interface
####################################################

####################################################
#INPUT PARAMETERS
####################################################

base = '/projects/'																	#projects directory
fin  = base+'MESA/Kepler11/d/'
fin  = fin+'newprofile.data'												#input profile
fdir = base+'Planet_Profiles/Kepler11/d2/'
fout = fdir+'Run68_Kepler11d2.txt'         					#where to save output profile

h    = -10e-6*RE                                    #Step size (in Earth Radii); negative means integrating backwards

#Interface parameters
gam_c= 5.                                           #equation of state for the core


####################################################
#LEGEND
####################################################

#EoS - equations of state:
#Seager(component,T0,P0), 0-Fe, 1-MgSiO3, 2-H2O
#Salpeter
#Tillotson
#Polytrope(gamma,P0,d0)

####################################################
#READ INTERFACE PROFILE
####################################################
if not os.path.exists(fdir):
  os.makedirs(fdir)

M_E,R_E,P_E,d_E,T_E,mu_E,u_E = np.loadtxt(fin,usecols=(0,1,2,3,4,5,6),unpack=True)

#Assign boundary conditions for the core
M0 = M_E[len(M_E)-1]
R0 = R_E[len(M_E)-1]
P0 = P_E[len(M_E)-1]
d0 = d_E[len(M_E)-1]
T0 = 10000.
dmin = d_E[len(M_E)-1]

####################################################
#Equation of State of the Core
####################################################
dmax = 20.                                          #Guess for the core density at the interface
Ps = np.zeros(1)                                      #arrays containing mass fraction to switch EoS


####################################################
#Searches for Solution
####################################################
EoS = [EOS.Seager(1,P0),EOS.Seager(0,P0)]

#Check if heaviest cases is valid
print 'checking heaviest case'
Ps[0] = 1.
PROFILE, density_high = EVOLVE.EVOLVE(h,EoS,Ps,P0,R0,M0,T0)
print min(PROFILE.r)/RE, min(PROFILE.m)/ME, density_high, density_high/M0
if (density_high < 0.):
  print 'pure iron is too light'
  quit()

#Check if lightest case is valid
print 'checking lightest case'
Ps[0] = 0.
PROFILE, density_high = EVOLVE.EVOLVE(h,EoS,Ps,P0,R0,M0,T0)
print min(PROFILE.r)/RE, min(PROFILE.m)/ME, density_high, density_high/M0
if (density_high > 0.):
  print 'pure silicate is too low'
  quit()

#Solve for core solution and output final plots and profile
efficiency = math.pow(10.,-4.)*100.
count = 0
Ps_min = 0.
Ps_max = 1.
while efficiency > 5.*math.pow(10.,-6.):
  Ps[0] = 0.5*(Ps_min+Ps_max)
  PROFILE, dm = EVOLVE.EVOLVE(h,EoS,Ps,P0,R0,M0,T0)

  efficiency = math.fabs(dm/M0)
  count += 1
  print count, efficiency, Ps[0], min(PROFILE.r)/RE, min(PROFILE.m)/ME, dm/M0

  if min(PROFILE.r)/RE >=math.fabs(h)*10.:
    Ps_max = Ps[0]
    efficiency = math.pow(10.,-4)*100.
  else:
    print 'radius is less than 10xh'
    if (dm < 0):
      Ps_min = Ps[0]
      efficiency = math.pow(10.,-4)*100.
    if (dm >= 0):
      print 'dm > 0', efficiency
      Ps_max = Ps[0]
  #dm is the density*volume - mass left of innermost point; cannot fit excess mass

PLOTS.TWO_COMPONENT_DIFF(math.fabs(h),fout,fdir,M_E,R_E,P_E,d_E,T_E,mu_E,u_E,gam_c,PROFILE,Ps[0])







