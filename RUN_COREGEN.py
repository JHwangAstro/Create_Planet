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
fin  = base+'Kepler11/d/'
fin  = fin+'newprofile.data'												#input profile
fdir = base+'Planet_Profiles/Kepler11/d/'           #where to save output profile
fout = fdir+'Run54_Kepler11d.txt'          					#where to save output profile

h    = -10e-8*RE                                    #Step size (in Earth Radii); negative means integrating backwards

#Interface parameters
q    = 100.                                         #parameter in Sigmoid Function
rhop = 0.4                                          #parameter in Sigmoid Function
rhoc = 8.0                                          #guess for density of core
gam_c= 5.                                           #equation of state for the core
gam_e= 4./3.                                        #equation of state for the core


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
Ps   = [2./3.]                                      #arrays containing mass fraction to switch EoS


####################################################
#Searches for Solution
####################################################
#Check if dmax is valid
EoS = [EOS.Polytrope(gam_c,P0,dmax)]              #Flag for Equation of State
#EoS = [EOS.Seager(1,dmax),EOS.Seager(0,dmax)]
PROFILE, density_high = EVOLVE.EVOLVE(h,EoS,Ps,P0,R0,M0,T0)
print min(PROFILE.r)/RE, min(PROFILE.m)/ME, density_high
if (density_high < 0.):
  print 'maximum density guess is too low'
  quit()

#Solve for core solution and output final plots and profile
efficiency = math.pow(10.,-4.)*100.
count = 0
while efficiency > math.pow(10.,-5.):
  d0 = 0.5*(dmin+dmax)
  EoS  = [EOS.Polytrope(gam_c,P0,d0)]              #Flag for Equation of State
  PROFILE, dm = EVOLVE.EVOLVE(h,EoS,Ps,P0,R0,M0,T0)

  efficiency = math.sqrt(math.pow(min(PROFILE.r)/RE,2.)+math.pow(min(PROFILE.m)/ME,2.))
  count += 1
  print count, efficiency, d0, min(PROFILE.r)/RE, min(PROFILE.m)/ME, dm/M0
  if (dm < 0):
    dmin = d0
    efficiency = math.pow(10.,-4)*100.
  if (dm >= 0):dmax = d0
  #dm is the density*volume - mass left of innermost point; cannot fit excess mass

  #PROFILE.fname = PROFILE.fname + '_' + str(count)
  #PLOTS.TWO_COMPONENT(math.fabs(h),fdir,M_E,R_E,P_E,d_E,PROFILE)

PLOTS.TWO_COMPONENT(math.fabs(h),fout,M_E,R_E,P_E,d_E,T_E,mu_E,u_E,gam_c,PROFILE)







