
import numpy as np
import math
from CONST import PI,G
from PROFILE import PROFILE

#Integrates radius, pressure, density, and mass
#Sets initial conditions at the center of the star
def EVOLVE(h,EoS,Ps,P0,R0,M0,T0):

  #x = radius
  #y = pressure
  #z = dPdr
  n = 0

  #Change in Pressure/Change in radius
  def dPdr(x,y,mass,tmp):
    rho = EoS[n].rho(y,tmp)
    return -G*mass*rho/x**2.

  #Second derivative of pressure
  #def d2pdr(x,y,z,tmp):
  #  rho = EoS.rho(y,tmp)
  #  drho= EoS.drho_dP(y,tmp)
  #  return (1./rho*drho*z**2.-2.*z/x-4.*PI*G*rho**2.)

  #Change in mass/change in radius
  def dmdr(x,y,tmp):
    rho = EoS[n].rho(y,tmp)
    return 4.*PI*x**2.*rho

  #Initial conditions
  r  = [R0]
  P  = [P0]
  m  = [M0]
  T  = [T0]
  #print n
  d  = [EoS[n].rho(P0,T0)]

  i = 1

  #Check boundary conditions
  while (P[i-1] > math.pow(10.,-10) and T[i-1] > 0. and d[i-1] > 0. and m[i-1] >= 0 and r[i-1] >= 0):

    k1 =  dPdr(r[i-1],P[i-1],m[i-1],T[i-1])
    #l1 = d2pdr(r[i-1],P[i-1],k1,T[i-1])
    m1 =  dmdr(r[i-1],P[i-1],T[i-1])
    t1 = EoS[n].dt_dr(P[i-1],k1,T[i-1])

    if P[i-1]+h*k1/2. < 0: break
    if T[i-1]+h*t1/2. < 0: break
    if m[i-1]+h*m1/2. < 0: break

    m2 =  dmdr(r[i-1]+h/2.,P[i-1]+h*k1/2.,T[i-1]+h*t1/2.)
    k2 =  dPdr(r[i-1]+h/2.,P[i-1]+h*k1/2.,m[i-1]+h*m1/2.,T[i-1]+h*t1/2.)
    #l2 = d2pdr(r[i-1]+h/2.,P[i-1]+h*k1/2.,k2,T[i-1]+h*t1/2.)
    t2 = EoS[n].dt_dr(P[i-1]+h*k1/2.,k2,T[i-1]+h*t1/2.)

    if P[i-1]+h*k2/2. < 0: break
    if T[i-1]+h*t2/2. < 0: break
    if m[i-1]+h*m2/2. < 0: break

    m3 =  dmdr(r[i-1]+h/2.,P[i-1]+h*k2/2.,T[i-1]+h*t2/2.)
    k3 =  dPdr(r[i-1]+h/2.,P[i-1]+h*k2/2.,m[i-1]+h*m2/2.,T[i-1]+h*t2/2.)
    #l3 = d2pdr(r[i-1]+h/2.,P[i-1]+h*k2/2.,k3,T[i-1]+h*t2/2.)
    t3 = EoS[n].dt_dr(P[i-1]+h*k2/2.,k3,T[i-1]+h*t2/2.)

    if P[i-1]+h*k3 < 0: break
    if T[i-1]+h*t3 < 0: break
    if m[i-1]+h*m3 < 0: break
    if r[i-1]+h <= 0: break

    m4 =  dmdr(r[i-1]+h,P[i-1]+h*k3,T[i-1]+h*t3)
    k4 =  dPdr(r[i-1]+h,P[i-1]+h*k3,m[i-1]+h*m3,T[i-1]+h*t3)
    #l4 = d2pdr(r[i-1]+h,P[i-1]+h*k3,k4,T[i-1]+h*t3)
    t4 = EoS[n].dt_dr(P[i-1]+h*k3,k4,T[i-1]+h*t3)

    if P[i-1] + (h/6.)*(k1 + 2.*k2 + 2.*k3 + k4) < 0.: break

    #Update radius, pressure, dP/dr, density, mass
    r.append(r[i-1] + h)
    P.append(P[i-1] + (h/6.)*(k1 + 2.*k2 + 2.*k3 + k4))
    m.append(m[i-1] + (h/6.)*(m1 + 2.*m2 + 2.*m3 + m4))
    T.append(T[i-1] + (h/6.)*(t1 + 2.*t2 + 2.*t3 + t4))
    d.append(EoS[n].rho(P[i],T[i]))
    #print i, P[i], T[i], d[i], r[i], m[i]#, dPdr[i]
    i = i + 1

    #check to see if condition is met to switch to different EoS
    if n < len(Ps):
      #print m[len(m)-1]/M0, Ps[n]
      if m[len(m)-1]/M0<Ps[n]:
        n = n+1
        print 'switching to iron'

  #density_high = 0
  k = len(r)-1
  #if 4.*PI*r[n]**3.*d[n]/3. > m[n]: density_high=1

  return PROFILE(EoS,np.array(r),np.array(m),np.array(d),np.array(P),np.array(T)), 4.*PI*r[k]**3.*d[k]/3.-m[k]















