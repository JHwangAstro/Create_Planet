
import numpy as np
import math
from CONST import PI, G, kb
from scipy.optimize import fsolve

class Salpeter:

  def rho(y):

    gamma = np.zeros([4,5])
    gamma[0,:] = [ .01512,  .08955, .109,  5.089,-5.980]
    gamma[1,:] = [ .002181,-.4015, 1.698, -9.566, 9.873]
    gamma[2,:] = [-.0003328,.5167,-2.369, 13.49, -14.27]
    gamma[3,:] = [-.01384, -.652,  3.529,-20.95,  22.64]

    yo = 9.524*math.pow(10.,13)/10.
    xi = math.pow(y/yo,1./5.)
    epsilon = math.pow(3./(32.*PI**2.),1./3.)

    beta = np.zeros(4)
    for i in range(0,5):
      for j in range(0,5):
        beta[i] += gamma[i,j]*math.pow(epsilon,j*.5)
      beta[i] = math.pow(beta[i],i+1)

    betasum = 0.
    for i in range(1,6):
      betasum = betasum + beta[i-1]*math.pow(xi,i)

    alpha = 1.941*.01-math.sqrt(epsilon)*6.277*.01+epsilon*1.076
    phi = math.pow(3.,1./3.)/20.+epsilon/math.pow(4.3,1./3.)
    xo = (1.+math.exp(-alpha*xi))*betasum/(xi+phi)
    rho = 3.886*math.pow(xo,-3)*1000.

    return rho


  def drho_dP(y):

    gamma = np.zeros([4,4])
    gamma[0,:] = [ .01512,  .08955, .109,  5.089,-5.980]
    gamma[1,:] = [ .002181,-.4015, 1.698, -9.566, 9.873]
    gamma[2,:] = [-.0003328,.5167,-2.369, 13.49, -14.27]
    gamma[3,:] = [-.01384, -.652,  3.529,-20.95,  22.64]

    yo = 9.524*math.pow(10.,13)/10.
    xi = math.pow(y/yo,1./5.)
    epsilon = math.pow(3./(32.*PI**2.),1./3.)

    beta = np.zeros(4)
    for i in range(0,5):
      for j in range(0,5):
        beta[i] += gamma[i,j]*math.pow(epsilon,j*.5)
      beta[i] = math.pow(beta[i],i+1)

    alpha = 1.941*.01-math.sqrt(epsilon)*6.277*.01+epsilon*1.076
    phi = math.pow(3.,1./3.)/20.+epsilon/math.pow(4.3,1./3.)

    xo = (1.+math.exp(-alpha*xi))*betasum/(xi+phi)

    betasum = 0.
    for i in range(1,6):
      betasum = betasum + beta[i-1]*math.pow(xi,i)*(float(i)/xi-alpha-math.pow(xi+phi,-1.))

    drho = -3.*3.866*math.pow(xo,-4.)
    drho = drho*(1.+math.exp(-alpha*xi))*betasum/(xi+phi)
    drho = drho*0.2*math.pow(y/yo,-.8)/yo
    return drho


class Seager:

  def __init__(self,c,P0):
    #Fe
    if c == 0:
      self.n = 0.528
      self.rho_c = 8.3
      self.co = 3.49*math.pow(10.,-(self.n+6.))
      #self.co = 3.49*math.pow(10.,-6.*self.n)

    #MgSiO3
    if c == 1:
      self.n = 0.549#0.541
      self.rho_c = 4.26#4.1
      self.co = 1.27*math.pow(10.,-(self.n+6.))#1.61*math.pow(10.,-(self.n+6.))
      #self.co = 1.27*math.pow(10.,-6.*self.n)#1.61*math.pow(10.,-(self.n+6.))

    #H2O
    if c == 2:
#      self.rho_c = 1460.
#      self.co = 0.00311
      self.n = 0.513


  def rho(self,y,T):

    return self.rho_c + self.co*y**self.n


  def drho_dP(self,y,T):

    return self.co*self.n*y**(self.n-1.)

  def dt_dr(self,y,z,T):

    return 0.
#    return (self.K*self.alpha*rho**(self.alpha-1.)*drho*z)



class Polytrope:

  def __init__(self,gamma,P0,d0):
    self.gamma = gamma
    self.K = P0/math.pow(d0,gamma)

  def rho(self,y,T):
    return (y/self.K)**(1./self.gamma)

  def drho_dP(self,y,T):
    return (1./(self.gamma*self.K))*(y/self.K)**(1./self.gamma-1.)

  def dt_dr(self,y,z,T):
    return 0.


class Tillotson:

  def rho(c,y):
    #Fe
    if c == 0:
      def f(x):
        a = 0.5
        b = 1.5
        A = 1.270e11
        B = 1.05e11
        u = 0.
        u0 = 0.
        eta = x/x0
        mu = eta - 1.
        return (a+(b*u0*math.pow(eta,2))/(mu*u0*math.pow(eta,2)))*u*x+A*mu+B*math.pow(mu,2)



#Woolfson fitted to Stevenson & Salpeter 1976 EoS
class Woolfson:

  def __init__(self):
    self.c = 0.003 #m3/kg
    self.mu = 3.*math.pow(10.,-27) #kg

  def rho(self,y,T):
    return (math.sqrt(1.+4.*mu*c/(kb*T))-1.)/(2.*c)


  def drho_dP(self,y,T):
    return mu/(kb*T(1.+2.*c*self.rho(y,T)))


  def dt_dr(self,y,z,T):
    drhodP = self.drho_dP(y,T)
    rho = self.rho(y,T)
    return mu/kb/(rho+c*math.pow(rho,2.))*(1.-rho*(1.+2.*c*rho)/(rho+c*math.pow(rho,2.))*drhodP)*z







