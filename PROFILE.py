import numpy as np
import math
from CONST import PI, G, mu, NA, kb, RE, ME

class PROFILE:
  def __init__(self,EoS,r,m,rho,P,T):
    N = len(m)
    self.N = N

    self.r = r
    self.m = m
    self.rho = rho
    self.P = P
    self.T = T

    self.Pc = P[0]

    self.fname = 'Polytrope'




















