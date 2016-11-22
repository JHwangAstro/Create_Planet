import numpy as np
import math
import PLOTS
import matplotlib.pyplot as plt
from CONST import PI, G, mu, NA, kb, RE, ME
from scipy.special import hyp2f1


def TWO_COMPONENT_DIFF(h,fout,fdir,M_E,R_E,P_E,d_E,T_E,mu_E,u_E,gam_c,PROFILE,mi):
  xx = [0.9207007352,
  0.0783643418,
  0.0002478098,
  0.000062247,
  0.0004509397,
  7.83643417535699E-005,
  3.66537564619843E-005,
  2.97932917843479E-005,
  2.91151136664617E-005
  ]

  M_C = PROFILE.m[PROFILE.r>h*100.]
  R_C = PROFILE.r[PROFILE.r>h*100.]
  d_C = PROFILE.rho[PROFILE.r>h*100.]
  P_C = PROFILE.P[PROFILE.r>h*100.]
  T_C = T_E[len(M_E)-1]
  mu_C= 2.0
  print 'check: ', M_C/max(M_C)
  print 'check: ', mi, len(P_C)
#  u_C = P_C+1.
  #for i in range(len(u_C)):
  #  print i
  #  u_C[i] = P_C[i]/(d_C[i]-4.26)/(1./0.549-1.)
  #  if (M_C[i]/max(M_C)) >= mi:
  #    u_C[i] = P_C[i]/(d_C[i]-8.3)/(1./0.528-1.)
  #u_C = np.array([P_C[i]/(d_C[i]-4.26)/(1./0.549-1.) if (M_C[i]/max(M_C))>=mi else P_C[i]/(d_C[i]-8.3)/(1./0.528-1.) for i in range(len(P_C))])
  u_C = P_C/(d_C-4.26)/(1./0.549-1.)
  gam1 = 1./0.528
  gam2 = 1./0.549
  u_C[M_C/max(M_C)<=mi] = P_C[M_C/max(M_C)<=mi]/d_C[M_C/max(M_C)<=mi]/(gam1-1.)*(1.-8.30/d_C[M_C/max(M_C)<=mi])**(-gam1)*hyp2f1(1.-gam1,-gam1,2.-gam1,8.30/d_C[M_C/max(M_C)<=mi])#(-0.71640366*d_C[M_C/max(M_C)<=mi]**2+15.98296227*d_C[M_C/max(M_C)<=mi]+0.99663563)
  u_C[M_C/max(M_C)>=mi] = P_C[M_C/max(M_C)>=mi]/d_C[M_C/max(M_C)>=mi]/(gam2-1.)*(1.-4.26/d_C[M_C/max(M_C)>=mi])**(-gam2)*hyp2f1(1.-gam2,-gam2,2.-gam2,4.26/d_C[M_C/max(M_C)>=mi])#(-0.71640366*d_C[M_C/max(M_C)<=mi]**2+15.98296227*d_C[M_C/max(M_C)<=mi]+0.99663563)
  p2 = u_C + 1.
  p2[M_C/max(M_C)<=mi] = u_C[M_C/max(M_C)<=mi]*(gam1-1.)*d_C[M_C/max(M_C)<=mi]*(1.-8.30/d_C[M_C/max(M_C)<=mi])**gam1/(-0.71640366*(8.3/d_C[M_C/max(M_C)<=mi])**2.+15.98296227*(8.3/d_C[M_C/max(M_C)<=mi])+0.99663563)
  p2[M_C/max(M_C)>=mi] = u_C[M_C/max(M_C)>=mi]*(gam2-1.)*d_C[M_C/max(M_C)>=mi]*(1.-4.26/d_C[M_C/max(M_C)>=mi])**gam2/(-0.5654*(4.26/d_C[M_C/max(M_C)>=mi])**2.+8.40976*(4.26/d_C[M_C/max(M_C)>=mi])+0.99535)
  N = len(PROFILE.m)

  fname = fout#fdir + PROFILE.fname
  #f = open(fname + '.log', 'w')
  f = open(fname, 'w')

  print fname

  #f.write('Radius\tDensity\tPressure\tEnc Mass\tTemperature')
  skip = len(M_C)/750
  f.write(str(M_C[len(M_C)-1])+'\t'+str(R_C[len(M_C)-1])+'\t'+str(P_C[len(M_C)-1])+'\t'+str(d_C[len(M_C)-1])+'\t'+str(T_C)+'\t'+str(mu_C)+'\t'+str(u_C[len(u_C)-1]))
  for j in range(skip,len(M_C),skip):
    k = len(M_C)-j-1
    f.write('\n'+str(M_C[k])+'\t'+str(R_C[k])+'\t'+str(P_C[k])+'\t'+str(d_C[k])+'\t'+str(T_C)+'\t'+str(mu_C)+'\t'+str(u_C[k]))
  for j in range(0,len(M_E)):
    k = len(M_E)-j-1
    f.write('\n'+str(M_E[k])+'\t'+str(R_E[k])+'\t'+str(P_E[k])+'\t'+str(d_E[k])+'\t'+str(T_E[k])+'\t'+str(mu_E[k])+'\t'+str(u_E[k]))

  f.close()


  plt.plot(R_C, u_C)
  plt.show()
  plt.clf()


  M_C = PROFILE.m
  R_C = PROFILE.r
  d_C = PROFILE.rho
  P_C = PROFILE.P

  #2 panel figure for paper
  if 1:
    plt.plot(R_C/max(R_E),d_C,color='b')
    plt.plot(R_E/max(R_E),d_E,color='g')
    #plt.xscale('log')
    plt.yscale('log')
    #plt.yticks([])
    plt.ylabel(r'$\rho\ [\mathrm{g}\ \mathrm{cm}^{-3}]$', fontsize=20)
    plt.xlabel(r'$r/R$', fontsize=20)
    plt.savefig('/home/jah752/Kepler11dRun68_Profile_rho.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.show()
    plt.clf()

    plt.plot(R_C/max(R_E),M_C/max(M_E),color='b')
    plt.plot(R_E/max(R_E),M_E/max(M_E),color='g')
    plt.axis([0.,1.,0.,1.1])
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.yticks([])
    plt.ylabel(r'$m/M$', fontsize=20)
    plt.xlabel(r'$r/R$', fontsize=20)

    plt.savefig('/home/jah752/Kepler11dRun68_Profile_m.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.show()
    plt.clf()

  plt.figure(1)

  plt.subplot(221)
  plt.plot(R_C,P_C,color='b')
  plt.plot(R_C,p2,color='r')
  plt.plot(R_E,P_E,color='g')
  #plt.axis([min(R)*0.8,max(R),0.,max(P_Degen)*1.1])
  #plt.xscale('log')
  plt.yscale('log')
  plt.yticks([])
  plt.xlabel(r'$R/R_\mathrm{E}$')
  plt.ylabel(r'$P$')

  plt.subplot(222)
  plt.plot(R_C/RE,M_C,color='b')
  plt.plot(R_E/RE,M_E,color='g')
  #plt.xscale('log')
  #plt.yscale('log')
  plt.yticks([])
  plt.ylabel(r'$M$', fontsize=20)
  plt.xlabel(r'$R$', fontsize=20)

  plt.subplot(223)
  plt.plot(R_C/RE,d_C,color='b')
  plt.plot(R_E/RE,d_E,color='g')
  #plt.xscale('log')
  #plt.yscale('log')
  #plt.yticks([])
  plt.ylabel(r'$\rho$', fontsize=20)
  plt.xlabel(r'$R$', fontsize=20)

  if 0:
    plt.subplot(224)
    plt.plot(P_C,d_C,color='b')
    plt.plot(P_E,d_E,color='g')
    plt.xscale('log')
    plt.yscale('log')
    plt.yticks([])
    plt.ylabel(r'$\rho$', fontsize=20)
    plt.xlabel(r'$P$', fontsize=20)

    plt.show()
    plt.clf()

  plt.subplot(224)
  d_C2 = [(d_C[i+1]-d_C[i])/d_C[i] for i in range(0,len(d_C)-1)]
  d_E2 = [(d_E[i+1]-d_E[i])/d_E[i] for i in range(0,len(d_E)-1)]
  #print R_E[0], R_E[len(R_E)-1],R_C[len(R_C)-1],R_C[0]
  print (d_C[0]-d_E[len(R_E)-1])/d_E[len(d_E)-1]
  plt.plot(R_C[1:],d_C2,color='b')
  plt.plot(R_E[1:],d_E2,color='g')
  #plt.xscale('log')
  #plt.yscale('log')
  #plt.yticks([])
  plt.xlabel(r'$r$', fontsize=20)
  plt.ylabel(r'$d\rho$', fontsize=20)

  plt.show()
  plt.clf()

  if 1:
    plt.plot(R_C/RE,M_C)
    plt.plot(R_E/RE,M_E)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$M$', fontsize=20)
    plt.xlabel(r'$R/R_\mathrm{E}$', fontsize=20)
    plt.savefig(fdir + PROFILE.fname + '_log_m_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    plt.plot(R_C/max(R_E),M_C/max(M_E))
    plt.plot(R_E/max(R_E),M_E/max(M_E))
    plt.ylabel(r'$m/M$', fontsize=20)
    plt.xlabel(r'$r/R$', fontsize=20)
    plt.savefig(fdir + PROFILE.fname + '_m_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.savefig(fdir + PROFILE.fname + '_m_r.ps', facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()

    plt.plot(R_C/RE,d_C)
    plt.plot(R_E/RE,d_E)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$\rho$', fontsize=20)
    plt.xlabel(r'$R$', fontsize=20)
    plt.savefig(fdir + PROFILE.fname + '_log_density_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    plt.plot(R_C/RE,d_C)
    plt.plot(R_E/RE,d_E)
    plt.yscale('log')
    plt.ylabel(r'$\rho\ \mathrm{[g\ cm}^{-3}]$', fontsize=20)
    plt.xlabel(r'$R/R_\mathrm{E}$', fontsize=20)
    plt.savefig(fdir + PROFILE.fname + 'density_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    plt.plot(d_C,P_C)
    #plt.axis([10000.,max(R_E),min(d_E),max(PROFILE.rho)])
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$\rho$', fontsize=20)
    plt.xlabel(r'$P$', fontsize=20)
    plt.savefig(fdir + PROFILE.fname + 'P_density.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    #plt.plot(PROFILE.r,PROFILE.T)
    #plt.ylabel(r'$T$', fontsize=20)
    #plt.xlabel(r'$R/R_\mathrm{E}$', fontsize=20)
    #plt.savefig(fdir + PROFILE.fname + 'T_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    #plt.clf()

def TWO_COMPONENT(h,fout,M_E,R_E,P_E,d_E,T_E,mu_E,u_E,gam_c,PROFILE):
  xx = [0.9207007352,
  0.0783643418,
  0.0002478098,
  0.000062247,
  0.0004509397,
  7.83643417535699E-005,
  3.66537564619843E-005,
  2.97932917843479E-005,
  2.91151136664617E-005
  ]

  M_C = PROFILE.m[PROFILE.r>h*100.]
  R_C = PROFILE.r[PROFILE.r>h*100.]
  d_C = PROFILE.rho[PROFILE.r>h*100.]
  P_C = PROFILE.P[PROFILE.r>h*100.]
  T_C = T_E[len(M_E)-1]
  mu_C= 2.0
  u_C = P_C/d_C/(gam_c-1.)

  N = len(PROFILE.m)

  fname = fout#fdir + PROFILE.fname
  #f = open(fname + '.log', 'w')
  f = open(fname, 'w')

  print fname

  #f.write('Radius\tDensity\tPressure\tEnc Mass\tTemperature')
  skip = len(M_C)/750
  f.write(str(M_C[len(M_C)-1])+'\t'+str(R_C[len(M_C)-1])+'\t'+str(P_C[len(M_C)-1])+'\t'+str(d_C[len(M_C)-1])+'\t'+str(T_C)+'\t'+str(mu_C)+'\t'+str(u_C[len(u_C)-1]))
  for j in range(skip,len(M_C),skip):
    k = len(M_C)-j-1
    f.write('\n'+str(M_C[k])+'\t'+str(R_C[k])+'\t'+str(P_C[k])+'\t'+str(d_C[k])+'\t'+str(T_C)+'\t'+str(mu_C)+'\t'+str(u_C[k]))
  for j in range(0,len(M_E)):
    k = len(M_E)-j-1
    f.write('\n'+str(M_E[k])+'\t'+str(R_E[k])+'\t'+str(P_E[k])+'\t'+str(d_E[k])+'\t'+str(T_E[k])+'\t'+str(mu_E[k])+'\t'+str(u_E[k]))

  f.close()



  M_C = PROFILE.m
  R_C = PROFILE.r
  d_C = PROFILE.rho
  P_C = PROFILE.P

  plt.figure(1)

  plt.subplot(221)
  plt.plot(R_C,P_C,color='b')
  plt.plot(R_E,P_E,color='g')
  #plt.axis([min(R)*0.8,max(R),0.,max(P_Degen)*1.1])
  #plt.xscale('log')
  plt.yscale('log')
  plt.yticks([])
  plt.xlabel(r'$R/R_\mathrm{E}$')
  plt.ylabel(r'$P$')

  plt.subplot(222)
  plt.plot(R_C,M_C,color='b')
  plt.plot(R_E,M_E,color='g')
  #plt.xscale('log')
  #plt.yscale('log')
  plt.yticks([])
  plt.ylabel(r'$M$', fontsize=20)
  plt.xlabel(r'$R$', fontsize=20)

  plt.subplot(223)
  plt.plot(R_C,d_C,color='b')
  plt.plot(R_E,d_E,color='g')
  #plt.xscale('log')
  #plt.yscale('log')
  #plt.yticks([])
  plt.ylabel(r'$\rho$', fontsize=20)
  plt.xlabel(r'$R$', fontsize=20)

  if 0:
    plt.subplot(224)
    plt.plot(P_C,d_C,color='b')
    plt.plot(P_E,d_E,color='g')
    plt.xscale('log')
    plt.yscale('log')
    plt.yticks([])
    plt.ylabel(r'$\rho$', fontsize=20)
    plt.xlabel(r'$P$', fontsize=20)

    plt.show()
    plt.clf()

  plt.subplot(224)
  d_C2 = [(d_C[i+1]-d_C[i])/d_C[i] for i in range(0,len(d_C)-1)]
  d_E2 = [(d_E[i+1]-d_E[i])/d_E[i] for i in range(0,len(d_E)-1)]
  #print R_E[0], R_E[len(R_E)-1],R_C[len(R_C)-1],R_C[0]
  print (d_C[0]-d_E[len(R_E)-1])/d_E[len(d_E)-1]
  plt.plot(R_C[1:],d_C2,color='b')
  plt.plot(R_E[1:],d_E2,color='g')
  #plt.xscale('log')
  #plt.yscale('log')
  #plt.yticks([])
  plt.xlabel(r'$r$', fontsize=20)
  plt.ylabel(r'$d\rho$', fontsize=20)

  plt.show()
  plt.clf()

  if 0:
    plt.plot(R_C,M_C)
    plt.plot(R_E,M_E)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$M$', fontsize=20)
    plt.xlabel(r'$R$', fontsize=20)
    plt.savefig(fdir + PROFILE.fname + '_log_m_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    plt.plot(R_C/max(R_E),M_C/max(M_E))
    plt.plot(R_E/max(R_E),M_E/max(M_E))
    plt.ylabel(r'$m/M$', fontsize=20)
    plt.xlabel(r'$r/R$', fontsize=20)
    plt.savefig(fdir + PROFILE.fname + '_m_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.savefig(fdir + PROFILE.fname + '_m_r.ps', facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()

    plt.plot(R_C,d_C)
    plt.plot(R_E,d_E)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$\rho$', fontsize=20)
    plt.xlabel(r'$R$', fontsize=20)
    plt.savefig(fdir + PROFILE.fname + '_log_density_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    plt.plot(R_C,d_C)
    plt.plot(R_E,d_E)
    plt.ylabel(r'$\rho$', fontsize=20)
    plt.xlabel(r'$R$', fontsize=20)
    plt.savefig(fdir + PROFILE.fname + 'density_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    plt.plot(d_C,P_C)
    #plt.axis([10000.,max(R_E),min(d_E),max(PROFILE.rho)])
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$\rho$', fontsize=20)
    plt.xlabel(r'$P$', fontsize=20)
    plt.savefig(fdir + PROFILE.fname + 'P_density.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    #plt.plot(PROFILE.r,PROFILE.T)
    #plt.ylabel(r'$T$', fontsize=20)
    #plt.xlabel(r'$R/R_\mathrm{E}$', fontsize=20)
    #plt.savefig(fdir + PROFILE.fname + 'T_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
    #plt.clf()




def GENLOGS(fdir,PROFILE):

  N = len(PROFILE.m)

  fname = fdir + PROFILE.fname
  f = open(fname + '.log', 'w')

  f.write('Radius\tDensity\tPressure\tEnc Mass\tTemperature')
  for i in range(0,N):
    f.write('\n' + str(PROFILE.r[i])+'\t'+str(PROFILE.rho[i])+'\t'+str(PROFILE.P[i])+'\t'+str(PROFILE.m[i])+'\t'+str(PROFILE.T[i]))
  f.close()

  #index = np.zeros(n_c-1)
  #index[:] = 1.
  #for i in range(0,n_c-1):
  #  if Ps[i] > 0.:
  #    index[i] = min(np.compress(np.where(psr2 < Ps[i],1,0),np.arange(n+1)))




def PLOT_PROFILE(fdir,PROFILE):
  plt.plot(PROFILE.r,PROFILE.P)
  #plt.xscale('log')
  #plt.yscale('log')
  plt.ylabel('$P\ \mathrm{[GPa]}$', fontsize=20)
  plt.xlabel(r'$R/R_\mathrm{E}$', fontsize=20)
  #plt.savefig(fdir + 'P_M.ps',  facecolor='w', edgecolor='w', format='ps', dpi=200)
  plt.savefig(fdir + PROFILE.fname + 'P_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()

  plt.plot(PROFILE.r,PROFILE.m)
  #plt.xscale('log')
  #plt.yscale('log')
  plt.ylabel(r'$M/M_\mathrm{E}$', fontsize=20)
  plt.xlabel(r'$R/R_\mathrm{E}$', fontsize=20)
  #plt.savefig(fdir + 'P_M.ps',  facecolor='w', edgecolor='w', format='ps', dpi=200)
  plt.savefig(fdir + PROFILE.fname + 'm_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()

  plt.plot(PROFILE.r,PROFILE.rho)
  #plt.xscale('log')
  #plt.yscale('log')
  plt.ylabel(r'$\rho\ [M_\mathrm{E}/R_\mathrm{E}^3]$', fontsize=20)
  plt.xlabel(r'$R/R_\mathrm{E}$', fontsize=20)
  #plt.savefig(fdir + 'P_M.ps',  facecolor='w', edgecolor='w', format='ps', dpi=200)
  plt.savefig(fdir + PROFILE.fname + 'density_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()

  plt.plot(PROFILE.r,PROFILE.T)
  #plt.xscale('log')
  #plt.yscale('log')
  plt.ylabel(r'$T$', fontsize=20)
  plt.xlabel(r'$R/R_\mathrm{E}$', fontsize=20)
  #plt.savefig(fdir + 'P_M.ps',  facecolor='w', edgecolor='w', format='ps', dpi=200)
  plt.savefig(fdir + PROFILE.fname + 'T_r.png', facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()




def PLOT_SUMMARY(fdir,PROFILES):

  #fdir = fdir + 'Summary/'
  #f = open(fdir + 'Summary.txt', 'w')
  #f.write('\tPressure\tRadius\tMass')

  #for i in range(0,m):
  #  for j in range(0,n):
  #    for k in range(0,o):
  #      f.write('\n' + str(P[i,j,k]) + '\t' + str(P[i,j,k]*Ps[1]) + '\t'+str(R[i,j])+'\t'+str(M[i,j]))#+'\t'+str(C[j,i]))
  #      for k in range(0,len(c)-1):
  #        f.write('\t'+str(C[i,j,k]))
  #f.close()

  N = len(PROFILES)

  Pc = np.array([PROFILES[i].Pc for i in range(N)])
  print Pc
  M  = np.array([PROFILES[i].m[PROFILES[i].N-1] for i in range(N)])
  R  = np.array([PROFILES[i].r[PROFILES[i].N-1] for i in range(N)])

  plt.plot(Pc, M)
  plt.xscale('log')
  plt.xlabel('$P\ \mathrm{[GPa]}$', fontsize=20)
  plt.ylabel(r'$m_\mathrm{Fe}/M$', fontsize=20)
  #plt.savefig(fdir + 'Summary_P_M.ps',  facecolor='w', edgecolor='w', format='ps', dpi=200)
  plt.savefig(fdir + 'Summary_P_M.png', facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()

  plt.plot(Pc, R)
  plt.xscale('log')
  plt.xlabel('$P\ \mathrm{[GPa]}$', fontsize=20)
  plt.ylabel(r'$R/R_\mathrm{E}$', fontsize=20)
  #plt.savefig(fdir + 'Summary_P_R.ps',  facecolor='w', edgecolor='w', format='ps', dpi=200)
  plt.savefig(fdir + 'Summary_P_R.png', facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()

  plt.plot(M, R)
  #plt.xscale('log')
  plt.xlabel(r'$M/M_\mathrm{E}$', fontsize=20)
  plt.ylabel(r'$R/R_\mathrm{E}$', fontsize=20)
  #plt.savefig(fdir + 'Summary_M_R.ps',  facecolor='w', edgecolor='w', format='ps', dpi=200)
  plt.savefig(fdir + 'Summary_R_M.png', facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()













