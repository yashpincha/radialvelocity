import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit,fsolve

font = 11 

# the keplerian solver

def keplerian_fit(t,K,P,e,omega,tau,kappa):
  e_anomaly = solve_kepler((t-tau)*2*np.pi/P,e) # solve_kepler takes returns eanom
  nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(0.5*e_anomaly), np.sqrt(1-e)*np.cos(0.5*e_anomaly)) # np.arctan2[x, y] returns arctan[x/y]
  return K*(np.cos(nu+omega)+e*np.cos(omega))+kappa # the cosine fit to the keplerian problem

def solve_kepler(M,e):
  eanom = np.zeros(M.shape) 
  for i, me in enumerate(M):
    # do iterative root solve with e=0 giving E=M as guess
    temp, =fsolve(lambda E: E-e*np.sin(E)-me,me)
    eanom[i] = temp
  return eanom

'''
The code is based on the following three equations
M = E - esinE 
nu = 2arctan(sqrt(1+e/1-e)tan(E/2))
v = K(cos*nu + omega) + ecos(omega)) -- K = na(sini)/sqrt(1-e2), the semi-amplitude
'''

# mass of star with its error, obtain from the NASA exoplanet archive
ms       = 1.44 
ms_error = 0.050

data = pd.read_csv('')
jd = np.array(data.t) 
rv = np.array(data.vel)
er = np.array(data.errvel)

# make plot of the RV data
fig1=plt.figure(figsize=(4,3))
axes = fig1.add_axes([0.14,0.14,0.85,0.78])
axes.plot(jd,rv,'.k')
axes.tick_params(which='both',direction='in',top=True,right=True)
axes.tick_params(which='major',length=5)
axes.tick_params(which='minor',length=3)
axes.tick_params(axis='x',pad=5)
plt.title('RV data for 51 Peg b',fontsize=font)
plt.xlabel('HJD (days)',fontsize=font)
plt.ylabel('radial velocity (m/s)',fontsize=font)
plt.show()
frequency, power = LombScargle(jd,rv,er).autopower(maximum_frequency=1,minimum_frequency=0.01) # the period should be within this frame, otherwise adjust it accordingly

fig2=plt.figure(figsize=(4,3))
axes = fig2.add_axes([0.14,0.14,0.85,0.78])
axes.semilogx(1./frequency,power,'k')
axes.tick_params(which='both',direction='in',top=True,right=True)
axes.tick_params(which='major',length=5)
axes.tick_params(which='minor',length=3)
axes.tick_params(axis='x',pad=5)
plt.title('Lomb-Scargle periodogram for 51 Peg b',fontsize=font)
plt.xlabel('period (days)',fontsize=font)
plt.ylabel('power',fontsize=font)
plt.show()

period = 1/frequency[np.argmax(power)]
print("The fit period is %10.5f days" % period)

tfit = np.linspace(0,period,1000)
rvfit = LombScargle(jd,rv,er).model(tfit,1/period)
semi_amplitude = 0.5*(np.max(rvfit)-np.min(rvfit)) # averaging out to find
print("The fit semi-amplitude is %10.5f m/s" % semi_amplitude)
phase = (jd % period)

voffset = np.mean(rvfit)
print("The velocity offset is %10.5f m/s" % voffset)


fig3=plt.figure(figsize=(4,3))
axes = fig3.add_axes([0.14,0.14,0.85,0.78])
axes.errorbar(phase,rv,er,fmt='.k')
axes.plot(tfit,rvfit,'-k')
axes.tick_params(which='both',direction='in',top=True,right=True)
axes.tick_params(which='major',length=5)
axes.tick_params(which='minor',length=3)
axes.tick_params(axis='x',pad=5)
plt.title('phased RV data for 51 Peg b',fontsize=font)
plt.xlabel('time (days)',fontsize=font)
plt.ylabel('radial velocity (m/s)',fontsize=font)
plt.yticks([-75,-50,-25,0,25,50,75])
plt.show()

K = semi_amplitude
P = period
e = 0.
omega = 0.
tau = jd[np.argmax(rv)]
kappa = voffset
guess = (K,P,e,omega,tau,kappa)

rvfit = keplerian_fit(jd,K,P,e,omega,tau,kappa)
chisq = np.sum(((rv-rvfit)/er)**2) #minimize leastqs
print("Chi-squared of initial guess is %10.5f" % chisq)

popt, pcov = curve_fit(keplerian_fit,jd,rv,
                      sigma=er,absolute_sigma=True,
                      p0=guess)

(K,P,e,omega,tau,kappa) = popt
print(popt)
print(f'semi amplitude:{K}')
print(f'period:{P}')
print(f'eccentricity:{e}')
print(f'periastron:{omega}')
print('epoch of periastron %10.5f' % tau)
print(f'vel_offset:{kappa}')

rvfit = keplerian_fit(jd,K,P,e,omega,tau,kappa)
chisq = np.sum(((rv-rvfit)/er)**2)
print("Chi-squared of least-squares fit is %10.5f" % chisq)

tfit = np.linspace(0,P,1000)
rvfit = keplerian_fit(tfit,K,P,e,omega,tau,kappa)
phase = (jd % P)

fig4=plt.figure(figsize=(4,3))
axes = fig4.add_axes([0.14,0.14,0.85,0.78])
axes.errorbar(phase,rv,er,fmt='.k')
axes.plot(tfit,rvfit,'-k')
axes.tick_params(which='both',direction='in',top=True,right=True)
axes.tick_params(which='major',length=5)
axes.tick_params(which='minor',length=3)
axes.tick_params(axis='x',pad=5)
plt.title('phased RV data for 51 Peg b',fontsize=font)
plt.xlabel('time (days)',fontsize=font)
plt.ylabel('radial velocity (m/s)',fontsize=font)
plt.yticks([-75,-50,-25,0,25,50,75])
plt.savefig('51Peg-RVphased-best.pdf')

if e<0:
  omega -= np.pi
  e *= -1

if K<0:
  K *= -1
  omega += np.pi


P_yr   = P/365.2422                    # period in years
a_au   = (ms*P_yr**2)**(1./3)          # semi-major axis in au
K_auyr = K*2.1096256684e-4             # K in au/yr

mp      = (2*np.pi)**(-1)*K_auyr*np.sqrt(1-e**2)*(ms**2*P_yr)**(1/3)
mp_mjup = mp*1047.59421               # convert to Jupiter mass

w_deg = omega*180/np.pi

print(f'semi major: {a_au}, mass in mjup: {mp_mjup}')

K_error   = np.sqrt(pcov[0,0])
P_error   = np.sqrt(pcov[1,1])
e_error   = np.sqrt(pcov[2,2])
w_error   = np.sqrt(pcov[3,3])
tau_error = np.sqrt(pcov[4,4])
a_error   = a_au * np.sqrt( (2*P_error/(3*P))**2 + (ms_error/(3*ms))**2 )
mp_error  = mp_mjup * np.sqrt( (2*ms_error/(3*ms))**2 + (K_error/K)**2
                             + (P_error/(3*P))**2 + (e*e_error/np.sqrt(1-e*e))**2 )

print(K_error,P_error,e_error,w_error,tau_error,a_error,mp_error)
#rvfit = np.array([keplerian_fit(t,*popt) for t in tfit])
#resid = rv - np.array
