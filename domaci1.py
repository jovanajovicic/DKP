# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:10:16 2026
DKP
DOBRO SI KRENULA SAMO NASTAVIIII :)
@author: Jocke
"""

import numpy as np
import matplotlib.pyplot as plt


q_e = 1.602176634e-19 #C
e0 = 8.8541878188e-12 #C**2 kg-1 m-3 s**2
m_e = 9.1093837139e-31 #kg
kB = 1.380649e-23  #J/K

h   = 6.62607015e-34 
c   = 2.99792458e8 

n   = np.logspace(5, 35, 1000)    # gustina 
T_e = np.logspace(1, 10, 1000)    # temperatura


kT   = kB * T_e                              # termalna energija 
w_pe = np.sqrt(n * q_e**2 / (e0 * m_e))     # ugaona ucestanost plazme 
nu_pe = w_pe / (2 * np.pi * 1e6)               # ucestanost plazme 

#debajev radijus
def T_Debajev (n, r_d):
    return (r_d**2 * n * q_e**2) /(e0 * kB)

#granica kvantnih efekata: T = h^2 * n^(2/3) / (kB * me)
def T_q_ef(n):
    return ((h**2/(8 * m_e))*((3 * n)/np.pi)**(2/3) )/kB

#granica N_De = 1
def T_NDE1(n):
    return (3 / (4 * np.pi * n))**(2/3) * n * q_e**2 / (e0 * kB)

#kvantne plazme
def T_q (n):
    return ((h**2/(8 * m_e))*((3 * n)/np.pi)**(2/3) )/kB

#relativ plazme
rel = m_e * c**2 / kB


#pomocne funkcije za sekundarne ose
def m3_to_cm3(n): return np.asarray(n) * 1e-6
def cm3_to_m3(n): return np.asarray(n) * 1e6
def K_to_eV(T):   return np.asarray(T) * kB / q_e
def eV_to_K(E):   return np.asarray(E) * q_e / kB
def n_to_MHz(n_m3):
    n_m3 = np.asarray(n_m3)
    return np.sqrt(n_m3 * q_e**2 / (e0 * m_e)) / (2 * np.pi * 1e6)
def MHz_to_n(f_MHz):
    f_MHz = np.asarray(f_MHz)
    w = f_MHz * 2 * np.pi * 1e6
    return w**2 * e0 * m_e / q_e**2

#plotovanje
fig, ax1 = plt.subplots(figsize=(12, 9))
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(1e5, 1e35)
ax1.set_ylim(1e1, 1e10)
ax1.set_xlabel(r'$n_e\ [\mathrm{m}^{-3}]$', fontsize=13)
ax1.set_ylabel(r'$T_e\ [\mathrm{K}]$', fontsize=13)
ax1.set_title(r'$n_e - T_e$ dijagram za različite tipove plazmi', fontsize=13)
ax1.grid(True, which='major', alpha=0.1)


#kvantne plazme
T_kvantni = T_q(n)
ax1.plot(n, T_kvantni, 'k:', linewidth=2)
ax1.text(1e31, 1e4, 'KVANTNI\nEFEKTI', fontsize=10)

#granice nde=1
T_nde1 = T_NDE1(n)
ax1.plot(n, T_nde1, 'k--', linewidth=2)
ax1.text(5e24, 2e4 , r'$N_{De} = 1$',  fontsize=10)

#relativisticka plazma
ax1.axhline(rel, color='k', lw=2)
ax1.text(1e18, rel * 1.4, 'REL. EFEKTI', fontsize=10)

#dodatna horizontalna linija na 1e7 K
ax1.axhline(1e7, color='black', lw=3)

#postavljanje osa
ax2 = ax1.secondary_xaxis('bottom', functions=(m3_to_cm3, cm3_to_m3))
ax2.spines['bottom'].set_position(('outward', 50))
ax2.set_xlabel(r'$n_e \ [cm^{-3}]$', fontsize=12)

ax3 = ax1.secondary_xaxis('top', functions=(n_to_MHz, MHz_to_n))
ax3.set_xlabel(r'$\omega_{pe} \ [MHz]$', fontsize=12)

ax4 = ax1.secondary_yaxis('right', functions=(K_to_eV, eV_to_K))
ax4.set_ylabel(r'$kT_e \ [eV]$', fontsize=12)

#linije samo iznad n_de=1
ax1.text(2e5, 3e7, r'$r_{De}\ \mathrm{[m]}$', fontsize=10)
for r_d in np.logspace(4, -10, 15):
    T_r = T_Debajev(n, r_d)
    mask = (T_r >= T_nde1) & (T_r >= 1e1) & (T_r <= 1e10)
    if np.any(mask):
        ax1.plot(n[mask], T_r[mask], 'k-.', lw=0.8, alpha=0.6)
        idx = np.where(mask)[0][-1]
        exp = int(round(np.log10(r_d)))
        if exp in [-10, -8, -6, -4, -2, 0, 2, 4]:
            lbl = '$1$' if exp == 0 else f'$10^{{{exp}}}$'
            ax1.text(n[idx]*1.1, T_r[idx], lbl, fontsize=6, color='gray')
            
 #karakteristicne tacke
plazme = [
    (5e5,    1.2e8,  2e4,   2.2e5, "S.V."),
    (2.5e6,  5e8,    9e2,   5e3,   "HDR"),
    (6e10,   2e12,   5e2,   3e3,   "Z.J."),
    (1e10,   8e13,   8e5,   1.6e6, "S.K."),
    (1e15,   8e16,   4.5e5, 1.5e6, "S.K.P."),
    (1e19,   5e20,   2e3,   1.2e4, "S.F."),
    (8e18,   1.5e21, 2e7,   7e7,   "TLAB"),
    (1.5e31, 4e33,   8e6,   3.1e7, "S.J.")
]  
         
for n_min, n_max, t_min, t_max, label in plazme:
    n_c = np.sqrt(n_min * n_max)
    t_c = np.sqrt(t_min * t_max)
    ax1.scatter(n_c, t_c, s=50, color='black', zorder=6)
    ax1.text(n_c, t_c * 2.0, label, ha='center', va='bottom',
             fontsize=9, fontweight='bold', zorder=6)
            
plt.savefig('domaci1.png')
plt.show()