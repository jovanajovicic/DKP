# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:41:46 2026

Kalkulator parametara plazme
Ulaz:  Te [K], ne [cm^-3], B [G]
Izlaz: ne [m^-3], B [T], kT [keV], rDe [m], omega_pe [MHz],
       NDe, omega_ce [MHz], rc_th [m], beta_e

@author: Jocke
"""
import numpy as np

kB   = 1.380649e-23      # J/K
e    = 1.602176634e-19   # C
me   = 9.1093837015e-31  # kg
eps0 = 8.8541878128e-12  # F/m
mu0  = 4 * np.pi * 1e-7  # H/m


def parametri(Te_K, ne, B):
    """
    Te_K   : elektronska temperatura [K]
    ne_cm3 : elektronska gustina [cm^-3]
    B_G    : magnetno polje [G]
    """
    # konverzija u SI
    #ne  = ne * 1e6          # [m^-3]
    #B   = B * 1e-4            # [T]

    # termalna energija
    kT_J   = kB * Te_K          # [J]
    kT_keV = kT_J / (1e3 * e)  # [keV]

    # Debajeva duzina
    rDe = np.sqrt(eps0 * kT_J / (ne * e**2))   # [m]

    # broj cestica u Debajevoj sferi
    NDe = (4/3) * np.pi * ne * rDe**3

    # elektronska plazmena frekvencija
    omega_pe_rad = np.sqrt(ne * e**2 / (eps0 * me))   # [rad/s]
    omega_pe_MHz = omega_pe_rad / (2 * np.pi * 1e6)   # [MHz]

    # elektronska ciklotronska frekvencija
    omega_ce_rad = e * B / me                          # [rad/s]
    omega_ce_MHz = omega_ce_rad / (2 * np.pi * 1e6)   # [MHz]

    # termalna brzina elektrona (v_th = sqrt(kT/me))
    v_th = np.sqrt(kT_J / me)                         # [m/s]

    # termalni Larmorov poluprecnik
    rc_th = v_th / omega_ce_rad    #[m]

    # beta_e = pritisak plazme / magnetski pritisak
    beta_e = (ne * kT_J) / (B**2 / (2 * mu0)) 

    return  ne,B,kT_keV,rDe,omega_pe_MHz,NDe,omega_ce_MHz,rc_th,beta_e
    

ne = np.array([1e20,1e16,1e14,1e11,1e8,1e7,1e3])
Te = np.array([6000,1e6,1e6,1200,5000,1e5,1e6])
B = np.array([0.1,3e-2,1e-4,3e-5,5e-10,1e-8,5e-10])

provera = parametri(Te, ne, B)
print('Kada proverimo funkciju sa vrednostima iz tabele dobijamo (za prvi red):',
      '\n DEBAJEV RADIJUS', provera[3][0], '[m]',
      '\n PLAZMENA FREKVENCIJA' , provera[4][0], '[MHz]',
      '\n CIKLOTRONSKA FREKVENCIJA' , provera[6][0], '[MHz]',
      '\n LARMOROV POLUPRECNIK' , provera[7][0], '[m]',
      '\n BETA' , provera[8][0],
      '\nSto se poklapa sa vrednostima u tabeli.')













