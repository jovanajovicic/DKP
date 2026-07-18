# -*- coding: utf-8 -*-
"""
Relativisticka jednacina kretanja naelektrisanja u dipolnom magnetnom polju Zemlje.

Integracija se vrsi Boris-C metodom

Struktura programa (moduli PVREDNOSTI, DEFPOLJA, IZLAZ) je ista kao u
primeru sa vezbi (REL_PUTANJE_RK4.py), samo je algoritam integracije
zamenjen.
"""

import numpy as np
import matplotlib.pyplot as plt
import PVREDNOSTI, IZLAZ, DEFPOLJA  #U istom direktorijumu kao i glavni program
from mpl_toolkits.mplot3d import Axes3D

N = 1  # Broj cestica u razmatranoj plazmi - npr. jedna cestica

T_sim = 6  #trajanje simulacije
dt = 0.0001  #vremenski korak
T_smp = 0.0001  #interval uzorkovanja za grafik

c = 299792458.0  #[m/s]

#inicijalizacije
t = 0.0
q, m, x, y, z, vx, vy, vz = PVREDNOSTI.puslovi(N, np)
#Pocetnu (obicnu) brzinu pretvaramo u prostorni deo 4-impulsa u = gamma*v,
#jer se u Boris algoritmu radi sa u, a ne direktno sa v (jedn. 1 i 2)
v2 = vx ** 2.0 + vy ** 2.0 + vz ** 2.0
gama0 = 1.0 / np.sqrt(1.0 - (v2 / (c ** 2.0)))
ux = gama0 * vx
uy = gama0 * vy
uz = gama0 * vz

fig, sp1, sp2 = IZLAZ.init_plots(plt)

xx = np.zeros(60001)
yy = np.zeros(60001)
zz = np.zeros(60001)

ii = 0

xx[ii] = x[0]
yy[ii] = y[0]
zz[ii] = z[0]

l = 0
#Boris algoritam radi u "leapfrog" rasporedu: pozicija x se racuna na
#polu-celobrojnim vremenskim koracima (n-1/2, n+1/2, ...), dok se u
#(impuls) racuna na celobrojnim koracima (n, n+1, ...) - j-na 1 u radu
#zato prvo pomeramo pocietnu poziciju za pola koraka unapred.

gama_n = np.sqrt(1.0 + (ux ** 2.0 + uy ** 2.0 + uz ** 2.0) / (c ** 2.0))
x = x + (0.5 * dt * ux / gama_n)
y = y + (0.5 * dt * uy / gama_n)
z = z + (0.5 * dt * uz / gama_n)

# Racun putanje
#Boris-C algoritam (Zenitani & Umeda, 2018, Phys. Plasmas 25, 112110)
#j-ne 3,6,11,12,5,1 

while t < (T_sim - dt):

    print(t)

    #polje na trenutnoj poziciji cestice
    BPx, BPy, BPz, EPx, EPy, EPz = DEFPOLJA.magdipol(x, y, z, t, np)

    
    #Prva polovina ubrzanja usled E polja - j-na 3
    epsx = (q / (2.0 * m)) * EPx
    epsy = (q / (2.0 * m)) * EPy
    epsz = (q / (2.0 * m)) * EPz

    umx = ux + (epsx * dt)
    umy = uy + (epsy * dt)
    umz = uz + (epsz * dt)

    #Boris-C rotacija usled B polja
    #j-ne 6,11,12   

    Bmag2 = (BPx ** 2.0) + (BPy ** 2.0) + (BPz ** 2.0)
    Bmag2 = np.maximum(Bmag2, 1.0e-30)  # da se izbegne deljenje nulom
    Bmag = np.sqrt(Bmag2)

    gama_m = np.sqrt(1.0 + (umx ** 2.0 + umy ** 2.0 + umz ** 2.0) / (c ** 2.0))

    #tacan ugao rotacije za dati korak j-na 6
    theta = (q * dt * Bmag) / (m * gama_m)

    bx = BPx / Bmag
    by = BPy / Bmag
    bz = BPz / Bmag

    #komponenta u paralelna sa B j-na 11
    u_dot_b = (umx * bx) + (umy * by) + (umz * bz)
    uparx = u_dot_b * bx
    upary = u_dot_b * by
    uparz = u_dot_b * bz

    uperpx = umx - uparx
    uperpy = umy - upary
    uperpz = umz - uparz

    #vektorski proizvod u x b - za j-nu 12
    ucrossx = (umy * bz) - (umz * by)
    ucrossy = (umz * bx) - (umx * bz)
    ucrossz = (umx * by) - (umy * bx)

    #rotacija - j-na 12
    upx = uparx + (uperpx * np.cos(theta)) + (ucrossx * np.sin(theta))
    upy = upary + (uperpy * np.cos(theta)) + (ucrossy * np.sin(theta))
    upz = uparz + (uperpz * np.cos(theta)) + (ucrossz * np.sin(theta))

    #druga polovina ubrzanja usled E polja - j-na 5
    
    ux = upx + (epsx * dt)
    uy = upy + (epsy * dt)
    uz = upz + (epsz * dt)

    t = t + dt
    
    #pomeranje pozicije za ceo korak koriscenjem novog 4-impulsa
    #j-na 1

    gama_np1 = np.sqrt(1.0 + (ux ** 2.0 + uy ** 2.0 + uz ** 2.0) / (c ** 2.0))
    x = x + (dt * ux / gama_np1)
    y = y + (dt * uy / gama_np1)
    z = z + (dt * uz / gama_np1)

    ii = ii + 1
    xx[ii] = x[0]
    yy[ii] = y[0]
    zz[ii] = z[0]

    #graficko predstavljanje rezultata

    l += 1

    if l * dt >= T_smp:

        l = 0.0
        IZLAZ.write_plots(sp1, sp2, plt, x, y, z, t, T_sim)  #prikaz u realnom vremenu


print("x range:", np.min(xx), np.max(xx))
print("y range:", np.min(yy), np.max(yy))
print("z range:", np.min(zz), np.max(zz))
print("any NaN in xx?", np.isnan(xx).any())
print("any NaN in yy?", np.isnan(yy).any())
print("any NaN in zz?", np.isnan(zz).any())

d = 6378137.0
fig = plt.figure()  #graficki prikaz
ax = fig.add_subplot(111, projection="3d")
ax.set_aspect("auto")
ax.grid(False)
plt.xlabel("$x[R_T]$")
plt.ylabel("$y[R_T]$")
ax.set_zlabel("$z[R_T]$")
plt.axis("on")
plt.plot(xx / d, yy / d, zz / d, color="green")
plt.savefig('grafik.png')
plt.show()
