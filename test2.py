import numpy as numpy


import numpy as np



T_i = 150*1.602e-13
r_e = 2.8179e-15
m_e = 9.1094e-31
c = 3e8
m_p = 1.6726e-27
Énergie_moyenne_excitation_eau = 75*1.602e-19
Énergie_moyenne_excitation_os = 91.9*1.602e-19
rho_eau = 997
rho_os = 1850
def densité_électronique(composition_atomique, masse_volumique):
    nbr_électrons_volumique = 0
    avogadro = 6.022e23
    masse_atomique = {1: 0.001007975, 6: 0.0120106, 7: 0.014006855, 8: 0.0159940, 12: 0.0243055, 15: 0.03097396200, 16: 0.0320675, 20: 0.040078} 
    for x in composition_atomique:
        nbr_électrons_volumique += masse_volumique * x[1] * avogadro * x[0] / masse_atomique[x[0]]
    return nbr_électrons_volumique
eau = [(1, 0.111894), (8, 0.888106)]#
densité_électronique_eau = densité_électronique(eau, 997)
os = [(1,0.063984), (6,0.278000), (7,0.027000), (8,0.410016), (12,0.002), (15,0.07), (16,0.002), (20,0.147)]
densité_électronique_os = densité_électronique(os, 1850)
def f_eau(T):
    gamma = T/(m_p*c**2) + 1
    Beta = np.sqrt((gamma**2-1)/gamma**2)
    a = 2*m_e*c**2
    b = 1 + (m_e/m_p)**2
    delta = 2*m_e/m_p
    T_emax = (a*(gamma**2-1))/(b+(delta*gamma))
    S_col_eau = 2 * np.pi * (r_e ** 2)* m_e * (c ** 2) * (densité_électronique_eau) * (1/Beta**2) * (np.log((2*m_e*(c**2)*(Beta**2)*(gamma**2)*T_emax)/((Énergie_moyenne_excitation_eau)**2))-2*(Beta**2))
    return rho_eau/S_col_eau
def f_os(T):
    gamma = T/(m_p*c**2) + 1
    Beta = np.sqrt((gamma**2-1)/gamma**2)
    a = 2*m_e*c**2
    b = 1 + (m_e/m_p)**2
    delta = 2*m_e/m_p
    T_emax = (a*(gamma**2-1))/(b+(delta*gamma))
    S_col_os = 2 * np.pi * (r_e ** 2)* m_e * (c ** 2) * (densité_électronique_os) * (1/Beta**2) * (np.log((2*m_e*(c**2)*(Beta**2)*(gamma**2)*T_emax)/((Énergie_moyenne_excitation_os)**2))-2*(Beta**2))
    return rho_os/S_col_os
def trap_eau(N):
    a = 0
    b = T_i
    h = (b-a)/N
    s = 0.5*f_eau(a) + 0.5*f_eau(b)
    for k in range(1,N):
        s += f_eau(a+k*h)
    return (h*s)
def trap_os(N):
    a = 0
    b = T_i
    h = (b-a)/N
    s = 0.5*f_os(a) + 0.5*f_os(b)
    for k in range(1,N):
        s += f_os(a+k*h)
    return (h*s)
print('''La portée des protons dans l'os est de '''+str(trap_os(1000))+''' [Kg/J*m^2]''')
print('''La portée des protons dans l'os est de '''+str(trap_os(2000))+''' [Kg/J*m^2]''')
print('''La portée des protons dans l'eau est de '''+str(trap_eau(1000))+''' [Kg/J*m^2]''')
print('''La portée des protons dans l'eau est de '''+str(trap_eau(2000))+''' [Kg/J*m^2]''')
print('''La portée des protons dans l'eau est de '''+str(trap_eau(2000))+''' [Kg/J*m^2]''')
print('''La portée des protons dans l'eau est de '''+str(trap_eau(10000))+''' [Kg/J*m^2]''')
tranches_eau = 2
tranches_os = 2
#while True is True:
    #I_ii = trap_os(tranches_os*2)
    #I_i = trap_os(tranches_os)
    #eps = (1/3)*(I_ii-I_i)
    #tranches_os *= 2
    #print(eps)
    #print(tranches_os)
    #if abs(eps) < 2.2e-16:
    #    print('''Le nombre de tranches os est de '''+ str(2*(tranches_os-1)))
    #    break

#while True is True:
    #I_ii = trap_eau(tranches_eau*2)
    #I_i = trap_eau(tranches_eau)
    #eps = (1/3)*(I_ii-I_i)
    #tranches_eau *= 2
    #print(eps)
    #print(tranches_eau)
    #if abs(eps) < 2.2e-16:
    #    print('''Le nombre de tranches eau est de '''+ str(2*(tranches_eau-1)))
    #    print(trap_eau((tranches_eau-1)*2))
    #    print(trap_eau(4355))
    #    break
#while True is True:
    #I_ii = trap_os(tranches_os*2)
    #I_i = trap_os(tranches_os)
    #eps = (1/3)*(I_ii-I_i)
    #tranches_os *= 2
    #print(eps)
    #print(tranches_os)
    #if abs(eps) < 2.2e-16:
    #    print('''Le nombre de tranches os est de '''+ str(2*(tranches_os-1)))
    #    break


print(trap_eau(1), trap_eau(2))

'''Romberg eau'''
def romberg_eau(t):
    R = np.zeros((t, t))
    for x in range(0, t):
        R[x, 0] = trap_eau(2**x)
        for j in range(0, x):
            R[x, j+1] = R[x, j] + (1/((4**(j+1))-1)) * (R[x,j]-R[x-1,j])
    return R
'''Romberg os'''
def romberg_os(t):
    R = np.zeros((t, t))
    for x in range(0, t):
        R[x, 0] = trap_os(2**x)
        for j in range(0, x):
            R[x, j+1] = R[x, j] + (1/((4**(j+1))-1)) * (R[x,j]-R[x-1,j])
    return R
print(romberg_eau(3))
tranches_romberg_eau = 3
# while True is True:
#     end = False
#     n = tranches_romberg_eau
#     tranches_romberg_eau += 1
#     R = romberg_eau(n)
#     #print(R)
#     for i in range(n):
#         for m in range(1,i):
#             eps = (1/(4**m-1))*(R[i,m]-R[i-1,m])
#             print(eps, tranches_romberg_eau)
#             if abs(eps) < 2.2e-16:
#                 print('''Le nombre de tranche avec romberg est '''+str(tranches_romberg_eau))
#                 end = True
#     if end:
#         break

while True is True:
    end = False
    n = tranches_romberg_eau
    tranches_romberg_eau += 1
    R = romberg_eau(n)
    print(R)
    print(R[n-1])
    print(len(R[n-1]))
    i = R[n-1]
    print(i)
    for m in range(len(i)):
        eps = (1/(4**(m+1)-1))*(i[m]-R[n-2][m])
        print(eps, n)
        if abs(eps) < 2.2e-16:
            print('''Le nombre de tranche avec romberg est '''+str(tranches_romberg_eau))
            end = True
    if end:
        break
    