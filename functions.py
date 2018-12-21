import matplotlib.pyplot as plt
import numpy as np

# CONSTANTS
mm2m = 1e-3
P_ref = 101.325     # kPa reference pressure = 1 atm
T_kelvin = 273.15
TOL = 1e-10
M_air = 28.96     #molecular weight of air

def CtoF(T_C):
    return (T_C)*9/5+32


#Standing, Marshall Burton. Volumetric and phase behavior of oil field hydrocarbon systems, 1977
def mu1_uncorrected(G, T_C):
    T_F = CtoF(T_C)
    return (1.709e-5-2.062e-6*G)*T_F + 8.188e-3 -6.15e-3*np.log10(G)

def N2_correction(moleN2,G):
    return moleN2*(8.48e-3*np.log10(G) + 9.59e-3)

def CO2_correction(moleCO2,G):
    return moleCO2*(9.08e-3*np.log10(G) + 6.24e-3)

def H2S_correction(moleH2S,G):
    return moleH2S*(8.49e-3*np.log10(G) + 3.73e-3)

#Dempsey's correlation
def S_correction(Pr, Tr):

    A = np.array([-2.46211820,  2.97054714,  -0.28626405, 0.008054205, 2.80860949,
         -3.49803305, 0.360373020, -0.01044324, -0.793385684, 1.39643306,
         -0.149144925, 0.00441015512, 0.0839387178 , -0.186408848, 0.0203367881,
         -0.000609579263])
    Pr2 = np.power(Pr, 2)
    Pr3 = np.power(Pr, 3)
    Tr2 = np.power(Tr, 2)
    Tr3 = np.power(Tr, 3)
    temp =          A[0]  + A[1]  * Pr + A[2]  * Pr2 + A[3]  * Pr3  \
           + Tr  * (A[4]  + A[5]  * Pr + A[6]  * Pr2 + A[7]  * Pr3) \
           + Tr2 * (A[8]  + A[9]  * Pr + A[10] * Pr2 + A[11] * Pr3) \
           + Tr3 * (A[12] + A[13] * Pr + A[14] * Pr2 + A[15] * Pr3)
    return np.exp(temp)/Tr


def calc_mu(Pr, Tr, G, T_C ,moleN2,moleCO2, moleH2S):
    mu1 = mu1_uncorrected(G, T_C) + N2_correction(moleN2,G) + CO2_correction(moleCO2,G) + H2S_correction(moleH2S,G)
    S = S_correction(Pr, Tr)
    return S*mu1

def mu_test():

    M = np.linspace(15,60, 20)
    G = M/M_air
    temp = np.linspace(-10,100,12)
    plt.figure(1)
    plt.subplot(221)

    for t in temp:
        mu1_u = []
        for g in G:
            mu1_u.append(mu1_uncorrected(g, t))

        plt.plot(M,mu1_u, label= str('%.0f' % t) + ' C$\degree$')

    plt.xlabel('Molecular weight')
    plt.ylabel('Viscosity (cp)')
    plt.title('$\mu_{1 uncorrected}$')
    plt.grid(True)
    plt.legend(fontsize='xx-small')
    plt.xlim(15, 60)
    plt.ylim(0.0060, 0.0140)

    plt.subplot(222)
    N2 = np.linspace(0, 20, 30)
    N2_r = N2/100
    M = np.array([15., 20., 25., 30., 40., 50., 60.])
    G = M / M_air
    for g in G:
        N2_c = []
        for n in N2_r:
            N2_c.append(N2_correction(n,g))

        plt.plot(N2, N2_c, label='M= ' + str('%.0f' % (g*M_air)))

    plt.xlabel('$N_2$ mole fraction')
    plt.ylabel('$N_2$ correction (cp)')
    plt.title('$N_{2 correction}$')
    plt.grid(True)
    plt.legend(fontsize='xx-small')
    plt.xlim(0, 20)
    plt.ylim(0.00, 0.0025)

    plt.subplot(223)
    CO2 = np.linspace(0, 20, 30)
    CO2_r = CO2 / 100
    M = np.array([15., 20., 25., 30., 40., 50., 60.])
    G = M / M_air
    for g in G:
        CO2_c = []
        for n in CO2_r:
            CO2_c.append(CO2_correction(n, g))

        plt.plot(CO2, CO2_c, label='M= ' + str('%.0f' % (g * M_air)))

    plt.xlabel('$CO_2$ mole fraction')
    plt.ylabel('$CO_2$ correction (cp)')
    plt.title('$CO_{2 correction}$')
    plt.grid(True)
    plt.legend(fontsize='xx-small')
    plt.xlim(0, 20)
    plt.ylim(0.00, 0.0025)

    plt.subplot(224)
    H2S = np.linspace(0, 20, 30)
    H2S_r = H2S / 100
    M = np.array([15., 20., 25., 30., 40., 50., 60.])
    G = M / M_air
    for g in G:
        H2S_c = []
        for n in H2S_r:
            H2S_c.append(H2S_correction(n, g))

        plt.plot(H2S, H2S_c, label='M= ' + str('%.0f' % (g * M_air)))

    plt.xlabel('$H_2S$ mole fraction')
    plt.ylabel('$H_2S$ correction (cp)')
    plt.title('$H_2S_{ correction}$')
    plt.grid(True)
    plt.legend(fontsize='xx-small')
    plt.xlim(0, 20)
    plt.ylim(0.00, 0.0025)

    plt.show()
    plt.gcf().clear()

    original = True

    plt.figure(2)
    if original:
        Pr = np.linspace(0, 4, 30)
        Tr = np.array([1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2., 2.2, 2.5, 3.])
        for t in Tr:
            S = []
            for p in Pr:
                S.append(S_correction(p, t))

            plt.plot(Pr, S, label='$T_R$ = ' + str('%.1f' % (t)))

        plt.xlabel('Reduced pressure $P_R$')
        plt.xlim(0, 4)
        plt.ylim(1, 3)
    else:
        Tr = np.linspace(1, 3, 20)
        Pr = np.array([1, 2, 3, 4, 5, 6, 8, 10, 15, 20])
        for p in Pr:
            S = []
            for t in Tr:
                S.append(S_correction(p, t))

            plt.plot(Tr, S, label='$P_R$ = ' + str('%.1f' % (p)))

        plt.xlabel('Reduced temperature $T_R$')
        plt.xlim(0.8, 3.2)
        plt.ylim(1, 6)

    plt.ylabel('$S$ correction')
    plt.title('$S$ correction')
    plt.grid(True)
    plt.legend()
    plt.show()
    plt.gcf().clear()


#Check John Lee's book
#Lee, W. John, and Robert A. Wattenbarger. "Gas reservoir engineering." (1996): 214-226.
# Z factor (Hall-Yarborough)
def Z_fac_HY(Pr, Tr):
    t = 1/Tr
    y = 0.001

    while True:
        F = -0.06125*Pr*t*np.exp(-1.2*(1-t)*(1-t)) + (y+y*y+y*y*y-y*y*y*y)/np.power(1.-y,3.) - (14.76*t-9.76*t*t+4.58*t*t*t)*y*y \
            + (90.7*t - 242.2*t*t + 42.4*t*t*t)*np.power(y,2.18+2.82*t)

        if abs(F) < TOL:
            break

        dFdy = (1 + 4 * y + 4 * y * y - 4 * y * y * y + y * y * y * y) / pow(1 - y, 4) - (29.52 * t -19.52*t*t + 9.16*t*t*t)*y \
               + (2.18 + 2.82*t)*(90.7*t - 242.2*t*t + 42.4*t*t*t)*pow(y, 1.18+2.82*t)
        y = y - F/dFdy
        if y < 0:
            y = 0.

    return  0.06125*Pr*t*np.exp(-1.2*(1-t)*(1-t))/y

# Z factor (Dranchuk)
def Z_fac_Dran(Pr, Tr):
    rhor = 0.0001
    A = [0.3265, -1.0700/Tr, -0.5339/pow(Tr,3), 0.01569/pow(Tr,4), -0.05165/pow(Tr,5),
         0.5475, -0.7361/Tr, 0.1844/pow(Tr,2),
         0.1056, 0.6134, 0.7210]
    B1 = A[0] + A[1] + A[2] + A[3] + A[4]
    B2 = A[5] + A[6] + A[7]
    B3 = A[8] * (A[6] + A[7])
    B4 = A[9]/pow(Tr,3)

    while True:

        EXP = np.exp(-A[10] * pow(rhor, 2))
        F = -0.27*Pr/rhor/Tr + 1 + B1*rhor + B2*pow(rhor,2) \
            - B3*pow(rhor,5) + B4*(pow(rhor,2)+A[10]*pow(rhor,4))*EXP

        dFdy = 0.27*Pr/rhor/rhor/Tr + B1 + 2*B2*rhor - 5*B3*pow(rhor,4) \
               + B4*(2*rhor + 4*A[10]*pow(rhor,3))*EXP \
               + B4*(pow(rhor,2)+A[10]*pow(rhor,4))*EXP*(-2*A[10]*rhor)

        if abs(F / dFdy) < TOL:
            rhor = rhor - F / dFdy
            break

        rhor = rhor - F/dFdy

    Z = 0.27*Pr/rhor/Tr
    return Z


def Zfac_test():

    Pr_high =  np.linspace(0.01,5, 30)
    Pr_low = np.linspace(0.01, 0.5, 30)
    Tr = np.array([1.2, 1.25, 1.3, 1.35, 1.4,1.45,1.5,1.6,1.7,1.8,1.9,2.0, 2.2, 2.4,2.6,2.8,3])

    for t in Tr:
        Z_high = []
        Z_low = []
        for p in Pr_high:
            Z_high.append(Z_fac_HY(p, t)) if t >= 1.6 else  Z_high.append(Z_fac_Dran(p, t))

        for p in Pr_low:
            Z_low.append(Z_fac_HY(p, t)) if t >= 1.6 else Z_low.append(Z_fac_Dran(p, t))

        plt.figure(1)
        plt.subplot(211)
        plt.plot(Pr_high,Z_high, label= 'Tr = '+ str(t))
        plt.xlabel('Reduced pressure $P_R$')
        plt.ylabel('$Z$ factor')
        plt.grid(True)
        plt.xlim(0, 5)
        plt.ylim(0.5, 1.1)
        plt.legend(fontsize='xx-small')

        plt.subplot(212)
        plt.plot(Pr_low, Z_low, label='Tr = ' + str(t))
        plt.xlabel('Reduced pressure $P_R$')
        plt.ylabel('$Z$ factor')
        plt.grid(True)
        plt.xlim(0, 0.5)
        plt.ylim(0.89, 1.01)
        plt.legend(fontsize='xx-small')

    plt.show()
    plt.close()


# Saturation gas vapor pressure
# Goff–Gratch formula
def PD_goff(T): # T is temp in kelvin
    Tst = T_kelvin + 100
    logPD = -7.90298*(Tst/T-1.) + 5.02808*np.log10(Tst/T)-1.3816e-7*(pow(10,11.344*(1-T/Tst))-1.) \
             + 8.1328e-3*(pow(10,-3.49149*(Tst/T-1.))-1.) + np.log10(1013.25)
    PD = pow(10,logPD)/10

    return PD


mu_test()
Zfac_test()