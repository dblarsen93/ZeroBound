'''
Modeldefs and ModelDyn functions for Simple ILA Model
'''
import numpy as np


def Steady(params):
    [delta, beta, gamma, chi, theta, alpha, xi, phiy, phipi, pibar, rho, psi, sigma, omega] = params
    
    # find steady state values
    rbar = beta**(1/gamma) + delta - 1
    kbar = (alpha/rbar)**(1/(1-alpha))
    GDPbar = kbar**alpha
    wbar = (1-alpha)*GDPbar
    cbar = wbar + (rbar-delta)*kbar
    ibar = (1+rbar)*(1+pibar) - 1
    mbar = (cbar**-gamma/chi*(1+pibar/beta - 1))**(1/theta)
    Pbar = 1.
    return kbar, mbar, ibar, Pbar

def Modeldefs(Kp, K, M, P, i, Z, X, params, t):
    '''
    This function takes vectors of endogenous and exogenous state variables
    along with a vector of 'jump' variables and returns explicitly defined
    values for consumption, gdp, wages, real interest rates, and transfers
    
    Inputs are:
        Kp: value of capital in next period
        K: value of capital this period
        Z: value of productivity this period
        X: Money supply shock
        params: list of parameter values
    
    Outputs are:
        GDP: GDP
        w: wage rate
        r: rental rate on capital
        c: consumption
        i: investment
        u: utiity
    '''
    
    # unpack input vectors
    kp = Kp
    k = K
    m = M
    z = Z
    x = X
    
    # unpack params
    [delta, beta, gamma, chi, theta, alpha, xi, phiy, phipi, pibar, rho, psi, sigma, omega] = params
    
    # find steady state values
    rbar = beta**(1/gamma) + delta - 1
    kbar = (alpha/rbar)**(1/(1-alpha))
    GDPbar = kbar**alpha
    # find definition values
    GDP = k**alpha*(np.exp(z)**(1-alpha))
    w = (1-alpha)*GDP
    r = alpha*GDP/k
    #pi = ((1 + r - delta)/(1+i)) - 1
    pi = (1+i)/(1+r) - 1
    Gam = (rbar + pibar) + phiy*np.log(GDP/GDPbar) + phipi*(pi-pibar) + x
    ip = np.sqrt(1/4*Gam**2 + xi**2) + 1/2*Gam
    c = w + (1+r-delta)*k - kp
    mhat = m*(1+pibar)**(-t)
    phat = P*(1+pibar)**(-t)
    #phat = (1+pi)/(1+pibar)
    u = c**(1-gamma)-1/(1-gamma) - chi*(m**(1-theta)-1)/(1-theta)
    return GDP, w, r, c, ip, mhat, phat, u, pi


def Modeldyn(theta0, params):
    '''
    This function takes vectors of endogenous and exogenous state variables
    along with a vector of 'jump' variables and returns values from the
    characterizing Euler equations.
    
    Inputs are:
        theta: a vector containng (Kpp, Kp, K, Mp, M, Pp, P, ip, i, Zp, Z, Xp, X) where:
            Kpp: value of capital in two periods
            Kp: value of capital in next period
            K: value of capital this period
            Mp: value of money next period
            M: value of money this period
            Pp: prices next period
            P: prices this period
            ip: interest rates next period
            i: interest rates this period
            Zp: value of productivity in next period
            Z: value of productivity this period
            Xp: shock to money supply next period
            X: shock to money supply this period
        params: list of parameter values
    
    Output are:
        Euler: a vector of Euler equations written so that they are zero at the
            steady state values of X, Y & Z.  This is a 2x1 numpy array. 
    '''
    
    # unpack theat0
    (Kpp, Kp, K, mp, m, Pp, P, ip, i, Zp, Z, Xp, X) = theta0
    
    # unpack params
    [delta, beta, gamma, chi, theta, alpha, xi, phiy, phipi, pibar, rho, psi, sigma, omega] = params
    
    # find definitions for now and next period
    GDP, w, r, c, ip, mhat, phat, u, pi = Modeldefs(Kp, K, m, P, i, Z, X, params, 101)
    GDPp, wp, rp, cp, ipp, mhatp, phatp, up, pip = Modeldefs(Kpp, Kp, mp, Pp, ip, Zp, Xp, params, 102)
    
    # Evaluate Euler equations
    E1 = (c**(-gamma)*(1 + ip)) / (beta*(cp**(-gamma) + chi*(mp/Pp)**(-theta))*(1 + rp)) - 1
    E2 = (c**(-gamma)) / (beta*cp**(-gamma)*(1 + rp - delta)) - 1
    
    return np.array([E1, E2])