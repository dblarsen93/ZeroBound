'''
Modeldefs and ModelDyn functions for Simple ILA Model
'''
import numpy as np


def Modeldefs(Xp, X, Z, params):
    '''
    This function takes vectors of endogenous and exogenous state variables
    along with a vector of 'jump' variables and returns explicitly defined
    values for consumption, gdp, wages, real interest rates, and transfers
    
    Inputs are:
        Xp: value of capital in next period
        X: value of capital this period
        Y: value of labor this period
        Z: value of productivity this period
        params: list of parameter values
    
    Outputs are:
        GDP: GDP
        w: wage rate
        r: rental rate on capital
        T: transfer payments
        c: consumption
        i: investment
        u: utiity
    '''
    
    # unpack input vectors
    kp = Xp
    k = X
    z = Z
    
    # unpack params
    [delta, beta, gamma, chi, theta, alpha, phiy, phipi, pibar, rho, psi, sigma, omega] = params
    
    # find steady state values
    rbar = beta**(1/gamma) + delta - 1
    kbar = (alpha/rbar)**(1/(1-alpha))
    GDPbar = kbar**alpha
    wbar = (1-alpha)*GDPbar
    cbar = wbar + (rbar-delta)*kbar
    ibar = (1+rbar)*(1+pibar) - 1
    mbar = (cbar**-gamma/chi*(1+pibar/beta - 1))**(1/theta)
    # find definition values
    GDP = k**alpha*(np.exp(z))**(1-alpha)
    w = (1-alpha)*GDP
    r = alpha*GDP/k
    Gam = (rbar + pibar) + phiy*np.log(GDP/GDPbar) + phipi*(pi-pibar) + x
    ip = np.sqrt(1/4*Gam + xi**2) + 1/2*Gam
    #i = np.sqrt(1/4*)
    pi = (1 + r - delta)/(1+i) - 1
    c = w + (1+r-delta)*k - kp

    u = c**(1-gamma)/(1-gamma) - chi*ell**(1+theta)/(1+theta)
    return Y, w, r, T, c, i, u


def Modeldyn(theta0, params):
    '''
    This function takes vectors of endogenous and exogenous state variables
    along with a vector of 'jump' variables and returns values from the
    characterizing Euler equations.
    
    Inputs are:
        theta: a vector containng (Xpp, Xp, X, Yp, Y, Zp, Z) where:
            Xpp: value of capital in two periods
            Xp: value of capital in next period
            X: value of capital this period
            Yp: value of labor in next period
            Y: value of labor this period
            Zp: value of productivity in next period
            Z: value of productivity this period
        params: list of parameter values
    
    Output are:
        Euler: a vector of Euler equations written so that they are zero at the
            steady state values of X, Y & Z.  This is a 2x1 numpy array. 
    '''
    
    # unpack theat0
    (Xpp, Xp, X, Yp, Y, Zp, Z) = theta0
    
    # unpack params
    [alpha, beta, gamma, delta, chi, theta, tau, rho, sigma] = params
    
    # find definitions for now and next period
    ell = Y
    if ell > 0.9999:
        ell = 0.9999
    elif ell < 0.0001:
        ell = 0.0001
    GDP, w, r, T, c, i, u = Modeldefs(Xp, X, Y, Z, params)
    GDPp, wp, rp, Tp, cp, ip, up = Modeldefs(Xpp, Xp, Yp, Zp, params)
    
    # Evaluate Euler equations
    E1 = (c**(-gamma)*(1-tau)*w) / (chi*ell**theta) - 1
    E2 = (c**(-gamma)) / (beta*cp**(-gamma)*(1 + (1-tau)*(rp - delta))) - 1
    
    return np.array([E1, E2])