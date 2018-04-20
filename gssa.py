# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 16:55:35 2017
@author: Daryl Larsen
"""
'''
This program implements the Generalized Stochastic Simulation Algorthim from 
Judd, Maliar and Mailar (2011) "Numerically stable and accurate stochastic 
simulation approaches for solving dynamic economic models", Quantitative
Economics vol. 2, pp. 173-210.
'''
import numpy as np
import matplotlib.pyplot as plt
from LinApp_FindSS import LinApp_FindSS
'''
We test the algorithm with a simple DSGE model with endogenous labor.
'''


def poly1(Xin, XYparams):
    '''
    Includes polynomial terms up to order 'pord' for each element and quadratic 
    cross terms  One observation (row) at a time
    '''
    nX = nx + nz
    Xbasis = np.ones((1, 1))
    # generate polynomial terms for each element
    for i in range(1, 3):
        Xbasis = np.append(Xbasis, Xin**i)
    # generate cross terms
    for i in range (0, nX):
        for j in range(i+1, nX):
            temp = Xin[i]*Xin[j]
            Xbasis = np.append(Xbasis, temp)
    return Xbasis

def XYfunc(Xm, Zn, XYparams, coeffs):
    An = np.exp(Zn)
    XZin = np.append(Xm, An)
    XYbasis = np.append(1., XZin)
    for i in range(1, 3):
        XYbasis = poly1(XZin, XYparams)
    XYout = np.dot(XYbasis, coeffs)
    Xn = XYout[0:nx]
    Y = XYout[nx:nx+ny]
    if Xn < 0:
        Xn = np.array([0])
    if Y > 0.9999:
        Y = np.array([0.9999])
    elif Y < 0:
        Y = np.array([0])
    return Xn, Y
    
def MVOLS(Y, X):
    '''
    OLS regression with observations in rows
    '''
    XX = np.dot(np.transpose(X), X)
    XY = np.dot(np.transpose(X), Y)
    coeffs = np.linalg.solve(XX, XY)
    return coeffs
 
def GSSA(params, kbar, ellbar):
    T = 10000
    regtype = 'poly1' # functional form for X & Y functions 
    fittype = 'MVOLS'   # regression fitting method
    pord = 3  # order of polynomial for fitting function
    ccrit = 1.0E-8  # convergence criteria for XY change
    nx = 1
    ny = 1
    nz = 1
    damp = 0.01  # damping paramter for fixed point algorithm
    
    [delta, beta, gamma, chi, theta, alpha, phiy, phipi, pibar, rho, psi, sigma, omega] = params
    XYparams = (pord, nx, ny, nz)

    Xstart = kbar
    
    #create history of Z's and X's
    Z = np.zeros([T,nz])
    x = np.zeros([T,nz])
    for t in range(1,T):
        Z[t,:] = rho*Z[t-1] + np.random.randn(1)*sigma
        x[t,:] = psi*x[t-1] + np.random.randn(1)*omega
        
    if regtype == 'poly1':
        coeffs = np.array([[  1.14283341e+03,   2.48767872e+01], \
                           [ -9.97029153e+01,  -3.18909650e+00], \
                           [ -3.13420276e+02,   1.38021793e+01], \
                           [  2.15133629e+00,   1.01468495e-01], \
                           [ -1.89909459e+00,   7.37298161e-01], \
                           [  1.63137901e+01,  -7.87431895e-01]])
    
    dist = 1.
    distold = 2.
    count = 0
    damp = .01
    XYold = np.ones((T-1, nx+ny))

    while dist > ccrit:
        count = count + 1
        X = np.zeros((T+1, nx))
        Xin = np.zeros((T, nx+nz))
        A = np.exp(Z)
        x = np.zeros((T,(pord*2)))
        X[0], Y[0] = XYfunc(Xstart, Z[0], XYparams, coeffs)
        for t in range(1,T+1):
            X[t] = XYfunc(X[t-1], Z[t-1], XYparams, coeffs)
            Xin[t-1,:] = np.concatenate((X[t-1], A[t-1]))
            #x[t-1,:] = poly1(Xin[t-1,:], XYparams)
            X1 = X[0:T]
            Y1 = Y[0:T]
            # plot time series
        if count % 10 == 0:
            timeperiods = np.asarray(range(0,T))
            plt.subplot(2,1,1)
            plt.plot(timeperiods, X1, label='X')
            #plt.axhline(y=kbar, color='r')
            plt.subplot(2,1,2)
            plt.plot(timeperiods, Y1, label='Y')
            #plt.axhline(y=ellbar, color='g')
            plt.xlabel('time')
            # plt.legend(loc=9, ncol=(nx+ny))
            plt.show()    
    
        # Generate consumption, lambda, and gamma series
        r = alpha*X[0:T]**(alpha-1)*(A[0:T]*Y[0:T])**(1-alpha)
        w = (1-alpha)*X[0:T]**alpha*(A[0:T]*Y[0:T])**(-alpha)
        c = (1-tau)*(w[0:T]*Y[0:T] + (r[0:T] - delta)*X[0:T]) + X[0:T] \
            + tau*(w[0:T]*Y[0:T] + (r[0:T] - delta)*X[0:T]) - X[1:T+1]
        # T-by-1
        Lam = (c[0:T-1]**(-gamma)*(1-tau)*w[0:T-1]) / (chi*Y[0:T-1]**theta)
        # (T-1)-by-1
        Gam = (beta*c[1:T]**(-gamma)*(1 + (1-tau)*(r[1:T] - delta))) \
            / (c[0:T-1]**(-gamma))
        # (T-1)-by-1
    
        # update values for X and Y
        Xnew = (Gam)*X[1:T]
        Ynew = (Lam)*Y[1:T]
        XY = np.append(Xnew, Ynew, axis = 1)
        x = x[0:T-1,:]
        
        if fittype == 'MVOLS':
            coeffsnew = MVOLS(XY, x)
        
        dist = np.mean(np.abs(1-XY/XYold))
        print('count ', count, 'distance', dist, 'distold', distold, \
              'damp', damp)
        
        if dist < distold:
            damp = damp*1.05
            if damp > 1.:
                damp = 1.
        else:
            damp = damp*.8

        distold = 1.*dist
    
        # update coeffs
        XYold = XY
        coeffs = (1-damp)*coeffs + damp*coeffsnew
        if count % 10 == 0:
            print('coeffs', coeffs)
    return coeffs