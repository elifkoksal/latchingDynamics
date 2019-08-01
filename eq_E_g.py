# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 14:42:47 2018

@author: Carlos
"""

import numpy as np


def eq_E_g(t, x, flag, I, U, tau, mu, J, Lambda):
	
	N = J.shape
	F = np.zeros(2*N)
	S = x[N:2*N]*x[0:N]
	X = sum(x[0:N])

    # Excitatiory neuron
	for i in range(N):
    F[i] = x[i]*(1-x[i])*(-mu[0]*x[i] -I[0] - Lambda*X+ J[i,0:N]*S.flatten(1))
	
  F[N:2*N] = (1-x[N:2*N])/tau[0] - U[0]*x[0:N]*x[N:2*N]

	return F
 
if __name__ == "__main__":
	pass
