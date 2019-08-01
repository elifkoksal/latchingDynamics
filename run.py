"""
@author: Elif Koksal Ersoz

This file is for running simulations in simu_E_g.py file for a given combination.  
    """
import numpy as np
import simu_E_g
import time
import os

def runSimulation(arguments):
    simulationspercombination = arguments[-1]
    for n in range(int(simulationspercombination)):
        simu_E_g.simulation(N, simTime, eta, tau, u, mu, G, I, n)


if __name__ == '__main__':
    simulationspercombination = 1
    I=0.0
    N=8
    eta = 0.04 
    rho = 2.4
    tau = 300
    u = rho/tau
    S = 1/(1+rho)
    G = 0.501 
    mu = 0.2501
    time_start1 = time.time()
    print('Combination : ',[eta, tau, u,G, mu])
    simTime = 3000 
    argumentsToRun = [N, simTime, eta, tau, u, mu, G, I, simulationspercombination]                                                
    runSimulation(argumentsToRun)
    time_end1 = time.time()
    print("Completed in {} sec".format(time_end1 - time_start1))

