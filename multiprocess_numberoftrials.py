"""
@author: Elif Koksal Ersoz

This file is for running simulations on multiple processors (cpuNumber).  
    """
import multiprocessing
import numpy as np
import simu_E_g
import time
import os

def printArgumentsToRun(arguments):
    simulationspercombination = arguments[-2]
    stepSize = arguments[-1]
    for n in range(simulationspercombination, stepSize + simulationspercombination, 1):
        
        print(n) 

def runSimulation(arguments):

    simulationspercombination = arguments[-2]
    stepSize = arguments[-1]
    for n in range(simulationspercombination, stepSize + simulationspercombination, 1):
        simu_E_g.simulation(N, simTime, eta, tau, u, mu, G, I, n)


if __name__ == '__main__':

    liste_eta = np.arange(0.02, 0.041, 0.02)
    liste_rho = np.array([2.4, 1.2])
    liste_tau = np.array([300, 900])
    liste_G = np.arange(0.501, 0.601, 0.05)
    liste_mu = np.arange(0.0501, 0.501, 0.05)
    liste_simTime = np.array([[3000, 8000],[6500, 10000]])
    I=0.0
    N=8
    count = len(liste_eta)*len(liste_rho)*len(liste_tau)*len(liste_G)
    simulationspercombination = 100
    for eta in liste_eta:
        for rho in liste_rho:
            indexRho = liste_rho.tolist().index(rho) 
            for tau in liste_tau:
                indexTau = liste_tau.tolist().index(tau) 
                u = rho/tau
                S = 1/(1+rho)
                for G in liste_G:
                    time_start1 = time.time()
                    print('Combination : ',[eta, tau, u,G])
                    for mu in liste_mu:
                        simTime = liste_simTime[indexRho][indexTau]
                        cpuNumber = 1
                        stepSize = int(np.floor(simulationspercombination/cpuNumber))
                        totalSimulationsPerCombination = stepSize * cpuNumber
                        argumentsToRun = [(N, simTime, eta, tau, u, mu, G, I, simulationspercombinationpercpu, stepSize) for simulationspercombinationpercpu in range(0, totalSimulationsPerCombination, stepSize)]                                                
                        p = multiprocessing.Pool(processes = cpuNumber)
                        async_result = p.map_async(runSimulation, argumentsToRun)
                        p.close()
                        p.join()

                    time_end1 = time.time()
                    count -= 1 
                    print("Completed in {} sec".format(time_end1 - time_start1), 'Number of the remaining folders = ', count)

