# -*- coding: utf-8 -*-
"""
@author: Elif Koksal Ersoz
The system of N units is simulated using the Euler-Maruyama method for dt =0.01.
Units above the threshold 0.5 are considered to be active. To faciitate the analysis, solutions are filtered (smoothened).
Passages from the threshold are recorded for furhter analysis.

"""
import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.signal import savgol_filter
from matplotlib import colors as mcolors

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

################################################
##### Synaptic plasticity in E-E matrices ######
################################################

def create_folder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def writeActivityPattern(working_txt_patternActivity, activityPattern,activityTime): 
    with open(working_txt_patternActivity, "w") as output:
        indexActivityTime = 0
        for p in activityPattern[1:]:
            output.write('{!s};{!s}'.format(p, activityTime[indexActivityTime]))
            output.write("\n")
            indexActivityTime += 1
        output.close()

def connection(N, pattern):
	wCons = np.zeros((N, N))

	for h in range(N):
		for j in range(N):
			for p in range(np.size(pattern, axis=1)):
				wCons[h, j] = wCons[h, j]+(pattern[h, p])*(pattern[j, p])
				
	return wCons

def find_in_iterable(x, iterable):
    for i, item in enumerate(iterable):
        if item == x:
            return i
    return None

def model(t, x, I, U, tau, mu, J, G):
    # inhibition by the excitatory population 
    # Excitatiory neuron
    N = J.shape[0]
    dxdt = np.zeros((1, 2*N))

    S = x[N:2*N]*x[0:N]
    sigma_x = sum(x[0:N])
    i=1
    
    for i in range(N):
	    dxdt[0][i] = x[i]*(1-x[i])*(-mu*x[i] -I - G*sigma_x+ np.dot(J[i,0:N], S.flatten("F")))
	
	# S
    dxdt[0][N:2*N] = (1 - x[N:2*N])/tau - U*x[N:2*N]*x[0:N]
    #for i in range(N, 2*N):
    #    dxdt[0][i] = (1 - x[i])/tau - U*x[i]*x[i-N]
    dxdt = dxdt.reshape((1, 2*N))
    return dxdt

def dW(delta_t, N):
    return np.sqrt(delta_t)*(1-2*np.random.random_sample((1,N)))

def simulation(N, simTime, eta, Tau, u, mu, G, I, n):
	
    current_folder = os.path.dirname(os.path.realpath(__file__))

    # define folders, file names	

    folder_name = "scaleG_N=%d"%(int(N))+ "_G=%d"%(int(1000*G))+ "_eta=%d"%(int(1000*eta))+\
            "_tau=%d"%(int(Tau)) + "_U=%d"%(int(1000*u)) 
    working_folder = os.path.join(current_folder, folder_name)
    create_folder(working_folder)
        
    name_file = folder_name + "_lmu=%d"%(int(1000*mu)) +"_%d"%(int(n))
    
    working_png_file = os.path.join(working_folder, name_file + '.png')
    working_txt_patternActivity= os.path.join(working_folder, name_file+'.txt')

    # asymétrique
    M=2 # number of overlaps
    pt = N-2 # N-1: cycle, N-2:chain
    pattern1 = np.concatenate((np.ones((M,1)), np.zeros((N-M,1))), axis=0)
    pattern = pattern1.copy()
		
    for pi in range(pt):
	    pattern1= np.concatenate(([pattern1[-1]], pattern1[0:-1]), axis=0)
	    pattern=np.concatenate((pattern, pattern1), axis=1)

    # Jmax-Hebb rule
    Jmax = connection(N, pattern)
        
    dt = 0.01
    t0 = 0.0
    ts  = np.arange(t0, simTime, dt)
    t = 1.0
    
	# Initial conditions
    x0 = np.empty(2*N)
    x0[0:N] = pattern[:,0].conj().transpose() # E neurons
    x0[N:2*N] = np.ones((1, N)) # S values
	
    xComplete = np.zeros((int(ts.size),2*N))	
    xComplete[0] = x0.copy()	
    # List to record neuron activities
    liste_neurones = [[] for _ in range(N)]
    activityPattern = [[0,1]]
    # pattern durations
    activityTime = []
    patternList=[ ]
    for i in range(N):
        patternList.append([i, i+1])
    patternList.append([ ])
    #####################
    # Simulate the system
    # ##################		
    Winc = np.empty(2*N)
    #noise = np.zeros((int(ts.size),1))

    for simCount in range(1, ts.size):
        
        t = (simCount - 1)*dt

        Winc[0:N]  = dW(dt, N)
        Winc[N:2*N] = np.zeros(N)
        
        xComplete[simCount] = xComplete[simCount-1] + model(t,xComplete[simCount-1], I, u, Tau, mu, Jmax, G) * dt + eta * Winc
        
        # xi values are bounded by [0,1]
        for i in range(N):
            if xComplete[simCount][i]<0:
                xComplete[simCount][i]=-xComplete[simCount][i]
            elif xComplete[simCount][i]>1:
                xComplete[simCount][i]=2-xComplete[simCount][i]
    
    # find active neurons which are above 0.5 after applying a filter

    xCompleteFilterred = np.zeros((int(ts.size),N))
    for i in range(N):
        xCompleteFilterred[:,i] = savgol_filter(xComplete[:,i], 1001, 2)
    
    for simCount in range(ts.size):
        liste_neurones.append([])
        for  countActive, x_i in np.ndenumerate(xCompleteFilterred[simCount]):
            if x_i > 0.5:
                liste_neurones[-1].extend(list(countActive))

        if liste_neurones[-1] != liste_neurones[-2]:
            if len(liste_neurones[-1]) < len(liste_neurones[-2]): 
                activityTime.append([simCount*dt,''])
                activityPattern.append(liste_neurones[-1])
            elif len(liste_neurones[-1]) > len(liste_neurones[-2]): 
                activityTime.append(['',simCount*dt])
                activityPattern.append(liste_neurones[-1])

    writeActivityPattern(working_txt_patternActivity, activityPattern, activityTime)
###########################
###########################
###########################


### plot figure ####
    plt.figure(n)
    #plt.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.965, wspace=0.2, hspace=0.2)
    colorList = ['blue','red', 'orange','purple','g','cyan','firebrick','magenta']
    plt.subplot(211)
    for i in range(N):
      plt.plot(ts, xComplete[:,i], linewidth=1.0, color = colorList[i], label="{}".format(i+1))
    plt.xlabel('Time (ms)')
    plt.ylabel('x(t)')
    plt.legend()

    plt.subplot(212)
    for i in range(N, 2*N):
      plt.plot(ts, xComplete[:,i], linewidth=1.0, label="{}".format(i+1-N))
    plt.xlabel('Time (ms)')
    plt.ylabel('s(t)')
    plt.legend()
    plt.savefig(working_png_file)


