import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path
import subprocess

# Python script to plot the outputs
# directory
graphdir ='../graphs/'
datadir  ='../data/'
pardir   ='../par/'
figformat = 'png'

# Program to be run
program = "./main"

# test case
tc = 2

# advection scheme
hords = (0,8)

#grid types (0-equiedge, 2-equiangular)
gtypes = (2,2)

# dp scheme
dps = (1,2)

# inner adv scheme
iadvs = (1,2)

# mass fixers
mfs = (1,1)

# N
N = 48
Ns=[]

# number of grids
ngrids = 3

# error arrays
errors_linf = np.zeros((ngrids, len(dps), len(hords)))
errors_l1   = np.zeros((ngrids, len(dps), len(hords)))
errors_l2   = np.zeros((ngrids, len(dps), len(hords)))

# compile the code
for n in range(0, ngrids):
   for m in range(0, len(hords)):
      hord = hords[m]
      for k in range(0, len(dps)):
         iadv  = iadvs[k]
         dp = dps[k]
         mf = mfs[k]
         gtype = gtypes[k]

         # error filename
         filename = "g"+str(gtype)+"_tc"+str(tc)+"_N"+str(N)+"_hord"+str(hord)+"_iadv"+str(iadv)+"_dp"+str(dp)+"_mf"+str(mf)+"_errors.txt"
   
         # load the errors
         errors = np.loadtxt(datadir+filename)
         errors_linf[n,k,m] = errors[0]
         errors_l1  [n,k,m] = errors[1]
         errors_l2  [n,k,m] = errors[2]

   # update N and dt
   Ns.append(N)
   N = 2*N
   

# Plotting parameters
colors = ('green', 'red', 'blue','orange','green', 'red', 'blue')
markers = ('o','+','x','*', '+', 'x', '*', '+')
lines_style = ('-','--','-','--')
lines_style = ('-','-','-','-','--','--','--')

# Plot error graph 
errors = [errors_linf, errors_l1, errors_l2]
names = [r'$L_{\infty}$',r'$L_1$',r'$L_2$']
enames = ['linf','l1','l2']

for e in range(0, len(errors)):
   Error = errors[e]
   emin, emax = np.amin(Error), np.amax(Error)
   emax = 10**(np.floor(np.log10(emax)+1))
   emin = 10**(np.floor(np.log10(emin)-1))
   fig, axs = plt.subplots(1, 2, figsize=(12, 6))  # Creating subplots, 1 row, 2 columns
   title = names[e] + " error - TC" + str(tc)
   fig.suptitle(title)
   for m in range(0, len(hords)):
      error = Error[:,:,m]
      hord = hords[m]
      for k in range(0, len(dps)):
         iadv = iadvs[k]
         dp = dps[k]         

         if iadv==1:
           sp = 'PL'
         elif iadv==2:
           sp = 'LT'
         scheme = sp+'-dp'+str(dp)+'-hord'+str(hord)

         # convergence rate
         n = len(Ns)-1
         CR = (np.log(error[n-1,k])-np.log(error[n,k]))/np.log(2.0)
         CR = str("{:2.1f}".format(CR))

         # plot in the respective subplot
         axs[m].loglog(Ns, error[:, k], lines_style[k], color=colors[k], marker=markers[k], \
                               label=str(scheme) + " - order " + CR)

      # plot reference lines
      eref = np.amin(Error)
      order1 = [eref, eref/2.0]
      order2 = [eref, eref/4.0]
      order3 = [eref, eref/8.0]
      if m==0:
         axs[m].loglog(Ns[ngrids - 2:ngrids], order1, '-',  color='black', label='1st order')
         axs[m].loglog(Ns[ngrids - 2:ngrids], order2, '--', color='black', label='2nd order')
         axs[m].loglog(Ns[ngrids - 2:ngrids], order3, '-.', color='black', label='3rd order')
      else:
         axs[m].loglog(Ns[ngrids - 2:ngrids], order1, '-',  color='black')
         axs[m].loglog(Ns[ngrids - 2:ngrids], order2, '--', color='black')
         axs[m].loglog(Ns[ngrids - 2:ngrids], order3, '-.', color='black')
 

      # Label
      title = 'hord'+str(hord)
      axs[m].set_xlabel('$N$')
      axs[m].set_ylabel('Error')
      axs[m].set_xlim(0, 1000)
      axs[m].set_ylim(emin, emax)
      axs[m].set_title(title)
      axs[m].legend()
      axs[m].grid(True, which="both")
      
   plt.tight_layout() 
   figname = graphdir + 'adv2d_tc' + str(tc) + '_' + enames[e] + "_error"
   plt.savefig(figname+'.'+figformat, format=figformat)
   plt.close()
