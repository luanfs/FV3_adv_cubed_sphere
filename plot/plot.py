import numpy as np
import matplotlib.pyplot as plt
from plotting_routines  import plot_scalarfield


# Python script to plot the outputs
# directory
graphdir ='../graphs/' # must exist
datadir  ='../data/' # must exist
nbfaces = 6
figformat='png'

# some constants
N    = 96 # number of cells
tc   = 1  # test case
hord = 8 # 1d adv scheme
dp   = 2  #2d adv scheme
iadv = 2
gtype = 2
mf = 0
nplots = 12


if tc==1:
   qmin, qmax = 1000.0, 3000.0
elif tc==2 or tc==3:
   qmin =  0.0
   qmax =  1.0
elif tc==4:
   qmin =  0.0
   qmax =  3.7

map_projection='mercator'
dpi=100
figformat='png'
colormap='jet'


basename = "g"+str(gtype)+"_tc"+str(tc)+"_N"+str(N)+"_hord"+str(hord)+"_iadv"+str(iadv)+"_dp"+str(dp)+"_mf"+str(mf)+"_t"
# Get scalar field
q = np.zeros((N,N,6,nplots+1))
for t in range(0, nplots+1):
   for p in range(0, nbfaces):
      # basename for plotting
      input_name  = datadir+basename+str(t)+'.txt'
      output_name = graphdir+'adv_cs_'+basename+str(t)+'.'+figformat
      data_info = np.loadtxt(input_name)

      # get binary data
      input_name_bin  = datadir+basename+str(t)+'_face'+str(p+1)+'.dat'

      z = open(input_name_bin, 'rb')
      z = np.fromfile(z, dtype=np.float64)
      z = np.reshape(z, (N,N))
      q[:,:,p,t] = z

      # get info
      time = data_info[0]
      time = str("{:.2e}".format(time))

      massvar = data_info[1]
      massvar = str("{:.2e}".format(massvar))

      cfl = data_info[2]
      cfl = str("{:.2e}".format(cfl))

   if iadv==1:
     sp = 'PL'
   elif iadv==2:
     sp = 'LT'
   
   output_name = graphdir+'advcs_'+basename+str(t)+'.'+figformat
   title = "N="+str(N)+", time = "+time+" days, CFL="+cfl+ '\n'+ \
   sp +'.hord'+str(hord)+'.dp'+str(dp)

   # plot the graph
   #plot_scalarfield(q[:,:,:,t], map_projection, title, output_name, colormap, qmin, qmax, dpi, figformat)


#------------------------------------------------------------------------------------------------
# Plot the error
if tc<=1:
   ts = 0
   te = nplots+1
elif tc>=2:
   ts = nplots
   te = nplots+1

q_error = np.zeros((N,N,6,te))
for t in range(ts,te):
   # basename for plotting
   input_name  = datadir+basename+str(t)+'.txt'
   output_name = graphdir+'adv_cs_'+basename+str(t)+'.'+figformat
   data_info = np.loadtxt(input_name)

   # get info
   time = data_info[0]

   q_error[:,:,:,t] = q[:,:,:,t] - q[:,:,:,0]

emax = np.amax(abs(q_error))
emin = -emax
colormap='seismic'
for t in range(ts,te):
   # basename for plotting
   input_name  = datadir+basename+str(t)+'.txt'
   output_name = graphdir+'adv_cs_'+basename+str(t)+'.'+figformat
   data_info = np.loadtxt(input_name)

   # get info
   time = data_info[0]
   time
   time = str("{:.2e}".format(time))

   massvar = data_info[1]
   massvar = str("{:.2e}".format(massvar))

   cfl = data_info[2]
   cfl = str("{:.2e}".format(cfl))
   title = "N="+str(N)+", time = "+time+" days, CFL="+cfl+ '\n'+ \
   sp +'.hord'+str(hord)+'.dp'+str(dp)
   output_name = graphdir+'advcs_error_'+basename+str(t)+'.'+figformat
   plot_scalarfield(q_error[:,:,:,t], map_projection, title, output_name, colormap, emin, emax, dpi, figformat)
