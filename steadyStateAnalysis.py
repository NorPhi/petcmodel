# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 11:42:21 2014

Calculate the steady state of the system when state transition are switched off


Copyright (C) 2014-2015  Anna Matuszyńska, Oliver Ebenhöh, Philipp Norf

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (license.txt).  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np

import sys
import dill
from functools import partial

import petcModel
import parametersPETC
import simulate

from time import time, sleep

p = parametersPETC.ParametersPETC()
# ---------------------------------------- #
# switch off the state transitions #
p.staticAntI = 0
p.staticAntII = 0
p.kStt7 = 0
p.kPph1 = 0

m = petcModel.PETCModel(p)
s = simulate.Sim(m)

# Length of PFDrange and STrange set through a command line argument:
if len(sys.argv) == 3 and sys.argv[2].isdigit() and int(sys.argv[2]) >= 3:
    N = int(sys.argv[2])
else:
    N = 20 # default. This also prevent input errors.

PFDrange = np.linspace(75,1500,N)
STrange = np.linspace(0,1,N)

Ys = np.zeros([N,N,8])

# dark adapted state. Not important         #[[]] see line 51
y0=np.array([[p.PQtot, 0.0202, 5.000, 0.0000, 0.0000, 0.0001, 0.9, 0.0000]])

# This replaces the part of y0[6] = STrange[i] in the previous loop.
y0 = y0.repeat(N,0) # this works because of the double brackets [[]] in y0
y0[:,6] = STrange

# Fix PFDrange parameter
# This might be unnecassary, pathos.multiprocessing can deal with functions
# which require multiple arguments. However, I haven't tested if it can deal
# with different numbers of arguments, yet.
steadyState = partial(s.steadyStateLightScan,PFDrange)

# --------------------------------------- #
#  This approach uses only a single core  #
'''
clock = time() # starting time
for i in range(N):
    print('Teraz liczymy dla ' + str(i)) # I don't get this line...
    Y = steadyState(y0[i])
    Ys[i,:] = Y
clock = time()-clock # finishing time for N iterations

# Serilisation of results and stats: 
ss = {'STrange': STrange, 'PFDrange': PFDrange, 'Ys': Ys, 'Sec': clock}
output = open('steadyStateAnalysisFixedST_SC_N' + str(N) + '.pkl', 'wb')
dill.dump(ss,output,2)
output.close()

print(clock)
'''
# --------------------------------------- #
#    This approach uses multiple cores    #
from pathos.multiprocessing import ProcessingPool, cpu_count

if __name__ == '__main__': # This is essential if Used on windows!
    clock = time() # starting time
    
    # creates a worker pool from given comand line parameter. If the given
    # parameter is to large all detectable CPUs will be utilised. If the given
    # parameter is nonsense only 1 core will be utilized.
    workers = 1
    if len(sys.argv) >= 2 and sys.argv[1].isdigit() and int(sys.argv[1]) > 0:
        workers = cpu_count()
        if int(sys.argv[1]) <= workers:
            workers = int(sys.argv[1])
    
    print 'N:  ' + str(N)
    print 'PW: ' + str(workers)
    sleep(3) # just 3 seconds pause to read the input again.

    # All the magic happens here:
    pool = ProcessingPool(workers)
    Ys = pool.map(steadyState,y0)   

    clock = time()-clock # elapsed time
    print 'Seconds: ' + str(clock) # Not essential but useful.

    # Serilisation of results and stats:
    ss = {'STrange': STrange, 'PFDrange': PFDrange, 'Ys': Ys, 'Sec': clock, 'PoolWorkers': workers}
    output = open('steadyStateAnalysisFixedST_MC_N' + str(N) + '.pkl', 'wb')
    dill.dump(ss,output,2)
    output.close()

else:
    print('Well, something went wrong.')

#================================================================= #
# 3 D plotting routine to obtain figure as in Ebenhoeh et al. 2014 #
'''
import dill
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

file = open('steadyStateAnalysisFixedST_MC_N20.pkl', 'rb')
data = dill.load(file)


ST = data['STrange']
PFD = data['PFDrange']
X,Y = np.meshgrid(PFD, ST)

Yss = np.zeros([len(ST),len(data['Ys'])])

for i in range(len(ST)):
        for j in range(len(PFD)):
         Yss[i,j] = 1 - data['Ys'][i][j][0] / 17.5

cm = matplotlib.cm.get_cmap('jet')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X,Y,Yss, rstride=1, cstride=1, cmap=cm,
                       linewidth=1, antialiased=True)

y_formatter = matplotlib.ticker.ScalarFormatter(useOffset = False)
ax.yaxis.set_major_formatter(y_formatter)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

ax.set_zlim3d(0, 1)
fig.colorbar(surf) #, shrink=0.5, aspect=5)
plt.title('steady state of reduced plastoquinon pool')
plt.xlabel('PFD')
plt.ylabel('PSII cross section')
plt.show()
'''