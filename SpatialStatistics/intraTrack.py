
#import packages
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import uproot
import random
import os
#One proton
def intratrack(TrackMatrix):
        mean_sum=0
        sum=0
        TrackMatrix=np.reshape(TrackMatrix,newshape=(3,14000))
        N=len(TrackMatrix[0,TrackMatrix[0,:]!=0])
        #print('N =',N,'resulting in ',N**N,' calculations')
        for i in range(N):
            #if TrackMatrix[2,i]==0:
            #    continue
            x_i=np.array([TrackMatrix[0,i],TrackMatrix[1,i],TrackMatrix[2,i]])
            #mean_sum+=np.sqrt((x_i[0]-x_mean[0])**2+(x_i[1]-x_mean[1])**2+(x_i[2]-x_mean[2])**2)
            for j in range(N):
                if j==i:
                    continue
                x_j=np.array([TrackMatrix[0,j],TrackMatrix[1,j],TrackMatrix[2,j]])
                dist_squared=np.sqrt((x_i[0]-x_j[0])**2+(x_i[1]-x_j[1])**2+(x_i[2]-x_j[2])**2)
                sum+=dist_squared

        intertrackvalue=sum/(N**2)
        return intertrackvalue
