#import packages
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import uproot
import random
import os
import math


def intertrack(TrackMatrix):
    length= len(TrackMatrix[:,0,0])
    mean=np.zeros(3*length)
    mean=mean.reshape(3,length)
    mean_sum=0
    sum=0

    for i in range(length):
        mean[0,i]=np.mean(TrackMatrix[i,0,TrackMatrix[i,0,:]!=0])
        mean[1,i]=np.mean(TrackMatrix[i,1,TrackMatrix[i,1,:]!=0])
        mean[2,i]=np.mean(TrackMatrix[i,2,TrackMatrix[i,2,:]!=0])

        if math.isnan(mean[0,i])==True:
            mean[0,i]=mean[0,i-1]
            mean[1,i]=mean[1,i-1]
            mean[2,i]=mean[2,i-1]
    for i in range(length):
        x_i=np.array([mean[0,i],mean[1,i],mean[2,i]])
        plt.plot(mean[0,i],mean[1,i],'.')
        for j in range(length):
            if j==i:
                continue
            x_j=np.array([mean[0,j],mean[1,j],mean[2,j]])
            dist_squared=np.sqrt((x_i[0]-x_j[0])**2+(x_i[1]-x_j[1])**2+(x_i[2]-x_j[2])**2)
            sum+=dist_squared
    plt.xlabel('Dist [nm]')
    plt.ylabel('Dist [nm]')
    plt.grid()
    plt.plot(np.cos(np.linspace(0,2*np.pi,1000))*6000,np.sin(np.linspace(0,2*np.pi,1000))*6000,'--')
    #plt.show()
    print(sum,length)
    return sum/(length-1)**2
