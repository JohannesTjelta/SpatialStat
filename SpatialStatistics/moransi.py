
#import packages
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import uproot
import random
import os
from timeit import default_timer as timer

from gearysC import gearysC
from intraTrack import intratrack
from interTrack import intertrack

# the libraries from line 12-15 needs to be installed.
# pysal espesially needs the conda enviroment to be installed
# to avoid headake drop these librarries and comment out

#import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pysal.explore as ps
import libpysal

print('whawt do you want to calculate?')
print('Moran\'s I=0,Geary\'s C=1, intraTrack=2, interTrack=3')
Method =int(input())
def prot_pos(N):
    """
    This function is distrebuting protons within the
    cellular core, N is the number of protons based on
    simulations done in the main.py program.
    """
    i = 0
    x_prot =np.zeros(N)
    y_prot=np.zeros(N)
    while i<N:
        x_temp=random.uniform(-6,6)  # edge values is for a 6um diameter cell
        y_temp=random.uniform(-6,6)
        dist_temp=np.sqrt(x_temp**2+y_temp**2)
        if dist_temp>6:
            x_temp=None
            y_temp=None
        else:
            x_prot[i]=x_temp
            y_prot[i]=y_temp
            x_temp=None
            y_temp=None
            i+=1

    return x_prot,y_prot

def protonChoice(N,path):
    """
    Extract the proton tracs from the proton tracks generated
    in Geant4 for different energies, the path is the directory
    containing the proton data.
    """
    root_list=[]
    root_files=os.listdir(path)
    if 'Res'in root_files:
        root_files.remove('Res')
    if 'analasis' in root_files:
        root_files.remove('analasis')
    i=0
    while i<N:
        root_temp=np.random.choice(root_files)
        file = uproot.open(path+'/'+root_temp)['microdosimetry']
        x=file['x'].array()
        if len(x)==0:
            None
        else:
            root_list.append(root_temp)
            i+=1
    return root_list

def ion_pos(N,path,rootfiles,x_prot,y_prot):
    """
    extract info from the proton files and assign the x and y
    possition from prot pos to the proton from the root files
    """
    dz=4000

    x_prot=x_prot*1e3
    y_prot=y_prot*1e3
    max_ion=14000  # assumed maximal ion per proton
    ion_tot=np.zeros((N,3,max_ion))
    for i in range(N):
        file = uproot.open(path+'/'+rootfiles[i])['microdosimetry']
        x=file['x'].array()
        y=file['y'].array()
        z=file['z'].array()
        type = file['flagProcess'].array()
        if len(z)==0:
            print('you did a poopy')
        else:
            x=x[np.where((z>z[0])&(z<z[0]+dz)&((type==13)|(type==23)))]
            y=y[np.where((z>z[0])&(z<z[0]+dz)&((type==13)|(type==23)))]
            z=z[np.where((z>z[0])&(z<z[0]+dz)&((type==13)|(type==23)))]
            x0=x[0]
            y0=y[0]
            z0=z[0]
            x=x-x0
            y=y-y0
            z=z-z0

            x=x+x_prot[i]
            y=y+y_prot[i]
            if np.amax(x)>6000:
                x[np.where(x>6000)]=0
                y[np.where(x>6000)]=0
                z[np.where(x>6000)]=0
            elif np.min(x)< -6000:
                x[np.where(x<-6000)]=0
                y[np.where(x<-6000)]=0
                z[np.where(x<-6000)]=0
            elif np.amax(y)>6000:
                x[np.where(y>6000)]=0
                y[np.where(y>6000)]=0
                z[np.where(y>6000)]=0
            elif np.min(y)< -6000:
                x[np.where(y<-6000)]=0
                y[np.where(y<-6000)]=0
                z[np.where(y<-6000)]=0
            if np.mean(x)>6000 or np.mean(x)<-6000 or np.mean(y)>6000 or np.mean(y)<-6000:

                x[np.where(y<-6000)]=0
                y[np.where(y<-6000)]=0
                z[np.where(y<-6000)]=0
                x[np.where(y>6000)]=0
                y[np.where(y>6000)]=0
                z[np.where(y>6000)]=0
            ion_tot[i,0,0:len(x)]=x
            ion_tot[i,1,0:len(y)]=y
            ion_tot[i,2,0:len(z)]=z
    return ion_tot,max_ion


def SpatialDataMining(ion_tot,grid_size,voxel_size):
    """
    SpatialDataMining count the number of ionizations in a voxel and sums
    the total ionizations in that voxel
    """
    x_grid=np.arange(-6000,6000,step=voxel_size)  # define grid in x 0.002
    y_grid=np.arange(-6000,6000,step=voxel_size)   # define grid in y
    print(len(x_grid))
    z_grid=np.arange(0,4000,step=voxel_size)  # define grid in z
    hit_matrix=np.zeros((len(x_grid),len(y_grid),len(z_grid)))
    ion_tot_x = 2
    for i in range(len(x_grid)-1):
        if i==0:
            start = timer()
        for j in range(len(y_grid)-1):
            #print(j)
            for k in range(len(z_grid)-1):
                hits=np.where(((ion_tot[:,0,:]>x_grid[i])&(ion_tot[:,0,:]<x_grid[i+1])&(ion_tot[:,1,:]>y_grid[j])&(ion_tot[:,1,:]<y_grid[j+1])&(ion_tot[:,2,:]>z_grid[k])&(ion_tot[:,2,:]<z_grid[k+1])))
                hit_matrix[i,j,k]=len(hits[0])
        if i==0:
            print('Time in s for one layer to be calculated:',timer()-start)
            print('Total time will approximatly be:',(timer()-start)*len(x_grid),'s')
    return hit_matrix




if __name__=='__main__':
    InputName=input('Enter energy of ptotons:')
    InputPath=input('Enter folder:')
    #InputName='1.2MeV'
    #InputPath='Data1_2'
    print('Enter 6 values for protons:')
    N_array=np.array([10,20,30,51,81,102])
    N_array=np.array([int(input('N1:')),int(input('N2:')),int(input('N3:')),
                    int(input('N4:')),int(input('N5:')),int(input('N6:'))])
    #print('First value for Geary\'s C is a test to compile the code via njit.')

    if Method==1:
        print('What is the parameter cut value?')
        parameter_limit=int(input('value [nm]:'))
    if Method==0:
        voxel=int(input('input voxelsize [nm]:'))
    meanint=np.zeros(len(N_array))
    temp=0
    #fig, axs = plt.subplots(2,3, sharex=True, sharey=True)
    fig = plt.figure(figsize=(10,6.3))
    gs = fig.add_gridspec(2, 3, hspace=0, wspace=0)
    (ax1, ax2, ax3), (ax4, ax5, ax6) = gs.subplots(sharex='col', sharey='row')
    for lol,i in enumerate(N_array):
        N=i
        path='../'+InputPath
        e = InputName
        x_prot,y_prot=prot_pos(N)
        rootfiles=protonChoice(N,path)


        ion_tot,max_ion=ion_pos(N,path,rootfiles,x_prot,y_prot)
        IonMatrix=ion_tot.reshape(3,N*max_ion)

        if Method==2:
            print(np.shape(ion_tot))
            intratrackvalue=intratrack(ion_tot)
            print('The intratrackvalue is {}'.format(intratrackvalue))
            meanint[temp]=intratrackvalue
            print(np.mean(meanint),np.std(meanint),meanint)
            temp+=1
        elif Method==1:
            print('Beginning calculating geary\'s C. This might take some time')
            C=gearysC(IonMatrix,parameter_limit)
            print('Geary\'s C for',N,e+' protons is:',C)
            print('This is with the parameter for the weight matrix set at {}nm'.format(parameter_limit))

        elif Method==0:



            hit_matrix=SpatialDataMining(ion_tot,3,voxel)
            Mean_I=np.zeros(len(hit_matrix[1,1,:]))
            for ii in range(len(hit_matrix[1,1,:])):
                w = libpysal.weights.lat2W(hit_matrix.shape[0], hit_matrix.shape[1])
                mi = ps.esda.Moran(hit_matrix[:,:,ii],w)
                #print('\nMoran\'s I for the {}\'th layer:'.format(i),mi.I)
                #print('And standard deviation for the {}\'th layer:'.format(i),mi.p_norm)
                Mean_I[ii]=mi.I

            print('Mean Moran=',np.mean(Mean_I[:-1]),'for {} protons'.format(i))
            #print('Moran\'s I for the jth layer:',mi.I)
            #np.save('export_celle_kjerne/test{}nm'.format(voxel),matlab_hit_matrix)
        elif Method==3:
            interTrackvalue=intertrack(ion_tot)
            print('The intertrackvalue is {} for {} Protons'.format(interTrackvalue,N))
        elif Method==9:
            if lol==0:
                ax1.plot(x_prot,y_prot,'.',label='1Gy')
                ax1.legend(loc='upper right')
                ax1.plot(np.cos(np.linspace(0,2*np.pi,1000))*6,np.sin(np.linspace(0,2*np.pi,1000))*6,'--')
            elif lol==1:
                ax2.plot(x_prot,y_prot,'.',label='2Gy')
                ax2.legend(loc='upper right')
                ax2.plot(np.cos(np.linspace(0,2*np.pi,1000))*6,np.sin(np.linspace(0,2*np.pi,1000))*6,'--')
            elif lol==2:
                ax3.plot(x_prot,y_prot,'.',label='3Gy')
                ax3.legend(loc='upper right')
                ax3.plot(np.cos(np.linspace(0,2*np.pi,1000))*6,np.sin(np.linspace(0,2*np.pi,1000))*6,'--')
            elif lol==3:
                ax4.plot(x_prot,y_prot,'.',label='5Gy')
                ax4.legend(loc='upper right')
                ax4.plot(np.cos(np.linspace(0,2*np.pi,1000))*6,np.sin(np.linspace(0,2*np.pi,1000))*6,'--')
            elif lol==4:
                ax5.plot(x_prot,y_prot,'.',label='8Gy')
                ax5.legend(loc='upper right')
                ax5.plot(np.cos(np.linspace(0,2*np.pi,1000))*6,np.sin(np.linspace(0,2*np.pi,1000))*6,'--')
            elif lol==5:
                ax6.plot(x_prot,y_prot,'.',label='10Gy')
                ax6.legend(loc='upper right')
                ax6.plot(np.cos(np.linspace(0,2*np.pi,1000))*6,np.sin(np.linspace(0,2*np.pi,1000))*6,'--')

    ax1.set(ylabel='Dist [$\mu m$]')

    ax4.set(xlabel='Dist [$\mu m$]', ylabel='Dist [$\mu m$]')
    ax5.set(xlabel='Dist [$\mu m$]')
    ax6.set(xlabel='Dist [$\mu m$]')
    plt.savefig('../Plots/incidentprotonexample.png',bbox_inches='tight')
