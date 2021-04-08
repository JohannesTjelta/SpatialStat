import numpy as np
import matplotlib.pyplot as plt
import random
def gearysC(IonMatrix,parameter):


    x_mean=np.array([np.mean(IonMatrix[0,IonMatrix[0,:]!=0]),np.mean(IonMatrix[1,IonMatrix[1,:]!=0]),np.mean(IonMatrix[2,IonMatrix[2,:]!=0])])
    W=0
    mean_sum=0
    sum=0
    N=len(IonMatrix[0,IonMatrix[0,:]!=0])
    print('N=',N)

    for i in range(N):
        if IonMatrix[2,i]==0:
            continue
        x_i=np.array([IonMatrix[0,i],IonMatrix[1,i],IonMatrix[2,i]])
        mean_sum+=(np.sqrt((x_i[0]-x_mean[0])**2+(x_i[1]-x_mean[1])**2+(x_i[2]-x_mean[2])**2))**2
        for j in range(N):
            if i>j:
                continue
            x_j=np.array([IonMatrix[0,j],IonMatrix[1,j],IonMatrix[2,j]])
            dist=(np.sqrt((x_i[0]-x_j[0])**2+(x_i[1]-x_j[1])**2+(x_i[2]-x_j[2])**2))**2
            dist2=np.sqrt(dist)
            if dist2>parameter:
                #print('lol')#,dist2,x_i[0],x_j[0],x_i[1],x_j[1],x_i[2],x_j[2])
                w_ij=0

            elif dist2<parameter:
                 sum+=dist
                 W+=1
            else:
                print('no')
    C=(N-1)*sum/(2*W*mean_sum)

    print(C,W)
    return C


if __name__=='__main__':
    parameterArray=np.array((0.3,0.7,10))
    for parameter in parameterArray:
        points =20
        #parameter=10
        points3D=points*3
        #rm=random matrix. rmc= random matrix clusterd
        rm=np.random.uniform(0,1,points3D).reshape(3,points)
        rmc=np.random.uniform(0.3,1,points3D).reshape(3,points)
        rmc[1]=rmc[1]*-1
        rmc2=np.random.uniform(0.3,0.2,points3D).reshape(3,points)
        def xy(n):
            x=np.ones(n)
            y=np.ones(n)
            i=0
            while i<n:
                x_temp = random.uniform(-3,3)
                y_temp = random.uniform(-3,3)
                dist=np.sqrt(x_temp**2+y_temp**2)
                if dist>3:
                    continue
                else:
                    x[i]=x_temp
                    y[i]=y_temp
                    i+=1
            return x,y
        def plot_circle():
            n=2000
            whole_pi=np.linspace(0,2*np.pi,n)
            c1=np.zeros(n)
            c2=np.cos(whole_pi)*3
            c3=np.sin(whole_pi)*3
            return c1,c2,c3
        qwe=plot_circle()
        xy=xy(points)

        rm[0]=xy[0];rm[1]=xy[1]
        c1=gearysC(rmc2,parameter)
        c2=gearysC(rmc,parameter)
        c3=gearysC(rm,parameter)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        #ax=fig.add_subplot(211,projection='3d')
        ax.plot3D(rm[0],rm[1],rm[2],'.',label='Randomly distrebuted C={:1.3f}'.format(c3))
        ax.plot3D(rmc[0],rmc[1],rmc[2],'.',label='Semi clustering C={:1.3f}'.format(c2))
        ax.plot3D(rmc2[0],rmc2[1],rmc2[2],'.',label='Max clustering C={:1.3f}'.format(c1))
        #ax.plot3D(qwe[0],qwe[1],qwe[2],'--',color='gray')
        ax.plot3D(qwe[1],qwe[2],qwe[0],'--',color='gray')
        ax.plot3D(qwe[1],qwe[2],np.ones(2000),'--',color='gray')
        #ax.plot3D(qwe[2],qwe[0],qwe[1])
        plt.title('Geary\'s C with a cut parameter={}'.format(parameter))
        for label in (ax.get_xticklabels() + ax.get_yticklabels() +ax.get_zticklabels()):
    	       label.set_fontsize(13)
        ax.legend()
        plt.savefig('Plots/qualitativeC{}.png'.format(parameter),bbox_inches='tight')
        plt.show()
