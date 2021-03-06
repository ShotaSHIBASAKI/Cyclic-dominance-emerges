#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last update, 16/01/2018
Discrete Replicator Dynamic in Asexual model
@author: Shota Shibasaki
"""

import numpy as np
import matplotlib.pyplot as plt

#parameter
a1=1
b1=1.5
e1=0.0001
a2=0.5
b2=1.2
e2=0.0001
e0=0.000001
h=0.5#effect of forming smaller fruiting bodies
theta=0.5#prob of macrocyst formation
def Matrix(a2,b2,theta):  
    """
   Defining patoff matrixes
   A1: payoff matrix for fruiting body formation
   A2: payoff matrix for macrocyst formation
   A: patoff matrix generated by the mixture of A1 and A2.
      see Eq [2] in the main text for more detail.
   Parameters
   a2: benefit of cmutual ooperation in macrocyst formation
   b2: benefit of exploitation for defectors in macrocyst formation
   theta: prob. of macrocyst fromation
      
   """
    A1=np.array([[a1,h*a1,0],\
                 [e0,e0,e0],\
                 [b1,h*e1,e1]])
    
    A2=np.array([[a2,0,a2],\
                 [b2,e2,b2],\
                 [a2,0,a2]])
    A=(1-theta)*A1+(theta)*A2
    return A

def func(x, t, a2,b2,theta):
    """
    Discrete replicator dynamics
    x: array of the frequency of the three strategy
    t: time
    a2, b2, theta: defined in function Matrix
    Note the background fitness = 0 here
    """
    #x is array and A is matrix
    A=Matrix(a2,b2,theta)
    phi=np.dot(x,A)
    phi=np.dot(phi,x.T)#mean fitness in the population
    #calculate replicator dynamics
    f0=np.dot(A[0,:],x.T)/phi
    f1=np.dot(A[1,:],x.T)/phi
    f2=np.dot(A[2,:],x.T)/phi     
    f=np.array([x[0]*f0,x[1]*f1,x[2]*f2])
    return(f)

def triplot(sol,col):
    """
    Mapping the phase space into the triangl plain
    """
    TP=np.array([1/2,np.sqrt(3)/2])#C
    LS=np.array([0,0])#D_M
    RS=np.array([1,0])#D_F
    orbX=sol[:,0]*TP[0]+sol[:,1]*LS[0]+sol[:,2]*RS[0]
    orbY=sol[:,0]*TP[1]+sol[:,1]*LS[1]+sol[:,2]*RS[1]
    #plot the dynamics
    plt.xlim(-0.2,1.2)
    plt.ylim(-0.1,1.2*TP[1])
    plt.plot([LS[0],TP[0]],[LS[1],TP[1]],linewidth=4,color="k")
    plt.plot([RS[0],TP[0]],[RS[1],TP[1]],linewidth=4,color="k")
    plt.plot([RS[0],LS[0]],[RS[1],LS[1]],linewidth=4,color="k")
    plt.text(-0.15,0,'$D^M$',fontsize=20,)
    plt.text(1.05,0,'$D^F$',fontsize=20,)
    plt.text(0.475,0.9,'$C$',fontsize=20,)
    plt.plot(orbX,orbY,linewidth=2,color=col)
    plt.tick_params(labelbottom='off',labelleft="off")
    
def Stability(A,p):
    """
    Check the stability of the interior point G 
    in the discrete replicator dynamics
    A is payoff matrix
    p is the frequency of three strategies
    """
    Ave=(p.dot(A)).dot(p.T)#mean fitness given p
    ev1=np.array([1,0,0]).T
    ev2=np.array([0,1,0]).T
    ev3=np.array([0,0,1]).T
    """
    L: Linearization of payoff A given p
    Then, check the stability of the interipor point G.
    Note: there is a proof of  the necessary and sufficient condition 
    that for a cubic equation  all the absolute values of its solutions are
    less than or equal to 1.
    See Sato (2013), "On the necessary and sufficient condition that 
    for a cubic equation  all the absolute values of its solutions are
    less than or equal to 1." , Network and Information, 21, pp87-92.
    But, this paper is written in Japanese.
    """
    L=np.zeros([3,3])
    L[0,0]=1+p[0]*(A[0,0]-(p.dot(A)).dot(ev1))/Ave
    L[0,1]=0+p[0]*(A[0,1]-(p.dot(A)).dot(ev2))/Ave
    L[0,2]=0+p[0]*(A[0,2]-(p.dot(A)).dot(ev3))/Ave
    L[1,0]=0+p[1]*(A[1,0]-(p.dot(A)).dot(ev1))/Ave
    L[1,1]=1+p[1]*(A[1,1]-(p.dot(A)).dot(ev2))/Ave
    L[1,2]=0+p[1]*(A[1,2]-(p.dot(A)).dot(ev3))/Ave
    L[2,0]=0+p[2]*(A[2,0]-(p.dot(A)).dot(ev1))/Ave
    L[2,1]=0+p[2]*(A[2,1]-(p.dot(A)).dot(ev2))/Ave
    L[2,2]=1+p[2]*(A[2,2]-(p.dot(A)).dot(ev3))/Ave
    trL=L.trace()
    detL=np.linalg.det(L)
    ML=(L[0,0]*L[1,1]-L[0,1]*L[1,0])\
       +(L[1,1]*L[2,2]-L[1,2]*L[2,1])+(L[2,2]*L[0,0]-L[0,2]*L[2,0])
    if abs(-trL)<3 and abs(-trL-detL)-1<ML and ML<1+(-trL)*(-detL)-(-detL)**2:
        stability="stable"
    else:
        stability="unstable"
    #numerical calculation of the eigenvalue of matrix L for double check
    (eig,v)=np.linalg.eig(L)
    return [stability,eig]
                


def main():
    """
    Simulate the continuous replicator dymnamics in asexual model
    from the three different initial conditions
    """
    A=Matrix(a2,b2,theta) #define payoff matrix
    T=500# end the simulation at T
    init=np.array([10,1,1])#initial condition 1
    init=init/sum(init)
    p=np.empty((0,3))
    p=np.append(p, [init], axis=0)
    #Note that this code is for Discrete replicator dynamics
    for t in range(T):
        q=p[t]
        pnext=func(q,t,a2,b2,theta)
        p=np.append(p, [pnext], axis=0)
    triplot(p,"c")
    
    init=np.array([1,10,1])#initial condition 2
    init=init/sum(init)
    p=np.empty((0,3))
    p=np.append(p, [init], axis=0)
    for t in range(T):
        q=p[t]
        pnext=func(q,t,a2,b2,theta)
        p=np.append(p, [pnext], axis=0)
    triplot(p,"m")
    
    init=np.array([1,1,12])#initial condition 3
    init=init/sum(init)
    p=np.empty((0,3))
    p=np.append(p, [init], axis=0)
    for t in range(T):
        q=p[t]
        pnext=func(q,t,a2,b2,theta)
        p=np.append(p, [pnext], axis=0)
    triplot(p,"y")
    """
    Check the inteirior equilibrium G
    """
    GX=p[-1,0]*0.5+p[-1,1]*0+p[-1,2]*1
    GY=p[-1,0]*np.sqrt(3)/2
    plt.plot(GX,GY,'k*',markersize=15)
    plt.text(GX+0.05,GY,'G',fontsize=20)
    print(GX,GY)
    plt.savefig("Dic_rep_asexual.pdf")
    plt.show()
    #check the stability of the interior equilibrium G
    [stability,eig]=Stability(A,p[-1,:])
    print(stability)
    meig=abs(eig)
    #check the maximum absolute value of eigenvalue is smaller than 1
    if max(meig)<1:
        print("true")
    else:
        print(max(meig))
if __name__ == '__main__':
    main()    
    