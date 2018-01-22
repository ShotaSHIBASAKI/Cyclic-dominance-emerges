#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 14:46:56 2017

@author: shibasakishota
"""

import matplotlib.pyplot as plt
import numpy as np
import csv
from sklearn import svm, cross_validation
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

def original_cmp(y):
    if y==0:
        return "y"
    elif y==1:
        return "m"
    elif y==2:
        return "c"
    else:
        return "g"
        

def make_meshgrid(x, y, zeroX, zeroY, h=.001):
    """Create a mesh of points to plot in

    Parameters
    ----------
    x: data to base x-axis meshgrid on
    y: data to base y-axis meshgrid on
    h: stepsize for meshgrid, optional

    Returns
    -------
    xx, yy : ndarray
    """
    x_min, x_max = min(x.min(),zeroX) - 1, x.max() + 1
    y_min, y_max = min(y.min(),zeroY) - 1, y.max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))
    return xx, yy



def plot_contours(clf, xx, yy, **params):
    """Plot the inverese-transformed decision boundaries for a classifier.

    Parameters
    ----------
    ax: matplotlib axes object
    clf: a classifier
    xx: meshgrid ndarray
    yy: meshgrid ndarray
    params: dictionary of params to pass to contourf, optional
    """
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    out = plt.contourf(xx, yy, Z, **params)
    return out

def plot_contours_tr(clf, xx, yy, xx_tr, yy_tr, **params):
    """Plot the decision boundaries for a classifier.

    Parameters
    ----------
    ax: matplotlib axes object
    clf: a classifier
    xx: meshgrid ndarray
    yy: meshgrid ndarray
    params: dictionary of params to pass to contourf, optional
    """
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    out = plt.contourf(xx_tr, yy_tr, Z, **params)
    return out

with open('exp_result.csv','r') as f:
    data=csv.reader(f)
    header=next(data)
    counter=0
    c0x=[]
    c0y=[]
    c1y=[]
    c1x=[]
    c2y=[]
    c2x=[]
    c3x=[]
    c3y=[]
    X=np.empty((0,2),float) #input data
    y=np.array([])
    C0X=np.array([])
    C0Y=np.array([])
    C1X=np.array([])
    C1Y=np.array([])
    C2X=np.array([])
    C2Y=np.array([])
    C3X=np.array([])
    C3Y=np.array([])
    for row in data:
        
        counter+=1
        if row[2]=="0": #neither FBs nor MCs are formed
            c0x.append(row[0])
            c0y.append(row[1])
            
            X=np.append(X, np.array([[float(row[0])*0.001,float(row[1])]]), axis=0)
            y=np.append(y, np.array([[int(row[2])]]))
            C0X=np.append(C0X, np.array([[float(row[0])*0.001]]))
            C0Y=np.append(C0Y, np.array([[float(row[1])]]))
            
        
        elif row[2]=="1": #only FBs are formed
            c1x.append(row[0])
            c1y.append(row[1])
            C1X=np.append(C1X, np.array([[float(row[0])*0.001]]))
            C1Y=np.append(C1Y, np.array([[float(row[1])]]))
            X=np.append(X, np.array([[float(row[0])*0.001,float(row[1])]]), axis=0)
            y=np.append(y, np.array([[int(row[2])]]))
                
        elif row[2]=="2": # only MCs are formed  
            c2x.append(row[0])
            c2y.append(row[1])
            X=np.append(X, np.array([[float(row[0])*0.001,float(row[1])]]), axis=0)
            y=np.append(y, np.array([[int(row[2])]]))
            C2X=np.append(C2X, np.array([[float(row[0])*0.001]]))
            C2Y=np.append(C2Y, np.array([[float(row[1])]]))
                
        else:
            #both FBs and MCs are formed  
            c3x.append(row[0])
            c3y.append(row[1])
            X=np.append(X, np.array([[float(row[0])*0.001,float(row[1])]]), axis=0)
            y=np.append(y, np.array([[int(row[2])]]))
            C3X=np.append(C3X, np.array([[float(row[0])*0.001]]))
            C3Y=np.append(C3Y, np.array([[float(row[1])]]))

# Scaling
sc = StandardScaler()
sc.fit(X)
X=sc.transform(X)

#grid research of parameter C
C_list=np.array([0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000])
CV=np.zeros(np.size(C_list))
MaxScore=-10000
cand=-1000
for i in range(np.size(C_list)):
    # data since we want to plot the support vectors
    C = C_list[i]  # SVM regularization parameter
    clf=svm.LinearSVC(C=C, dual=False, multi_class='ovr')
    clf.fit(X, y)
    # tcross validation
    scores = cross_val_score(clf, X, y, cv=10) # cv number of splits
    mscore=scores.mean()#mean of score
    CV[i]=mscore
    if MaxScore<mscore:
        MaxScore=mscore
        cand=C
print(str("parameter value %.4f" %(cand)))
print(str("cross validation %.4f" %(MaxScore)))

#plot the cross validation
plt.xscale("log")
plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.3f'))
plt.plot(C_list,CV)
plt.xlabel("parameter value",fontsize=16)
plt.ylabel("accuracy",fontsize=16)
plt.savefig("CrossValidation.pdf")
plt.show()

#after find the optimal parameter value
clf=svm.LinearSVC(C=cand, dual=False, multi_class='ovr')
clf.fit(X, y)
#inverse transform
X_tr=sc.inverse_transform(X)
#print(clf)
# Set-up 2x2 grid for plotting.
X0, X1 = X[:, 0], X[:, 1]
X0_tr, X1_tr=X_tr[:,0], X_tr[:,1]
zeroX,zeroY=sc.transform([0,0])
xx, yy = make_meshgrid(X0, X1, zeroX, zeroY)
xx_tr, yy_tr=np.meshgrid(np.linspace(0, np.max(X0_tr), np.size(xx,1)),
                         np.linspace(0, np.max(X1_tr), np.size(yy,0)))
#plot_contours(clf, xx, yy, cmap=plt.cm.cool, alpha=0.8)
plot_contours_tr(clf, xx, yy, xx_tr, yy_tr, cmap=plt.cm.cool, alpha=0.8)
#plt.scatter(C1X,C1Y,color=plt.cm.ocean,s=50,edgecolors='k')
#plt.scatter(X0, X1, c=y, cmap=plt.cm.cool, s=50, edgecolors='k')
plt.scatter(X0_tr, X1_tr, c=y, cmap=plt.cm.cool, s=50, edgecolors='k')
plt.xlim(0, max(X0_tr))
plt.ylim(0, max(X1_tr))
plt.xlabel("luminosity x 1000 (lux)",fontsize=16)
plt.ylabel("amount of liquid (ml)",fontsize=16)
plt.savefig("ResultClassification_dot.pdf")
plt.show()


print(mscore)
