# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 19:26:31 2019

@author: Subrata
"""
#this code is for figure_2 and figure_3
import numpy as np
from numpy import *
from scipy import stats
from scipy.integrate import ode
import pylab as pl
import matplotlib as mpl
import random
import math
#import scipy.io
import networkx as nx
from collections import Counter
#from scipy.signal import find_peaks

A= np.loadtxt('adjacency_erdos_renyi(500,0.01).txt',float)
G=nx.DiGraph(A)
N = 500
row_sums= A.sum(axis=1)
A1=A/row_sums[:,np.newaxis]
degree=np.sum(A,axis=0)

I=np.full(N,10.)
for i in range(350,500):
    I[i]=3.0
    
def f(t,y,p):
    a=p[0]; b=p[1]; c1=p[2]; d1=p[3]; g=p[4]
    dy = np.zeros(2*N,float)
    S1=np.dot(A1,y[1:2*N+1:2])
    for i in range(0,2*N,2):
    	dy[i] = a*(b*y[i+1]-y[i])
    	dy[i+1] = c1*(y[i+1]**2)+ d1*y[i+1]+ g -y[i]+I[int(i/2)]+K*(S1[int(i/2)]-y[i+1])
    return dy
#choose coupling strength for create data     
K=.3

y0=np.zeros(2*N,float)
ep = np.random.uniform(0.1,0.9,N)
for i in range(0,N):
    y0[2*i]= -13+ep[i]
    y0[2*i+1]= -65+ep[i]

tf = 2000
t = np.linspace(0, tf, tf*100)
y = np.zeros((len(t), len(y0)))
y[0, :] = y0
T=[] 
p = [0.02,0.2,0.04,5.0,140.0]
r = ode(f).set_integrator('vode',method='Adams')
r.set_f_params(p).set_initial_value(y0,t[0])


for j in range(1, len(t)):
    print(j)
    y[j, :] = r.integrate(t[j])   
    y0 = y[j,:]
    for k in range(0,N):
    	if y0[2*k+1] >= 30:
    		y0[2*k] = y0[2*k]+8
    		y0[2*k+1] = -65
    r.set_initial_value(y0, t[j])
    
################finding nearest neighbours
distance=[]
for j in range(0,N):
    distance.append([])
    for i in range(0,N):
        distance[j].append(nx.shortest_path_length(G,j,i))

ne=[]
for i in range(350,N):
    ne.append([])
    for j in range(0,N):
        ne[i-350].append(distance[i][j])
        
neigh=ne
for i in range(0,len(ne)):
    for j in range(0,N):
        if ne[i][j]==1.0:
            neigh[i][j]=1.
        elif ne[i][j]==2.0:
            neigh[i][j]=2.
        else:
            neigh[i][j]=0.
#finding nearest neighbours            
x1=[]
for j in range(0,len(ne)):
    x1.append([])
    for i in range(0,N):
        if ne[j][i]==1.0:
            x1[j].append(i)

#finding nearest neighbour exitable nodes            
x3=[]
for i in range(0,len(x1)):
    x3.append([])
    for j in range(0,len(x1[i])):
        if x1[i][j]<350:
            x3[i].append(x1[i][j])
#calculate r            
x4=[]
for i in range(0,len(x1)):
    x4.append([])
    for j in range(0,1):
        x4[i].append(float(len(x3[i][:]))/float(len(x1[i][:])))
        
###############calculate ISI of all quiesent nodes###############
i1=[]
for i in range(350,500):
    i1.append([])
    for j in range(50000,199999):
        if abs(y[j+1,2*i+1]-y[j,2*i+1])>65:
            i1[i-350].append(t[j+1])

is1=[]
isi=[]
for j in range(0,len(i1)):
    is1.append([])
    isi.append([])
    for i in range(0,len(i1[j])-1):
        is1[j].append(i1[j][i+1]-i1[j][i])
    isi[j].append(np.mean(is1[j]))

#data of r and ISI stored for all quisent nodes            
data=[]
for i in range(0,150):
    data.append([])
    for j in range(0,1):
        data[i].append([x4[i],isi[i]])
data_f=[]
for i in range(0,150):
    data_f.append([])
    for j in range(0,2):
        data_f[i].append(data[i][0][j][0])

#save data
np.savetxt('erdos_renyi_fraction_(k0.3).txt',data_f)
#using data we can draw r vs mean ISI
#linear average of the data
b=np.loadtxt('erdos_renyi_fraction_(k0.3).txt',float)

b11=[]
b12=[]
for i in range(0,150):
    b11.append(b[i,0])
    b12.append(b[i,1])
   
data_2=np.linspace(min(b11)-0.02,max(b11)+0.02,8)
L3= len(data_2)
r=[]
r_av=[]
for j in range(1,L3):
    for i in range(0,150):
        if data_2[j-1] < b[i][0] < data_2[j]:
            r.append(b[i][1])
    r_av.append(np.mean(r))
    del r[:]    
   
d0=[]
for j in range(1,L3):
    d0.append((data_2[j-1]+data_2[j])/2)        
#plot the average data
pl.figure(1)
pl.plot(d0,r_av,'r')

    
    
    
    