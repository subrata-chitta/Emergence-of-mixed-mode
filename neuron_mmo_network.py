# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 12:42:43 2019

@author: subrata
"""

#this code is for figure_1
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
from scipy.signal import find_peaks
#load the adjacency matrix
A=np.loadtxt('adjacency_erdos_renyi(500,0.01).txt',float)
N = 500
row_sums= A.sum(axis=1)
A1=A/row_sums[:,np.newaxis]
degree=np.sum(A,axis=0)
#create array I
I=np.full(N,10.)
for i in range(350,500):
    I[i]=3.0

#define the network function    
def f(t,y,p):
    a=p[0]; b=p[1]; c1=p[2]; d1=p[3]; g=p[4]
    dy = np.zeros(2*N,float)
    S1=np.dot(A1,y[1:2*N+1:2])
    for i in range(0,2*N,2):
    	dy[i] = a*(b*y[i+1]-y[i])
    	dy[i+1] = c1*(y[i+1]**2)+ d1*y[i+1]+ g -y[i]+I[int(i/2)]+K*(S1[int(i/2)]-y[i+1])
    return dy
#K is the coupling strength    
K=0.3
#initial condition
y0=np.zeros(2*N,float)
ep = np.random.uniform(0.1,0.9,N)
for i in range(0,N):
    y0[2*i]= -13+ep[i]
    y0[2*i+1]= -65+ep[i]

tf = 1000
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
    
#plotting figures
    
ind = (0,50000,100000)
ind2=(-80,-40,0,40)
#plotting time series of three different quiesent nodes 
pl.figure(1)
pl.plot(y[0:100000,731],'r')
pl.plot(y[0:100000,803],'b')
pl.plot(y[0:100000,997],'g')
pl.xlabel('$t$',fontsize=35)
pl.ylabel(r'$v$',fontsize=35)
pl.xticks(ind,('0','50','100'),fontsize=35)
pl.yticks(ind2,('-80','-40','0','40'),fontsize=35)
pl.yticks(fontsize=35)
pl.tight_layout()
#plotting time series of three different excitable nodes
pl.figure(2)
pl.plot(y[0:100000,101],'b')
pl.plot(y[0:100000,391],'r')
pl.plot(y[0:100000,393],'g')
pl.xlabel('$t$',fontsize=35)
pl.ylabel(r'$v$',fontsize=35)
pl.xticks(ind,('0','50','100'),fontsize=35)
pl.yticks(ind2,('-80','-40','0','40'),fontsize=35)
pl.yticks(fontsize=35)  
pl.tight_layout()  
#plotting the spatio-temporal figure
pl.figure(3)
pl.imshow(y[0:10000,1:1001:2],origin='lower',aspect='auto',cmap='cool')
pl.xlabel('$N$',fontsize=35)
pl.ylabel(r'$time$',fontsize=35)
pl.xticks(ind3,('0','250','500'),fontsize=35)
pl.yticks(ind,('0','50','100'),fontsize=35)
pl.colorbar()
pl.tight_layout()

