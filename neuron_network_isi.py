# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 21:56:10 2019

@author: Subrata
"""
#this code is for network_isi
import numpy as np
from numpy import *
from scipy import stats
from scipy.integrate import ode
import pylab as pl
#import matplotlib.pyplot as plt
import random
import math
#import scipy.io
import networkx as nx
from collections import Counter
import scipy.signal
from scipy.signal import find_peaks_cwt

A=np.loadtxt('adjacency_erdos_renyi(500,0.01).txt',float)
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
    
    
inter_spike1=[]
inter_spike2=[]
std1=[]
std2=[]
#create array of different coupling strength(0-2)
L=np.arange(0.3,2,.01)
for q in range(0,len(L)):
    K=L[q]
    y0 = np.zeros(2*N,float)
    ep = np.random.uniform(0.1,0.9,N)
    for i in range(0,N):
        y0[2*i] = -13+ep[i]
        y0[2*i+1] = -65.0+ep[i]
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
    # we take a threshold absolute difference between two consecutive points
    # must be greater than 65,then we call it a spike   
    i1=[]
    for i in range(0,350):
        i1.append([])
        for j in range(20000,199999):
            if abs(y[j+1,2*i+1]-y[j,2*i+1])>65:
                i1[i].append(t[j+1])
    is1=[]
    isi1=[] 
    for j in range(0,len(i1)):
        is1.append([])
        isi1.append([])
        for i in range(0,len(i1[j])-1):
            is1[j].append(i1[j][i+1]-i1[j][i])
        isi1[j].append(np.mean(is1[j]))
    inter_spike1.append(np.mean(isi1)) 
    std1.append(np.std(isi1))
    
    i2=[]
    for i in range(350,500):
        i2.append([])
        for j in range(20000,199999):
            if abs(y[j+1,2*i+1]-y[j,2*i+1])>65:
                i2[i-350].append(t[j+1])
    is2=[]
    isi2=[] 
    for j in range(0,len(i2)):
        is2.append([])
        isi2.append([])
        for i in range(0,len(i2[j])-1):
            is2[j].append(i2[j][i+1]-i2[j][i])
        isi2[j].append(np.mean(is2[j]))
    inter_spike2.append(np.mean(isi2)) 
    std2.append(np.std(isi2))  

       
#inter_spike1 will stored the mean isi of all oscillatory nodes with coupling strenth(0-2)
#inter_spike1 will stored the mean isi of all quiscent nodes with coupling strenth(0-2)
#figure of mean ISI of all exitable nodes with coupling strength
pl.figure(1)
pl.plot(L,inter_spike1,'b.')
#figure of mean ISI of all exitable nodes with coupling strength
pl.figure(2)  
pl.plot(L,inter_spike2,'r.')     

        
             
              
          
                   
