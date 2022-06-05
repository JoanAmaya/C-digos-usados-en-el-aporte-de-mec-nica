# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 02:29:49 2022

@author: soyyo
"""

import numpy as np
import matplotlib.pyplot as plt


##Constantes usadas##

R=0.0375
g=9.807
p=1.5
tetha=0.261799

#Ecuacion diferencial de cambio de variable#
def f_0(G,b):
    
    return G

#Ecuación diferencial de movimiento
def f_1(G,b):
    nume=(1+p)*g*R*np.sin(tetha)-((R**2)*np.sin(b)*(G**2))-(g*R*np.sin(b+tetha))
    den=2*(R**2)*(p+((1-np.cos(b))))
    fnc=nume/den
    return fnc
    



def Rungekuta(f_0,f_1,r0,t):
    
    h = (t[-1] - t[0])/(len(t)-1)
    
    db = np.zeros(len(t))
    dmomento_b = np.zeros(len(t))
    
    db[0] = r0[0]
    dmomento_b[0] = r0[1]
    
    
    K1 = np.zeros(2)
    K2 = np.zeros(2)
    K3 = np.zeros(2)
    K4 = np.zeros(2)
    
    
    for i in range(1,len(t)):
        
        K1[0] = f_0(dmomento_b[i-1],db[i-1])
        K1[1] = f_1(dmomento_b[i-1],db[i-1])
      
        
        K2[0] = f_0(dmomento_b[i-1]+(h/2)*K1[1],db[i-1]+(h/2)*K1[0])
        K2[1] = f_1(dmomento_b[i-1]+(h/2)*K1[1],db[i-1]+(h/2)*K1[0])
       
        
        K3[0] = f_0(dmomento_b[i-1]+(h/2)*K2[1],db[i-1]+(h/2)*K2[0])
        K3[1] = f_1(dmomento_b[i-1]+(h/2)*K2[1],db[i-1]+(h/2)*K2[0])
       
        
        K4[0] = f_0(dmomento_b[i-1]+(h)*K3[1],db[i-1]+(h)*K3[0])
        K4[1] = f_1(dmomento_b[i-1]+(h)*K3[1],db[i-1]+(h)*K3[0])
        
        
        db[i] = db[i-1] + (h/6)*(K1[0]+2*K2[0]+2*K3[0]+K4[0])
        
        dmomento_b[i] = dmomento_b[i-1] + (h/6)*(K1[1]+2*K2[1]+2*K3[1]+K4[1])
        
    return db,dmomento_b

p0=0
b0=0.523599
t=np.arange(0,2,0.0001)
r0_=np.array([b0,p0])
db,dmomento_b=Rungekuta(f_0,f_1,r0_,t)
l=0.2
x=(l-(R*db))*np.cos(tetha)+2*R*np.sin(db/2)*np.cos(db/2)*np.cos(tetha)
y=(l-(R*db))*np.sin(tetha)+2*R*np.sin(db/2)*np.sin((db/2)+tetha)
plt.plot(t,db)
