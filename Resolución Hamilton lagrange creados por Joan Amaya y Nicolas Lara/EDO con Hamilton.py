# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 00:44:46 2022

@author: soyyo
"""
import numpy as np
import matplotlib.pyplot as plt


##Constantes usadas##
m0=0.075
M=0.15
R=0.02
g=9.807
tetha=0.261799
#Ecuacion diferencial de la posición#
def f_0(P_b,b):
    sn=np.sin(b/2)
    den=(2*M*(R**2)+(4*m0*(R**2)*(sn**2)))
    fnc=P_b/den
    
    return fnc

#Ecuación diferencial del momento
def f_1(P_b,b):
    sn1=np.sin(b/2)
    den1=(2*M*(R**2)+(4*m0*(R**2)*(sn1**2)))
    pr1=((P_b**2)*(R**2)*(m0*np.sin(b)))/(den1**2)
    pr2=M*g*R*np.sin(tetha)
    pr3=m0*g*R*np.sin(b+tetha)
    pr4=m0*g*R*np.sin(tetha)
    fnc=pr1+pr2-pr3+pr4
    return fnc
    

#Metodo Runge-Kutta de 4 orden

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

#condiciones iniciales 
p0=0
b0=1.04467329258475
t=np.arange(0,2,0.001)
r0_=np.array([b0,p0])
#Resolucion del sistema

db,dmomento_b=Rungekuta(f_0,f_1,r0_,t)
l=0.2
#Coordenadas cartesianas
lb=list(db)
mx=-1000000000000
n=0
m=0
mn=1000000000000
for i in range(0,len(db)):
    if db[i]>mx:
        mx=db[i]
        n=i
    if db[i]>db[i+1]:    
        break

for i in range(n,len(db)):
    if db[i]<mn:
        mn=db[i]
        m=i
    if db[i]<db[i+1]:    
        break

print(m) 
periodo=(t[m]-t[n])*2


x=((R*db))*np.cos(tetha)+2*R*np.sin(db/2)*np.cos(db/2)*np.cos(tetha)
y=((R*db))*np.sin(tetha)+2*R*np.sin(db/2)*np.sin((db/2)+tetha)
plt.title(label="B vs T"+"         Periodo de " +str(periodo)+"s")
plt.ylabel("B(rad)")
plt.xlabel("t(s)")
plt.plot(t,db, label="B_0 = "+ str(75))




