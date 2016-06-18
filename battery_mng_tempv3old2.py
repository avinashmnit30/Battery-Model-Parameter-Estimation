# -*- coding: utf-8 -*-
"""
Created on Tue May 31 17:52:39 2016

@author: Avinash
"""

import numpy
from math import *
from numpy import *

def BM_data_fun(state=1):
    ''' 
    
    1. Choose State=0 for charging and 1 for discharging     
    '''    
    
    global Vl,Vh
    global SOCl,SOCh
    global Tl,Th
    global I,Qr,V,tc_total,Cr,data_text,n_samples
    Tl=[10,0,-10]
    n_samples=[168,163,172]
    Th=25
    if state==0:
        nstr='chr'
    elif state==1:
        nstr='dis'
    data_text=''
    text1='C:/Users/Avinash/Desktop/temp_data/temp_data v2/batterydata/'
    datal1=numpy.genfromtxt(text1+nstr+str(10)+'.csv', delimiter=',')
    datal2=numpy.genfromtxt(text1+nstr+str(0)+'.csv', delimiter=',')
    datal3=numpy.genfromtxt(text1+nstr+str(-10)+'.csv', delimiter=',')
    datah=numpy.genfromtxt(text1+nstr+str(25)+'.csv', delimiter=',')
    
    I=1.75 # In Amps (Current)
    Qr=3.35 # In Ah (Rated Capacity)
    V=3.6 # in Volts (Normal Voltage)
    tc_total=1  # In hr    

    SOCl=[]    
    SOCl1=datal1.T[0]
    SOCl2=datal2.T[0]
    SOCl3=datal3.T[0]
    for i in range(len(SOCl1)):
        SOCl.append(SOCl1[i])
    for i in range(len(SOCl2)):
        SOCl.append(SOCl2[i])
    for i in range(len(SOCl3)):
        SOCl.append(SOCl3[i])
    SOCl=numpy.array(SOCl)
    
    SOCh=datah.T[0]
    
    Vl=[]    
    Vl1=datal1.T[1]
    Vl2=datal2.T[1]
    Vl3=datal3.T[1]
    for i in range(len(Vl1)):
        Vl.append(Vl1[i])
    for i in range(len(Vl2)):
        Vl.append(Vl2[i])
    for i in range(len(Vl3)):
        Vl.append(Vl3[i])
    Vl=numpy.array(Vl)
    Vh=datah.T[1]
    Cr=1
    #print Vh
    #print SOCh

def battery_fun_low(al,D,lb,ub,data_extract=0,state=1):
    N=shape(SOCh)[0]

    V_c=numpy.zeros(shape=(1,N))[0]
    Rsl=numpy.zeros(shape=(1,N))[0]
    Rll=numpy.zeros(shape=(1,N))[0]
    Rol=numpy.zeros(shape=(1,N))[0]    
    Csl=numpy.zeros(shape=(1,N))[0]
    Cll=numpy.zeros(shape=(1,N))[0]    
    V0=numpy.zeros(shape=(1,N))[0]
    SOC=numpy.zeros(shape=(1,N))[0]
    V1=numpy.zeros(shape=(1,N))[0]
    V2=numpy.zeros(shape=(1,N))[0] 
    
    fitness=0
    for i in range(N):
        if state==0:
            SOC[i]=SOCh[i]
        elif state==1:
            SOC[i]=1-SOCh[i]
    #try:
    for i in range(N):    
        if state==0:
            tc=(SOC[i]/SOC[N-1])*tc_total
        elif state==1:
            tc=(SOC[N-1]/SOC[i])*tc_total
        Rol[i]=((al[0]+al[1]*Cr+al[2]*Cr*Cr)*exp(-al[3]*SOC[i])+(al[4]+al[5]*Cr+al[6]*Cr*Cr))
        Rsl[i]=((al[7]+al[8]*Cr+al[9]*Cr*Cr)*exp(-al[10]*SOC[i])+(al[11]+al[12]*Cr+al[13]*Cr*Cr))
        Rll[i]=((al[14]+al[15]*Cr+al[16]*Cr*Cr)*exp(-al[17]*SOC[i])+(al[18]+al[19]*Cr+al[20]*Cr*Cr))
        Csl[i]=(-(al[21]+al[22]*Cr+al[23]*Cr*Cr)*exp(-al[24]*SOC[i])+(al[25]+al[26]*Cr+al[27]*Cr*Cr))
        Cll[i]=(-(al[28]+al[29]*Cr+al[30]*Cr*Cr)*exp(-al[31]*SOC[i])+(al[32]+al[33]*Cr+al[34]*Cr*Cr))
        V0[i]=((al[35]+al[36]*Cr+al[37]*Cr*Cr)*exp(-al[38]*SOC[i])+(al[39]+al[40]*SOC[i]+al[41]*SOC[i]*SOC[i]+al[42]*(pow(SOC[i],3)))-al[43]*Cr+al[44]*Cr*Cr)
        V1[i]=(Vh[i]-V0[i]-I*Rol[i])/(Csl[i]*(1/Csl[i]+1/Cll[i]))
        V2[i]=(Vh[i]-V0[i]-I*Rol[i])/(Cll[i]*(1/Csl[i]+1/Cll[i]))
        V_c[i]=V0[i]+I*Rol[i]+I*Rsl[i]*(1-exp(-tc/(Rsl[i]*Csl[i])))+I*Rll[i]*(1-exp(-tc/(Rll[i]*Cll[i])))+V1[i]*(exp(-tc/(Rsl[i]*Csl[i])))+V2[i]*(exp(-tc/(Rll[i]*Cll[i])))
        fitness+=abs(V_c[i]-Vh[i])
    fitness=fitness
    #except(OverflowError):
     #   return(10000) 
    if data_extract==1:
        numpy.savetxt("Vc"+data_text+".csv",numpy.array(V_c),delimiter=",")
        numpy.savetxt("Rol"+data_text+".csv",numpy.array(Rol),delimiter=",")
        numpy.savetxt("Rsl"+data_text+".csv",numpy.array(Rsl),delimiter=",")
        numpy.savetxt("Csl"+data_text+".csv",numpy.array(Csl),delimiter=",")
        numpy.savetxt("Rll"+data_text+".csv",numpy.array(Rll),delimiter=",")
        numpy.savetxt("Cll"+data_text+".csv",numpy.array(Cll),delimiter=",")
        numpy.savetxt("V"+data_text+".csv",numpy.array(V_c),delimiter=",")
        #numpy.savetxt("V0"+data_text+".csv",numpy.array(V0),delimiter=",")       
    for i in range(D):
        if al[i]>ub:
            fitness+=10*(al[i]-ub)**2              
        if al[i]<lb:
            fitness+=10*(al[i]-lb)**2          
    return(fitness)
    
def battery_fun_high(a,al,D,lb,ub,data_extract=0,state=1):  
    fitness=0    

    N=shape(SOCl)[0]

    Roh=numpy.zeros(shape=(1,N))[0]
    Rsh=numpy.zeros(shape=(1,N))[0]
    Csh=numpy.zeros(shape=(1,N))[0]    
    Rlh=numpy.zeros(shape=(1,N))[0]
    Clh=numpy.zeros(shape=(1,N))[0]    
    
    V_c=numpy.zeros(shape=(1,N))[0]
    Rsl=numpy.zeros(shape=(1,N))[0]
    Rll=numpy.zeros(shape=(1,N))[0]
    Rol=numpy.zeros(shape=(1,N))[0]    
    Csl=numpy.zeros(shape=(1,N))[0]
    Cll=numpy.zeros(shape=(1,N))[0]    
    V0=numpy.zeros(shape=(1,N))[0]
    SOC=numpy.zeros(shape=(1,N))[0]
    V1=numpy.zeros(shape=(1,N))[0]
    V2=numpy.zeros(shape=(1,N))[0] 

    for i in range(N):
        if state==0:
            SOC[i]=SOCl[i]
        elif state==1:
            SOC[i]=1-SOCl[i]
    n_temps=len(Tl)
    
    for k in range(n_temps):
        add=0
        i=0
        while i<n_samples[k]:
            if state==0:
                tc=(SOC[add+i]/SOC[add+n_samples[k]-1])*tc_total
            elif state==1:
                tc=(SOC[add+n_samples[k]-1]/SOC[add+i])*tc_total
            Rol[i]=((al[0]+al[1]*Cr+al[2]*Cr*Cr)*exp(-al[3]*SOC[i])+(al[4]+al[5]*Cr+al[6]*Cr*Cr))
            Rsl[i]=((al[7]+al[8]*Cr+al[9]*Cr*Cr)*exp(-al[10]*SOC[i])+(al[11]+al[12]*Cr+al[13]*Cr*Cr))
            Rll[i]=((al[14]+al[15]*Cr+al[16]*Cr*Cr)*exp(-al[17]*SOC[i])+(al[18]+al[19]*Cr+al[20]*Cr*Cr))
            Csl[i]=(-(al[21]+al[22]*Cr+al[23]*Cr*Cr)*exp(-al[24]*SOC[i])+(al[25]+al[26]*Cr+al[27]*Cr*Cr))
            Cll[i]=(-(al[28]+al[29]*Cr+al[30]*Cr*Cr)*exp(-al[31]*SOC[i])+(al[32]+al[33]*Cr+al[34]*Cr*Cr))
            dT=Th-Tl[k]
            Roh[i]=Rol[i]*a[0]*exp(a[1]/(Tl[k]-a[2]))
            Rsh[i]=Rsl[i]+a[3]*dT+a[4]*dT*SOC[add+i]
            Rlh[i]=Rll[i]+(a[5]*dT)*exp(-al[15]*SOC[add+i])+a[6]*dT
            Csh[i]=Csl[i]+a[7]*dT*SOC[add+i]+a[8]*dT
            Clh[i]=Cll[i]*a[9]*exp(a[10]/Tl[k])
#            elif state==1:
#                Roh[i]=Rol[i]*a[0]*exp(a[1]/(Tl[k]-a[2]))
#                Rsh[i]=a[3]*exp(a[4]*Tl[k]-(a[5]*Tl[k]+a[6])*SOC[add+i])+a[7]*Tl[k]+a[8]
#                Rlh[i]=Rll[i]*(a[9]*Tl[k])
#                Csh[i]=Csl[i]+a[7]*dT*SOC[add+i]+a[8]*dT
#                Clh[i]=Cll[i]*a[9]*exp(a[10]/Tl[k])
            V0[i]=((al[35]+al[36]*Cr+al[37]*Cr*Cr)*exp(-al[38]*SOC[i])+(al[39]+al[40]*SOC[i]+al[41]*SOC[i]*SOC[i]+al[42]*(pow(SOC[i],3)))-al[43]*Cr+al[44]*Cr*Cr)
            V1[i]=(Vl[i]-V0[i]-I*Roh[i])/(Csh[i]*(1/Csh[i]+1/Clh[i]))
            V2[i]=(Vl[i]-V0[i]-I*Roh[i])/(Clh[i]*(1/Csh[i]+1/Clh[i]))
            V_c[i]=V0[i]+I*Roh[i]+I*Rsh[i]*(1-exp(-tc/(Rsh[i]*Csh[i])))+I*Rlh[i]*(1-exp(-tc/(Rll[i]*Cll[i])))+V1[i]*(exp(-tc/(Rsl[i]*Csl[i])))+V2[i]*(exp(-tc/(Rll[i]*Cll[i])))        
            temp=str(abs(V_c[i]-Vl[i]))
            if temp=='nan':
                fitness+=10
            else:
                fitness+=abs(V_c[i]-Vl[i])
            i+=3
        add+=n_samples[k]
        

            
                    
    #V_temp=V_c-Vl
    #print V_temp

        #Tsh=Rsh[i]*Csh[i]
        #Tlh=Rlh[i]*Clh[i]
        #Vsh=Rsh[i]*I
        #Vlh=Rlh[i]*I
        #Vth=Vsh*(1-exp(-t/Tsh))+Vlh*(1-exp(-t/Tlh))+I*Roh
        #fitness+=abs(Vh[i]-Vth)
    #print V_c
    if data_extract==1:
        numpy.savetxt("Vc"+data_text+'2'+".csv",numpy.array(V_c),delimiter=",")
        numpy.savetxt("Roh"+data_text+'2'+".csv",numpy.array(Roh),delimiter=",")
        numpy.savetxt("Rsh"+data_text+'2'+".csv",numpy.array(Rsh),delimiter=",")
        numpy.savetxt("Csh"+data_text+'2'+".csv",numpy.array(Csh),delimiter=",")
        numpy.savetxt("Rlh"+data_text+'2'+".csv",numpy.array(Rlh),delimiter=",")
        numpy.savetxt("Clh"+data_text+'2'+".csv",numpy.array(Clh),delimiter=",")
        numpy.savetxt("V"+data_text+'2'+".csv",numpy.array(V_c),delimiter=",")
    for i in range(D):
        if a[i]>ub:
            fitness+=10*(a[i]-ub)**2              
        if a[i]<lb:
            fitness+=10*(a[i]-lb)**2          
    return(fitness)