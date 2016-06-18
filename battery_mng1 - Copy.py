# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 11:29:07 2015

@author: Avinash
"""
import numpy
from numpy import *
from math import *
from numba import double
from numba.decorators import jit, autojit

def BM_data_fun(state,CDr):
    ''' 
    
    1. Choose State=0 for charging and 1 for discharging
    2. CDr is used for both both Cr (charging rate) and Dr (discharging rate) depeding on state
     
    '''

    
    if state==0:
        V_m=numpy.genfromtxt('Cr'+str(CDr)+'.csv', delimiter=',')
    elif state==1:
        V_m=numpy.genfromtxt('Dr'+str(CDr)+'.csv', delimiter=',')        
    return(V_m)
    
#@autojit
def BM_fun1_arr(a,D,battery_no,state,Cr,V_m,SOC_m,lb=0,ub=50,data_extract=0,data_text=" "):
    '''    
    
    1. Choose battery_no=1 for EIG, 2 for sony_us18650, 3 for panasonic, 4 for sanyo
    2. Choose State=0 for charging and 1 for discharging    
    
    '''
    if battery_no==1:
        V=2.5 # Voltage
        Qr=8 # Capacity
        SOC_max=1 # For discharge it acts as DOD_max
        C=8 # Ampere
        if Cr==1:
            tc_total=1  # In hr
        elif Cr==4:
            tc_total=0.333
        
    Ic=C*Cr
    N=shape(SOC_m)[0]

    V_c=numpy.zeros(shape=(1,N))[0]
    R1=numpy.zeros(shape=(1,N))[0]
    R2=numpy.zeros(shape=(1,N))[0]
    C=numpy.zeros(shape=(1,N))[0]
    V0=numpy.zeros(shape=(1,N))[0]
    SOC=numpy.zeros(shape=(1,N))[0]

    fitness=0
    for i in range(N):
        if state==0:
            SOC[i]=SOC_m[i]
        elif state==1:
            SOC[i]=1-SOC_m[i]
#    try:
    for i in range(N):    
        if state==0:
            tc=(SOC[i]/SOC[N-1])*tc_total
        elif state==1:
            tc=(SOC[N-1]/SOC[i])*tc_total
        R1[i]=((a[0]+a[1]*Cr+a[2]*Cr*Cr)*exp(-a[3]*SOC[i])+(a[4]+a[5]*Cr+a[6]*Cr*Cr))
        R2[i]=((a[7]+a[8]*Cr+a[9]*Cr*Cr)*exp(-a[10]*SOC[i])+(a[11]+a[12]*Cr+a[13]*Cr*Cr))
        C[i]=(-(a[14]+a[15]*Cr+a[16]*Cr*Cr)*exp(-a[17]*SOC[i])+(a[18]+a[19]*Cr+a[20]*Cr*Cr))
        V0[i]=((a[21]+a[22]*Cr+a[23]*Cr*Cr)*exp(-a[24]*SOC[i])+(a[25]+a[26]*SOC[i]+a[27]*SOC[i]*SOC[i]+a[28]*(pow(SOC[i],3)))-a[29]*Cr+a[30]*Cr*Cr)
        V_c[i]=(((Qr/C[i])+Ic*R2[i])*exp(-tc/(R2[i]*C[i]))+V0[i]-Ic*(R1[i]+R2[i]))
        fitness+=abs(V_c[i]-V_m[i])
#    except(OverflowError):
#        return(10000) 
    if data_extract==1:
        numpy.savetxt("Vc"+data_text+".csv",numpy.array(V_c),delimiter=",")
        numpy.savetxt("R1"+data_text+".csv",numpy.array(R1),delimiter=",")
        numpy.savetxt("R2"+data_text+".csv",numpy.array(R2),delimiter=",")
        numpy.savetxt("C"+data_text+".csv",numpy.array(C),delimiter=",")
        numpy.savetxt("V0"+data_text+".csv",numpy.array(V0),delimiter=",")       
    for i in range(D):
        if a[i]>ub:
            fitness+=100*(a[i]-ub)**2              
        if a[i]<lb:
            fitness+=100*(a[i]-lb)**2          
#    print fitness
    return(fitness)

#@autojit
def BM_fun1_list(a,D,battery_no,state,Cr,V_m,SOC_m,lb=0,ub=50,data_extract=0,data_text=" "):
    '''    
    
    1. Choose battery_no=1 for EIG, 2 for sony_us18650, 3 for panasonic, 4 for sanyo
    2. Choose State=0 for charging and 1 for discharging    
    
    '''
    if battery_no==1:
        V=2.5 # Voltage
        Qr=8 # Capacity
        SOC_max=1 # For discharge it acts as DOD_max
        C=8 # Ampere
        if Cr==1:
            tc_total=1  # In hr
        elif Cr==4:
            tc_total=0.333
        
    Ic=C*Cr
    N=shape(SOC_m)[0]

    V_c=[]
    R1=[]
    R2=[]
    C=[]
    V0=[]
    SOC=[]

    fitness=0
    for i in range(N):
        if state==0:
            SOC.append(SOC_m[i])
        elif state==1:
            SOC.append(1-SOC_m[i])
    try:
        for i in range(N):    
            if state==0:
                tc=(SOC[i]/SOC[N-1])*tc_total
            elif state==1:
                tc=(SOC[N-1]/SOC[i])*tc_total
            R1.append(((a[0]+a[1]*Cr+a[2]*Cr*Cr)*exp(-a[3]*SOC[i])+(a[4]+a[5]*Cr+a[6]*Cr*Cr)))
            R2.append(((a[7]+a[8]*Cr+a[9]*Cr*Cr)*exp(-a[10]*SOC[i])+(a[11]+a[12]*Cr+a[13]*Cr*Cr)))
            C.append((-(a[14]+a[15]*Cr+a[16]*Cr*Cr)*exp(-a[17]*SOC[i])+(a[18]+a[19]*Cr+a[20]*Cr*Cr)))
            V0.append(((a[21]+a[22]*Cr+a[23]*Cr*Cr)*exp(-a[24]*SOC[i])+(a[25]+a[26]*SOC[i]+a[27]*SOC[i]*SOC[i]+a[28]*(pow(SOC[i],3)))-a[29]*Cr+a[30]*Cr*Cr))
            V_c.append((((Qr/C[i])+Ic*R2[i])*exp(-tc/(R2[i]*C[i]))+V0[i]-Ic*(R1[i]+R2[i])))
            fitness+=abs(V_c[i]-V_m[i])
    except(OverflowError):
        return(10000) 
    if data_extract==1:
        numpy.savetxt("Vc"+data_text+".csv",numpy.array(V_c),delimiter=",")
        numpy.savetxt("R1"+data_text+".csv",numpy.array(R1),delimiter=",")
        numpy.savetxt("R2"+data_text+".csv",numpy.array(R2),delimiter=",")
        numpy.savetxt("C"+data_text+".csv",numpy.array(C),delimiter=",")
        numpy.savetxt("V0"+data_text+".csv",numpy.array(V0),delimiter=",")       
    for i in range(D):
        if a[i]>ub:
            fitness+=100*(a[i]-ub)**2              
        if a[i]<lb:
            fitness+=100*(a[i]-lb)**2          
    return(fitness/N)
                


    
def BM_fun1_fast(a,D,battery_no,state,Cr,V_m,SOC_m,lb=0,ub=50,data_extract=0,data_text=" "):
    '''    
    
    1. Choose battery_no=1 for EIG, 2 for sony_us18650, 3 for panasonic, 4 for sanyo
    2. Choose State=0 for charging and 1 for discharging    
    
    '''
    if battery_no==1:
        V=2.5 # Voltage
        Qr=8 # Capacity
        SOC_max=1 # For discharge it acts as DOD_max
        C=8 # Ampere
        if Cr==1:
            tc_total=1  # In hr
        elif Cr==4:
            tc_total=0.333
        
    Ic=C*Cr
    N=shape(SOC_m)[0]

    V_c=numpy.zeros(shape=(1,N))[0]
    R1=numpy.zeros(shape=(1,N))[0]
    R2=numpy.zeros(shape=(1,N))[0]
    C=numpy.zeros(shape=(1,N))[0]
    V0=numpy.zeros(shape=(1,N))[0]
    SOC=numpy.zeros(shape=(1,N))[0]

    fitness=0
    for i in range(N):
        if state==0:
            SOC[i]=SOC_m[i]
        elif state==1:
            SOC[i]=1-SOC_m[i]
    for i in range(N):    
        tc=(SOC[i]/SOC[N-1])*tc_total
        R1[i]=((a[0]+a[1]*Cr+a[2]*Cr*Cr)*exp(-a[3]*SOC[i])+(a[4]+a[5]*Cr+a[6]*Cr*Cr))
        R2[i]=((a[7]+a[8]*Cr+a[9]*Cr*Cr)*exp(-a[10]*SOC[i])+(a[11]+a[12]*Cr+a[13]*Cr*Cr))
        C[i]=(-(a[14]+a[15]*Cr+a[16]*Cr*Cr)*exp(-a[17]*SOC[i])+(a[18]+a[19]*Cr+a[20]*Cr*Cr))
        V0[i]=((a[21]+a[22]*Cr+a[23]*Cr*Cr)*exp(-a[24]*SOC[i])+(a[25]+a[26]*SOC[i]+a[27]*SOC[i]*SOC[i]+a[28]*(pow(SOC[i],3)))-a[29]*Cr+a[30]*Cr*Cr)
        V_c[i]=(((Qr/C[i])+Ic*R2[i])*exp(-tc/(R2[i]*C[i]))+V0[i]-Ic*(R1[i]+R2[i]))
        fitness+=abs(V_c[i]-V_m[i])
    if data_extract==1:
        numpy.savetxt("Vc"+data_text+".csv",numpy.array(V_c),delimiter=",")
        numpy.savetxt("R1"+data_text+".csv",numpy.array(R1),delimiter=",")
        numpy.savetxt("R2"+data_text+".csv",numpy.array(R2),delimiter=",")
        numpy.savetxt("C"+data_text+".csv",numpy.array(C),delimiter=",")
        numpy.savetxt("V0"+data_text+".csv",numpy.array(V0),delimiter=",")
    for i in range(D):
        if a[i]>ub:
            fitness=100*(a[i]-ub)**2              
        if a[i]<lb:
            fitness=100*(a[i]-lb)**2  
    return(fitness)
                
def BM_fun1_fast2(a,D,battery_no,state,Cr,V_m,SOC_m,data_extract=0,data_text=" "):
    '''    
    
    1. Choose battery_no=1 for EIG, 2 for sony_us18650, 3 for panasonic, 4 for sanyo
    2. Choose State=0 for charging and 1 for discharging    
    
    '''
    if battery_no==1:
        V=2.5 # Voltage
        Qr=8 # Capacity
        SOC_max=1 # For discharge it acts as DOD_max
        C=8 # Ampere
        if Cr==1:
            tc_total=1  # In hr
        elif Cr==4:
            tc_total=0.333
        
    Ic=C*Cr
    N=shape(SOC_m)[0]
    V_c=numpy.zeros(shape=(1,N))[0]
    SOC=numpy.zeros(shape=(1,N))[0]

    fitness=0
    for i in range(N):
        if state==0:
            SOC[i]=SOC_m[i]
        elif state==1:
            SOC[i]=1-SOC_m[i]

    for i in range(N):    
        tc=(SOC[i]/SOC[N-1])*tc_total
        R1=((a[0]+a[1]*Cr+a[2]*Cr*Cr)*exp(-a[3]*SOC[i])+(a[4]+a[5]*Cr+a[6]*Cr*Cr))
        R2=((a[7]+a[8]*Cr+a[9]*Cr*Cr)*exp(-a[10]*SOC[i])+(a[11]+a[12]*Cr+a[13]*Cr*Cr))
        C=(-(a[14]+a[15]*Cr+a[16]*Cr*Cr)*exp(-a[17]*SOC[i])+(a[18]+a[19]*Cr+a[20]*Cr*Cr))
        V0=((a[21]+a[22]*Cr+a[23]*Cr*Cr)*exp(-a[24]*SOC[i])+(a[25]+a[26]*SOC[i]+a[27]*SOC[i]*SOC[i]+a[28]*(pow(SOC[i],3)))-a[29]*Cr+a[30]*Cr*Cr)
        V_c[i]=(((Qr/C)+Ic*R2)*exp(-tc/(R2*C))+V0-Ic*(R1+R2))
        fitness+=abs(V_c[i]-V_m[i])
             
    return(fitness)


def BM_fun2(a,D,battery_no,state,Cr,V_m,SOC_m,lb=0,ub=50,data_extract=0,data_text=" "):
    '''    
    
    1. Choose battery_no=1 for EIG, 2 for sony_us18650, 3 for panasonic, 4 for sanyo
    2. Choose State=0 for charging and 1 for discharging    
    
    '''
    if battery_no==1:
        V=2.5 # Voltage
        Qr=8 # Capacity
        SOC_max=1 # For discharge it acts as DOD_max
        C=8 # Ampere
        if Cr==1:
            tc_total=1  # In hr
        elif Cr==4:
            tc_total=0.333
        
    Ic=C*Cr
    N=shape(SOC_m)[0]

    V_c=[]
    R1=[]
    R2=[]
    R=[]
    C1=[]
    C2=[]
    V0=[]
    SOC=[]
    V1=[]
    V2=[]

    fitness=0
    for i in range(N):
        if state==0:
            SOC.append(SOC_m[i])
        elif state==1:
            SOC.append(1-SOC_m[i])
    #try:
    for i in range(N):    
        tc=(SOC[i]/SOC[N-1])*tc_total
        R.append(((a[0]+a[1]*Cr+a[2]*Cr*Cr)*exp(-a[3]*SOC[i])+(a[4]+a[5]*Cr+a[6]*Cr*Cr)))
        R1.append(((a[7]+a[8]*Cr+a[9]*Cr*Cr)*exp(-a[10]*SOC[i])+(a[11]+a[12]*Cr+a[13]*Cr*Cr)))
        R2.append(((a[14]+a[15]*Cr+a[16]*Cr*Cr)*exp(-a[17]*SOC[i])+(a[18]+a[19]*Cr+a[20]*Cr*Cr)))
        C1.append((-(a[21]+a[22]*Cr+a[23]*Cr*Cr)*exp(-a[24]*SOC[i])+(a[25]+a[26]*Cr+a[27]*Cr*Cr)))
        C2.append((-(a[28]+a[29]*Cr+a[30]*Cr*Cr)*exp(-a[31]*SOC[i])+(a[32]+a[33]*Cr+a[34]*Cr*Cr)))
        V0.append(((a[35]+a[36]*Cr+a[37]*Cr*Cr)*exp(-a[38]*SOC[i])+(a[39]+a[40]*SOC[i]+a[41]*SOC[i]*SOC[i]+a[42]*(pow(SOC[i],3)))-a[43]*Cr+a[44]*Cr*Cr))
        V1.append((V_m[0]-V0[0]-Ic*R[i])/(C1[i]*(1/C1[i]+1/C2[i])))            
        V2.append((V_m[0]-V0[0]-Ic*R[i])/(C2[i]*(1/C1[i]+1/C2[i])))          
        V_c.append(V0[i]+Ic*R[i]+Ic*R1[i]*(1-exp(-tc/(R1[i]*C1[i])))+Ic*R2[i]*(1-exp(-tc/(R2[i]*C2[i])))+V1[i]*(exp(-tc/(R1[i]*C1[i])))+V2[i]*(exp(-tc/(R2[i]*C2[i]))))
        fitness+=abs(V_c[i]-V_m[i])
    #except(OverflowError):
     #   return(10000) 
    if data_extract==1:
        numpy.savetxt("Vc"+data_text+".csv",numpy.array(V_c),delimiter=",")
        numpy.savetxt("R1"+data_text+".csv",numpy.array(R1),delimiter=",")
        numpy.savetxt("R2"+data_text+".csv",numpy.array(R2),delimiter=",")
        numpy.savetxt("C"+data_text+".csv",numpy.array(C),delimiter=",")
        numpy.savetxt("V0"+data_text+".csv",numpy.array(V0),delimiter=",")       
    for i in range(D):
        if a[i]>ub:
            fitness+=100*(a[i]-ub)**2              
        if a[i]<lb:
            fitness+=100*(a[i]-lb)**2          
    return(fitness)

def BM_fun3(a,D,battery_no1,state1,Cr1,V_m1,SOC_m1,battery_no2,state2,Cr2,V_m2,SOC_m2,lb=0,ub=50,data_extract=0,data_text=" "):
    # Multi-objective variant of fun1
    '''    
    
    1. Choose battery_no=1 for EIG, 2 for sony_us18650, 3 for panasonic, 4 for sanyo
    2. Choose State=0 for charging and 1 for discharging    
    
    '''
    if battery_no==1:
        V=2.5 # Voltage
        Qr=8 # Capacity
        SOC_max=1 # For discharge it acts as DOD_max
        C=8 # Ampere
        if Cr==1:
            tc_total=1  # In hr
        elif Cr==4:
            tc_total=0.333
        
    Ic=C*Cr
    N=shape(SOC_m)[0]

    V_c=[]
    R1=[]
    R2=[]
    C=[]
    V0=[]
    SOC=[]

    fitness=0
    for i in range(N):
        if state==0:
            SOC.append(SOC_m[i])
        elif state==1:
            SOC.append(1-SOC_m[i])
    try:
        for i in range(N):    
            tc=(SOC[i]/SOC[N-1])*tc_total
            R1.append(((a[0]+a[1]*Cr+a[2]*Cr*Cr)*exp(-a[3]*SOC[i])+(a[4]+a[5]*Cr+a[6]*Cr*Cr)))
            R2.append(((a[7]+a[8]*Cr+a[9]*Cr*Cr)*exp(-a[10]*SOC[i])+(a[11]+a[12]*Cr+a[13]*Cr*Cr)))
            C.append((-(a[14]+a[15]*Cr+a[16]*Cr*Cr)*exp(-a[17]*SOC[i])+(a[18]+a[19]*Cr+a[20]*Cr*Cr)))
            V0.append(((a[21]+a[22]*Cr+a[23]*Cr*Cr)*exp(-a[24]*SOC[i])+(a[25]+a[26]*SOC[i]+a[27]*SOC[i]*SOC[i]+a[28]*(pow(SOC[i],3)))-a[29]*Cr+a[30]*Cr*Cr))
            V_c.append((((Qr/C[i])+Ic*R2[i])*exp(-tc/(R2[i]*C[i]))+V0[i]-Ic*(R1[i]+R2[i])))
            fitness+=abs(V_c[i]-V_m[i])
    except(OverflowError):
        return(10000) 
    if data_extract==1:
        numpy.savetxt("Vc"+data_text+".csv",numpy.array(V_c),delimiter=",")
        numpy.savetxt("R1"+data_text+".csv",numpy.array(R1),delimiter=",")
        numpy.savetxt("R2"+data_text+".csv",numpy.array(R2),delimiter=",")
        numpy.savetxt("C"+data_text+".csv",numpy.array(C),delimiter=",")
        numpy.savetxt("V0"+data_text+".csv",numpy.array(V0),delimiter=",")       
    for i in range(D):
        if a[i]>ub:
            fitness+=100*(a[i]-ub)**2              
        if a[i]<lb:
            fitness+=100*(a[i]-lb)**2          
    return(fitness)
                


def BM_fun4(a,D,battery_noa,statea,Cra,V_ma,SOC_ma,battery_nob,stateb,Crb,V_mb,SOC_mb,lb=0,ub=50,data_extract=0,data_text=" "):
    # Multi-objective variant of fun2    
    '''    
    
    1. Choose battery_no=1 for EIG, 2 for sony_us18650, 3 for panasonic, 4 for sanyo
    2. Choose State=0 for charging and 1 for discharging    
    
    '''
    if battery_noa==1:
        Va=2.5 # Voltage
        Qra=8 # Capacity
        SOC_maxa=1 # For discharge it acts as DOD_max
        Ca=8 # Ampere
        if Cra==1:
            tc_totala=1  # In hr
        elif Cra==4:
            tc_totala=0.333
    if battery_nob==1:
        Vb=2.5 # Voltage
        Qrb=8 # Capacity
        SOC_maxb=1 # For discharge it acts as DOD_max
        Cb=8 # Ampere
        if Crb==1:
            tc_totalb=1  # In hr
        elif Crb==4:
            tc_totalb=0.333
        
    Ica=Ca*Cra  
    Icb=Cb*Crb    
    Na=shape(SOC_ma)[0]
    Nb=shape(SOC_mb)[0]

    V_ca=[]
    R1a=[]
    R2a=[]
    Ra=[]
    C1a=[]
    C2a=[]
    V0a=[]
    SOCa=[]
    V1a=[]
    V2a=[]
    V_cb=[]
    R1b=[]
    R2b=[]
    Rb=[]
    C1b=[]
    C2b=[]
    V0b=[]
    SOCb=[]
    V1b=[]
    V2b=[]

    fitness=0
    for i in range(Na):
        if statea==0:
            SOCa.append(SOC_ma[i])
        elif statea==1:
            SOCa.append(1-SOC_ma[i])
    for i in range(Nb):
        if stateb==0:
            SOCb.append(SOC_mb[i])
        elif stateb==1:
            SOCb.append(1-SOC_mb[i])
    try:
        for i in range(Na):    
            tca=(SOCa[i]/SOCa[Na-1])*tc_totala
            R1a.append(((a[0]+a[1]*Cra+a[2]*Cra*Cra)*exp(-a[3]*SOCa[i])+(a[4]+a[5]*Cra+a[6]*Cra*Cra)))
            R2a.append(((a[7]+a[8]*Cra+a[9]*Cra*Cra)*exp(-a[10]*SOCa[i])+(a[11]+a[12]*Cra+a[13]*Cra*Cra)))
            Ra.append(((a[31]+a[32]*Cra+a[33]*Cra*Cra)*exp(-a[34]*SOCa[i])+(a[35]+a[36]*Cra+a[37]*Cra*Cra)))
            C1a.append((-(a[14]+a[15]*Cra+a[16]*Cra*Cra)*exp(-a[17]*SOCa[i])+(a[18]+a[19]*Cra+a[20]*Cra*Cra)))
            C2a.append((-(a[38]+a[39]*Cra+a[40]*Cra*Cra)*exp(-a[41]*SOCa[i])+(a[42]+a[43]*Cra+a[44]*Cra*Cra)))
            V0a.append(((a[21]+a[22]*Cra+a[23]*Cra*Cra)*exp(-a[24]*SOCa[i])+(a[25]+a[26]*SOCa[i]+a[27]*SOCa[i]*SOCa[i]+a[28]*(pow(SOCa[i],3)))-a[29]*Cra+a[30]*Cra*Cra))
            V1a.append((V_ma[0]-V0a[0]-Ica*Ra[i])/(C1a[i]*(1/C1a[i]+1/C2a[i])))            
            V2a.append((V_ma[0]-V0a[0]-Ica*Ra[i])/(C2a[i]*(1/C1a[i]+1/C2a[i])))          
            V_ca.append(V0a[i]+Ica*Ra[i]+Ica*R1a[i]*(1-exp(-tca/(R1a[i]*C1a[i])))+Ica*R2a[i]*(1-exp(-tca/(R2a[i]*C2a[i])))+V1a[i]*(exp(-tca/(R1a[i]*C1a[i])))+V2a[i]*(exp(-tca/(R2a[i]*C2a[i]))))
            fitness+=abs(V_ca[i]-V_ma[i])
        for i in range(Nb):    
            tcb=(SOCb[i]/SOCb[Nb-1])*tc_totalb
            R1b.append(((a[0]+a[1]*Crb+a[2]*Crb*Crb)*exp(-a[3]*SOCb[i])+(a[4]+a[5]*Crb+a[6]*Crb*Crb)))
            R2b.append(((a[7]+a[8]*Crb+a[9]*Crb*Crb)*exp(-a[10]*SOCb[i])+(a[11]+a[12]*Crb+a[13]*Crb*Crb)))
            Rb.append(((a[31]+a[32]*Crb+a[33]*Crb*Crb)*exp(-a[34]*SOCb[i])+(a[35]+a[36]*Crb+a[37]*Crb*Crb)))
            C1b.append((-(a[14]+a[15]*Crb+a[16]*Crb*Crb)*exp(-a[17]*SOCb[i])+(a[18]+a[19]*Crb+a[20]*Crb*Crb)))
            C2b.append((-(a[38]+a[39]*Crb+a[40]*Crb*Crb)*exp(-a[41]*SOCb[i])+(a[42]+a[43]*Crb+a[44]*Crb*Crb)))
            V0b.append(((a[21]+a[22]*Crb+a[23]*Crb*Crb)*exp(-a[24]*SOCb[i])+(a[25]+a[26]*SOCb[i]+a[27]*SOCb[i]*SOCb[i]+a[28]*(pow(SOCb[i],3)))-a[29]*Crb+a[30]*Crb*Crb))
            V1b.append((V_mb[0]-V0b[0]-Icb*Rb[i])/(C1b[i]*(1/C1b[i]+1/C2b[i])))            
            V2b.append((V_mb[0]-V0b[0]-Icb*Rb[i])/(C2b[i]*(1/C1b[i]+1/C2b[i])))          
            V_cb.append(V0b[i]+Icb*Rb[i]+Icb*R1b[i]*(1-exp(-tcb/(R1b[i]*C1b[i])))+Icb*R2b[i]*(1-exp(-tcb/(R2b[i]*C2b[i])))+V1b[i]*(exp(-tcb/(R1b[i]*C1b[i])))+V2b[i]*(exp(-tcb/(R2b[i]*C2b[i]))))
            fitness+=abs(V_cb[i]-V_mb[i])
    except(OverflowError):
        return(10000) 
    if data_extract==1:
        numpy.savetxt("Vc"+data_text+".csv",numpy.array(V_c),delimiter=",")
        numpy.savetxt("R1"+data_text+".csv",numpy.array(R1),delimiter=",")
        numpy.savetxt("R2"+data_text+".csv",numpy.array(R2),delimiter=",")
        numpy.savetxt("C"+data_text+".csv",numpy.array(C),delimiter=",")
        numpy.savetxt("V0"+data_text+".csv",numpy.array(V0),delimiter=",")       
    for i in range(D):
        if a[i]>ub:
            fitness+=100*(a[i]-ub)**2              
        if a[i]<lb:
            fitness+=100*(a[i]-lb)**2          
    return(fitness)                

def ga_model1(guy):
    D=31
    battery_no=1
    state=0
    CDr=1
    data=BM_data_fun(state=state,CDr=CDr).T
    V_m=data[1]
    SOC_m=data[0]

    ub=50
    lb=0
    def fit(guy):
        return BM_fun1(guy.genes,D,battery_no=battery_no,state=state,Cr=CDr,V_m=V_m,SOC_m=SOC_m,lb=lb[0],ub=ub[0])              
    return fit

    
    