# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 09:33:19 2016

@author: Avinash
"""
import numpy
from math import *

def BM_data_fun(state=1,data_set=1):
    ''' 
    
    1. Choose State=0 for charging and 1 for discharging
    2. CDr is used for both both Cr (charging rate) and Dr (discharging rate) depeding on state
     
    '''
    global Va,Vb,Vc
    global SOCa,SOCb,SOCc
    global V, SOC
    global I,t_total
    global Ta,Tb,Tc

    dataa=numpy.genfromtxt('tempset'+str(data_set)+'dis'+str(20)+'.csv', delimiter=',')
    datab=numpy.genfromtxt('tempset'+str(data_set)+'dis'+str(40)+'.csv', delimiter=',')
    datac=numpy.genfromtxt('tempset'+str(data_set)+'dis'+str(60)+'.csv', delimiter=',')
    
    s1=[0,4,9,15,19,23,27,32,39,49,59,69,79,89,99,109,119,129,139,149,159,169,179,185]
    s2=[0,4,9,15,19,23,27,32,39,49,59,69,79,89,99,109,119,129,139,149,159,169,174]
    s3=[0,4,9,15,19,23,27,32,39,49,59,69,79,89,99,109,119,129,139,149,159,169,179,184]
    
    Va=numpy.zeros(shape=(len(s1),))
    Vb=numpy.zeros(shape=(len(s2),))
    Vc=numpy.zeros(shape=(len(s3),))
    SOCa=numpy.zeros(shape=(len(s1),))
    SOCb=numpy.zeros(shape=(len(s2),))
    SOCc=numpy.zeros(shape=(len(s3),))
    
    for i in range(len(s1)):
        Va[i]=dataa[s1[i]][2]
        SOCa[i]=dataa[s1[i]][1]
        
    for i in range(len(s2)):
        Vb[i]=datab[s2[i]][2]
        SOCb[i]=datab[s2[i]][1]
    
    for i in range(len(s3)):
        Vc[i]=datac[s3[i]][2]
        SOCc[i]=datac[s3[i]][1]
        
    V=numpy.zeros(shape=(len(s1)+len(s2)+len(s3),))
    SOC=numpy.zeros(shape=(len(s1)+len(s2)+len(s3),))
    
    for i in range(len(s1)):
        V[i]=dataa[s1[i]][2]
        SOC[i]=dataa[s1[i]][1]
        
    for i in range(len(s2)):
        V[len(s1)+i]=datab[s2[i]][2]
        SOC[len(s1)+i]=datab[s2[i]][1]
    
    for i in range(len(s3)):
        V[len(s1)+len(s2)+i]=datac[s3[i]][2]
        SOC[len(s1)+len(s2)+i]=datac[s3[i]][1]
    I=0.156
    t_total=1
    Ta=20
    Tb=40
    Tc=60    
    
def BM_temp_low(al,mode):
    fitness=0
    for i in range(len(Va)):
        t=(SOCa[len(Va)]/SOC[i])*t_total
        Rol=al[0]*(SOCa[i]**4)+al[1]*(SOCa[i]**3)+al[2]*(SOCa[i]**2)+al[3]*(SOCa[i])+al[4]
        Rsl=al[5]*exp(-al[6]*SOCa[i])+al[7]+al[8]*SOCa[i]
        Csl=al[9]*(SOCa[i]**3)+al[10]*(SOCa[i]**2)+al[11]*(SOCa[i])+al[12]
        Rll=al[13]*exp(-al[14]*SOCa[i])+al[15]+al[16]*SOCa[i]
        Cll=al[17]*(SOCa[i]**6)+al[18]*(SOCa[i]**5)+al[19]*(SOCa[i]**4)+al[20]*(SOCa[i]**3)+al[21]*(SOCa[i]**2)+al[22]*(SOCa[i])+al[23]
        Tsl=Rsl*Csl
        Tll=Rll*Cll
        Vsl=Rsl*I
        Vll=Rll*I
        Vta=Vsl(1-exp(-t/Tsl))+Vll(1-exp(-t/Tll))+I*Rol
        fitness+=abs(Va[i]-Vta)
    return(fitness)

def BM_temp_high(a,al,mode):
    fitness=0
    for i in range(len(Vb)):
        t=(SOCb[len(Vb)]/SOC[i])*t_total 
        Rol=al[0]*(SOCb[i]**4)+al[1]*(SOCb[i]**3)+al[2]*(SOCb[i]**2)+al[3]*(SOCb[i])+al[4]
        Rsl=al[5]*exp(-al[6]*SOCb[i])+al[7]+al[8]*SOCb[i]
        Csl=al[9]*(SOCb[i]**3)+al[10]*(SOCb[i]**2)+al[11]*(SOCb[i])+al[12]
        Rll=al[13]*exp(-al[14]*SOCb[i])+al[15]+al[16]*SOCb[i]
        Cll=al[17]*(SOCb[i]**6)+al[18]*(SOCb[i]**5)+al[19]*(SOCb[i]**4)+al[20]*(SOCb[i]**3)+al[21]*(SOCb[i]**2)+al[22]*(SOCb[i])+al[23]
        dT=Tb-Ta        
        Roh=Rol*a[0]*exp(a[1]/(Tb-a[2]))
        Rsh=Rsl+a[3]*dT+a[4]*dT*SOCb[i]
        Csh=Csl+a[5]*dT*SOCb[i]+a[6]*dT
        Rlh=Rll+(a[7]*dT)*exp(-al[15]*SOCb[i])+a[8]*dT
        Clh=Cll*a[9]*exp(a[10]/T)
        Tsh=Rsh*Csh
        Tlh=Rlh*Clh
        Vsh=Rsh*I
        Vlh=Rlh*I
        Vtb=Vsh(1-exp(-t/Tsh))+Vlh(1-exp(-t/Tlh))+I*Roh
        fitness+=abs(Vb[i]-Vtb)
    for i in range(len(Vc)):
        t=(SOCc[len(Vc)]/SOC[i])*t_total 
        Rol=al[0]*(SOCc[i]**4)+al[1]*(SOCc[i]**3)+al[2]*(SOCc[i]**2)+al[3]*(SOCc[i])+al[4]
        Rsl=al[5]*exp(-al[6]*SOCc[i])+al[7]+al[8]*SOCc[i]
        Csl=al[9]*(SOCc[i]**3)+al[10]*(SOCc[i]**2)+al[11]*(SOCc[i])+al[12]
        Rll=al[13]*exp(-al[14]*SOCc[i])+al[15]+al[16]*SOCc[i]
        Cll=al[17]*(SOCc[i]**6)+al[18]*(SOCc[i]**5)+al[19]*(SOCc[i]**4)+al[20]*(SOCc[i]**3)+al[21]*(SOCc[i]**2)+al[22]*(SOCc[i])+al[23]
        dT=Tc-Ta        
        Roh=Rol*a[0]*exp(a[1]/(Tc-a[2]))
        Rsh=Rsl+a[3]*dT+a[4]*dT*SOCc[i]
        Csh=Csl+a[5]*dT*SOCc[i]+a[6]*dT
        Rlh=Rll+(a[7]*dT)*exp(-al[15]*SOCc[i])+a[8]*dT
        Clh=Cll*a[9]*exp(a[10]/T)
        Tsh=Rsh*Csh
        Tlh=Rlh*Clh
        Vsh=Rsh*I
        Vlh=Rlh*I
        Vtb=Vsh(1-exp(-t/Tsh))+Vlh(1-exp(-t/Tlh))+I*Roh
        fitness+=abs(Vb[i]-Vtb)
    return(fitness)     