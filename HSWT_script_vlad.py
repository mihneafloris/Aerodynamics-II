# -*- coding: utf-8 -*-

from flowtools_vlad import *
import numpy as np
import matplotlib.pyplot as plt

plt.close()
plt.style.use("ggplot")

#constants
pi    = 3.14159
gamma = 1.4

#input files
flist1=["1A2019-01-18_15-35-57.txt","2A2019-01-18_15-38-20.txt"]
geometry_file="geometry.txt"

#tunnel geometry

geometry=np.genfromtxt(geometry_file,delimiter=" ",comments="%")
xtunnel=geometry[:,0]
htab=geometry[:,1]
Atab=geometry[:,2]

def part1(file,Atab):
    name=str(file)
    data=np.genfromtxt(file,delimiter="\t",comments="%",skip_header=6)
    xtab=data[:,0]
    ptab=data[:,1]
    
    Ax0_Adot=flowisentropic2(gamma,ptab[0],'pres')[4]
    coef=Ax0_Adot/Atab[0]
    coef=1
    
    
    p_theory=[]
    M_theory=[]
    T_theory=[]
    for i in range(len(Atab)):
        Aratio=Atab[i]
        x=xtunnel[i]
        if x<65.0:    #since the throat is at x=65.00 [mm]
            Aratio=Aratio*coef
            output=flowisentropic2(gamma,Aratio,"sub")
            
        else:
            if name[0]=='1':
                output=flowisentropic2(gamma,Aratio,"sup")
            else:
                
                output=flowisentropic2(gamma,Aratio,"sub")
        p_theory.append(output[2])
        M_theory.append(output[0])
        T_theory.append(output[1])
    p_theory=np.array(p_theory)
    M_theory=np.array(M_theory)
    T_theory=np.array(T_theory)
    
    plt.figure(name)
    plt.subplot(3,1,1)
    plt.grid(True)
    plt.plot(xtunnel,Atab,'-r',linewidth=0.5)
    plt.xlabel(r"x [mm]")
    plt.ylabel(r"$A/A_t$ [-]")
    
    plt.subplot(3,1,2)
    #p measured
    plt.plot(xtab[xtab<=200.0],ptab[xtab<=200.0],'-r',marker='+',markevery=2,label=r"Measured $p/p_t$",linewidth=0.5)
    # p calculated 
    plt.plot(xtunnel[xtunnel>=44.80],p_theory[xtunnel>=44.80],'-b',marker='^',markevery=2,label=r"Theoretical $p/p_t$",linewidth=0.5)
    plt.xlabel(r"x [mm]")
    plt.ylabel(r"$p/p_t$ [-]")
    plt.grid(True)
    plt.legend()
    
    plt.subplot(3,1,3)
    #compute the M in the tunnel from the pressure
    M_calculated=[]
    for i in range(len(xtab)):
        pratio=ptab[i]
        x=xtab[i]
        output=flowisentropic2(gamma,pratio,"pres")
        M_calculated.append(output[0])
    M_calculated=np.array(M_calculated)

    #M in tunnel (calculated based on real pressure)
    plt.plot(xtab[xtab<=200.0],M_calculated[xtab<=200.0],'-r',marker='+',markevery=2,label="Measured Mach",linewidth=0.5)
    #M theoretical 
    plt.plot(xtunnel[xtunnel<=200.0],M_theory[xtunnel<=200.00],'-b',marker='^',markevery=2,label="Theoretical Mach",linewidth=0.5)
    plt.xlabel(r"x [mm]")
    plt.ylabel(r"M [-]" )
    plt.grid(True)
    plt.legend()

for file in flist1:
    part1(file,Atab)
    plt.show()
    
    
#Part 2
def tunnel_height(x):
    h0=8
    a=0.01692
    return a*x + h0
flist2=["3A2019-01-18_15-44-58.txt","4A2019-01-18_15-46-06.txt","5A2019-01-18_15-51-52.txt"]
i=0
graphlist=[['r','+'],['b','^'],['y','*']]
plt.figure("Experiment 3,4 and 5")
for file in flist2:
    data=np.genfromtxt(file,delimiter="\t",comments="%",skip_header=6)

    graph=graphlist[i]
    color=graph[0]
    symbol=graph[1]
    xtab=data[:,0]
    ptab=data[:,1]
    name="Experiment " + str(i+3)
    plt.plot(xtab,ptab,color,marker=symbol,markevery=2,label=name,linewidth=0.5)
    plt.axvline(x=405)
    plt.axvline(x=530)
    plt.axvline(x=630)
    plt.legend()
    plt.xlabel(r"x [mm]")
    plt.ylabel(r"$p/p_t$ [-]")
    plt.grid(True)
    i=i+1
    
hthroat=tunnel_height(630)
A_As=hthroat/14.6

output=flowisentropic2(gamma,A_As,'sub')
M_sub=output[0]
output1=flowisentropic2(gamma,A_As,'sup')
M_sup=output1[0]

output2=flownormalshock2(gamma,M_sup,'mach')
p_pt6=output2[5]*output2[6]
print(output2[2],output2[6])
p_pt5=output2[2]*output2[6]


output3=flowisentropic2(gamma,M_sub,'mach')
p_pt3=output3[2]
print(p_pt3,p_pt5,p_pt6)
plt.plot(530,0.77,"gx")
plt.plot(530,0.72,"gx")
plt.plot(530,0.28,"gx")
plt.plot(630,0.82,"gx")
plt.plot(630,0.75,"gx")
plt.plot(630,0.22,"gx")
plt.annotate(0.77,(540,0.77))
plt.annotate(0.72,(540,0.72))
plt.annotate(0.28,(540,0.28))
plt.annotate(0.82,(640,0.82))
plt.annotate(0.75,(640,0.75))
plt.annotate(0.22,(640,0.22))
plt.figure("Experiment 3,4 and 5")
plt.show()



