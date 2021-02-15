
# coding: utf-8

# # Stability Driver
# Runs the stability calculation with varing parameters

# In[ ]:

import numpy as np
import sys, subprocess;
sys.path.append('~/');
import driveClass 
from glob import glob       
import os


# In[163]:
c=driveClass.Case();


path = glob("*/");
path.remove("__pycache__/")

cwd = os.getcwd();

for alpha in path:
    os.chdir(alpha);


    subprocess.call("mkdir 0.7",shell=True);
    subprocess.call("cp geom.xml geom.bse al-* 0.7/",shell=True);
    path1=glob("al-*");
    
    a = np.loadtxt(path1[0], delimiter=",");
    

    a = a[a[:,5].argmax()]  

    #alphas
    alStart= a[0] ;
    #S
    sStart=0.7; 
    #sigmas
    siStart= a[5]; srStart= a[6];
    #Reynolds
    ReStart = a[4]; ReStep= 0.1*a[4];
    #betas
    lzStart= a[3];

    clist=[];

    Slist=[ round(i*0.1,1) for i in range(1, 7)]

    Slist.reverse()
    #Slist = Slist + [0.09, 0.08, 0.07 , 0.06, 0.05] 


    driveClass.SContin(alStart, Slist, ReStart, lzStart, siStart, srStart, clist, 40, ReStep)

    os.chdir(cwd);








