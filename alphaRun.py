
# coding: utf-8

# # Stability Driver
# Runs the stability calculation with varing parameters

# In[ ]:

import numpy as np
import sys;
sys.path.append('~/');
import driveClass2 
        


# In[163]:
c=driveClass2.Case();

#alphas
alStart= 2 ; alEnd=5; alStep=0.25;
#S
sStart=0.7; sEnd=0.4; sStep=0.4;
#sigmas
siStart= -0.00711398; srStart= 0.512488;
#Reynolds
ReStart = 110; ReStep=30;
#betas
lzStart= 8.45;

clist=[];

#ReStart, lzStart, siStart, srStart=unstableLz(ReStart, lzStart, siStart, srStart, sStart, alStart, clist, 200);





driveClass2.alContin(alStart, alEnd, alStep, sStart, ReStart, lzStart, siStart, srStart, clist, 40 , ReStep);











