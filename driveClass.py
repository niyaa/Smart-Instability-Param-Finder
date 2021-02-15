
# coding: utf-8

# # Stability Driver
# Runs the stability calculation with varing parameters

# In[ ]:

import os
import fileinput

import numpy as np

from scipy import interpolate
from scipy import  optimize
from scipy.optimize import curve_fit

import time
import signal
import subprocess

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from IPython.core.display import display, HTML
import time
display(HTML("<style>.container { width:100% !important; }</style>"))


exe='mpirun -np 1 '
#exe='/home/sgepner/Projects/nektar/build_hyp/dist/bin/'
exe=''

# In[ ]:

'''Stores the result'''
class Case:
    def __init__(self):
        self.S = 0.0
        self.alpha = 0.0
        self.beta = 0.0
        self.lz = 0.0
        self.Re = 0.0
        self.sig_i = 0.0
        self.sig_r = 0.0
        self.sig_r1 = 0.0
    def setupxml(self):
        filexml = "stab.xml"
        sr = -1.0 * self.sig_r
        si = self.sig_i
        #print sr
        for line in fileinput.input(filexml, inplace=1):
            if line.startswith("      <P> Re"):
                print ("      <P> Re            = "+str(self.Re)+"   </P>")
            elif line.startswith("      <P> LZ"):
                print ("      <P> LZ            = "+str(self.lz)+"   </P>")
            elif line.startswith("      <P> realShift"):
                print ("      <P> realShift     = "+str(si)+"   </P>")
            elif line.startswith("      <P> imagShift"):
                print ("      <P> imagShift     = "+str(sr)+"   </P>")
            else:
                print (line, end='')
                
    def setupscr(self):
        filescr = "geom.scr"
        #print sr
        #self.alpha=c.alpha
        #self.S=c.S
        for line in fileinput.input(filescr, inplace=1):
            if line.startswith("a="):
                print ("a="+str(self.alpha)+";")
            elif line.startswith("S="):
                print ("S="+str(self.S)+";")
            else:
                print (line, end='');

    def run(self):
        global p
        p=subprocess.Popen(exe+'IncNavierStokesSolver geom.xml stab.xml', shell=True)
        
        pid = os.getpid()
        print ("Running case for Re="+str(self.Re)+" Lz="+str(self.lz)+" sr= ", str(-1.0 * self.sig_r));
        #print "my pid", pid
        #print "process runnng as ", p.pid
        #poll = p.poll()
        #if poll == None:
        #    print "Process alive"
        
        p.wait()
        print ("Finished!");
        
        #poll = p.poll()
        #if poll != None:
        #    print poll, "Process dead"
        
        #command = "rm *.fld *.chk"
        #os.system(command)
             

    
    def readevl(self):
        fileevl = "geom.evl"
        readnext = False
        readline = ""
        si = []
        sr = []
        with open(fileevl, "r") as ins:
            for line in ins:
                if line.startswith("         Real        Imaginary "):
                    readnext = True
                    continue
                if readnext:
                    readline = line
                    readline.split()
                    si.append(float(readline.split()[2]))
                    sr.append(float(readline.split()[3]))
        si, sr = zip(*sorted(zip(si, sr)))
        #print si, sr
        
        self.sig_i = si[-1]
        self.sig_r = sr[-1]
        self.sig_r1 = sr[-2]
        print ("Re="+str(self.Re)+" Lz="+str(self.lz)+" (si, sr): ("+str(self.sig_i)+","+str(self.sig_r)+","+ str(self.sig_r1)+ ")"+ "delta r:", self.sig_r - self.sig_r1)
        #command = "rm *.evl"
        #os.system(command)
        
    def readscr(self):
        readline = ""
        with open("geom.scr", "r") as ins:
            for line in ins:
                if line.startswith("S="):
                    readline = line
                    self.S=float(readline.split("=")[1][:-2])
                if line.startswith("a="):
                    readline = line
                    self.alpha=float(readline.split("=")[1][:-2])
                            
    def cleanLast(self):
        subprocess.call("rm geom.msh geom.xml *bak* *.evl",shell=True);

            


# In[162]:

def criticalReBeta3(d,prnt=0):
    
    y=d[:,4];
    Re=list(set(y));
    Re.sort();
    betaCr=[];
    sigmaMax=[];
    crsigmaR=[];
    ReNew=[];
    for i in Re:
        xx=d[d[:,4]==i][:,2];
        yy=d[d[:,4]==i][:,5];
        zz=d[d[:,4]==i][:,6];
        l1, l2 = zip(*sorted(zip(xx, yy)))
        xx=list(l1);yy=list(l2);
        if(prnt==1):print(i);
        xx=np.asarray(xx);yy=np.asarray(yy);
        if len(xx)>3:
            if(prnt==1):print(xx); print(yy);
            f=interpolate.interp1d(xx,-yy,kind='linear');
            f1=interpolate.interp1d(xx,zz,kind='linear')
            try:
                betaTemp=optimize.fmin(f,xx[np.argmax(yy)],maxiter=100000);
                betaCr.append(betaTemp[0]);
                sigmaMax.append(-f(betaTemp).flatten()[0]);
                ReNew.append(i);
                crsigmaR.append(f1(betaTemp).flatten()[0]);

            except:
                if(prnt==1):print(str(i)+" No \n");
                continue;

    l3, l2, l1, l0 = zip(*sorted(zip(ReNew,sigmaMax,betaCr, crsigmaR)));
    sigmaMax=list(l2); betaCr=list(l1);ReNew=list(l3); crsigmaR=list(l0);
    if(prnt==1):print('The most important 1. RE 2. sigmaMax 3. betacr 4. crsigmaR \n');
    if(prnt==1):print(ReNew); print(sigmaMax); print(betaCr); print(crsigmaR);
    Re=ReNew;

    f=interpolate.interp1d(Re,sigmaMax,kind='linear');
    f1=interpolate.interp1d(Re,betaCr,kind='linear');
    f2=interpolate.interp1d(Re,crsigmaR,kind='linear');

    s=np.where(np.diff(np.sign(sigmaMax)));
    if(prnt==1):print("the sign change "+str(s))
    recr=optimize.brentq(f,Re[s[0][0]],Re[s[0][0]+1]);
    betacr=f1(recr);
    crsigmaR=f2(recr)
    if(prnt==1):print(Re);print('\n'); print(sigmaMax); print('\n'); print(betaCr);
    print("The Converged values of the critical Re and Beta are " +str(recr) +" ------------ "+str(betacr));
    return recr, betacr.tolist(), crsigmaR.tolist();


########################################################################################################################################################
def unstableLz(Res, lzStart, siStart, srStart, sStart, alStart, clist, Maxitr, ReStep, fileWr=1):
    clist=[]
    filename =  "al-"+str(alStart)+"-.dat" 
    Lz = lzStart; si = siStart; sr = srStart
    ResRet=0; LzRet=Lz; siRet=0; srRet=0

    dLz = 0.2
    #dLz=0.2;
    maxList=[]; reList=[]; lzList=[]; sigList=[]; srList=[]
    sigMax=-1000;  countItr = 0
    
    numberOfSigmaPositive = 0
    while((sigMax < 0) and (Lz > 1) and (Lz < 1000) and (Res > 0) and (Res < 6000) and (countItr < Maxitr)  ):
        run = True
        lzAtSiMax=0
        srAtSiMax=0
        sigMax=-1000
        sigLowLimit = -5e-4

        while(run):
            c = Case()
            c.S = sStart
            c.alpha = alStart
            c.Re=Res
            c.lz = Lz
            c.beta = 2.0*np.pi / Lz
            c.sig_i = si
            c.sig_r = sr-sr*0.1
            c.setupxml()
            c.run()
            c.readevl()
            
            if(fileWr):
                text_file=open(filename,"a");
                text_file.write("%2.2f, %2.2f, %f, %3.2f, %d, %0.12f, %f, %f, %f, %f \n" % (c.alpha, c.S, c.beta, c.lz, c.Re, c.sig_i, c.sig_r, c.sig_r1, si, sr-sr*0.1))
                print("print 1 \n");
                text_file.close()

            si = c.sig_i
            sr = c.sig_r
            
            if ( (si < sigLowLimit) or (countItr > Maxitr) or (Lz < 1) or (Lz > 1000) or (Res <10) or (Res > 6000) ):
                run = False
            
            if( c.sig_i >sigMax): 
                sigMax = c.sig_i; lzAtSiMax = c.lz; srAtSiMax = c.sig_r; numberOfSigmaPositive = numberOfSigmaPositive + 1;

            Lz = Lz + dLz*Lz
            si = c.sig_i
            sr = c.sig_r
                    
        run = True
        Lz = lzStart-dLz*lzStart
        si = siStart
        sr = srStart
        while((run) and  (Lz > 1) and (Lz < 1000) ):
            c = Case()
            c.S = sStart
            c.alpha = alStart
            c.Re = Res
            c.lz = Lz
            c.beta = 2.0*np.pi / Lz
            c.sig_i = si
            c.sig_r = sr-sr*0.1;
            c.setupxml()
            c.run()
            c.readevl()

            if(fileWr):
                text_file=open(filename,"a");
                text_file.write("%2.2f, %2.2f, %f, %3.2f, %d, %0.12f, %f, %f, %f, %f \n" % (c.alpha, c.S, c.beta, c.lz, c.Re, c.sig_i, c.sig_r, c.sig_r1, si, sr-sr*0.1))
                print("print 2 \n");
                text_file.close()

            si = c.sig_i
            sr = c.sig_r
            Lz = Lz - dLz*Lz

            if ( (si < sigLowLimit) or (countItr > Maxitr) or (Lz < 1) or (Lz > 1000) or (Res <10) or (Res > 6000) ):
                run = False
            
            if( c.sig_i >sigMax): 
                sigMax = c.sig_i; lzAtSiMax = c.lz; srAtSiMax = c.sig_r;
                

        Res = Res+ReStep; Lz = lzStart = round(lzAtSiMax,4); si = siStart = sigMax; sr = srStart = srAtSiMax;
        countItr = countItr + 1;

    sigLowLimit = -5e-4  
    if(countItr < 2):
        Res = Res - 2*ReStep
        while((sigMax>0)  and  (Lz > 1) and (Lz < 1000) and (Res >0) and (Res < 6000) and (countItr < Maxitr) ):
            sigMax=-1000;
            run = True;
            
            while(run):
                c = Case()
                c.S = sStart
                c.alpha = alStart
                c.Re= Res;
                c.lz = Lz
                c.beta = 2.0*np.pi / Lz
                c.sig_i = si
                c.sig_r = sr-sr*0.1
                c.setupxml()
                c.run()
                c.readevl()

                if(fileWr):
                    text_file=open(filename,"a");
                    text_file.write("%2.2f, %2.2f, %f, %3.2f, %d, %0.12f, %f, %f, %f, %f \n" % (c.alpha, c.S, c.beta, c.lz, c.Re, c.sig_i, c.sig_r, c.sig_r1, si, sr-sr*0.1))
                    print("print 8 \n");
                    text_file.close()
                
                si = c.sig_i
                sr = c.sig_r
                Lz = Lz + dLz*Lz

                if ( (si < sigLowLimit) or (countItr > Maxitr) or (Lz < 1) or (Lz > 1000) or (Res <10) or (Res > 6000) ):
                    run = False
                
                if( c.sig_i > sigMax): 
                    sigMax = c.sig_i; lzAtSiMax = c.lz; srAtSiMax = c.sig_r;
            
            run = True
            Lz = lzStart-dLz*lzStart
            
            while(run):
                c = Case()
                c.S = sStart
                c.alpha = alStart
                c.Re = Res
                c.lz = Lz
                c.beta = 2.0*np.pi / Lz
                c.sig_i = si
                c.sig_r = sr-sr*0.1
                c.setupxml()
                c.run()
                c.readevl()
                                    
                if(fileWr):
                    text_file=open(filename,"a");
                    text_file.write("%2.2f, %2.2f, %f, %3.2f, %d, %0.12f, %f, %f, %f, %f \n" % (c.alpha, c.S, c.beta, c.lz, c.Re, c.sig_i, c.sig_r, c.sig_r1, si, sr-sr*0.1))
                    print("print 9 \n");
                    text_file.close()

                si = c.sig_i
                sr = c.sig_r
                Lz = Lz - dLz*Lz

                if ( (si < sigLowLimit) or (countItr > Maxitr) or (Lz < 1) or (Lz > 1000) or (Res <10) or (Res > 6000) ):
                    run = False

                if( c.sig_i  > sigMax):
                    sigMax = c.sig_i; lzAtSiMax = c.lz; srAtSiMax = c.sig_r;

            Res = Res - ReStep; Lz = lzStart = round(lzAtSiMax,4); si = siStart = sigMax; sr = srStart = srAtSiMax;
            countItr = countItr + 1;
        
        Res = Res + 2*ReStep



    return  Res, lzStart, siStart, srStart;	


####################################################################################################################

        
def alContin(alStart, alEnd, alStep, sStart, ReStart, lzStart, siStart, srStart, clist,maxItr, ReStep):
    nStep=int(abs(alEnd - alStart)/abs(alStep))
    for i in range(0,nStep):
        al=alStart+i*alStep;
        al=round(al,2);
        pwd=os.getcwd();
        if not os.path.exists(os.getcwd()+'/'+str(al)):
            os.makedirs(str(al));
        os.chdir(str(al));
        subprocess.call('cp ../geom.scr .',shell=True);
        subprocess.call('cp ../stab.xml .',shell=True);
        subprocess.call('cp ../base.xml .',shell=True);
        subprocess.call('cp ../base3d.xml .',shell=True);
        c=Case();
        c.S=sStart;
        c.alpha=al;
        c.setupscr();
        global p
        p=subprocess.Popen('gmsh -2 -order 9 geom.scr', shell=True)
        (output, err) = p.communicate()  

        #This makes the wait possible
        p_status = p.wait()

        p=subprocess.Popen('NekMesh geom.msh geom.xml', shell=True)
        
        (output, err) = p.communicate()  
        p_status = p.wait()
        
        p=subprocess.Popen('ADRSolver geom.xml base.xml', shell=True)

        (output, err) = p.communicate()  
        p_status = p.wait()
    
        p=subprocess.Popen('mv geom.fld geomB2D.fld', shell=True)
        (output, err) = p.communicate()  
        p_status = p.wait()
        
        p=subprocess.Popen('IncNavierStokesSolver geom.xml base3d.xml', shell=True)
        (output, err) = p.communicate()  
        p_status = p.wait()

        p=subprocess.Popen('mv geom.fld geom.bse', shell=True)
        (output, err) = p.communicate()  

        ReStart, lzStart, siStart, srStart = unstableLz(ReStart, lzStart, siStart, srStart, sStart, al, clist, maxItr, ReStep, 1);
        
        os.chdir(pwd);

        text_file=open('al-cont.dat',"a");
        text_file.write("%2.2f, %2.2f, %f, %3.2f, %d, %0.12f, %f \n" % (round(c.alpha,2), c.S, (2*np.pi)/lzStart, round(lzStart,2), ReStart, siStart, srStart))
        text_file.close()







def SContin(alStart, Slist, ReStart, lzStart, siStart, srStart, clist, maxItr, ReStep):
    for ss in Slist:
        #ss=sStart+i*sStep;
        #ss=round(ss,2);
        pwd=os.getcwd();
        if not os.path.exists(os.getcwd()+'/'+str(ss)):
            os.makedirs(str(ss));
        os.chdir(str(ss));
        subprocess.call('cp ../geom.scr .',shell=True);
        subprocess.call('cp ../stab.xml .',shell=True);
        subprocess.call('cp ../base.xml .',shell=True);
        subprocess.call('cp ../base3d.xml .',shell=True);
        c=Case();
        c.S=ss;
        c.alpha=alStart;
        c.setupscr();
        global p
        p=subprocess.Popen('gmsh -2 -order 9 geom.scr', shell=True)
        (output, err) = p.communicate()  

        #This makes the wait possible
        p_status = p.wait()

        p=subprocess.Popen('NekMesh geom.msh geom.xml', shell=True)
        
        (output, err) = p.communicate()  
        p_status = p.wait()

        
        p=subprocess.Popen('ADRSolver geom.xml base.xml', shell=True)

        (output, err) = p.communicate()  
        p_status = p.wait()
    
        p=subprocess.Popen('mv geom.fld geomB2D.fld', shell=True)
        (output, err) = p.communicate()  
        p_status = p.wait()
        
        p=subprocess.Popen('IncNavierStokesSolver geom.xml base3d.xml', shell=True)
        (output, err) = p.communicate()  
        p_status = p.wait()

        p=subprocess.Popen('mv geom.fld geom.bse', shell=True)
        (output, err) = p.communicate()  

        ReStart, lzStart, siStart, srStart = unstableLz(ReStart, lzStart, siStart, srStart, ss, alStart, clist, maxItr, ReStep, 1);

        os.chdir(pwd);

        text_file=open('al-cont.dat',"a");
        text_file.write("%2.2f, %2.2f, %f, %3.2f, %d, %0.12f, %f \n" % (round(c.alpha,2), c.S, (2*np.pi)/lzStart, round(lzStart,2), ReStart, siStart, srStart))
        text_file.close()




