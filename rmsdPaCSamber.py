from __future__ import print_function
import sys, math, numpy, os 

mdcmd="pmemd.cuda "    #command for pmemd 
anlcmd="cpptraj "      #command for CPPTRAJ
topfn="complex.parm7"  #the topology file
initfn="prod_1.rst7"   #the initial configuration file (should be well-relaxed)
confinitfn="md.mdin"   #MD setting file for the cyle 0
conffn="mdpacs.mdin"   #MD setting file for the cycle 1++
logfn="pacs.log"       #the log file for PaCS-MD simulation
centergroup=":9@CA"    #centering the atoms in the group so that protein does not broken in the main box
fittinggroup="@CA"     #fitting group for structural comparison
rmsdgroup=":C,O,N"     #group for RMSD calculation
posresfn=""  #leave it blank as "" if you don't use positional restraints
skipdata=2 #skip the datapoint in the files
cycle=200   #number of cycle to perform
nrep=30    #number of replicas
endchar="\n"
rankcrit=-1  # -1 for descending and 1 for ascending

global CVmat

def mdrun(n,m):
    for x in range(1,m+1):
        wd="pacs-"+str(n)+"-"+str(x)
        if n==0:
            conffile=confinitfn
        else:
            conffile=conffn
        if len(posresfn)==0:
            os.system(mdcmd+" -O -p "+topfn+" -i "+conffile+" -c "+wd+"/input.rst7 -o "+wd+"/pacs.mdout -r "+wd+"/pacs.rst7 -inf "+wd+"/pacs.mdinfo -x "+wd+"/pacs.nc ")
        else:
            os.system(mdcmd+" -O -p "+topfn+" -i "+conffile+" -c "+wd+"/input.rst7 -o "+wd+"/pacs.mdout -r "+wd+"/pacs.rst7 -inf "+wd+"/pacs.mdinfo -x "+wd+"/pacs.nc -ref "+wd+"/"+posresfn) 
    return

def writecalCV(n,mdx):
    wd="pacs-"+str(n)+"-"+str(mdx)
    f=open("calCV.cpptraj","w")
    f.write("parm "+topfn+endchar)
    f.write("trajin "+wd+"/pacs.nc "+endchar)
    f.write("center "+centergroup+endchar)
    f.write("image "+endchar)
    f.write("reference "+initfn+" [refstr]"+endchar)
    f.write("rms fit ref [refstr] @CA,O,N "+endchar)
    f.write("rms cal ref [refstr] @CA nofit out "+wd+"/pacs-"+str(n)+"-"+str(mdx)+".xvg"+endchar)
    f.write("run "+endchar)
    f.write("quit")
    f.close() 
    return

def calCV(n,m,rankcrit):
    global CVmat
    CVmat=[]
    for x in range(1,m+1):
        wd="pacs-"+str(n)+"-"+str(x)
        writecalCV(n,x)
        os.system(anlcmd+" -i calCV.cpptraj ") 
        CVmattmp1=numpy.loadtxt(wd+"/pacs-"+str(n)+"-"+str(x)+".xvg",skiprows=skipdata) 
        CVmattmp2=[]
        CVmattmp2=numpy.concatenate((CVmattmp1,numpy.ones((len(CVmattmp1[:,1]),1))*(x)),1)
        if x==1:
            CVmattmp=CVmattmp2
        else:
            CVmattmp=numpy.concatenate((CVmattmp,CVmattmp2),0)
    CVmat=CVmattmp[numpy.lexsort((CVmattmp[:,0],numpy.array(CVmattmp[:,1])*rankcrit))]

    return

def writeexportframe(n,rep):
    global CVmat
    if rep>len(CVmat):
        print("Number of replicas is larger than the available data points.")
        quit
    else:
        f=open("export.cpptraj","w")
        f.write("parm "+topfn+endchar)
        for x in range(1,rep+1):
            wd="pacs-"+str(n-1)+"-"+str(int(CVmat[x,2])) 
            f.write("trajin "+wd+"/pacs.nc "+endchar)
            wd="pacs-"+str(n)+"-"+str(x) 
            f.write("trajout "+wd+"/input.rst7 onlyframes "+str(int(CVmat[x,0]))+endchar  )
            f.write("run "+endchar)
            f.write("clear trajin "+endchar) 
            f.write("clear trajout "+endchar)
        f.write("quit"+endchar)
        f.close()
    return

def cycleinit(n,m):
    global CVmat
    if len(CVmat[:])==0:
        print("There is no frame to export.")
        exit
    else:
        writeexportframe(n,m)
        for x in range(1,m+1):
            wd="pacs-"+str(n)+"-"+str(x) 
            os.system("mkdir "+wd)
        os.system(anlcmd+" -i export.cpptraj")
    return

def writeCV(n,m,fn):
    global CVmat
    os.system("echo '++++++++++++++++'  >> "+fn)
    os.system("echo 'ROUND "+str(n)+". '  >> "+fn)
    os.system("echo '1111 Max to Min Ranking 1111'  >> "+fn)
    os.system("echo '+ Rep +++ Step +++ RMSD +'  >> "+fn)
    for x in range(0,m):
        os.system("echo '"+str(CVmat[x,2])+"   "+str(CVmat[x,0])+"   "+str(CVmat[x,1])+"'  >> "+fn)
    return

def main():
    global CVmat
    #initiating the first cycle
    n=0
    m=1
    os.system("mkdir pacs-"+str(n)+"-"+str(m)) 
    os.system("cp "+initfn+" pacs-"+str(n)+"-"+str(m)+"/input.rst7") 
    mdrun(0,1)
    calCV(0,1,rankcrit)
    writeCV(n,nrep,logfn)
    for n in range(1,cycle+1):
        cycleinit(n,nrep)
        mdrun(n,nrep)
        calCV(n,nrep,rankcrit)
        writeCV(n,nrep,logfn)
    return

main()

