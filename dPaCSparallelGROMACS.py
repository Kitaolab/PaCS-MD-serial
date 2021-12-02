from __future__ import print_function
import sys, math, numpy, os, re
from subprocess import call
from multiprocessing import Pool

#################################################
# PaCS MD for GROMACS
# Number of data should to be the same with
#        number of frames
# DATE updated: Apr 7th 2015
# Apr. 8 2015: change for loop for m to for m in range(1,nbin+1)
# Apr. 24 2015: adding if there is itp file, copy that file to working directory
# May 18 2015: adding restart to code
# June 16 2015: adding cutoff stop for distance in z direction
# Jan 24 2021: adding support for parallel multidir option in GROMACS:
#				to run in serial mode, set runmode to 1
#				to run in parallel mode, set the runmode to number of parallel process
# Oct 27 2021: fixing the mixing multiprocessing of python and denial executibility of MPI process on Fugaku
# Dec 02 2021: rearrange the code
#################################################

nbin=30          #number of replicas of PaCS-MD to run 
nround=100       #number of cycles of PaCS-MD to run
rest=0          #set this < 0 to begin the new PaCSMD simulation, == 0 for automatic restart
restep=2         #this parameter only use when rest is set greater than 0, this give the specific step to run
nroundadd=100    #this parameter only use when rest is set greater than 0, this is the number of cycle to add to
comdistmax=7.0   #stop PaCS-MD when reaching this value 

######################CV setting#############################

groupA="Protein"
groupB="DNA"

#####################input files setting#####################

grofn="input.gro"
mdfn="mdpacs.mdp"
mdinitfn="md.mdp"
vfn=""
ndxfn="index.ndx"
topolfn="topol.top"
itpfn="*.itp"
comdistfn="comdist"
outfnpf="pacs"
logfn="pacs.log"
ndxfn="index.ndx"
selfn="seltext.txt"
skipdata=2 

#######################clean trajectory option###############

clntrj=True
keepgroup="non-Water"
clnimdtraj=True

##################GROMACS run-setting########################

runmode=30 #using for multidir option in GROMACS
runmode2=4  #using for multiprocessing package of python within a node.
gmxcommand=" gmx_mpi "
mpicommand="mpirun -np " 
usemp=True #use Multiprocessing library from Python
if usemp==True:
	gmxcmd="gmx_serial "   #change it to exact non-MPI version of GROMACS
else:
	gmxcmd=mpicommand+" 1 "+gmxcommand  #gmx serial calling
gmxcmd2=mpicommand+str(runmode)+" "+gmxcommand #call MPI task in this variable
gpu=-1 #use negative number for system choice of MD running conf.
gpuid="0011"
ntomp=50   
wdir=os.getcwd()
outfn=wdir+"/"+outfnpf

#################################################
# Subroutine is below:
#################################################

def gmxcmd2lastloop(runmodelastloop):
	tmpstr=mpicommand+str(runmodelastloop)+" gmx_mpi "   
	return tmpstr

def GROcheck():
	os.system(gmxcmd+" --version | grep 'GROMACS version' > gmxversion ")
	f=open("gmxversion","r")
	ln=f.readlines()
	vers=ln[0].split()[-1].split(".")[0]
	return vers
	
def writesel(fn,txt):
	f=open(fn,"w")
	f.write(txt)
	f.close()
	return fn

def checkcycle():
	n=1
	chkfile=True
	while chkfile:
		print("Currently check cycle "+str(n)+"!")
		failcyc=0
		for cnt in range(1,nbin+1):
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/topol.tpr")):
				chkfile=False 
				failcyc=cnt 
				break
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/confout.gro")):
				chkfile=False 
				failcyc=cnt 
				break
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/traj_comp.xtc")):
				chkfile=False 
				failcyc=cnt 
				break
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/"+outfnpf+"-"+str(n)+"-"+str(cnt)+".xvg")):
				chkfile=False 
				failcyc=cnt 
				break
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/input.gro")):
				chkfile=False 
				failcyc=cnt 
				break
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/ener.edr")):
				chkfile=False 
				failcyc=cnt 
				break
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/md.log")):
				chkfile=False 
				failcyc=cnt 
				break
		n=n+1
	print("Simulation was terminated at CYCLE "+str(n-1)+" and REPLICA "+str(cnt)+". Simulation will restart at this cycle.")
	return n-1

def fn(m):
    global n,comdistcp
	#create directory
    os.system("mkdir "+outfn+"-"+str(n)+"-"+str(m))
	#dump frame to folder for MD run
    if usemp==True:
        os.system("echo 'System' | "+gmxcmd+" trjconv -f "+outfn+"-"+str(n-1)+"-"+str(int(comdistcp[len(comdistcp[:,0])-m,2]))+"/traj_comp.xtc  -s "+outfn+"-"+str(n-1)+"-"+str(int(comdistcp[len(comdistcp[:,0])-m,2]))+"/topol.tpr -o "+outfn+"-"+str(n)+"-"+str(m)+"/input.gro -dump "+str(float(comdistcp[len(comdistcp[:,0])-m,0]))  )
    else:
        os.system(mpicommand+" 1 --stdin "+wdir+"/"+selfn+".system"+ gmxcommand+" trjconv -f "+outfn+"-"+str(n-1)+"-"+str(int(comdistcp[len(comdistcp[:,0])-m,2]))+"/traj_comp.xtc  -s "+outfn+"-"+str(n-1)+"-"+str(int(comdistcp[len(comdistcp[:,0])-m,2]))+"/topol.tpr -o "+outfn+"-"+str(n)+"-"+str(m)+"/input.gro -dump "+str(float(comdistcp[len(comdistcp[:,0])-m,0]))  )
    #copy needed file for MD run
    os.system("cp -r "+wdir+"/"+mdfn+" "+outfn+"-"+str(n)+"-"+str(m)+"/")
    os.system("cp -r "+wdir+"/"+topolfn+" "+outfn+"-"+str(n)+"-"+str(m)+"/")
    if len(itpfn)>0:
        os.system("cp -r "+wdir+"/"+itpfn+" "+outfn+"-"+str(n)+"-"+str(m)+"/")
    os.system("cp -r "+wdir+"/index.ndx"+" "+outfn+"-"+str(n)+"-"+str(m)+"/")
    #running simulation with random vel (set in mdp file)
    if vers <=2016:
        os.system(gmxcmd+" grompp -f "+outfn+"-"+str(n)+"-"+str(m)+"/"+mdfn+" -c "+outfn+"-"+str(n)+"-"+str(m)+"/input.gro -n "+ndxfn+" -o "+outfn+"-"+str(n)+"-"+str(m)+"/topol.tpr -maxwarn 10")
    else:
        os.system(gmxcmd+" grompp -f "+outfn+"-"+str(n)+"-"+str(m)+"/"+mdfn+" -c "+outfn+"-"+str(n)+"-"+str(m)+"/input.gro -n "+ndxfn+" -o "+outfn+"-"+str(n)+"-"+str(m)+"/topol.tpr -r "+grofn+" -maxwarn 10")
    return 

def f2(m):
    global n
    if usemp==True:
        #apply nopbc for calculating the CV
        os.system("echo 'non-Water' | "+gmxcmd+" trjconv -s "+outfn+"-"+str(n)+"-"+str(m)+"/topol.tpr -f "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp.xtc -o "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-nowat.xtc -pbc mol -ur compact")
        os.system("mv "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-nowat.xtc "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp.xtc ")
        #remove the noPBC trajectories
        os.system("rm -r "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-noPBC.xtc")
    else:
        #apply nopbc for calculating the CV
        os.system(mpicommand+" 1 --stdin "+wdir+"/"+selfn+".keepgroup "+ gmxcommand+" trjconv -s "+outfn+"-"+str(n)+"-"+str(m)+"/topol.tpr -f "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp.xtc -o "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-nowat.xtc -pbc mol -ur compact")
        os.system("mv "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-nowat.xtc "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp.xtc ")
        #remove the noPBC trajectories
        os.system("rm -r "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-noPBC.xtc")
    return

def f1a(m):
    global n
    #apply nopbc for calculating the CV
    if usemp==True:
        os.system("echo 'non-Water' | "+gmxcmd+" trjconv -s "+outfn+"-"+str(n)+"-"+str(m)+"/topol.tpr -f "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp.xtc -o "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-nowat.xtc -pbc mol -ur compact")
        os.system("mv "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-nowat.xtc "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp.xtc ")
    else:
        os.system(mpicommand+" 1 --stdin "+wdir+"/"+selfn+".keepgroup "+ gmxcommand+" trjconv -s "+outfn+"-"+str(n)+"-"+str(m)+"/topol.tpr -f "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp.xtc -o "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-nowat.xtc -pbc mol -ur compact")
        os.system("mv "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-nowat.xtc "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp.xtc ")
    return

def f1b(m):
    global n
    #remove the noPBC trajectories
    os.system("rm -r "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-noPBC.xtc")
    return

def f(m):
    global n
    #apply nopbc for calculating the CV
    if usemp==True:
        os.system("echo 'System' | "+gmxcmd+" trjconv -s "+outfn+"-"+str(n)+"-"+str(m)+"/topol.tpr -f "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp.xtc -o "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-noPBC.xtc -pbc mol -ur compact")
        #calculating the CV of the run
        os.system(gmxcmd+" distance -f "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-noPBC.xtc -s "+outfn+"-"+str(n)+"-"+str(m)+"/topol.tpr -n "+ndxfn+" -oall "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".xvg -xvg none -tu ps -sf "+wdir+"/sel.dat")
    else:
        os.system(mpicommand+" 1 --stdin "+wdir+"/"+selfn+".system "+ gmxcommand+" trjconv -s "+outfn+"-"+str(n)+"-"+str(m)+"/topol.tpr -f "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp.xtc -o "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-noPBC.xtc -pbc mol -ur compact")
        #calculating the CV of the run
        os.system(gmxcmd+" distance -f "+outfn+"-"+str(n)+"-"+str(m)+"/traj_comp-noPBC.xtc -s "+outfn+"-"+str(n)+"-"+str(m)+"/topol.tpr -n "+ndxfn+" -oall "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".xvg -xvg none -tu ps -sf "+wdir+"/sel.dat")
    return

#################################################
# Main program is below:
#################################################

#check version of GROMACS
#vers=GROcheck()
vers=2021.3 
#vers=re.sub("[^\d\.]", "", GROcheck())

######write selection file#####
writesel(wdir+"/sel.dat",'com of group "'+groupA+'" plus com of group "'+groupB+'"\n')
#system
writesel(wdir+"/"+selfn+".system",'System\n')
#keepgroup
writesel(wdir+"/"+selfn+".keepgroup",keepgroup+"\n")

if rest<0:
	#save the old log file
	if os.path.exists(wdir+"/"+logfn):
		i=0
		while(os.path.exists(wdir+"/"+logfn+"."+str(i))):
			i=i+1 
		os.system("mv "+wdir+"/"+logfn+" "+wdir+"/"+logfn+"."+str(i))
	os.system("echo 'PaCS MD is going to run ' > "+wdir+"/"+logfn)
	print("PaCS MD is going to run \n")
	os.system("echo 'Number of BIN is  "+str(nbin)+" ' >> "+wdir+"/"+logfn)
	print("echo 'Number of BIN is  "+str(nbin)+"\n")
	os.system("echo 'Number of ROUND is "+str(nround)+" ' >> "+wdir+"/"+logfn)
	print("echo 'Number of ROUND is "+str(nround)+"\n")
	os.system("echo 'Simulation is running from the beginning. '  >> "+wdir+"/"+logfn)
	print("Simulation is running from the beginning. \n")
	#call gromacs to run for the first cycle
	os.system("echo '++++++++++++++++'  >> "+wdir+"/"+logfn)
	print('++++++++++++++++\n')
	os.system("echo 'ROUND 0. '  >> "+wdir+"/"+logfn)
	print('ROUND 0. \n')
	os.system("echo '1111 Max to Min Ranking 1111'  >> "+wdir+"/"+logfn)
	print('+ Ranking +\n' )
	os.system("echo '+ Bin +++ Step +'  >> "+wdir+"/"+logfn)
	print('+ Bin +++ Step +\n')
	os.system("mkdir "+outfn+"-0-0")
	#copy needed file for MD run
	os.system("cp -r "+wdir+"/"+mdinitfn+" "+outfn+"-0-0/")
	os.system("cp -r "+wdir+"/"+topolfn+" "+outfn+"-0-0/")
	if len(itpfn)>0:
		os.system("cp -r "+wdir+"/"+itpfn+" "+outfn+"-0-0/")
	os.system("cp -r "+wdir+"/"+grofn+" "+outfn+"-0-0/") 
	if len(vfn)>0:
		os.system("cp -r "+wdir+"/"+vfn+" "+outfn+"-0-0/")
		if vers<=2016:
			os.system(gmxcmd+" grompp -f "+outfn+"-0-0/"+mdinitfn+" -c "+outfn+"-0-0/"+grofn+" -t "+outfn+"-0-0/"+vfn+" -p "+outfn+"-0-0/"+topolfn+" -n "+ndxfn+" -o "+outfn+"-0-0/topol.tpr -maxwarn 10")
		else:
			os.system(gmxcmd+" grompp -f "+outfn+"-0-0/"+mdinitfn+" -c "+outfn+"-0-0/"+grofn+" -t "+outfn+"-0-0/"+vfn+" -p "+outfn+"-0-0/"+topolfn+" -n "+ndxfn+" -o "+outfn+"-0-0/topol.tpr -r "+grofn+" -maxwarn 10")
	else:
		if vers<=2016:
					os.system(gmxcmd+" grompp -f "+outfn+"-0-0/"+mdinitfn+" -c "+outfn+"-0-0/"+grofn+" -p "+outfn+"-0-0/"+topolfn+" -n "+ndxfn+" -o "+outfn+"-0-0/topol.tpr -maxwarn 10")
		else:
					os.system(gmxcmd+" grompp -f "+outfn+"-0-0/"+mdinitfn+" -c "+outfn+"-0-0/"+grofn+" -p "+outfn+"-0-0/"+topolfn+" -n "+ndxfn+" -o "+outfn+"-0-0/topol.tpr -r "+grofn+" -maxwarn 10")		
	print("################################")
	if gpu>=0 and ntomp>0:
		os.system(gmxcmd2+" mdrun -s "+outfn+"-0-0/topol.tpr -o "+outfn+"-0-0/traj.trr -x "+outfn+"-0-0/traj_comp.xtc  -e "+outfn+"-0-0/ener.edr  -g "+outfn+"-0-0/md.log  -c "+outfn+"-0-0/confout.gro  -cpo "+outfn+"-0-0/state.cpt -pme gpu -npme 1 -v -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
	else:
		os.system(gmxcmd2+" mdrun -s "+outfn+"-0-0/topol.tpr -o "+outfn+"-0-0/traj.trr -x "+outfn+"-0-0/traj_comp.xtc  -e "+outfn+"-0-0/ener.edr  -g "+outfn+"-0-0/md.log  -c "+outfn+"-0-0/confout.gro  -cpo "+outfn+"-0-0/state.cpt -pme gpu -npme 1 -v -ntomp "+str(ntomp))
	#apply nopbc for calculating the rmsd
	os.system(mpicommand+" 1 --stdin "+wdir+"/"+selfn+".system "+ gmxcommand+" trjconv -s "+outfn+"-0-0/topol.tpr -f "+outfn+"-0-0/traj_comp.xtc -o "+outfn+"-0-0/traj_comp-noPBC.xtc -pbc mol -ur compact")
	print(mpicommand+" 1 --stdin "+wdir+"/"+selfn+".system "+ gmxcommand+" trjconv -s "+outfn+"-0-0/topol.tpr -f "+outfn+"-0-0/traj_comp.xtc -o "+outfn+"-0-0/traj_comp-noPBC.xtc -pbc mol -ur compact")
	#calculating the rmsd of the first cycle and pick up the 10 best one
	os.system(gmxcmd+" distance -f "+outfn+"-0-0/traj_comp-noPBC.xtc -s "+outfn+"-0-0/topol.tpr  -n "+ndxfn+" -oall "+outfn+"-0-0/"+outfnpf+"-0-0.xvg -xvg none -tu ps -sf "+wdir+"/sel.dat")
	print(gmxcmd+" distance -f "+outfn+"-0-0/traj_comp-noPBC.xtc -s "+outfn+"-0-0/topol.tpr  -n "+ndxfn+" -oall "+outfn+"-0-0/"+outfnpf+"-0-0.xvg -xvg none -tu ps -sf "+wdir+"/sel.dat")
	if clnimdtraj==True:
		os.system("rm -r "+outfn+"-0-0/traj_comp-noPBC.xtc")
	#reading the rmsd
	comdist=numpy.loadtxt(outfn+"-0-0/"+outfnpf+"-0-0.xvg",skiprows=skipdata)
	comdistcp=numpy.concatenate((comdist,numpy.zeros((len(comdist[:,1]),1))),1)
	comdistcp=comdistcp[numpy.lexsort((comdistcp[:,0],comdistcp[:,1]))]
	#writing the ranking to log file
	for m in range(1,nbin+1):
		os.system("echo '"+str(comdistcp[len(comdistcp[:,0])-m,2])+"   "+str(comdistcp[len(comdistcp[:,0])-m,0])+"   "+str(comdistcp[len(comdistcp[:,0])-m,1])+"'  >> "+wdir+"/"+logfn)
		print(str(comdistcp[len(comdistcp[:,0])-m,2])+"   "+str(comdistcp[len(comdistcp[:,0])-m,0])+"   "+str(comdistcp[len(comdistcp[:,0])-m,1])+"\n")
	
else:
	print("Simulation is going restart at ROUND "+str(nround)+ "\n" )
	#checking the need of restart file or quit the program
	if not(os.path.exists(wdir+"/"+logfn)):
		print(logfn+" file is not found! Please check the file.")
		exit()
	f=open(wdir+"/"+logfn,'r')
	linesfn=f.readlines()
	f.close()
	for m in range(0,len(linesfn)-1):
		linesfn[m]=linesfn[m].rstrip('\n')
	bn=int(linesfn[1][17:(len(linesfn[1]))])
	print("BIN of current PACS MD is ",bn)
	if rest!=0:
		rnd=int(linesfn[2][19:(len(linesfn[2]))])
		print("ROUND of current PACS MD is ",rnd)
		if restep>rnd-1:
			print("Restart round is larger than the number of round that we have. Please recheck")
			exit()
	elif rest==0:
		rnd=checkcycle()
		restep=rnd-1
		print("ROUND of current PACS MD is ",rnd)
	#re-executing the MD code for regenerating the last cycle
	n=restep
	if runmode==1:
		for m in range(1,nbin+1):
			if gpu>=0 and ntomp>0:
				os.system(gmxcmd2+" mdrun -deffnm "+outfn+"-"+str(n)+"-"+str(m)+"/topol -v -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
			else:
				os.system(gmxcmd2+" mdrun -deffnm "+outfn+"-"+str(n)+"-"+str(m)+"/topol -v -ntomp "+str(ntomp))
	elif runmode>1: 
		mdloop=nbin//runmode
		print("Expect number of mdloop "+str(mdloop))
		if (mdloop*runmode)<nbin:
			lastloop=nbin%runmode
			mdloop=mdloop
		else:
			lastloop=0
		print("Expect number of lastloop "+str(lastloop))
		for x in range(0,mdloop):
			multidir=" "
			for m in range(x*runmode+1,(x+1)*runmode+1):
				multidir=multidir+outfn+"-"+str(n)+"-"+str(m)+"  "
			if gpu>=0 and ntomp>0:
				os.system(gmxcmd2+" mdrun -multidir "+multidir+" -s topol -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
				print(gmxcmd2+" mdrun -multidir "+multidir+" -s topol -v -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
			else:
				os.system(gmxcmd2+" mdrun -multidir "+multidir+" -s topol -ntomp "+str(ntomp))		
				print(gmxcmd2+" mdrun -multidir "+multidir+" -s topol -v -ntomp "+str(ntomp))
		if lastloop>0:
			multidir=" "
			for m in range(mdloop*runmode+1,nbin+1):
				multidir=multidir+outfn+"-"+str(n)+"-"+str(m)+"  "
			print(gmxcmd2lastloop(lastloop))
			if gpu>=0 and ntomp>0:
				os.system(gmxcmd2lastloop(lastloop)+" mdrun -multidir "+multidir+" -s topol -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
			else:
				os.system(gmxcmd2lastloop(lastloop)+" mdrun -multidir "+multidir+" -s topol -ntomp "+str(ntomp))		
	# rewrite the PACS MD log:
	if os.path.exists(wdir+"/"+logfn):
		i=0
		while(os.path.exists(wdir+"/"+logfn+"."+str(i))):
			i=i+1
		os.system("mv "+wdir+"/"+logfn+" "+wdir+"/"+logfn+"."+str(i))
	f=open(wdir+"/"+logfn,'w')
	wrtstepneed=nbin*(restep+1)+4*(restep+1)+4
	while wrtstepneed>len(linesfn):
		rnd=rnd-1
		restep=rnd
		wrtstepneed=nbin*(restep+1)+4*(restep+1)+4
		print("Due to the lack of data in "+wdir+"/"+logfn+". Resetting ROUND of current PACS MD is ",rnd)
	for m in range(0,wrtstepneed):
		if m==2:
			f.write(linesfn[2][0:18]+" "+str(int(restep+nroundadd))+"\n")
		elif (restep!=0) and (m==wrtstepneed-1-nbin-1):
			f.write(linesfn[m]+"\n")
		elif m==(wrtstepneed-1):
			f.write(linesfn[m]+" ")
		else:
			f.write(linesfn[m]+"\n")
	f.close()
	#rebuild the comdist array
	cdcptemp=[] #clean buffer array
	for m in range(1,nbin+1):
		if not(os.path.exists(outfn+"-"+str(restep)+"-"+str(m)+"/"+outfnpf+"-"+str(restep)+"-"+str(m)+".xvg")):
			print("Distance file is not found! Please check the file!")
			exit()
		cdtemp=numpy.loadtxt(outfn+"-"+str(restep)+"-"+str(m)+"/"+outfnpf+"-"+str(restep)+"-"+str(m)+".xvg",skiprows=skipdata)
		cdcptemp1=[]
		cdcptemp1=numpy.concatenate((cdtemp,numpy.ones((len(cdtemp[:,1]),1))*(m)),1)
		if m>1:
			cdcptemp =numpy.concatenate((cdcptemp,cdcptemp1),0)
		else:
			cdcptemp=cdcptemp1
	comdistcp=cdcptemp[numpy.lexsort((cdcptemp[:,0],cdcptemp[:,1]))]
#check for restart flag
if rest<0:
	nbloop=1
else:
	nbloop=restep+1
	nround=restep+nroundadd+1


n=nbloop
while n<nround:
#for n in range(nbloop,nround):
	cdcptemp=[] #clean buffer array
	os.system("echo '++++++++++++++++'  >> "+wdir+"/"+logfn)
	print('++++++++++++++++\n')
	os.system("echo 'ROUND "+str(n)+". '  >> "+wdir+"/"+logfn)
	print("ROUND "+str(n)+"\n" )
	os.system("echo '1111 Max to Min Ranking 1111'  >> "+wdir+"/"+logfn)
	print("++++ Ranking +++")
	os.system("echo '+ Bin +++ Step +'  >> "+wdir+"/"+logfn)
	print("+ Bin +++ Step +")
	#preparing the tpr files:
	if usemp==True:
		p=Pool(runmode2)
		with p:
			p.map(fn,[tmpvar for tmpvar in range(1,nbin+1)])  
		p.join()  
		p.close()   
	else:
		for tmpvar in range(1,nbin+1):
			fn(tmpvar)
	#cleanning trajectories if being selected
	if ((clntrj==True) and (clnimdtraj==True)):
		if usemp==True:
			p=Pool(runmode2)
			with p:
				p.map(f2,[tmpvar for tmpvar in range(1,nbin+1)] )    
			p.join()  
			p.close() 
		else:
			for tmpvar in range(1,nbin+1):
				f2(tmpvar)
	elif ((clntrj==True) and (clnimdtraj==False)):
		if usemp==True:
			p=Pool(runmode2)
			with p:
				p.map(f1a,[tmpvar for tmpvar in range(1,nbin+1)] )    
			p.join()  
			p.close() 
		else:
			for tmpvar in range(1,nbin+1):
				f1a(tmpvar)
	elif ((clntrj==False) and (clnimdtraj==True)):
		if usemp==True:
			p=Pool(runmode2)
			with p:
				p.map(f1b,[tmpvar for tmpvar in range(1,nbin+1)] )    
			p.join()  
			p.close() 		
		else:	
			for tmpvar in range(1,nbin+1):
				f1b(tmpvar)
	#executing the MD code
	if runmode==1:
		for m in range(1,nbin+1):
			if gpu>=0 and ntomp>0:
				os.system(gmxcmd2+" mdrun -deffnm "+outfn+"-"+str(n)+"-"+str(m)+"/topol -v -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
			else:
				os.system(gmxcmd2+" mdrun -deffnm "+outfn+"-"+str(n)+"-"+str(m)+"/topol -v -ntomp "+str(ntomp))
	elif runmode>1: 
		mdloop=nbin//runmode
		print("Expect number of mdloop "+str(mdloop))
		if (mdloop*runmode)<nbin:
			lastloop=nbin%runmode
			mdloop=mdloop
		else:
			lastloop=0
		print("Expect number of lastloop "+str(lastloop))
		for x in range(0,mdloop):
			multidir=" "
			for m in range(x*runmode+1,(x+1)*runmode+1):
				multidir=multidir+outfn+"-"+str(n)+"-"+str(m)+"  "
			if gpu>=0 and ntomp>0:
				os.system(gmxcmd2+" mdrun -multidir "+multidir+" -s topol -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
				print(gmxcmd2+" mdrun -multidir "+multidir+" -s topol -v -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
			else:
				os.system(gmxcmd2+" mdrun -multidir "+multidir+" -s topol -ntomp "+str(ntomp))		
				print(gmxcmd2+" mdrun -multidir "+multidir+" -s topol -v -ntomp "+str(ntomp))
		if lastloop>0:
			multidir=" "
			for m in range(mdloop*runmode+1,nbin+1):
				multidir=multidir+outfn+"-"+str(n)+"-"+str(m)+"  "
			print(gmxcmd2lastloop(lastloop))
			if gpu>=0 and ntomp>0:
				os.system(gmxcmd2lastloop(lastloop)+" mdrun -multidir "+multidir+" -s topol -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
			else:
				os.system(gmxcmd2lastloop(lastloop)+" mdrun -multidir "+multidir+" -s topol -ntomp "+str(ntomp))		
	#check the distribution and add to temperary array
	if usemp==True:
		p=Pool(runmode2)
		with p:
			p.map(f,[tmpvar for tmpvar in range(1,nbin+1)] )    
		p.join()  
		p.close()
	else:
		for tmpvar in range(1,nbin+1):
			f(tmpvar)   
	for m in range(1,nbin+1):	
		#reading the CV
		cdtemp=numpy.loadtxt(outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".xvg",skiprows=skipdata)
		cdcptemp1=[]
		cdcptemp1=numpy.concatenate((cdtemp,numpy.ones((len(cdtemp[:,1]),1))*(m)),1)
		if m>1:
			cdcptemp =numpy.concatenate((cdcptemp,cdcptemp1),0)
		else:
			cdcptemp=cdcptemp1
	#Preparing for the next PaCS MD step
	#Ranking the trajectory
	comdistcp=cdcptemp[numpy.lexsort((cdcptemp[:,0],cdcptemp[:,1]))]
	#writing the ranking to log file
	for l in range(1,nbin+1):
		os.system("echo '"+str(comdistcp[len(comdistcp[:,0])-l,2])+"   "+str(comdistcp[len(comdistcp[:,0])-l,0])+"   "+str(comdistcp[len(comdistcp[:,0])-l,1])+"'  >> "+wdir+"/"+logfn)
		print(str(comdistcp[len(comdistcp[:,0])-l,2])+"   "+str(comdistcp[len(comdistcp[:,0])-l,0])+"   "+str(comdistcp[len(comdistcp[:,0])-l,1])+"\n")
	if ((comdistmax>0.0) and (comdistcp[len(comdistcp[:,2])-nbin-1,1]<comdistmax) and (n==(nround-1))):
		nround=nround+1
	if ((comdistmax>0.0) and (comdistcp[len(comdistcp[:,2])-nbin-1,1]>comdistmax) and (n<(nround-1))):		
		break
	n=n+1    #this command is for replacing "for" loop with "while" loop



