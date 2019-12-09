#!/usr/bin/python
import sys, math, numpy, os
from subprocess import call

#################################################
# PaCS MD for GROMACS
# Number of data should to be the same with
#        number of frames
# DATE updated: Apr 7th 2015
# Apr. 8 2015: change for loop for m to for m in range(1,nbin+1)
# Apr. 24 2015: adding if there is itp file, copy that file to working directory
# May 18 2015: adding restart to code
# June 16 2015: adding cutoff stop for distance in z direction
#################################################

nbin=30
nround=100
rest=-1      #set this < 0 to begin the new PaCSMD simulation 
restep=2 
nroundadd=100 
comdistmax=7.0

groupA="Protein"
groupB="triNAG"
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



gmxcmd="mpirun -np 1 gmx "  #gmx serial calling
gmxcmd2="mpirun -np 1 gmx "#call MPI task in this variable
#use negative number for system choice of MD running conf.
gpu=-1
gpuid="0011"
ntomp=4  

wdir=os.getcwd()
outfn=wdir+"/"+outfnpf


#################################################
# Subroutine is below:
#################################################

def checkcycle():
	n=1
	chkfile=True
	while chkfile:
		print("Currently check cycle "+str(n)+"!")
		failcyc=0
		for cnt in range(1,nbin+1):
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/"+outfnpf+"-"+str(n)+"-"+str(cnt)+".tpr")):
				chkfile=False 
				failcyc=cnt 
				break
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/"+outfnpf+"-"+str(n)+"-"+str(cnt)+".gro")):
				chkfile=False 
				failcyc=cnt 
				break
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/"+outfnpf+"-"+str(n)+"-"+str(cnt)+".xtc")):
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
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/"+outfnpf+"-"+str(n)+"-"+str(cnt)+".edr")):
				chkfile=False 
				failcyc=cnt 
				break
			if not(os.path.exists(outfn+"-"+str(n)+"-"+str(cnt)+"/"+outfnpf+"-"+str(n)+"-"+str(cnt)+".log")):
				chkfile=False 
				failcyc=cnt 
				break
		n=n+1
	print("Simulation was terminated at CYCLE "+str(n-1)+" and REPLICA "+str(cnt)+". Simulation will restart at this cycle.")
	return n-1
			


#################################################
# Main program is below:
#################################################

#write selection file
f=open(wdir+"/sel.dat","w")
f.write('com of group "'+groupA+'" plus com of group "'+groupB+'"')
f.close()

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
		os.system(gmxcmd+" grompp -f "+outfn+"-0-0/"+mdinitfn+" -c "+outfn+"-0-0/"+grofn+" -t "+outfn+"-0-0/"+vfn+" -p "+outfn+"-0-0/"+topolfn+" -o "+outfn+"-0-0/"+outfnpf+"-0-0.tpr -maxwarn 10")
	else:
                os.system(gmxcmd+" grompp -f "+outfn+"-0-0/"+mdinitfn+" -c "+outfn+"-0-0/"+grofn+" -p "+outfn+"-0-0/"+topolfn+" -o "+outfn+"-0-0/"+outfnpf+"-0-0.tpr -maxwarn 10")
	print("################################")
	if gpu>=0 and ntomp>0:
		os.system(gmxcmd2+" mdrun -deffnm "+outfn+"-0-0/"+outfnpf+"-0-0 -v -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
	else:
		os.system(gmxcmd2+" mdrun -deffnm "+outfn+"-0-0/"+outfnpf+"-0-0 -v -ntomp "+str(ntomp))
	#apply nopbc for calculating the rmsd
	os.system("echo 'System' | "+gmxcmd+" trjconv -s "+outfn+"-0-0/"+outfnpf+"-0-0.tpr -f "+outfn+"-0-0/"+outfnpf+"-0-0.xtc -o "+outfn+"-0-0/"+outfnpf+"-0-0-noPBC.xtc -pbc mol -ur compact")
	print("echo 'System' | "+gmxcmd+" trjconv -s "+outfn+"-0-0/"+outfnpf+"-0-0.tpr -f "+outfn+"-0-0/"+outfnpf+"-0-0.xtc -o "+outfn+"-0-0/"+outfnpf+"-0-0-noPBC.xtc -pbc mol -ur compact")
	#calculating the rmsd of the first cycle and pick up the 10 best one
	os.system(gmxcmd+" distance -f "+outfn+"-0-0/"+outfnpf+"-0-0-noPBC.xtc -s "+outfn+"-0-0/"+outfnpf+"-0-0.tpr  -n "+ndxfn+" -oall "+outfn+"-0-0/"+outfnpf+"-0-0.xvg -xvg none -tu ps -sf "+wdir+"/sel.dat")
	print(gmxcmd+" distance -f "+outfn+"-0-0/"+outfnpf+"-0-0-noPBC.xtc -s "+outfn+"-0-0/"+outfnpf+"-0-0.tpr  -n "+ndxfn+" -oall "+outfn+"-0-0/"+outfnpf+"-0-0.xvg -xvg none -tu ps -sf "+wdir+"/sel.dat")
	#reading the rmsd
	comdist=numpy.loadtxt(outfn+"-0-0/"+outfnpf+"-0-0.xvg")
	comdistcp=numpy.concatenate((comdist,numpy.zeros((len(comdist[:,1]),1))),1)
	comdistcp=comdistcp[numpy.lexsort((comdistcp[:,0],comdistcp[:,1]))]
	#writing the ranking to log file
	for m in range(1,nbin+1):
		os.system("echo '"+str(comdistcp[len(comdistcp[:,0])-m,2])+"   "+str(comdistcp[len(comdistcp[:,0])-m,0])+"   "+str(comdistcp[len(comdistcp[:,0])-m,1])+"'  >> "+wdir+"/"+logfn)
		print(str(comdistcp[len(comdistcp[:,0])-m,2])+"   "+str(comdistcp[len(comdistcp[:,0])-m,0])+"   "+str(comdistcp[len(comdistcp[:,0])-m,1])+"\n")
	
else:
	#os.system("echo 'Simulation is going restart at ROUND "+str(nround)+" '  >> "+wdir+"/"+logfn)
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
	#print(linesfn[1][17:(len(linesfn[1]))]+".")
	bn=int(linesfn[1][17:(len(linesfn[1]))])
	print("BIN of current PACS MD is ",bn)
	#print(linesfn[2][19:(len(linesfn[2]))]+".")
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
		cdtemp=numpy.loadtxt(outfn+"-"+str(restep)+"-"+str(m)+"/"+outfnpf+"-"+str(restep)+"-"+str(m)+".xvg")
		cdcptemp1=[]
		cdcptemp1=numpy.concatenate((cdtemp,numpy.ones((len(cdtemp[:,1]),1))*(m)),1)
		#cdcptemp1=cdcptemp1[numpy.lexsort((cdcptemp1[:,0],cdcptemp1[:,1]))]
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
	for m in range(1,nbin+1):
		
		#create directory
		os.system("mkdir "+outfn+"-"+str(n)+"-"+str(m))
		
		#dump frame to folder for MD run
		os.system("echo 'System' | "+gmxcmd+" trjconv -f "+outfn+"-"+str(n-1)+"-"+str(int(comdistcp[len(comdistcp[:,0])-m,2]))+"/"+outfnpf+"-"+str(n-1)+"-"+str(int(comdistcp[len(comdistcp[:,0])-m,2]))+".xtc  -s "+outfn+"-"+str(n-1)+"-"+str(int(comdistcp[len(comdistcp[:,0])-m,2]))+"/"+outfnpf+"-"+str(n-1)+"-"+str(int(comdistcp[len(comdistcp[:,0])-m,2]))+".tpr -o "+outfn+"-"+str(n)+"-"+str(m)+"/input.gro -dump "+str(float(comdistcp[len(comdistcp[:,0])-m,0]))  )
		
		#copy needed file for MD run
		os.system("cp -r "+wdir+"/"+mdfn+" "+outfn+"-"+str(n)+"-"+str(m)+"/")
		os.system("cp -r "+wdir+"/"+topolfn+" "+outfn+"-"+str(n)+"-"+str(m)+"/")
		if len(itpfn)>0:
			os.system("cp -r "+wdir+"/"+itpfn+" "+outfn+"-"+str(n)+"-"+str(m)+"/")
		os.system("cp -r "+wdir+"/index.ndx"+" "+outfn+"-"+str(n)+"-"+str(m)+"/")

		#running simulation with random vel (set in mdp file)
		os.system(gmxcmd+" grompp -f "+outfn+"-"+str(n)+"-"+str(m)+"/"+mdfn+" -c "+outfn+"-"+str(n)+"-"+str(m)+"/input.gro -o "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".tpr -maxwarn 10")
		if gpu>=0 and ntomp>0:
			os.system(gmxcmd2+" mdrun -deffnm "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+" -v -ntomp "+str(ntomp)+" -gpu_id "+str(gpuid))
		else:
			os.system(gmxcmd2+" mdrun -deffnm "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+" -v -ntomp "+str(ntomp))
		
		#check the distribution and add to temperary array

		#apply nopbc for calculating the CV
		
		os.system("echo 'System' | "+gmxcmd+" trjconv -s "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".tpr -f "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".xtc -o "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+"-noPBC.xtc -pbc mol -ur compact")
		
		print("echo 'System' | "+gmxcmd+" trjconv -s "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".tpr -f "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".xtc -o "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+"-noPBC.xtc -pbc mol -ur compact")


		#calculating the CV of the run

		os.system(gmxcmd+" distance -f "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+"-noPBC.xtc -s "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".tpr -n "+ndxfn+" -oall "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".xvg -xvg none -tu ps -sf "+wdir+"/sel.dat")
		print(gmxcmd+" distance -f "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+"-noPBC.xtc -s "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".tpr -n "+ndxfn+" -oall "+outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".xvg -xvg none -tu ps -sf "+wdir+"/sel.dat")

		#reading the CV
		cdtemp=numpy.loadtxt(outfn+"-"+str(n)+"-"+str(m)+"/"+outfnpf+"-"+str(n)+"-"+str(m)+".xvg")
		cdcptemp1=[]
		cdcptemp1=numpy.concatenate((cdtemp,numpy.ones((len(cdtemp[:,1]),1))*(m)),1)
		#cdcptemp1=cdcptemp1[numpy.lexsort((cdcptemp1[:,0],cdcptemp1[:,1]))]
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

