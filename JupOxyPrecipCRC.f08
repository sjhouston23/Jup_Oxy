program precipitation
!*******************************************************************************
!* Created by Stephen J. Houston 2.8.18
!*******************************************************************************
!* This program simulates the path of an energetic oxygen ion as it penetrates
!* into the Jovian atmosphere. Different initial energies are considered.
!* I use 1-25,000 keV/u as initial energies for this exploratory work.
!* A random pitch angle is considered for the precipitating ion.
!* A Monte Carlo simulation is used to determine the type of collision
!* and where the collision occurs. To calculate the collision path dN
!* I use 1-Prob = exp{-sigtot*dN}. After each collision I track the secondary
!* electrons produced, as well as the new ion energy. The results are binned by
!* column density which correspond to a particular altitude.
!******************************
!* External files required:
!* The secondary electron distributions are calculated by other codes based on
!* cross sections calculated by Dave Schultz.
!* The singly-differential cross sections for the ion+H2 processes are read in
!* at the beginning of the model. These are as function of energy.
!* The distribution functions for ejected electrons are read as eprobfunc and
!* aprobfunc. These are as function of energy and angle, respectively.
!* Ejected electron energy is calculated in this model by using the cross
!* sections and energy loss model presented in Schultz et al., 2018.
!* See input files for the files that are needed.
!******************************
!* Goals:
!* 1. Read in the distribution functions for all available collision types,
!* charge states and energies and store this info in matrices
!* 2. Set a matrix with the angular distribution to determine wether the
!* electron will be scattered forward or backward in a collision. The incident
!* ion pitch angle is added to the ejected electron angle.
!* 3. Read in all the total xs calculated by Dave.
!* 4. Create altitude bins and find the corresponding column density.
!* 5. Set the initial conditions for the ion -> charge state and initial energy,
!* incident angle (normally kept at 0) and initial pitch angle
!* 6. Follow ion as it penetrates the atmosphere and has collisions determined
!* by the MC.
!* 7. Track the charge state of the ion, number and energy of electrons
!* produced, and ion energy at each altitude bin in the atmosphere, until the
!* ion runs out of energy (E<1keV/u)
!*******************************************************************************

use, intrinsic :: ISO_FORTRAN_ENV
use formatting
implicit none

!**************************** Variable Declaration *****************************
integer i,j,k,l,m,n,run,ion !Do-loop variables
!* Parameter Variables:
integer number_of_energies !Number of initial ion energies
integer atmosLen !"Length" of the atmosphere (2998-200 km) with 2 km steps
integer k1,k2,in,lux !Seed and luxury level for RNG
integer nE2strBins,nOxEngBins !Number of 2-Stream and oxygen
integer nStopPowerEBins !Number of stopping power energy bins
integer nProc,nChS !Number of processes and charge states
integer ne,na !Number of electron energy and angle bins
integer nOutputFiles
integer stephen

real*8 oxEngBinSize !Size of oxygen bins
real*8 stopPowerEBinSize !Size of stopping power bins
real*8 mass !Mass of species (oxygen=16) used in energy loss equation

parameter(number_of_energies=12,atmosLen=1544,nE2strBins=260,nProc=36,nChS=10)
parameter(ne=2600,na=1800) !Number of electron energy and angle bins
parameter(k1=0,k2=0,lux=3) !lux set to 3 for optimal randomness and timeliness
parameter(nOxEngBins=5000,oxEngBinSize=5.0,mass=16.0) !Oxygen binning/mass
parameter(nStopPowerEBins=295,stopPowerEBinSize=10.0) !Stopping power bins
parameter(nOutputFiles=21) !Number of output files

integer trial,number_of_ions,energy
integer t1,t2,clock_maxTotal,clock_rateTotal !Used to calculate comp. time
integer hrs,min
integer t3,t4,clock_max,clock_rate !Used to calculate comp. time
integer numSim,dpt,maxDpt !Depth integers
integer initQ,tempQ,tempQold,excite !Charge states and excitation
integer addElect,bin,disso,process,processE,processC(36),PID(2),ds
integer processEnergies(12),binNo
integer CXProc(14) !Charge exchange processes
!* For large integer counts
integer(kind=int64) :: elect,totalElect,tElectFwd,tElectBwd !Total Electrons
integer(kind=int64) :: eCounts(nE2strBins),electFwdA(atmosLen) !Electrons fwd
integer(kind=int64) :: electFwdAE(atmosLen,nE2strBins) !Electrons fwd
integer(kind=int64) :: electBwdA(atmosLen),electBwdAE(atmosLen,nE2strBins) !bwd
integer(kind=int64) :: oxygen(nProc,atmosLen,nChS) !Oxygen counter
integer(kind=int64) :: oxygenCX(atmosLen,nChS) !Oxygen charge exchange
integer(kind=int64) :: H2Ex(atmosLen) !Excited H2 counter
integer(kind=int64) :: OxyVsEng(nChS,nOxEngBins) !Charge state fraction variable
integer(kind=int64) :: totO(nOxEngBins) !Total oxygen counter
integer(kind=int64) :: nSPions(nStopPowerEBins) !Normalization for stop power
integer(kind=int64) :: pnSPions(nStopPowerEBins,36,10)
integer(kind=int64),dimension(nProc,atmosLen) :: totHp,totH2p !H+ and H2+ counts

real*8 sec,sigTotOld,dEold !Time
real*8 E,dE,Eion(number_of_energies) !Array of initial input ion energies
real*8 es,del(nE2strBins),E2str(nE2strBins) !2-Stream energy bin creation
real*8 engBins(nOxEngBins),oxEngBins(nOxEngBins) !Oxygen energy bins
real*8 dEsp,delSP(nStopPowerEBins) !Stopping power dE
real*8,dimension(nStopPowerEBins) :: SigTotvsEng,dEvsEng,dNvsEng
real*8 stopPowerEBins(nStopPowerEBins) !Stopping power energy bins
real*8 SPvsEng(nStopPowerEBins) !Stopping power vs energy
real*8 dum,norm !Dummy variable, normalization
real*8 incB,kappa,avg(1000)
real*8 sigTot(nChS,25000),dN,dNTot,dZ,dZTot
real*8 eEnergy,eAngle,eEnergyTmp
real*8 eAngleSS,eEnergySS,eAngleDS(2),eEnergyDS(2)
real*8 ProcessdE(nStopPowerEBins,36,10)
real*8,dimension(6,9,13,ne) :: eProbFunc !Ejected electron probability function
real*8,dimension(6,9,13,na) :: aProbFunc
real*8,dimension(atmosLen,nE2strBins) :: prode2stF,prode2stB
!*Total column density, array of altitude in km, altitude bin size, scale height
real*8,dimension(atmosLen) :: totalCD,altitude,altDelta,totalDens,H

!* Outputs
integer(kind=int64) :: collisions(8,5),NSIM(8),SIM(5)
integer(kind=int64),dimension(atmosLen) :: totalHp,totalH2p
integer(kind=int64),dimension(nProc) :: pHp,npHp,pH2p,npH2p

character(len=100) filename,filenames(nOutputFiles) !Output file names
character(len=100) arg

real ranVecA(1000002) !Want randome vector to be plenty big
real pangle !Pitch angle variables
real dissRan !Random number to determine dissociation probability
real stpnuc
real,allocatable :: angle(:) !Want as many random number angles as ions later

real*8 delAVGe(790),elecEbins(790),eAngleAVG(790)
integer eAngleCounts(790),elecAbins(180),EAC(180)

logical open

!****************************** Data Declaration *******************************
!* Initial ion enegy input:
data Eion/1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0,5000.0,10000.0,&
          25000.0/
!dE for each 2-Stream energy bin. Must match two stream code binning
data del/20*0.5,70*1.0,10*2.0,20*5.0,10*10.0,20*10.0,10*50.0,10*100.0,40*200.0,&
         10*400,10*1000,10*2000,10*5000,10*10000.0/
data engBins/nOxEngBins*oxEngBinSize/ !Used for oxygen binning
data delSP/100*1.0,90*10,90*100,15*1000/ !Binning for stopping power (n=295)
data delAVGe/200*10,400*20,140*100,50*1000/
data processC/28,30,18,20,13,15,23,25,27,29,22,24,12,14,33,35,3,5,22,24,8,10,&
              16,26,11,23,24,2,4,21,7,9,31,1,6,36/
data processEnergies/1,10,50,75,100,110,140,190,200,230,280,295/!nStopPowerEBins
!* Charge exchange process numbers
data CXProc/2,4,6,8,9,10,12,14,20,23,24,25,30,36/
!* Output data file names: (nOutputFiles)
data filenames/'H+_Prod','H2+_Prod','H2_Excite_Prod','Oxy_Vs_Energy',&
'Stopping_Power','Processes','Oxy_Neg','Oxy0_','Oxy1_','Oxy2_','Oxy3_','Oxy4_',&
'Oxy5_','Oxy6_','Oxy7_','Oxy8_','Oxy_CX','XRay_DE','XRay_CX','2Str_Elect_Fwd',&
'2Str_Elect_Bwd'/
!********************************** Run Time ***********************************
!Calculate the total computational run time of the model:
call system_clock (t1,clock_rateTotal,clock_maxTotal)
!**************************** Initialize Variables *****************************
totalCD  =0.0;oxEngBins=0.0;stopPowerEBins=0.0;
altitude =0.0;sigTot   =0.0
altDelta =0.0;elect    =0
energy   =0  ;nSPions  =0
es       =0.0;
E2str    =0.0;
!**************************** Create the Atmosphere ****************************
open(unit=200,file='./Atmosphere/Input/JunoColumnDensity_2km.dat',status='old')
open(unit=201,file='./Atmosphere/Input/JunoAtmosphere_2km.dat',status='old')
read(200,*);read(201,*) !Skip header lines
do i=1,atmosLen
  read(200,*)altitude(i),dum,dum,dum,dum,totalCD(i)!Read in the atmosphere
  read(201,*)dum,dum,dum,dum,dum,totalDens(atmosLen-i+1),dum,dum,H(atmosLen-i+1)
  altDelta(i)=2.0
end do
close(200) !Close column density file
close(201) !Close atmosphere file
!*********************** Ejected Electron Probabilities ************************
!* Created by ReadElectDist.f08 and ProbDist.f08
open(unit=202,file='./NewElectronDist/ProbDistFunc/eprobfunc.dat',status='old')
open(unit=203,file='./NewElectronDist/ProbDistFunc/aprobfunc.dat',status='old')
read(202,*) eProbFunc
read(203,*) aProbFunc
close(202)
close(203)
!**************************** Various Bin Creation *****************************
!2-Stream energy bins:
do i=1,nE2strBins
  es=es+del(i)
  E2str(i)=Es
end do
!Oxygen bins for charge state fractions:
oxEngBins(1)=oxEngBinSize
do i=2,nOxEngBins
  oxEngBins(i)=oxEngBins(i-1)+engBins(i) !5-5000 keV/u.May need to go to 25MeV
end do
!Stopping power bins:
es=1.0
do i=1,nStopPowerEBins
  es=es+delSP(i)
  stopPowerEBins(i)=es
end do
!Bins for ejected electron energy:
es=0.0
do i=1,790
  es=es+delAVGe(i)
  elecEbins(i)=es
end do
!Bins for ejected electron angle:
do i=1,180
  elecAbins(i)=i
end do
!*******************************************************************************
!******************************** MAIN PROGRAM *********************************
!*******************************************************************************
!* The following run number corresponds to the energy in keV/u
!* Juno:
!* 1=10, 2=15, 3=20, 4=30, 5=45, 6=60, 7=75, 8=120, 9=220, 10=450, 11=500,
!* 12=750, 13=1000, 14=1250, 15=1500, 16=1750, 17=2000, 18=2500, 19=3000,
!* 20=4000, 21=5000, 22=10000, 23=25000
!* Regular:
!* 1=1, 2=10, 3=50, 4=75, 5=100, 6=200, 7=500, 8=1000, 9=2000, 10=5000,
!* 11=10000, 12=25000
!*******************************************************************************
number_of_ions=10!000

call get_command_argument(1,arg)
read(arg,'(I100)') trial
do run=9,9!1,number_of_energies
  call system_clock(t3,clock_rate,clock_max) !Comp. time of each run
  energy=int(Eion(run))
  write(filename,'("./Output/",I0,"keV/Seeds.dat")') energy
  open(unit=205,file=filename,status='unknown',access='append',action='write')
  write(205,*) trial
  close(205)
  write(filename,'("./Output/",I0,"keV/Overview",I0,".dat")') energy,trial
  open(unit=206,file=filename,status='unknown')
  write(206,*) "Number of ions:         ", number_of_ions
  write(206,*) "Initial energy:         ", energy, 'keV'
  write(206,*) "Trial number (RNG Seed):", trial
  write(206,*) "***************************************************************&
                ****************************"
!*************************** Random Number Generator ***************************
  !k1=0,k2=0 Should be set to zero unless restarting at a break (See ranlux.f08)
  in=trial !RNG seed
  write(206,*) "RNG Seed = ", in !~
  call rluxgo(lux,in,k1,k2)
  allocate(angle(number_of_ions))
  call ranlux(angle,number_of_ions) !Calculate all the angles to be used here
!********************** Reset Counters For New Ion Energy **********************
  totHp =0;totalElect=0;tElectFwd =0;tElectBwd =0;SPvsEng=0.0    ;nSPions=0
  totH2p=0;eCounts   =0;electFwdA =0;electBwdA =0;SigTotvsEng=0.0;maxDpt=0
  H2Ex  =0;oxygen    =0;electFwdAE=0;electBwdAE=0;dEvsEng=0.0
!************************ Ion Precipitation Begins Here ************************
  write(206,*) 'Starting Ion Precipitiaton: ', energy,'keV/u' !Double check eng
  flush(206)
  do ion=1,number_of_ions !Each ion starts here
    !*****************************
    !Reset Variables:

    !*****************************
    !Initial Conditions:
    pangle=0.0         !Reset the pitch angle for every run
    incB  =0.0         !Incident B-field
    kappa =0.0         !Used to account for pitch angle
    numSim=energy*1000 !Number of simulations for a single ion. Must be great !~
                       !enough to allow the ion to lose all energy
    E=Eion(run)        !Start with initial ion energy
    dE=0.0             !Energy loss
    initQ=3            !1 is an initial charge state of -1, 3 is +1
    tempQ=initQ        !Set the charge state variable that will be changed
    tempQold=initQ     !Need another charge state variable for energyLoss.f08
    dNTot=0.0          !Reset the column density to the top of the atm.
    dZTot=3000.0       !Start from the top of the atmosphere
    dpt=4              !Depth of penetration for bins. (integer value)
    !Beginning scale height (H) at 4 or 5 seems to be more accurate than 1-3
    l=0                !Used as index for dN calculation (ranVecA(l))
    process=0;excite=0 !CollisionSim outputs
    PID=0              !Process identification numbers
    sigTotOld=0.0
    !*****************************
    pangle=(2.0*atan(1.0))-acos(angle(ion)) !Pitch angle calculation has a
    !cosine dist. Straight down is pitch angle of 0, random number must be 0
    write(206,*) 'Ion Number: ', ion, 'Pitch angle: ', pangle*90/acos(0.0)
    flush(206)
    kappa=1.0/(cos(pangle)*cos(incB)) !Used to convert from ds to dz
    call ranlux(ranVecA,1000002) !Get a random vector for collisions
    do i=1,numSim !This loop repeats after each collision until E < 1 keV/u
      !*****************************
      !Reset Variables:
      eEnergy=0.0;eEnergyTmp=0.0 !Ejected electron energy
      eAngle =0.0;ds        =1   !Double stripping electron angle variable
      addElect=0 ;processE  =0   !Ejected electron integers
      eAngleSS=0.0;eEnergySS=0.0
      eAngleDS=0.0;eEnergyDS=0.0
      !*****************************
      call CollisionSim(E,tempQ,sigTot,process,excite,elect,disso,PID)
      collisions(PID(1)+1,PID(2)+1)=collisions(PID(1)+1,PID(2)+1)+1!Counting
7000 continue
      l=l+1
      if(l.ge.1000000)then
        !Filling ranVecA with a huge amount of numbers is a big time waster
        call ranlux(ranVecA,1000002) !Only get ranVec as needed
        l=1 !Reset l back to 1 (Start at 1 because ranVecA(l) is called next)
      end if
      !Calculate how far ion moves before a collision (dN)
      dN=-log(1-ranVecA(l))/sigTot(tempQold,nint(E))
      !Sometimes ranVecA is small enough to make DN 0
      if(dN.lt.1.0)goto 7000 !Get a new dN
      sigTotOld=sigTot(tempQold,nint(E))
      dNTot=dNTot+dN !Total change in column density
      do j=1,atmosLen !Loop through all of the atmosphere
        if(dNTot.le.totalCD(j+1))then !Move to proper CD of atmosphere
          !Calculate change in z based on the movement through the column dens.
          dZ=log((cos(pangle)*dN/(totalDens(dpt)*H(dpt)))+1)*H(dpt)
          dZTot=dZTot-dZ*1e-5 !Convert to km and keep subtracting from alt.
          do k=1,atmosLen
            if(dZTot.gt.altitude(k))then
              dNTot=totalCD(k)
              dpt=k !dpt is now the bin corresponding to depth
              if(dpt.gt.maxDpt) maxDpt=dpt !Used to see how deep we go
              goto 1000 !Get out of the do-loop that finds depth of penetration
            end if !Altitude if-statemet
          end do !Altitude do-loop
          !If we get here, then the ion has went through the entire atmosphere
          write(206,*)"JupOxyPrecip.f08: WARNING: Ion exited the bottom of the &
                     atmosphere, proceeding to next ion."
          goto 4000 !Continue on to the next ion
        end if !Column density if-statement
      end do !Column density do-loop
!*********************** Secondary Electron Calculations ***********************
      1000 continue
      if(PID(1).eq.4.or.PID(1).eq.5.or.PID(1).eq.7)then !Other processes
        addElect=addElect+11
      end if
      do j=1,elect !Do loop through all of the electrons ejected
        if(PID(1).eq.1.and.addElect.le.10)then !Single Ionization
          processE=1 !Process that goes into the EjectedElectron subroutine
          addElect=addElect+10 !Once addElect >10, go to look at PID(2)
        elseif(PID(1).eq.2.and.addElect.le.10)then !Double Ionization
          processE=2 !Process that goes into the EjectedElectron subroutine
          addElect=addElect+5 !Once addElect >10, go to look at PID(2)
        elseif(PID(1).eq.3.and.addElect.le.10)then !Transfer Ionization
          processE=3
          addElect=addElect+10
        elseif(PID(1).eq.6.and.addElect.le.10)then !DoubleCapture Autoionization
          processE=4
          addElect=addElect+10
        elseif(PID(2).eq.1.and.addElect.gt.10)then !Single Stripping
          processE=5
        elseif(PID(2).eq.2.and.addElect.gt.10)then !Double Stripping
          processE=6
        else
          write(206,*) "JupOxyPrecip.f08: WARNING: processE for &
          EjectedElectron.f08 not initialized."
          write(206,*) "JupOxyPrecip.f08: PID(1):",PID(1),"PID(2):",PID(2),&
                     "Electron:",j
        end if
        call EjectedElectron(E,processE,tempQold,eProbFunc,aProbFunc,&
                             eEnergyTmp,eAngle,bin)
        if(processE.eq.5)then
          eAngleSS=eAngle !Need the ejection angle for energy transformation
          eEnergySS=eEnergyTmp
        end if
        if(processE.eq.6)then
          eAngleDS(ds)=eAngle
          eEnergyDS(ds)=eEnergyTmp
          ds=ds+1
        end if
        totalElect=totalElect+1 !Total number of electrons produced
        eCounts(bin)=eCounts(bin)+1 !Total number of electrons produced vs. eng.
        eAngle=eAngle+(pangle*90/acos(0.0))
        !Must add the pitch angle to ejected e angle. pangle = [0,acos(0.0)]
        if(eAngle.le.90.0)then !Counting electrons going foward (downward)
          tElectFwd=tElectFwd+1 !Total electrons forward
          electFwdA(dpt)=electFwdA(dpt)+1 !Electrons forward vs. alt.
          electFwdAE(dpt,bin)=electFwdAE(dpt,bin)+1 !Elect fwd vs. alt. and eng.
        elseif(eAngle.le.270.0)then !Electrons going backward (0 is down)
          tElectBwd=tElectBwd+1 !Total electrons backward
          electBwdA(dpt)=electBwdA(dpt)+1 !Electrons backward vs. alt.
          electBwdAE(dpt,bin)=electBwdAE(dpt,bin)+1 !Elect bwd vs. alt. and eng.
        else !If the electron is ejected so far backward it's going fwd again
          write(206,*) "JupOxyPrecip.f08: WARNING: Elect ejection angle &
                        greater than 270 degrees."
        end if
        !Only want to add the electron energies for the NSIM process since SS
        !and DS have to be transformed into a different reference frame
        if(processE.le.4)eEnergy=eEnergy+eEnergyTmp
        addElect=addElect+1
      end do !j=1,elect
!************************** Counting X-Ray Production **************************
!* Note:
!*  An X-Ray count at a specific altitude and charge state means that there was
!*  an X-Ray producing collision at that specific altitude and the resultant ion
!*  was at the recorded charge state. That means, a collision that goes from
!*  O^8+ to O^7+ will be recorded as O^7+; therefore, the O^8+ bin will never
!*  produce an X-Ray. The last bin should ALWAYS be 0. Each processes can only
!*  create one X-Ray, if it's an X-Ray producing collision.
!*  Direct excitation (X-RayDE) producing collisions:
!*    TEX+SPEX(27) ,SI+SPEX(29), DI+SPEX(32)
!*  Charge exchange (X-RayCX) producing collisions:
!*    SC+SS(19), TI(25), SC(30)
!*  This is all done in the writing of the outputfiles. X-Ray productions are
!*  files 118 and 119 using the oxygen variable.
!*******************************************************************************
!********************** Counting Oxygen & H/H2 Production **********************
      oxygen(process,dpt,tempQ)=oxygen(process,dpt,tempQ)+1
      !If the process is SI, SC, or NEG, then there's a chance of dissociation
      if(disso.eq.2)then
        totHp(process,dpt)=totHp(process,dpt)+2
!*** Need to check which processes will ionize and dissociate
      elseif(disso.eq.1)then
        call ranlux(dissRan,1) !Random number to determine dissociation
        if(dissRan.le.0.1)then !10% chance of dissociation (H + H^+)
          totHp(process,dpt)=totHp(process,dpt)+1
        else !90% chance of no dissociation
          totH2p(process,dpt)=totH2p(process,dpt)+1
        end if
      elseif(disso.eq.0)then !TEX never dissociates, result is H2*
        H2Ex(dpt)=H2Ex(dpt)+1
      end if
!************************** Energy Loss Calculations ***************************
      call energyloss(E,tempQold,eEnergy,PID,dE,eAngleSS,eEnergySS,&
                      eAngleDS,eEnergyDS)
      dEsp=(dE)/dN !stopping power (calc before dE is recalculated)
      dEold=dE
      dE=(1/mass)*(1.0e-3)*(dE+stpnuc(E)*dN)*kappa !Total dE function
      if(dN.lt.0.0)then !Change in column density should never be less than 0
        write(206,10001) E,dEsp,dE,dN,dEold,process,PID(1),PID(2),tempQold
      end if
!********************** Oxygen Charge State Distribution ***********************
      do j=1,nOxEngBins
        if(E.le.oxEngBins(j))then
          OxyVsEng(tempQ,j)=OxyVsEng(tempQ,j)+1
          goto 2000
        end if
      end do
2000 continue
      do j=1,nStopPowerEBins
        if(E.le.stopPowerEBins(j))then
          SPvsEng(j)=SPvsEng(j)+dEsp
          SigTotvsEng(j)=SigTotvsEng(j)+SigTotOld
          dEvsEng(j)=dEvsEng(j)+dEold !Times 16 to get rid of eV/u
          dNvsEng(j)=dNvsEng(j)+dN
          ProcessdE(j,processC(process),tempQold)=&
            ProcessdE(j,processC(process),tempQold)+dEold
          nSPions(j)=nSPions(j)+1
          pnSPions(j,processC(process),tempQold)=&
            pnSPions(j,processC(process),tempQold)+1
          goto 3000
        end if
      end do
3000 continue
      E=E-dE !Find the new energy !~
      tempQold=tempQ !Assign newly acquired charge state to old variable
      if(E.lt.1.0) goto 4000 !Stop once the energy is less than 1 keV/u
      if(i.eq.numSim)then !~
        write(206,*) 'JupOxyPrecip.f08: ERROR: numSim not large enough.' !~
        write(206,*) 'JupOxyPrecip.f08: Ion energy was: ',E
        goto 4000
      end if
    end do !i=1,numSim
4000 continue !Continue to go on to the next ion
  end do !ion=1,number_of_ions
!****************************** Calculate Oxy CX *******************************
  do i=1,atmosLen
    do j=1,nChS
      do k=1,14
        oxygenCX(i,j)=oxygenCX(i,j)+oxygen(CXProc(k),i,j)
      end do
    end do
  end do
!******************************** Output Header ********************************
  write(206,*)"------------------------------------------NEW RUN---------------&
              ---------------------------"
  write(206,*)"Number of ions: ", number_of_ions
  write(206,*)"Initial energy: ", energy, 'keV'
  write(206,*)"Trial number:   ", trial
  write(206,*)"****************************************************************&
              ***************************"
  !******* Check various electron counters
  write(206,*)'Sum of total electrons foward:         ',tElectFwd,sum(electFwdA)
  write(206,*)'Sum of total electrons backward:       ',tElectBwd,sum(electBwdA)
  write(206,*)'Sum of total electrons foward+backward:',sum(electFwdA+electBwdA)
  write(206,*)'Sum of total electrons:                ',totalElect,sum(eCounts)
  write(206,*)'Max Depth:                             ',altitude(maxDpt)
  flush(206)
!********** Open output data files for each set of initial energies ************
  do i=1,nOutputFiles
    write(filename,'("./Output/",I0,"keV/",A,I0,".dat")') &
          energy,trim(filenames(i)),trial
    filename=trim(filename)
    open(unit=100+i,file=filename,status='unknown')
  end do
!***************************** Write out to files ******************************
  norm=number_of_ions*2e5 !Normalization condition to per ion per km
  totalHp=sum(totHp,dim=1)
  totalH2p=sum(totH2p,dim=1)
  pHp=sum(totHp,dim=2)
  pH2p=sum(totH2p,dim=2)
  l=0
  k=0
  do i=1,nProc
    if(pHp(i).gt.0)then !Pick out which processes cause H+
      l=l+1
      npHp(l)=i
    end if
    if(pH2p(i).gt.0)then !Pick out which processes cause H2+
      k=k+1
      npH2p(k)=i
    end if
  end do
  write(101,H01) !H+ header
  write(101,*) "Alt [km] ", (HProc(npHp(i)),i=1,l), "Total"
  write(102,H02) !H2+ header
  write(102,*) "Alt [km] ", (HProc(npH2P(i)),i=1,k), "Total"
  write(103,H03) !H2* header
  do i=1,atmosLen !Ionization/Excitation vs. altitude
    write(101,F01) altitude(i),(real(totHp(npHp(j),i))/norm,j=1,l),&
                   real(totalHp(i))/norm
    write(102,F01) altitude(i),(real(totH2p(npH2p(j),i))/norm,j=1,k),&
                   real(totalH2p(i))/norm
    write(103,F02) altitude(i),real(H2Ex(i))/norm
  end do
  totO=sum(OxyVsEng,dim=1)
  write(104,H04) !Oxy vs energy header
  do i=1,nOxEngBins !Oxygen charge state distribution
    write(104,F03) oxEngBins(i)-(oxEngBinSize/2.0), &
                   (real(OxyVsEng(j,i))/real(totO(i)),j=1,nChS)
  end do
  write(105,H05) !Stopping power header
  do i=1,nStopPowerEBins !Stopping power vs. ion energy
    write(105,F04) stopPowerEBins(i)-(delSP(i)/2.0), &
                   SPvsEng(i)/real(nSPions(i)**2), &
                   SigTotvsEng(i)/real(nSPions(i)), &
                   dEvsEng(i)/real(nSPions(i)), &
                   dNvsEng(i)/real(nSPions(i)), &
                   (SigTotvsEng(i)*dEvsEng(i))/(real(nSPions(i))**2), &
                   nSPions(i)
  end do
  NSIM=sum(collisions,dim=2)
  SIM=sum(collisions,dim=1)
  write(106,H07) !Collisions header
  do i=1,8 !Total number of each type of collision
    write(106,F06) Coll(i), (collisions(i,j),j=1,5),NSIM(i)
  end do
  write(106,*) '------------------------------------------------------------&
                --------------------------'
  write(106,F06) 'Sum     ',SIM(1),SIM(2),SIM(3),SIM(4),SIM(5)
  write(106,*) ''
  write(106,H07) !Collisions percentage header
  do i=1,8 !Total percentage of each type of collision
    write(106,F07) Coll(i), &
    (real(collisions(i,j))/real(sum(collisions))*100.0,j=1,5),&
    real(NSIM(i))/real(sum(NSIM))*100
  end do
  write(106,*) '------------------------------------------------------------&
                --------------------------'
  write(106,F07) 'Sum     ',real(SIM(1))/real(sum(SIM))*100,&
   real(SIM(2))/real(sum(SIM))*100,real(SIM(3))/real(sum(SIM))*100,&
   real(SIM(4))/real(sum(SIM))*100,real(SIM(5))/real(sum(SIM))*100
  do i=1,nChS !Oxygen production
    write(106+i,*) "Alt [km] ", (HProc(k),k=1,nProc)
    do j=1,atmosLen
      write(106+i,F01) altitude(j),(real(oxygen(k,j,i))/norm,k=1,nProc)
    end do
  end do
  do i=1,3
    write(116+i,H06)
  end do
  do i=1,atmosLen !Oxygen production from charge exchange
    write(117,F05) altitude(i),(real(oxygenCX(i,j))/norm,j=1,nChS)
  end do
  do i=1,atmosLen !DE - TEX+SPEX,SI+SPEX,DI+SPEX, CX - SC+SS,TI,SC
    write(118,F05) altitude(i),& !X-Ray production from direct excitation
      (real(oxygen(27,i,j)+oxygen(29,i,j)+oxygen(32,i,j))/norm,j=1,nChs)
    write(119,F05) altitude(i),& !X-Ray production from charge exchange
      (real(oxygen(19,i,j)+oxygen(25,i,j)+oxygen(30,i,j))/norm,j=1,nChs)
  end do
!***************************** Secondary Electrons *****************************
  do i=1,atmosLen
    do j=1,nE2strBins
      prode2stF(i,j)=real(electFwdAE(i,j))/norm
      prode2stB(i,j)=real(electBwdAE(i,j))/norm
    end do
  end do
  do j=1,nE2strBins !2-Stream electrons, forward and backward
    write(120,F2Str) (prode2stF(i,j),i=atmosLen,1,-1)
    write(121,F2Str) (prode2stB(i,j),i=atmosLen,1,-1)
  end do
!**************** Close all of the files that have been opened *****************
  do i=1,nOutputFiles
    close(100+i)
  end do
!****************** Write to screen some general information *******************
  write(206,*) '--------------------------------------------------'
  write(206,*) ' (NEG)  =',NSIM(1)
  write(206,*) ' (SI)   =',NSIM(2)
  write(206,*) ' (DI)   =',NSIM(3)
  write(206,*) ' (TI)   =',NSIM(4)
  write(206,*) ' (SC)   =',NSIM(5)
  write(206,*) ' (DC)   =',NSIM(6)
  write(206,*) ' (DCAI) =',NSIM(7)
  write(206,*) ' (TEX)  =',NSIM(8)
  write(206,*) ' NSIM only processes: ',SIM(1)
  write(206,*) '--------------------------------------------------'
  write(206,*) ' (SS)   =',SIM(2)
  write(206,*) ' (DS)   =',SIM(3)
  write(206,*) ' (SPEX) =',SIM(4)
  write(206,*) ' (DPEX) =',SIM(5)
  write(206,*) ' SIM processes:      ',sum(SIM)-SIM(1)
  write(206,*) '--------------------------------------------------'
  write(206,*) ' Total collisions:   ',sum(SIM)

  call system_clock(t4,clock_rate,clock_max) !Elapsed time for a single energy
  hrs=int(real(t4-t3)/clock_rate/3600.0)
  min=int(((real(t4-t3)/clock_rate)-hrs*3600)/60)
  sec=mod(real(t4-t3)/clock_rate,60.0)
  write(206,*) 'Individual run elapsed real time = ',hrs,':',min,':',sec
  deallocate(angle)
end do !run=1,number_of_energies
call system_clock (t2,clock_rateTotal,clock_maxTotal) !Total elapsed time
write(206,*) 'Total elapsed real time = ', real(t2-t1)/clock_rateTotal
close(206)
write(filename,"('./Output/',I0,'keV/Elapsed_Times.dat')") energy
1003 continue
inquire(file=filename,opened=open)
if(open)goto 1003
open(unit=207,file=filename,status='unknown',access='append',action='write')
write(207,*) trial,hrs,':',min,':',sec
close(207)
!**************************** Formatting conditions ****************************
!* Located in formatting.08. First compile with:
!*   gfortran -c formatting.08
!* this creates formatting.o, which then needs to be included in the command to
!* run this code.
!*******************************************************************************
10001 format(F8.2,3(2x,ES10.3E2),2x,F9.2,4(2x,I2))
end program
