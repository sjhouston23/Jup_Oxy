subroutine CollisionSim(Energy,tempQ,sigTot,process,excite,elect,disso,PID)
!*******************************************************************************
!* Created by Stephen J. Houston 2.20.18
!*******************************************************************************
!* This subroutine reads in each individual cross-section and the total cross-
!* section for all the processes for each energy and charge state for 35
!* different collision types (see processes below). Then uses a random
!* number from 0-1 to and compare it to the probability for each
!* collision process (individual cross-section divided by total cross-section).
!* This probability will be calcultated from the cross-sections for the
!* different collision processes. The probability that is larger than the random
!* number will be the selected outcome process of the collision.
!*******************************************************************************
!*
!* 	Input:
!*		Energy --> Energy of the ion
!*		 Type: Real
!*		 Units: keV/u
!*
!*		tempQ --> Charge state
!*		 Type: Integer
!*		 Units: None
!*
!*   Returns:
!*    tempQ --> Replaced with new charge state
!*		 Type: Integer
!*		 Units: None
!*
!*    sigTot --> Total cross-section
!*		 Type: Real
!*		 Units: cm^-2
!*
!*		process --> Type of collision that occured (ordered by # of charge states)
!*		 Type: Integer
!*		 Units: None
!*       1  = DC+DS
!*       2  = DC+DPEX
!*       3  = DCAI+DS
!*       4  = DCAI+DPEX
!*       5  = TI+DS
!*       6  = TI+DPEX
!*       7  = SC+DS
!*       8  = SC+DPEX
!*       9  = DC+SS
!*       10 = DC+SPEX
!*       11 = DCAI+SS
!*       12 = DCAI+SPEX
!*       13 = TI+SS
!*       14 = TI+SPEX
!*       15 = TEX+DS
!*       16 = TEX+DPEX
!*       17 = SI+DS
!*       18 = SI+DPEX
!*       19 = SC+SS
!*       20 = SC+SPEX
!*       21 = DI+DS
!*       22 = DI+DPEX
!*       23 = DCAI
!*       24 = DC
!*       25 = TI
!*       26 = TEX+SS
!*       27 = TEX+SPEX
!*       28 = SI+SS
!*       29 = SI+SPEX
!*       30 = SC
!*       31 = DI+SS
!*       32 = DI+SPEX
!*       33 = TEX
!*       34 = SI
!*       35 = DI
!*
!*    excite --> Number of excitations of target (H2)
!*     Type: Integer
!*     Units: None
!*
!*    elect --> Number of Electrons Ejected
!*     Type: Integer
!*     Units: None
!*
!*    disso --> Number indicating whether dissociation is possible (0-2)
!*     Type: Integer
!*     Units: None
!*
!*    PID --> Collision Type (1-7,0-4)
!*     Type: 2-Component Array, Integer
!*     Units: None
!*
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************

integer, intent(inout) :: tempQ
integer, intent(out) :: process,excite,elect,disso
integer, intent(out) :: PID(2)
!Process ID (first number is NSIM, second is SIM process)
real*8, intent(in) :: Energy
real*8, intent(out) :: sigTot(10,25000)

integer ChS1,ChS5,ChS6,ChS7,ChS8,ChS9,ChS10 !Number of charge states
integer nEng  !Number of energies.
integer nFiles !Number of cross-section files.

parameter (ChS1=1,ChS5=5,ChS6=6,ChS7=7,ChS8=8,ChS9=9,ChS10=10)
parameter (nEng=25000,nFiles=37)

!integer tempQ
!integer process,excite
!real*8 Energy
!real*8 sigTot(ChS9,nEng)

integer t1,t2,clock_maxTotal,clock_rateTotal !Used to calculate comp. time
integer error,errorP(nFiles),E !errorP - processes error
integer Q !Used to calculated number of electrons produced (tempQ-Q)

character(len=100) filename !Will be open() filename.
character(len=20) XSFile(nFiles) !Cross-section filenames

real*8,dimension(ChS1,nEng),save :: neg
real*8,dimension(ChS5,nEng),save :: dcds,dcdpex,dcaids,dcaidpex
real*8,dimension(ChS6,nEng),save :: tids,tidpex,scds,scdpex,dcss,dcspex,dcaiss,&
                                    dcaispex
real*8,dimension(ChS7,nEng),save :: tiss,tispex,texds,texdpex,sids,sidpex,scss,&
                                    scspex,dids,didpex,dcai,dc
real*8,dimension(ChS8,nEng),save :: ti,texspex,sispex,sc,dispex
real*8,dimension(ChS9,nEng),save :: texss,siss,diss
real*8,dimension(ChS10,nEng),save:: tex,si,di
real*8,save :: dum
real*8 prob(36) !dummy variable, energy, and probability for each coll.
real*8 sumProb !The sum of probabilities

real ranVecB(1)

external ranlux
!****************************** Data Declaration *******************************

data XSFile /'dcds','dcdpex','dcaids','dcaidpex','tids','tidpex','scds',&
             'scdpex','dcss','dcspex','dcaiss','dcaispex','tiss','tispex',&
             'texds','texdpex','sids','sidpex','scss','scspex','dids','didpex',&
             'dcai','dc','ti','texss','texspex','siss','sispex','sc','diss',&
             'dispex','tex','si','di','neg','TotalXS'/

!******************************** Main Program *********************************
!Calculate the total computational run time of the collision read in:
if(dum.gt.24999)goto 2000
call system_clock (t1,clock_rateTotal,clock_maxTotal)
write(206,*) 'CollisionSim.f08: Opening collision process files...'
do i=1,nFiles
  write(filename,'("./SIMXSInterp/",A,".dat")') trim(XSfile(i))
  open(unit=99+i,file=filename,status='old',iostat=error) !Open all files
  if(error.gt.0)then !If there is an error opening one of the files
    write(206,*) 'CollisionSim.f08: Error opening file: ',filename
    STOP 'CollisionSim.f08: Stopping program...'
  end if
end do
write(206,*) 'CollisionSim.f08: Reading collision process files...'
do i=1,25000 !Read in all the cross-sections including the total cross-section
  read(100,*,iostat=errorP(1))  dum, (dcds(j,i),j=1,ChS5)
  read(101,*,iostat=errorP(2))  dum, (dcdpex(j,i),j=1,ChS5)
  read(102,*,iostat=errorP(3))  dum, (dcaids(j,i),j=1,ChS5)
  read(103,*,iostat=errorP(4))  dum, (dcaidpex(j,i),j=1,ChS5)
  read(104,*,iostat=errorP(5))  dum, (tids(j,i),j=1,ChS6)
  read(105,*,iostat=errorP(6))  dum, (tidpex(j,i),j=1,ChS6)
  read(106,*,iostat=errorP(7))  dum, (scds(j,i),j=1,ChS6)
  read(107,*,iostat=errorP(8))  dum, (scdpex(j,i),j=1,ChS6)
  read(108,*,iostat=errorP(9))  dum, (dcss(j,i),j=1,ChS6)
  read(109,*,iostat=errorP(10)) dum, (dcspex(j,i),j=1,ChS6)
  read(110,*,iostat=errorP(11)) dum, (dcaiss(j,i),j=1,ChS6)
  read(111,*,iostat=errorP(12)) dum, (dcaispex(j,i),j=1,ChS6)
  read(112,*,iostat=errorP(13)) dum, (tiss(j,i),j=1,ChS7)
  read(113,*,iostat=errorP(14)) dum, (tispex(j,i),j=1,ChS7)
  read(114,*,iostat=errorP(15)) dum, (texds(j,i),j=1,ChS7)
  read(115,*,iostat=errorP(16)) dum, (texdpex(j,i),j=1,ChS7)
  read(116,*,iostat=errorP(17)) dum, (sids(j,i),j=1,ChS7)
  read(117,*,iostat=errorP(18)) dum, (sidpex(j,i),j=1,ChS7)
  read(118,*,iostat=errorP(19)) dum, (scss(j,i),j=1,ChS7)
  read(119,*,iostat=errorP(20)) dum, (scspex(j,i),j=1,ChS7)
  read(120,*,iostat=errorP(21)) dum, (dids(j,i),j=1,ChS7)
  read(121,*,iostat=errorP(22)) dum, (didpex(j,i),j=1,ChS7)
  read(122,*,iostat=errorP(23)) dum, (dcai(j,i),j=1,ChS7)
  read(123,*,iostat=errorP(24)) dum, (dc(j,i),j=1,ChS7)
  read(124,*,iostat=errorP(25)) dum, (ti(j,i),j=1,ChS8)
  read(125,*,iostat=errorP(26)) dum, (texss(j,i),j=1,ChS9)
  read(126,*,iostat=errorP(27)) dum, (texspex(j,i),j=1,ChS8)
  read(127,*,iostat=errorP(28)) dum, (siss(j,i),j=1,ChS9)
  read(128,*,iostat=errorP(29)) dum, (sispex(j,i),j=1,ChS8)
  read(129,*,iostat=errorP(30)) dum, (sc(j,i),j=1,ChS8)
  read(130,*,iostat=errorP(31)) dum, (diss(j,i),j=1,ChS9)
  read(131,*,iostat=errorP(32)) dum, (dispex(j,i),j=1,ChS8)
  read(132,*,iostat=errorP(33)) dum, (tex(j,i),j=1,ChS10)
  read(133,*,iostat=errorP(34)) dum, (si(j,i),j=1,ChS10)
  read(134,*,iostat=errorP(35)) dum, (di(j,i),j=1,ChS10)
  read(135,*,iostat=errorP(36)) dum, (neg(j,i),j=1,ChS1)
  read(136,*,iostat=errorP(37)) dum, (sigTot(j,i),j=1,ChS10)
  !Make sure all the files are read correctly and if not, where was the problem
  do j=1,nFiles
    if(errorP(j).gt.0)then
      write(206,'(A,I3.2,A,I6)') 'CollisionSim.f08: Error reading file&
                               corresponding to no.',j,' at an energy of',i
      STOP 'CollisionSim.f08: Stopping program...'
    end if
  end do
end do
do i=1,nFiles !Close all the files we've opened
  close(99+i)
end do
call system_clock (t2,clock_rateTotal,clock_maxTotal)
write(206,*) 'CollisionSim.f08: Collision process files opened, read, and closed &
            in', real(t2-t1)/clock_rateTotal, 'seconds.'
2000 continue
!*******************************************************************************
!******************** Collision-Type Probability Calculation *******************
!*******************************************************************************
!* Calculate the transition probabilities by taking the collision process XS and
!* dividing it by the total cross-section. Since each process has different
!* charge state range the IQQ must match the actual charge state desired!
!* Must skip all of the processes that cannot occur for a given charge state.
!*******************************************************************************
!call system_clock (t1,clock_rateTotal,clock_maxTotal)
!write(*,*) 'CollisionSim.f08: Calculating process probabilities...'
!Energy=2000.0
E=nint(Energy)
!tempQ=1
prob=0.0
if(tempQ.eq.1)then !**************************** O- ****************************
  prob(26)=texss(tempQ,E)/sigTot(tempQ,E)
  prob(28)=siss(tempQ,E)/sigTot(tempQ,E)
  prob(31)=diss(tempQ,E)/sigTot(tempQ,E)
  prob(33)=tex(tempQ,E)/sigTot(tempQ,E)
  prob(34)=si(tempQ,E)/sigTot(tempQ,E)
  prob(35)=di(tempQ,E)/sigTot(tempQ,E)
elseif(tempQ.eq.2)then !************************** O ***************************
  prob(15)=texds(tempQ-1,E)/sigTot(tempQ,E)
  prob(16)=texdpex(tempQ-1,E)/sigTot(tempQ,E)
  prob(17)=sids(tempQ-1,E)/sigTot(tempQ,E)
  prob(18)=sidpex(tempQ-1,E)/sigTot(tempQ,E)
  prob(21)=dids(tempQ-1,E)/sigTot(tempQ,E)
  prob(22)=didpex(tempQ-1,E)/sigTot(tempQ,E)
  prob(26)=texss(tempQ,E)/sigTot(tempQ,E)
  prob(27)=texspex(tempQ-1,E)/sigTot(tempQ,E)
  prob(28)=siss(tempQ,E)/sigTot(tempQ,E)
  prob(29)=sispex(tempQ-1,E)/sigTot(tempQ,E)
  prob(31)=diss(tempQ,E)/sigTot(tempQ,E)
  prob(32)=dispex(tempQ-1,E)/sigTot(tempQ,E)
  prob(33)=tex(tempQ,E)/sigTot(tempQ,E)
  prob(34)=si(tempQ,E)/sigTot(tempQ,E)
  prob(35)=di(tempQ,E)/sigTot(tempQ,E)
  prob(36)=neg(tempQ-1,E)/sigTot(tempQ,E)
elseif(tempQ.eq.3)then !************************** O+ **************************
  prob(5) =tids(tempQ-2,E)/sigTot(tempQ,E)
  prob(6) =tidpex(tempQ-2,E)/sigTot(tempQ,E)
  prob(7) =scds(tempQ-2,E)/sigTot(tempQ,E)
  prob(8) =scdpex(tempQ-2,E)/sigTot(tempQ,E)
  prob(13)=tiss(tempQ-2,E)/sigTot(tempQ,E)
  prob(14)=tispex(tempQ-2,E)/sigTot(tempQ,E)
  prob(15)=texds(tempQ-1,E)/sigTot(tempQ,E)
  prob(16)=texdpex(tempQ-1,E)/sigTot(tempQ,E)
  prob(17)=sids(tempQ-1,E)/sigTot(tempQ,E)
  prob(18)=sidpex(tempQ-1,E)/sigTot(tempQ,E)
  prob(19)=scss(tempQ-2,E)/sigTot(tempQ,E)
  prob(20)=scspex(tempQ-2,E)/sigTot(tempQ,E)
  prob(21)=dids(tempQ-1,E)/sigTot(tempQ,E)
  prob(22)=didpex(tempQ-1,E)/sigTot(tempQ,E)
  prob(25)=ti(tempQ-2,E)/sigTot(tempQ,E)
  prob(26)=texss(tempQ,E)/sigTot(tempQ,E)
  prob(27)=texspex(tempQ-1,E)/sigTot(tempQ,E)
  prob(28)=siss(tempQ,E)/sigTot(tempQ,E)
  prob(29)=sispex(tempQ-1,E)/sigTot(tempQ,E)
  prob(30)=sc(tempQ-2,E)/sigTot(tempQ,E)
  prob(31)=diss(tempQ,E)/sigTot(tempQ,E)
  prob(32)=dispex(tempQ-1,E)/sigTot(tempQ,E)
  prob(33)=tex(tempQ,E)/sigTot(tempQ,E)
  prob(34)=si(tempQ,E)/sigTot(tempQ,E)
  prob(35)=di(tempQ,E)/sigTot(tempQ,E)
elseif(tempQ.ge.4.and.tempQ.le.8)then !*************** O++ - O6+ ***************
  prob(1) =dcds(tempQ-3,E)/sigTot(tempQ,E)
  prob(2) =dcdpex(tempQ-3,E)/sigTot(tempQ,E)
  prob(3) =dcaids(tempQ-3,E)/sigTot(tempQ,E)
  prob(4) =dcaidpex(tempQ-3,E)/sigTot(tempQ,E)
  prob(5) =tids(tempQ-2,E)/sigTot(tempQ,E)
  prob(6) =tidpex(tempQ-2,E)/sigTot(tempQ,E)
  prob(7) =scds(tempQ-2,E)/sigTot(tempQ,E)
  prob(8) =scdpex(tempQ-2,E)/sigTot(tempQ,E)
  prob(9) =dcss(tempQ-3,E)/sigTot(tempQ,E)
  prob(10)=dcspex(tempQ-3,E)/sigTot(tempQ,E)
  prob(11)=dcaiss(tempQ-3,E)/sigTot(tempQ,E)
  prob(12)=dcaispex(tempQ-3,E)/sigTot(tempQ,E)
  prob(13)=tiss(tempQ-2,E)/sigTot(tempQ,E)
  prob(14)=tispex(tempQ-2,E)/sigTot(tempQ,E)
  prob(15)=texds(tempQ-1,E)/sigTot(tempQ,E)
  prob(16)=texdpex(tempQ-1,E)/sigTot(tempQ,E)
  prob(17)=sids(tempQ-1,E)/sigTot(tempQ,E)
  prob(18)=sidpex(tempQ-1,E)/sigTot(tempQ,E)
  prob(19)=scss(tempQ-2,E)/sigTot(tempQ,E)
  prob(20)=scspex(tempQ-2,E)/sigTot(tempQ,E)
  prob(21)=dids(tempQ-1,E)/sigTot(tempQ,E)
  prob(22)=didpex(tempQ-1,E)/sigTot(tempQ,E)
  prob(23)=dcai(tempQ-3,E)/sigTot(tempQ,E)
  prob(24)=dc(tempQ-3,E)/sigTot(tempQ,E)
  prob(25)=ti(tempQ-2,E)/sigTot(tempQ,E)
  prob(26)=texss(tempQ,E)/sigTot(tempQ,E)
  prob(27)=texspex(tempQ-1,E)/sigTot(tempQ,E)
  prob(28)=siss(tempQ,E)/sigTot(tempQ,E)
  prob(29)=sispex(tempQ-1,E)/sigTot(tempQ,E)
  prob(30)=sc(tempQ-2,E)/sigTot(tempQ,E)
  prob(31)=diss(tempQ,E)/sigTot(tempQ,E)
  prob(32)=dispex(tempQ-1,E)/sigTot(tempQ,E)
  prob(33)=tex(tempQ,E)/sigTot(tempQ,E)
  prob(34)=si(tempQ,E)/sigTot(tempQ,E)
  prob(35)=di(tempQ,E)/sigTot(tempQ,E)
elseif(tempQ.eq.9)then !************************* O7+ **************************
  prob(9) =dcss(tempQ-3,E)/sigTot(tempQ,E)
  prob(10)=dcspex(tempQ-3,E)/sigTot(tempQ,E)
  prob(11)=dcaiss(tempQ-3,E)/sigTot(tempQ,E)
  prob(12)=dcaispex(tempQ-3,E)/sigTot(tempQ,E)
  prob(13)=tiss(tempQ-2,E)/sigTot(tempQ,E)
  prob(14)=tispex(tempQ-2,E)/sigTot(tempQ,E)
  prob(19)=scss(tempQ-2,E)/sigTot(tempQ,E)
  prob(20)=scspex(tempQ-2,E)/sigTot(tempQ,E)
  prob(23)=dcai(tempQ-3,E)/sigTot(tempQ,E)
  prob(24)=dc(tempQ-3,E)/sigTot(tempQ,E)
  prob(25)=ti(tempQ-2,E)/sigTot(tempQ,E)
  prob(26)=texss(tempQ,E)/sigTot(tempQ,E)
  prob(27)=texspex(tempQ-1,E)/sigTot(tempQ,E)
  prob(28)=siss(tempQ,E)/sigTot(tempQ,E)
  prob(29)=sispex(tempQ-1,E)/sigTot(tempQ,E)
  prob(30)=sc(tempQ-2,E)/sigTot(tempQ,E)
  prob(31)=diss(tempQ,E)/sigTot(tempQ,E)
  prob(32)=dispex(tempQ-1,E)/sigTot(tempQ,E)
  prob(33)=tex(tempQ,E)/sigTot(tempQ,E)
  prob(34)=si(tempQ,E)/sigTot(tempQ,E)
  prob(35)=di(tempQ,E)/sigTot(tempQ,E)
elseif(tempQ.eq.10)then !************************* O8+ **************************
  prob(23)=dcai(tempQ-3,E)/sigTot(tempQ,E)
  prob(24)=dc(tempQ-3,E)/sigTot(tempQ,E)
  prob(25)=ti(tempQ-2,E)/sigTot(tempQ,E)
  prob(30)=sc(tempQ-2,E)/sigTot(tempQ,E)
  prob(33)=tex(tempQ,E)/sigTot(tempQ,E)
  prob(34)=si(tempQ,E)/sigTot(tempQ,E)
  prob(35)=di(tempQ,E)/sigTot(tempQ,E)
else !Error readout
  write(206,*) "CollisionSim.f08: WARNING: Collision probability is uninitialized"
  write(206,*) "CollisionSim.f08: The charge state (tempQ) does not lie between &
              1 and 9."
  STOP "CollisionSim.f08: Stopping program..."
endif
if(sum(prob).ge.1.01.or.sum(prob).le.0.9999)then
  write(206,*) "CollisionSim.f08: WARNING: Normalized collision probability is &
              not close enough to 1. The value is: ", sum(prob)
  write(206,*) "tempQ: ", tempQ, "Energy: ", E
end if
!*******************************************************************************
!************************* Collision-Type Determination ************************
!*******************************************************************************
!* Now that the collision probability is determined, a simple monte carlo
!* technique is used to determine the collision type. After collision-type is
!* determined the final charge state and collision-type are recorded as
!* tempQ and process, respectively.
!*
!* Algorithm:
!*  Step 1) Generate a random number on the interval [0,1]
!*  Step 2) Parce a sub-interval to cover the probability of a single
!*          collision type; i.e. for Single Ionization from 0->P(SI) [P(SI)<1]
!*  Step 3) Ask if the randomly generated number is in the interval [0,P(SI)],
!*          if so write out the final charge state and process
!*  Step 4) if not, parce the next sub-interval of [0,1] from [P(SI),P(DI)]
!*  Step 5) Repeat step 3.
!*  Step 6) Repreat step 4&5 through all 35 collision types until collision is
!*          determined.
!*
!* PID(1) NSIM Processes             Abbreviation  Charge Change  Oxygen Excite
!*   0    Negative Ion Channel           (NEG)          +1              0
!*   1    Single Ionization              (SI)            0              0
!*   2    Double Ionization              (DI)            0              0
!*   3    Transfer Ionization            (TI)           -1              0
!*   4    Single Capture                 (SC)           -1              1
!*   5    Double Capture                 (DC)           -2              2
!*   6    Double-Capture Autoionization  (DCAI)         -1            2->0
!*   7    Target Excitation              (TEX)           0              0
!* PID(2) SIM Processes
!*   1    Single Stripping               (SS)           +1              0
!*   2    Double Stripping               (DS)           +2              0
!*   3    Single Projectile Excitation   (SPEX)          0              1
!*   4    Double Projectile Excitation   (DPEX)          0              2
!*
!* Some processes will dissociate H2 always, sometimes, or never. This is
!* followed with the "disso" variable.
!* disso = 0, never dissociates
!* disso = 1, 10% chance of dissociation (determined in the main program)
!* disso = 2, always dissociates
!*
!* PID is used so the same NSIM and SIM processes can be categorized for easier
!* use later on in the main program. If an SIM process is 0, then that means the
!* process that occured was solely NSIM.
!*
!*******************************************************************************
call ranlux(ranVecB,1) !Only need 1 random number every time there's a collision
sumProb=0
elect=0
Q=tempQ
if(ranVecB(1).gt.0.99999)ranVecB(1)=ranVecB(1)-0.00001
!write(*,*) ranVecB(1),SigTot(tempQ,E)
!Go through the probabilities until the random number falls within a range
do i=1,36
  sumProb=sumProb+prob(i)
  if(ranVecB(1).le.sumProb)then
    process=i !Assign the process variable the correct index
    !DC+DS
    if(i.eq.1)then
      tempQ=tempQ !The charge state may change if an electron is gained or lost
      excite=2 !The number of excited electrons bound to the oxygen ion
      elect=2 !How many electrons are produced
      disso=2 !Is there dissociation
      PID(1)=5
      PID(2)=2
    !DC+DPEX
    elseif(i.eq.2)then
      tempQ=tempQ-2
      excite=4
      elect=0
      disso=2
      PID(1)=5
      PID(2)=4
    !DCAI+DS
    elseif(i.eq.3)then
      tempQ=tempQ+1
      excite=0
      elect=3
      disso=2
      PID(1)=6
      PID(2)=2
    !DCAI+DPEX
    elseif(i.eq.4)then
      tempQ=tempQ-1
      excite=2
      elect=1
      disso=2
      PID(1)=6
      PID(2)=4
    !TI+DS
    elseif(i.eq.5)then
      tempQ=tempQ+1
      excite=0
      elect=3
      disso=2
      PID(1)=3
      PID(2)=2
    !TI+DPEX
    elseif(i.eq.6)then
      tempQ=tempQ-1
      excite=2
      elect=1
      disso=2
      PID(1)=3
      PID(2)=4
    !SC+DS
    elseif(i.eq.7)then
      tempQ=tempQ+1
      excite=1
      elect=2
      disso=1
      PID(1)=4
      PID(2)=2
    !SC+DPEX
    elseif(i.eq.8)then
      tempQ=tempQ-1
      excite=3
      elect=0
      disso=1
      PID(1)=4
      PID(2)=4
    !DC+SS
    elseif(i.eq.9)then
      tempQ=tempQ-1
      excite=2
      elect=1
      disso=2
      PID(1)=5
      PID(2)=1
    !DC+SPEX
    elseif(i.eq.10)then
      tempQ=tempQ-2
      excite=3
      elect=0
      disso=2
      PID(1)=5
      PID(2)=3
    !DCAI+SS
    elseif(i.eq.11)then
      tempQ=tempQ
      excite=0
      elect=2
      disso=2
      PID(1)=6
      PID(2)=1
    !DCAI+SPEX
    elseif(i.eq.12)then
      tempQ=tempQ-1
      excite=1
      elect=1
      disso=2
      PID(1)=6
      PID(2)=3
    !TI+SS
    elseif(i.eq.13)then
      tempQ=tempQ
      excite=0
      elect=2
      disso=2
      PID(1)=3
      PID(2)=1
    !TI+SPEX
    elseif(i.eq.14)then
      tempQ=tempQ-1
      excite=1
      elect=1
      disso=2
      PID(1)=3
      PID(2)=3
    !TEX+DS
    elseif(i.eq.15)then
      tempQ=tempQ+2
      excite=0
      elect=2
      disso=0
      PID(1)=7
      PID(2)=2
    !TEX+DPEX
    elseif(i.eq.16)then
      tempQ=tempQ
      excite=2
      elect=0
      disso=0
      PID(1)=7
      PID(2)=4
    !SI+DS
    elseif(i.eq.17)then
      tempQ=tempQ+2
      excite=0
      elect=3
      disso=1
      PID(1)=1
      PID(2)=2
    !SI+DPEX
    elseif(i.eq.18)then
      tempQ=tempQ
      excite=2
      elect=1
      disso=1
      PID(1)=1
      PID(2)=4
    !SC+SS
    elseif(i.eq.19)then
      tempQ=tempQ
      excite=1
      elect=1
      disso=1
      PID(1)=4
      PID(2)=1
    !SC+SPEX
    elseif(i.eq.20)then
      tempQ=tempQ-1
      excite=2
      elect=0
      disso=1
      PID(1)=4
      PID(2)=3
    !DI+DS
    elseif(i.eq.21)then
      tempQ=tempQ+2
      excite=0
      elect=4
      disso=2
      PID(1)=2
      PID(2)=2
    !DI+DPEX
    elseif(i.eq.22)then
      tempQ=tempQ
      excite=2
      elect=2
      disso=2
      PID(1)=2
      PID(2)=4
    !DCAI
    elseif(i.eq.23)then
      tempQ=tempQ-1
      excite=0
      elect=1
      disso=2
      PID(1)=6
      PID(2)=0
    !DC
    elseif(i.eq.24)then
      tempQ=tempQ-2
      excite=2
      elect=0
      disso=2
      PID(1)=5
      PID(2)=0
    !TI
    elseif(i.eq.25)then
      tempQ=tempQ-1
      excite=0
      elect=1
      disso=2
      PID(1)=3
      PID(2)=0
    !TEX+SS
    elseif(i.eq.26)then
      tempQ=tempQ+1
      excite=0
      elect=1
      disso=0
      PID(1)=7
      PID(2)=1
    !TEX+SPEX
    elseif(i.eq.27)then
      tempQ=tempQ
      excite=1
      elect=0
      disso=0
      PID(1)=7
      PID(2)=3
    !SI+SS
    elseif(i.eq.28)then
      tempQ=tempQ+1
      excite=0
      elect=2
      disso=1
      PID(1)=1
      PID(2)=1
    !SI+SPEX
    elseif(i.eq.29)then
      tempQ=tempQ
      excite=1
      elect=1
      disso=1
      PID(1)=1
      PID(2)=3
    !SC
    elseif(i.eq.30)then
      tempQ=tempQ-1
      excite=1
      elect=0
      disso=1
      PID(1)=4
      PID(2)=0
    !DI+SS
    elseif(i.eq.31)then
      tempQ=tempQ+1
      excite=0
      elect=3
      disso=2
      PID(1)=2
      PID(2)=1
    !DI+SPEX
    elseif(i.eq.32)then
      tempQ=tempQ
      excite=1
      elect=2
      disso=2
      PID(1)=2
      PID(2)=3
    !TEX
    elseif(i.eq.33)then
      tempQ=tempQ
      excite=0
      elect=0
      disso=0
      PID(1)=7
      PID(2)=0
    !SI
    elseif(i.eq.34)then
      tempQ=tempQ
      excite=0
      elect=1
      disso=1
      PID(1)=1
      PID(2)=0
    !DI
    elseif(i.eq.35)then
      tempQ=tempQ
      excite=0
      elect=2
      disso=2
      PID(1)=2
      PID(2)=0
    !NEG
    elseif(i.eq.36)then
      tempQ=tempQ-1
      excite=0
      elect=0
      disso=1
      PID(1)=0
      PID(2)=0
    end if
    goto 1000
  end if
end do
if(ranVecB(1).ge.sumProb)then
  write(206,*)'CollisionSim.f08: ERROR: Random number greater than normalized &
  probability:', ranVecB(1), sumProb
  STOP 'CollisionSim.f08: Stopping program...'
end if
1000 continue
!call system_clock (t2,clock_rateTotal,clock_maxTotal)
!write(*,*) 'CollisionSim.f08: Collision process calculated &
!            in', real(t2-t1)/clock_rateTotal, 'seconds.'

end subroutine
