subroutine EnergyLoss(E,Q,electEnergy,PID,dE,eAngleSS,eEnergySS,eAngleDS,&
  eEnergyDS)!,eE,eA)
!*******************************************************************************
!* Created by Stephen J. Houston 2.23.18
!*******************************************************************************
!* This routine calculates the energy loss of a precpitating ion.
!* The energy loss obtained by energy loss models based on models
!* calculated by Schultz et. al. 2018 in: Data for secondary electron production
!* from ion precipitation at Jupiter II: Simultaneous and non-simultaneous
!* target and projectile processes in collisions of O^(q+)+H_2 (q=0-8)
!* Table A, pg. 10, Tables 1 and 2, pg. 31 and 32.
!* Also from Schultz et. al. 2016 Ionization of molecular hydrogen and stripping
!* of oxygen atoms and ions... Table 4, pg. 37.
!*******************************************************************************
!*
!* Input:
!*    E --> Total Energy/nucleon
!*     Type: Real
!*     Units: keV/u
!*
!*    tempQ --> Number of Ion Charges
!*     Type: Integer
!*	   Units: None
!*
!*    electEnergy --> Energy of Secondary Electron(s)
!*     Type: Real*8
!*     Units: eV
!*
!*    PID --> Collision Type (1-7,0-4) See CollisionSim.f08
!*     Type: 2-Component Array, Integer
!*     Units: None
!*
!* Returns:
!*    dE --> Delta Energy
!*     Type: Real
!*     Units: eV
!*
!*******************************************************************************
!*
!* FYI: For single stripping and double stripping, the electron energies need to
!* be boosted into the projectile frame rather than the target frame, produced
!* by the SDXS. This calculation can be found in the appendix of Schultz et. al.
!* 2018.
!*
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************

integer,intent(in) :: Q,PID(2) !Input charge state and process ID
real*8,intent(in)  :: E,electEnergy !Input ion energy and sec. electron energy
real*8,intent(in)  :: eAngleSS,eEnergySS !Single stripping conditions
real*8,intent(in)  :: eAngleDS(2),eEnergyDS(2) !Double stripping conditions
real*8,intent(out) :: dE!,eE,eA !Output the change in energy of the ion

integer ChS,tempQ !Charge states
real*8 IP1,IP2 !First two ionization potentials for H2
real*8 pi,c,oMass !Speed of light and oxygen mass (14.8952e6 GeV/c^2)
real*8 auCon !eV to a.u. conversion (1 a.u. = 27.2116 eV)
real*8 alpha !Fine-structure constant (1/137)
real*8 eMass !Electron mass (0.5109989 MeV/c^2)
real*8 ThetaP !Theta in the projectile frame

parameter (n=12,ChS=9) !Number of data points and charge states
parameter (c=137.036) !Speed of light in a.u.
parameter (oMass=14.903e6,auCon=27.2116) !oMass conversion to keV
parameter (pi=4.0*atan(1.0d0),eMass=0.5109989e3)
parameter (IP1=15.4254,IP2=16.4287)
!IP Data from J. Liu et al, Determination of the ionization and dissociation
!energy of the hydrogen molecule, J. Chem. Phys. 130, 174306 (2009).

integer bin,Ebin(n) !Bin variable and energy bins

real*8 Vproj !Projectile velocity (precipitating ion)
real*8 Vz !Electron velocity in projectile frame
real*8 Vsquared !The square of the ejected electron's vel. in projectile frame
real*8 electEnergySS,electEnergyDS1,electEnergyDS2,electEnergyAU1,electEnergyAU2
real*8 IPO(8),dENEG(n) !Oxy IP and dE for negative ion channel
real*8,dimension(n,ChS) :: dETI,dESC,dEDC,dESPEX,dEDPEX,dESS,dEDS
!dE for process that don't eject electrons (Table)

!****************************** Data Declaration *******************************

!* Energy bins determined by Schultz et al. 2018 in Table 1 (keV/u)
data Ebin/1,10,50,75,100,200,500,1000,2000,5000,10000,25000/

!* Ionization potentials of oxygen (eV)
data IPO/13.62,35.12,54.94,77.41,113.9,138.12,739.28,871.41/

!* Energy losses for negative ion channel (Schultz et al. 2018 Table ) (eV)
data dENEG/14.5,19.4,41.2,54.8,68.4,123,286,559,1103,2737,5460,13629/

!* New energy loss/gain is in Schultz et al. 2018 Tables 1-11 (eV)
!* Energy loss/gain for Transfer Ionization (Schultz et al. 2016 Table 4) (eV)
data dETI/&
  0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&   !O
  -7.47,-2.57,19.2,32.8,46.4,101,264,537,1081,2715,5438,13607,&   !O+
  -21.1,-16.2,26.2,39.8,32.8,87.3,251,523,1068,2702,5425,13594,&  !O++
  -7.92,-3.02,18.8,32.4,46.0,101,234,506,1051,2685,5408,13577,&   !O3+
  -16.7,-11.3,10.0,23.6,37.2,91.7,214,487,1031,2665,5388,13557,&  !O4+
  -26.9,-22.0,-0.21,13.4,27.0,81.5,245,464,1008,2642,5365,13534,& !O5+
  -16.2,-11.3,10.5,24.1,37.7,67.2,159,431,976,2610,5333,13502,&   !O6+
  -26.6,-21.7,15.4,13.7,27.3,81.8,212,390,934,2568,4781,12950,&   !O7+
  -38.5,-33.6,18.5,32.1,45.7,69.9,191,463,1008,2521,4591,12760/   !O8+

!* Energy loss/gain for Single Capture Charge Transfer (Schultz 2016 Table 4)
data dESC/&
  0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
  -7.47,-2.57,19.2,32.8,46.4,101,264,537,1081,2715,5438,13607,&
  -21.1,-16.2,26.2,39.8,53.4,87.3,251,523,1068,2701,5425,13594,&
  -7.92,-3.02,18.8,32.4,46.0,101,234,506,1051,2685,5408,13547,&
  -16.7,-11.3,10.0,23.6,37.2,91.7,214,487,1031,2665,5388,13557,&
  -8.14,-3.24,18.5,32.2,45.8,81.5,245,464,1008,2642,5365,13534,&
  -16.2,-11.3,10.5,35.7,49.3,104,231,431,976,2610,5333,13502,&
  -26.6,-21.7,15.4,29.0,42.6,97.1,212,484,1029,2568,4781,12950,&
  -18.9,-14.0,18.5,32.1,45.7,89.5,233,463,1008,2520,4591,12760/

!* Energy loss/gain for Double Capture Charge Transfer (Schultz 2016 Table 4)
data dEDC/&
  0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
  0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
  -20.4,-10.6,62.6,89.8,117,196,523,1068,2157,5425,10871,27209,&
  -3.60,6.21,49.8,77.0,104,213,494,1039,2128,5396,10842,27181,&
  -16.7,-6.88,36.7,63.9,91.2,200,465,1010,2099,5466,10813,27104,&
  -6.68,3.12,46.7,73.9,101,179,506,963,2052,5320,10766,27104,&
  -21.4,-11.6,32.0,78.8,106,215,480,904,1993,5261,10630,26969,&
  -40.6,-30.8,39.2,66.5,93.7,203,446,990,2080,5347,10630,26969,&
  -28.3,-18.5,43.8,71.0,98.2,188,481,951,2040,5095,9394,25732/

!* Energy loss for Single Projectile Excitation (Schultz 2018 Table 2)
data dESPEX/&
  12.3,12.5,12.5,12.5,12.5,12.5,12.4,12.4,12.4,12.4,12.4,12.4,&
  29.2,30.3,30.4,30.4,30.3,30.3,30.3,30.2,30.2,30.1,30.1,30.1,&
  41.0,42.5,44.1,43.9,43.8,43.6,43.4,43.3,43.2,43.0,43.0,42.9,&
  47.8,52.3,56.7,57.3,57.0,56.4,55.9,55.6,55.5,55.1,55.0,54.9,&
  63.1,70.3,78.2,80.5,81.5,80.8,79.9,79.2,78.9,78.3,78.1,77.9,&
  63.6,71.1,82.8,86.2,88.5,88.5,86.5,85.3,84.9,84.1,83.3,83.3,&
  420.,469.,479.,490.,497.,515.,543.,550.,546.,544.,542.,540.,&
  586.,508.,526.,533.,538.,566.,601.,620.,614.,608.,606.,603.,&
  0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/

!* Energy loss for Double Projectile Excitation (Schultz 2018 Table 2)
data dEDPEX/&
  48.0,54.9,61.7,64.2,64.9,63.2,65.9,63.7,58.5,65.1,66.1,62.3,&
  41.0,46.7,47.4,47.3,46.8,46.6,46.7,46.1,45.8,45.7,45.5,44.8,&
  97.0,108.,111.,111.,111.,110.,110.,110.,110.,110.,109.,110.,&
  120.,147.,162.,164.,168.,170.,169.,170.,169.,167.,169.,171.,&
  140.,168.,181.,184.,186.,188.,186.,185.,184.,183.,183.,0.00,&
  550.,620.,661.,662.,655.,685.,758.,789.,794.,779.,790.,0.00,&
  800.,900.,1000,1302,1115,1183,1211,1230,1242,1214,1237,0.00,&
  0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
  0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/

!* Energy loss for Single Stripping
data dESS/&
  21.1,24.0,44.9,50.5,53.8,58.8,60.7,60.5,59.0,57.8,58.6,58.0,&
  46.9,46.1,65.2,75.0,83.0,100.,116.,126.,131.,138.,143.,145.,&
  68.4,67.2,82.3,93.3,102.,128.,159.,177.,189.,209.,216.,245.,&
  92.2,91.6,105.,114.,123.,154.,199.,231.,259.,286.,321.,346.,&
  133.,131.,145.,153.,159.,193.,257.,303.,350.,399.,444.,531.,&
  149.,160.,172.,179.,185.,218.,293.,354.,408.,470.,503.,552.,&
  829.,839.,858.,864.,897.,915.,955.,1084,1339,1780,2135,2396,&
  971.,991.,1011,1019,1045,1080,1103,1227,1490,2018,2377,2698,&
  0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/

!* Energy loss for Double Stripping
data dEDS/&
  60.7,75.3,110.,127.,139.,175.,219.,241.,244.,234.,244.,243.,&
  107.,116.,151.,172.,189.,230.,303.,335.,351.,334.,349.,384.,&
  142.,161.,192.,213.,228.,284.,380.,421.,446.,452.,530.,551.,&
  201.,221.,263.,275.,291.,352.,462.,549.,617.,696.,735.,864.,&
  267.,292.,334.,342.,352.,413.,543.,637.,720.,823.,912.,1280,&
  897.,942.,1027,1077,1111,1214,1264,1506,1883,2440,2838,3344,&
  1691,1811,2011,2061,2109,2160,2134,2318,2894,3481,3972,5389,&
  0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&
  0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/

!******************************** Main Program *********************************

tempQ=Q
if(tempQ.gt.1)then !Treat negative ion electron ejection same as neutral
  tempQ=tempQ-1!Need to adjust the arrays since they don't include negative ions
end if !But if q equals 1, then treat it as a neutral oxy energy loss.

!* Determine which bin is to be used based on the oxygen energy
do i=12,1,-1 !This will only apply to the dE(process) variables
  if(E.ge.real(Ebin(i-1)+(Ebin(i)-Ebin(i-1))/2.0))then
    bin=i
    k=1
    goto 10
  elseif(E.ge.real(Ebin(i-1)))then
    bin=i-1
    k=2
    goto 10
  end if
end do
10 continue
if (bin.eq.0) write(*,*) 'EnergyLoss.f08: Error! i=0'

!* Want to use f to linearly interpolate energy losses between big bins for the
!* dE(process) data points.
if (k.eq.1) f=(E-Ebin(bin-1))/(Ebin(bin)-Ebin(bin-1))
if (k.eq.2) f=(E-Ebin(bin))/(Ebin(bin+1)-Ebin(bin))


dE=0.0 !Initialize

!* Find the dE for the given process based on Schultz et al. 2018 Table 1
!* Go through the NSIM processes
if(PID(1).eq.0)then !Negative Ion Charge State
  if(f.ge.0.5)dE=(f*dENEG(bin)+(1-f)*dENEG(bin-1))
  if(f.lt.0.5)dE=(f*dENEG(bin+1)+(1-f)*dENEG(bin))
elseif(PID(1).eq.1)then !Single Ionization
  dE=IP1+electEnergy
elseif(PID(1).eq.2)then !Double Ionization
  dE=IP1+IP2+electEnergy !Both electron energies have been added together
elseif(PID(1).eq.3)then !Transfer Ionization
!  dE=IP1+electenergy+dETI(bin,tempQ)
  if(f.ge.0.5)dE=IP1+electenergy+(f*dETI(bin,tempQ)+(1-f)*dETI(bin-1,tempQ))
  if(f.lt.0.5)dE=IP1+electenergy+(f*dETI(bin+1,tempQ)+(1-f)*dETI(bin,tempQ))
elseif(PID(1).eq.4)then !Single Capture
  if(f.ge.0.5)dE=(f*dESC(bin,tempQ)+(1-f)*dESC(bin-1,tempQ))
  if(f.lt.0.5)dE=(f*dESC(bin+1,tempQ)+(1-f)*dESC(bin,tempQ))
elseif(PID(1).eq.5)then !Double Capture
  if(f.ge.0.5)dE=(f*dEDC(bin,tempQ)+(1-f)*dEDC(bin-1,tempQ))
  if(f.lt.0.5)dE=(f*dEDC(bin+1,tempQ)+(1-f)*dEDC(bin,tempQ))
elseif(PID(1).eq.6)then !Double-Capture Autoionization
  dE=electEnergy
elseif(PID(1).eq.7)then !Target Excitation
  dE=7.7
end if
 dE=0.0 !~
!Go through the SIM processes
! if(PID(2).eq.1)then !Single Stripping Average Eloss
!   if(f.ge.0.5)dE=(f*dESS(bin,tempQ)+(1-f)*dESS(bin-1,tempQ))
!   if(f.lt.0.5)dE=(f*dESS(bin+1,tempQ)+(1-f)*dESS(bin,tempQ))
if(PID(2).eq.1)then !Single Stripping
  electEnergyAU1=eEnergySS/auCon !Convert to a.u.
  Vproj=sqrt(2.0*E*16.0/oMass)*c !In a.u.
  Vz=sqrt(2.0*electEnergyAU1)*cos(eAngleSS*pi/180.0)-Vproj !In a.u.
  Vsquared=(2*electEnergyAU1)-(2*Vz*Vproj)-Vproj**2 !Electron velocity squared
  electEnergySS=(0.5)*Vsquared*auCon !Convert back to eV
  ThetaP=acos(Vz/sqrt(Vsquared))*180/pi
  dE=dE+IPO(tempQ)+electEnergySS
  ! eE=electEnergySS
  ! eA=ThetaP
!  write(*,*) eEnergySS,electEnergySS,eAngleSS,ThetaP,IPO(tempQ)
!**  write(*,10000) eEnergySS,eAngleSS,electEnergyAU1,Vproj,Vz,Vsquared,electEnergySS,dE
! elseif(PID(2).eq.2)then !Double Stripping Average Eloss
!   if(f.ge.0.5)dE=(f*dEDS(bin,tempQ)+(1-f)*dEDS(bin-1,tempQ))
!   if(f.lt.0.5)dE=(f*dEDS(bin+1,tempQ)+(1-f)*dEDS(bin,tempQ))
elseif(PID(2).eq.2)then !Double Stripping
  electEnergyAU1=eEnergyDS(1)/auCon !Convert to a.u.
  Vproj=sqrt(2.0*E*16.0/oMass)*c !In a.u.
  Vz=sqrt(2.0*electEnergyAU1)*cos(eAngleDS(1)*pi/180.0)-Vproj !In a.u.
  Vsquared=(2*electEnergyAU1)-(2*Vz*Vproj)-Vproj**2 !Electron velocity squared
  electEnergyDS1=(0.5)*Vsquared*auCon !Convert back to eV
!  if(E.lt.360.5)write(*,*) Vproj,Vz,Vsquared,electenergyDS1,eEnergyDS(1),eAngleDS(1),pi
  electEnergyAU2=eEnergyDS(2)/auCon !Convert to a.u.
  Vproj=sqrt(2.0*E*16.0/oMass)*c !In a.u.
  Vz=sqrt(2.0*electEnergyAU2)*cos(eAngleDS(2)*pi/180.0)-Vproj !In a.u.
  Vsquared=(2*electEnergyAU2)-(2*Vz*Vproj)-Vproj**2 !Electron velocity squared
  electenergyDS2=(0.5)*Vsquared*auCon !Convert back to eV
!  if(E.lt.360.5)write(*,*) Vproj,Vz,Vsquared,electenergyDS2,eEnergyDS(2),eAngleDS(2),pi
  dE=dE+IPO(tempQ)+IPO(tempQ+1)+electenergyDS1+electenergyDS2
!  if(E.lt.360.5)write(*,*) IPO(tempQ)+IPO(tempQ+1)+electenergyDS1+electenergyDS2
elseif(PID(2).eq.3)then !Single Projectile Excitation
  if(f.ge.0.5)dE=dE+(f*dESPEX(bin,tempQ)+(1-f)*dESPEX(bin-1,tempQ))
  if(f.lt.0.5)dE=dE+(f*dESPEX(bin+1,tempQ)+(1-f)*dESPEX(bin,tempQ))
elseif(PID(2).eq.4)then !Double Projectile Excitation
  if(f.ge.0.5)dE=dE+(f*dEDPEX(bin,tempQ)+(1-f)*dEDPEX(bin-1,tempQ))
  if(f.lt.0.5)dE=dE+(f*dEDPEX(bin+1,tempQ)+(1-f)*dEDPEX(bin,tempQ))
end if

!**10000 format(5x,F9.2,20x,F9.2,11x,3(F9.3,5x),2x,F9.2,10x,F9.2,10x,F9.2)!,F9.2,2x,F8.2,2x,F8.2,)
!**if(E.lt.360.5)write(*,*) dE, IP1, IP2, electenergyDS1, electenergyDS2
end subroutine
