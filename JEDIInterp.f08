program JEDIInterpolator
!*******************************************************************************
!* Created by Stephen J. Houston 9.5.18
!*******************************************************************************
!* This program interpolates JEDI ion spectra into a finer energy grid. It
!* outputs the intensity, energy bin width, and then a renormalized ion flux
!* which is too be input into the oxygen ion precipitation model. Alternatively,
!* the flux can instead be multiplied by a normalized (1 ion/cm^2/s) output from
!* JupOxyPrecip.
!*******************************************************************************

implicit real*8(a-h,o-z)

!*******************************************************************************
integer run

parameter(number_of_energies=10) !Number of JEDI energy bins
!* JEDI variables
real*8,dimension(number_of_energies) :: Jenergy,Jintensity,Jebins,Jflux
!*   Jenergy - JEDI energy bin [keV]
!*   Jintensity - JEDI ion flux [c/s/ster/cm^2/keV]
!*   Jebins - Size of JEDI energy bins [keV]
!*   Jflux - Jintensity values converted to [counts/cm^2/s]

character(len=10) date,version
character(len=12) time
character(len=100) filename

!****************************** Data Declaration *******************************
!* Width of JEDI energy bins: (May eventually need to be adjusted)
data Jebins/66.0,71.0,105.0,216.0,346.0,251.0,300.0,880.0,2280.0,5340.0/
!********************************* Initialize **********************************
pi=4.0d0*atan(1.0d0);Jenergy=0.0;Jintensity=0.0;Jbins=0.0;Jflux=0.0
!*************************** Open JEDI Ion Spectrum ****************************
write(version,'("v1")') !Filename of a JEDI spectrum (.d2s file)
write(filename,'("./JunoData/Spectra/",A,".d2s")') trim(version)
open(unit=100,file=trim(filename),status='old')
write(*,*)
do i=1,25 !Reading in the data measured by JEDI
  if(i.le.2.or.i.ge.4.and.i.le.15)read(100,*)
  if(i.eq.3)read(100,1001) date, time !Read the date and time of the flyby
  if(i.eq.3)write(*,1000) trim(version),date,time !Write to screen
  if(i.ge.16)read(100,1002) Jenergy(i-15),Jintensity(i-15)
end do
write(*,*)
write(*,1003)'Energy Bin:','JEDI Intensity:','Energy Bin Width:',&
             'Normalized Flux:' !Write out general information
do run=1,number_of_energies !Convert to [counts/cm^2/s]
!* The first 3 energy bins include both sulfur and oxygen. I'm assuming a 2:1
!* ratio of oxygen:sulfur (from SO_2)
  if(run.le.3)Jflux(run)=Jintensity(run)*2*pi*Jebins(run)*2/3
  if(run.ge.4)Jflux(run)=Jintensity(run)*2*pi*Jebins(run)
  write(*,1004)Jenergy(run),Jintensity(run),Jebins(run),Jflux(run)
end do
close(100) !Close JEDI measurement file

1000 format('FILE: ',A5,'.d2s',/,'DATE: ',A10,/,'TIME: ',A12)
1001 format(39X,A10,1X,A12)
1002 format(5X,ES12.9,1X,ES13.10)
1003 format(3x,A11,2x,A15,2x,A17,2x,A16)
1004 format(F11.3,1x,F15.4,2x,F17.2,2x,F16.3)

end program
