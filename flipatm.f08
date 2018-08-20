program flipatm

implicit none

integer length,i

parameter(length=1544)

real,dimension(length) :: alt,H2,He,CH4,H,Tot,Temp,Press,ScaleHeight

open(unit=100,file='./Atmosphere/Input/JunoAtmosphere_2km.dat',status='old')
open(unit=101,file='./Atmosphere/Input/JunoAtmosphere_2kmFlip.dat',status='unknown')

read(100,*)
do i=1,length
  read(100,200) alt(i),H2(i),He(i),CH4(i),H(i),Tot(i),Temp(i),Press(i),ScaleHeight(i)
end do
write(101,*) '# Alt [km]  H2 [cm^-3]   He [cm^-3]  CH4 [cm^-3]   H [cm^-3]   &
              Tot [cm^-3]  Temp [K]  Press [mBar] Scale height [km]'
do i=length,1,-1
  write(101,200) alt(i),H2(i),He(i),CH4(i),H(i),Tot(i),Temp(i),Press(i),ScaleHeight(i)
end do

200 format(1x,F8.2,2X,5(ES11.3E2,2X),F8.2,2x,ES11.3E2,5x,ES11.3E2)

end program
