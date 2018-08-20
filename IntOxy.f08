program IntegrateCX

use formatting
implicit real*8(a-h,o-z)

integer run,energy

real*8 altitude(1544),oxygenCX(1544,10),Ioxy(10)
real*8 Eion(12)
character(len=100) filename

data Eion/1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0,5000.0,10000.0,&
          25000.0/

do run=2,10
  energy=int(Eion(run))
  write(filename,'("./Output/",I0,"keV/Oxy_CX_Comb.dat")') energy
  open(unit=100,file=filename,status='old')
  read(100,*)

  do i=1,1544
    read(100,F05) altitude(i),(oxygenCX(i,j),j=1,10)
  end do
  Ioxy=sum(oxygenCX,dim=1)*2e5
  write(*,*) '--------------------------------------------------'
  write(*,*) 'Energy: ',energy
  write(*,*) 'Altitude integrated oxygen production from charge exchage: '
  write(*,"('   ONeg',7x,'O',10x,'O+',9x,'O2+',8x,'O3+',8x,'O4+',8x,'O5+',8x,&
        'O6+',8x,'O7+',8x,'O8+')")
  write(*,1000) (Ioxy(i),i=1,10)
end do


1000 format(10(2x,ES9.2))


end program
