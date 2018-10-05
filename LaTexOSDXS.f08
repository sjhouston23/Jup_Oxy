program LaTexOSDXS
!*******************************************************************************
!* Created by Stephen J Houston 10.4.18
!*******************************************************************************
!* This program is designed to read in oxygen singly differential cross-section
!* data and reformat it into LaTex tables.
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
parameter(neProc=6)

real*8, dimension(6,9,13,149) :: edist
real*8, dimension(6,9,13,55) :: adist

character(len=30),dimension(neProc) :: LaTexProc
character(len=6),dimension(2) :: engang

data LaTexProc/'single ionization','double ionization','transfer ionization',&
'double-capture autoionization','single stripping','double stripping'/
data engang/'energy','angle'/
!*******************************************************************************
call ReadElectDist(edist,adist) !Read in all electron distribution information

open(unit=100,file='./LaTex.txt')
no=1
do i=1,6 !Collision type
  do m=1,2 !Energy then angle tables
    do j=1,9 !Charge state
      if(m.eq.1)then
        if(edist(i,j,1,1).eq.0.00)goto 20000
      elseif(m.eq.2)then
        if(adist(i,j,1,1).eq.0.00)goto 20000
      end if
      write(100,*)
      write(100,1002)
      write(100,1001)
      if(j.eq.1)write(100,1003) no,trim(engang(m)),trim(LaTexProc(i))
      if(j.eq.2)write(100,1004) no,trim(engang(m)),trim(LaTexProc(i))
      if(j.ge.3)write(100,1005) no,trim(engang(m)),trim(LaTexProc(i)),j-1
      no=no+1
      write(100,1006)
      write(100,1007)
      write(100,1008)
      if(m.eq.1)write(100,1009)
      if(m.eq.2)write(100,1010)
      write(100,1007)
      if(m.eq.1)then
        do l=1,149 !Ejected electron energy (k=1) AND ion energy (k=2-13)
          if(edist(i,j,1,l).eq.0.00)goto 10000
          write(100,1000) (edist(i,j,k,l),k=1,13)
        end do
      elseif(m.eq.2)then
        do l=1,55
          if(adist(i,j,1,l).eq.0.00)goto 10000
          write(100,1000) (adist(i,j,k,l),k=1,13)
        end do
      end if
      10000 continue
      write(100,1011)
      20000 continue
    end do
  end do
end do

1000 format(F9.3,12(1x,ES8.2E2))
1001 format("\begin{singlespace}")
1002 format("\scriptsize")
1003 format("Table B.",I0,": The singly differential cross-section as a &
function of the ejected electron ",A,/," for ",A," in 1 - 25000 keV/u &
O + H$_2$ collisions. (From \cite{schultz2018})")
1004 format("Table B.",I0,": The singly differential cross-section as a &
function of the ejected electron ",A,/," for ",A," in 1 - 25000 keV/u &
O$^{\rm +}$ + H$_2$ collisions. (From \cite{schultz2018})")
1005 format("Table B.",I0,": The singly differential cross-section as a &
function of the ejected electron ",A,/," for ",A," in 1 - 25000 keV/u &
O$^{\rm ",I0,"+}$ + H$_2$ collisions. (From \cite{schultz2018})")
1006 format("\begin{Verbatim}[baselinestretch=0.85]")
1007 format("------------------------------------------------------------------&
             ---------------------------------------------------")
1008 format("Electron    Ion Energy (keV/u)")
1009 format("Energy (eV)  1        10       50       75      100      200      &
             500      1000     2000     5000    10000    25000")
1010 format("Angle (deg) 1        10       50       75       100      200      &
            500      1000     2000     5000    10000    25000")
1011 format("\end{Verbatim}",/,"\end{singlespace}",/,/,"\newpage")

end program
