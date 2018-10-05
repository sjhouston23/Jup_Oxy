subroutine ReadElectDist(edist,adist)
!*******************************************************************************
!* Created by Stephen J. Houston 3.7.18 from previous work
!*******************************************************************************
!* This subroutine will calculate input all electron energy and angle
!* SDXS files and output them as two 4D arrays. The arrays are
!* indexed by collision-type and ion charge state so that all information
!* can be accessed from 2 variables.
!**********************************************************************
!**********************************************************************
!C
!C 	Input:
!C		Nothing
!C
!C  Returns:
!C		edist --> single differential cross section as a function of ejected electron energy
!C		Type: Real
!C		Units: cm^-2
!C    *See comments under definition of edist for more information.
!C
!C		adist --> single differential cross section as a function of ejected electron angle
!C		Type: Real
!C		Units: cm^-2
!C    *See comments under definition of edist for more information.
!
!**********************************************************************
!C
implicit none
!C
!********************** DECLARATION OF VARIABLES **********************
integer c,q,i,j !* Indexing variables (See edist and adist commenting below)
character(len=100) energyfile !* For writing electron energy filename for opening and reading
character(len=100) anglefile  !* For writing electron angle filename for opening and reading
!* Define a name array for oppening the files
character(len=4) collision(6)
data collision/'SI','DI','TI','DCAI','SS','DS'/
!* Charge state limits per each collison type
integer, dimension(6,2) :: q_limit
DATA q_limit/&
1,1,2,3,1,1,& ! q_limit(c,1) is minimum charge-state per collision type
9,9,9,9,8,7/  ! q_limit(c,2) is maximum charge-state per collision type
!* nrowse & nrowsa are the number of rows per collision type in either
!* the energy files (nrowse, hence "e") and angle files (nrowsa, hence "a")
integer, dimension(6) :: nrowse
DATA nrowse/133,115,72,64,149,138/
integer, dimension(6) :: nrowsa
DATA nrowsa/46,46,55,46,48,43/
!****************** Electron Energy SDXS File Sizes *******************
!* Size of arrays for inputting electron energy file data (eng+nrows)
real*8, dimension(13,133) :: eng133 !* Sized for SI
real*8, dimension(13,115) :: eng115 !* Sized for DI
real*8, dimension(13,72)  :: eng72  !* Sized for TI
real*8, dimension(13,64)  :: eng64  !* Sized for DCAI
real*8, dimension(13,149) :: eng149 !* Sized for SS
real*8, dimension(13,138) :: eng138 !* Sized for DS
!**********************************************************************
!****************** Electron Angle SDXS File Sizes ********************
!* Size of arrays for inputting electron angle file data
real*8, dimension(13,46) :: ang46 !* Sized for SI,DI,DCAI
real*8, dimension(13,55) :: ang55 !* Sized for TI
!real*8, dimension(13,46) :: ang46 !* Sized for DCAI
real*8, dimension(13,48) :: ang48 !* Sized for SS
real*8, dimension(13,43) :: ang43 !* Sized for DS
!**********************************************************************
!* Set matricies big enough for any of the cross-section distributions
real*8, dimension(6,9,13,149) :: edist
real*8, dimension(6,9,13,55) :: adist
! edist(c,q,i,j) & adist(c,q,i,j):
! - Indexing:
!                c - Collision Type:
!                                   1) Single Ionization
!                                   2) Double Ionization
!                                   3) Transfer Ionization
!                                   4) Double-Capture Auto Ionization
!                                   5) Single Stripping
!                                   6) Double stripping
!
!                 q - Charge state + 1 (O4+ -> q=5)
!
!                 i - distribution column
!                                   1) Ejected electron (energy list (eV)/ angle list (deg))
!                                   2-13) Ion Energy [1,10,50,75,100,200,500,1000,2000,5000,10000,25000](keV/u)
!
!                 j - distribution row (i.e. which ejected electron energy(angle) in eV(deg))
!
!******* Start Looping Through All 96 Differential Cross Section Files *********
do c=1,6 !Collision types, see indexing comments above
    do q=q_limit(c,1),q_limit(c,2) !Loop through charge states
!* Open energy file for a given collision type (c) and charge-state (q)
        write(energyfile, '("./NewElectronDist/",A,"_",I0,"_eng.dat")')& !* energy file name
        trim(collision(c)),q-1 !* Specifies collision type (c) & charge-state (q)
        open(1,file=energyfile,status='old') !* Open energy file

        write(anglefile, '("./NewElectronDist/",A,"_",I0,"_ang.dat")')& !* angle file name
        trim(collision(c)),q-1 !* Specifies collision type (c) & charge-state (q)
        open(2,file=anglefile,status='old') !* Open angle file

!* Adjust input 2D array size for corresponding size of electron energy file
      if(c.eq.1)then         !* SI......................................
          read(1,*) eng133          !* Input electron energy sdxs distribution
          read(2,*) ang46           !* Input electron angle sdxs distribution
      elseif(c.eq.2)then     !* DI...........................................
          read(1,*) eng115          !* Input electron energy sdxs distribution
          read(2,*) ang46           !* Input electron angle sdxs distribution
      elseif(c.eq.3)then     !* TI.........................................
          read(1,*) eng72           !* Input electron energy sdxs distribution
          read(2,*) ang55           !* Input electron angle sdxs distribution
      elseif(c.eq.4)then     !* DCAI......................................
          read(1,*) eng64           !* Input electron energy sdxs distribution
          read(2,*) ang46           !* Input electron angle sdxs distribution
      elseif(c.eq.5)then     !* SS......................................
          read(1,*) eng149          !* Input electron energy sdxs distribution
          read(2,*) ang48           !* Input electron angle sdxs distribution
      elseif(c.eq.6)then     !* DS......................................
          read(1,*) eng138          !* Input electron energy sdxs distribution
          read(2,*) ang43           !* Input electron angle sdxs distribution
      endif
      1000 continue
!* Add individually inputted energy cross-sections into larger arrays
      do i=1,13 !* Loop through columns of individual files/arrays
          do j=1,nrowse(c) !* Loop through rows of individual files/arrays
            if(c.eq.1)then            !* SI............................
              edist(c,q,i,j) = eng133(i,j)
              eng133(i,j) = 0.0 !Clear temporary array after filling edist()
            elseif(c.eq.2)then        !* DI.................................
              edist(c,q,i,j) = eng115(i,j)
              eng115(i,j) = 0.0 !Clear temporary array after filling edist()
            elseif(c.eq.3)then        !* TI...............................
              edist(c,q,i,j) = eng72(i,j)
              eng72(i,j) = 0.0 !Clear temporary array after filling edist()
            elseif(c.eq.4)then        !* DCAI............................
              edist(c,q,i,j) = eng64(i,j)
              eng64(i,j) = 0.0 !Clear temporary array after filling edist()
            elseif(c.eq.5)then        !* SS............................
              edist(c,q,i,j) = eng149(i,j)
              eng149(i,j) = 0.0 !Clear temporary array after filling edist()
            elseif(c.eq.6)then        !* DS............................
              edist(c,q,i,j) = eng138(i,j)
              eng138(i,j) = 0.0 !Clear temporary array after filling edist()
            endif
          enddo
      enddo
!* Add individually inputted angle cross-sectons into larger arrays
      do i=1,13 !* Loop through columns of individual files/arrays
          do j=1,nrowsa(c) !* Loop through rows of individual files/arrays
            if(c.eq.1)then            !* SI.................
              adist(c,q,i,j) = ang46(i,j)
              ang46(i,j) = 0.0 !Clear temporary array after filling edist()
            elseif(c.eq.2)then        !* DI............................
              adist(c,q,i,j) = ang46(i,j)
              ang46(i,j) = 0.0 !Clear temporary array after filling edist()
            elseif(c.eq.3)then        !* TI............................
              adist(c,q,i,j) = ang55(i,j)
              ang55(i,j) = 0.0 !Clear temporary array after filling edist()
            elseif(c.eq.4)then        !* DCAI............................
              adist(c,q,i,j) = ang46(i,j)
              ang46(i,j) = 0.0 !Clear temporary array after filling edist()
            elseif(c.eq.5)then        !* SS............................
              adist(c,q,i,j) = ang48(i,j)
              ang48(i,j) = 0.0 !Clear temporary array after filling edist()
            elseif(c.eq.6)then        !* DS............................
              adist(c,q,i,j) = ang43(i,j)
              ang43(i,j) = 0.0 !Clear temporary array after filling edist()
            endif
          enddo
      enddo
    enddo
enddo
close(1)
close(2)
return
end
