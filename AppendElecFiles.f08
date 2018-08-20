program append
!*******************************************************************************
!* recreated by Stephen J. Houston 4.20.18
!*******************************************************************************
!*******************************************************************************
!*
!* This program combines the electron ejection energy and angle cross-sections
!* from Schultz et. al, 2016(1-2000 keV/u) and Schultz et. al, 2018(5-25 MeV/u).
!* The energy bins vary between the two papers. To prevent gaps in the data,
!* there is a linear interpolation between data points if there would normally
!* be a gap.
!*
!* Collision Type:
!*    1) Single Ionization
!*    2) Double Ionization
!*    3) Transfer Ionization
!*    4) Double-Capture Auto Ionization
!*    5) Single Stripping
!*    6) Double stripping
!*
!* Possible Ion Energies:
!*    1, 10, 50, 75, 100, 200, 500, 1000, 2000, 5000, 10000, 25000
!*
!*******************************************************************************

implicit real*8(a-h,o-z)

!*******************************************************************************

integer err
integer,dimension(6,2) :: q_bounds !Range of charge states for each collision

real,allocatable,dimension(:) :: Ee1,Ee2 !Electron energy arrays
real,allocatable,dimension(:) :: E1,E10,E50,E75,E100,E200,E500,E1000,E2000,&
                                 E5000,E10000,E25000 !Ion Energies for XS's

character(len=100) energyFile !Full path to energy file that gets trimmed
character(len=100) angleFile  !Full path to angle file that gets trimmed
character(len=2) collision(6) !Abbreviation of collision in filename
character(len=4) collisionNew(6) !Abbreviation of collision in new filename

data collision /'si','di','ti','da','ss','ds'/ !File abbreviations
data collisionNew /'SI','DI','TI','DCAI','SS','DS'/
data q_bounds/1,1,2,3,1,1,& !Lower charge state bounds for collisions 1-6
              9,9,9,9,8,7/  !Upper charge state bounds for collisions 1-6

do i=1,6 !Do-loop for every collision type
  do j=q_bounds(i,1),q_bounds(i,2) !Do-loop through each charge state
!* Open energy file for a given collision type (i) and charge-state (j)
    write(energyFile, '("./OldElectronDist/sp-",A,"-E-O",I0)')& !1-2000 keV/u
    trim(collision(i)),j-1 !Specifies collision type (i) & charge-state (j)
    open(100,file=trim(energyfile),status='old') !Open file

    if(i.eq.4)goto 3000 !There are no new DCAI (da) ejection energies/angles
    write(energyfile, '("./UpdatedElectronDist/sp-",A,"E-O",I0,".dat")')& !5-25
    trim(collision(i)),j-1 !Specifies collision type (i) & charge-state (j)
    open(200,file=trim(energyfile),status='old') !Open file
    3000 continue
    write(energyfile, '("./NewElectronDist/",A,"_",I0,"_eng.dat")')& !energy file
    trim(collisionNew(i)),j-1 !Specifies collision type (i) & charge-state (j)
    open(300,file=trim(energyfile),status='unknown') !Open file

    maxLines=1000 !Initialize, this can be made larger without any problems
    nLines1=1 !Make the arrays 1 bigger than needed, but last index won't be
    nLines2=1 !filled with a value
    nLinesTot=0
    do k=1,maxLines !Determine how many energy/angle bins there are
      read(100,*,end=1000)
      if(k.eq.maxLines) write(*,*)'Error: maximum number of lines exceeded. &
      Increase maxlines for unit=100.'
      nLines1=nLines1+1 !Number of lines in the file
    end do !End do for k=1,maxLines
    1000 continue !Exit the do loop once the end of the file is reached
    do k=1,maxLines
      read(200,*,end=2000)
      if(k.eq.maxLines) write(*,*)'Error: maximum number of lines exceeded. &
      Increase maxlines for unit=200.'
      nLines2=nLines2+1
    end do
    2000 continue
    rewind(100)
    rewind(200)
    nLinesTot=nLines1+nLines2 !Make all arrays the size of both combined.
    write(*,*) 'Collision: ', i, '  Charge State: ', j-1 !Output check
    allocate(Ee1(nLines1)) !Electron energy bins
    allocate(E1(nLines1)) !Allocate all energies/angles to the size of the file
    allocate(E10(nLines1))
    allocate(E50(nLines1))
    allocate(E75(nLines1))
    allocate(E100(nLines1))
    allocate(E200(nLines1))
    allocate(E500(nLines1))
    allocate(E1000(nLines1))
    allocate(E2000(nLines1))
    if(i.ne.4)then
      allocate(Ee2(nLines2)) !5-25 MeV/u values
      allocate(E5000(nLines2))
      allocate(E10000(nLines2))
      allocate(E25000(nLines2))
    end if
    do k=1,nLines1-1 !Read in all of the cross-section values for 1-2000 keV/u
      read(100,*) Ee1(k),E1(k),E10(k),E50(k),E75(k),E100(k),E200(k),E500(k),&
                  E1000(k),E2000(k)
    end do !End do for k=1,nLines1
    if(i.eq.4)then !Set the 5-25 MeV/u to same length as 1-2000 for DCAI
      allocate(Ee2(nLines1))
      allocate(E5000(nLines1))
      allocate(E10000(nLines1))
      allocate(E25000(nLines1))
      Ee2=Ee1 !Set energy bins to the same
      goto 5000 !For DCAI, can skip everything and just write out the XS's
    end if
    do k=1,nLines2-1 !Read in all of the cross-section values for 5-25 MeV/u
      read(200,*) Ee2(k),E5000(k),E10000(k),E25000(k)
    end do
    close(200)
    5000 continue !Skip everything for DCAI
    close(100)
    n1=1
    n2=1
    err=0
    if(i.eq.4)err=4
    do k=1,nLinesTot
      if(Ee1(n1).eq.Ee2(n2).and.err.eq.0)then
        write(300,501) Ee1(n1),E1(n1),E10(n1),E50(n1),E75(n1),E100(n1),E200(n1),&
        E500(n1),E1000(n1),E2000(n1),E5000(n2),E10000(n2),E25000(n2)
        n1=n1+1 !Advance the Ee1 array
        n2=n2+1 !Advance the Ee2 array
        if(n1.eq.nLines1.and.n2.eq.nLines2) goto 4000 !Completed both arrays
        if(n1.eq.nLines1)err=1 !Completed the Ee1 array
        if(n2.eq.nLines2)err=2 !Completed the Ee2 array
      elseif(Ee1(n1).lt.Ee2(n2).and.err.eq.0)then
        !Linearly interpolate missing data points
        if(n2.eq.1)then
          write(300,503) Ee1(n1),E1(n1),E10(n1),E50(n1),E75(n1),E100(n1),E200(n1),&
          E500(n1),E1000(n1),E2000(n1)," 1.00E-30"," 1.00E-30"," 1.00E-30"
        else
          f=(Ee1(n1)-Ee2(n2-1))/(Ee2(n2)-Ee2(n2-1))
          write(300,501) Ee1(n1),E1(n1),E10(n1),E50(n1),E75(n1),E100(n1),E200(n1),&
          E500(n1),E1000(n1),E2000(n1),(1-f)*E5000(n2-1)+f*E5000(n2),&
          (1-f)*E10000(n2-1)+f*E10000(n2),(1-f)*E25000(n2-1)+f*E25000(n2)
        end if
        n1=n1+1 !Advance the Ee1 array
        if(n1.eq.nLines1)err=1 !Completed the Ee1 array
      elseif(Ee1(n1).gt.Ee2(n2).and.err.eq.0)then
        if(n1.eq.1)then
          write(300,502) Ee2(n2)," 1.00E-30"," 1.00E-30"," 1.00E-30"," 1.00E-30",&
          " 1.00E-30"," 1.00E-30"," 1.00E-30"," 1.00E-30"," 1.00E-30",&
          E5000(n2),E10000(n2),E25000(n2)
        else
          f=(Ee2(n2)-Ee1(n1-1))/(Ee1(n1)-Ee1(n1-1))
          write(300,501) Ee2(n2),(1-f)*E1(n1-1)+f*E1(n1),(1-f)*E10(n1-1)+f*E10(n1),&
          (1-f)*E50(n1-1)+f*E50(n1),(1-f)*E75(n1-1)+f*E75(n1),&
          (1-f)*E100(n1-1)+f*E100(n1),(1-f)*E200(n1-1)+f*E200(n1),&
          (1-f)*E500(n1-1)+f*E500(n1),(1-f)*E1000(n1-1)+f*E1000(n1),&
          (1-f)*E2000(n1-1)+f*E2000(n1),E5000(n2),E10000(n2),E25000(n2)
        end if
        n2=n2+1
        if(n2.eq.nLines2)err=2 !Completed the Ee2 array
      elseif(err.eq.1)then
        write(300,502) Ee2(n2)," 1.00E-30"," 1.00E-30"," 1.00E-30"," 1.00E-30",&
        " 1.00E-30"," 1.00E-30"," 1.00E-30"," 1.00E-30"," 1.00E-30",&
        E5000(n2),E10000(n2),E25000(n2)
        n2=n2+1
        if(n2.eq.nLines2)goto 4000
      elseif(err.eq.2.or.err.eq.4)then !4 is DCAI
        write(300,503) Ee1(n1),E1(n1),E10(n1),E50(n1),E75(n1),E100(n1),E200(n1),&
        E500(n1),E1000(n1),E2000(n1)," 1.00E-30"," 1.00E-30"," 1.00E-30"
        n1=n1+1
        if(n1.eq.nLines1)goto 4000
      end if
    end do !End do for k=1,nLinesTot
    !write(300,501)
    4000 continue
    deallocate(Ee1,E1,E10,E50,E75,E100,E200,E500,E1000,E2000,Ee2,E5000,E10000,&
    E25000)
  end do !End do for i=1,6
end do !End do for j=q_bounds(i,1),q_bounds(i,2)
close(300)
501 format (F10.3,12(ES9.2E2))
502 format (F10.3,9(A9),3(ES9.2E2))
503 format (F10.3,9(ES9.2E2),3(A9))


end program
