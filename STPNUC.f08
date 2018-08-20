FUNCTION STPNUC(EP)
!C
!C  Updated 6-30-92 by Edward J. Howell
!C
!C  Purpose: To calculate the nuclear stopping power of Oxygen in
!C           Hydrogen based on the Thomas-Fermi Int. Mod.
!C
!C  Argument List: EP
!C
!C  Input:
!C        EP = Proton Energy
!C        Type: REAL
!C        Units: KeV/amu
!C
!C  Returns:
!C        STPNUC = Nuclear Stopping Power
!C        Type: REAL
!C        Units: eV*cm**2
!C
!C  External Procedures: None
!C
!C  Limitations: None
!C
!C *********************** Declare Variables *********************
!C
IMPLICIT NONE
INTEGER I
REAL EP, GP, F, X, A, STPNUC
!C
!C ********************** Dimension Variables *******************
!C
DIMENSION GP(11), F(11)
!C
!C ********************** Data Statements ***********************
!C
DATA GP /0.000, 0.001, 0.002, 0.003, 0.004, 0.005, 0.010, 0.015, 0.020, 0.022, 0.500/

DATA F /0.0, 1.5, 2.0, 2.1, 2.0, 1.9, 1.2, 0.8, 0.6, 0.5, 0.0/

!C
!C ******************** Main Program Center *********************
!C
stpnuc = 0.0
X = 0.23*SQRT(EP*1.E-3)
IF (EP.LE.(4725.9)) THEN
  DO I = 11, 2, -1
    IF (X.LE.GP(I)) THEN
      A = F(I-1) + (F(I)-F(I-1))/(GP(I)-GP(I-1))*(X-GP(I-1))
      STPNUC = A*0.497E-15
    ENDIF
  ENDDO
ELSE
  STPNUC = 0.00000
ENDIF
RETURN
END
