******************************************************************************
*                             Subroutine GRIDPAR
*                            ~~~~~~~~~~~~~~~~~~~~
*
*  Purpose:
*
*     Compute parameters of the grid for a given number of rings for
*     the SREAG pixelization
*
*  Usage:
*
*     CALL GRIDPAR ( NRING, NCELL, CA, B, CL, NCR, NC1R, DL, IER )
*
*  Input arguments:
*
*     NRING - number of rings, INTEGER(4);
*             must be even in the range [4:41068]
*
*  Output arguments:
*
*     NCELL - number of cells, INTEGER(4)
*        CA - cell area, sq. deg, REAL(8)
*         B - array of dimension NRING+1 of latitudinal ring boundaries, deg,
*             REAL(8)
*        CL - array of dimension NRING of central latitudes by rings, deg,
*             REAL(8)
*       NCR - array of dimension NRING of number of cells by rings, INTEGER(4)
*      NC1R - array of dimension NRING of number of the first cell in each ring,
*             INTEGER(4)
*        DL - array of dimension NRING of longitudinal cell span, deg, REAL(8)
*       IER - status code, INTEGER(4):
*             IER=0  - normal exit
*             IER=-1 - ilvalid NRING
*
*  Modified arguments:
*
*     none
*
*  Called routines:
*
*     none
*
*  Method:
*
*     Malkin Z. A new equal-area isolatitudinal grid on a spherical surface.
*     AJ, Vol. 158, No. 4, id. 158, 2019. DOI: 10.3847/1538-3881/ab3a44
*
*  Notes:
*
*     1. The rings are numbered beginning from the North pole.
*     2. Grid resolution can be computed as the square root of CA.
*
*
*                                                 (C) Zinovy Malkin, 2019-2020
******************************************************************************

      SUBROUTINE GRIDPAR ( NRING, NCELL, CA, B, CL, NCR, NC1R, DL, IER )

      IMPLICIT NONE

      SAVE

* Input variables

      INTEGER(4)  NRING          ! number of rings

* Output variables

      INTEGER(4)  NCELL          ! total # of cells in grid
      REAL(8)  CA                ! cell area, sq. deg
      REAL(8)  B(NRING+1)        ! latitudinal ring boundaries, deg
      REAL(8)  CL(NRING)         ! central latitude, deg
      INTEGER(4)  NCR(NRING)     ! # of cells by rings
      INTEGER(4)  NC1R(NRING)    ! number of first cells in rings
      REAL(8)  DL(NRING)         ! longitudinal cell span, deg
      INTEGER(4)  IER

* Modified variables



* Local variables

      INTEGER(4), PARAMETER :: NRINGMAX = 41068    ! INTEGER(4) limit
      INTEGER(4)  I, N2

      REAL(8), PARAMETER :: DEGRAD = 1.7453292519943295769D-2
      REAL(8), PARAMETER :: HALFPI = 1.5707963267948966192D0    ! pi/2
      REAL(8), PARAMETER :: FOURPI = 1.2566370614359172954D1    ! 4pi
* Square degrees on sphere.
      REAL(8), PARAMETER :: SDEGS = 41252.961249419277010D0
* Square degrees on hemisphere.
      REAL(8), PARAMETER :: SDEGHS = SDEGS / 2D0

      REAL(8)  DB          ! initial step in latitude, deg
      REAL(8)  CS          ! initial (non-integer) # of cells in ring
      REAL(8)  B2, CL0, CCL0

* ----------------------------------------------------------------------------

      IER = 0

* Check if the given NRING is valid.
      IF ( NRING < 4 .OR. NRING > NRINGMAX .OR. NRING/2*2 /= NRING )
     &     THEN
         IER = -1
         RETURN
      END IF

      N2 = NRING/2
      DB = 90D0 / DBLE(N2)                ! nominal ring width
      NCELL = 0

* --- Compute data for the North hemisphere.

      DO I=1,N2         ! rings from North pole
         CL0 = DB/2D0 + DB*DBLE(N2-I)     ! nominal central latitude
         CCL0 = DCOS(CL0*DEGRAD)
         CS = DB/CCL0      ! initial (non-integer) # of cells in ring
         NCR(I) = NINT(360D0/CS)    ! final (integer) # of cells in ring
         NCELL = NCELL + NCR(I)           ! total # of cells in grid
         DL(I) = 360D0/DBLE(NCR(I))       ! longitudinal cell span, deg
      END DO

* Cell area, sq. deg..
      CA = SDEGHS / DBLE(NCELL)

* --- Compute final latitudinal ring boundaries.

* Upper boundary, deg.
      B(1) = 90D0

      DO I=1,N2
* Lower boundary, deg.
         B2 = DASIN ( DSIN(B(I)*DEGRAD) -
     &      (CA/SDEGS*FOURPI)/(DL(I)*DEGRAD) ) / DEGRAD
* Central latutude, deg.
         CL(I) = (B(I)+B2) / 2D0
         B(I+1) = B2
      END DO

* Complete data for South hemisphere.
      DO I=N2+1,NRING
         NCR(I) = NCR(NRING-I+1)
         CL(I) = -CL(NRING-I+1)
         B(I+1) = -B(NRING-I+1)
         DL(I) = DL(NRING-I+1)
      END DO

* Compute the number of the first cell in each ring.

      NC1R(1) = 1

      DO I=2,NRING
         NC1R(I) = NC1R(I-1) + NCR(I-1)
      END DO

      NCELL = NCELL * 2

      END   ! GRIDPAR
