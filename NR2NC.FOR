******************************************************************************
*                              Function NR2NC
*                             ~~~~~~~~~~~~~~~~~
*
*  Purpose:
*
*     Compute the number of cells for a given number of rings for the SREAG
*     pixelization
*
*  Usage:
*
*     NR2NC ( NRING )
*
*  Input arguments:
*
*     NRING - number of rings, INTEGER(4);
*             must be even in the range [4:41068]
*
*  Returned value:
*
*     Number of cells, INTEGER(4)
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
*     1. If invalid NRING is given, zero value returned.
*
*
*                                                 (C) Zinovy Malkin, 2019-2020
******************************************************************************

      FUNCTION NR2NC ( NRING )

      IMPLICIT NONE

* Function

      INTEGER(4)  NR2NC

* Input variables

      INTEGER(4)  NRING          ! number of rings

* Local variables

      INTEGER(4), PARAMETER :: NRINGMAX = 41068    ! INTEGER(4) limit
      INTEGER(4)  I, N2, NC, NCR, IER

      REAL(8), PARAMETER :: DEGRAD = 1.7453292519943295769D-2
      REAL(8)  DB          ! initial step in latitude, deg
      REAL(8)  CS          ! initial (non-integer) # of cells in ring
      REAL(8)  CL0, CCL0

* ----------------------------------------------------------------------------

      IER = 0

* Check if given NRING is valid.
      IF ( NRING < 4 .OR. NRING > NRINGMAX .OR. NRING/2*2 /= NRING )
     &     THEN
         NR2NC = 0
         RETURN
      END IF

      N2 = NRING/2
      DB = 90D0 / DBLE(N2)      ! nominal ring width, deg
      NC = 0

* --- Compute data for the North hemisphere.

      DO I=1,N2                 ! rings from North pole
         CL0 = DB/2D0 + DB*DBLE(N2-I)    ! nominal central latitude, deg
         CCL0 = DCOS(CL0*DEGRAD)
         CS = DB/CCL0           ! initial (real) # of cells in the ring
         NCR = NINT(360D0/CS)   ! final (integer) # of cells in the ring
         NC = NC + NCR          ! total # of cells in the grid
      END DO

      NR2NC = NC + NC

      END   ! NR2NC
