******************************************************************************
*                             Subroutine POS2CN2
*                            ~~~~~~~~~~~~~~~~~~~~
*
*  Purpose:
*
*     Compute the cell number for a given point position for the SREAG
*     pixelization
*
*  Usage:
*
*     CALL POS2CN2 ( NRING, LON, LAT, ICELL, IRING, IER )
*
*  Input arguments:
*
*     NRING - number of rings, INTEGER(4);
*             must be even in the range [4:41068]
*       LON - longitude of the point, deg, REAL(8);
*             must be in the range [0:360)
*       LAT - latitude of the point, deg, REAL(8)
*
*  Output arguments:
*
*     ICELL - cell number starting from 1, INTEGER(4)
*     IRING - ring number starting from 1 and N, INTEGER(4);
*             point on the boundary will be assigned to upper ring
*      IER - status code, INTEGER*4:
*            IER=0  - normal exit
*            IER=-1 - invalid NRING
*            IER=-2 - invalid LON
*            IER=-3 - invalid LAT
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
*     1. If IER<0, ICELL=0 and IRING=0 returned.
*
*                                                 (C) Zinovy Malkin, 2019-2020
******************************************************************************

      SUBROUTINE POS2CN2 ( NRING, LON, LAT, ICELL, IRING, IER )

      IMPLICIT NONE

      SAVE

* Input variables

      INTEGER(4)  NRING       ! # of rings in selected subdivision

      REAL(8)  LON, LAT       ! longitude and latitude of point, deg

* Output variables

      INTEGER(4)  ICELL       ! cell number
      INTEGER(4)  IRING       ! ring number
      INTEGER(4)  IER

* Modified variables



* Local variables

      INTEGER(4), PARAMETER :: NRINGMAX = 41068    ! INTEGER(4) limit
      INTEGER(4)  I, IC, IER1
      INTEGER(4)  NRINGL / 0 /      ! last NRING value

* Subdivison data

      INTEGER(4)  NCELL             ! total # of cells in grid
      INTEGER(4)  NCR(NRINGMAX)     ! number of cells by rings
      INTEGER(4)  NC1R(NRINGMAX)    ! number of first cells in rings

      REAL(8)  CA                   ! cell area, sq. deg
      REAL(8)  CB(NRINGMAX+1)       ! latitudinal boundaries, deg
      REAL(8)  CL(NRINGMAX)         ! central latitude, deg
      REAL(8)  DL(NRINGMAX)         ! longitudinal cell span, deg

* ----------------------------------------------------------------------------

      IER = 0
      ICELL = 0
      IRING = 0

* Check if given NRING is valid.
      IF ( NRING < 4 .OR. NRING > NRINGMAX .OR. NRING/2*2 /= NRING )
     &     THEN
         IER = -1
         RETURN
      END IF

* Check if given longitude is in the range.
      IF ( LON < 0D0 .OR. LON >= 360D0 )  THEN
         IER = -2
         RETURN
      END IF

* Check if given latitude is in the range.
      IF ( LAT > 90D0 .OR. LAT < -90D0 )  THEN
         IER = -3
         RETURN
      END IF

* Get grid parameters during the first call and if NRING has changed since
* the previous call.
      IF ( NRING /= NRINGL )  THEN
         CALL GRIDPAR ( NRING, NCELL, CA, CB, CL, NCR, NC1R, DL, IER1 )
         IF ( IER1 /= 0 )  THEN
            IER = -1
            RETURN
         END IF
      END IF

* Ring number IRING.
      DO I=1,NRING
         IF ( LAT <= CB(I) .AND. LAT >= CB(I+1) )  THEN
            IRING = I
            EXIT
         END IF
      END DO

* Cell number in the IRING-th ring.
      IC = INT( LON / DL(IRING) ) + 1

* Cell number in the whole grid.
      ICELL = NC1R(IRING) + IC - 1

* Save the last given ring number.
      NRINGL = NRING

      END   ! POS2CN2
