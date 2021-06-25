******************************************************************************
*                             Subroutine CN2POS2
*                            ~~~~~~~~~~~~~~~~~~~~
*
*  Purpose:
*
*     Compute the cell center coordinates for a given cell number for
*     the SREAG pixelization
*
*  Usage:
*
*     CALL CN2POS2 ( NRING, ICELL, LON, LAT, IRING, IER )
*
*  Input arguments:
*
*     NRING - number of rings, INTEGER(4);
*             must be even in the range [4:41068]
*     ICELL - cell number, INTEGER(4)
*
*  Output arguments:
*
*       LON - longitude of the center of the cell, deg, REAL(8)
*       LAT - latitude of the center of the cell, deg, REAL(8)
*     IRING - ring number, INTEGER*4
*       IER - status code, INTEGER*4:
*             IER=0  - normal exit
*             IER=-1 - invalid NRING
*             IER=-2 - invalid ICELL
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
*
*  Notes:
*
*     1. Zero values are returned in case of IER /= 0.
*
*                                                 (C) Zinovy Malkin, 2019-2020
******************************************************************************

      SUBROUTINE CN2POS2 ( NRING, ICELL, LON, LAT, IRING, IER )

      IMPLICIT NONE

      SAVE

* Input variables

      INTEGER(4)  NRING       ! # of rings in selected subdivision
      INTEGER(4)  ICELL       ! cell number

* Output variables

      INTEGER(4)  IER, IRING

      REAL(8)  LON, LAT       ! coordinates of the cell center, deg

* Modified variables



* Local variables

      INTEGER(4), PARAMETER :: NRINGMAX = 41068    ! INTEGER(4) limit
      INTEGER(4)  NRINGL / 0 /      ! last NRING value
      INTEGER(4)  I, IC, IER1

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
      LON = 0D0
      LAT = 0D0
      IRING = 0

* Check if given NRING is valid.
      IF ( NRING < 4 .OR. NRING > NRINGMAX .OR. NRING/2*2 /= NRING )
     &     THEN
         IER = -1
         RETURN
      END IF

* Get grid parameters if NRING changed after previous call.
      IF ( NRING /= NRINGL )  THEN

         CALL GRIDPAR ( NRING, NCELL, CA, CB, CL, NCR, NC1R, DL, IER1 )
         IF ( IER1 /= 0 )  THEN
            IER = -1
            RETURN
         END IF

      END IF

      IF ( ICELL < 1 .OR. ICELL > NCELL )  THEN
         IER = -2
         RETURN
      END IF

* Ring number IRING containing cell ICELL.
      IF ( ICELL >= NC1R(NRING) )  THEN
         IRING = NRING
      ELSE
         DO I=1,NRING-1
            IF ( ICELL >= NC1R(I) .AND. ICELL < NC1R(I+1) )  THEN
               IRING = I
               EXIT
            END IF
         END DO
      END IF

* Central latitude of the cell ICELL.
      LAT = CL(IRING)

* Cell number in the IRING-th ring.
      IC = ICELL - NC1R(IRING) + 1

* Central longitude of the cell ICELL.
      LON = DL(IRING) * (DBLE(IC)-0.5D0)

* Save the last given ring number.
      NRINGL = NRING

      END   ! CN2POS2
