******************************************************************************
*                              Function NC2NR
*                             ~~~~~~~~~~~~~~~~
*
*  Purpose:
*
*     Compute the nearest number of rings for a given number of cells for
*     the SREAG pixelization
*
*  Usage:
*
*     NC2NR ( NCELL)
*
*  Input arguments:
*
*     NCELL - number of cells, INTEGER(4);
*             must be in the range [20:2147421180]
*
*  Returned value:
*
*     The nearest even number of rings corresponding to NCELL, INTEGER(4)
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
*     1. If NCELL < 20 or NCELL > 2147421180, zero value is returned.
*
*
*                                                 (C) Zinovy Malkin, 2019-2020
******************************************************************************

      FUNCTION NC2NR ( NCELL )

      IMPLICIT NONE

* Function

      INTEGER(4)  NC2NR

* Input variables

      INTEGER(4)  NCELL             ! number of cells

* Local variables

      REAL(8)  NR1

* ----------------------------------------------------------------------------

      IF ( NCELL >= 20 .AND. NCELL <= 2147421180 ) THEN
* Compute the real number of rings.
         NR1 = DSQRT(DBLE(NCELL)) * 0.8862269225D0
* Round to the nearest even number and convert to integer.
         NC2NR = NINT( DNINT(NR1/2D0) * 2D0 )
      ELSE
         NC2NR = 0
      END IF

      END   ! NC2NR
