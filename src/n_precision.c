
#include <stdio.h>
#include <string.h>

#define   FORTRAN_FUNCTION
#define   MATCHN( s1, s2, n)   (memcmp((s1), (s2), (n)) == 0)
#define   DIFFERN(s1, s2, n)   (memcmp((s1), (s2), (n)) != 0)

#define        INTEGER         long

/*
   ----------------------------------------------------------------
   Define the number of BYTES and IEEE representations of 1 and 1.0
   ----------------------------------------------------------------
*/

/*
 *===============================================*
 *             H_  => HALF   precision           *
 *             S_  => SINGLE precision           *
 *             D_  => DOUBLE precision           *
 *             Q_  => QUAD   precision           *
 *              _I => INTEGER*n                  *
 *              _R => REAL*n                     *
 *===============================================*
*/

#define        H_ILEN          2
#define        S_ILEN          4
#define        D_ILEN          8
#define        S_RLEN          4
#define        D_RLEN          8
#define        Q_RLEN          16

#if defined(LITTLE_ENDIAN) || defined(LITTLE)
#define        H_IONE          "\001\000"
#define        S_IONE          "\001\000\000\000"
#define        D_IONE          "\001\000\000\000\000\000\000\000"
#define        S_RONE          "\000\000\200\077"
#define        D_RONE          "\000\000\000\000\000\000\360\077"
#define        Q_RONE          "\000\000\000\000\000\000\000\000\000\000\000\000\000\000\377\077"
#else
#define        H_IONE          "\000\001"
#define        S_IONE          "\000\000\000\001"
#define        D_IONE          "\000\000\000\000\000\000\000\001"
#define        S_RONE          "\077\200\000\000"
#define        D_RONE          "\077\360\000\000\000\000\000\000"
#if defined(VPP)
#define        Q_RONE          "\077\377\000\000\000\000\000\000\000\000\000\000\000\000\000\000"
#else
#define        Q_RONE          "\077\360\000\000\000\000\000\000\074\260\000\000\000\000\000\000"
#endif
#endif

/*
   --------------------------------------
   Modify for specific systems & machines
   --------------------------------------
*/

#define        NPREC           n_precision_

#ifdef    __uxp__
#define        NPREC           n_precision_
#endif /* __uxp__ */

#ifdef    __sgi
#define        NPREC           n_precision_
#endif /* __sgi */

#ifdef      sun
#define        NPREC           n_precision_
#endif /*   sun */

#ifdef    __hpux
#define        NPREC           n_precision
#endif /* __hpux */

#ifdef    _CRAY
#define        NPREC           N_PRECISION

#ifndef   _CRAYIEEE              /*  Not IEEE, (not Cray T3d or Cray T90)  */

#define        H_IONE          "\000\001"
#define        S_IONE          "\000\000\000\001"
#define        D_IONE          "\000\000\000\000\000\000\000\001"
#define        S_RONE          "\100\001\200\000\000\000\000\000"
#define        D_RONE          "\100\001\200\000\000\000\000\000\000\000\000\000\000\000\000\000"
#define        Q_RONE          "\100\001\200\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000"
#endif /* _CRAYIEEE */
#endif /* _CRAY */

#ifndef ECMWF_PROTOTYPING
#define ECMWF_PROTOTYPING
#ifdef __STDC__

/* Use ANSI C protyping and declaration */

#ifndef   __
#define   __(_S)      _S         /* use ANSI prototype */
#endif /* __ */

#ifndef   ___
#define   ___(_T,_V)  _T _V      /* use ANSI formal argument declarations */
#endif /* ___ */

#else  /* __STDC__ */

/* Use traditional (K & R) C prototyping and declaration */

#ifndef   __
#define   __(_S)      ()         /* use K & R prototype */
#endif /* __ */

#ifndef   ___
#define   ___(_T,_V)  _V         /* use K & R formal argument decalaration */
#endif /* ___ */

#define   const
#define   volatile
#define   void      char
#endif /* __STDC__ */
#endif /* ECMWF_PROTOTYPING */

INTEGER  NPREC __((char *));

/*                        =========================    */
INTEGER FORTRAN_FUNCTION  NPREC(___(char *, value))
/*                        =========================    */

#ifndef __STDC__
                                    char *  value;
#endif /* !__STDC__ */

{
/*
******************************************************************************
**
*  NAME
*       N_PRECISION - Return the no. of bytes in a REAL/DOUBLE_PREC/INTEGER
*
*  SYNOPSIS
*       INTEGER FUNCTION N_PRECISION(VAL_ARRAY)
*
*  STANDARDS
*       ECMWF extension. A FORTRAN-CALLABLE routine.
*
*  DESCRIPTION
*       The N_PRECISION function returns the number of bytes of storage that
*       are occupied by a variable, whose type is given by the values in the
*       VAL_ARRAY. The amount of storage is not necessarily the number of
*       bytes used to represent the quantity. For example, on the C90 the
*       number of bytes of storage for an INTEGER is 8 (64 bits), while the
*       internal representation uses only the lower 46 bits!
*
*       VAL_ARRAY is a 2-element array, both of which are set to the value 1.
*       This value, and the spacing between the 2 values, can then be used to
*       calculate the number of bytes of storage occupied by a variable of the
*       same type.
*
*   EXAMPLE
*
*            PROGRAM testit
*            INTEGER*4 n_precision
*
*            REAL*4            x4(2)
*            REAL*8            x8(2)
*            REAL*16           x16(2)  ! DOUBLE    NOT ALLOWED ON THE T3D
*            REAL              xr(2)
*            DOUBLE PRECISION  xd(2)   ! DOUBLE    NOT ALLOWED ON THE T3D
*            INTEGER*2         i2(2)
*            INTEGER*4         i4(2)
*            INTEGER*8         i8(2)   ! INTEGER*8 NOT ALLOWED ON HP or SUN
*            INTEGER           ii(2)
*
*            DO 10 i = 1, 2
*            x4 (i) = 1.0
*            x8 (i) = 1.0
*            x16(i) = 1.0
*            xr (i) = 1.0
*            xd (i) = 1.0
*            i2 (i) = 1
*            i4 (i) = 1
*            i8 (i) = 1
*            ii (i) = 1
*        10  CONTINUE
*
*            print *, 'n_precision(i2 ) = ', n_precision(i2 )
*            print *, 'n_precision(i4 ) = ', n_precision(i4 )
*            print *, 'n_precision(i8 ) = ', n_precision(i8 )
*            print *, 'n_precision(ii ) = ', n_precision(ii )
*            print *, 'n_precision(x4 ) = ', n_precision(x4 )
*            print *, 'n_precision(x8 ) = ', n_precision(x8 )
*            print *, 'n_precision(x16) = ', n_precision(x16)
*            print *, 'n_precision(xr ) = ', n_precision(xr )
*            print *, 'n_precision(xd ) = ', n_precision(xd )
*
*            print *, 'n_precision(1.0) = ', n_precision(1.0)
*            print *, 'n_precision(0.0) = ', n_precision(0.0)
*
*            STOP
*            END
*
*      On the VPP300, with the default compiler setting, the output is:
*
*            n_precision(i2 ) = 2
*            n_precision(i4 ) = 4
*            n_precision(i8 ) = 8
*            n_precision(ii ) = 4
*            n_precision(x4 ) = 4
*            n_precision(x8 ) = 8
*            n_precision(x16) = 16
*            n_precision(xr ) = 4
*            n_precision(xd ) = 8
*            n_precision(1.0) = 0    # BAD INPUT
*            n_precision(0.0) = 0    # BAD INPUT
*
*      On the SGI, with the default compiler settings, the output is:
*
*            n_precision(i2 ) = 2
*            n_precision(i4 ) = 4
*            n_precision(i8 ) = 8
*            n_precision(ii ) = 4
*            n_precision(x4 ) = 4
*            n_precision(x8 ) = 8
*            n_precision(x16) = 8    # REAL*16 treated as REAL*8
*            n_precision(xr ) = 4
*            n_precision(xd ) = 8
*            n_precision(1.0) = 0    # BAD INPUT
*            n_precision(0.0) = 0    # BAD INPUT
*
*      On the HP, with the default compiler settings, the output is:
*
*            n_precision(i2 ) = 2
*            n_precision(i4 ) = 4
*            n_precision(i8 ) = 4    # HP DOES NOT ALLOW 8-BYTE INTEGERS
*            n_precision(ii ) = 4
*            n_precision(x4 ) = 4
*            n_precision(x8 ) = 8
*            n_precision(x16) = 16
*            n_precision(xr ) = 4
*            n_precision(xd ) = 8
*            n_precision(1.0) = 0    # BAD INPUT
*            n_precision(0.0) = 0    # BAD INPUT
*
*      On the SUN, with the default compiler settings, the output is:
*
*            n_precision(i2 ) = 2
*            n_precision(i4 ) = 4
*            n_precision(i8 ) = 4    # SUN DOES NOT ALLOW 8-BYTE INTEGERS
*            n_precision(ii ) = 4
*            n_precision(x4 ) = 4
*            n_precision(x8 ) = 8
*            n_precision(x16) = 16
*            n_precision(xr ) = 4
*            n_precision(xd ) = 8
*            n_precision(1.0) = 0    # BAD INPUT
*            n_precision(0.0) = 0    # BAD INPUT
*
*      On the Cray C90, with the default compiler settings, the output is:
*
*            n_precision(i2 ) = 8    # NOT 2, Storage is 1 x 64-bit word
*            n_precision(i4 ) = 8    # NOT 4, Storage is 1 x 64-bit word
*            n_precision(i8 ) = 8
*            n_precision(ii ) = 8
*            n_precision(x4 ) = 8    # NOT 4, REAL*4 => REAL*8
*            n_precision(x8 ) = 8
*            n_precision(x16) = 16
*            n_precision(xr ) = 8
*            n_precision(xd ) = 16
*            n_precision(1.0) = 0    # BAD INPUT
*            n_precision(0.0) = 0    # BAD INPUT
*
*      On the Cray T3D, with the default compiler settings, the output is:
*
*            n_precision(i2 ) = 8    # NOT 2, Storage is 1 x 64-bit word
*            n_precision(i4 ) = 8    # NOT 4, Storage is 1 x 64-bit word
*            n_precision(i8 ) = 8
*            n_precision(ii ) = 8
*            n_precision(x4 ) = 8    # NOT 4, REAL*4 => REAL*8
*            n_precision(x8 ) = 8
*            n_precision(x16) = 8    # T3D DOES NOT ALLOW DOUBLE-PRECISION
*            n_precision(xr ) = 8
*            n_precision(xd ) = 8    # T3D DOES NOT ALLOW DOUBLE-PRECISION
*            n_precision(1.0) = 0    # BAD INPUT
*            n_precision(0.0) = 0    # BAD INPUT
*
*  RETURN VALUE
*      Successful completion is indicated by a non-zero value.
*
*      0  -  The input array was probably incorrect, i.e. it was not a
*            2-element array, or one or both of the elements were not set to
*            the value 1. (See above example)
*
*      n  -  The number of BYTEs that the array element, and hence the
*            variable TYPE, occupies. This can be one of: 2, 4, 8, 16.
*
*  SEE ALSO
*      NR_COMPAT_EC, ND_COMPAT_EC, NI_COMPAT_EC
**
******************************************************************************
*/

/* =========  DEBUG ========
   unsigned char  *val = (unsigned char *) value;
   int             i, j;

   printf("INPUT VALUE=");
   for (i = 0; i < 16 ; i++)
   {
      j = val[i];
      printf("/%3.3o", j);
   }
   printf("\n");
   ========== DEBUG ======== */

/*
   ----------------------
   Half-precision INTEGER
   ----------------------
*/

   if (MATCHN(value, H_IONE, H_ILEN))
   {
      if (MATCHN(value + H_ILEN, H_IONE, H_ILEN))
         return H_ILEN;
   }

/*
   ------------------------
   Single-precision INTEGER
   ------------------------
*/

   if (MATCHN(value, S_IONE, S_ILEN))
   {
      if (MATCHN(value + S_ILEN, S_IONE, S_ILEN))
         return S_ILEN;
   }

/*
   ---------------------
   Single-precision REAL
   ---------------------
*/

   if (MATCHN(value, S_RONE, S_RLEN))
   {
      if (MATCHN(value + S_RLEN, S_RONE, S_RLEN))
         return S_RLEN;
   }

/*
   ------------------------
   Double-precision INTEGER
   ------------------------
*/

   if (MATCHN(value, D_IONE, D_ILEN))
   {
      if (MATCHN(value + D_ILEN, D_IONE, D_ILEN))
         return D_ILEN;
   }

/*
   ---------------------
   Double-precision REAL
   ---------------------
*/

   if (MATCHN(value, D_RONE, D_RLEN))
   {
      if (MATCHN(value + D_RLEN, D_RONE, D_RLEN))
         return D_RLEN;
   }

/*
   ------------------------
   Quadruple-precision REAL
   ------------------------
*/

   if (MATCHN(value, Q_RONE, Q_RLEN))
   {
      if (MATCHN(value + Q_RLEN, Q_RONE, Q_RLEN))
         return Q_RLEN;
   }

/*
   ------------------=============
   NONE OF THE ABOVE - ERROR ERROR
   ------------------=============
*/

   return 0;
}

