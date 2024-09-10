      module spher_harm
      contains 
      SUBROUTINE SPHAR(NMAX, GZON, GLAT, SH)

C***********************************************************************
C
C
C     This subroutine computes the surface spherical harmonic functions
C     (V and W), their latitudinal derivatives (V_N and W_N) and their
C     zonal derivatives (V_E and W_E) as a function of the geocentric
C     latitude GLAT and the geocentric zonal coordinate GZON.  GZON may   
C     be longitude, Right Ascension (RA), Local Solar Time (LST) or any
C     other zonal coordinate, depending on the application.  SPHAR uses
C     the convention adopted for gravity models which does NOT use the
C     Condon-Shortley phase.  Therefore, the even surface spherical
C     harmonic functions are always positive at the northern end of the
C     Greenwich (prime) meridian.  This subroutine also uses the
C     convention of normalizing the orthogonality conditions to 4 * PI
C     instead of unity.  The user must specify the degree NMAX of the
C     truncation needed.  
C
C               INPUTS
C
C     NMAX    : Degree of truncation for spherical harmonic series
C     GZON    : Geocentric (east) longitude (or any zonal coordinate)
C               in radians
C     GLAT    : Geocentric latitude in radians
C
C               OUTPUTS
C
C     V       : Array for even surface spherical harmonic functions
C     W       : Array for odd surface spherical harmonic functions
C     V_E     : Array for partial derivatives dV/dE
C     W_E     : Array for partial derivatives dW/dE
C     V_N     : Array for partial derivatives dV/dN
C     W_N     : Array for partial derivatives dW/dN
C
C     These output arrays contain the values associated with the input
C     position.
C
C     E is normalized zonal coordinate measured eastward (radians)
C     N is latitude direction measured northward (radians)
C
C     Normalized zonal coordinate means the zonal coordinate is divided 
C     by DCOS(GLAT) so that the derivatives are evaluated per radian 
C     along a great circle tangent to the local latitude circle in the 
C     easterly direction.
C
C     References:  This algorithm was developed through combining
C     the recursion formulas in TP SCC 008 "Mathematical Foundation
C     for SCC Astrodynamic Theory" (1982) (pp 105-108) with the  
C     recursion formulas in "Satellite Orbits - Models, Methods and
C     Applications", (2000) by Oliver Montenbruck and Eberhard Gill
C     (pp 66-67).
C
C
C     Written by:    Mark Storz
C
C     Last Updated:  2001 April 27
C
C***********************************************************************

      IMPLICIT NONE   
            
C     Double precision declarations for arrays
      REAL*8 SH(1:(NMAX+1)*(NMAX+1))
      REAL*8 V(0:NMAX,0:NMAX)     ! Even surface spherical harmonics
      REAL*8 W(NMAX+1,NMAX+1)         ! Odd surface spherical harmonics
      REAL*8 V_E(0:NMAX,0:NMAX)   ! Partial derivatives dV/dE
      REAL*8 W_E(NMAX,NMAX)       ! Partial derivatives dW/dE
      REAL*8 V_N(0:NMAX,0:NMAX)   ! Partial derivatives dV/dN
      REAL*8 W_N(NMAX,NMAX)       ! Partial derivatives dW/dN
      REAL*8 V0(-1:NMAX)       ! Zeroth order spherical harmonics
      REAL*8 V0_N(0:NMAX)      ! Partial derivatives dV0/dN
      REAL*8 VOP(0:NMAX+1,NMAX+1) ! Ratio V over P (axial distance)
      REAL*8 WOP(0:NMAX+1,NMAX+1) ! Ratio W over P (axial distance)

C     Double precision declarations for scalars

      REAL*8 CLON       ! Cosine of GZON 
      REAL*8 EM         ! M (Order of spherical harmonic term)
      REAL*8 EMAM       ! 2M
      REAL*8 EMAMA1     ! 2M + 1
      REAL*8 EMAMS1     ! 2M - 1
      REAL*8 EMS1       ! M - N
      REAL*8 EN         ! N (Degree of spherical harmonic term)
      REAL*8 ENAM       ! N + M
      REAL*8 ENAMS1     ! N + M - 1
      REAL*8 ENANA1     ! N + N + 1
      REAL*8 ENANS1     ! N + N - 1
      REAL*8 ENSM       ! N - M
      REAL*8 FACTOR     ! Normalization factor
      REAL*8 GLAT       ! Geocentric latitude (rad)
      REAL*8 GZON       ! Geocentric longitude or any zonal coord. (rad)
      REAL*8 P          ! Axial distance to the surface of unit sphere
      REAL*8 SECTOR     ! Sectorial factor
      REAL*8 SLON       ! Sine of longitude or any zonal coord. (rad)
      REAL*8 TESSER     ! Tesseral factor
      REAL*8 X          ! Cartesian X-coord. on surface of unit sphere 
      REAL*8 Y          ! Cartesian Y-coord. on surface of unit sphere
      REAL*8 Z          ! Cartesian Z-ccord. on surface of unit sphere          
      
C     Integer declarations

      INTEGER M         ! Order of spherical harmonic term
      INTEGER MAM       ! M + M
      INTEGER N         ! N
      INTEGER NAM       ! N + M
      INTEGER NAN       ! N + N
      INTEGER NMAX      ! Order of spherical harmonic truncation
      INTEGER NI
      INTEGER IFUNC
      INTEGER MI
C     V0 contains the zeroth order spherical harmonic functions.
C     V0_N contains the latitude derivatives of V0.
C     VOP contains the even spherical harmonic functions divided by P.
C     WOP contains the odd spherical harmonic functions divided by P.

C     Compute P: axial distance of a point on surface of unit sphere
C     Compute Z: coordinate on the unit sphere perpendicular to the
C     equatorial plane measured parallel to the Earth's axis in the 
C     northward direction

      P = DCOS(GLAT)
      Z = DSIN(GLAT)

      CLON = DCOS(GZON) ! Compute cosine of zonal coordinate.
      SLON = DSIN(GZON) ! Compute sine of zonal coordinate.

C     Compute X:  Coordinate on unit sphere perpendicular to the +/-90 
C     meridian plane (positive direction is toward 0 deg longitude).
C     Compute Y:  Coordinate on unit sphere perpendicular to the prime 
C     meridian plane (positive direction is toward +90 deg longitude).
C     Longitude may be replaced by any zonal coordinate.

      X = P * CLON
      Y = P * SLON

C     Load the values for zeroth degree vectors and arrays.

      V0(0)    = 1.D0
      V0_N(0)  = 0.D0
      V(0,0)   = 1.D0
      V_E(0,0) = 0.D0
      V_N(0,0) = 0.D0

C     Set the value for fictitious -1 degree term of zeroth order to 0.
C     This is needed to 'jump start' the recursion for V0 values.

      V0(-1) = 0.D0

C     Load the (1,1) values for the ratio V/P and W/P, where P is the
C     distance of point on the unit sphere from the axis.  P represents
C     the Greek letter "rho" which is a common symbol for this          
C     cylindrical coordinate.

      VOP(1,1) = CLON
      WOP(1,1) = SLON

C     Initialize the sectorial factor to 2.

      SECTOR = 2.D0

C     The following are two nested LOOPs that evaluate all recursions.

      DO M = 1,NMAX

C       Compute the single precision index factors for zeroth order
C       terms and sectorial terms.  In constructing these variable
C       names, "A" means "Add" and "S" means "Subtract".

        MAM    = M + M

        EM     = DBLE(M)       ! M
        EMAM   = DBLE(MAM)     ! 2M
        EMS1   = DBLE(M-1)     ! M-1
        EMAMA1 = DBLE(MAM+1)   ! 2M+1
        EMAMS1 = DBLE(MAM-1)   ! 2M-1

C       Evaluate recursions for the zeroth order V and dV/dN vectors.

        V0(M) = (EMAMS1 * Z * V0(M-1) - EMS1 * V0(M-2)) / EM
        V0_N(M) = Z * V0_N(M-1) + EM * P * V0(M-1)

C       Evaluate normalization factor (FACTOR) for zeroth order
C       spherical harmonics and use it to produce the normalized output
C       arrays for V and dV/dN.  dV/dE is always zero for the zeroth
C       order (zonal) spherical harmonics.

        FACTOR = DSQRT(EMAMA1)
        V(M, 0) = V0(M) * FACTOR
        V_E(M, 0) = 0.D0
        V_N(M, 0) = V0_N(M) * FACTOR

C       Set the fictitious sub-sectorial values for V/P and W/P to 0.
C       These are needed to 'jump start' the recursion for V/P and W/P.

        VOP(M-1, M) = 0.D0
        WOP(M-1, M) = 0.D0

C       Evaluate normalization factor (FACTOR) for sectorial spherical
C       harmonics and use it to produce the normalized output arrays
C       for V, W, V_N, W_N, V_E and W_E corresponding to the sectorial
C       spherical harmonic functions.  Begin by computing sectorial
C       factor SECTOR and by initializing TESSER to SECTOR.

        SECTOR =  SECTOR / (EMAMS1 * EMAM)
        TESSER =  SECTOR
        FACTOR   =  DSQRT(SECTOR * EMAMA1)
        V(M,M)   =  P * VOP(M,M) * FACTOR
        W(M,M)   =  P * WOP(M,M) * FACTOR
        V_E(M,M) = -EM * WOP(M, M) * FACTOR
        W_E(M,M) =  EM * VOP(M, M) * FACTOR
        V_N(M,M) = (EMAM * VOP(M-1,M) - EM * Z * VOP(M,M)) * FACTOR
        W_N(M,M) = (EMAM * WOP(M-1,M) - EM * Z * WOP(M,M)) * FACTOR

        DO N = M+1,NMAX

C         Compute the single precision index factors for tesseral terms.
C         In constructing these variable names, "A" means "Add" and "S"
C         means "Subtract".

          NAN    = N + N
          NAM    = N + M

          EN     = DBLE(N)       ! N
          ENAM   = DBLE(NAM)     ! N+M
          ENSM   = DBLE(N-M)     ! N-M
          ENAMS1 = DBLE(NAM-1)   ! N+M-1
          ENANA1 = DBLE(NAN+1)   ! 2N+1
          ENANS1 = DBLE(NAN-1)   ! 2N-1

C         Evaluate recursions for the tesseral V/P and W/P terms.

          VOP(N,M) = (ENANS1*Z*VOP(N-1,M) - ENAMS1*VOP(N-2,M)) / ENSM
          WOP(N,M) = (ENANS1*Z*WOP(N-1,M) - ENAMS1*WOP(N-2,M)) / ENSM

C         Evaluate the tesseral factor and use it to evaluate the
C         normalization factor (FACTOR) for tesseral spherical
C         harmonics.  Then use FACTOR to produce the normalized output
C         arrays for V, W, V_N, W_N, V_E and W_E corresponding to the
C         tesseral spherical harmonic functions.

          TESSER   =  TESSER * ENSM / ENAM
          FACTOR   =  DSQRT(TESSER * ENANA1)
          V(N,M)   =  P * VOP(N,M) * FACTOR
          W(N,M)   =  P * WOP(N,M) * FACTOR
          V_E(N,M) = -EM * WOP(N,M) * FACTOR
          W_E(N,M) =  EM * VOP(N,M) * FACTOR
          V_N(N,M) = (ENAM * VOP(N-1,M) - EN * Z * VOP(N,M)) * FACTOR
          W_N(N,M) = (ENAM * WOP(N-1,M) - EN * Z * WOP(N,M)) * FACTOR

        END DO

C       Evaluate the recursions for the sectorial V/P and W/P terms.

        VOP(M+1,M+1) = EMAMA1 * (X * VOP(M,M) - Y * WOP(M,M))
        WOP(M+1,M+1) = EMAMA1 * (X * WOP(M,M) + Y * VOP(M,M))

      END DO

      ! Re-normalize from geodesy to fully-normalized convention
      ! and load SH in order V00,V10,V11,W11,V20,V21,W21,V22,W22,...
      IFUNC = 1
      
      DO NI=0, NMAX
        ! print *, 'NI = ', NI
        SH(IFUNC) = (V(NI,0) * 0.5*DSQRT(DBLE(2.0))); ! MI=0
        IFUNC = IFUNC + 1
        
        DO MI = 1,NI
          ! print *, 'MI = ', MI
          SH(IFUNC) = (V(NI,MI) * 0.5);
          IFUNC = IFUNC + 1
          SH(IFUNC) = (W(NI,MI) * 0.5);
          ! print *, 'W', NI, MI
          IFUNC = IFUNC + 1
          
        ENDDO
      ENDDO

      RETURN
      END
      end module
