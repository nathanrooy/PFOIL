C--- THIS SINGLE SUBROUTINE RUNS EVERYTHING WE NEED ----------------------------

      
      SUBROUTINE PFOIL_VISC(X_IN,Y_IN,ADEG_IN,RE_IN,N_IN,ITMAX_IN,
     & VERB,CL_OUT, CM_OUT, CD_OUT, CDF_OUT, CDP_OUT)

Cf2py intent(in) X_IN, Y_IN, ADEG_IN, RE_IN, N_IN , ITMAX_IN, VERB
Cf2py intent(out) CL_OUT, CM_OUT, CD_OUT, CDF_OUT, CDP_OUT
      
      INCLUDE 'XFOIL.INC'
      INTEGER N_IN
      LOGICAL VERB
      REAL X_IN(IBX), Y_IN(IBX)
      REAL ADEG_IN, RE_IN
      
      IF(VERB) THEN
       WRITE(*,*) '> PFOIL VERBOSE MODE: ON'
      ENDIF

C     max panel angle threshold for warning
      DATA ANGTOL / 40.0 /
      
      CALL INIT
      
C     transfer input variables over to xfoil
      NB     = N_IN
      REINF1 = RE_IN
      ADEG   = ADEG_IN
      ITMAX  = ITMAX_IN
      
      DO I=1, IBX
       XB(I) = X_IN(I)
       YB(I) = Y_IN(I)
      ENDDO
      
C     run remaining subroutines      
      CALL LOAD
      CALL ABCOPY(.TRUE.)
      CALL CANG(X,Y,N,0, IMAX,AMAX)

      IF(ABS(AMAX).GT.ANGTOL .AND. VERB)THEN
         WRITE(*,1081) AMAX, IMAX
 1081    FORMAT(
     &  /' WARNING: Poor input coordinate distribution'
     &  /'          Excessive panel angle', F7.1,'  at i =', I4
     &  /'          Repaneling with PANE and/or PPAR suggested'
     &  /'           (doing GDES,CADD before repaneling _may_'
     &  /'            improve excessively coarse LE spacing' )
      ENDIF
      
      IF(VERB) THEN
        CALL MRSHOW(.TRUE.,.TRUE.)
      ENDIF
      
      LALFA  = .TRUE.
      LVISC  = .TRUE.
      LVCONV = .FALSE.
      
      ALFA = DTOR * ADEG
      QINF = 1.0
     
      CALL SPECAL
      
      IF(ABS(ALFA-AWAKE) .GT. 1.0E-5) LWAKE  = .FALSE.
      IF(ABS(ALFA-AVISC) .GT. 1.0E-5) LVCONV = .FALSE.
      IF(ABS(MINF-MVISC) .GT. 1.0E-5) LVCONV = .FALSE.
      
      CALL VISCAL(ITMAX)
C      CALL CPX
C      CALL FCPMIN
  
      CDPDIF = CD - CDF
     
C     transfer final values for output
      CL_OUT = CL
      CM_OUT = CM
      CD_OUT = CD
      CDF_OUT = CDF
      CDP_OUT = CDPDIF  
      
      WRITE(*,*) '> FINISHED!'
      
      END
      
C--- END -----------------------------------------------------------------------

C ### AREAD.f ##################################################################

      SUBROUTINE AREAD(LU,FNAME,NMAX,X,Y,N,NAME,ISPARS,ITYPE,INFO)
      
      DIMENSION X(NMAX), Y(NMAX)
      CHARACTER*(*) FNAME
      CHARACTER*(*) NAME
      CHARACTER*(*) ISPARS
   
C--------------------------------------------------------
C     Reads in several types of airfoil coordinate file.
C
C  Input:
C       LU      logical unit to use for reading
C       FNAME   name of coordinate file to be read,
C               if FNAME(1:1).eq.' ', unit LU is assumed 
C               to be already open
C       INFO   0 keep quiet
C              1 print info on airfoil
C  Output:
C       X,Y     coordinates
C       N       number of X,Y coordinates
C       NAME    character name string        (if ITYPE > 1)
C       ISPARS  ISES/MSES domain-size string (if ITYPE > 2)
C       ITYPE returns type of file:
C           0  None.  Read error occurred.
C           1  Generic.
C           2  Labeled generic.
C           3  MSES single element.
C           4  MSES multi-element.
C--------------------------------------------------------
      CHARACTER*80 LINE1,LINE2,LINE
      LOGICAL LOPEN, ERROR
      DIMENSION A(10)
C
      IEL = 0
      NEL = 0
C
C---- assume read error will occur
      ITYPE = 0
C      
      LOPEN = FNAME(1:1) .NE. ' '
      IF(LOPEN) OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=98)
C
 11   READ(LU,1000,END=99,ERR=98) LINE1
      IF(INDEX('#!',LINE1(1:1)) .NE. 0) GO TO 11
C
 12   READ(LU,1000,END=99) LINE2
      IF(INDEX('#!',LINE2(1:1)) .NE. 0) GO TO 12
C
      I = 1
C
C---- try to read two numbers from first line
      NA = 2
      CALL GETFLT(LINE1,A,NA,ERROR)
      IF(ERROR .OR. NA.LT.2) THEN
C------ must be a name string
        NAME = LINE1
      ELSE
C------ no name, just two valid numbers... must be plain airfoil file
        NAME = ' '
        IF(INFO.GT.0) THEN
         WRITE(*,*)
         WRITE(*,*) 'Plain airfoil file'
        ENDIF
        ITYPE = 1
        REWIND(LU)
        GO TO 50
      ENDIF
C
C---- if we got here, there's a name line,
C-    so now try to read four MSES domain numbers from second line
      NA = 4
      CALL GETFLT(LINE2,A,NA,ERROR)
      IF(ERROR .OR. NA.LT.2) THEN
C------ less than two valid numbers... not a valid format
        GO TO 99
C
      ELSEIF(NA.LT.4) THEN
C------ less than four numbers... usual .dat labeled file
        NAME = LINE1
        IF(INFO.GT.0) THEN
         WRITE(*,*)
         WRITE(*,*) 'Labeled airfoil file.  Name:  ', NAME
        ENDIF
        ITYPE = 2
        REWIND(LU)
        READ(LU,1000,END=99) LINE1
        GO TO 50
C
      ELSE
C------ four or more numbers... MSES or MISES file
        IF(INFO.GT.0) THEN
         WRITE(*,*)
         WRITE(*,*) 'MSES airfoil file.  Name:  ', NAME
        ENDIF
        ITYPE = 3
        ISPARS = LINE2
      ENDIF
C
C---- read each element until 999.0 or end of file is encountered
   50 NEL = NEL + 1
      DO 55 I=1, NMAX
 51     READ(LU,1000,END=60) LINE
C
C------ skip comment line
        IF(INDEX('#!',LINE(1:1)) .NE. 0) GO TO 51
C
        NA = 2
        CALL GETFLT(LINE,A,NA,ERROR)
        IF(ERROR) GO TO 99
C
C------ skip line without at least two numbers
        IF(NA.LT.2) GO TO 51
C
        X(I) = A(1)
        Y(I) = A(2)
C
        IF (X(I) .EQ. 999.0 .AND. Y(I) .EQ. 999.0) THEN
C-------- if this is the element we want, just exit
          IF(IEL .EQ. NEL) GO TO 60
C
          IF(IEL.EQ.0) THEN
           CALL ASKI('Enter element number^',IEL)
           ITYPE = 4
          ENDIF
C
C-------- if this is the specified element, exit.
          IF(IEL .EQ. NEL) GO TO 60
          GO TO 50
        ENDIF
   55 CONTINUE
      WRITE(*,5030) NMAX
      WRITE(*,5900)
      IF(LOPEN) CLOSE(LU)
      ITYPE = 0
      RETURN
C
   60 N = I-1
      IF(LOPEN) CLOSE(LU)
      RETURN
C
   98 CONTINUE
      NFN = INDEX(FNAME,' ') + 1
      WRITE(*,5050) FNAME(1:NFN)
      WRITE(*,5900)
      ITYPE = 0
      RETURN
C
   99 CONTINUE
      IF(LOPEN) CLOSE(LU)
      WRITE(*,5100)
      WRITE(*,5900)
      ITYPE = 0
      RETURN
C...............................................................
 1000 FORMAT(A)
 5030 FORMAT(/' Buffer array size exceeded'
     &       /' Maximum number of points: ', I4 )
 5050 FORMAT(/' File OPEN error.  Nonexistent file:  ', A)
 5100 FORMAT(/' File READ error.  Unrecognizable file format')
 5900 FORMAT( ' *** LOAD NOT COMPLETED ***' )
      END ! AREAD
      
      
C***********************************************************************
C    Module:  spline.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************


      FUNCTION SEVAL(SS,X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)
C--------------------------------------------------
C     Calculates X(SS)                             |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      SEVAL = T*X(I) + (1.0-T)*X(I-1) + (T-T*T)*((1.0-T)*CX1 - T*CX2)
      RETURN
      END ! SEVAL
      

      SUBROUTINE SCALC(X,Y,S,N)
      DIMENSION X(N), Y(N), S(N)
C----------------------------------------
C     Calculates the arc length array S  |
C     for a 2-D array of points (X,Y).   |
C----------------------------------------
C
      S(1) = 0.
      DO 10 I=2, N
        S(I) = S(I-1) + SQRT((X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2)
   10 CONTINUE
C
      RETURN
      END ! SCALC


      SUBROUTINE SEGSPL(X,XS,S,N)
C-----------------------------------------------
C     Splines X(S) array just like SPLINE,      |
C     but allows derivative discontinuities     |
C     at segment joints.  Segment joints are    |
C     defined by identical successive S values. |
C-----------------------------------------------
      DIMENSION X(N), XS(N), S(N)
C
      IF(S(1).EQ.S(2)  ) STOP 'SEGSPL:  First input point duplicated'
      IF(S(N).EQ.S(N-1)) STOP 'SEGSPL:  Last  input point duplicated'
C
      ISEG0 = 1
      DO 10 ISEG=2, N-2
        IF(S(ISEG).EQ.S(ISEG+1)) THEN
         NSEG = ISEG - ISEG0 + 1
         CALL SPLIND(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG,-999.0,-999.0)
         ISEG0 = ISEG+1
        ENDIF
   10 CONTINUE
C
      NSEG = N - ISEG0 + 1
      CALL SPLIND(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG,-999.0,-999.0)
C
      RETURN
      END ! SEGSPL
      

      SUBROUTINE SPLIND(X,XS,S,N,XS1,XS2)
      DIMENSION X(N),XS(N),S(N)
      PARAMETER (NMAX=1000)
      DIMENSION  A(NMAX),B(NMAX),C(NMAX)
C-------------------------------------------------------
C     Calculates spline coefficients for X(S).          |
C     Specified 1st derivative and/or usual zero 2nd    |
C     derivative end conditions are used.               |
C     To evaluate the spline at some value of S,        |
C     use SEVAL and/or DEVAL.                           |
C                                                       |
C     S        independent variable array (input)       |
C     X        dependent variable array   (input)       |
C     XS       dX/dS array                (calculated)  |
C     N        number of points           (input)       |
C     XS1,XS2  endpoint derivatives       (input)       |
C              If = 999.0, then usual zero second       |
C              derivative end condition(s) are used     |
C              If = -999.0, then zero third             |
C              derivative end condition(s) are used     |
C                                                       |
C-------------------------------------------------------
      IF(N.GT.NMAX) STOP 'SPLIND: array overflow, increase NMAX'
C     
      DO 1 I=2, N-1
        DSM = S(I) - S(I-1)
        DSP = S(I+1) - S(I)
        B(I) = DSP
        A(I) = 2.0*(DSM+DSP)
        C(I) = DSM
        XS(I) = 3.0*((X(I+1)-X(I))*DSM/DSP + (X(I)-X(I-1))*DSP/DSM)
    1 CONTINUE
C
      IF(XS1.EQ.999.0) THEN
C----- set zero second derivative end condition
       A(1) = 2.0
       C(1) = 1.0
       XS(1) = 3.0*(X(2)-X(1)) / (S(2)-S(1))
      ELSE IF(XS1.EQ.-999.0) THEN
C----- set zero third derivative end condition
       A(1) = 1.0
       C(1) = 1.0
       XS(1) = 2.0*(X(2)-X(1)) / (S(2)-S(1))
      ELSE
C----- set specified first derivative end condition
       A(1) = 1.0
       C(1) = 0.
       XS(1) = XS1
      ENDIF
C
      IF(XS2.EQ.999.0) THEN
       B(N) = 1.0
       A(N) = 2.0
       XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ELSE IF(XS2.EQ.-999.0) THEN
       B(N) = 1.0
       A(N) = 1.0
       XS(N) = 2.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ELSE
       A(N) = 1.0
       B(N) = 0.
       XS(N) = XS2
      ENDIF
C
      IF(N.EQ.2 .AND. XS1.EQ.-999.0 .AND. XS2.EQ.-999.0) THEN
       B(N) = 1.0
       A(N) = 2.0
       XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ENDIF
C
C---- solve for derivative array XS
      CALL TRISOL(A,B,C,XS,N)
C
      RETURN
      END ! SPLIND
      
      
      SUBROUTINE TRISOL(A,B,C,D,KK)
      DIMENSION A(KK),B(KK),C(KK),D(KK)
C-----------------------------------------
C     Solves KK long, tri-diagonal system |
C                                         |
C             A C          D              |
C             B A C        D              |
C               B A .      .              |
C                 . . C    .              |
C                   B A    D              |
C                                         |
C     The righthand side D is replaced by |
C     the solution.  A, C are destroyed.  |
C-----------------------------------------
C
      DO 1 K=2, KK
        KM = K-1
        C(KM) = C(KM) / A(KM)
        D(KM) = D(KM) / A(KM)
        A(K) = A(K) - B(K)*C(KM)
        D(K) = D(K) - B(K)*D(KM)
    1 CONTINUE
C
      D(KK) = D(KK)/A(KK)
C
      DO 2 K=KK-1, 1, -1
        D(K) = D(K) - C(K)*D(K+1)
    2 CONTINUE
C
      RETURN
      END ! TRISOL
      
      
      SUBROUTINE SINVRT(SI,XI,X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)
C-------------------------------------------------------
C     Calculates the "inverse" spline function S(X).    |
C     Since S(X) can be multi-valued or not defined,    |
C     this is not a "black-box" routine.  The calling   |
C     program must pass via SI a sufficiently good      |
C     initial guess for S(XI).                          |
C                                                       |
C     XI      specified X value       (input)           |
C     SI      calculated S(XI) value  (input,output)    |
C     X,XS,S  usual spline arrays     (input)           |
C                                                       |
C-------------------------------------------------------
C
      SISAV = SI
C
      DO 10 ITER=1, 10
        RES  = SEVAL(SI,X,XS,S,N) - XI
        RESP = DEVAL(SI,X,XS,S,N)
        DS = -RES/RESP
        SI = SI + DS
        IF(ABS(DS/(S(N)-S(1))) .LT. 1.0E-5) RETURN
   10 CONTINUE
      WRITE(*,*)
     &  'SINVRT: spline inversion failed. Input value returned.'
      SI = SISAV
C
      RETURN
      END ! SINVRT
      

      FUNCTION DEVAL(SS,X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)
C--------------------------------------------------
C     Calculates dX/dS(SS)                         |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      DEVAL = X(I) - X(I-1) + (1.-4.0*T+3.0*T*T)*CX1 + T*(3.0*T-2.)*CX2
      DEVAL = DEVAL/DS
      RETURN
      END ! DEVAL
      
      
      FUNCTION CURV(SS,X,XS,Y,YS,S,N)
      DIMENSION X(N), XS(N), Y(N), YS(N), S(N)
C-----------------------------------------------
C     Calculates curvature of splined 2-D curve |
C     at S = SS                                 |
C                                               |
C     S        arc length array of curve        |
C     X, Y     coordinate arrays of curve       |
C     XS,YS    derivative arrays                |
C              (calculated earlier by SPLINE)   |
C-----------------------------------------------
C     
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
C
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      XD = X(I) - X(I-1) + (1.0-4.0*T+3.0*T*T)*CX1 + T*(3.0*T-2.0)*CX2
      XDD = (6.0*T-4.0)*CX1 + (6.0*T-2.0)*CX2
C
      CY1 = DS*YS(I-1) - Y(I) + Y(I-1)
      CY2 = DS*YS(I)   - Y(I) + Y(I-1)
      YD = Y(I) - Y(I-1) + (1.0-4.0*T+3.0*T*T)*CY1 + T*(3.0*T-2.0)*CY2
      YDD = (6.0*T-4.0)*CY1 + (6.0*T-2.0)*CY2
C 
      SD = SQRT(XD*XD + YD*YD)
      SD = MAX(SD,0.001*DS)
C
      CURV = (XD*YDD - YD*XDD) / SD**3
C
      RETURN
      END ! CURV
      
      
      FUNCTION D2VAL(SS,X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)
C--------------------------------------------------
C     Calculates d2X/dS2(SS)                       |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      D2VAL = (6.*T-4.)*CX1 + (6.*T-2.0)*CX2
      D2VAL = D2VAL/DS**2
      RETURN
      END ! D2VAL
      
      
C***********************************************************************
C    Module:  userio.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************
C
C
C==== user input routines with prompting and error trapping
C
C
      SUBROUTINE ASKI(PROMPT,IINPUT)
C
C---- integer input
C
      CHARACTER*(*) PROMPT
      INTEGER IINPUT
      CHARACTER LINE*80
C
      NP = INDEX(PROMPT,'^') - 1
      IF(NP.LE.0) NP = LEN(PROMPT)
C
 10   WRITE(*,1000) PROMPT(1:NP)
C
      READ (*,1001,ERR=10) LINE
      IF(LINE.NE.' ') THEN
        READ (LINE,*,ERR=10) IINPUT
      ENDIF  
      RETURN
C
 1000 FORMAT(/A,'   i>  ',$)
 1001 FORMAT(A)
      END ! ASKI


      SUBROUTINE ASKR(PROMPT,RINPUT)
C
C---- real input
C
      CHARACTER*(*) PROMPT
      REAL RINPUT
      CHARACTER LINE*80
C
      NP = INDEX(PROMPT,'^') - 1
      IF(NP.LE.0) NP = LEN(PROMPT)
C
 10   WRITE(*,1000) PROMPT(1:NP)
C
      READ (*,1001,ERR=10) LINE
      IF(LINE.NE.' ') THEN
        READ (LINE,*,ERR=10) RINPUT
      ENDIF  
      RETURN
C
 1000 FORMAT(/A,'   r>  ',$)
 1001 FORMAT(A)
      END ! ASKR


      SUBROUTINE ASKL(PROMPT,LINPUT)
C
C---- logical input
C
      CHARACTER*(*) PROMPT
      LOGICAL LINPUT
      CHARACTER*1 CHAR
C
      NP = INDEX(PROMPT,'^') - 1
      IF(NP.LE.0) NP = LEN(PROMPT)
C
 10   WRITE(*,1000) PROMPT(1:NP)
      READ (*,1010) CHAR
      IF(CHAR.EQ.'y') CHAR = 'Y'
      IF(CHAR.EQ.'n') CHAR = 'N'
      IF(CHAR.NE.'Y' .AND. CHAR.NE.'N') GO TO 10
C
      LINPUT = CHAR .EQ. 'Y'
      RETURN
C
 1000 FORMAT(/A,' y/n>  ',$)
 1010 FORMAT(A)
      END ! ASKL


      SUBROUTINE ASKS(PROMPT,INPUT)
C
C---- string of arbitrary length input
C
      CHARACTER*(*) PROMPT
      CHARACTER*(*) INPUT
C
      NP = INDEX(PROMPT,'^') - 1
      IF(NP.LE.0) NP = LEN(PROMPT)
C
      WRITE(*,1000) PROMPT(1:NP)
      READ (*,1010) INPUT
C
      RETURN
C
 1000 FORMAT(/A,'   s>  ',$)
 1010 FORMAT(A)
      END ! ASKS


      SUBROUTINE ASKC(PROMPT,COMAND,CARGS)
C
C---- returns 4-byte character string input converted to uppercase
C---- also returns rest of input characters in CARGS string
C
      CHARACTER*(*) PROMPT
      CHARACTER*(*) COMAND, CARGS
C
      CHARACTER*128 LINE
      LOGICAL ERROR
C
      IZERO = ICHAR('0')
C
      NP = INDEX(PROMPT,'^') - 1
      IF(NP.LE.0) NP = LEN(PROMPT)
C
      WRITE(*,1000) PROMPT(1:NP)
      READ (*,1020) LINE
C
C---- strip off leading blanks
      DO K=1, 128
        IF(LINE(1:1) .EQ. ' ') THEN
         LINE = LINE(2:128)
        ELSE
         GO TO 5
        ENDIF
      ENDDO
 5    CONTINUE
C
C---- find position of first blank, "+", "-", ".", ",", or numeral
      K = INDEX(LINE,' ')
      KI = INDEX(LINE,'-')
      IF(KI.NE.0) K = MIN(K,KI)
      KI = INDEX(LINE,'+')
      IF(KI.NE.0) K = MIN(K,KI)
      KI = INDEX(LINE,'.')
      IF(KI.NE.0) K = MIN(K,KI)
      KI = INDEX(LINE,',')
      IF(KI.NE.0) K = MIN(K,KI)
      DO I=0, 9
        KI = INDEX(LINE,CHAR(IZERO+I))
        IF(KI.NE.0) K = MIN(K,KI)
      ENDDO
C
C---- there is no blank between command and argument... use first 4 characters
      IF(K.LE.0) K = 5
C
      IF(K.EQ.1) THEN
C------ the "command" is a number... set entire COMAND string with it
        COMAND = LINE
      ELSE
C------ the "command" is some string... just use the part up to the argument
        COMAND = LINE(1:K-1)
      ENDIF
C
C---- convert it to uppercase
      CALL LC2UC(COMAND)
C
      CARGS = LINE(K:128)
      CALL STRIP(CARGS,NCARGS)
      RETURN
C
 1000 FORMAT(/A,'   c>  ',$)
 1020 FORMAT(A)
      END ! ASKC


      SUBROUTINE LC2UC(INPUT)
      CHARACTER*(*) INPUT
C
      CHARACTER*26 LCASE, UCASE
      DATA LCASE / 'abcdefghijklmnopqrstuvwxyz' /
      DATA UCASE / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
C
      N = LEN(INPUT)
C
      DO 10 I=1, N
        K = INDEX( LCASE , INPUT(I:I) )
        IF(K.GT.0) INPUT(I:I) = UCASE(K:K)
 10   CONTINUE
C
      RETURN
      END ! LC2UC



      SUBROUTINE READI(N,IVAR,ERROR)
      DIMENSION IVAR(N)
      LOGICAL ERROR
C--------------------------------------------------
C     Reads N integer variables, leaving unchanged 
C     if only <return> is entered.
C--------------------------------------------------
      DIMENSION IVTMP(40)
      CHARACTER*80 LINE
C
      READ(*,1000) LINE
 1000 FORMAT(A80)
C
      DO 10 I=1, N
        IVTMP(I) = IVAR(I)
 10   CONTINUE
C
      NTMP = 40
      CALL GETINT(LINE,IVTMP,NTMP,ERROR)
C
      IF(ERROR) RETURN
C
      DO 20 I=1, N
        IVAR(I) = IVTMP(I)
 20   CONTINUE
C
      RETURN
      END ! READI



      SUBROUTINE READR(N,VAR,ERROR)
      DIMENSION VAR(N)
      LOGICAL ERROR
C-------------------------------------------------
C     Reads N real variables, leaving unchanged 
C     if only <return> is entered.
C-------------------------------------------------
      DIMENSION VTMP(40)
      CHARACTER*80 LINE
C
      READ(*,1000) LINE
 1000 FORMAT(A80)
C
      DO 10 I=1, N
        VTMP(I) = VAR(I)
 10   CONTINUE
C
      NTMP = 40
      CALL GETFLT(LINE,VTMP,NTMP,ERROR)
C
      IF(ERROR) RETURN
C
      DO 20 I=1, N
        VAR(I) = VTMP(I)
 20   CONTINUE
C
      RETURN
      END ! READR




      SUBROUTINE GETINT(INPUT,A,N,ERROR)
      CHARACTER*(*) INPUT
      INTEGER A(*)
      LOGICAL ERROR
C----------------------------------------------------------
C     Parses character string INPUT into an array
C     of integer numbers returned in A(1...N)
C
C     Will attempt to extract no more than N numbers, 
C     unless N = 0, in which case all numbers present
C     in INPUT will be extracted.
C
C     N returns how many numbers were actually extracted.
C----------------------------------------------------------
      CHARACTER*130 REC
      CHARACTER*1 TAB
C
      TAB = CHAR(9)
C
C---- only first 128 characters in INPUT will be parsed
      ILEN = MIN( LEN(INPUT) , 128 )
      ILENP = ILEN + 2
C
C---- put input into local work string (which will be munched)
      REC(1:ILENP) = INPUT(1:ILEN) // ' ,'
C
C---- ignore everything after a "!" character
      K = INDEX(REC,'!')
      IF(K.GT.0) REC(1:ILEN) = REC(1:K-1)
C
C---- change tabs to spaces
 5    K = INDEX(REC(1:ILEN),TAB)
      IF(K.GT.0) THEN
       REC(K:K) = ' '
       GO TO 5
      ENDIF
C
      NINP = N
C
C---- count up how many numbers are to be extracted
      N = 0
      K = 1
      DO 10 IPASS=1, ILEN
C------ search for next space or comma starting with current index K
        KSPACE = INDEX(REC(K:ILENP),' ') + K - 1
        KCOMMA = INDEX(REC(K:ILENP),',') + K - 1
C
        IF(K.EQ.KSPACE) THEN
C------- just skip this space
         K = K+1
         GO TO 9
        ENDIF
C
        IF(K.EQ.KCOMMA) THEN
C------- comma found.. increment number count and keep looking
         N = N+1
         K = K+1
         GO TO 9
        ENDIF
C
C------ neither space nor comma found, so we ran into a number...
C-    ...increment number counter and keep looking after next space or comma
        N = N+1
        K = MIN(KSPACE,KCOMMA) + 1
C
  9     IF(K.GE.ILEN) GO TO 11
 10   CONTINUE
C
C---- decide on how many numbers to read, and go ahead and read them
 11   IF(NINP.GT.0) N = MIN( N, NINP )
      READ(REC(1:ILEN),*,ERR=20) (A(I),I=1,N)
      ERROR = .FALSE.
      RETURN
C
C---- bzzzt !!!
 20   CONTINUE
ccc   WRITE(*,*) 'GETINT: String-to-integer conversion error.'
      N = 0
      ERROR = .TRUE.
      RETURN
      END ! GETINT


      SUBROUTINE GETFLT(INPUT,A,N,ERROR)
      CHARACTER*(*) INPUT
      REAL A(*)
      LOGICAL ERROR
C----------------------------------------------------------
C     Parses character string INPUT into an array
C     of real numbers returned in A(1...N)
C
C     Will attempt to extract no more than N numbers, 
C     unless N = 0, in which case all numbers present
C     in INPUT will be extracted.
C
C     N returns how many numbers were actually extracted.
C----------------------------------------------------------
      CHARACTER*130 REC
      CHARACTER*1 TAB
C
      TAB = CHAR(9)
C
C---- only first 128 characters in INPUT will be parsed
      ILEN = MIN( LEN(INPUT) , 128 )
      ILENP = ILEN + 2
C
C---- put input into local work string (which will be munched)
      REC(1:ILENP) = INPUT(1:ILEN) // ' ,'
C
C---- ignore everything after a "!" character
      K = INDEX(REC,'!')
      IF(K.GT.0) REC(1:ILEN) = REC(1:K-1)
C
C---- change tabs to spaces
 5    K = INDEX(REC(1:ILEN),TAB)
      IF(K.GT.0) THEN
       REC(K:K) = ' '
       GO TO 5
      ENDIF
C
      NINP = N
C
C---- count up how many numbers are to be extracted
      N = 0
      K = 1
      DO 10 IPASS=1, ILEN
C------ search for next space or comma starting with current index K
        KSPACE = INDEX(REC(K:ILENP),' ') + K - 1
        KCOMMA = INDEX(REC(K:ILENP),',') + K - 1
C
        IF(K.EQ.KSPACE) THEN
C------- just skip this space
         K = K+1
         GO TO 9
        ENDIF
C
        IF(K.EQ.KCOMMA) THEN
C------- comma found.. increment number count and keep looking
         N = N+1
         K = K+1
         GO TO 9
        ENDIF
C
C------ neither space nor comma found, so we ran into a number...
C-    ...increment number counter and keep looking after next space or comma
        N = N+1
        K = MIN(KSPACE,KCOMMA) + 1
C
  9     IF(K.GE.ILEN) GO TO 11
 10   CONTINUE
C
C---- decide on how many numbers to read, and go ahead and read them
 11   IF(NINP.GT.0) N = MIN( N, NINP )
      READ(REC(1:ILEN),*,ERR=20) (A(I),I=1,N)
      ERROR = .FALSE.
      RETURN
C
C---- bzzzt !!!
 20   CONTINUE
ccc   WRITE(*,*) 'GETFLT: String-to-integer conversion error.'
      N = 0
      ERROR = .TRUE.
      RETURN
      END ! GETFLT



      SUBROUTINE STRIP(STRING,NS)
      CHARACTER*(*) STRING
C----------------------------------------------------
C     Strips leading blanks off STRING and returns 
C     length NS of non-blank part.
C----------------------------------------------------
      NLEN = LEN(STRING)
C
C---- find last non-blank character
      DO K2 = NLEN, 1, -1
        IF(STRING(K2:K2).NE.' ') GO TO 11
      ENDDO
      K2 = 0
   11 CONTINUE
C
C---- find first non-blank character
      DO K1 = 1, K2
        IF(STRING(K1:K1).NE.' ') GO TO 21
      ENDDO
   21 CONTINUE
C
C---- number of non-blank characters
      NS = K2 - K1 + 1
      IF(NS.EQ.0) RETURN
C
C---- shift STRING so first character is non-blank
      STRING(1:NS) = STRING(K1:K2)
C
C---- pad tail of STRING with blanks
      DO K = NS+1, NLEN
        STRING(K:K) = ' '
      ENDDO
C
      RETURN
      END





      SUBROUTINE BSTRIP(STRING,NS)
      CHARACTER*(*) STRING
C--------------------------------------------------
C     Strips all blanks from STRING and returns 
C     length NS of non-blank part.
C     If STRING is all blanks, just returns NS=0
C--------------------------------------------------
C
C---- first remove any leading blanks and get length to be processed
      CALL STRIP(STRING,NS)
C
C---- pass over STRING and strip out all interior blanks
      K = 1
C
 10   CONTINUE
      IF(K.GE.NS) THEN
       RETURN
C
      ELSEIF(STRING(K:K) .EQ. ' ') THEN
       STRING(K:NS-1) = STRING(K+1:NS)
       NS = NS - 1
C
      ELSE
       K = K + 1
C
      ENDIF
C
      GO TO 10
C
      END




      SUBROUTINE GETARG0(IARG,ARG)
C------------------------------------------------
C     Same as GETARG, but...
C
C     ...in the case of Intel Fortran, this one
C     doesn't barf if there's no Unix argument 
C      (just returns blank string instead)
C------------------------------------------------
      CHARACTER*(*) ARG
C
      NARG = IARGC()
      IF(NARG.GE.IARG) THEN
       CALL GETARG(IARG,ARG)
      ELSE
       ARG = ' '
      ENDIF
C
      RETURN
      END ! GETARG0


C***********************************************************************
C    Module:  xbl.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************
C
      SUBROUTINE SETBL
C-------------------------------------------------
C     Sets up the BL Newton system coefficients
C     for the current BL variables and the edge
C     velocities received from SETUP. The local
C     BL system coefficients are then
C     incorporated into the global Newton system.  
C-------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      REAL USAV(IVX,2)
      REAL U1_M(2*IVX), U2_M(2*IVX)
      REAL D1_M(2*IVX), D2_M(2*IVX)
      REAL ULE1_M(2*IVX), ULE2_M(2*IVX)
      REAL UTE1_M(2*IVX), UTE2_M(2*IVX)
      REAL MA_CLMR, MSQ_CLMR, MDI
C
C---- set the CL used to define Mach, Reynolds numbers
      IF(LALFA) THEN
       CLMR = CL
      ELSE
       CLMR = CLSPEC
      ENDIF
C
C---- set current MINF(CL)
      CALL MRCL(CLMR,MA_CLMR,RE_CLMR)
      MSQ_CLMR = 2.0*MINF*MA_CLMR
C
C---- set compressibility parameter TKLAM and derivative TK_MSQ
      CALL COMSET
C
C---- set gas constant (= Cp/Cv)
      GAMBL = GAMMA
      GM1BL = GAMM1
C
C---- set parameters for compressibility correction
      QINFBL = QINF
      TKBL    = TKLAM
      TKBL_MS = TKL_MSQ
C
C---- stagnation density and 1/enthalpy
      RSTBL    = (1.0 + 0.5*GM1BL*MINF**2) ** (1.0/GM1BL)
      RSTBL_MS = 0.5*RSTBL/(1.0 + 0.5*GM1BL*MINF**2)
C
      HSTINV    = GM1BL*(MINF/QINFBL)**2 / (1.0 + 0.5*GM1BL*MINF**2)
      HSTINV_MS = GM1BL*( 1.0/QINFBL)**2 / (1.0 + 0.5*GM1BL*MINF**2)
     &                - 0.5*GM1BL*HSTINV / (1.0 + 0.5*GM1BL*MINF**2)
C
C---- set Reynolds number based on freestream density, velocity, viscosity
      HERAT    = 1.0 - 0.5*QINFBL**2*HSTINV
      HERAT_MS =     - 0.5*QINFBL**2*HSTINV_MS
C
      REYBL    = REINF * SQRT(HERAT**3) * (1.0+HVRAT)/(HERAT+HVRAT)
      REYBL_RE =         SQRT(HERAT**3) * (1.0+HVRAT)/(HERAT+HVRAT)
      REYBL_MS = REYBL * (1.5/HERAT - 1.0/(HERAT+HVRAT))*HERAT_MS
C
      IDAMPV = IDAMP
C
C---- save TE thickness
      DWTE = WGAP(1)
C
      IF(.NOT.LBLINI) THEN
C----- initialize BL by marching with Ue (fudge at separation)
       WRITE(*,*)
       WRITE(*,*) 'Initializing BL ...'
       CALL MRCHUE
       LBLINI = .TRUE.
      ENDIF
C
      WRITE(*,*)
C
C---- march BL with current Ue and Ds to establish transition
      CALL MRCHDU
C
      DO 5 IS=1, 2
        DO 6 IBL=2, NBL(IS)
          USAV(IBL,IS) = UEDG(IBL,IS)
    6   CONTINUE
    5 CONTINUE
C
      CALL UESET
C
      DO 7 IS=1, 2
        DO 8 IBL=2, NBL(IS)
          TEMP = USAV(IBL,IS)
          USAV(IBL,IS) = UEDG(IBL,IS)
          UEDG(IBL,IS) = TEMP
    8   CONTINUE
    7 CONTINUE
C
      ILE1 = IPAN(2,1)
      ILE2 = IPAN(2,2)
      ITE1 = IPAN(IBLTE(1),1)
      ITE2 = IPAN(IBLTE(2),2)
C
      JVTE1 = ISYS(IBLTE(1),1)
      JVTE2 = ISYS(IBLTE(2),2)
C
      DULE1 = UEDG(2,1) - USAV(2,1)
      DULE2 = UEDG(2,2) - USAV(2,2)
C
C---- set LE and TE Ue sensitivities wrt all m values
      DO 10 JS=1, 2
        DO 110 JBL=2, NBL(JS)
          J  = IPAN(JBL,JS)
          JV = ISYS(JBL,JS)
          ULE1_M(JV) = -VTI(       2,1)*VTI(JBL,JS)*DIJ(ILE1,J)
          ULE2_M(JV) = -VTI(       2,2)*VTI(JBL,JS)*DIJ(ILE2,J)
          UTE1_M(JV) = -VTI(IBLTE(1),1)*VTI(JBL,JS)*DIJ(ITE1,J)
          UTE2_M(JV) = -VTI(IBLTE(2),2)*VTI(JBL,JS)*DIJ(ITE2,J)
  110   CONTINUE
   10 CONTINUE
C
      ULE1_A = UINV_A(2,1)
      ULE2_A = UINV_A(2,2)
C
      TINDEX(1) = 0.0
      TINDEX(2) = 0.0
C
C**** Go over each boundary layer/wake
      DO 2000 IS=1, 2
C
C---- there is no station "1" at similarity, so zero everything out
      DO 20 JS=1, 2
        DO 210 JBL=2, NBL(JS)
          JV = ISYS(JBL,JS)
          U1_M(JV) = 0.
          D1_M(JV) = 0.
  210   CONTINUE
   20 CONTINUE
      U1_A = 0.
      D1_A = 0.
C
      DUE1 = 0.
      DDS1 = 0.
C
C---- similarity station pressure gradient parameter  x/u du/dx
      IBL = 2
      BULE = 1.0
C
      AMCRIT = ACRIT(IS)
C
C---- set forced transition arc length position
      CALL XIFSET(IS)
C
      TRAN = .FALSE.
      TURB = .FALSE.
C
C**** Sweep downstream setting up BL equation linearizations
      DO 1000 IBL=2, NBL(IS)
C
      IV  = ISYS(IBL,IS)
C
      SIMI = IBL.EQ.2
      WAKE = IBL.GT.IBLTE(IS)
      TRAN = IBL.EQ.ITRAN(IS)
      TURB = IBL.GT.ITRAN(IS)
C
      I = IPAN(IBL,IS)
C
C---- set primary variables for current station
      XSI = XSSI(IBL,IS)
      IF(IBL.LT.ITRAN(IS)) AMI = CTAU(IBL,IS)
      IF(IBL.GE.ITRAN(IS)) CTI = CTAU(IBL,IS)
      UEI = UEDG(IBL,IS)
      THI = THET(IBL,IS)
      MDI = MASS(IBL,IS)
C
      DSI = MDI/UEI
C
      IF(WAKE) THEN
       IW = IBL - IBLTE(IS)
       DSWAKI = WGAP(IW)
      ELSE
       DSWAKI = 0.
      ENDIF
C
C---- set derivatives of DSI (= D2)
      D2_M2 =  1.0/UEI
      D2_U2 = -DSI/UEI
C
      DO 30 JS=1, 2
        DO 310 JBL=2, NBL(JS)
          J  = IPAN(JBL,JS)
          JV = ISYS(JBL,JS)
          U2_M(JV) = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
          D2_M(JV) = D2_U2*U2_M(JV)
  310   CONTINUE
   30 CONTINUE
      D2_M(IV) = D2_M(IV) + D2_M2
C
      U2_A = UINV_A(IBL,IS)
      D2_A = D2_U2*U2_A
C
C---- "forced" changes due to mismatch between UEDG and USAV=UINV+dij*MASS
      DUE2 = UEDG(IBL,IS) - USAV(IBL,IS)
      DDS2 = D2_U2*DUE2
C
      CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
      CALL BLKIN
C
C---- check for transition and set TRAN, XT, etc. if found
      IF(TRAN) THEN
        CALL TRCHEK
        AMI = AMPL2
      ENDIF
      IF(IBL.EQ.ITRAN(IS) .AND. .NOT.TRAN) THEN
       WRITE(*,*) 'SETBL: Xtr???  n1 n2: ', AMPL1, AMPL2
      ENDIF
C
C---- assemble 10x4 linearized system for dCtau, dTh, dDs, dUe, dXi
C     at the previous "1" station and the current "2" station
C
      IF(IBL.EQ.IBLTE(IS)+1) THEN
C
C----- define quantities at start of wake, adding TE base thickness to Dstar
       TTE = THET(IBLTE(1),1) + THET(IBLTE(2),2)
       DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
       CTE = ( CTAU(IBLTE(1),1)*THET(IBLTE(1),1)
     &       + CTAU(IBLTE(2),2)*THET(IBLTE(2),2) ) / TTE
       CALL TESYS(CTE,TTE,DTE)
C
       TTE_TTE1 = 1.0
       TTE_TTE2 = 1.0
       DTE_MTE1 =               1.0 / UEDG(IBLTE(1),1)
       DTE_UTE1 = -DSTR(IBLTE(1),1) / UEDG(IBLTE(1),1)
       DTE_MTE2 =               1.0 / UEDG(IBLTE(2),2)
       DTE_UTE2 = -DSTR(IBLTE(2),2) / UEDG(IBLTE(2),2)
       CTE_CTE1 = THET(IBLTE(1),1)/TTE
       CTE_CTE2 = THET(IBLTE(2),2)/TTE
       CTE_TTE1 = (CTAU(IBLTE(1),1) - CTE)/TTE
       CTE_TTE2 = (CTAU(IBLTE(2),2) - CTE)/TTE
C
C----- re-define D1 sensitivities wrt m since D1 depends on both TE Ds values
       DO 35 JS=1, 2
         DO 350 JBL=2, NBL(JS)
           J  = IPAN(JBL,JS)
           JV = ISYS(JBL,JS)
           D1_M(JV) = DTE_UTE1*UTE1_M(JV) + DTE_UTE2*UTE2_M(JV)
  350    CONTINUE
   35  CONTINUE
       D1_M(JVTE1) = D1_M(JVTE1) + DTE_MTE1
       D1_M(JVTE2) = D1_M(JVTE2) + DTE_MTE2
C
C----- "forced" changes from  UEDG --- USAV=UINV+dij*MASS  mismatch
       DUE1 = 0.
       DDS1 = DTE_UTE1*(UEDG(IBLTE(1),1) - USAV(IBLTE(1),1))
     &      + DTE_UTE2*(UEDG(IBLTE(2),2) - USAV(IBLTE(2),2))
C
      ELSE
C
       CALL BLSYS
C
      ENDIF
C
C
C---- Save wall shear and equil. max shear coefficient for plotting output
      TAU(IBL,IS) = 0.5*R2*U2*U2*CF2
      DIS(IBL,IS) =     R2*U2*U2*U2*DI2*HS2*0.5
      CTQ(IBL,IS) = CQ2
      DELT(IBL,IS) = DE2
      USLP(IBL,IS) = 1.60/(1.0+US2)
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c      IF(WAKE) THEN
c        ALD = DLCON
c      ELSE
c       ALD = 1.0
c      ENDIF
cC
c      IF(TURB .AND. .NOT.WAKE) THEN
c        GCC = GCCON
c        HKC     = HK2 - 1.0 - GCC/RT2
c        IF(HKC .LT. 0.01) THEN
c         HKC = 0.01
c        ENDIF
c       ELSE
c        HKC = HK2 - 1.0
c       ENDIF
cC
c       HR = HKC     / (GACON*ALD*HK2)
c       UQ = (0.5*CF2 - HR**2) / (GBCON*D2)
cC
c       IF(TURB) THEN
c        IBLP = MIN(IBL+1,NBL(IS))
c        IBLM = MAX(IBL-1,2      )
c        DXSSI = XSSI(IBLP,IS) - XSSI(IBLM,IS)
c        IF(DXXSI.EQ.0.0) DXSSI = 1.0
c        GUXD(IBL,IS) = -LOG(UEDG(IBLP,IS)/UEDG(IBLM,IS)) / DXSSI
c        GUXQ(IBL,IS) = -UQ
c       ELSE
c        GUXD(IBL,IS) = 0.0
c        GUXQ(IBL,IS) = 0.0
c       ENDIF
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C---- set XI sensitivities wrt LE Ue changes
      IF(IS.EQ.1) THEN
       XI_ULE1 =  SST_GO
       XI_ULE2 = -SST_GP
      ELSE
       XI_ULE1 = -SST_GO
       XI_ULE2 =  SST_GP
      ENDIF
C
C---- stuff BL system coefficients into main Jacobian matrix
C
      DO 40 JV=1, NSYS
        VM(1,JV,IV) = VS1(1,3)*D1_M(JV) + VS1(1,4)*U1_M(JV)
     &              + VS2(1,3)*D2_M(JV) + VS2(1,4)*U2_M(JV)
     &              + (VS1(1,5) + VS2(1,5) + VSX(1))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   40 CONTINUE
C
      VB(1,1,IV) = VS1(1,1)
      VB(1,2,IV) = VS1(1,2)
C
      VA(1,1,IV) = VS2(1,1)
      VA(1,2,IV) = VS2(1,2)
C
      IF(LALFA) THEN
       VDEL(1,2,IV) = VSR(1)*RE_CLMR + VSM(1)*MSQ_CLMR
      ELSE
       VDEL(1,2,IV) = 
     &       (VS1(1,4)*U1_A + VS1(1,3)*D1_A)
     &     + (VS2(1,4)*U2_A + VS2(1,3)*D2_A)
     &     + (VS1(1,5) + VS2(1,5) + VSX(1))
     &      *(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      VDEL(1,1,IV) = VSREZ(1)
     &   + (VS1(1,4)*DUE1 + VS1(1,3)*DDS1)
     &   + (VS2(1,4)*DUE2 + VS2(1,3)*DDS2)
     &   + (VS1(1,5) + VS2(1,5) + VSX(1))
     &    *(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      DO 50 JV=1, NSYS
        VM(2,JV,IV) = VS1(2,3)*D1_M(JV) + VS1(2,4)*U1_M(JV)
     &              + VS2(2,3)*D2_M(JV) + VS2(2,4)*U2_M(JV)
     &              + (VS1(2,5) + VS2(2,5) + VSX(2))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   50 CONTINUE
C
      VB(2,1,IV)  = VS1(2,1)
      VB(2,2,IV)  = VS1(2,2)
C
      VA(2,1,IV) = VS2(2,1)
      VA(2,2,IV) = VS2(2,2)
C
      IF(LALFA) THEN
       VDEL(2,2,IV) = VSR(2)*RE_CLMR + VSM(2)*MSQ_CLMR
      ELSE
       VDEL(2,2,IV) = 
     &       (VS1(2,4)*U1_A + VS1(2,3)*D1_A)
     &     + (VS2(2,4)*U2_A + VS2(2,3)*D2_A)
     &     + (VS1(2,5) + VS2(2,5) + VSX(2))
     &      *(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      VDEL(2,1,IV) = VSREZ(2)
     &   + (VS1(2,4)*DUE1 + VS1(2,3)*DDS1)
     &   + (VS2(2,4)*DUE2 + VS2(2,3)*DDS2)
     &   + (VS1(2,5) + VS2(2,5) + VSX(2))
     &    *(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      DO 60 JV=1, NSYS
        VM(3,JV,IV) = VS1(3,3)*D1_M(JV) + VS1(3,4)*U1_M(JV)
     &              + VS2(3,3)*D2_M(JV) + VS2(3,4)*U2_M(JV)
     &              + (VS1(3,5) + VS2(3,5) + VSX(3))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   60 CONTINUE
C
      VB(3,1,IV) = VS1(3,1)
      VB(3,2,IV) = VS1(3,2)
C
      VA(3,1,IV) = VS2(3,1)
      VA(3,2,IV) = VS2(3,2)
C
      IF(LALFA) THEN
       VDEL(3,2,IV) = VSR(3)*RE_CLMR + VSM(3)*MSQ_CLMR
      ELSE
       VDEL(3,2,IV) = 
     &       (VS1(3,4)*U1_A + VS1(3,3)*D1_A)
     &     + (VS2(3,4)*U2_A + VS2(3,3)*D2_A)
     &     + (VS1(3,5) + VS2(3,5) + VSX(3))
     &      *(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      VDEL(3,1,IV) = VSREZ(3)
     &   + (VS1(3,4)*DUE1 + VS1(3,3)*DDS1)
     &   + (VS2(3,4)*DUE2 + VS2(3,3)*DDS2)
     &   + (VS1(3,5) + VS2(3,5) + VSX(3))
     &    *(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      IF(IBL.EQ.IBLTE(IS)+1) THEN
C
C----- redefine coefficients for TTE, DTE, etc
       VZ(1,1)    = VS1(1,1)*CTE_CTE1
       VZ(1,2)    = VS1(1,1)*CTE_TTE1 + VS1(1,2)*TTE_TTE1
       VB(1,1,IV) = VS1(1,1)*CTE_CTE2
       VB(1,2,IV) = VS1(1,1)*CTE_TTE2 + VS1(1,2)*TTE_TTE2
C
       VZ(2,1)    = VS1(2,1)*CTE_CTE1
       VZ(2,2)    = VS1(2,1)*CTE_TTE1 + VS1(2,2)*TTE_TTE1
       VB(2,1,IV) = VS1(2,1)*CTE_CTE2
       VB(2,2,IV) = VS1(2,1)*CTE_TTE2 + VS1(2,2)*TTE_TTE2
C
       VZ(3,1)    = VS1(3,1)*CTE_CTE1
       VZ(3,2)    = VS1(3,1)*CTE_TTE1 + VS1(3,2)*TTE_TTE1
       VB(3,1,IV) = VS1(3,1)*CTE_CTE2
       VB(3,2,IV) = VS1(3,1)*CTE_TTE2 + VS1(3,2)*TTE_TTE2
C
      ENDIF
C
C---- turbulent intervals will follow if currently at transition interval
      IF(TRAN) THEN
        TURB = .TRUE.
C
C------ save transition location
        ITRAN(IS) = IBL
        TFORCE(IS) = TRFORC
        XSSITR(IS) = XT
C
C------ interpolate airfoil geometry to find transition x/c
C-      (for user output)
        IF(IS.EQ.1) THEN
         STR = SST - XT
        ELSE
         STR = SST + XT
        ENDIF
        CHX = XTE - XLE
        CHY = YTE - YLE
        CHSQ = CHX**2 + CHY**2
        XTR = SEVAL(STR,X,XP,S,N)
        YTR = SEVAL(STR,Y,YP,S,N)
        XOCTR(IS) = ((XTR-XLE)*CHX + (YTR-YLE)*CHY)/CHSQ
        YOCTR(IS) = ((YTR-YLE)*CHX - (XTR-XLE)*CHY)/CHSQ
      ENDIF
C
      TRAN = .FALSE.
C
      IF(IBL.EQ.IBLTE(IS)) THEN
C----- set "2" variables at TE to wake correlations for next station
C
       TURB = .TRUE.
       WAKE = .TRUE.
       CALL BLVAR(3)
       CALL BLMID(3)
      ENDIF
C
      DO 80 JS=1, 2
        DO 810 JBL=2, NBL(JS)
          JV = ISYS(JBL,JS)
          U1_M(JV) = U2_M(JV)
          D1_M(JV) = D2_M(JV)
  810   CONTINUE
   80 CONTINUE
C
      U1_A = U2_A
      D1_A = D2_A
C
      DUE1 = DUE2
      DDS1 = DDS2
C      
      IF(IBL .EQ. ITRAN(IS) .AND. X2 .GT. X1) THEN
       IF(IS.EQ.1) THEN
        TINDEX(IS) = FLOAT(IST-ITRAN(IS)+3) - (XT-X1)/(X2-X1)
       ELSE
        TINDEX(IS) = FLOAT(IST+ITRAN(IS)-2) + (XT-X1)/(X2-X1)
       ENDIF
      ENDIF
C
C---- set BL variables for next station
      DO 190 ICOM=1, NCOM
        COM1(ICOM) = COM2(ICOM)
  190 CONTINUE
C
C---- next streamwise station
 1000 CONTINUE
C
      IF(TFORCE(IS)) THEN
       WRITE(*,9100) IS,XOCTR(IS),ITRAN(IS)
 9100  FORMAT(1X,'Side',I2,' forced transition at x/c = ',F7.4,I5)
      ELSE
       WRITE(*,9200) IS,XOCTR(IS),ITRAN(IS)
 9200  FORMAT(1X,'Side',I2,'  free  transition at x/c = ',F7.4,I5)
      ENDIF
C
C---- next airfoil side
 2000 CONTINUE
C
      RETURN
      END


      SUBROUTINE IBLSYS
C---------------------------------------------
C     Sets the BL Newton system line number
C     corresponding to each BL station.
C---------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
C
      IV = 0
      DO 10 IS=1, 2
        DO 110 IBL=2, NBL(IS)
          IV = IV+1
          ISYS(IBL,IS) = IV
  110   CONTINUE
   10 CONTINUE
C
      NSYS = IV
      IF(NSYS.GT.2*IVX) STOP '*** IBLSYS: BL system array overflow. ***'
C
      RETURN
      END


      SUBROUTINE MRCHUE
C----------------------------------------------------
C     Marches the BLs and wake in direct mode using
C     the UEDG array. If separation is encountered,
C     a plausible value of Hk extrapolated from
C     upstream is prescribed instead.  Continuous
C     checking of transition onset is performed.
C----------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      LOGICAL DIRECT
      REAL MSQ
C
C---- shape parameters for separation criteria
      HLMAX = 3.8
      HTMAX = 2.5
C
      DO 2000 IS = 1, 2
C
      WRITE(*,*) '   side ', IS, ' ...'
C
      AMCRIT = ACRIT(IS)
C
C---- set forced transition arc length position
      CALL XIFSET(IS)
C
C---- initialize similarity station with Thwaites' formula
      IBL = 2
      XSI = XSSI(IBL,IS)
      UEI = UEDG(IBL,IS)
C      BULE = LOG(UEDG(IBL+1,IS)/UEI) / LOG(XSSI(IBL+1,IS)/XSI)
C      BULE = MAX( -.08 , BULE )
      BULE = 1.0
      UCON = UEI/XSI**BULE
      TSQ = 0.45/(UCON*(5.0*BULE+1.0)*REYBL) * XSI**(1.0-BULE)
      THI = SQRT(TSQ)
      DSI = 2.2*THI
      AMI = 0.0
C
C---- initialize Ctau for first turbulent station
      CTI = 0.03
C
      TRAN = .FALSE.
      TURB = .FALSE.
      ITRAN(IS) = IBLTE(IS)
C
C---- march downstream
      DO 1000 IBL=2, NBL(IS)
        IBM = IBL-1
C
        IW = IBL - IBLTE(IS)
C
        SIMI = IBL.EQ.2
        WAKE = IBL.GT.IBLTE(IS)
C
C------ prescribed quantities
        XSI = XSSI(IBL,IS)
        UEI = UEDG(IBL,IS)
C
        IF(WAKE) THEN
         IW = IBL - IBLTE(IS)
         DSWAKI = WGAP(IW)
        ELSE
         DSWAKI = 0.
        ENDIF
C
        DIRECT = .TRUE.
C
C------ Newton iteration loop for current station
        DO 100 ITBL=1, 25
C
C-------- assemble 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi
C         at the previous "1" station and the current "2" station
C         (the "1" station coefficients will be ignored)
C
C
          CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
          CALL BLKIN
C
C-------- check for transition and set appropriate flags and things
          IF((.NOT.SIMI) .AND. (.NOT.TURB)) THEN
           CALL TRCHEK
           AMI = AMPL2
C
           IF(TRAN) THEN
            ITRAN(IS) = IBL
            IF(CTI.LE.0.0) THEN
             CTI = 0.03
             S2 = CTI
            ENDIF
           ELSE
            ITRAN(IS) = IBL+2
           ENDIF
C
C
          ENDIF
C
          IF(IBL.EQ.IBLTE(IS)+1) THEN
           TTE = THET(IBLTE(1),1) + THET(IBLTE(2),2)
           DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
           CTE = ( CTAU(IBLTE(1),1)*THET(IBLTE(1),1)
     &           + CTAU(IBLTE(2),2)*THET(IBLTE(2),2) ) / TTE
           CALL TESYS(CTE,TTE,DTE)
          ELSE
           CALL BLSYS
          ENDIF
C

          IF(DIRECT) THEN
C
C--------- try direct mode (set dUe = 0 in currently empty 4th line)
           VS2(4,1) = 0.
           VS2(4,2) = 0.
           VS2(4,3) = 0.
           VS2(4,4) = 1.0
           VSREZ(4) = 0.
C
C--------- solve Newton system for current "2" station
           CALL GAUSS(4,4,VS2,VSREZ,1)
C
C--------- determine max changes and underrelax if necessary
           DMAX = MAX( ABS(VSREZ(2)/THI),
     &                 ABS(VSREZ(3)/DSI)  )
           IF(IBL.LT.ITRAN(IS)) DMAX = MAX(DMAX,ABS(VSREZ(1)/10.0))
           IF(IBL.GE.ITRAN(IS)) DMAX = MAX(DMAX,ABS(VSREZ(1)/CTI ))
C
           RLX = 1.0
           IF(DMAX.GT.0.3) RLX = 0.3/DMAX
C
C--------- see if direct mode is not applicable
           IF(IBL .NE. IBLTE(IS)+1) THEN
C
C---------- calculate resulting kinematic shape parameter Hk
            MSQ = UEI*UEI*HSTINV / (GM1BL*(1.0 - 0.5*UEI*UEI*HSTINV))
            HTEST = (DSI + RLX*VSREZ(3)) / (THI + RLX*VSREZ(2))
            CALL HKIN( HTEST, MSQ, HKTEST, DUMMY, DUMMY)
C
C---------- decide whether to do direct or inverse problem based on Hk
            IF(IBL.LT.ITRAN(IS)) HMAX = HLMAX
            IF(IBL.GE.ITRAN(IS)) HMAX = HTMAX
            DIRECT = HKTEST.LT.HMAX
           ENDIF
C
           IF(DIRECT) THEN
C---------- update as usual
ccc            IF(IBL.LT.ITRAN(IS)) AMI = AMI + RLX*VSREZ(1)
            IF(IBL.GE.ITRAN(IS)) CTI = CTI + RLX*VSREZ(1)
            THI = THI + RLX*VSREZ(2)
            DSI = DSI + RLX*VSREZ(3)
           ELSE
C---------- set prescribed Hk for inverse calculation at the current station
            IF(IBL.LT.ITRAN(IS)) THEN
C----------- laminar case: relatively slow increase in Hk downstream
             HTARG = HK1 + 0.03*(X2-X1)/T1
            ELSE IF(IBL.EQ.ITRAN(IS)) THEN
C----------- transition interval: weighted laminar and turbulent case
             HTARG = HK1 + (0.03*(XT-X1) - 0.15*(X2-XT))/T1
            ELSE IF(WAKE) THEN
C----------- turbulent wake case:
C-           asymptotic wake behavior with approximate Backward Euler
             CONST = 0.03*(X2-X1)/T1
             HK2 = HK1
             HK2 = HK2 - (HK2 +     CONST*(HK2-1.0)**3 - HK1)
     &                  /(1.0 + 3.0*CONST*(HK2-1.0)**2)
             HK2 = HK2 - (HK2 +     CONST*(HK2-1.0)**3 - HK1)
     &                  /(1.0 + 3.0*CONST*(HK2-1.0)**2)
             HK2 = HK2 - (HK2 +     CONST*(HK2-1.0)**3 - HK1)
     &                  /(1.0 + 3.0*CONST*(HK2-1.0)**2)
             HTARG = HK2
            ELSE
C----------- turbulent case: relatively fast decrease in Hk downstream
             HTARG = HK1 - 0.15*(X2-X1)/T1
            ENDIF
C
C---------- limit specified Hk to something reasonable
            IF(WAKE) THEN
             HTARG = MAX( HTARG , 1.01 )
            ELSE
             HTARG = MAX( HTARG , HMAX )
            ENDIF
C
            WRITE(*,1300) IBL, HTARG
 1300       FORMAT(' MRCHUE: Inverse mode at', I4, '     Hk =', F8.3)
C
C---------- try again with prescribed Hk
            GO TO 100
C
           ENDIF
C
          ELSE
C
C-------- inverse mode (force Hk to prescribed value HTARG)
           VS2(4,1) = 0.
           VS2(4,2) = HK2_T2
           VS2(4,3) = HK2_D2
           VS2(4,4) = HK2_U2
           VSREZ(4) = HTARG - HK2
C
           CALL GAUSS(4,4,VS2,VSREZ,1)
C
C--------- added Ue clamp   MD  3 Apr 03
           DMAX = MAX( ABS(VSREZ(2)/THI),
     &                 ABS(VSREZ(3)/DSI),
     &                 ABS(VSREZ(4)/UEI)  )
           IF(IBL.GE.ITRAN(IS)) DMAX = MAX( DMAX , ABS(VSREZ(1)/CTI))
C
           RLX = 1.0
           IF(DMAX.GT.0.3) RLX = 0.3/DMAX
C
C--------- update variables
ccc           IF(IBL.LT.ITRAN(IS)) AMI = AMI + RLX*VSREZ(1)
           IF(IBL.GE.ITRAN(IS)) CTI = CTI + RLX*VSREZ(1)
           THI = THI + RLX*VSREZ(2)
           DSI = DSI + RLX*VSREZ(3)
           UEI = UEI + RLX*VSREZ(4)
C
          ENDIF
C
C-------- eliminate absurd transients
          IF(IBL.GE.ITRAN(IS)) THEN
           CTI = MIN(CTI , 0.30 )
           CTI = MAX(CTI , 0.0000001 )
          ENDIF
C
          IF(IBL.LE.IBLTE(IS)) THEN
            HKLIM = 1.02
          ELSE
            HKLIM = 1.00005
          ENDIF
          MSQ = UEI*UEI*HSTINV / (GM1BL*(1.0 - 0.5*UEI*UEI*HSTINV))
          DSW = DSI - DSWAKI
          CALL DSLIM(DSW,THI,UEI,MSQ,HKLIM)
          DSI = DSW + DSWAKI
C
          IF(DMAX.LE.1.0E-5) GO TO 110
C
  100   CONTINUE
        WRITE(*,1350) IBL, IS, DMAX 
 1350   FORMAT(' MRCHUE: Convergence failed at',I4,'  side',I2,
     &         '    Res =', E12.4)
C
C------ the current unconverged solution might still be reasonable...
CCC        IF(DMAX .LE. 0.1) GO TO 110
        IF(DMAX .LE. 0.1) GO TO 109
C
C------- the current solution is garbage --> extrapolate values instead
         IF(IBL.GT.3) THEN 
          IF(IBL.LE.IBLTE(IS)) THEN
           THI = THET(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5
           DSI = DSTR(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5
          ELSE IF(IBL.EQ.IBLTE(IS)+1) THEN
           CTI = CTE
           THI = TTE
           DSI = DTE
          ELSE
           THI = THET(IBM,IS)
           RATLEN = (XSSI(IBL,IS)-XSSI(IBM,IS)) / (10.0*DSTR(IBM,IS))
           DSI = (DSTR(IBM,IS) + THI*RATLEN) / (1.0 + RATLEN)
          ENDIF
          IF(IBL.EQ.ITRAN(IS)) CTI = 0.05
          IF(IBL.GT.ITRAN(IS)) CTI = CTAU(IBM,IS)
C
          UEI = UEDG(IBL,IS)
          IF(IBL.GT.2 .AND. IBL.LT.NBL(IS))
     &     UEI = 0.5*(UEDG(IBL-1,IS) + UEDG(IBL+1,IS))
         ENDIF
C
 109     CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
         CALL BLKIN
C
C------- check for transition and set appropriate flags and things
         IF((.NOT.SIMI) .AND. (.NOT.TURB)) THEN
          CALL TRCHEK
          AMI = AMPL2
          IF(     TRAN) ITRAN(IS) = IBL
          IF(.NOT.TRAN) ITRAN(IS) = IBL+2
         ENDIF
C
C------- set all other extrapolated values for current station
         IF(IBL.LT.ITRAN(IS)) CALL BLVAR(1)
         IF(IBL.GE.ITRAN(IS)) CALL BLVAR(2)
         IF(WAKE) CALL BLVAR(3)
C
         IF(IBL.LT.ITRAN(IS)) CALL BLMID(1)
         IF(IBL.GE.ITRAN(IS)) CALL BLMID(2)
         IF(WAKE) CALL BLMID(3)
C
C------ pick up here after the Newton iterations
  110   CONTINUE
C
C------ store primary variables
        IF(IBL.LT.ITRAN(IS)) CTAU(IBL,IS) = AMI
        IF(IBL.GE.ITRAN(IS)) CTAU(IBL,IS) = CTI
        THET(IBL,IS) = THI
        DSTR(IBL,IS) = DSI
        UEDG(IBL,IS) = UEI
        MASS(IBL,IS) = DSI*UEI
        TAU(IBL,IS)  = 0.5*R2*U2*U2*CF2
        DIS(IBL,IS)  =     R2*U2*U2*U2*DI2*HS2*0.5
        CTQ(IBL,IS)  = CQ2
        DELT(IBL,IS) = DE2
        TSTR(IBL,IS) = HS2*T2
C
C------ set "1" variables to "2" variables for next streamwise station
        CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
        CALL BLKIN
        DO 310 ICOM=1, NCOM
          COM1(ICOM) = COM2(ICOM)
  310   CONTINUE
C
C------ turbulent intervals will follow transition interval or TE
        IF(TRAN .OR. IBL.EQ.IBLTE(IS)) THEN
         TURB = .TRUE.
C
C------- save transition location
         TFORCE(IS) = TRFORC
         XSSITR(IS) = XT
        ENDIF
C
        TRAN = .FALSE.
C
        IF(IBL.EQ.IBLTE(IS)) THEN
         THI = THET(IBLTE(1),1) + THET(IBLTE(2),2)
         DSI = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
        ENDIF
C
 1000 CONTINUE
 2000 CONTINUE
C
      RETURN
      END
  
 
      SUBROUTINE MRCHDU
C----------------------------------------------------
C     Marches the BLs and wake in mixed mode using
C     the current Ue and Hk.  The calculated Ue
C     and Hk lie along a line quasi-normal to the
C     natural Ue-Hk characteristic line of the
C     current BL so that the Goldstein or Levy-Lees
C     singularity is never encountered.  Continuous
C     checking of transition onset is performed.
C----------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      REAL VTMP(4,5), VZTMP(4)
      REAL MSQ
ccc   REAL MDI
C
      DATA DEPS / 5.0E-6 /
C
C---- constant controlling how far Hk is allowed to deviate
C-    from the specified value.
      SENSWT = 1000.0
C
      DO 2000 IS = 1, 2
C
      AMCRIT = ACRIT(IS)
C
C---- set forced transition arc length position
      CALL XIFSET(IS)
C
C---- set leading edge pressure gradient parameter  x/u du/dx
      IBL = 2
      XSI = XSSI(IBL,IS)
      UEI = UEDG(IBL,IS)
CCC      BULE = LOG(UEDG(IBL+1,IS)/UEI) / LOG(XSSI(IBL+1,IS)/XSI)
CCC      BULE = MAX( -.08 , BULE )
      BULE = 1.0
C
C---- old transition station
      ITROLD = ITRAN(IS)
C
      TRAN = .FALSE.
      TURB = .FALSE.
      ITRAN(IS) = IBLTE(IS)
C
C---- march downstream
      DO 1000 IBL=2, NBL(IS)
        IBM = IBL-1
C
        SIMI = IBL.EQ.2
        WAKE = IBL.GT.IBLTE(IS)
C
C------ initialize current station to existing variables
        XSI = XSSI(IBL,IS)
        UEI = UEDG(IBL,IS)
        THI = THET(IBL,IS)
        DSI = DSTR(IBL,IS)

CCC        MDI = MASS(IBL,IS)
C
C------ fixed BUG   MD 7 June 99
        IF(IBL.LT.ITROLD) THEN
         AMI = CTAU(IBL,IS)
         CTI = 0.03
        ELSE
         CTI = CTAU(IBL,IS)
         IF(CTI.LE.0.0) CTI = 0.03
        ENDIF
C
CCC        DSI = MDI/UEI
C
        IF(WAKE) THEN
         IW = IBL - IBLTE(IS)
         DSWAKI = WGAP(IW)
        ELSE
         DSWAKI = 0.
        ENDIF
C
        IF(IBL.LE.IBLTE(IS)) DSI = MAX(DSI-DSWAKI,1.02000*THI) + DSWAKI
        IF(IBL.GT.IBLTE(IS)) DSI = MAX(DSI-DSWAKI,1.00005*THI) + DSWAKI
C
C------ Newton iteration loop for current station
        DO 100 ITBL=1, 25
C
C-------- assemble 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi
C         at the previous "1" station and the current "2" station
C         (the "1" station coefficients will be ignored)
C
          CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
          CALL BLKIN
C
C-------- check for transition and set appropriate flags and things
          IF((.NOT.SIMI) .AND. (.NOT.TURB)) THEN
           CALL TRCHEK
           AMI = AMPL2
           IF(     TRAN) ITRAN(IS) = IBL
           IF(.NOT.TRAN) ITRAN(IS) = IBL+2
          ENDIF
C
          IF(IBL.EQ.IBLTE(IS)+1) THEN
           TTE = THET(IBLTE(1),1) + THET(IBLTE(2),2)
           DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
           CTE = ( CTAU(IBLTE(1),1)*THET(IBLTE(1),1)
     &           + CTAU(IBLTE(2),2)*THET(IBLTE(2),2) ) / TTE
           CALL TESYS(CTE,TTE,DTE)
          ELSE
           CALL BLSYS
          ENDIF
C
C-------- set stuff at first iteration...
          IF(ITBL.EQ.1) THEN
C
C--------- set "baseline" Ue and Hk for forming  Ue(Hk)  relation
           UEREF = U2
           HKREF = HK2
C
C--------- if current point IBL was turbulent and is now laminar, then...
           IF(IBL.LT.ITRAN(IS) .AND. IBL.GE.ITROLD ) THEN
C---------- extrapolate baseline Hk
            UEM = UEDG(IBL-1,IS)
            DSM = DSTR(IBL-1,IS)
            THM = THET(IBL-1,IS)
            MSQ = UEM*UEM*HSTINV / (GM1BL*(1.0 - 0.5*UEM*UEM*HSTINV))
            CALL HKIN( DSM/THM, MSQ, HKREF, DUMMY, DUMMY )
           ENDIF
C
C--------- if current point IBL was laminar, then...
           IF(IBL.LT.ITROLD) THEN
C---------- reinitialize or extrapolate Ctau if it's now turbulent
            IF(TRAN) CTAU(IBL,IS) = 0.03
            IF(TURB) CTAU(IBL,IS) = CTAU(IBL-1,IS)
            IF(TRAN .OR. TURB) THEN
             CTI = CTAU(IBL,IS)
             S2 = CTI
            ENDIF
           ENDIF
C
          ENDIF
C
C
          IF(SIMI .OR. IBL.EQ.IBLTE(IS)+1) THEN
C
C--------- for similarity station or first wake point, prescribe Ue
           VS2(4,1) = 0.
           VS2(4,2) = 0.
           VS2(4,3) = 0.
           VS2(4,4) = U2_UEI
           VSREZ(4) = UEREF - U2
C
          ELSE
C
C********* calculate Ue-Hk characteristic slope
C
           DO 20 K=1, 4
             VZTMP(K) = VSREZ(K)
             DO 201 L=1, 5
               VTMP(K,L) = VS2(K,L)
  201        CONTINUE
   20      CONTINUE
C
C--------- set unit dHk
           VTMP(4,1) = 0.
           VTMP(4,2) = HK2_T2
           VTMP(4,3) = HK2_D2
           VTMP(4,4) = HK2_U2*U2_UEI
           VZTMP(4)  = 1.0
C
C--------- calculate dUe response
           CALL GAUSS(4,4,VTMP,VZTMP,1)
C
C--------- set  SENSWT * (normalized dUe/dHk)
           SENNEW = SENSWT * VZTMP(4) * HKREF/UEREF
           IF(ITBL.LE.5) THEN
            SENS = SENNEW
           ELSE IF(ITBL.LE.15) THEN
            SENS = 0.5*(SENS + SENNEW)
           ENDIF
C
C--------- set prescribed Ue-Hk combination
           VS2(4,1) = 0.
           VS2(4,2) =  HK2_T2 * HKREF
           VS2(4,3) =  HK2_D2 * HKREF
           VS2(4,4) =( HK2_U2 * HKREF  +  SENS/UEREF )*U2_UEI
           VSREZ(4) = -(HKREF**2)*(HK2 / HKREF - 1.0)
     &                     - SENS*(U2  / UEREF - 1.0)
C
          ENDIF
C
C-------- solve Newton system for current "2" station
          CALL GAUSS(4,4,VS2,VSREZ,1)
C
C-------- determine max changes and underrelax if necessary
C-------- (added Ue clamp   MD  3 Apr 03)
          DMAX = MAX( ABS(VSREZ(2)/THI),
     &                ABS(VSREZ(3)/DSI),
     &                ABS(VSREZ(4)/UEI)  )
          IF(IBL.GE.ITRAN(IS)) DMAX = MAX(DMAX,ABS(VSREZ(1)/(10.0*CTI)))
C
          RLX = 1.0
          IF(DMAX.GT.0.3) RLX = 0.3/DMAX
C
C-------- update as usual
          IF(IBL.LT.ITRAN(IS)) AMI = AMI + RLX*VSREZ(1)
          IF(IBL.GE.ITRAN(IS)) CTI = CTI + RLX*VSREZ(1)
          THI = THI + RLX*VSREZ(2)
          DSI = DSI + RLX*VSREZ(3)
          UEI = UEI + RLX*VSREZ(4)
C
C-------- eliminate absurd transients
          IF(IBL.GE.ITRAN(IS)) THEN
           CTI = MIN(CTI , 0.30 )
           CTI = MAX(CTI , 0.0000001 )
          ENDIF
C
          IF(IBL.LE.IBLTE(IS)) THEN
            HKLIM = 1.02
          ELSE
            HKLIM = 1.00005
          ENDIF
          MSQ = UEI*UEI*HSTINV / (GM1BL*(1.0 - 0.5*UEI*UEI*HSTINV))
          DSW = DSI - DSWAKI
          CALL DSLIM(DSW,THI,UEI,MSQ,HKLIM)
          DSI = DSW + DSWAKI
C
          IF(DMAX.LE.DEPS) GO TO 110
C
  100   CONTINUE
C
        WRITE(*,1350) IBL, IS, DMAX 
 1350   FORMAT(' MRCHDU: Convergence failed at',I4,'  side',I2,
     &         '    Res =', E12.4)
C
C------ the current unconverged solution might still be reasonable...
CCC        IF(DMAX .LE. 0.1) GO TO 110
        IF(DMAX .LE. 0.1) GO TO 109
C
C------- the current solution is garbage --> extrapolate values instead
         IF(IBL.GT.3) THEN
          IF(IBL.LE.IBLTE(IS)) THEN
           THI = THET(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5
           DSI = DSTR(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5
           UEI = UEDG(IBM,IS)
          ELSE IF(IBL.EQ.IBLTE(IS)+1) THEN
           CTI = CTE
           THI = TTE
           DSI = DTE
           UEI = UEDG(IBM,IS)
          ELSE
           THI = THET(IBM,IS)
           RATLEN = (XSSI(IBL,IS)-XSSI(IBM,IS)) / (10.0*DSTR(IBM,IS))
           DSI = (DSTR(IBM,IS) + THI*RATLEN) / (1.0 + RATLEN)
           UEI = UEDG(IBM,IS)
          ENDIF
          IF(IBL.EQ.ITRAN(IS)) CTI = 0.05
          IF(IBL.GT.ITRAN(IS)) CTI = CTAU(IBM,IS)
         ENDIF
C
 109     CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
         CALL BLKIN
C
C------- check for transition and set appropriate flags and things
         IF((.NOT.SIMI) .AND. (.NOT.TURB)) THEN
          CALL TRCHEK
          AMI = AMPL2
          IF(     TRAN) ITRAN(IS) = IBL
          IF(.NOT.TRAN) ITRAN(IS) = IBL+2
         ENDIF
C
C------- set all other extrapolated values for current station
         IF(IBL.LT.ITRAN(IS)) CALL BLVAR(1)
         IF(IBL.GE.ITRAN(IS)) CALL BLVAR(2)
         IF(WAKE) CALL BLVAR(3)
C
         IF(IBL.LT.ITRAN(IS)) CALL BLMID(1)
         IF(IBL.GE.ITRAN(IS)) CALL BLMID(2)
         IF(WAKE) CALL BLMID(3)
C
C------ pick up here after the Newton iterations
  110   CONTINUE
C
        SENS = SENNEW
C
C------ store primary variables
        IF(IBL.LT.ITRAN(IS)) CTAU(IBL,IS) = AMI
        IF(IBL.GE.ITRAN(IS)) CTAU(IBL,IS) = CTI
        THET(IBL,IS) = THI
        DSTR(IBL,IS) = DSI
        UEDG(IBL,IS) = UEI
        MASS(IBL,IS) = DSI*UEI
        TAU(IBL,IS)  = 0.5*R2*U2*U2*CF2
        DIS(IBL,IS)  =     R2*U2*U2*U2*DI2*HS2*0.5
        CTQ(IBL,IS)  = CQ2
        DELT(IBL,IS) = DE2
        TSTR(IBL,IS) = HS2*T2
C
C------ set "1" variables to "2" variables for next streamwise station
        CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
        CALL BLKIN
        DO 310 ICOM=1, NCOM
          COM1(ICOM) = COM2(ICOM)
  310   CONTINUE
C
C
C------ turbulent intervals will follow transition interval or TE
        IF(TRAN .OR. IBL.EQ.IBLTE(IS)) THEN
         TURB = .TRUE.
C
C------- save transition location
         TFORCE(IS) = TRFORC
         XSSITR(IS) = XT
        ENDIF
C
        TRAN = .FALSE.
C
 1000 CONTINUE
C
 2000 CONTINUE
C
      RETURN
      END
  
 
      SUBROUTINE XIFSET(IS)
C-----------------------------------------------------
C     Sets forced-transition BL coordinate locations.
C-----------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
C
      IF(XSTRIP(IS).GE.1.0) THEN
       XIFORC = XSSI(IBLTE(IS),IS)
       RETURN
      ENDIF
C
      CHX = XTE - XLE
      CHY = YTE - YLE
      CHSQ = CHX**2 + CHY**2
C
C---- calculate chord-based x/c, y/c
      DO 10 I=1, N
        W1(I) = ((X(I)-XLE)*CHX + (Y(I)-YLE)*CHY) / CHSQ
        W2(I) = ((Y(I)-YLE)*CHX - (X(I)-XLE)*CHY) / CHSQ
 10   CONTINUE
C
      CALL SPLIND(W1,W3,S,N,-999.0,-999.0)
      CALL SPLIND(W2,W4,S,N,-999.0,-999.0)
C
      IF(IS.EQ.1) THEN
C
C----- set approximate arc length of forced transition point for SINVRT
       STR = SLE + (S(1)-SLE)*XSTRIP(IS)
C
C----- calculate actual arc length
       CALL SINVRT(STR,XSTRIP(IS),W1,W3,S,N)
C
C----- set BL coordinate value
       XIFORC = MIN( (SST - STR) , XSSI(IBLTE(IS),IS) )
C
      ELSE
C----- same for bottom side
C
       STR = SLE + (S(N)-SLE)*XSTRIP(IS)
       CALL SINVRT(STR,XSTRIP(IS),W1,W3,S,N)
       XIFORC = MIN( (STR - SST) , XSSI(IBLTE(IS),IS) )
C
      ENDIF
C
      IF(XIFORC .LT. 0.0) THEN
       WRITE(*,1000) IS
 1000  FORMAT(/' ***  Stagnation point is past trip on side',I2,'  ***')
       XIFORC = XSSI(IBLTE(IS),IS)
      ENDIF
C
      RETURN
      END




      SUBROUTINE UPDATE
C------------------------------------------------------------------
C      Adds on Newton deltas to boundary layer variables.
C      Checks for excessive changes and underrelaxes if necessary.
C      Calculates max and rms changes.
C      Also calculates the change in the global variable "AC".
C        If LALFA=.TRUE. , "AC" is CL
C        If LALFA=.FALSE., "AC" is alpha
C------------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL UNEW(IVX,2), U_AC(IVX,2)
      REAL QNEW(IQX),   Q_AC(IQX)
      EQUIVALENCE (VA(1,1,1), UNEW(1,1)) ,
     &            (VB(1,1,1), QNEW(1)  )
      EQUIVALENCE (VA(1,1,IVX), U_AC(1,1)) ,
     &            (VB(1,1,IVX), Q_AC(1)  )
      REAL MSQ
C
C---- max allowable alpha changes per iteration
      DALMAX =  0.5*DTOR
      DALMIN = -0.5*DTOR
C
C---- max allowable CL change per iteration
      DCLMAX =  0.5
      DCLMIN = -0.5
      IF(MATYP.NE.1) DCLMIN = MAX(-0.5 , -0.9*CL)
C
      HSTINV = GAMM1*(MINF/QINF)**2 / (1.0 + 0.5*GAMM1*MINF**2)
C
C---- calculate new Ue distribution assuming no under-relaxation
C-    also set the sensitivity of Ue wrt to alpha or Re
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
C
          DUI    = 0.
          DUI_AC = 0.
          DO 100 JS=1, 2
            DO 1000 JBL=2, NBL(JS)
              J  = IPAN(JBL,JS)
              JV = ISYS(JBL,JS)
              UE_M = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
              DUI    = DUI    + UE_M*(MASS(JBL,JS)+VDEL(3,1,JV))
              DUI_AC = DUI_AC + UE_M*(            -VDEL(3,2,JV))
 1000       CONTINUE
  100     CONTINUE
C
C-------- UINV depends on "AC" only if "AC" is alpha
          IF(LALFA) THEN
           UINV_AC = 0.
          ELSE
           UINV_AC = UINV_A(IBL,IS)
          ENDIF
C
          UNEW(IBL,IS) = UINV(IBL,IS) + DUI
          U_AC(IBL,IS) = UINV_AC      + DUI_AC
C
   10   CONTINUE
    1 CONTINUE
C
C---- set new Qtan from new Ue with appropriate sign change
      DO 2 IS=1, 2
        DO 20 IBL=2, IBLTE(IS)
          I = IPAN(IBL,IS)
          QNEW(I) = VTI(IBL,IS)*UNEW(IBL,IS)
          Q_AC(I) = VTI(IBL,IS)*U_AC(IBL,IS)
   20   CONTINUE
    2 CONTINUE
C
C---- calculate new CL from this new Qtan
      SA = SIN(ALFA)
      CA = COS(ALFA)
C
      BETA = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
C
      BFAC     = 0.5*MINF**2 / (1.0 + BETA)
      BFAC_MSQ = 0.5         / (1.0 + BETA)
     &         - BFAC        / (1.0 + BETA) * BETA_MSQ
C
      CLNEW = 0.
      CL_A  = 0.
      CL_MS = 0.
      CL_AC = 0.
C
      I = 1
      CGINC = 1.0 - (QNEW(I)/QINF)**2
      CPG1  = CGINC / (BETA + BFAC*CGINC)
      CPG1_MS = -CPG1/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
      CPI_Q = -2.0*QNEW(I)/QINF**2
      CPC_CPI = (1.0 - BFAC*CPG1)/ (BETA + BFAC*CGINC)
      CPG1_AC = CPC_CPI*CPI_Q*Q_AC(I)
C
      DO 3 I=1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
C
        CGINC = 1.0 - (QNEW(IP)/QINF)**2
        CPG2  = CGINC / (BETA + BFAC*CGINC)
        CPG2_MS = -CPG2/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
        CPI_Q = -2.0*QNEW(IP)/QINF**2
        CPC_CPI = (1.0 - BFAC*CPG2)/ (BETA + BFAC*CGINC)
        CPG2_AC = CPC_CPI*CPI_Q*Q_AC(IP)
C
        DX   =  (X(IP) - X(I))*CA + (Y(IP) - Y(I))*SA
        DX_A = -(X(IP) - X(I))*SA + (Y(IP) - Y(I))*CA
C
        AG    = 0.5*(CPG2    + CPG1   )
        AG_MS = 0.5*(CPG2_MS + CPG1_MS)
        AG_AC = 0.5*(CPG2_AC + CPG1_AC)
C
        CLNEW = CLNEW + DX  *AG
        CL_A  = CL_A  + DX_A*AG
        CL_MS = CL_MS + DX  *AG_MS
        CL_AC = CL_AC + DX  *AG_AC
C
        CPG1    = CPG2
        CPG1_MS = CPG2_MS
        CPG1_AC = CPG2_AC
    3 CONTINUE
C
C---- initialize under-relaxation factor
      RLX = 1.0
C
      IF(LALFA) THEN
C===== alpha is prescribed: AC is CL
C
C----- set change in Re to account for CL changing, since Re = Re(CL)
       DAC = (CLNEW - CL) / (1.0 - CL_AC - CL_MS*2.0*MINF*MINF_CL)
C
C----- set under-relaxation factor if Re change is too large
       IF(RLX*DAC .GT. DCLMAX) RLX = DCLMAX/DAC
       IF(RLX*DAC .LT. DCLMIN) RLX = DCLMIN/DAC
C
      ELSE
C===== CL is prescribed: AC is alpha
C
C----- set change in alpha to drive CL to prescribed value
       DAC = (CLNEW - CLSPEC) / (0.0 - CL_AC - CL_A)
C
C----- set under-relaxation factor if alpha change is too large
       IF(RLX*DAC .GT. DALMAX) RLX = DALMAX/DAC
       IF(RLX*DAC .LT. DALMIN) RLX = DALMIN/DAC
C
      ENDIF
C
      RMSBL = 0.
      RMXBL = 0.
C
      DHI = 1.5
      DLO = -.5
C
C---- calculate changes in BL variables and under-relaxation if needed
      DO 4 IS=1, 2
        DO 40 IBL=2, NBL(IS)
          IV = ISYS(IBL,IS)
C


C-------- set changes without underrelaxation
          DCTAU = VDEL(1,1,IV) - DAC*VDEL(1,2,IV)
          DTHET = VDEL(2,1,IV) - DAC*VDEL(2,2,IV)
          DMASS = VDEL(3,1,IV) - DAC*VDEL(3,2,IV)
          DUEDG = UNEW(IBL,IS) + DAC*U_AC(IBL,IS)  -  UEDG(IBL,IS)
          DDSTR = (DMASS - DSTR(IBL,IS)*DUEDG)/UEDG(IBL,IS)
C
C-------- normalize changes
          IF(IBL.LT.ITRAN(IS)) DN1 = DCTAU / 10.0
          IF(IBL.GE.ITRAN(IS)) DN1 = DCTAU / CTAU(IBL,IS)
          DN2 = DTHET / THET(IBL,IS)
          DN3 = DDSTR / DSTR(IBL,IS)
          DN4 = ABS(DUEDG)/0.25
C
C-------- accumulate for rms change
          RMSBL = RMSBL + DN1**2 + DN2**2 + DN3**2 + DN4**2
C          
C-------- see if Ctau needs underrelaxation
          RDN1 = RLX*DN1
          IF(ABS(DN1) .GT. ABS(RMXBL)) THEN
           RMXBL = DN1
           IF(IBL.LT.ITRAN(IS)) VMXBL = 'n'
           IF(IBL.GE.ITRAN(IS)) VMXBL = 'C'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN1 .GT. DHI) RLX = DHI/DN1
          IF(RDN1 .LT. DLO) RLX = DLO/DN1
C
C-------- see if Theta needs underrelaxation
          RDN2 = RLX*DN2
          IF(ABS(DN2) .GT. ABS(RMXBL)) THEN
           RMXBL = DN2
           VMXBL = 'T'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN2 .GT. DHI) RLX = DHI/DN2
          IF(RDN2 .LT. DLO) RLX = DLO/DN2
C
C-------- see if Dstar needs underrelaxation
          RDN3 = RLX*DN3
          IF(ABS(DN3) .GT. ABS(RMXBL)) THEN
           RMXBL = DN3
           VMXBL = 'D'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN3 .GT. DHI) RLX = DHI/DN3
          IF(RDN3 .LT. DLO) RLX = DLO/DN3
C
C-------- see if Ue needs underrelaxation
          RDN4 = RLX*DN4
          IF(ABS(DN4) .GT. ABS(RMXBL)) THEN
           RMXBL = DUEDG
           VMXBL = 'U'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN4 .GT. DHI) RLX = DHI/DN4
          IF(RDN4 .LT. DLO) RLX = DLO/DN4
C
   40   CONTINUE
    4 CONTINUE
C
C---- set true rms change
      RMSBL = SQRT( RMSBL / (4.0*FLOAT( NBL(1)+NBL(2) )) )
C
C
      IF(LALFA) THEN
C----- set underrelaxed change in Reynolds number from change in lift
       CL = CL + RLX*DAC
      ELSE
C----- set underrelaxed change in alpha
       ALFA = ALFA + RLX*DAC
       ADEG = ALFA/DTOR
      ENDIF
C
C---- update BL variables with underrelaxed changes
      DO 5 IS=1, 2
        DO 50 IBL=2, NBL(IS)
          IV = ISYS(IBL,IS)
C
          DCTAU = VDEL(1,1,IV) - DAC*VDEL(1,2,IV)
          DTHET = VDEL(2,1,IV) - DAC*VDEL(2,2,IV)
          DMASS = VDEL(3,1,IV) - DAC*VDEL(3,2,IV)
          DUEDG = UNEW(IBL,IS) + DAC*U_AC(IBL,IS)  -  UEDG(IBL,IS)
          DDSTR = (DMASS - DSTR(IBL,IS)*DUEDG)/UEDG(IBL,IS)
C
          CTAU(IBL,IS) = CTAU(IBL,IS) + RLX*DCTAU
          THET(IBL,IS) = THET(IBL,IS) + RLX*DTHET
          DSTR(IBL,IS) = DSTR(IBL,IS) + RLX*DDSTR
          UEDG(IBL,IS) = UEDG(IBL,IS) + RLX*DUEDG
C
          IF(IBL.GT.IBLTE(IS)) THEN
           IW = IBL - IBLTE(IS)
           DSWAKI = WGAP(IW)
          ELSE
           DSWAKI = 0.
          ENDIF
C
C-------- eliminate absurd transients
          IF(IBL.GE.ITRAN(IS))
     &      CTAU(IBL,IS) = MIN( CTAU(IBL,IS) , 0.25 )
C
          IF(IBL.LE.IBLTE(IS)) THEN
            HKLIM = 1.02
          ELSE
            HKLIM = 1.00005
          ENDIF
          MSQ = UEDG(IBL,IS)**2*HSTINV
     &        / (GAMM1*(1.0 - 0.5*UEDG(IBL,IS)**2*HSTINV))
          DSW = DSTR(IBL,IS) - DSWAKI
          CALL DSLIM(DSW,THET(IBL,IS),UEDG(IBL,IS),MSQ,HKLIM)
          DSTR(IBL,IS) = DSW + DSWAKI
C
C-------- set new mass defect (nonlinear update)
          MASS(IBL,IS) = DSTR(IBL,IS) * UEDG(IBL,IS)
C
   50   CONTINUE
C
C------ make sure there are no "islands" of negative Ue
        DO IBL = 3, IBLTE(IS)
          IF(UEDG(IBL-1,IS) .GT. 0.0 .AND.
     &       UEDG(IBL  ,IS) .LE. 0.0       ) THEN
           UEDG(IBL,IS) = UEDG(IBL-1,IS)
           MASS(IBL,IS) = DSTR(IBL,IS) * UEDG(IBL,IS)
          ENDIF
        ENDDO
    5 CONTINUE
C
C
C---- equate upper wake arrays to lower wake arrays
      DO 6 KBL=1, NBL(2)-IBLTE(2)
        CTAU(IBLTE(1)+KBL,1) = CTAU(IBLTE(2)+KBL,2)
        THET(IBLTE(1)+KBL,1) = THET(IBLTE(2)+KBL,2)
        DSTR(IBLTE(1)+KBL,1) = DSTR(IBLTE(2)+KBL,2)
        UEDG(IBLTE(1)+KBL,1) = UEDG(IBLTE(2)+KBL,2)
         TAU(IBLTE(1)+KBL,1) =  TAU(IBLTE(2)+KBL,2)
         DIS(IBLTE(1)+KBL,1) =  DIS(IBLTE(2)+KBL,2)
         CTQ(IBLTE(1)+KBL,1) =  CTQ(IBLTE(2)+KBL,2)
        DELT(IBLTE(1)+KBL,1) = DELT(IBLTE(2)+KBL,2)
        TSTR(IBLTE(1)+KBL,1) = TSTR(IBLTE(2)+KBL,2)
    6 CONTINUE
C
      RETURN
      END



      SUBROUTINE DSLIM(DSTR,THET,UEDG,MSQ,HKLIM)
      IMPLICIT REAL (A-H,M,O-Z)
C
      H = DSTR/THET
      CALL HKIN(H,MSQ,HK,HK_H,HK_M)
C
      DH = MAX( 0.0 , HKLIM-HK ) / HK_H
      DSTR = DSTR + DH*THET
C
      RETURN
      END



      SUBROUTINE BLPINI
      INCLUDE 'BLPAR.INC'
C
      SCCON = 5.6
      GACON = 6.70
      GBCON = 0.75
      GCCON = 18.0
      DLCON =  0.9
C
      CTRCON = 1.8
      CTRCEX = 3.3
C
      DUXCON = 1.0
C
      CTCON = 0.5/(GACON**2 * GBCON)
C
      CFFAC = 1.0
C
      RETURN
      END


C***********************************************************************
C    Module:  xblsys.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************


      SUBROUTINE TRCHEK
C
C---- 1st-order amplification equation
cc      CALL TRCHEK1
C
C---- 2nd-order amplification equation
      CALL TRCHEK2
C
      RETURN
      END



      SUBROUTINE AXSET( HK1,    T1,    RT1,    A1,
     &                  HK2,    T2,    RT2,    A2,  ACRIT, IDAMPV,
     &           AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &               AX_HK2, AX_T2, AX_RT2, AX_A2 )
C----------------------------------------------------------
C     Returns average amplification AX over interval 1..2
C----------------------------------------------------------
C
cC==========================
cC---- 1st-order -- based on "1" quantities only
c      CALL DAMPL( HK1, T1, RT1, AX1, AX1_HK1, AX1_T1, AX1_RT1 )
c      AX2_HK2 = 0.0
c      AX2_T2  = 0.0
c      AX2_RT2 = 0.0
cC
c      AX1_A1 = 0.0
c      AX2_A2 = 0.0
cC
c      AX     = AX1
c      AX_AX1 = 1.0
c      AX_AX2 = 0.0
cC
c      ARG = MIN( 20.0*(ACRIT-A1) , 20.0 )
c      EXN    = EXP(-ARG)
c      EXN_A1 = 20.0*EXN
c      EXN_A2 = 0.
cC
c      DAX    = EXN   * 0.0004/T1
c      DAX_A1 = EXN_A1* 0.0004/T1
c      DAX_A2 = 0.
c      DAX_T1 = -DAX/T1
c      DAX_T2 = 0.
C
C==========================
C---- 2nd-order
      IF(IDAMPV.EQ.0) THEN
       CALL DAMPL( HK1, T1, RT1, AX1, AX1_HK1, AX1_T1, AX1_RT1 )
       CALL DAMPL( HK2, T2, RT2, AX2, AX2_HK2, AX2_T2, AX2_RT2 )
      ELSE
       CALL DAMPL2( HK1, T1, RT1, AX1, AX1_HK1, AX1_T1, AX1_RT1 )
       CALL DAMPL2( HK2, T2, RT2, AX2, AX2_HK2, AX2_T2, AX2_RT2 )
      ENDIF
C
CC---- simple-average version
C      AXA = 0.5*(AX1 + AX2)
C      IF(AXA .LE. 0.0) THEN
C       AXA = 0.0
C       AXA_AX1 = 0.0
C       AXA_AX2 = 0.0
C      ELSE
C       AXA_AX1 = 0.5
C       AXA_AX2 = 0.5
C      ENDIF
C
C---- rms-average version (seems a little better on coarse grids)
      AXSQ = 0.5*(AX1**2 + AX2**2)
      IF(AXSQ .LE. 0.0) THEN
       AXA = 0.0
       AXA_AX1 = 0.0
       AXA_AX2 = 0.0
      ELSE
       AXA = SQRT(AXSQ)
       AXA_AX1 = 0.5*AX1/AXA
       AXA_AX2 = 0.5*AX2/AXA
      ENDIF
C
C----- small additional term to ensure  dN/dx > 0  near  N = Ncrit
       ARG = MIN( 20.0*(ACRIT-0.5*(A1+A2)) , 20.0 )
       IF(ARG.LE.0.0) THEN
        EXN    = 1.0
CC      EXN_AC = 0.
        EXN_A1 = 0.
        EXN_A2 = 0.
       ELSE
        EXN    = EXP(-ARG)
CC      EXN_AC = -20.0    *EXN
        EXN_A1 =  20.0*0.5*EXN
        EXN_A2 =  20.0*0.5*EXN
       ENDIF
C
       DAX    = EXN    * 0.002/(T1+T2)
CC     DAX_AC = EXN_AC * 0.002/(T1+T2)
       DAX_A1 = EXN_A1 * 0.002/(T1+T2)
       DAX_A2 = EXN_A2 * 0.002/(T1+T2)
       DAX_T1 = -DAX/(T1+T2)
       DAX_T2 = -DAX/(T1+T2)
C
c
c        DAX    = 0.
c        DAX_A1 = 0.
c        DAX_A2 = 0.
c        DAX_AC = 0.
c        DAX_T1 = 0.
c        DAX_T2 = 0.
C==========================
C
      AX     = AXA             + DAX
C
      AX_HK1 = AXA_AX1*AX1_HK1
      AX_T1  = AXA_AX1*AX1_T1  + DAX_T1
      AX_RT1 = AXA_AX1*AX1_RT1
      AX_A1  =                   DAX_A1
C
      AX_HK2 = AXA_AX2*AX2_HK2
      AX_T2  = AXA_AX2*AX2_T2  + DAX_T2
      AX_RT2 = AXA_AX2*AX2_RT2
      AX_A2  =                   DAX_A2
C
      RETURN
      END


c      SUBROUTINE TRCHEK1
cC-------------------------------------------------
cC     Checks if transition occurs in the current
cC     interval 1..2  (IBL-1...IBL) on side IS.
cC
cC     Old first-order version. 
cC
cC     Growth rate is evaluated at the upstream 
cC     point "1". The discrete amplification 
cC     equation is
cC
cC       Ncrit - N(X1)     
cC       -------------  =  N'(X1)
cC          XT - X1        
cC
cC     which can be immediately solved for 
cC     the transition location XT.
cC-------------------------------------------------
c      INCLUDE 'XBL.INC'
cC
cC---- calculate AMPL2 value
c      CALL AXSET( HK1,    T1,    RT1, AMPL1,
c     &            HK2,    T2,    RT2, AMPL2,  AMCRIT, IDAMPV,
c     &     AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
c     &         AX_HK2, AX_T2, AX_RT2, AX_A2 )
c      AMPL2 = AMPL1 + AX*(X2-X1)
cC
cC---- test for free or forced transition
c      TRFREE = AMPL2.GE.AMCRIT
c      TRFORC = XIFORC.GT.X1 .AND. XIFORC.LE.X2
cC
cC---- set transition interval flag
c      TRAN = TRFORC .OR. TRFREE
cC
cC---- if no transition yet, just return
c      IF(.NOT.TRAN) RETURN
cC
cC---- resolve if both forced and free transition
c      IF(TRFREE .AND. TRFORC) THEN
c       XT = (AMCRIT-AMPL1)/AX  +  X1
c       TRFORC = XIFORC .LT. XT
c       TRFREE = XIFORC .GE. XT
c      ENDIF
cC
c      IF(TRFORC) THEN
cC----- if forced transition, then XT is prescribed
c       XT = XIFORC
c       XT_A1 = 0.
c       XT_X1 = 0.
c       XT_T1 = 0.
c       XT_D1 = 0.
c       XT_U1 = 0.
c       XT_X2 = 0.
c       XT_T2 = 0.
c       XT_D2 = 0.
c       XT_U2 = 0.
c       XT_MS = 0.
c       XT_RE = 0.
c       XT_XF = 1.0
c      ELSE
cC----- if free transition, XT is related to BL variables
cC-     by the amplification equation
cC
c       XT    =  (AMCRIT-AMPL1)/AX     + X1
c       XT_AX = -(AMCRIT-AMPL1)/AX**2
cC
c       XT_A1 = -1.0/AX - (AMCRIT-AMPL1)/AX**2 * AX_A1
c       XT_X1 = 1.0
c       XT_T1 = XT_AX*(AX_HK1*HK1_T1 + AX_T1 + AX_RT1*RT1_T1)
c       XT_D1 = XT_AX*(AX_HK1*HK1_D1                        )
c       XT_U1 = XT_AX*(AX_HK1*HK1_U1         + AX_RT1*RT1_U1)
c       XT_X2 = 0.
c       XT_T2 = 0.
c       XT_D2 = 0.
c       XT_U2 = 0.
c       XT_MS = XT_AX*(AX_HK1*HK1_MS         + AX_RT1*RT1_MS)
c       XT_RE = XT_AX*(                        AX_RT1*RT1_RE)
c       XT_XF = 0.0
c      ENDIF
cC
c      RETURN
c      END
 
 
      SUBROUTINE TRCHEK2
C----------------------------------------------------------------
C     New second-order version:  December 1994.
C
C     Checks if transition occurs in the current interval X1..X2.
C     If transition occurs, then set transition location XT, and 
C     its sensitivities to "1" and "2" variables.  If no transition, 
C     set amplification AMPL2.
C
C
C     Solves the implicit amplification equation for N2:
C
C       N2 - N1     N'(XT,NT) + N'(X1,N1)
C       -------  =  ---------------------
C       X2 - X1               2
C
C     In effect, a 2-point central difference is used between
C     X1..X2 (no transition), or X1..XT (transition).  The switch
C     is done by defining XT,NT in the equation above depending
C     on whether N2 exceeds Ncrit.
C
C  If N2<Ncrit:  NT=N2    , XT=X2                  (no transition)
C
C  If N2>Ncrit:  NT=Ncrit , XT=(Ncrit-N1)/(N2-N1)  (transition)
C
C
C----------------------------------------------------------------
      INCLUDE 'XBL.INC'
      DATA DAEPS / 5.0E-5 /
CCC   DATA DAEPS / 1.0D-12 /
C
C---- save variables and sensitivities at IBL ("2") for future restoration
      DO 5 ICOM=1, NCOM
        C2SAV(ICOM) = COM2(ICOM)
    5 CONTINUE
C
C---- calculate average amplification rate AX over X1..X2 interval
      CALL AXSET( HK1,    T1,    RT1, AMPL1,
     &            HK2,    T2,    RT2, AMPL2,  AMCRIT, IDAMPV,
     &     AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &         AX_HK2, AX_T2, AX_RT2, AX_A2 )
C
C---- set initial guess for iterate N2 (AMPL2) at X2
      AMPL2 = AMPL1 + AX*(X2-X1)
C
C---- solve implicit system for amplification AMPL2
      DO 100 ITAM=1, 30
C
C---- define weighting factors WF1,WF2 for defining "T" quantities from 1,2
C
      IF(AMPL2 .LE. AMCRIT) THEN
C------ there is no transition yet,  "T" is the same as "2"
        AMPLT    = AMPL2
        AMPLT_A2 = 1.0
        SFA    = 1.0
        SFA_A1 = 0.
        SFA_A2 = 0.
      ELSE
C------ there is transition in X1..X2, "T" is set from N1, N2
        AMPLT    = AMCRIT
        AMPLT_A2 = 0.
        SFA    = (AMPLT - AMPL1)/(AMPL2-AMPL1)
        SFA_A1 = ( SFA  - 1.0  )/(AMPL2-AMPL1)
        SFA_A2 = (      - SFA  )/(AMPL2-AMPL1)
      ENDIF
C
      IF(XIFORC.LT.X2) THEN
        SFX    = (XIFORC - X1 )/(X2-X1)
        SFX_X1 = (SFX    - 1.0)/(X2-X1)
        SFX_X2 = (       - SFX)/(X2-X1)
        SFX_XF =  1.0          /(X2-X1)
      ELSE
        SFX    = 1.0
        SFX_X1 = 0.
        SFX_X2 = 0.
        SFX_XF = 0.
      ENDIF
C
C---- set weighting factor from free or forced transition
      IF(SFA.LT.SFX) THEN
        WF2    = SFA
        WF2_A1 = SFA_A1
        WF2_A2 = SFA_A2
        WF2_X1 = 0.
        WF2_X2 = 0.
        WF2_XF = 0.
      ELSE
        WF2    = SFX
        WF2_A1 = 0.
        WF2_A2 = 0.
        WF2_X1 = SFX_X1
        WF2_X2 = SFX_X2
        WF2_XF = SFX_XF
      ENDIF
C
C
C=====================
CC---- 1st-order (based on "1" quantites only, for testing)
C      WF2    = 0.0
C      WF2_A1 = 0.0
C      WF2_A2 = 0.0
C      WF2_X1 = 0.0
C      WF2_X2 = 0.0
C      WF2_XF = 0.0
C=====================
C
      WF1    = 1.0 - WF2
      WF1_A1 =     - WF2_A1
      WF1_A2 =     - WF2_A2
      WF1_X1 =     - WF2_X1
      WF1_X2 =     - WF2_X2
      WF1_XF =     - WF2_XF
C
C---- interpolate BL variables to XT
      XT    = X1*WF1    + X2*WF2
      TT    = T1*WF1    + T2*WF2
      DT    = D1*WF1    + D2*WF2
      UT    = U1*WF1    + U2*WF2
C
      XT_A2 = X1*WF1_A2 + X2*WF2_A2
      TT_A2 = T1*WF1_A2 + T2*WF2_A2
      DT_A2 = D1*WF1_A2 + D2*WF2_A2
      UT_A2 = U1*WF1_A2 + U2*WF2_A2
C
C---- temporarily set "2" variables from "T" for BLKIN
      X2 = XT
      T2 = TT
      D2 = DT
      U2 = UT
C
C---- calculate laminar secondary "T" variables HKT, RTT
      CALL BLKIN
C
      HKT    = HK2
      HKT_TT = HK2_T2
      HKT_DT = HK2_D2
      HKT_UT = HK2_U2
      HKT_MS = HK2_MS
C
      RTT    = RT2
      RTT_TT = RT2_T2
      RTT_UT = RT2_U2
      RTT_MS = RT2_MS
      RTT_RE = RT2_RE
C
C---- restore clobbered "2" variables, except for AMPL2
      AMSAVE = AMPL2
      DO 8 ICOM=1, NCOM
        COM2(ICOM) = C2SAV(ICOM)
 8    CONTINUE
      AMPL2 = AMSAVE
C
C---- calculate amplification rate AX over current X1-XT interval
      CALL AXSET( HK1,    T1,    RT1, AMPL1,
     &            HKT,    TT,    RTT, AMPLT,  AMCRIT, IDAMPV,
     &     AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &         AX_HKT, AX_TT, AX_RTT, AX_AT )
C
C---- punch out early if there is no amplification here
      IF(AX .LE. 0.0) GO TO 101
C
C---- set sensitivity of AX(A2)
      AX_A2 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_A2
     &      + (AX_HKT*HKT_DT                        )*DT_A2
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_A2
     &      +  AX_AT                                 *AMPLT_A2
C
C---- residual for implicit AMPL2 definition (amplification equation)
      RES    = AMPL2 - AMPL1 - AX   *(X2-X1) 
      RES_A2 = 1.0           - AX_A2*(X2-X1)
C
      DA2 = -RES/RES_A2
C
      RLX = 1.0
      DXT = XT_A2*DA2
C
      IF(RLX*ABS(DXT/(X2-X1)) .GT. 0.05) RLX = 0.05*ABS((X2-X1)/DXT)
      IF(RLX*ABS(DA2)         .GT. 1.0 ) RLX = 1.0 *ABS(   1.0 /DA2)
C
C---- check if converged
      IF(ABS(DA2) .LT. DAEPS) GO TO 101
C
      IF((AMPL2.GT.AMCRIT .AND. AMPL2+RLX*DA2.LT.AMCRIT).OR.
     &   (AMPL2.LT.AMCRIT .AND. AMPL2+RLX*DA2.GT.AMCRIT)    ) THEN
C------ limited Newton step so AMPL2 doesn't step across AMCRIT either way
        AMPL2 = AMCRIT
      ELSE
C------ regular Newton step
        AMPL2 = AMPL2 + RLX*DA2
      ENDIF
C
 100  CONTINUE
      WRITE(*,*) 'TRCHEK2: N2 convergence failed.'
      WRITE(*,6700) X1, XT, X2, AMPL1, AMPLT, AMPL2, AX, DA2
 6700 FORMAT(1X,'x:', 3F9.5,'  N:',3F7.3,'  Nx:',F8.3,'   dN:',E10.3)
C
 101  CONTINUE
C
C
C---- test for free or forced transition
      TRFREE = AMPL2 .GE. AMCRIT
      TRFORC = XIFORC.GT.X1 .AND. XIFORC.LE.X2
C
C---- set transition interval flag
      TRAN = TRFORC .OR. TRFREE
C
      IF(.NOT.TRAN) RETURN
C
C---- resolve if both forced and free transition
      IF(TRFREE .AND. TRFORC) THEN
       TRFORC = XIFORC .LT. XT
       TRFREE = XIFORC .GE. XT
      ENDIF
C
      IF(TRFORC) THEN
C----- if forced transition, then XT is prescribed,
C-     no sense calculating the sensitivities, since we know them...
       XT = XIFORC
       XT_A1 = 0.
       XT_X1 = 0.
       XT_T1 = 0.
       XT_D1 = 0.
       XT_U1 = 0.
       XT_X2 = 0.
       XT_T2 = 0.
       XT_D2 = 0.
       XT_U2 = 0.
       XT_MS = 0.
       XT_RE = 0.
       XT_XF = 1.0
       RETURN
      ENDIF
C
C---- free transition ... set sensitivities of XT
C
C---- XT( X1 X2 A1 A2 XF ),  TT( T1 T2 A1 A2 X1 X2 XF),   DT( ...
CC    XT    = X1*WF1    + X2*WF2
CC    TT    = T1*WF1    + T2*WF2
CC    DT    = D1*WF1    + D2*WF2
CC    UT    = U1*WF1    + U2*WF2
C
      XT_X1 =    WF1
      TT_T1 =    WF1
      DT_D1 =    WF1
      UT_U1 =    WF1
C
      XT_X2 =                WF2
      TT_T2 =                WF2
      DT_D2 =                WF2
      UT_U2 =                WF2
C
      XT_A1 = X1*WF1_A1 + X2*WF2_A1
      TT_A1 = T1*WF1_A1 + T2*WF2_A1
      DT_A1 = D1*WF1_A1 + D2*WF2_A1
      UT_A1 = U1*WF1_A1 + U2*WF2_A1
C
CC    XT_A2 = X1*WF1_A2 + X2*WF2_A2
CC    TT_A2 = T1*WF1_A2 + T2*WF2_A2
CC    DT_A2 = D1*WF1_A2 + D2*WF2_A2
CC    UT_A2 = U1*WF1_A2 + U2*WF2_A2
C
      XT_X1 = X1*WF1_X1 + X2*WF2_X1 + XT_X1
      TT_X1 = T1*WF1_X1 + T2*WF2_X1
      DT_X1 = D1*WF1_X1 + D2*WF2_X1
      UT_X1 = U1*WF1_X1 + U2*WF2_X1
C
      XT_X2 = X1*WF1_X2 + X2*WF2_X2 + XT_X2
      TT_X2 = T1*WF1_X2 + T2*WF2_X2
      DT_X2 = D1*WF1_X2 + D2*WF2_X2
      UT_X2 = U1*WF1_X2 + U2*WF2_X2
C
      XT_XF = X1*WF1_XF + X2*WF2_XF
      TT_XF = T1*WF1_XF + T2*WF2_XF
      DT_XF = D1*WF1_XF + D2*WF2_XF
      UT_XF = U1*WF1_XF + U2*WF2_XF
C
C---- at this point, AX = AX( HK1, T1, RT1, A1, HKT, TT, RTT, AT )
C
C---- set sensitivities of AX( T1 D1 U1 A1 T2 D2 U2 A2 MS RE )
      AX_T1 =  AX_HK1*HK1_T1 + AX_T1 + AX_RT1*RT1_T1
     &      + (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_T1
      AX_D1 =  AX_HK1*HK1_D1
     &      + (AX_HKT*HKT_DT                        )*DT_D1
      AX_U1 =  AX_HK1*HK1_U1         + AX_RT1*RT1_U1
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_U1
      AX_A1 =  AX_A1
     &      + (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_A1
     &      + (AX_HKT*HKT_DT                        )*DT_A1
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_A1
      AX_X1 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_X1
     &      + (AX_HKT*HKT_DT                        )*DT_X1
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_X1
C
      AX_T2 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_T2
      AX_D2 = (AX_HKT*HKT_DT                        )*DT_D2
      AX_U2 = (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_U2
      AX_A2 =  AX_AT                                 *AMPLT_A2
     &      + (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_A2
     &      + (AX_HKT*HKT_DT                        )*DT_A2
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_A2
      AX_X2 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_X2
     &      + (AX_HKT*HKT_DT                        )*DT_X2
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_X2
C
      AX_XF = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_XF
     &      + (AX_HKT*HKT_DT                        )*DT_XF
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_XF
C
      AX_MS =  AX_HKT*HKT_MS         + AX_RTT*RTT_MS
     &      +  AX_HK1*HK1_MS         + AX_RT1*RT1_MS
      AX_RE =                          AX_RTT*RTT_RE
     &                               + AX_RT1*RT1_RE
C
C
C---- set sensitivities of residual RES
CCC   RES  = AMPL2 - AMPL1 - AX*(X2-X1)
      Z_AX =               -    (X2-X1)
C
      Z_A1 = Z_AX*AX_A1 - 1.0
      Z_T1 = Z_AX*AX_T1
      Z_D1 = Z_AX*AX_D1
      Z_U1 = Z_AX*AX_U1
      Z_X1 = Z_AX*AX_X1 + AX
C
      Z_A2 = Z_AX*AX_A2 + 1.0
      Z_T2 = Z_AX*AX_T2
      Z_D2 = Z_AX*AX_D2
      Z_U2 = Z_AX*AX_U2
      Z_X2 = Z_AX*AX_X2 - AX
C
      Z_XF = Z_AX*AX_XF
      Z_MS = Z_AX*AX_MS
      Z_RE = Z_AX*AX_RE
C
C---- set sensitivities of XT, with RES being stationary for A2 constraint
      XT_A1 = XT_A1 - (XT_A2/Z_A2)*Z_A1
      XT_T1 =       - (XT_A2/Z_A2)*Z_T1
      XT_D1 =       - (XT_A2/Z_A2)*Z_D1
      XT_U1 =       - (XT_A2/Z_A2)*Z_U1
      XT_X1 = XT_X1 - (XT_A2/Z_A2)*Z_X1
      XT_T2 =       - (XT_A2/Z_A2)*Z_T2
      XT_D2 =       - (XT_A2/Z_A2)*Z_D2
      XT_U2 =       - (XT_A2/Z_A2)*Z_U2
      XT_X2 = XT_X2 - (XT_A2/Z_A2)*Z_X2
      XT_MS =       - (XT_A2/Z_A2)*Z_MS
      XT_RE =       - (XT_A2/Z_A2)*Z_RE
      XT_XF = 0.0
C
      RETURN
      END


      SUBROUTINE BLSYS
C------------------------------------------------------------------
C
C     Sets up the BL Newton system governing the current interval:
C
C     |       ||dA1|     |       ||dA2|       |     |
C     |  VS1  ||dT1|  +  |  VS2  ||dT2|   =   |VSREZ|
C     |       ||dD1|     |       ||dD2|       |     |
C              |dU1|              |dU2|
C              |dX1|              |dX2|
C
C        3x5    5x1         3x5    5x1          3x1
C
C     The system as shown corresponds to a laminar station
C     If TRAN, then  dS2  replaces  dA2
C     If TURB, then  dS1, dS2  replace  dA1, dA2
C
C------------------------------------------------------------------
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
C
C---- calculate secondary BL variables and their sensitivities
      IF(WAKE) THEN
       CALL BLVAR(3)
       CALL BLMID(3)
      ELSE IF(TURB.OR.TRAN) THEN
       CALL BLVAR(2)
       CALL BLMID(2)
      ELSE
       CALL BLVAR(1)
       CALL BLMID(1)
      ENDIF
C
C---- for the similarity station, "1" and "2" variables are the same
      IF(SIMI) THEN
       DO 3 ICOM=1, NCOM
         COM1(ICOM) = COM2(ICOM)
    3  CONTINUE
      ENDIF
C
C---- set up appropriate finite difference system for current interval
      IF(TRAN) THEN
       CALL TRDIF
      ELSE IF(SIMI) THEN
       CALL BLDIF(0)
      ELSE IF(.NOT.TURB) THEN
       CALL BLDIF(1)
      ELSE IF(WAKE) THEN
       CALL BLDIF(3)
      ELSE IF(TURB) THEN
       CALL BLDIF(2)
      ENDIF
C
      IF(SIMI) THEN
C----- at similarity station, "1" variables are really "2" variables
       DO 10 K=1, 4
         DO 101 L=1, 5
           VS2(K,L) = VS1(K,L) + VS2(K,L)
           VS1(K,L) = 0.
  101    CONTINUE
   10  CONTINUE
      ENDIF
C
C---- change system over into incompressible Uei and Mach
      DO 20 K=1, 4
C
C------ residual derivatives wrt compressible Uec
        RES_U1 = VS1(K,4)
        RES_U2 = VS2(K,4)
        RES_MS = VSM(K)
C
C------ combine with derivatives of compressible  U1,U2 = Uec(Uei M)
        VS1(K,4) = RES_U1*U1_UEI
        VS2(K,4) =                RES_U2*U2_UEI
        VSM(K)   = RES_U1*U1_MS + RES_U2*U2_MS  + RES_MS
   20 CONTINUE
C
      RETURN
      END
 

      SUBROUTINE TESYS(CTE,TTE,DTE)
C--------------------------------------------------------
C     Sets up "dummy" BL system between airfoil TE point 
C     and first wake point infinitesimally behind TE.
C--------------------------------------------------------
      IMPLICIT REAL (M)
      INCLUDE 'XBL.INC'
C
      DO 55 K=1, 4
        VSREZ(K) = 0.
        VSM(K)   = 0.
        VSR(K)   = 0.
        VSX(K)   = 0.
        DO 551 L=1, 5
          VS1(K,L) = 0.
          VS2(K,L) = 0.
  551   CONTINUE
   55 CONTINUE
C
      CALL BLVAR(3)
C
      VS1(1,1) = -1.0
      VS2(1,1) = 1.0
      VSREZ(1) = CTE - S2      
C
      VS1(2,2) = -1.0
      VS2(2,2) = 1.0
      VSREZ(2) = TTE - T2
C
      VS1(3,3) = -1.0
      VS2(3,3) = 1.0
      VSREZ(3) = DTE - D2 - DW2
C
      RETURN
      END


      SUBROUTINE BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
C----------------------------------------------------------
C     Set BL primary "2" variables from parameter list
C----------------------------------------------------------
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
C
      X2 = XSI
      AMPL2 = AMI
      S2  = CTI
      T2  = THI
      D2  = DSI - DSWAKI
      DW2 = DSWAKI
C
      U2 = UEI*(1.0-TKBL) / (1.0 - TKBL*(UEI/QINFBL)**2)
      U2_UEI = (1.0 + TKBL*(2.0*U2*UEI/QINFBL**2 - 1.0))
     &       / (1.0 - TKBL*(UEI/QINFBL)**2)
      U2_MS  = (U2*(UEI/QINFBL)**2  -  UEI)*TKBL_MS
     &                    / (1.0 - TKBL*(UEI/QINFBL)**2)
C
      RETURN
      END ! BLPRV

 
      SUBROUTINE BLKIN
C----------------------------------------------------------
C     Calculates turbulence-independent secondary "2" 
C     variables from the primary "2" variables.
C----------------------------------------------------------
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
C
C---- set edge Mach number ** 2
      M2    = U2*U2*HSTINV / (GM1BL*(1.0 - 0.5*U2*U2*HSTINV))
      TR2   = 1.0 + 0.5*GM1BL*M2
      M2_U2 = 2.0*M2*TR2/U2
      M2_MS = U2*U2*TR2    / (GM1BL*(1.0 - 0.5*U2*U2*HSTINV))
     &      * HSTINV_MS
C
C---- set edge static density (isentropic relation)
      R2    = RSTBL   *TR2**(-1.0/GM1BL)
      R2_U2 = -R2/TR2 * 0.5*M2_U2
      R2_MS = -R2/TR2 * 0.5*M2_MS
     &      + RSTBL_MS*TR2**(-1.0/GM1BL)
C
C---- set shape parameter
      H2    =  D2/T2
      H2_D2 = 1.0/T2
      H2_T2 = -H2/T2
C
C---- set edge static/stagnation enthalpy
      HERAT = 1.0 - 0.5*U2*U2*HSTINV
      HE_U2 =     -        U2*HSTINV
      HE_MS =     - 0.5*U2*U2*HSTINV_MS
C
C---- set molecular viscosity
      V2 = SQRT((HERAT)**3) * (1.0+HVRAT)/(HERAT+HVRAT)/REYBL
      V2_HE = V2*(1.5/HERAT - 1.0/(HERAT+HVRAT))
C
      V2_U2 =                        V2_HE*HE_U2
      V2_MS = -V2/REYBL * REYBL_MS + V2_HE*HE_MS
      V2_RE = -V2/REYBL * REYBL_RE
C
C---- set kinematic shape parameter
      CALL HKIN( H2, M2, HK2, HK2_H2, HK2_M2 )
C
      HK2_U2 =                HK2_M2*M2_U2
      HK2_T2 = HK2_H2*H2_T2
      HK2_D2 = HK2_H2*H2_D2
      HK2_MS =                HK2_M2*M2_MS
C
C---- set momentum thickness Reynolds number
      RT2    = R2*U2*T2/V2
      RT2_U2 = RT2*(1.0/U2 + R2_U2/R2 - V2_U2/V2)
      RT2_T2 = RT2/T2
      RT2_MS = RT2*(         R2_MS/R2 - V2_MS/V2)
      RT2_RE = RT2*(                  - V2_RE/V2)
C
      RETURN
      END ! BLKIN


 
      SUBROUTINE BLVAR(ITYP)
C----------------------------------------------------
C     Calculates all secondary "2" variables from
C     the primary "2" variables X2, U2, T2, D2, S2.
C     Also calculates the sensitivities of the
C     secondary variables wrt the primary variables.
C
C      ITYP = 1 :  laminar
C      ITYP = 2 :  turbulent
C      ITYP = 3 :  turbulent wake
C----------------------------------------------------
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
C
      IF(ITYP.EQ.3) HK2 = MAX(HK2,1.00005)
      IF(ITYP.NE.3) HK2 = MAX(HK2,1.05000)
C
C---- density thickness shape parameter     ( H** )
      CALL HCT( HK2, M2, HC2, HC2_HK2, HC2_M2 )
      HC2_U2 = HC2_HK2*HK2_U2 + HC2_M2*M2_U2
      HC2_T2 = HC2_HK2*HK2_T2
      HC2_D2 = HC2_HK2*HK2_D2
      HC2_MS = HC2_HK2*HK2_MS + HC2_M2*M2_MS
C
C---- set KE thickness shape parameter from  H - H*  correlations
      IF(ITYP.EQ.1) THEN
       CALL HSL( HK2, RT2, M2, HS2, HS2_HK2, HS2_RT2, HS2_M2 )
      ELSE
       CALL HST( HK2, RT2, M2, HS2, HS2_HK2, HS2_RT2, HS2_M2 )
      ENDIF
C
      HS2_U2 = HS2_HK2*HK2_U2 + HS2_RT2*RT2_U2 + HS2_M2*M2_U2
      HS2_T2 = HS2_HK2*HK2_T2 + HS2_RT2*RT2_T2
      HS2_D2 = HS2_HK2*HK2_D2
      HS2_MS = HS2_HK2*HK2_MS + HS2_RT2*RT2_MS + HS2_M2*M2_MS
      HS2_RE =                  HS2_RT2*RT2_RE
C
C---- normalized slip velocity  Us
      US2     = 0.5*HS2*( 1.0 - (HK2-1.0)/(GBCON*H2) )
      US2_HS2 = 0.5  *  ( 1.0 - (HK2-1.0)/(GBCON*H2) )
      US2_HK2 = 0.5*HS2*(     -  1.0     /(GBCON*H2) )
      US2_H2  = 0.5*HS2*        (HK2-1.0)/(GBCON*H2**2)
C
      US2_U2 = US2_HS2*HS2_U2 + US2_HK2*HK2_U2
      US2_T2 = US2_HS2*HS2_T2 + US2_HK2*HK2_T2 + US2_H2*H2_T2
      US2_D2 = US2_HS2*HS2_D2 + US2_HK2*HK2_D2 + US2_H2*H2_D2
      US2_MS = US2_HS2*HS2_MS + US2_HK2*HK2_MS
      US2_RE = US2_HS2*HS2_RE
C
      IF(ITYP.LE.2 .AND. US2.GT.0.95) THEN
CCC       WRITE(*,*) 'BLVAR: Us clamped:', US2
       US2 = 0.98
       US2_U2 = 0.
       US2_T2 = 0.
       US2_D2 = 0.
       US2_MS = 0.
       US2_RE = 0.
      ENDIF
C
      IF(ITYP.EQ.3 .AND. US2.GT.0.99995) THEN
CCC       WRITE(*,*) 'BLVAR: Wake Us clamped:', US2
       US2 = 0.99995
       US2_U2 = 0.
       US2_T2 = 0.
       US2_D2 = 0.
       US2_MS = 0.
       US2_RE = 0.
      ENDIF
C
C---- equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2
C   ...  NEW  12 Oct 94
      GCC = 0.0
      HKC = HK2 - 1.0
      HKC_HK2 = 1.0
      HKC_RT2 = 0.0
      IF(ITYP.EQ.2) THEN
       GCC = GCCON
       HKC     = HK2 - 1.0 - GCC/RT2
       HKC_HK2 = 1.0
       HKC_RT2 =             GCC/RT2**2
       IF(HKC .LT. 0.01) THEN
        HKC = 0.01
        HKC_HK2 = 0.0
        HKC_RT2 = 0.0
       ENDIF
      ENDIF
C
      HKB = HK2 - 1.0
      USB = 1.0 - US2
      CQ2     =
     &    SQRT( CTCON*HS2*HKB*HKC**2 / (USB*H2*HK2**2) )
      CQ2_HS2 = CTCON    *HKB*HKC**2 / (USB*H2*HK2**2)       * 0.5/CQ2
      CQ2_US2 = CTCON*HS2*HKB*HKC**2 / (USB*H2*HK2**2) / USB * 0.5/CQ2
      CQ2_HK2 = CTCON*HS2    *HKC**2 / (USB*H2*HK2**2)       * 0.5/CQ2
     &        - CTCON*HS2*HKB*HKC**2 / (USB*H2*HK2**3) * 2.0 * 0.5/CQ2
     &        + CTCON*HS2*HKB*HKC    / (USB*H2*HK2**2) * 2.0 * 0.5/CQ2
     &         *HKC_HK2
      CQ2_RT2 = CTCON*HS2*HKB*HKC    / (USB*H2*HK2**2) * 2.0 * 0.5/CQ2
     &         *HKC_RT2
      CQ2_H2  =-CTCON*HS2*HKB*HKC**2 / (USB*H2*HK2**2) / H2  * 0.5/CQ2
C
      CQ2_U2 = CQ2_HS2*HS2_U2 + CQ2_US2*US2_U2 + CQ2_HK2*HK2_U2
      CQ2_T2 = CQ2_HS2*HS2_T2 + CQ2_US2*US2_T2 + CQ2_HK2*HK2_T2
      CQ2_D2 = CQ2_HS2*HS2_D2 + CQ2_US2*US2_D2 + CQ2_HK2*HK2_D2
      CQ2_MS = CQ2_HS2*HS2_MS + CQ2_US2*US2_MS + CQ2_HK2*HK2_MS
      CQ2_RE = CQ2_HS2*HS2_RE + CQ2_US2*US2_RE
C
      CQ2_U2 = CQ2_U2                + CQ2_RT2*RT2_U2
      CQ2_T2 = CQ2_T2 + CQ2_H2*H2_T2 + CQ2_RT2*RT2_T2
      CQ2_D2 = CQ2_D2 + CQ2_H2*H2_D2
      CQ2_MS = CQ2_MS                + CQ2_RT2*RT2_MS
      CQ2_RE = CQ2_RE                + CQ2_RT2*RT2_RE
C
C
C---- set skin friction coefficient 
      IF(ITYP.EQ.3) THEN
C----- wake
       CF2     = 0.
       CF2_HK2 = 0.
       CF2_RT2 = 0.
       CF2_M2  = 0.
      ELSE IF(ITYP.EQ.1) THEN
C----- laminar
       CALL CFL( HK2, RT2, M2, CF2, CF2_HK2, CF2_RT2, CF2_M2 )
      ELSE
C----- turbulent
       CALL CFT( HK2, RT2, M2, CF2, CF2_HK2, CF2_RT2, CF2_M2 )
       CALL CFL( HK2, RT2, M2, CF2L,CF2L_HK2,CF2L_RT2,CF2L_M2)
       IF(CF2L.GT.CF2) THEN
C------- laminar Cf is greater than turbulent Cf -- use laminar
C-       (this will only occur for unreasonably small Rtheta)
ccc      write(*,*) 'Cft Cfl Rt Hk:', CF2, CF2L, RT2, HK2, X2
         CF2     = CF2L
         CF2_HK2 = CF2L_HK2
         CF2_RT2 = CF2L_RT2
         CF2_M2  = CF2L_M2
       ENDIF
      ENDIF
C
      CF2_U2 = CF2_HK2*HK2_U2 + CF2_RT2*RT2_U2 + CF2_M2*M2_U2
      CF2_T2 = CF2_HK2*HK2_T2 + CF2_RT2*RT2_T2
      CF2_D2 = CF2_HK2*HK2_D2
      CF2_MS = CF2_HK2*HK2_MS + CF2_RT2*RT2_MS + CF2_M2*M2_MS
      CF2_RE =                  CF2_RT2*RT2_RE
C
C---- dissipation function    2 CD / H*
      IF(ITYP.EQ.1) THEN
C
C----- laminar
       CALL DIL( HK2, RT2, DI2, DI2_HK2, DI2_RT2 )
C
       DI2_U2 = DI2_HK2*HK2_U2 + DI2_RT2*RT2_U2
       DI2_T2 = DI2_HK2*HK2_T2 + DI2_RT2*RT2_T2
       DI2_D2 = DI2_HK2*HK2_D2
       DI2_S2 = 0.
       DI2_MS = DI2_HK2*HK2_MS + DI2_RT2*RT2_MS
       DI2_RE =                  DI2_RT2*RT2_RE
C
      ELSE IF(ITYP.EQ.2) THEN
C
CCC       CALL DIT(     HS2,     US2,     CF2,     S2, DI2,
CCC     &           DI2_HS2, DI2_US2, DI2_CF2, DI2_S2      )
C
C----- turbulent wall contribution
       CALL CFT(HK2, RT2, M2, CF2T, CF2T_HK2, CF2T_RT2, CF2T_M2)
       CF2T_U2 = CF2T_HK2*HK2_U2 + CF2T_RT2*RT2_U2 + CF2T_M2*M2_U2
       CF2T_T2 = CF2T_HK2*HK2_T2 + CF2T_RT2*RT2_T2
       CF2T_D2 = CF2T_HK2*HK2_D2
       CF2T_MS = CF2T_HK2*HK2_MS + CF2T_RT2*RT2_MS + CF2T_M2*M2_MS
       CF2T_RE =                   CF2T_RT2*RT2_RE
C
       DI2      =  ( 0.5*CF2T*US2 ) * 2.0/HS2
       DI2_HS2  = -( 0.5*CF2T*US2 ) * 2.0/HS2**2
       DI2_US2  =  ( 0.5*CF2T     ) * 2.0/HS2
       DI2_CF2T =  ( 0.5     *US2 ) * 2.0/HS2
C
       DI2_S2 = 0.0
       DI2_U2 = DI2_HS2*HS2_U2 + DI2_US2*US2_U2 + DI2_CF2T*CF2T_U2
       DI2_T2 = DI2_HS2*HS2_T2 + DI2_US2*US2_T2 + DI2_CF2T*CF2T_T2
       DI2_D2 = DI2_HS2*HS2_D2 + DI2_US2*US2_D2 + DI2_CF2T*CF2T_D2
       DI2_MS = DI2_HS2*HS2_MS + DI2_US2*US2_MS + DI2_CF2T*CF2T_MS
       DI2_RE = DI2_HS2*HS2_RE + DI2_US2*US2_RE + DI2_CF2T*CF2T_RE
C
C
C----- set minimum Hk for wake layer to still exist
       GRT = LOG(RT2)
       HMIN = 1.0 + 2.1/GRT
       HM_RT2 = -(2.1/GRT**2) / RT2
C
C----- set factor DFAC for correcting wall dissipation for very low Hk
       FL = (HK2-1.0)/(HMIN-1.0)
       FL_HK2 =   1.0/(HMIN-1.0)
       FL_RT2 = ( -FL/(HMIN-1.0) ) * HM_RT2
C
       TFL = TANH(FL)
       DFAC  = 0.5 + 0.5* TFL
       DF_FL =       0.5*(1.0 - TFL**2)
C
       DF_HK2 = DF_FL*FL_HK2
       DF_RT2 = DF_FL*FL_RT2
C
       DI2_S2 = DI2_S2*DFAC
       DI2_U2 = DI2_U2*DFAC + DI2*(DF_HK2*HK2_U2 + DF_RT2*RT2_U2)
       DI2_T2 = DI2_T2*DFAC + DI2*(DF_HK2*HK2_T2 + DF_RT2*RT2_T2)
       DI2_D2 = DI2_D2*DFAC + DI2*(DF_HK2*HK2_D2                )
       DI2_MS = DI2_MS*DFAC + DI2*(DF_HK2*HK2_MS + DF_RT2*RT2_MS)
       DI2_RE = DI2_RE*DFAC + DI2*(                DF_RT2*RT2_RE)
       DI2    = DI2   *DFAC
C
      ELSE
C
C----- zero wall contribution for wake
       DI2    = 0.0
       DI2_S2 = 0.0
       DI2_U2 = 0.0
       DI2_T2 = 0.0
       DI2_D2 = 0.0
       DI2_MS = 0.0
       DI2_RE = 0.0
C
      ENDIF
C
C
C---- Add on turbulent outer layer contribution
      IF(ITYP.NE.1) THEN
C
       DD     =  S2**2 * (0.995-US2) * 2.0/HS2
       DD_HS2 = -S2**2 * (0.995-US2) * 2.0/HS2**2
       DD_US2 = -S2**2               * 2.0/HS2
       DD_S2  =  S2*2.0* (0.995-US2) * 2.0/HS2
C
       DI2    = DI2    + DD
       DI2_S2 =          DD_S2
       DI2_U2 = DI2_U2 + DD_HS2*HS2_U2 + DD_US2*US2_U2
       DI2_T2 = DI2_T2 + DD_HS2*HS2_T2 + DD_US2*US2_T2
       DI2_D2 = DI2_D2 + DD_HS2*HS2_D2 + DD_US2*US2_D2
       DI2_MS = DI2_MS + DD_HS2*HS2_MS + DD_US2*US2_MS
       DI2_RE = DI2_RE + DD_HS2*HS2_RE + DD_US2*US2_RE
C
C----- add laminar stress contribution to outer layer CD
c###
       DD     =  0.15*(0.995-US2)**2 / RT2  * 2.0/HS2
       DD_US2 = -0.15*(0.995-US2)*2. / RT2  * 2.0/HS2
       DD_HS2 = -DD/HS2
       DD_RT2 = -DD/RT2
C
       DI2    = DI2    + DD
       DI2_U2 = DI2_U2 + DD_HS2*HS2_U2 + DD_US2*US2_U2 + DD_RT2*RT2_U2
       DI2_T2 = DI2_T2 + DD_HS2*HS2_T2 + DD_US2*US2_T2 + DD_RT2*RT2_T2
       DI2_D2 = DI2_D2 + DD_HS2*HS2_D2 + DD_US2*US2_D2
       DI2_MS = DI2_MS + DD_HS2*HS2_MS + DD_US2*US2_MS + DD_RT2*RT2_MS
       DI2_RE = DI2_RE + DD_HS2*HS2_RE + DD_US2*US2_RE + DD_RT2*RT2_RE
C
      ENDIF
C
C
      IF(ITYP.EQ.2) THEN
        CALL DIL( HK2, RT2, DI2L, DI2L_HK2, DI2L_RT2 )
C
        IF(DI2L.GT.DI2) THEN
C------- laminar CD is greater than turbulent CD -- use laminar
C-       (this will only occur for unreasonably small Rtheta)
ccc       write(*,*) 'CDt CDl Rt Hk:', DI2, DI2L, RT2, HK2
          DI2    = DI2L
          DI2_S2 = 0.
          DI2_U2 = DI2L_HK2*HK2_U2 + DI2L_RT2*RT2_U2
          DI2_T2 = DI2L_HK2*HK2_T2 + DI2L_RT2*RT2_T2
          DI2_D2 = DI2L_HK2*HK2_D2
          DI2_MS = DI2L_HK2*HK2_MS + DI2L_RT2*RT2_MS
          DI2_RE =                   DI2L_RT2*RT2_RE
        ENDIF
      ENDIF
C
cC----- add on CD contribution of inner shear layer
c       IF(ITYP.EQ.3 .AND. DW2.GT.0.0) THEN
c        DKON = 0.03*0.75**3
c        DDI = DKON*US2**3
c        DDI_US2 = 3.0*DKON*US2**2
c        DI2 = DI2 + DDI * DW2/DWTE
c        DI2_U2 = DI2_U2 + DDI_US2*US2_U2 * DW2/DWTE
c        DI2_T2 = DI2_T2 + DDI_US2*US2_T2 * DW2/DWTE
c        DI2_D2 = DI2_D2 + DDI_US2*US2_D2 * DW2/DWTE
c        DI2_MS = DI2_MS + DDI_US2*US2_MS * DW2/DWTE
c        DI2_RE = DI2_RE + DDI_US2*US2_RE * DW2/DWTE
c       ENDIF
C
      IF(ITYP.EQ.3) THEN
C------ laminar wake CD
        CALL DILW( HK2, RT2, DI2L, DI2L_HK2, DI2L_RT2 )
        IF(DI2L .GT. DI2) THEN
cc        IF(.true.) THEN
C------- laminar wake CD is greater than turbulent CD -- use laminar
C-       (this will only occur for unreasonably small Rtheta)
ccc         write(*,*) 'CDt CDl Rt Hk:', DI2, DI2L, RT2, HK2
         DI2    = DI2L
         DI2_S2 = 0.
         DI2_U2 = DI2L_HK2*HK2_U2 + DI2L_RT2*RT2_U2
         DI2_T2 = DI2L_HK2*HK2_T2 + DI2L_RT2*RT2_T2
         DI2_D2 = DI2L_HK2*HK2_D2
         DI2_MS = DI2L_HK2*HK2_MS + DI2L_RT2*RT2_MS
         DI2_RE =                   DI2L_RT2*RT2_RE
        ENDIF
      ENDIF
C
C
      IF(ITYP.EQ.3) THEN
C----- double dissipation for the wake (two wake halves)
       DI2    = DI2   *2.0
       DI2_S2 = DI2_S2*2.0
       DI2_U2 = DI2_U2*2.0
       DI2_T2 = DI2_T2*2.0
       DI2_D2 = DI2_D2*2.0
       DI2_MS = DI2_MS*2.0
       DI2_RE = DI2_RE*2.0
      ENDIF
C
C---- BL thickness (Delta) from simplified Green's correlation
      DE2     = (3.15 + 1.72/(HK2-1.0)   )*T2  +  D2
      DE2_HK2 = (     - 1.72/(HK2-1.0)**2)*T2
C
      DE2_U2 = DE2_HK2*HK2_U2
      DE2_T2 = DE2_HK2*HK2_T2 + (3.15 + 1.72/(HK2-1.0))
      DE2_D2 = DE2_HK2*HK2_D2 + 1.0
      DE2_MS = DE2_HK2*HK2_MS
C
ccc      HDMAX = 15.0
      HDMAX = 12.0
      IF(DE2 .GT. HDMAX*T2) THEN
cccc      IF(DE2 .GT. HDMAX*T2 .AND. (HK2 .GT. 4.0 .OR. ITYP.EQ.3)) THEN
       DE2    = HDMAX*T2
       DE2_U2 =  0.0
       DE2_T2 = HDMAX
       DE2_D2 =  0.0
       DE2_MS =  0.0
      ENDIF
C
      RETURN
      END
 

      SUBROUTINE BLMID(ITYP)
C----------------------------------------------------
C     Calculates midpoint skin friction CFM
C
C      ITYP = 1 :  laminar
C      ITYP = 2 :  turbulent
C      ITYP = 3 :  turbulent wake
C----------------------------------------------------
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
C
C---- set similarity variables if not defined
      IF(SIMI) THEN
       HK1    = HK2
       HK1_T1 = HK2_T2
       HK1_D1 = HK2_D2
       HK1_U1 = HK2_U2
       HK1_MS = HK2_MS
       RT1    = RT2
       RT1_T1 = RT2_T2
       RT1_U1 = RT2_U2
       RT1_MS = RT2_MS
       RT1_RE = RT2_RE
       M1    = M2
       M1_U1 = M2_U2
       M1_MS = M2_MS
      ENDIF
C
C---- define stuff for midpoint CF
      HKA = 0.5*(HK1 + HK2)
      RTA = 0.5*(RT1 + RT2)
      MA  = 0.5*(M1  + M2 )
C
C---- midpoint skin friction coefficient  (zero in wake)
      IF(ITYP.EQ.3) THEN
       CFM     = 0.
       CFM_HKA = 0.
       CFM_RTA = 0.
       CFM_MA  = 0.
       CFM_MS  = 0.
      ELSE IF(ITYP.EQ.1) THEN
       CALL CFL( HKA, RTA, MA, CFM, CFM_HKA, CFM_RTA, CFM_MA )
      ELSE
       CALL CFT( HKA, RTA, MA, CFM, CFM_HKA, CFM_RTA, CFM_MA )
       CALL CFL( HKA, RTA, MA, CFML,CFML_HKA,CFML_RTA,CFML_MA)
       IF(CFML.GT.CFM) THEN
ccc      write(*,*) 'Cft Cfl Rt Hk:', CFM, CFML, RTA, HKA, 0.5*(X1+X2)
         CFM     = CFML
         CFM_HKA = CFML_HKA
         CFM_RTA = CFML_RTA
         CFM_MA  = CFML_MA
       ENDIF
      ENDIF
C
      CFM_U1 = 0.5*(CFM_HKA*HK1_U1 + CFM_MA*M1_U1 + CFM_RTA*RT1_U1)
      CFM_T1 = 0.5*(CFM_HKA*HK1_T1 +                CFM_RTA*RT1_T1)
      CFM_D1 = 0.5*(CFM_HKA*HK1_D1                                )
C
      CFM_U2 = 0.5*(CFM_HKA*HK2_U2 + CFM_MA*M2_U2 + CFM_RTA*RT2_U2)
      CFM_T2 = 0.5*(CFM_HKA*HK2_T2 +                CFM_RTA*RT2_T2)
      CFM_D2 = 0.5*(CFM_HKA*HK2_D2                                )
C
      CFM_MS = 0.5*(CFM_HKA*HK1_MS + CFM_MA*M1_MS + CFM_RTA*RT1_MS
     &            + CFM_HKA*HK2_MS + CFM_MA*M2_MS + CFM_RTA*RT2_MS)
      CFM_RE = 0.5*(                                CFM_RTA*RT1_RE
     &                                            + CFM_RTA*RT2_RE)
C
      RETURN
      END ! BLMID

 
      SUBROUTINE TRDIF
C-----------------------------------------------
C     Sets up the Newton system governing the
C     transition interval.  Equations governing
C     the  laminar  part  X1 < xi < XT  and
C     the turbulent part  XT < xi < X2
C     are simply summed.
C-----------------------------------------------
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
      REAL  BL1(4,5), BL2(4,5), BLREZ(4), BLM(4), BLR(4), BLX(4)
     &    , BT1(4,5), BT2(4,5), BTREZ(4), BTM(4), BTR(4), BTX(4)
C
C---- save variables and sensitivities for future restoration
      DO 5 ICOM=1, NCOM
        C1SAV(ICOM) = COM1(ICOM)
        C2SAV(ICOM) = COM2(ICOM)
    5 CONTINUE
C
C---- weighting factors for linear interpolation to transition point
      WF2    = (XT-X1)/(X2-X1)
      WF2_XT = 1.0/(X2-X1)
C
      WF2_A1 = WF2_XT*XT_A1
      WF2_X1 = WF2_XT*XT_X1 + (WF2-1.0)/(X2-X1)
      WF2_X2 = WF2_XT*XT_X2 -  WF2     /(X2-X1)
      WF2_T1 = WF2_XT*XT_T1
      WF2_T2 = WF2_XT*XT_T2
      WF2_D1 = WF2_XT*XT_D1
      WF2_D2 = WF2_XT*XT_D2
      WF2_U1 = WF2_XT*XT_U1
      WF2_U2 = WF2_XT*XT_U2
      WF2_MS = WF2_XT*XT_MS
      WF2_RE = WF2_XT*XT_RE
      WF2_XF = WF2_XT*XT_XF
C
      WF1    = 1.0 - WF2
      WF1_A1 = -WF2_A1
      WF1_X1 = -WF2_X1
      WF1_X2 = -WF2_X2
      WF1_T1 = -WF2_T1
      WF1_T2 = -WF2_T2
      WF1_D1 = -WF2_D1
      WF1_D2 = -WF2_D2
      WF1_U1 = -WF2_U1
      WF1_U2 = -WF2_U2
      WF1_MS = -WF2_MS
      WF1_RE = -WF2_RE
      WF1_XF = -WF2_XF
C
C
C**** FIRST,  do laminar part between X1 and XT
C
C-----interpolate primary variables to transition point
      TT    = T1*WF1    + T2*WF2
      TT_A1 = T1*WF1_A1 + T2*WF2_A1
      TT_X1 = T1*WF1_X1 + T2*WF2_X1
      TT_X2 = T1*WF1_X2 + T2*WF2_X2
      TT_T1 = T1*WF1_T1 + T2*WF2_T1 + WF1
      TT_T2 = T1*WF1_T2 + T2*WF2_T2 + WF2
      TT_D1 = T1*WF1_D1 + T2*WF2_D1
      TT_D2 = T1*WF1_D2 + T2*WF2_D2
      TT_U1 = T1*WF1_U1 + T2*WF2_U1
      TT_U2 = T1*WF1_U2 + T2*WF2_U2
      TT_MS = T1*WF1_MS + T2*WF2_MS
      TT_RE = T1*WF1_RE + T2*WF2_RE
      TT_XF = T1*WF1_XF + T2*WF2_XF
C
      DT    = D1*WF1    + D2*WF2
      DT_A1 = D1*WF1_A1 + D2*WF2_A1
      DT_X1 = D1*WF1_X1 + D2*WF2_X1
      DT_X2 = D1*WF1_X2 + D2*WF2_X2
      DT_T1 = D1*WF1_T1 + D2*WF2_T1
      DT_T2 = D1*WF1_T2 + D2*WF2_T2
      DT_D1 = D1*WF1_D1 + D2*WF2_D1 + WF1
      DT_D2 = D1*WF1_D2 + D2*WF2_D2 + WF2
      DT_U1 = D1*WF1_U1 + D2*WF2_U1
      DT_U2 = D1*WF1_U2 + D2*WF2_U2
      DT_MS = D1*WF1_MS + D2*WF2_MS
      DT_RE = D1*WF1_RE + D2*WF2_RE
      DT_XF = D1*WF1_XF + D2*WF2_XF
C
      UT    = U1*WF1    + U2*WF2
      UT_A1 = U1*WF1_A1 + U2*WF2_A1
      UT_X1 = U1*WF1_X1 + U2*WF2_X1
      UT_X2 = U1*WF1_X2 + U2*WF2_X2
      UT_T1 = U1*WF1_T1 + U2*WF2_T1
      UT_T2 = U1*WF1_T2 + U2*WF2_T2
      UT_D1 = U1*WF1_D1 + U2*WF2_D1
      UT_D2 = U1*WF1_D2 + U2*WF2_D2
      UT_U1 = U1*WF1_U1 + U2*WF2_U1 + WF1
      UT_U2 = U1*WF1_U2 + U2*WF2_U2 + WF2
      UT_MS = U1*WF1_MS + U2*WF2_MS
      UT_RE = U1*WF1_RE + U2*WF2_RE
      UT_XF = U1*WF1_XF + U2*WF2_XF
C
C---- set primary "T" variables at XT  (really placed into "2" variables)
      X2 = XT
      T2 = TT
      D2 = DT
      U2 = UT
C
      AMPL2 = AMCRIT
      S2 = 0.
C
C---- calculate laminar secondary "T" variables
      CALL BLKIN
      CALL BLVAR(1)
C
C---- calculate X1-XT midpoint CFM value
      CALL BLMID(1)
C=
C=    at this point, all "2" variables are really "T" variables at XT
C=
C
C---- set up Newton system for dAm, dTh, dDs, dUe, dXi  at  X1 and XT
      CALL BLDIF(1)
C
C---- The current Newton system is in terms of "1" and "T" variables,
C-    so calculate its equivalent in terms of "1" and "2" variables.
C-    In other words, convert residual sensitivities wrt "T" variables
C-    into sensitivities wrt "1" and "2" variables.  The amplification
C-    equation is unnecessary here, so the K=1 row is left empty.
      DO 10 K=2, 3
        BLREZ(K) = VSREZ(K)
        BLM(K)   = VSM(K)
     &           + VS2(K,2)*TT_MS
     &           + VS2(K,3)*DT_MS
     &           + VS2(K,4)*UT_MS
     &           + VS2(K,5)*XT_MS
        BLR(K)   = VSR(K)
     &           + VS2(K,2)*TT_RE
     &           + VS2(K,3)*DT_RE
     &           + VS2(K,4)*UT_RE
     &           + VS2(K,5)*XT_RE
        BLX(K)   = VSX(K)
     &           + VS2(K,2)*TT_XF
     &           + VS2(K,3)*DT_XF
     &           + VS2(K,4)*UT_XF
     &           + VS2(K,5)*XT_XF
C
        BL1(K,1) = VS1(K,1)
     &           + VS2(K,2)*TT_A1
     &           + VS2(K,3)*DT_A1
     &           + VS2(K,4)*UT_A1
     &           + VS2(K,5)*XT_A1
        BL1(K,2) = VS1(K,2)
     &           + VS2(K,2)*TT_T1
     &           + VS2(K,3)*DT_T1
     &           + VS2(K,4)*UT_T1
     &           + VS2(K,5)*XT_T1
        BL1(K,3) = VS1(K,3)
     &           + VS2(K,2)*TT_D1
     &           + VS2(K,3)*DT_D1
     &           + VS2(K,4)*UT_D1
     &           + VS2(K,5)*XT_D1
        BL1(K,4) = VS1(K,4)
     &           + VS2(K,2)*TT_U1
     &           + VS2(K,3)*DT_U1
     &           + VS2(K,4)*UT_U1
     &           + VS2(K,5)*XT_U1
        BL1(K,5) = VS1(K,5)
     &           + VS2(K,2)*TT_X1
     &           + VS2(K,3)*DT_X1
     &           + VS2(K,4)*UT_X1
     &           + VS2(K,5)*XT_X1
C
        BL2(K,1) = 0.
        BL2(K,2) = VS2(K,2)*TT_T2
     &           + VS2(K,3)*DT_T2
     &           + VS2(K,4)*UT_T2
     &           + VS2(K,5)*XT_T2
        BL2(K,3) = VS2(K,2)*TT_D2
     &           + VS2(K,3)*DT_D2
     &           + VS2(K,4)*UT_D2
     &           + VS2(K,5)*XT_D2
        BL2(K,4) = VS2(K,2)*TT_U2
     &           + VS2(K,3)*DT_U2
     &           + VS2(K,4)*UT_U2
     &           + VS2(K,5)*XT_U2
        BL2(K,5) = VS2(K,2)*TT_X2
     &           + VS2(K,3)*DT_X2
     &           + VS2(K,4)*UT_X2
     &           + VS2(K,5)*XT_X2
C
   10 CONTINUE
C
C
C**** SECOND, set up turbulent part between XT and X2  ****
C
C---- calculate equilibrium shear coefficient CQT at transition point
      CALL BLVAR(2)
C
C---- set initial shear coefficient value ST at transition point
C-    ( note that CQ2, CQ2_T2, etc. are really "CQT", "CQT_TT", etc.)
C
      CTR     = CTRCON*EXP(-CTRCEX/(HK2-1.0))
      CTR_HK2 = CTR * CTRCEX/(HK2-1.0)**2
C
c      CTR     = 1.1*EXP(-10.0/HK2**2)
c      CTR_HK2 = CTR * 10.0 * 2.0/HK2**3
C
CCC      CTR = 1.2
CCC      CTR = 0.7
CCC      CTR_HK2 = 0.0
C
      ST    = CTR*CQ2
      ST_TT = CTR*CQ2_T2 + CQ2*CTR_HK2*HK2_T2
      ST_DT = CTR*CQ2_D2 + CQ2*CTR_HK2*HK2_D2
      ST_UT = CTR*CQ2_U2 + CQ2*CTR_HK2*HK2_U2
      ST_MS = CTR*CQ2_MS + CQ2*CTR_HK2*HK2_MS
      ST_RE = CTR*CQ2_RE
C
C---- calculate ST sensitivities wrt the actual "1" and "2" variables
      ST_A1 = ST_TT*TT_A1 + ST_DT*DT_A1 + ST_UT*UT_A1
      ST_X1 = ST_TT*TT_X1 + ST_DT*DT_X1 + ST_UT*UT_X1
      ST_X2 = ST_TT*TT_X2 + ST_DT*DT_X2 + ST_UT*UT_X2
      ST_T1 = ST_TT*TT_T1 + ST_DT*DT_T1 + ST_UT*UT_T1
      ST_T2 = ST_TT*TT_T2 + ST_DT*DT_T2 + ST_UT*UT_T2
      ST_D1 = ST_TT*TT_D1 + ST_DT*DT_D1 + ST_UT*UT_D1
      ST_D2 = ST_TT*TT_D2 + ST_DT*DT_D2 + ST_UT*UT_D2
      ST_U1 = ST_TT*TT_U1 + ST_DT*DT_U1 + ST_UT*UT_U1
      ST_U2 = ST_TT*TT_U2 + ST_DT*DT_U2 + ST_UT*UT_U2
      ST_MS = ST_TT*TT_MS + ST_DT*DT_MS + ST_UT*UT_MS + ST_MS
      ST_RE = ST_TT*TT_RE + ST_DT*DT_RE + ST_UT*UT_RE + ST_RE
      ST_XF = ST_TT*TT_XF + ST_DT*DT_XF + ST_UT*UT_XF
C
      AMPL2 = 0.
      S2 = ST
C
C---- recalculate turbulent secondary "T" variables using proper CTI
      CALL BLVAR(2)
C
C---- set "1" variables to "T" variables and reset "2" variables
C-    to their saved turbulent values
      DO 30 ICOM=1, NCOM
        COM1(ICOM) = COM2(ICOM)
        COM2(ICOM) = C2SAV(ICOM)
   30 CONTINUE
C
C---- calculate XT-X2 midpoint CFM value
      CALL BLMID(2)
C
C---- set up Newton system for dCt, dTh, dDs, dUe, dXi  at  XT and X2
      CALL BLDIF(2)
C
C---- convert sensitivities wrt "T" variables into sensitivities
C-    wrt "1" and "2" variables as done before for the laminar part
      DO 40 K=1, 3
        BTREZ(K) = VSREZ(K)
        BTM(K)   = VSM(K) 
     &           + VS1(K,1)*ST_MS
     &           + VS1(K,2)*TT_MS
     &           + VS1(K,3)*DT_MS
     &           + VS1(K,4)*UT_MS
     &           + VS1(K,5)*XT_MS
        BTR(K)   = VSR(K) 
     &           + VS1(K,1)*ST_RE
     &           + VS1(K,2)*TT_RE
     &           + VS1(K,3)*DT_RE
     &           + VS1(K,4)*UT_RE
     &           + VS1(K,5)*XT_RE
        BTX(K)   = VSX(K)
     &           + VS1(K,1)*ST_XF
     &           + VS1(K,2)*TT_XF
     &           + VS1(K,3)*DT_XF
     &           + VS1(K,4)*UT_XF
     &           + VS1(K,5)*XT_XF
C
        BT1(K,1) = VS1(K,1)*ST_A1
     &           + VS1(K,2)*TT_A1
     &           + VS1(K,3)*DT_A1
     &           + VS1(K,4)*UT_A1
     &           + VS1(K,5)*XT_A1
        BT1(K,2) = VS1(K,1)*ST_T1
     &           + VS1(K,2)*TT_T1
     &           + VS1(K,3)*DT_T1
     &           + VS1(K,4)*UT_T1
     &           + VS1(K,5)*XT_T1
        BT1(K,3) = VS1(K,1)*ST_D1
     &           + VS1(K,2)*TT_D1
     &           + VS1(K,3)*DT_D1
     &           + VS1(K,4)*UT_D1
     &           + VS1(K,5)*XT_D1
        BT1(K,4) = VS1(K,1)*ST_U1
     &           + VS1(K,2)*TT_U1
     &           + VS1(K,3)*DT_U1
     &           + VS1(K,4)*UT_U1
     &           + VS1(K,5)*XT_U1
        BT1(K,5) = VS1(K,1)*ST_X1
     &           + VS1(K,2)*TT_X1
     &           + VS1(K,3)*DT_X1
     &           + VS1(K,4)*UT_X1
     &           + VS1(K,5)*XT_X1
C
        BT2(K,1) = VS2(K,1)
        BT2(K,2) = VS2(K,2)
     &           + VS1(K,1)*ST_T2
     &           + VS1(K,2)*TT_T2
     &           + VS1(K,3)*DT_T2
     &           + VS1(K,4)*UT_T2
     &           + VS1(K,5)*XT_T2
        BT2(K,3) = VS2(K,3)
     &           + VS1(K,1)*ST_D2
     &           + VS1(K,2)*TT_D2
     &           + VS1(K,3)*DT_D2
     &           + VS1(K,4)*UT_D2
     &           + VS1(K,5)*XT_D2
        BT2(K,4) = VS2(K,4)
     &           + VS1(K,1)*ST_U2
     &           + VS1(K,2)*TT_U2
     &           + VS1(K,3)*DT_U2
     &           + VS1(K,4)*UT_U2
     &           + VS1(K,5)*XT_U2
        BT2(K,5) = VS2(K,5)
     &           + VS1(K,1)*ST_X2
     &           + VS1(K,2)*TT_X2
     &           + VS1(K,3)*DT_X2
     &           + VS1(K,4)*UT_X2
     &           + VS1(K,5)*XT_X2
C
   40 CONTINUE
C
C---- Add up laminar and turbulent parts to get final system
C-    in terms of honest-to-God "1" and "2" variables.
      VSREZ(1) =            BTREZ(1)
      VSREZ(2) = BLREZ(2) + BTREZ(2)
      VSREZ(3) = BLREZ(3) + BTREZ(3)
      VSM(1)   =            BTM(1)
      VSM(2)   = BLM(2)   + BTM(2)
      VSM(3)   = BLM(3)   + BTM(3)
      VSR(1)   =            BTR(1)
      VSR(2)   = BLR(2)   + BTR(2)
      VSR(3)   = BLR(3)   + BTR(3)
      VSX(1)   =            BTX(1)
      VSX(2)   = BLX(2)   + BTX(2)
      VSX(3)   = BLX(3)   + BTX(3)
      DO 60 L=1, 5
        VS1(1,L) =            BT1(1,L)
        VS2(1,L) =            BT2(1,L)
        VS1(2,L) = BL1(2,L) + BT1(2,L)
        VS2(2,L) = BL2(2,L) + BT2(2,L)
        VS1(3,L) = BL1(3,L) + BT1(3,L)
        VS2(3,L) = BL2(3,L) + BT2(3,L)
   60 CONTINUE
C
C---- To be sanitary, restore "1" quantities which got clobbered
C-    in all of the numerical gymnastics above.  The "2" variables
C-    were already restored for the XT-X2 differencing part.
      DO 70 ICOM=1, NCOM
        COM1(ICOM) = C1SAV(ICOM)
   70 CONTINUE
C
      RETURN
      END
 
 
      SUBROUTINE BLDIF(ITYP)
C-----------------------------------------------------------
C     Sets up the Newton system coefficients and residuals
C
C        ITYP = 0 :  similarity station
C        ITYP = 1 :  laminar interval
C        ITYP = 2 :  turbulent interval
C        ITYP = 3 :  wake interval
C
C     This routine knows nothing about a transition interval,
C     which is taken care of by TRDIF.
C-----------------------------------------------------------
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
C
      IF(ITYP.EQ.0) THEN
C----- similarity logarithmic differences  (prescribed)
       XLOG = 1.0
       ULOG = BULE
       TLOG = 0.5*(1.0 - BULE)
       HLOG = 0.
       DDLOG = 0.
      ELSE
C----- usual logarithmic differences
       XLOG = LOG(X2/X1)
       ULOG = LOG(U2/U1)
       TLOG = LOG(T2/T1)
       HLOG = LOG(HS2/HS1)
C       XLOG = 2.0*(X2-X1)/(X2+X1)
C       ULOG = 2.0*(U2-U1)/(U2+U1)
C       TLOG = 2.0*(T2-T1)/(T2+T1)
C       HLOG = 2.0*(HS2-HS1)/(HS2+HS1)
       DDLOG = 1.0
      ENDIF
C
      DO 55 K=1, 4
        VSREZ(K) = 0.
        VSM(K) = 0.
        VSR(K) = 0.
        VSX(K) = 0.
        DO 551 L=1, 5
          VS1(K,L) = 0.
          VS2(K,L) = 0.
  551   CONTINUE
   55 CONTINUE
C
C---- set triggering constant for local upwinding
      HUPWT = 1.0
C
ccc      HDCON = 5.0*HUPWT
ccc      HD_HK1 = 0.0
ccc      HD_HK2 = 0.0
C
      HDCON  =  5.0*HUPWT/HK2**2
      HD_HK1 =  0.0
      HD_HK2 = -HDCON*2.0/HK2
C
C---- use less upwinding in the wake
      IF(ITYP.EQ.3) THEN
       HDCON  =  HUPWT/HK2**2
       HD_HK1 =  0.0
       HD_HK2 = -HDCON*2.0/HK2
      ENDIF
C
C---- local upwinding is based on local change in  log(Hk-1)
C-    (mainly kicks in at transition)
      ARG = ABS((HK2-1.0)/(HK1-1.0))
      HL = LOG(ARG)
      HL_HK1 = -1.0/(HK1-1.0)
      HL_HK2 =  1.0/(HK2-1.0)
C
C---- set local upwinding parameter UPW and linearize it
C
C       UPW = 0.5   Trapezoidal
C       UPW = 1.0   Backward Euler
C
      HLSQ = MIN( HL**2 , 15.0 )
      EHH = EXP(-HLSQ*HDCON)
      UPW = 1.0 - 0.5*EHH
      UPW_HL =        EHH * HL  *HDCON
      UPW_HD =    0.5*EHH * HLSQ
C
      UPW_HK1 = UPW_HL*HL_HK1 + UPW_HD*HD_HK1
      UPW_HK2 = UPW_HL*HL_HK2 + UPW_HD*HD_HK2
C
      UPW_U1 = UPW_HK1*HK1_U1
      UPW_T1 = UPW_HK1*HK1_T1
      UPW_D1 = UPW_HK1*HK1_D1
      UPW_U2 = UPW_HK2*HK2_U2
      UPW_T2 = UPW_HK2*HK2_T2
      UPW_D2 = UPW_HK2*HK2_D2
      UPW_MS = UPW_HK1*HK1_MS
     &       + UPW_HK2*HK2_MS
C
C
      IF(ITYP.EQ.0) THEN
C
C***** LE point -->  set zero amplification factor
       VS2(1,1) = 1.0
       VSR(1)   = 0.
       VSREZ(1) = -AMPL2
C
      ELSE IF(ITYP.EQ.1) THEN
C
C***** laminar part -->  set amplification equation
C
C----- set average amplification AX over interval X1..X2
       CALL AXSET( HK1,    T1,    RT1, AMPL1,  
     &             HK2,    T2,    RT2, AMPL2, AMCRIT, IDAMPV,
     &      AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &          AX_HK2, AX_T2, AX_RT2, AX_A2 )
C
       REZC = AMPL2 - AMPL1 - AX*(X2-X1)
       Z_AX = -(X2-X1)
C
       VS1(1,1) = Z_AX* AX_A1  -  1.0
       VS1(1,2) = Z_AX*(AX_HK1*HK1_T1 + AX_T1 + AX_RT1*RT1_T1)
       VS1(1,3) = Z_AX*(AX_HK1*HK1_D1                        )
       VS1(1,4) = Z_AX*(AX_HK1*HK1_U1         + AX_RT1*RT1_U1)
       VS1(1,5) =  AX
       VS2(1,1) = Z_AX* AX_A2  +  1.0
       VS2(1,2) = Z_AX*(AX_HK2*HK2_T2 + AX_T2 + AX_RT2*RT2_T2)
       VS2(1,3) = Z_AX*(AX_HK2*HK2_D2                        )         
       VS2(1,4) = Z_AX*(AX_HK2*HK2_U2         + AX_RT2*RT2_U2)
       VS2(1,5) = -AX
       VSM(1)   = Z_AX*(AX_HK1*HK1_MS         + AX_RT1*RT1_MS
     &                + AX_HK2*HK2_MS         + AX_RT2*RT2_MS)
       VSR(1)   = Z_AX*(                        AX_RT1*RT1_RE
     &                                        + AX_RT2*RT2_RE)
       VSX(1)   = 0.
       VSREZ(1) = -REZC
C
      ELSE
C
C***** turbulent part -->  set shear lag equation
C
       SA  = (1.0-UPW)*S1  + UPW*S2
       CQA = (1.0-UPW)*CQ1 + UPW*CQ2
       CFA = (1.0-UPW)*CF1 + UPW*CF2
       HKA = (1.0-UPW)*HK1 + UPW*HK2
C
       USA = 0.5*(US1 + US2)
       RTA = 0.5*(RT1 + RT2)
       DEA = 0.5*(DE1 + DE2)
       DA  = 0.5*(D1  + D2 )
C
C
       IF(ITYP.EQ.3) THEN
C------ increased dissipation length in wake (decrease its reciprocal)
        ALD = DLCON
       ELSE
        ALD = 1.0
       ENDIF
C
C----- set and linearize  equilibrium 1/Ue dUe/dx   ...  NEW  12 Oct 94
       IF(ITYP.EQ.2) THEN
        GCC = GCCON
        HKC     = HKA - 1.0 - GCC/RTA
        HKC_HKA = 1.0
        HKC_RTA =             GCC/RTA**2
        IF(HKC .LT. 0.01) THEN
         HKC = 0.01
         HKC_HKA = 0.0
         HKC_RTA = 0.0
        ENDIF
       ELSE
        GCC = 0.0
        HKC = HKA - 1.0
        HKC_HKA = 1.0
        HKC_RTA = 0.0
       ENDIF
C
       HR     = HKC     / (GACON*ALD*HKA)
       HR_HKA = HKC_HKA / (GACON*ALD*HKA) - HR / HKA
       HR_RTA = HKC_RTA / (GACON*ALD*HKA)
C
       UQ     = (0.5*CFA - HR**2) / (GBCON*DA)
       UQ_HKA =   -2.0*HR*HR_HKA  / (GBCON*DA)
       UQ_RTA =   -2.0*HR*HR_RTA  / (GBCON*DA)
       UQ_CFA =   0.5             / (GBCON*DA)
       UQ_DA  = -UQ/DA
       UQ_UPW = UQ_CFA*(CF2-CF1) + UQ_HKA*(HK2-HK1)
C
       UQ_T1 = (1.0-UPW)*(UQ_CFA*CF1_T1 + UQ_HKA*HK1_T1) + UQ_UPW*UPW_T1
       UQ_D1 = (1.0-UPW)*(UQ_CFA*CF1_D1 + UQ_HKA*HK1_D1) + UQ_UPW*UPW_D1
       UQ_U1 = (1.0-UPW)*(UQ_CFA*CF1_U1 + UQ_HKA*HK1_U1) + UQ_UPW*UPW_U1
       UQ_T2 =      UPW *(UQ_CFA*CF2_T2 + UQ_HKA*HK2_T2) + UQ_UPW*UPW_T2
       UQ_D2 =      UPW *(UQ_CFA*CF2_D2 + UQ_HKA*HK2_D2) + UQ_UPW*UPW_D2
       UQ_U2 =      UPW *(UQ_CFA*CF2_U2 + UQ_HKA*HK2_U2) + UQ_UPW*UPW_U2
       UQ_MS = (1.0-UPW)*(UQ_CFA*CF1_MS + UQ_HKA*HK1_MS) + UQ_UPW*UPW_MS
     &       +      UPW *(UQ_CFA*CF2_MS + UQ_HKA*HK2_MS)
       UQ_RE = (1.0-UPW)* UQ_CFA*CF1_RE
     &       +      UPW * UQ_CFA*CF2_RE
C
       UQ_T1 = UQ_T1             + 0.5*UQ_RTA*RT1_T1
       UQ_D1 = UQ_D1 + 0.5*UQ_DA
       UQ_U1 = UQ_U1             + 0.5*UQ_RTA*RT1_U1
       UQ_T2 = UQ_T2             + 0.5*UQ_RTA*RT2_T2
       UQ_D2 = UQ_D2 + 0.5*UQ_DA
       UQ_U2 = UQ_U2             + 0.5*UQ_RTA*RT2_U2
       UQ_MS = UQ_MS             + 0.5*UQ_RTA*RT1_MS
     &                           + 0.5*UQ_RTA*RT2_MS
       UQ_RE = UQ_RE             + 0.5*UQ_RTA*RT1_RE
     &                           + 0.5*UQ_RTA*RT2_RE
C
       SCC = SCCON*1.333/(1.0+USA)
       SCC_USA = -SCC/(1.0+USA)
C
       SCC_US1 = SCC_USA*0.5
       SCC_US2 = SCC_USA*0.5
C
C
       SLOG = LOG(S2/S1)
       DXI = X2 - X1
C
       REZC = SCC*(CQA - SA*ALD)*DXI
     &      - DEA*2.0*          SLOG
     &      + DEA*2.0*(UQ*DXI - ULOG)*DUXCON
C

c        if(  ! (rt2.gt.1.0e3 .and. rt1.le.1.0e3) .or.
c     &     (rt2.gt.1.0e4 .and. rt1.le.1.0e4) .or.
c     &     (rt2.gt.1.0e5 .and. rt1.le.1.0e5)        ) then
c           gga = (HKA-1.0-GCC/RTA)/HKA / sqrt(0.5*CFA)
c           write(*,4455) rta, hka, gga, cfa, cqa, sa, uq, ulog/dxi
c 4455      format(1x,f7.0, 2f9.4,f10.6,2f8.5,2f10.5)
c        endif


       Z_CFA = DEA*2.0*UQ_CFA*DXI * DUXCON
       Z_HKA = DEA*2.0*UQ_HKA*DXI * DUXCON
       Z_DA  = DEA*2.0*UQ_DA *DXI * DUXCON
       Z_SL = -DEA*2.0
       Z_UL = -DEA*2.0 * DUXCON
       Z_DXI = SCC    *(CQA - SA*ALD)     + DEA*2.0*UQ*DUXCON
       Z_USA = SCC_USA*(CQA - SA*ALD)*DXI
       Z_CQA = SCC*DXI
       Z_SA = -SCC*DXI*ALD
       Z_DEA = 2.0*((UQ*DXI - ULOG)*DUXCON - SLOG)
C
       Z_UPW = Z_CQA*(CQ2-CQ1) + Z_SA *(S2 -S1 )
     &       + Z_CFA*(CF2-CF1) + Z_HKA*(HK2-HK1)
       Z_DE1 = 0.5*Z_DEA
       Z_DE2 = 0.5*Z_DEA
       Z_US1 = 0.5*Z_USA
       Z_US2 = 0.5*Z_USA
       Z_D1  = 0.5*Z_DA
       Z_D2  = 0.5*Z_DA
       Z_U1  =                 - Z_UL/U1
       Z_U2  =                   Z_UL/U2
       Z_X1  = -Z_DXI
       Z_X2  =  Z_DXI
       Z_S1  = (1.0-UPW)*Z_SA  - Z_SL/S1
       Z_S2  =      UPW *Z_SA  + Z_SL/S2
       Z_CQ1 = (1.0-UPW)*Z_CQA
       Z_CQ2 =      UPW *Z_CQA
       Z_CF1 = (1.0-UPW)*Z_CFA
       Z_CF2 =      UPW *Z_CFA
       Z_HK1 = (1.0-UPW)*Z_HKA
       Z_HK2 =      UPW *Z_HKA
C
       VS1(1,1) = Z_S1
       VS1(1,2) =        Z_UPW*UPW_T1 + Z_DE1*DE1_T1 + Z_US1*US1_T1
       VS1(1,3) = Z_D1 + Z_UPW*UPW_D1 + Z_DE1*DE1_D1 + Z_US1*US1_D1
       VS1(1,4) = Z_U1 + Z_UPW*UPW_U1 + Z_DE1*DE1_U1 + Z_US1*US1_U1
       VS1(1,5) = Z_X1
       VS2(1,1) = Z_S2
       VS2(1,2) =        Z_UPW*UPW_T2 + Z_DE2*DE2_T2 + Z_US2*US2_T2
       VS2(1,3) = Z_D2 + Z_UPW*UPW_D2 + Z_DE2*DE2_D2 + Z_US2*US2_D2
       VS2(1,4) = Z_U2 + Z_UPW*UPW_U2 + Z_DE2*DE2_U2 + Z_US2*US2_U2
       VS2(1,5) = Z_X2
       VSM(1)   =        Z_UPW*UPW_MS + Z_DE1*DE1_MS + Z_US1*US1_MS
     &                                + Z_DE2*DE2_MS + Z_US2*US2_MS
C
       VS1(1,2) = VS1(1,2) + Z_CQ1*CQ1_T1 + Z_CF1*CF1_T1 + Z_HK1*HK1_T1
       VS1(1,3) = VS1(1,3) + Z_CQ1*CQ1_D1 + Z_CF1*CF1_D1 + Z_HK1*HK1_D1
       VS1(1,4) = VS1(1,4) + Z_CQ1*CQ1_U1 + Z_CF1*CF1_U1 + Z_HK1*HK1_U1
C
       VS2(1,2) = VS2(1,2) + Z_CQ2*CQ2_T2 + Z_CF2*CF2_T2 + Z_HK2*HK2_T2
       VS2(1,3) = VS2(1,3) + Z_CQ2*CQ2_D2 + Z_CF2*CF2_D2 + Z_HK2*HK2_D2
       VS2(1,4) = VS2(1,4) + Z_CQ2*CQ2_U2 + Z_CF2*CF2_U2 + Z_HK2*HK2_U2
C
       VSM(1)   = VSM(1)   + Z_CQ1*CQ1_MS + Z_CF1*CF1_MS + Z_HK1*HK1_MS
     &                     + Z_CQ2*CQ2_MS + Z_CF2*CF2_MS + Z_HK2*HK2_MS
       VSR(1)   =            Z_CQ1*CQ1_RE + Z_CF1*CF1_RE
     &                     + Z_CQ2*CQ2_RE + Z_CF2*CF2_RE
       VSX(1)   = 0.
       VSREZ(1) = -REZC
C
      ENDIF
C
C**** Set up momentum equation
      HA = 0.5*(H1 + H2)
      MA = 0.5*(M1 + M2)
      XA = 0.5*(X1 + X2)
      TA = 0.5*(T1 + T2)
      HWA = 0.5*(DW1/T1 + DW2/T2)
C
C---- set Cf term, using central value CFM for better accuracy in drag
      CFX     = 0.50*CFM*XA/TA  +  0.25*(CF1*X1/T1 + CF2*X2/T2)
      CFX_XA  = 0.50*CFM   /TA
      CFX_TA  = -.50*CFM*XA/TA**2
C
      CFX_X1  = 0.25*CF1   /T1     + CFX_XA*0.5
      CFX_X2  = 0.25*CF2   /T2     + CFX_XA*0.5
      CFX_T1  = -.25*CF1*X1/T1**2  + CFX_TA*0.5
      CFX_T2  = -.25*CF2*X2/T2**2  + CFX_TA*0.5
      CFX_CF1 = 0.25*    X1/T1
      CFX_CF2 = 0.25*    X2/T2
      CFX_CFM = 0.50*    XA/TA
C
      BTMP = HA + 2.0 - MA + HWA
C
      REZT  = TLOG + BTMP*ULOG - XLOG*0.5*CFX
      Z_CFX = -XLOG*0.5
      Z_HA  =  ULOG
      Z_HWA =  ULOG
      Z_MA  = -ULOG
      Z_XL  =-DDLOG * 0.5*CFX
      Z_UL  = DDLOG * BTMP
      Z_TL  = DDLOG
C
      Z_CFM = Z_CFX*CFX_CFM
      Z_CF1 = Z_CFX*CFX_CF1
      Z_CF2 = Z_CFX*CFX_CF2
C
      Z_T1 = -Z_TL/T1 + Z_CFX*CFX_T1 + Z_HWA*0.5*(-DW1/T1**2)
      Z_T2 =  Z_TL/T2 + Z_CFX*CFX_T2 + Z_HWA*0.5*(-DW2/T2**2)
      Z_X1 = -Z_XL/X1 + Z_CFX*CFX_X1
      Z_X2 =  Z_XL/X2 + Z_CFX*CFX_X2
      Z_U1 = -Z_UL/U1
      Z_U2 =  Z_UL/U2
C
      VS1(2,2) = 0.5*Z_HA*H1_T1 + Z_CFM*CFM_T1 + Z_CF1*CF1_T1 + Z_T1
      VS1(2,3) = 0.5*Z_HA*H1_D1 + Z_CFM*CFM_D1 + Z_CF1*CF1_D1
      VS1(2,4) = 0.5*Z_MA*M1_U1 + Z_CFM*CFM_U1 + Z_CF1*CF1_U1 + Z_U1
      VS1(2,5) =                                                Z_X1
      VS2(2,2) = 0.5*Z_HA*H2_T2 + Z_CFM*CFM_T2 + Z_CF2*CF2_T2 + Z_T2
      VS2(2,3) = 0.5*Z_HA*H2_D2 + Z_CFM*CFM_D2 + Z_CF2*CF2_D2
      VS2(2,4) = 0.5*Z_MA*M2_U2 + Z_CFM*CFM_U2 + Z_CF2*CF2_U2 + Z_U2
      VS2(2,5) =                                                Z_X2
C
      VSM(2)   = 0.5*Z_MA*M1_MS + Z_CFM*CFM_MS + Z_CF1*CF1_MS
     &         + 0.5*Z_MA*M2_MS                + Z_CF2*CF2_MS
      VSR(2)   =                  Z_CFM*CFM_RE + Z_CF1*CF1_RE
     &                                         + Z_CF2*CF2_RE
      VSX(2)   = 0.
      VSREZ(2) = -REZT
C
C**** Set up shape parameter equation
C
      XOT1 = X1/T1
      XOT2 = X2/T2
C
      HA  = 0.5*(H1  + H2 )
      HSA = 0.5*(HS1 + HS2)
      HCA = 0.5*(HC1 + HC2)
      HWA = 0.5*(DW1/T1 + DW2/T2)
C
      DIX = (1.0-UPW)*DI1*XOT1 + UPW*DI2*XOT2
      CFX = (1.0-UPW)*CF1*XOT1 + UPW*CF2*XOT2
      DIX_UPW = DI2*XOT2 - DI1*XOT1
      CFX_UPW = CF2*XOT2 - CF1*XOT1
C
      BTMP = 2.0*HCA/HSA + 1.0 - HA - HWA
C
      REZH  = HLOG + BTMP*ULOG + XLOG*(0.5*CFX-DIX)
      Z_CFX =  XLOG*0.5
      Z_DIX = -XLOG
      Z_HCA = 2.0*ULOG/HSA
      Z_HA  = -ULOG
      Z_HWA = -ULOG
      Z_XL  = DDLOG * (0.5*CFX-DIX)
      Z_UL  = DDLOG * BTMP
      Z_HL  = DDLOG
C
      Z_UPW = Z_CFX*CFX_UPW + Z_DIX*DIX_UPW
C
      Z_HS1 = -HCA*ULOG/HSA**2 - Z_HL/HS1
      Z_HS2 = -HCA*ULOG/HSA**2 + Z_HL/HS2
C
      Z_CF1 = (1.0-UPW)*Z_CFX*XOT1
      Z_CF2 =      UPW *Z_CFX*XOT2
      Z_DI1 = (1.0-UPW)*Z_DIX*XOT1
      Z_DI2 =      UPW *Z_DIX*XOT2
C
      Z_T1 = (1.0-UPW)*(Z_CFX*CF1 + Z_DIX*DI1)*(-XOT1/T1)
      Z_T2 =      UPW *(Z_CFX*CF2 + Z_DIX*DI2)*(-XOT2/T2)
      Z_X1 = (1.0-UPW)*(Z_CFX*CF1 + Z_DIX*DI1)/ T1        - Z_XL/X1
      Z_X2 =      UPW *(Z_CFX*CF2 + Z_DIX*DI2)/ T2        + Z_XL/X2
      Z_U1 =                                              - Z_UL/U1
      Z_U2 =                                                Z_UL/U2
C
      Z_T1 = Z_T1 + Z_HWA*0.5*(-DW1/T1**2)
      Z_T2 = Z_T2 + Z_HWA*0.5*(-DW2/T2**2)
C
      VS1(3,1) =                               Z_DI1*DI1_S1
      VS1(3,2) = Z_HS1*HS1_T1 + Z_CF1*CF1_T1 + Z_DI1*DI1_T1 + Z_T1
      VS1(3,3) = Z_HS1*HS1_D1 + Z_CF1*CF1_D1 + Z_DI1*DI1_D1
      VS1(3,4) = Z_HS1*HS1_U1 + Z_CF1*CF1_U1 + Z_DI1*DI1_U1 + Z_U1
      VS1(3,5) =                                              Z_X1
      VS2(3,1) =                               Z_DI2*DI2_S2
      VS2(3,2) = Z_HS2*HS2_T2 + Z_CF2*CF2_T2 + Z_DI2*DI2_T2 + Z_T2
      VS2(3,3) = Z_HS2*HS2_D2 + Z_CF2*CF2_D2 + Z_DI2*DI2_D2
      VS2(3,4) = Z_HS2*HS2_U2 + Z_CF2*CF2_U2 + Z_DI2*DI2_U2 + Z_U2
      VS2(3,5) =                                              Z_X2
      VSM(3)   = Z_HS1*HS1_MS + Z_CF1*CF1_MS + Z_DI1*DI1_MS
     &         + Z_HS2*HS2_MS + Z_CF2*CF2_MS + Z_DI2*DI2_MS
      VSR(3)   = Z_HS1*HS1_RE + Z_CF1*CF1_RE + Z_DI1*DI1_RE
     &         + Z_HS2*HS2_RE + Z_CF2*CF2_RE + Z_DI2*DI2_RE
C
      VS1(3,2) = VS1(3,2) + 0.5*(Z_HCA*HC1_T1+Z_HA*H1_T1) + Z_UPW*UPW_T1
      VS1(3,3) = VS1(3,3) + 0.5*(Z_HCA*HC1_D1+Z_HA*H1_D1) + Z_UPW*UPW_D1
      VS1(3,4) = VS1(3,4) + 0.5*(Z_HCA*HC1_U1           ) + Z_UPW*UPW_U1
      VS2(3,2) = VS2(3,2) + 0.5*(Z_HCA*HC2_T2+Z_HA*H2_T2) + Z_UPW*UPW_T2
      VS2(3,3) = VS2(3,3) + 0.5*(Z_HCA*HC2_D2+Z_HA*H2_D2) + Z_UPW*UPW_D2
      VS2(3,4) = VS2(3,4) + 0.5*(Z_HCA*HC2_U2           ) + Z_UPW*UPW_U2
C
      VSM(3)   = VSM(3)   + 0.5*(Z_HCA*HC1_MS           ) + Z_UPW*UPW_MS
     &                    + 0.5*(Z_HCA*HC2_MS           )
C
      VSX(3)   = 0.
      VSREZ(3) = -REZH
C
      RETURN
      END


 
      SUBROUTINE DAMPL( HK, TH, RT, AX, AX_HK, AX_TH, AX_RT )
C==============================================================
C     Amplification rate routine for envelope e^n method.
C     Reference:
C                Drela, M., Giles, M.,
C               "Viscous/Inviscid Analysis of Transonic and
C                Low Reynolds Number Airfoils",
C                AIAA Journal, Oct. 1987.
C
C     NEW VERSION.   March 1991       (latest bug fix  July 93)
C          - m(H) correlation made more accurate up to H=20
C          - for H > 5, non-similar profiles are used 
C            instead of Falkner-Skan profiles.  These 
C            non-similar profiles have smaller reverse 
C            velocities, are more representative of typical 
C            separation bubble profiles.
C--------------------------------------------------------------
C
C     input :   HK     kinematic shape parameter
C               TH     momentum thickness
C               RT     momentum-thickness Reynolds number
C
C     output:   AX     envelope spatial amplification rate
C               AX_(.) sensitivity of AX to parameter (.)
C
C
C     Usage: The log of the envelope amplitude N(x) is
C            calculated by integrating AX (= dN/dx) with
C            respect to the streamwise distance x.
C                      x
C                     /
C              N(x) = | AX(H(x),Th(x),Rth(x)) dx
C                     /
C                      0
C            The integration can be started from the leading
C            edge since AX will be returned as zero when RT
C            is below the critical Rtheta.  Transition occurs
C            when N(x) reaches Ncrit (Ncrit= 9 is "standard").
C==============================================================
      IMPLICIT REAL (A-H,M,O-Z)
ccc   DATA DGR / 0.04 /
      DATA DGR / 0.08 /
C
      HMI = 1.0/(HK - 1.0)
      HMI_HK = -HMI**2
C
C---- log10(Critical Rth) - H   correlation for Falkner-Skan profiles
      AA    = 2.492*HMI**0.43
      AA_HK =   (AA/HMI)*0.43 * HMI_HK
C
      BB    = TANH(14.0*HMI - 9.24)
      BB_HK = (1.0 - BB*BB) * 14.0 * HMI_HK
C
      GRCRIT = AA    + 0.7*(BB + 1.0)
      GRC_HK = AA_HK + 0.7* BB_HK
C
C
      GR = LOG10(RT)
      GR_RT = 1.0 / (2.3025851*RT)
C
      IF(GR .LT. GRCRIT-DGR) THEN
C
C----- no amplification for Rtheta < Rcrit
       AX    = 0.
       AX_HK = 0.
       AX_TH = 0.
       AX_RT = 0.
C
      ELSE
C
C----- Set steep cubic ramp used to turn on AX smoothly as Rtheta 
C-     exceeds Rcrit (previously, this was done discontinuously).
C-     The ramp goes between  -DGR < log10(Rtheta/Rcrit) < DGR
C
       RNORM = (GR - (GRCRIT-DGR)) / (2.0*DGR)
       RN_HK =     -  GRC_HK       / (2.0*DGR)
       RN_RT =  GR_RT              / (2.0*DGR)
C
       IF(RNORM .GE. 1.0) THEN
        RFAC    = 1.0
        RFAC_HK = 0.
        RFAC_RT = 0.
       ELSE
        RFAC    = 3.0*RNORM**2 - 2.0*RNORM**3
        RFAC_RN = 6.0*RNORM    - 6.0*RNORM**2
C
        RFAC_HK = RFAC_RN*RN_HK
        RFAC_RT = RFAC_RN*RN_RT
       ENDIF
C
C----- Amplification envelope slope correlation for Falkner-Skan
       ARG    = 3.87*HMI    - 2.52
       ARG_HK = 3.87*HMI_HK
C
       EX    = EXP(-ARG**2)
       EX_HK = EX * (-2.0*ARG*ARG_HK)
C
       DADR    = 0.028*(HK-1.0) - 0.0345*EX
       DADR_HK = 0.028          - 0.0345*EX_HK
C
C----- new m(H) correlation    1 March 91
       AF = -0.05 + 2.7*HMI -  5.5*HMI**2 + 3.0*HMI**3
       AF_HMI =     2.7     - 11.0*HMI    + 9.0*HMI**2
       AF_HK = AF_HMI*HMI_HK
C
       AX    = (AF   *DADR/TH                ) * RFAC
       AX_HK = (AF_HK*DADR/TH + AF*DADR_HK/TH) * RFAC
     &       + (AF   *DADR/TH                ) * RFAC_HK
       AX_TH = -AX/TH
       AX_RT = (AF   *DADR/TH                ) * RFAC_RT
C
      ENDIF
C
      RETURN
      END ! DAMPL


 
      SUBROUTINE DAMPL2( HK, TH, RT, AX, AX_HK, AX_TH, AX_RT )
C==============================================================
C     Amplification rate routine for modified envelope e^n method.
C     Reference: 
C                Drela, M., Giles, M.,
C               "Viscous/Inviscid Analysis of Transonic and 
C                Low Reynolds Number Airfoils", 
C                AIAA Journal, Oct. 1987.
C
C     NEWER VERSION.   Nov 1996
C          - Amplification rate changes to the Orr-Sommerfeld
C              maximum ai(H,Rt) function for H > 4 .
C          - This implicitly assumes that the frequency range
C              (around w = 0.09 Ue/theta) which experiences this 
C              maximum amplification rate contains the currently
C              most-amplified frequency.
C--------------------------------------------------------------
C
C     input :   HK     kinematic shape parameter
C               TH     momentum thickness
C               RT     momentum-thickness Reynolds number
C
C     output:   AX     envelope spatial amplification rate
C               AX_(.) sensitivity of AX to parameter (.)
C
C
C     Usage: The log of the envelope amplitude N(x) is 
C            calculated by integrating AX (= dN/dx) with 
C            respect to the streamwise distance x.
C                      x
C                     /
C              N(x) = | AX(H(x),Th(x),Rth(x)) dx
C                     /
C                      0
C            The integration can be started from the leading
C            edge since AX will be returned as zero when RT
C            is below the critical Rtheta.  Transition occurs
C            when N(x) reaches Ncrit (Ncrit= 9 is "standard").
C==============================================================
      IMPLICIT REAL (A-H,M,O-Z)
      DATA DGR / 0.08 /
      DATA HK1, HK2 / 3.5, 4.0 /
C
      HMI = 1.0/(HK - 1.0)
      HMI_HK = -HMI**2
C
C---- log10(Critical Rth) -- H   correlation for Falkner-Skan profiles
      AA    = 2.492*HMI**0.43
      AA_HK =   (AA/HMI)*0.43 * HMI_HK
C
      BB    = TANH(14.0*HMI - 9.24)
      BB_HK = (1.0 - BB*BB) * 14.0 * HMI_HK
C
      GRC = AA    + 0.7*(BB + 1.0)
      GRC_HK = AA_HK + 0.7* BB_HK
C
C
      GR = LOG10(RT)
      GR_RT = 1.0 / (2.3025851*RT)
C
      IF(GR .LT. GRC-DGR) THEN
C
C----- no amplification for Rtheta < Rcrit
       AX    = 0.
       AX_HK = 0.
       AX_TH = 0.
       AX_RT = 0.
C
      ELSE
C
C----- Set steep cubic ramp used to turn on AX smoothly as Rtheta 
C-     exceeds Rcrit (previously, this was done discontinuously).
C-     The ramp goes between  -DGR < log10(Rtheta/Rcrit) < DGR
C
       RNORM = (GR - (GRC-DGR)) / (2.0*DGR)
       RN_HK =     -  GRC_HK       / (2.0*DGR)
       RN_RT =  GR_RT              / (2.0*DGR)
C
       IF(RNORM .GE. 1.0) THEN
        RFAC    = 1.0
        RFAC_HK = 0.
        RFAC_RT = 0.
       ELSE
        RFAC    = 3.0*RNORM**2 - 2.0*RNORM**3
        RFAC_RN = 6.0*RNORM    - 6.0*RNORM**2
C
        RFAC_HK = RFAC_RN*RN_HK
        RFAC_RT = RFAC_RN*RN_RT
       ENDIF
C
C
C----- set envelope amplification rate with respect to Rtheta
C-       DADR = d(N)/d(Rtheta) = f(H)
C
       ARG    = 3.87*HMI    - 2.52
       ARG_HK = 3.87*HMI_HK
C
       EX    = EXP(-ARG**2)
       EX_HK = EX * (-2.0*ARG*ARG_HK)
C
       DADR    = 0.028*(HK-1.0) - 0.0345*EX
       DADR_HK = 0.028          - 0.0345*EX_HK
C
C
C----- set conversion factor from d/d(Rtheta) to d/dx
C-       AF = Theta d(Rtheta)/dx = f(H)
C
       BRG = -20.0*HMI
       AF = -0.05 + 2.7*HMI -  5.5*HMI**2 + 3.0*HMI**3 + 0.1*EXP(BRG)
       AF_HMI =     2.7     - 11.0*HMI    + 9.0*HMI**2 - 2.0*EXP(BRG)
       AF_HK = AF_HMI*HMI_HK
C
C
C----- set amplification rate with respect to x, 
C-     with RFAC shutting off amplification when below Rcrit
C
       AX    = (AF   *DADR/TH                ) * RFAC
       AX_HK = (AF_HK*DADR/TH + AF*DADR_HK/TH) * RFAC
     &       + (AF   *DADR/TH                ) * RFAC_HK
       AX_TH = -AX/TH
       AX_RT = (AF   *DADR/TH                ) * RFAC_RT
C
      ENDIF
C
      IF(HK .LT. HK1) RETURN
C
C---- non-envelope max-amplification correction for separated profiles
C
      HNORM = (HK - HK1) / (HK2 - HK1)
      HN_HK =       1.0  / (HK2 - HK1)
C
C---- set blending fraction HFAC = 0..1 over HK1 < HK < HK2
      IF(HNORM .GE. 1.0) THEN
       HFAC = 1.0
       HF_HK = 0.
      ELSE
       HFAC  =  3.0*HNORM**2 - 2.0*HNORM**3
       HF_HK = (6.0*HNORM    - 6.0*HNORM**2)*HN_HK
      ENDIF
C
C---- "normal" envelope amplification rate AX1
      AX1    = AX
      AX1_HK = AX_HK
      AX1_TH = AX_TH
      AX1_RT = AX_RT
C
C---- set modified amplification rate AX2
      GR0 = 0.30 + 0.35 * EXP(-0.15*(HK-5.0))
      GR0_HK =   - 0.35 * EXP(-0.15*(HK-5.0)) * 0.15
C
      TNR = TANH(1.2*(GR - GR0))
      TNR_RT =  (1.0 - TNR**2)*1.2*GR_RT
      TNR_HK = -(1.0 - TNR**2)*1.2*GR0_HK
C
      AX2    = (0.086*TNR    -     0.25/(HK-1.0)**1.5) / TH
      AX2_HK = (0.086*TNR_HK + 1.5*0.25/(HK-1.0)**2.5) / TH
      AX2_RT = (0.086*TNR_RT                         ) / TH
      AX2_TH = -AX2/TH
C
      IF(AX2 .LT. 0.0) THEN
       AX2    = 0.0
       AX2_HK = 0.
       AX2_RT = 0.
       AX2_TH = 0.
      ENDIF
C
C---- blend the two amplification rates
      AX    = HFAC*AX2    + (1.0 - HFAC)*AX1
      AX_HK = HFAC*AX2_HK + (1.0 - HFAC)*AX1_HK + HF_HK*(AX2-AX1)
      AX_RT = HFAC*AX2_RT + (1.0 - HFAC)*AX1_RT
      AX_TH = HFAC*AX2_TH + (1.0 - HFAC)*AX1_TH
C
      RETURN
      END ! DAMPL2

 
 
      SUBROUTINE HKIN( H, MSQ, HK, HK_H, HK_MSQ )
      REAL MSQ
C
C---- calculate kinematic shape parameter (assuming air)
C     (from Whitfield )
      HK     =    (H - 0.29*MSQ)/(1.0 + 0.113*MSQ)
      HK_H   =     1.0          /(1.0 + 0.113*MSQ)
      HK_MSQ = (-.29 - 0.113*HK)/(1.0 + 0.113*MSQ)
C
      RETURN
      END
 


      SUBROUTINE DIL( HK, RT, DI, DI_HK, DI_RT )
C
C---- Laminar dissipation function  ( 2 CD/H* )     (from Falkner-Skan)
      IF(HK.LT.4.0) THEN
       DI    = ( 0.00205  *  (4.0-HK)**5.5 + 0.207 ) / RT
       DI_HK = ( -.00205*5.5*(4.0-HK)**4.5         ) / RT
      ELSE
       HKB = HK - 4.0
       DEN = 1.0 + 0.02*HKB**2
       DI    = ( -.0016  *  HKB**2  /DEN   + 0.207             ) / RT
       DI_HK = ( -.0016*2.0*HKB*(1.0/DEN - 0.02*HKB**2/DEN**2) ) / RT
      ENDIF
      DI_RT = -DI/RT
C
      RETURN
      END


      SUBROUTINE DILW( HK, RT, DI, DI_HK, DI_RT )
      REAL MSQ
C
      MSQ = 0.
      CALL HSL( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
C
C---- Laminar wake dissipation function  ( 2 CD/H* )
      RCD    =  1.10 * (1.0 - 1.0/HK)**2  / HK
      RCD_HK = -1.10 * (1.0 - 1.0/HK)*2.0 / HK**3
     &       - RCD/HK
C
      DI    = 2.0*RCD   /(HS*RT)
      DI_HK = 2.0*RCD_HK/(HS*RT) - (DI/HS)*HS_HK
      DI_RT = -DI/RT             - (DI/HS)*HS_RT
C
      RETURN
      END


      SUBROUTINE HSL( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
      REAL MSQ
C
C---- Laminar HS correlation
      IF(HK.LT.4.35) THEN
       TMP = HK - 4.35
       HS    = 0.0111*TMP**2/(HK+1.0)
     &       - 0.0278*TMP**3/(HK+1.0)  + 1.528
     &       - 0.0002*(TMP*HK)**2
       HS_HK = 0.0111*(2.0*TMP    - TMP**2/(HK+1.0))/(HK+1.0)
     &       - 0.0278*(3.0*TMP**2 - TMP**3/(HK+1.0))/(HK+1.0)
     &       - 0.0002*2.0*TMP*HK * (TMP + HK)
      ELSE
       HS2 = 0.015
c       HS2 = 0.09
       HS    = HS2*    (HK-4.35)**2/HK + 1.528
       HS_HK = HS2*2.0*(HK-4.35)   /HK
     &       - HS2*    (HK-4.35)**2/HK**2
      ENDIF
C
      HS_RT  = 0.
      HS_MSQ = 0.
C
      RETURN
      END


      SUBROUTINE CFL( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
      REAL MSQ
C
C---- Laminar skin friction function  ( Cf )    ( from Falkner-Skan )
      IF(HK.LT.5.5) THEN
       TMP = (5.5-HK)**3 / (HK+1.0)
       CF    = ( 0.0727*TMP                      - 0.07       )/RT
       CF_HK = ( -.0727*TMP*3.0/(5.5-HK) - 0.0727*TMP/(HK+1.0))/RT
      ELSE
       TMP = 1.0 - 1.0/(HK-4.5)
       CF    = ( 0.015*TMP**2      - 0.07  ) / RT
       CF_HK = ( 0.015*TMP*2.0/(HK-4.5)**2 ) / RT
      ENDIF
      CF_RT = -CF/RT
      CF_MSQ = 0.0
C
      RETURN
      END



      SUBROUTINE DIT( HS, US, CF, ST, DI, DI_HS, DI_US, DI_CF, DI_ST )
C
C---- Turbulent dissipation function  ( 2 CD/H* )
      DI    =  ( 0.5*CF*US + ST*ST*(1.0-US) ) * 2.0/HS
      DI_HS = -( 0.5*CF*US + ST*ST*(1.0-US) ) * 2.0/HS**2
      DI_US =  ( 0.5*CF    - ST*ST          ) * 2.0/HS
      DI_CF =  ( 0.5   *US                  ) * 2.0/HS
      DI_ST =  (            2.0*ST*(1.0-US) ) * 2.0/HS
C
      RETURN
      END


      SUBROUTINE HST( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
      IMPLICIT REAL (A-H,M,O-Z)
C
C---- Turbulent HS correlation
C
      DATA HSMIN, DHSINF / 1.500, 0.015 /
C
C---- ###  12/4/94
C---- limited Rtheta dependence for Rtheta < 200
C
C
      IF(RT.GT.400.0) THEN
       HO    = 3.0 + 400.0/RT
       HO_RT =     - 400.0/RT**2
      ELSE
       HO    = 4.0
       HO_RT = 0.
      ENDIF
C
      IF(RT.GT.200.0) THEN
       RTZ    = RT
       RTZ_RT = 1.
      ELSE
       RTZ    = 200.0
       RTZ_RT = 0.
      ENDIF
C
      IF(HK.LT.HO) THEN
C----- attached branch
C=======================================================
C----- old correlation
C-     (from Swafford profiles)
c       SRT = SQRT(RT)
c       HEX = (HO-HK)**1.6
c       RTMP = 0.165 - 1.6/SRT
c       HS    = HSMIN + 4.0/RT + RTMP*HEX/HK
c       HS_HK = RTMP*HEX/HK*(-1.6/(HO-HK) - 1.0/HK)
c       HS_RT = -4.0/RT**2 + HEX/HK*0.8/SRT/RT
c     &             + RTMP*HEX/HK*1.6/(HO-HK)*HO_RT
C=======================================================
C----- new correlation  29 Nov 91
C-     (from  arctan(y+) + Schlichting  profiles)
       HR    = ( HO - HK)/(HO-1.0)
       HR_HK =      - 1.0/(HO-1.0)
       HR_RT = (1.0 - HR)/(HO-1.0) * HO_RT
       HS    = (2.0-HSMIN-4.0/RTZ)*HR**2  * 1.5/(HK+0.5) + HSMIN
     &       + 4.0/RTZ
       HS_HK =-(2.0-HSMIN-4.0/RTZ)*HR**2  * 1.5/(HK+0.5)**2
     &       + (2.0-HSMIN-4.0/RTZ)*HR*2.0 * 1.5/(HK+0.5) * HR_HK
       HS_RT = (2.0-HSMIN-4.0/RTZ)*HR*2.0 * 1.5/(HK+0.5) * HR_RT
     &       + (HR**2 * 1.5/(HK+0.5) - 1.0)*4.0/RTZ**2 * RTZ_RT
C
      ELSE
C
C----- separated branch
       GRT = LOG(RTZ)
       HDIF = HK - HO 
       RTMP = HK - HO + 4.0/GRT
       HTMP    = 0.007*GRT/RTMP**2 + DHSINF/HK
       HTMP_HK = -.014*GRT/RTMP**3 - DHSINF/HK**2
       HTMP_RT = -.014*GRT/RTMP**3 * (-HO_RT - 4.0/GRT**2/RTZ * RTZ_RT)
     &         + 0.007    /RTMP**2 / RTZ * RTZ_RT
       HS    = HDIF**2 * HTMP + HSMIN + 4.0/RTZ
       HS_HK = HDIF*2.0* HTMP
     &       + HDIF**2 * HTMP_HK
       HS_RT = HDIF**2 * HTMP_RT      - 4.0/RTZ**2 * RTZ_RT
     &       + HDIF*2.0* HTMP * (-HO_RT)
C
      ENDIF
C
C---- fudge HS slightly to make sure   HS -> 2   as   HK -> 1
C-    (unnecessary with new correlation)
c      HTF    = 0.485/9.0 * (HK-4.0)**2/HK  +  1.515
c      HTF_HK = 0.485/9.0 * (1.0-16.0/HK**2)
c      ARG = MAX( 10.0*(1.0 - HK) , -15.0 )
c      HXX = EXP(ARG)
c      HXX_HK = -10.0*HXX
cC
c      HS_HK  = (1.0-HXX)*HS_HK  +  HXX*HTF_HK
c     &       + (        -HS     +      HTF    )*HXX_HK
c      HS_RT  = (1.0-HXX)*HS_RT
c      HS     = (1.0-HXX)*HS     +  HXX*HTF
C
C---- Whitfield's minor additional compressibility correction
      FM = 1.0 + 0.014*MSQ
      HS     = ( HS + 0.028*MSQ ) / FM
      HS_HK  = ( HS_HK          ) / FM
      HS_RT  = ( HS_RT          ) / FM
      HS_MSQ = 0.028/FM  -  0.014*HS/FM
C
      RETURN
      END
 
 
 
      SUBROUTINE CFT( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
      IMPLICIT REAL (A-H,M,O-Z)
      INCLUDE 'BLPAR.INC'
C
      DATA GAM /1.4/
C
C---- Turbulent skin friction function  ( Cf )    (Coles)
      GM1 = GAM - 1.0
      FC = SQRT(1.0 + 0.5*GM1*MSQ)
      GRT = LOG(RT/FC)
      GRT = MAX(GRT,3.0)
C
      GEX = -1.74 - 0.31*HK
C
      ARG = -1.33*HK
      ARG = MAX(-20.0, ARG )
C
      THK = TANH(4.0 - HK/0.875)
C
      CFO =  CFFAC * 0.3*EXP(ARG) * (GRT/2.3026)**GEX
      CF     = ( CFO  +  1.1E-4*(THK-1.0) ) / FC
      CF_HK  = (-1.33*CFO - 0.31*LOG(GRT/2.3026)*CFO
     &         - 1.1E-4*(1.0-THK**2) / 0.875    ) / FC
      CF_RT  = GEX*CFO/(FC*GRT) / RT
      CF_MSQ = GEX*CFO/(FC*GRT) * (-0.25*GM1/FC**2) - 0.25*GM1*CF/FC**2
C
      RETURN
      END ! CFT


 
      SUBROUTINE HCT( HK, MSQ, HC, HC_HK, HC_MSQ )
      REAL MSQ
C
C---- density shape parameter    (from Whitfield)
      HC     = MSQ * (0.064/(HK-0.8) + 0.251)
      HC_HK  = MSQ * (-.064/(HK-0.8)**2     )
      HC_MSQ =        0.064/(HK-0.8) + 0.251
C
      RETURN
      END
 
 
C***********************************************************************
C    Module:  xfoil.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************


      SUBROUTINE INIT
C---------------------------------------------------
C     Variable initialization/default routine.
C     See file XFOIL.INC for variable description.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      PI = 4.0*ATAN(1.0)
      HOPI = 0.50/PI
      QOPI = 0.25/PI
      DTOR = PI/180.0
C
C---- default Cp/Cv (air)
      GAMMA = 1.4
      GAMM1 = GAMMA - 1.0
C
C---- set unity freestream speed
      QINF = 1.0
C
C---- initialize freestream Mach number to zero
      MATYP = 1
      MINF1 = 0.
C
      ALFA = 0.0
      COSA = 1.0
      SINA = 0.0
C
      DO 10 I=1, IQX
        GAMU(I,1) = 0.
        GAMU(I,2) = 0.
        GAM(I) = 0.
        GAM_A(I) = 0.
   10 CONTINUE
      PSIO = 0.
C
      CL = 0.
      CM = 0.
      CD = 0.
C
      SIGTE = 0.0
      GAMTE = 0.0
      SIGTE_A = 0.
      GAMTE_A = 0.
C
      DO 20 I=1, IZX
        SIG(I) = 0.
   20 CONTINUE
C
      NQSP = 0
      DO 30 K=1, IPX
        ALQSP(K) = 0.
        CLQSP(K) = 0.
        CMQSP(K) = 0.
        DO 302 I=1, IBX
          QSPEC(I,K) = 0.
 302    CONTINUE
 30   CONTINUE
C
      AWAKE = 0.0
      AVISC = 0.0
C
      KIMAGE = 1
      YIMAGE = -10.0
      LIMAGE = .FALSE.
C
C---- output coordinate file delimiters
C-    KDELIM =  0  blanks
C-           =  1  commas
C-           =  2  tabs
      KDELIM = 0
C
      LGAMU  = .FALSE.
      LQINU  = .FALSE.
      LVISC  = .FALSE.
      LWAKE  = .FALSE.
      LPACC  = .FALSE.
      LBLINI = .FALSE.
      LIPAN  = .FALSE.
      LQAIJ  = .FALSE.
      LADIJ  = .FALSE.
      LWDIJ  = .FALSE.
      LCPXX  = .FALSE.
      LQVDES = .FALSE.
      LQSPEC = .FALSE.
      LQREFL = .FALSE.
      LVCONV = .FALSE.
      LCPREF = .FALSE.
      LFOREF = .FALSE.
      LPFILE = .FALSE.
      LPFILX = .FALSE.
      LPPSHO = .FALSE.
      LBFLAP = .FALSE.
      LFLAP  = .FALSE.
      LEIW   = .FALSE.
      LSCINI = .FALSE.
      LPLOT  = .FALSE.
      LCLIP  = .FALSE.
      LVLAB  = .TRUE.
      LCMINP = .FALSE.
      LHMOMP = .FALSE.
      LFREQP = .TRUE.
C
      LCURS  = .TRUE.
      LLAND  = .TRUE.
      LGSAME = .FALSE.
C
      LGPARM = .TRUE.
      LPLCAM = .FALSE.
C
C---- input airfoil will not be normalized
      LNORM = .FALSE.
C
C---- airfoil will not be forced symmetric
      LQSYM = .FALSE.
      LGSYM = .FALSE.
C
C---- endpoint slopes will be matched
      LQSLOP = .TRUE.
      LGSLOP = .TRUE.
      LCSLOP = .TRUE.
C
C---- grids on Qspec(s) and buffer airfoil geometry plots will be plotted
      LQGRID = .TRUE.
      LGGRID = .TRUE.
      LGTICK = .TRUE.
C
C---- no grid on Cp plots
      LCPGRD = .FALSE.
C
C---- overlay inviscid Cp on viscous Cp(x) plot
      LCPINV = .TRUE.

C---- grid and no symbols are to be used on BL variable plots
      LBLGRD = .TRUE.
      LBLSYM = .FALSE.
C
C---- buffer and current airfoil flap hinge coordinates
      XBF = 0.0
      YBF = 0.0
      XOF = 0.0
      YOF = 0.0
C
      NCPREF = 0
C                                               n
C---- circle plane array size (257, or largest 2  + 1 that will fit array size)
      ANN = LOG(FLOAT((2*IQX)-1))/LOG(2.0)
      NN = INT( ANN + 0.00001 )
      NC1 = 2**NN + 1
      NC1 = MIN( NC1 , 257 )
C
C---- default paneling parameters
      NPAN = 160
      CVPAR = 1.0
      CTERAT = 0.15
      CTRRAT = 0.2
C
C---- default paneling refinement zone x/c endpoints
      XSREF1 = 1.0
      XSREF2 = 1.0
      XPREF1 = 1.0
      XPREF2 = 1.0
C
C---- no polars present to begin with
      NPOL = 0
      IPACT = 0
      DO IP = 1, NPX
        PFNAME(IP) = ' '
        PFNAMX(IP) = ' '
      ENDDO
C
C---- no reference polars
      NPOLREF = 0
C
C---- plot aspect ratio, character size
      PLOTAR = 0.55
      CH     = 0.015
C
C---- airfoil node tick-mark size (as fraction of arc length)
      GTICK = 0.0005
C
C---- Cp limits in  Cp vs x  plot
      CPMAX =  1.0
      CPMIN = -2.0
      CPDEL = -0.5
      PFAC = PLOTAR/(CPMAX-CPMIN)
C
C---- Ue limits in  Ue vs x  plot
      UEMAX =  1.8
      UEMIN = -1.0
      UEDEL =  0.2
      UFAC = PLOTAR/(UEMAX-UEMIN)
C
C---- DCp limits in CAMB loading plot
      YPMIN = -0.4
      YPMAX =  0.4
C
C---- scaling factor for Cp vector plot
      VFAC = 0.25
C
C---- offsets and scale factor for airfoil in  Cp vs x  plot
      XOFAIR = 0.09
      YOFAIR = -.01
      FACAIR = 0.70
C
C---- u/Qinf scale factor for profile plotting
      UPRWT = 0.02
C
C---- polar plot options, grid, list, legend, no CDW
      LPGRID = .TRUE.
      LPCDW  = .FALSE.
      LPLIST = .TRUE.
      LPLEGN = .TRUE.
      LAECEN = .FALSE.
      LPCDH   = .FALSE.
      LPCMDOT = .FALSE.
C
C---- axis limits and annotation deltas for polar plot
      CPOLPLF(1,ICD) = 0.0
      CPOLPLF(2,ICD) = 0.04
      CPOLPLF(3,ICD) = 0.01
C        
      CPOLPLF(1,ICL) = 0.
      CPOLPLF(2,ICL) = 1.5
      CPOLPLF(3,ICL) = 0.5
C        
      CPOLPLF(1,ICM) = -0.25
      CPOLPLF(2,ICM) =  0.0
      CPOLPLF(3,ICM) =  0.05
C        
      CPOLPLF(1,IAL) = -4.0
      CPOLPLF(2,IAL) = 10.0
      CPOLPLF(3,IAL) =  2.0
C
C---- widths of plot boxes in polar plot page
      XCDWID = 0.45
      XALWID = 0.25
      XOCWID = 0.20
C
C---- line style and color index for each polar
C
C     1  *****************************  SOLID
C     2  **** **** **** **** **** ****  LONG DASHED
C     3  ** ** ** ** ** ** ** ** ** **  SHORT DASHED
C     4  * * * * * * * * * * * * * * *  DOTTED
C     5  ***** * ***** * ***** * *****  DASH-DOT
C     6  ***** * * ***** * * ***** * *  DASH-DOT-DOT
C     7  ***** * * * ***** * * * *****  DASH-DOT-DOT-DOT
C     8  **** **** * * **** **** * *    DASH-DASH-DOT-DOT
C
C     3  red
C     4  orange
C     5  yellow
C     6  green
C     7  cyan
C     8  blue
C     9  violet
C    10  magenta
C
      DO IP=1, NPX
cc        ILINP(IP) = 1 + MOD(IP-1,8)
cc        ICOLP(IP) = 3 + MOD(IP-1,8)
C
C------ normally solid, going to dashed after IP=7
        ILINP(IP) = 1 + (IP-1)/7
C
C------ skip yellow (hard to see on white background)
        ICOLP(IP) = 3 + MOD(IP-1,7)
        IF(ICOLP(IP) .GE. 5) ICOLP(IP) = ICOLP(IP) + 1
      ENDDO
C
C---- polar variables to be written to polar save file
      IPOL(1) = IAL
      IPOL(2) = ICL
      IPOL(3) = ICD
      IPOL(4) = ICP
      IPOL(5) = ICM
      NIPOL = 5
      NIPOL0 = 5
C
      JPOL(1) = JTN
      JPOL(2) = JTI
      NJPOL = 2
C
C---- default Cm reference location
      XCMREF = 0.25
      YCMREF = 0.
C
C---- default viscous parameters
      RETYP = 1
      REINF1 = 0.
      ACRIT(1) = 9.0
      ACRIT(2) = 9.0
      XSTRIP(1) = 1.0
      XSTRIP(2) = 1.0
      XOCTR(1) = 1.0
      XOCTR(2) = 1.0
      YOCTR(1) = 0.
      YOCTR(2) = 0.
      WAKLEN = 1.0
C
      IDAMP = 0
C
C---- select default BL plotting coordinate (can be changed in VPLO)
      IXBLP = 1   !  x
ccc   IXBLP = 2   !  s
C
C---- set BL calibration parameters
      CALL BLPINI
C
C---- Newton iteration limit
      ITMAX = 20
C
C---- max number of unconverged sequence points for early exit
      NSEQEX = 4
C
C---- drop tolerance for BL system solver
      VACCEL = 0.01
C
C---- inverse-mapping auto-filter level
      FFILT = 0.0
C
C---- default overlay airfoil filename
      ONAME = ' '
C
C---- default filename prefix
      PREFIX = ' '
C
C---- Plotting flag
      IDEV = 1   ! X11 window only
c     IDEV = 2                  ! B&W PostScript output file only (no color)
c     IDEV = 3   ! both X11 and B&W PostScript file
c     IDEV = 4   ! Color PostScript output file only 
c     IDEV = 5   ! both X11 and Color PostScript file 
C
C---- Re-plotting flag (for hardcopy)
c     IDEVRP = 2   ! B&W PostScript
      IDEVRP = 4   ! Color PostScript
C
C---- PostScript output logical unit and file specification
      IPSLU = 0  ! output to file  plot.ps   on LU 4    (default case)
c     IPSLU = ?  ! output to file  plot?.ps  on LU 10+?
C
C---- screen fraction taken up by plot window upon opening
      SCRNFR = 0.80
C
C---- Default plot size in inches
C-    (Default plot window is 11.0 x 8.5)
C-   (Must be smaller than XPAGE if objects are to fit on paper page)
      SIZE = 10.0

C---- plot-window dimensions in inches for plot blowup calculations
C-    currently,  11.0 x 8.5  default window is hard-wired in libPlt
      XPAGE = 11.0
      YPAGE = 8.5
C
C---- page margins in inches
      XMARG = 0.0
      YMARG = 0.0
C
C---- set top and bottom-side colors
cc      ICOLS(1) = 5
cc      ICOLS(2) = 7
      ICOLS(1) = 8
      ICOLS(2) = 3
C
C   3  red
C   4  orange
C   5  yellow
C   6  green
C   7  cyan
C   8  blue
C   9  violet
C  10  magenta
C
C
C      CALL PLINITIALIZE
C
C---- set up color spectrum
      NCOLOR = 64
C      CALL COLORSPECTRUMHUES(NCOLOR,'RYGCBM')
C
C
      NNAME  = 32
      NAME   = '                                '
CCC             12345678901234567890123456789012
C
C---- MSES domain parameters (not used in XFOIL)
      ISPARS = ' -2.0  3.0  -2.5  3.5'
C
C---- set MINF, REINF, based on current CL-dependence
      CALL MRCL(1.0,MINF_CL,REINF_CL)
C
C---- set various compressibility parameters from MINF
      CALL COMSET
C
      RETURN
      END ! INIT
      
      
      SUBROUTINE MRCL(CLS,M_CLS,R_CLS)
C-------------------------------------------
C     Sets actual Mach, Reynolds numbers
C     from unit-CL values and specified CLS
C     depending on MATYP,RETYP flags.
C-------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL M_CLS
C
      CLA = MAX( CLS , 0.000001 )
C
      IF(RETYP.LT.1 .OR. RETYP.GT.3) THEN
        WRITE(*,*) 'MRCL:  Illegal Re(CL) dependence trigger.'
        WRITE(*,*) '       Setting fixed Re.'
        RETYP = 1
      ENDIF
      IF(MATYP.LT.1 .OR. MATYP.GT.3) THEN
        WRITE(*,*) 'MRCL:  Illegal Mach(CL) dependence trigger.'
        WRITE(*,*) '       Setting fixed Mach.'
        MATYP = 1
      ENDIF
C
C
      IF(MATYP.EQ.1) THEN
C
        MINF  = MINF1
        M_CLS = 0.
C
      ELSE IF(MATYP.EQ.2) THEN
C
        MINF  =  MINF1/SQRT(CLA)
        M_CLS = -0.5*MINF/CLA
C
      ELSE IF(MATYP.EQ.3) THEN
C
        MINF  = MINF1
        M_CLS = 0.
C
      ENDIF
C
C
      IF(RETYP.EQ.1) THEN
C
        REINF = REINF1
        R_CLS = 0.
C
      ELSE IF(RETYP.EQ.2) THEN
C
        REINF =  REINF1/SQRT(CLA)
        R_CLS = -0.5*REINF/CLA
C
      ELSE IF(RETYP.EQ.3) THEN
C
        REINF =  REINF1/CLA
        R_CLS = -REINF /CLA
C
      ENDIF
C
C
      IF(MINF .GE. 0.99) THEN
        WRITE(*,*)
        WRITE(*,*) 'MRCL: CL too low for chosen Mach(CL) dependence'
        WRITE(*,*) '      Aritificially limiting Mach to  0.99'
        MINF = 0.99
        M_CLS = 0.
      ENDIF
C
      RRAT = 1.0
      IF(REINF1 .GT. 0.0) RRAT = REINF/REINF1
C
      IF(RRAT .GT. 100.0) THEN
        WRITE(*,*)
        WRITE(*,*) 'MRCL: CL too low for chosen Re(CL) dependence'
        WRITE(*,*) '      Aritificially limiting Re to ',REINF1*100.0
        REINF = REINF1*100.0
        R_CLS = 0.
      ENDIF
C
      RETURN
      END ! MRCL


      SUBROUTINE COMSET
      INCLUDE 'XFOIL.INC'
C
C---- set Karman-Tsien parameter TKLAM
      BETA = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
C
      TKLAM   = MINF**2 / (1.0 + BETA)**2
      TKL_MSQ =     1.0 / (1.0 + BETA)**2
     &    - 2.0*TKLAM/ (1.0 + BETA) * BETA_MSQ
C
C---- set sonic Pressure coefficient and speed
      IF(MINF.EQ.0.0) THEN
       CPSTAR = -999.0
       QSTAR = 999.0
      ELSE
       CPSTAR = 2.0 / (GAMMA*MINF**2)
     &        * (( (1.0 + 0.5*GAMM1*MINF**2)
     &            /(1.0 + 0.5*GAMM1        ))**(GAMMA/GAMM1) - 1.0)
       QSTAR = QINF/MINF
     &       * SQRT( (1.0 + 0.5*GAMM1*MINF**2)
     &              /(1.0 + 0.5*GAMM1        ) )
      ENDIF
C
      RETURN
      END ! COMSET


      SUBROUTINE CPCALC(N,Q,QINF,MINF,CP)
C---------------------------------------------
C     Sets compressible Cp from speed.
C---------------------------------------------
      DIMENSION Q(N),CP(N)
      REAL MINF
C
      LOGICAL DENNEG
C
      BETA = SQRT(1.0 - MINF**2)
      BFAC = 0.5*MINF**2 / (1.0 + BETA)
C
      DENNEG = .FALSE.
C
      DO 20 I=1, N
        CPINC = 1.0 - (Q(I)/QINF)**2
        DEN = BETA + BFAC*CPINC
        CP(I) = CPINC / DEN
        IF(DEN .LE. 0.0) DENNEG = .TRUE.
  20  CONTINUE
C
      IF(DENNEG) THEN
       WRITE(*,*)
       WRITE(*,*) 'CPCALC: Local speed too large. ',
     &            'Compressibility corrections invalid.'
      ENDIF
C
      RETURN
      END ! CPCALC

 
      SUBROUTINE CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, 
     &                  XREF,YREF,
     &                  CL,CM,CDP, CL_ALF,CL_MSQ)
C-----------------------------------------------------------
C     Integrates surface pressures to get CL and CM.
C     Integrates skin friction to get CDF.
C     Calculates dCL/dAlpha for prescribed-CL routines.
C-----------------------------------------------------------
      DIMENSION X(N),Y(N), GAM(N), GAM_A(N)
      REAL MINF
C
ccC---- moment-reference coordinates
cc      XREF = 0.25
cc      YREF = 0.
C
      SA = SIN(ALFA)
      CA = COS(ALFA)
C
      BETA = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
C
      BFAC     = 0.5*MINF**2 / (1.0 + BETA)
      BFAC_MSQ = 0.5         / (1.0 + BETA)
     &         - BFAC        / (1.0 + BETA) * BETA_MSQ
C
      CL = 0.0
      CM = 0.0

      CDP = 0.0
C
      CL_ALF = 0.
      CL_MSQ = 0.
C
      I = 1
      CGINC = 1.0 - (GAM(I)/QINF)**2
      CPG1     = CGINC/(BETA + BFAC*CGINC)
      CPG1_MSQ = -CPG1/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
      CPI_GAM = -2.0*GAM(I)/QINF**2
      CPC_CPI = (1.0 - BFAC*CPG1)/ (BETA + BFAC*CGINC)
      CPG1_ALF = CPC_CPI*CPI_GAM*GAM_A(I)
C
      DO 10 I=1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
C
        CGINC = 1.0 - (GAM(IP)/QINF)**2
        CPG2     = CGINC/(BETA + BFAC*CGINC)
        CPG2_MSQ = -CPG2/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
        CPI_GAM = -2.0*GAM(IP)/QINF**2
        CPC_CPI = (1.0 - BFAC*CPG2)/ (BETA + BFAC*CGINC)
        CPG2_ALF = CPC_CPI*CPI_GAM*GAM_A(IP)
C
        DX = (X(IP) - X(I))*CA + (Y(IP) - Y(I))*SA
        DY = (Y(IP) - Y(I))*CA - (X(IP) - X(I))*SA
        DG = CPG2 - CPG1
C
        AX = (0.5*(X(IP)+X(I))-XREF)*CA + (0.5*(Y(IP)+Y(I))-YREF)*SA
        AY = (0.5*(Y(IP)+Y(I))-YREF)*CA - (0.5*(X(IP)+X(I))-XREF)*SA
        AG = 0.5*(CPG2 + CPG1)
C
        DX_ALF = -(X(IP) - X(I))*SA + (Y(IP) - Y(I))*CA
        AG_ALF = 0.5*(CPG2_ALF + CPG1_ALF)
        AG_MSQ = 0.5*(CPG2_MSQ + CPG1_MSQ)
C
        CL     = CL     + DX* AG
        CDP    = CDP    - DY* AG
        CM     = CM     - DX*(AG*AX + DG*DX/12.0)
     &                  - DY*(AG*AY + DG*DY/12.0)
C
        CL_ALF = CL_ALF + DX*AG_ALF + AG*DX_ALF
        CL_MSQ = CL_MSQ + DX*AG_MSQ
C
        CPG1 = CPG2
        CPG1_ALF = CPG2_ALF
        CPG1_MSQ = CPG2_MSQ
   10 CONTINUE
C
      RETURN
      END ! CLCALC



      SUBROUTINE CDCALC
      INCLUDE 'XFOIL.INC'
C
      SA = SIN(ALFA)
      CA = COS(ALFA)
C
      IF(LVISC .AND. LBLINI) THEN
C
C----- set variables at the end of the wake
       THWAKE = THET(NBL(2),2)
       URAT   = UEDG(NBL(2),2)/QINF
       UEWAKE = UEDG(NBL(2),2) * (1.0-TKLAM) / (1.0 - TKLAM*URAT**2)
       SHWAKE = DSTR(NBL(2),2)/THET(NBL(2),2)
C
C----- extrapolate wake to downstream infinity using Squire-Young relation
C      (reduces errors of the wake not being long enough)
       CD = 2.0*THWAKE * (UEWAKE/QINF)**(0.5*(5.0+SHWAKE))
C
      ELSE
C
       CD = 0.0
C
      ENDIF
C
C---- calculate friction drag coefficient
      CDF = 0.0
      DO 20 IS=1, 2
        DO 205 IBL=3, IBLTE(IS)
          I  = IPAN(IBL  ,IS)
          IM = IPAN(IBL-1,IS)
          DX = (X(I) - X(IM))*CA + (Y(I) - Y(IM))*SA
          CDF = CDF + 0.5*(TAU(IBL,IS)+TAU(IBL-1,IS))*DX * 2.0/QINF**2
 205    CONTINUE
 20   CONTINUE
C
      RETURN
      END ! CDCALC



      SUBROUTINE LOAD
C------------------------------------------------------
C     Reads airfoil file into buffer airfoil
C     and does various initial processesing on it.
C------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C      CHARACTER*(*) FILNAM
C
C      FNAME = FILNAM
C      IF(FNAME(1:1) .EQ. ' ') CALL ASKS('Enter filename^',FNAME)
C      
C      LU = 9
C      CALL AREAD(LU,FNAME,IBX,XB,YB,NB,NAME,ISPARS,ITYPE,1)
C      IF(ITYPE.EQ.0) RETURN
C
C      IF(ITYPE.EQ.1) CALL ASKS('Enter airfoil name^',NAME)
C      CALL STRIP(NAME,NNAME)
C
C---- set default prefix for other filenames
C      KDOT = INDEX(FNAME,'.')
C      IF(KDOT.EQ.0) THEN
C       PREFIX = FNAME
C      ELSE
C       PREFIX = FNAME(1:KDOT-1)
C      ENDIF
C      CALL STRIP(PREFIX,NPREFIX)
C
C---- calculate airfoil area assuming counterclockwise ordering
      AREA = 0.0
      DO 50 I=1, NB
        IP = I+1
        IF(I.EQ.NB) IP = 1
        AREA = AREA + 0.5*(YB(I)+YB(IP))*(XB(I)-XB(IP))
   50 CONTINUE
C
      IF(AREA.GE.0.0) THEN
       LCLOCK = .FALSE.
       WRITE(*,1010) NB
      ELSE
C----- if area is negative (clockwise order), reverse coordinate order
       LCLOCK = .TRUE.
       WRITE(*,1011) NB
       DO 55 I=1, NB/2
         XTMP = XB(NB-I+1)
         YTMP = YB(NB-I+1)
         XB(NB-I+1) = XB(I)
         YB(NB-I+1) = YB(I)
         XB(I) = XTMP
         YB(I) = YTMP
   55  CONTINUE
      ENDIF
C
      IF(LNORM) THEN
       CALL NORM(XB,XBP,YB,YBP,SB,NB)
       WRITE(*,1020)
      ENDIF
C
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB, W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
C
      XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XBTE = 0.5*(XB(1) + XB(NB))
      YBTE = 0.5*(YB(1) + YB(NB))
C
      WRITE(*,1050) XBLE,YBLE, CHORDB,
     &              XBTE,YBTE
C
C---- set reasonable MSES domain parameters for non-MSES coordinate file
      IF(ITYPE.LE.2 .AND. ISPARS.EQ.' ') THEN
        XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
        YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
        XINL = XBLE - 2.0*CHORDB
        XOUT = XBLE + 3.0*CHORDB
        YBOT = YBLE - 2.5*CHORDB
        YTOP = YBLE + 3.5*CHORDB
        XINL = AINT(20.0*ABS(XINL/CHORDB)+0.5)/20.0 * SIGN(CHORDB,XINL)
        XOUT = AINT(20.0*ABS(XOUT/CHORDB)+0.5)/20.0 * SIGN(CHORDB,XOUT)
        YBOT = AINT(20.0*ABS(YBOT/CHORDB)+0.5)/20.0 * SIGN(CHORDB,YBOT)
        YTOP = AINT(20.0*ABS(YTOP/CHORDB)+0.5)/20.0 * SIGN(CHORDB,YTOP)
        WRITE(ISPARS,1005) XINL, XOUT, YBOT, YTOP
 1005   FORMAT(1X, 4F8.2 )
      ENDIF
C
C---- wipe out old flap hinge location
      XBF = 0.0
      YBF = 0.0
      LBFLAP = .FALSE.
C
C---- wipe out off-design alphas, CLs
cc      NALOFF = 0
cc      NCLOFF = 0
C
      RETURN
C...............................................................
 1010 FORMAT(/' Number of input coordinate points:', I4
     &       /' Counterclockwise ordering')
 1011 FORMAT(/' Number of input coordinate points:', I4
     &       /' Clockwise ordering')
 1020 FORMAT(/' Airfoil has been normalized')
 1050 FORMAT(/'  LE  x,y  =', 2F10.5,'  |   Chord =',F10.5
     &       /'  TE  x,y  =', 2F10.5,'  |'                 )
      END ! LOAD


      SUBROUTINE TECALC
C-------------------------------------------
C     Calculates total and projected TE gap 
C     areas and TE panel strengths.
C-------------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- set TE base vector and TE bisector components
      DXTE = X(1) - X(N)
      DYTE = Y(1) - Y(N)
      DXS = 0.5*(-XP(1) + XP(N))
      DYS = 0.5*(-YP(1) + YP(N))
C
C---- normal and streamwise projected TE gap areas
      ANTE = DXS*DYTE - DYS*DXTE
      ASTE = DXS*DXTE + DYS*DYTE
C
C---- total TE gap area
      DSTE = SQRT(DXTE**2 + DYTE**2)
C
      SHARP = DSTE .LT. 0.0001*CHORD
C
      IF(SHARP) THEN
       SCS = 1.0
       SDS = 0.0
      ELSE
       SCS = ANTE/DSTE
       SDS = ASTE/DSTE
      ENDIF
C
C---- TE panel source and vorticity strengths
      SIGTE = 0.5*(GAM(1) - GAM(N))*SCS
      GAMTE = -.5*(GAM(1) - GAM(N))*SDS
C
      SIGTE_A = 0.5*(GAM_A(1) - GAM_A(N))*SCS
      GAMTE_A = -.5*(GAM_A(1) - GAM_A(N))*SDS
C
      RETURN
      END ! TECALC 
      
      
C***********************************************************************
C    Module:  xgdes.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************


      SUBROUTINE GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XF,YF)
      DIMENSION X(N),XP(N),Y(N),YP(N),S(N)
C
      IF(XF .EQ. -999.0)
     &  CALL ASKR('Enter flap hinge x location^',XF)
C
C---- find top and bottom y at hinge x location
      TOPS = S(1) + (X(1) - XF)
      BOTS = S(N) - (X(N) - XF)
      CALL SINVRT(TOPS,XF,X,XP,S,N)      
      CALL SINVRT(BOTS,XF,X,XP,S,N)      
      TOPY = SEVAL(TOPS,Y,YP,S,N)
      BOTY = SEVAL(BOTS,Y,YP,S,N)
C
      WRITE(*,1000) TOPY, BOTY
 1000 FORMAT(/'  Top    surface:  y =', F8.4,'     y/t = 1.0'
     &       /'  Bottom surface:  y =', F8.4,'     y/t = 0.0')
C
      IF(YF .EQ. -999.0)
     & CALL ASKR(
     &  'Enter flap hinge y location (or 999 to specify y/t)^',YF)
C
      IF(YF .EQ. 999.0) THEN
        CALL ASKR('Enter flap hinge relative y/t location^',YREL)
        YF = TOPY*YREL + BOTY*(1.0-YREL)
      ENDIF
C
      RETURN
      END ! GETXYF
      

      SUBROUTINE ABCOPY(LCONF)
      INCLUDE 'XFOIL.INC'
      LOGICAL LCONF
C
      IF(NB.LE.1) THEN
       WRITE(*,*) 'ABCOPY: Buffer airfoil not available.'
       RETURN
      ELSEIF(NB.GT.IQX-5) THEN
       WRITE(*,*) 'Maximum number of panel nodes  : ',IQX-5
       WRITE(*,*) 'Number of buffer airfoil points: ',NB
       WRITE(*,*) 'Current airfoil cannot be set.'
       WRITE(*,*) 'Try executing PANE at Top Level instead.'
       RETURN
      ENDIF
      IF(N.NE.NB) LBLINI = .FALSE.
C
      N = NB
      DO 101 I=1, N
        X(I) = XB(I)
        Y(I) = YB(I)
  101 CONTINUE
      LGSAME = .TRUE.
C
      IF(LBFLAP) THEN
       XOF = XBF
       YOF = YBF
       LFLAP = .TRUE.
      ENDIF
C
C---- strip out doubled points
      I = 1
 102  CONTINUE
      I = I+1
      IF(X(I-1).EQ.X(I) .AND. Y(I-1).EQ.Y(I)) THEN
        DO 104 J=I, N-1
          X(J) = X(J+1)
          Y(J) = Y(J+1)
 104    CONTINUE
        N = N-1
      ENDIF
      IF(I.LT.N) GO TO 102
C
      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XP,S,N)
      CALL SEGSPL(Y,YP,S,N)

      CALL NCALC(X,Y,S,N,NX,NY)

      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD  = SQRT( (XTE-XLE)**2 + (YTE-YLE)**2 )

      CALL TECALC
      CALL APCALC
C
      LGAMU = .FALSE.
      LQINU = .FALSE.
      LWAKE = .FALSE.
      LQAIJ = .FALSE.
      LADIJ = .FALSE.
      LWDIJ = .FALSE.
      LIPAN = .FALSE.
      LVCONV = .FALSE.
      LSCINI = .FALSE.
CCC      LBLINI = .FALSE.
C
      IF(LCONF) WRITE(*,1200) N
 1200 FORMAT(/' Current airfoil nodes set from buffer airfoil nodes (',
     &        I4,' )')
C
      RETURN
      END ! ABCOPY     
      
      
C***********************************************************************
C    Module:  xgeom.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************

 
      SUBROUTINE NORM(X,XP,Y,YP,S,N)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C-----------------------------------------------
C     Scales coordinates to get unit chord
C-----------------------------------------------
C
      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XP,S,N)
      CALL SEGSPL(Y,YP,S,N)
C
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
C
      XMAX = 0.5*(X(1) + X(N))
      XMIN = SEVAL(SLE,X,XP,S,N)
      YMIN = SEVAL(SLE,Y,YP,S,N)
C
      FUDGE = 1.0/(XMAX-XMIN)
      DO 40 I=1, N
        X(I) = (X(I)-XMIN)*FUDGE
        Y(I) = (Y(I)-YMIN)*FUDGE
        S(I) = S(I)*FUDGE
   40 CONTINUE
C
      RETURN
      END


      SUBROUTINE GEOPAR(X,XP,Y,YP,S,N, T,
     &             SLE,CHORD,AREA,RADLE,ANGTE,
     &             EI11A,EI22A,APX1A,APX2A,
     &             EI11T,EI22T,APX1T,APX2T,
     &             THICK,CAMBR)
      DIMENSION X(*), XP(*), Y(*), YP(*), S(*), T(*)
C
      PARAMETER (IBX=600)
      DIMENSION
     &     XCAM(2*IBX), YCAM(2*IBX), YCAMP(2*IBX),
     &     XTHK(2*IBX), YTHK(2*IBX), YTHKP(2*IBX)
C------------------------------------------------------
C     Sets geometric parameters for airfoil shape
C------------------------------------------------------
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
C
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
C
      CHSQ = (XTE-XLE)**2 + (YTE-YLE)**2
      CHORD = SQRT(CHSQ)
C
      CURVLE = CURV(SLE,X,XP,Y,YP,S,N)
C
      RADLE = 0.0
      IF(ABS(CURVLE) .GT. 0.001*(S(N)-S(1))) RADLE = 1.0 / CURVLE
C
      ANG1 = ATAN2( -YP(1) , -XP(1) )
      ANG2 = ATANC(  YP(N) ,  XP(N) , ANG1 )
      ANGTE = ANG2 - ANG1
C

      DO I=1, N
        T(I) = 1.0
      ENDDO
C
      CALL AECALC(N,X,Y,T, 1, 
     &            AREA,XCENA,YCENA,EI11A,EI22A,APX1A,APX2A)
C
      CALL AECALC(N,X,Y,T, 2, 
     &            SLEN,XCENT,YCENT,EI11T,EI22T,APX1T,APX2T)
C
C--- Old, approximate thickness,camber routine (on discrete points only)
      CALL TCCALC(X,XP,Y,YP,S,N, THICK,XTHICK, CAMBR,XCAMBR )
C
C--- More accurate thickness and camber estimates
cc      CALL GETCAM(XCAM,YCAM,NCAM,XTHK,YTHK,NTHK,
cc     &            X,XP,Y,YP,S,N )
cc      CALL GETMAX(XCAM,YCAM,YCAMP,NCAM,XCAMBR,CAMBR)
cc      CALL GETMAX(XTHK,YTHK,YTHKP,NTHK,XTHICK,THICK)
cc      THICK = 2.0*THICK
C
      WRITE(*,1000) THICK,XTHICK,CAMBR,XCAMBR
 1000 FORMAT( ' Max thickness = ',F12.6,'  at x = ',F7.3,
     &       /' Max camber    = ',F12.6,'  at x = ',F7.3)


C
      RETURN
      END ! GEOPAR
      

      SUBROUTINE LEFIND(SLE,X,XP,Y,YP,S,N)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C------------------------------------------------------
C     Locates leading edge spline-parameter value SLE
C
C     The defining condition is
C         
C      (X-XTE,Y-YTE) . (X',Y') = 0     at  S = SLE
C
C     i.e. the surface tangent is normal to the chord
C     line connecting X(SLE),Y(SLE) and the TE point.
C------------------------------------------------------
C
C---- convergence tolerance
      DSEPS = (S(N)-S(1)) * 1.0E-5
C
C---- set trailing edge point coordinates
      XTE = 0.5*(X(1) + X(N))
      YTE = 0.5*(Y(1) + Y(N))
C
C---- get first guess for SLE
      DO 10 I=3, N-2
        DXTE = X(I) - XTE
        DYTE = Y(I) - YTE
        DX = X(I+1) - X(I)
        DY = Y(I+1) - Y(I)
        DOTP = DXTE*DX + DYTE*DY
        IF(DOTP .LT. 0.0) GO TO 11
   10 CONTINUE
C
   11 SLE = S(I)
C
C---- check for sharp LE case
      IF(S(I) .EQ. S(I-1)) THEN
ccc        WRITE(*,*) 'Sharp LE found at ',I,SLE
        RETURN
      ENDIF
C
C---- Newton iteration to get exact SLE value
      DO 20 ITER=1, 50
        XLE  = SEVAL(SLE,X,XP,S,N)
        YLE  = SEVAL(SLE,Y,YP,S,N)
        DXDS = DEVAL(SLE,X,XP,S,N)
        DYDS = DEVAL(SLE,Y,YP,S,N)
        DXDD = D2VAL(SLE,X,XP,S,N)
        DYDD = D2VAL(SLE,Y,YP,S,N)
C
        XCHORD = XLE - XTE
        YCHORD = YLE - YTE
C
C------ drive dot product between chord line and LE tangent to zero
        RES  = XCHORD*DXDS + YCHORD*DYDS
        RESS = DXDS  *DXDS + DYDS  *DYDS
     &       + XCHORD*DXDD + YCHORD*DYDD
C
C------ Newton delta for SLE 
        DSLE = -RES/RESS
C
        DSLE = MAX( DSLE , -0.02*ABS(XCHORD+YCHORD) )
        DSLE = MIN( DSLE ,  0.02*ABS(XCHORD+YCHORD) )
        SLE = SLE + DSLE
        IF(ABS(DSLE) .LT. DSEPS) RETURN
   20 CONTINUE
      WRITE(*,*) 'LEFIND:  LE point not found.  Continuing...'
      SLE = S(I)
      RETURN
      END
      
      
      SUBROUTINE AECALC(N,X,Y,T, ITYPE, 
     &                  AREA,XCEN,YCEN,EI11,EI22,APX1,APX2)
      DIMENSION X(*),Y(*),T(*)
C---------------------------------------------------------------
C     Calculates geometric properties of shape X,Y
C
C     Input:
C       N      number of points
C       X(.)   shape coordinate point arrays
C       Y(.)
C       T(.)   skin-thickness array, used only if ITYPE = 2
C       ITYPE  = 1 ...   integration is over whole area  dx dy
C              = 2 ...   integration is over skin  area   t ds
C
C     Output:
C       XCEN,YCEN  centroid location
C       EI11,EI22  principal moments of inertia
C       APX1,APX2  principal-axis angles
C---------------------------------------------------------------
      DATA PI / 3.141592653589793238 /
C
      SINT  = 0.0
      AINT  = 0.0
      XINT  = 0.0
      YINT  = 0.0
      XXINT = 0.0
      XYINT = 0.0
      YYINT = 0.0
C
      DO 10 IO = 1, N
        IF(IO.EQ.N) THEN
          IP = 1
        ELSE
          IP = IO + 1
        ENDIF
C
        DX =  X(IO) - X(IP)
        DY =  Y(IO) - Y(IP)
        XA = (X(IO) + X(IP))*0.50
        YA = (Y(IO) + Y(IP))*0.50
        TA = (T(IO) + T(IP))*0.50
C
        DS = SQRT(DX*DX + DY*DY)
        SINT = SINT + DS

        IF(ITYPE.EQ.1) THEN
C-------- integrate over airfoil cross-section
          DA = YA*DX
          AINT  = AINT  +       DA
          XINT  = XINT  + XA   *DA
          YINT  = YINT  + YA   *DA/2.0
          XXINT = XXINT + XA*XA*DA
          XYINT = XYINT + XA*YA*DA/2.0
          YYINT = YYINT + YA*YA*DA/3.0
        ELSE
C-------- integrate over skin thickness
          DA = TA*DS
          AINT  = AINT  +       DA
          XINT  = XINT  + XA   *DA
          YINT  = YINT  + YA   *DA
          XXINT = XXINT + XA*XA*DA
          XYINT = XYINT + XA*YA*DA
          YYINT = YYINT + YA*YA*DA
        ENDIF
C
 10   CONTINUE
C
      AREA = AINT
C
      IF(AINT .EQ. 0.0) THEN
        XCEN  = 0.0
        YCEN  = 0.0
        EI11  = 0.0
        EI22  = 0.0
        APX1 = 0.0
        APX2 = ATAN2(1.0,0.0)
        RETURN
      ENDIF
C
C
C---- calculate centroid location
      XCEN = XINT/AINT
      YCEN = YINT/AINT
C
C---- calculate inertias
      EIXX = YYINT - YCEN*YCEN*AINT
      EIXY = XYINT - XCEN*YCEN*AINT
      EIYY = XXINT - XCEN*XCEN*AINT
C
C---- set principal-axis inertias, EI11 is closest to "up-down" bending inertia
      EISQ  = 0.25*(EIXX - EIYY)**2  + EIXY**2
      SGN = SIGN( 1.0 , EIYY-EIXX )
      EI11 = 0.5*(EIXX + EIYY) - SGN*SQRT(EISQ)
      EI22 = 0.5*(EIXX + EIYY) + SGN*SQRT(EISQ)
C
      IF(EI11.EQ.0.0 .OR. EI22.EQ.0.0) THEN
C----- vanishing section stiffness
       APX1 = 0.0
       APX2 = ATAN2(1.0,0.0)
C
      ELSEIF(EISQ/(EI11*EI22) .LT. (0.001*SINT)**4) THEN
C----- rotationally-invariant section (circle, square, etc.)
       APX1 = 0.0
       APX2 = ATAN2(1.0,0.0)
C
      ELSE
C----- normal airfoil section
       C1 = EIXY
       S1 = EIXX-EI11
C
       C2 = EIXY
       S2 = EIXX-EI22
C
       IF(ABS(S1).GT.ABS(S2)) THEN
         APX1 = ATAN2(S1,C1)
         APX2 = APX1 + 0.5*PI
       ELSE
         APX2 = ATAN2(S2,C2)
         APX1 = APX2 - 0.5*PI
       ENDIF

       IF(APX1.LT.-0.5*PI) APX1 = APX1 + PI
       IF(APX1.GT.+0.5*PI) APX1 = APX1 - PI
       IF(APX2.LT.-0.5*PI) APX2 = APX2 + PI
       IF(APX2.GT.+0.5*PI) APX2 = APX2 - PI
C
      ENDIF
C
      RETURN
      END ! AECALC
      
      
      SUBROUTINE TCCALC(X,XP,Y,YP,S,N, 
     &                  THICK,XTHICK, CAMBR,XCAMBR )
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C---------------------------------------------------------------
C     Calculates max thickness and camber at airfoil points
C
C     Note: this routine does not find the maximum camber or 
C           thickness exactly as it only looks at discrete points
C
C     Input:
C       N      number of points
C       X(.)   shape coordinate point arrays
C       Y(.)
C
C     Output:
C       THICK  max thickness
C       CAMBR  max camber
C---------------------------------------------------------------
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD = SQRT((XTE-XLE)**2 + (YTE-YLE)**2)
C
C---- set unit chord-line vector
      DXC = (XTE-XLE) / CHORD
      DYC = (YTE-YLE) / CHORD
C
      THICK = 0.
      XTHICK = 0.
      CAMBR = 0.
      XCAMBR = 0.
C
C---- go over each point, finding the y-thickness and camber
      DO 30 I=1, N
        XBAR = (X(I)-XLE)*DXC + (Y(I)-YLE)*DYC
        YBAR = (Y(I)-YLE)*DXC - (X(I)-XLE)*DYC
C
C------ set point on the opposite side with the same chord x value
        CALL SOPPS(SOPP, S(I), X,XP,Y,YP,S,N, SLE)
        XOPP = SEVAL(SOPP,X,XP,S,N)
        YOPP = SEVAL(SOPP,Y,YP,S,N)
C
        YBAROP = (YOPP-YLE)*DXC - (XOPP-XLE)*DYC
C
        YC = 0.5*(YBAR+YBAROP)
        YT =  ABS(YBAR-YBAROP)
C
        IF(ABS(YC) .GT. ABS(CAMBR)) THEN
         CAMBR = YC
         XCAMBR = XOPP
        ENDIF
        IF(ABS(YT) .GT. ABS(THICK)) THEN
         THICK = YT
         XTHICK = XOPP
        ENDIF
   30 CONTINUE
C
      RETURN
      END ! TCCALC
      
      
      SUBROUTINE SOPPS(SOPP, SI, X,XP,Y,YP,S,N, SLE)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C--------------------------------------------------
C     Calculates arc length SOPP of point 
C     which is opposite of point SI, on the 
C     other side of the airfoil baseline
C--------------------------------------------------
C
C---- reference length for testing convergence
      SLEN = S(N) - S(1)
C
C---- set chordline vector
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD = SQRT((XTE-XLE)**2 + (YTE-YLE)**2)
      DXC = (XTE-XLE) / CHORD
      DYC = (YTE-YLE) / CHORD
C
      IF(SI.LT.SLE) THEN
       IN = 1
       INOPP = N
      ELSE
       IN = N
       INOPP = 1
      ENDIF
      SFRAC = (SI-SLE)/(S(IN)-SLE)
      SOPP = SLE + SFRAC*(S(INOPP)-SLE)
C     
      IF(ABS(SFRAC) .LE. 1.0E-5) THEN
       SOPP = SLE
       RETURN
      ENDIF
C
C---- XBAR = x coordinate in chord-line axes
      XI  = SEVAL(SI , X,XP,S,N)
      YI  = SEVAL(SI , Y,YP,S,N)
      XLE = SEVAL(SLE, X,XP,S,N)
      YLE = SEVAL(SLE, Y,YP,S,N)
      XBAR = (XI-XLE)*DXC + (YI-YLE)*DYC
C
C---- converge on exact opposite point with same XBAR value
      DO 300 ITER=1, 12
        XOPP  = SEVAL(SOPP,X,XP,S,N)
        YOPP  = SEVAL(SOPP,Y,YP,S,N)
        XOPPD = DEVAL(SOPP,X,XP,S,N)
        YOPPD = DEVAL(SOPP,Y,YP,S,N)
C
        RES  = (XOPP -XLE)*DXC + (YOPP -YLE)*DYC - XBAR
        RESD =  XOPPD     *DXC +  YOPPD     *DYC
C
        IF(ABS(RES)/SLEN .LT. 1.0E-5) GO TO 305
        IF(RESD .EQ. 0.0) GO TO 303
C
        DSOPP = -RES/RESD
        SOPP = SOPP + DSOPP
C
        IF(ABS(DSOPP)/SLEN .LT. 1.0E-5) GO TO 305
 300  CONTINUE
 303  WRITE(*,*)
     &      'SOPPS: Opposite-point location failed. Continuing...'
      SOPP = SLE + SFRAC*(S(INOPP)-SLE)
C
 305  CONTINUE
      RETURN
      END ! SOPPS
      
      
      SUBROUTINE CANG(X,Y,N,IPRINT, IMAX,AMAX)
      DIMENSION X(*), Y(*)
C-------------------------------------------------------------------
C     IPRINT=2:   Displays all panel node corner angles
C     IPRINT=1:   Displays max panel node corner angle
C     IPRINT=0:   No display... just returns values
C-------------------------------------------------------------------
C
      AMAX = 0.0
      IMAX = 1
C
C---- go over each point, calculating corner angle
      IF(IPRINT.EQ.2) WRITE(*,1050)
      DO 30 I=2, N-1
        DX1 = X(I) - X(I-1)
        DY1 = Y(I) - Y(I-1)
        DX2 = X(I) - X(I+1)
        DY2 = Y(I) - Y(I+1)
C
C------ allow for doubled points
        IF(DX1.EQ.0.0 .AND. DY1.EQ.0.0) THEN
         DX1 = X(I) - X(I-2)
         DY1 = Y(I) - Y(I-2)
        ENDIF
        IF(DX2.EQ.0.0 .AND. DY2.EQ.0.0) THEN
         DX2 = X(I) - X(I+2)
         DY2 = Y(I) - Y(I+2)
        ENDIF
C
        CROSSP = (DX2*DY1 - DY2*DX1)
     &         / SQRT((DX1**2 + DY1**2) * (DX2**2 + DY2**2))
        ANGL = ASIN(CROSSP)*(180.0/3.1415926)
        IF(IPRINT.EQ.2) WRITE(*,1100) I, X(I), Y(I), ANGL
        IF(ABS(ANGL) .GT. ABS(AMAX)) THEN
         AMAX = ANGL
         IMAX = I
        ENDIF
   30 CONTINUE
C
      IF(IPRINT.GE.1) WRITE(*,1200) AMAX, IMAX, X(IMAX), Y(IMAX)
C
      RETURN
C
 1050 FORMAT(/'  i       x        y      angle')
CCC             120   0.2134  -0.0234   25.322
 1100 FORMAT(1X,I3, 2F9.4, F9.3)
 1200 FORMAT(/' Maximum panel corner angle =', F7.3,
     &        '   at  i,x,y  = ', I3, 2F9.4 )
      END ! CANG
      
      
C***********************************************************************
C    Module:  xoper.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************


      SUBROUTINE SPECAL
C-----------------------------------
C     Converges to specified alpha.
C-----------------------------------
      INCLUDE 'XFOIL.INC'
      REAL MINF_CLM, MSQ_CLM
C
C---- calculate surface vorticity distributions for alpha = 0, 90 degrees
      IF(.NOT.LGAMU .OR. .NOT.LQAIJ) CALL GGCALC
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
C---- superimpose suitably weighted  alpha = 0, 90  distributions
      DO 50 I=1, N
        GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
        GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   50 CONTINUE
      PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
C
      CALL TECALC
      CALL QISET
C
C---- set initial guess for the Newton variable CLM
      CLM = 1.0
C
C---- set corresponding  M(CLM), Re(CLM)
      CALL MRCL(CLM,MINF_CLM,REINF_CLM)
      CALL COMSET
C
C---- set corresponding CL(M)
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
C
C---- iterate on CLM
      DO 100 ITCL=1, 20
C
        MSQ_CLM = 2.0*MINF*MINF_CLM
        DCLM = (CL - CLM)/(1.0 - CL_MSQ*MSQ_CLM)
C
        CLM1 = CLM
        RLX = 1.0
C
C------ under-relaxation loop to avoid driving M(CL) above 1
        DO 90 IRLX=1, 12
C
          CLM = CLM1 + RLX*DCLM
C
C-------- set new freestream Mach M(CLM)
          CALL MRCL(CLM,MINF_CLM,REINF_CLM)
C
C-------- if Mach is OK, go do next Newton iteration
          IF(MATYP.EQ.1 .OR. MINF.EQ.0.0 .OR. MINF_CLM.NE.0.0) GO TO 91
C
          RLX = 0.5*RLX
   90   CONTINUE
   91   CONTINUE
C
C------ set new CL(M)
        CALL COMSET
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
C
        IF(ABS(DCLM).LE.1.0E-6) GO TO 110
C
  100 CONTINUE
      WRITE(*,*) 'SPECAL:  Minf convergence failed'
  110 CONTINUE
C
C---- set final Mach, CL, Cp distributions, and hinge moment
      CALL MRCL(CL,MINF_CL,REINF_CL)
      CALL COMSET
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
      CALL CPCALC(N,QINV,QINF,MINF,CPI)
      IF(LVISC) THEN
       CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
       CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
      ELSE
       CALL CPCALC(N,QINV,QINF,MINF,CPI)
      ENDIF
C      IF(LFLAP) CALL MHINGE
C
      RETURN
      END ! SPECAL
      

      SUBROUTINE VISCAL(NITER1)
C----------------------------------------
C     Converges viscous operating point
C----------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- convergence tolerance
      DATA EPS1 / 1.0E-4 /
C
      NITER = NITER1
C
C---- calculate wake trajectory from current inviscid solution if necessary
      IF(.NOT.LWAKE) THEN
       CALL XYWAKE
      ENDIF
C
C---- set velocities on wake from airfoil vorticity for alpha=0, 90
      CALL QWCALC
C
C---- set velocities on airfoil and wake for initial alpha
      CALL QISET
C
      IF(.NOT.LIPAN) THEN
C
       IF(LBLINI) CALL GAMQV
C
C----- locate stagnation point arc length position and panel index
       CALL STFIND
C
C----- set  BL position -> panel position  pointers
       CALL IBLPAN
C
C----- calculate surface arc length array for current stagnation point location
       CALL XICALC
C
C----- set  BL position -> system line  pointers
       CALL IBLSYS
C
      ENDIF
C
C---- set inviscid BL edge velocity UINV from QINV
      CALL UICALC
C
      IF(.NOT.LBLINI) THEN
C
C----- set initial Ue from inviscid Ue
       DO IBL=1, NBL(1)
         UEDG(IBL,1) = UINV(IBL,1)
       ENDDO
C
       DO IBL=1, NBL(2)
         UEDG(IBL,2) = UINV(IBL,2)
       ENDDO
C
      ENDIF
C
      IF(LVCONV) THEN
C----- set correct CL if converged point exists
       CALL QVFUE
       IF(LVISC) THEN
        CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
        CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
       ELSE
        CALL CPCALC(N,QINV,QINF,MINF,CPI)
       ENDIF
       CALL GAMQV
       CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &             CL,CM,CDP, CL_ALF,CL_MSQ)
       CALL CDCALC
      ENDIF
C
C---- set up source influence matrix if it doesn't exist
      IF(.NOT.LWDIJ .OR. .NOT.LADIJ) CALL QDCALC
C
C---- Newton iteration for entire BL solution
      IF(NITER.EQ.0) THEN
        NITER = 100
        WRITE(*,*) '> NO VALUE FOR NITER WAS SPECIFIED...'
      ENDIF
      
      WRITE(*,*)
      WRITE(*,*) 'Solving BL system ...'
      DO 1000 ITER=1, NITER
C
C------ fill Newton system for BL variables
        CALL SETBL
C
C------ solve Newton system with custom solver
        CALL BLSOLV
C
C------ update BL variables
        CALL UPDATE
C
        IF(LALFA) THEN
C------- set new freestream Mach, Re from new CL
         CALL MRCL(CL,MINF_CL,REINF_CL)
         CALL COMSET
        ELSE
C------- set new inviscid speeds QINV and UINV for new alpha
         CALL QISET
         CALL UICALC
        ENDIF
C
C------ calculate edge velocities QVIS(.) from UEDG(..)
        CALL QVFUE
C
C------ set GAM distribution from QVIS
        CALL GAMQV
C
C------ relocate stagnation point
        CALL STMOVE
C
C------ set updated CL,CD
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
        CALL CDCALC
C
C------ display changes and test for convergence
        IF(RLX.LT.1.0) 
     &   WRITE(*,2000) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL,RLX
        IF(RLX.EQ.1.0) 
     &   WRITE(*,2010) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL
         CDPDIF = CD - CDF
         WRITE(*,2020) ALFA/DTOR, CL, CM, CD, CDF, CDPDIF
c         CDSURF = CDP + CDF
c         WRITE(*,2025) CDSURF, CDF, CDP

        IF(RMSBL .LT. EPS1) THEN
         LVCONV = .TRUE.
         AVISC = ALFA
         MVISC = MINF
         GO TO 90
        ENDIF
C
 1000 CONTINUE
      WRITE(*,*) 'VISCAL:  Convergence failed'
C
   90 CONTINUE
      CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
      CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
C      IF(LFLAP) CALL MHINGE


        is = 1
        hkmax = 0.
        hkm = 0.0
        psep = 0.
        patt = 0.
        do ibl = 2, iblte(is)
          hki = dstr(ibl,is) / thet(ibl,is)
          hkmax = max(hki,hkmax)
          if(hkm .lt. 4.0 .and. 
     &       hki .ge. 4.0      ) then
           hfrac = (4.0 - hkm) / (hki - hkm )
           pdefm = uedg(ibl-1,is)**2 * thet(ibl-1,is)
           pdefi = uedg(ibl  ,is)**2 * thet(ibl  ,is)
           psep = pdefm*(1.0-hfrac) + pdefi*hfrac
          endif
          if(hkm .gt. 4.0 .and. 
     &       hki .lt. 4.0      ) then
           hfrac = (4.0 - hkm) / (hki - hkm )
           pdefm = uedg(ibl-1,is)**2 * thet(ibl-1,is)
           pdefi = uedg(ibl  ,is)**2 * thet(ibl  ,is)
           patt = pdefm*(1.0-hfrac) + pdefi*hfrac
          endif
          hkm = hki
        enddo
        delp = patt - psep

        write(*,9922) 
     &    acrit(is), hkmax, cd, 2.0*psep, 2.0*patt, 2.0*delp,
     &    xoctr(is)
 9922   format(1x, f10.3, f10.4, f11.6, 3f11.6, f10.4, '     #')

      RETURN
C....................................................................
 2000   FORMAT
     &   (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3,
     &     '   RLX:',F6.3)
 2010   FORMAT
     &   (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3)
 2020   FORMAT
     &   ( 1X,3X,'   a =', F7.3,'      CL =',F8.4  /
     &     1X,3X,'  Cm =', F8.4, '     CD =',F9.5,
     &           '   =>   CDf =',F9.5,'    CDp =',F9.5)
 2025   FORMAT
     &   ( 1X,3X, 6X     ,  8X , ' Int CD =',F9.5,
     &           '   =>   CDf =',F9.5,'    CDp =',F9.5)
      END ! VISCAL
      
      
      SUBROUTINE MRSHOW(LM,LR)
      INCLUDE 'XFOIL.INC'
      LOGICAL LM, LR
C
      IF(LM .OR. LR) WRITE(*,*)
C
      IF(LM) THEN
       IF(MATYP.EQ.1) WRITE(*,1100) MINF1
       IF(MATYP.EQ.2) WRITE(*,1100) MINF1, ' / sqrt(CL)'
       IF(MATYP.EQ.3) WRITE(*,1100) MINF1, ' / CL'
      ENDIF
C
      IF(LR) THEN
       IF(RETYP.EQ.1) WRITE(*,1200) REINF1
       IF(RETYP.EQ.2) WRITE(*,1200) REINF1, ' / sqrt(CL)'
       IF(RETYP.EQ.3) WRITE(*,1200) REINF1, ' / CL'
      ENDIF
C
      RETURN
C
 1100 FORMAT(1X,'M  =' , F10.4, A)
 1200 FORMAT(1X,'Re =' , G12.4, A)
      END ! MRSHOW 
      
      
C***********************************************************************
C    Module:  xpanel.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************


      SUBROUTINE APCALC
      INCLUDE 'XFOIL.INC'
C
C---- set angles of airfoil panels
      DO 10 I=1, N-1
        SX = X(I+1) - X(I)
        SY = Y(I+1) - Y(I)
        IF(SX.EQ.0.0 .AND. SY.EQ.0.0) THEN
          APANEL(I) = ATAN2( -NY(I) , -NX(I) )
        ELSE
          APANEL(I) = ATAN2( SX , -SY )
        ENDIF
   10 CONTINUE
C
C---- TE panel
      I = N
      IP = 1
      IF(SHARP) THEN
       APANEL(I) = PI
      ELSE
       SX = X(IP) - X(I)
       SY = Y(IP) - Y(I)
       APANEL(I) = ATAN2( -SX , SY ) + PI
      ENDIF
C
      RETURN
      END
 
 
      SUBROUTINE NCALC(X,Y,S,N,XN,YN)
C---------------------------------------
C     Calculates normal unit vector
C     components at airfoil panel nodes
C---------------------------------------
      DIMENSION X(N), Y(N), S(N), XN(N), YN(N)
C
      IF(N.LE.1) RETURN
C
      CALL SEGSPL(X,XN,S,N)
      CALL SEGSPL(Y,YN,S,N)
      DO 10 I=1, N
        SX =  YN(I)
        SY = -XN(I)
        SMOD = SQRT(SX*SX + SY*SY)
        IF(SMOD .EQ. 0.0) THEN
         XN(I) = -1.0
         YN(I) = 0.0
        ELSE
         XN(I) = SX/SMOD
         YN(I) = SY/SMOD
        ENDIF
   10 CONTINUE
C
C---- average normal vectors at corner points
      DO 20 I=1, N-1
        IF(S(I) .EQ. S(I+1)) THEN
          SX = 0.5*(XN(I) + XN(I+1))
          SY = 0.5*(YN(I) + YN(I+1))
          SMOD = SQRT(SX*SX + SY*SY)
          IF(SMOD .EQ. 0.0) THEN
           XN(I) = -1.0
           YN(I) = 0.0
           XN(I+1) = -1.0
           YN(I+1) = 0.0
          ELSE
           XN(I)   = SX/SMOD
           YN(I)   = SY/SMOD
           XN(I+1) = SX/SMOD
           YN(I+1) = SY/SMOD
          ENDIF
        ENDIF
 20   CONTINUE
C
      RETURN
      END

 
      SUBROUTINE PSILIN(I,XI,YI,NXI,NYI,PSI,PSI_NI,GEOLIN,SIGLIN)
C-----------------------------------------------------------------------
C     Calculates current streamfunction Psi at panel node or wake node
C     I due to freestream and all bound vorticity Gam on the airfoil. 
C     Sensitivities of Psi with respect to alpha (Z_ALFA) and inverse
C     Qspec DOFs (Z_QDOF0,Z_QDOF1) which influence Gam in inverse cases.
C     Also calculates the sensitivity vector dPsi/dGam (DZDG).
C
C     If SIGLIN=True, then Psi includes the effects of the viscous
C     source distribution Sig and the sensitivity vector dPsi/dSig
C     (DZDM) is calculated.
C
C     If GEOLIN=True, then the geometric sensitivity vector dPsi/dn
C     is calculated, where n is the normal motion of the jth node.
C
C          Airfoil:  1   < I < N
C          Wake:     N+1 < I < N+NW
C-----------------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL NXO, NYO, NXP, NYP, NXI, NYI
      LOGICAL GEOLIN,SIGLIN
C
C---- distance tolerance for determining if two points are the same
      SEPS = (S(N)-S(1)) * 1.0E-5
C
      IO = I
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
      DO 3 JO=1, N
        DZDG(JO) = 0.0
        DZDN(JO) = 0.0
        DQDG(JO) = 0.0
    3 CONTINUE
C
      DO 4 JO=1, N
        DZDM(JO) = 0.0
        DQDM(JO) = 0.0
    4 CONTINUE
C
      Z_QINF = 0.
      Z_ALFA = 0.
      Z_QDOF0 = 0.
      Z_QDOF1 = 0.
      Z_QDOF2 = 0.
      Z_QDOF3 = 0.
C
      PSI    = 0.
      PSI_NI = 0.
C
      QTAN1 = 0.
      QTAN2 = 0.
      QTANM = 0.
C
      IF(SHARP) THEN
       SCS = 1.0
       SDS = 0.0
      ELSE
       SCS = ANTE/DSTE
       SDS = ASTE/DSTE
      ENDIF
C
      DO 10 JO=1, N
        JP = JO+1
C
        JM = JO-1
        JQ = JP+1
C
        IF(JO.EQ.1) THEN
         JM = JO
        ELSE IF(JO.EQ.N-1) THEN
         JQ = JP
        ELSE IF(JO.EQ.N) THEN
         JP = 1
         IF((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2 .LT. SEPS**2) GO TO 12
        ENDIF
C
        DSO = SQRT((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2)
C
C------ skip null panel
        IF(DSO .EQ. 0.0) GO TO 10
C
        DSIO = 1.0 / DSO
C
        APAN = APANEL(JO)
C
        RX1 = XI - X(JO)
        RY1 = YI - Y(JO)
        RX2 = XI - X(JP)
        RY2 = YI - Y(JP)
C
        SX = (X(JP) - X(JO)) * DSIO
        SY = (Y(JP) - Y(JO)) * DSIO
C
        X1 = SX*RX1 + SY*RY1
        X2 = SX*RX2 + SY*RY2
        YY = SX*RY1 - SY*RX1
C
        RS1 = RX1*RX1 + RY1*RY1
        RS2 = RX2*RX2 + RY2*RY2
C
C------ set reflection flag SGN to avoid branch problems with arctan
        IF(IO.GE.1 .AND. IO.LE.N) THEN
C------- no problem on airfoil surface
         SGN = 1.0
        ELSE
C------- make sure arctan falls between  -/+  Pi/2
         SGN = SIGN(1.0,YY)
        ENDIF
C
C------ set log(r^2) and arctan(x/y), correcting for reflection if any
        IF(IO.NE.JO .AND. RS1.GT.0.0) THEN
         G1 = LOG(RS1)
         T1 = ATAN2(SGN*X1,SGN*YY) + (0.5 - 0.5*SGN)*PI
        ELSE
         G1 = 0.0
         T1 = 0.0
        ENDIF
C
        IF(IO.NE.JP .AND. RS2.GT.0.0) THEN
         G2 = LOG(RS2)
         T2 = ATAN2(SGN*X2,SGN*YY) + (0.5 - 0.5*SGN)*PI
        ELSE
         G2 = 0.0
         T2 = 0.0
        ENDIF
C
        X1I = SX*NXI + SY*NYI
        X2I = SX*NXI + SY*NYI
        YYI = SX*NYI - SY*NXI
C
        IF(GEOLIN) THEN
         NXO = NX(JO)
         NYO = NY(JO)
         NXP = NX(JP)
         NYP = NY(JP)
C
         X1O =-((RX1-X1*SX)*NXO + (RY1-X1*SY)*NYO)*DSIO-(SX*NXO+SY*NYO)
         X1P = ((RX1-X1*SX)*NXP + (RY1-X1*SY)*NYP)*DSIO
         X2O =-((RX2-X2*SX)*NXO + (RY2-X2*SY)*NYO)*DSIO
         X2P = ((RX2-X2*SX)*NXP + (RY2-X2*SY)*NYP)*DSIO-(SX*NXP+SY*NYP)
         YYO = ((RX1+X1*SY)*NYO - (RY1-X1*SX)*NXO)*DSIO-(SX*NYO-SY*NXO)
         YYP =-((RX1-X1*SY)*NYP - (RY1+X1*SX)*NXP)*DSIO
        ENDIF
C
        IF(JO.EQ.N) GO TO 11
C
        IF(SIGLIN) THEN
C
C------- set up midpoint quantities
         X0 = 0.5*(X1+X2)
         RS0 = X0*X0 + YY*YY
         G0 = LOG(RS0)
         T0 = ATAN2(SGN*X0,SGN*YY) + (0.5 - 0.5*SGN)*PI
C
C------- calculate source contribution to Psi  for  1-0  half-panel
         DXINV = 1.0/(X1-X0)
         PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0)
         PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)
     &        + (X0-X1)*YY) * DXINV
C
         PSX1 =  -(T1-APAN)
         PSX0 =    T0-APAN
         PSYY =  0.5*(G1-G0)
C
         PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV
         PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV
         PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV
C
         DSM = SQRT((X(JP)-X(JM))**2 + (Y(JP)-Y(JM))**2)
         DSIM = 1.0/DSM
C
CCC      SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC      SIG1 = (SIG(JP) - SIG(JM))*DSIM
CCC      SSUM = SIG0 + SIG1
CCC      SDIF = SIG0 - SIG1
C
         SSUM = (SIG(JP) - SIG(JO))*DSIO + (SIG(JP) - SIG(JM))*DSIM
         SDIF = (SIG(JP) - SIG(JO))*DSIO - (SIG(JP) - SIG(JM))*DSIM
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JM) = DZDM(JM) + QOPI*(-PSUM*DSIM + PDIF*DSIM)
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*DSIO - PDIF*DSIO)
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*(DSIO+DSIM)
     &                                          + PDIF*(DSIO-DSIM))
C
C------- dPsi/dni
         PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI
         PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JM) = DQDM(JM) + QOPI*(-PSNI*DSIM + PDNI*DSIM)
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*DSIO - PDNI*DSIO)
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*(DSIO+DSIM)
     &                                          + PDNI*(DSIO-DSIM))
C
C
C------- calculate source contribution to Psi  for  0-2  half-panel
         DXINV = 1.0/(X0-X2)
         PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2)
         PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN)
     &        + (X2-X0)*YY) * DXINV
C
         PSX0 =  -(T0-APAN)
         PSX2 =    T2-APAN
         PSYY =  0.5*(G0-G2)
C
         PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV
         PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV
         PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV
C
         DSP = SQRT((X(JQ)-X(JO))**2 + (Y(JQ)-Y(JO))**2)
         DSIP = 1.0/DSP
C
CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SSUM = SIG2 + SIG0
CCC         SDIF = SIG2 - SIG0
C
         SSUM = (SIG(JQ) - SIG(JO))*DSIP + (SIG(JP) - SIG(JO))*DSIO
         SDIF = (SIG(JQ) - SIG(JO))*DSIP - (SIG(JP) - SIG(JO))*DSIO
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*(DSIP+DSIO)
     &                                          - PDIF*(DSIP-DSIO))
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*DSIO - PDIF*DSIO)
         DZDM(JQ) = DZDM(JQ) + QOPI*( PSUM*DSIP + PDIF*DSIP)
C
C------- dPsi/dni
         PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI
         PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)
     &                                          - PDNI*(DSIP-DSIO))
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO)
         DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP)
C
        ENDIF
C
C------ calculate vortex panel contribution to Psi
        DXINV = 1.0/(X1-X2)
        PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
        PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2-RS1*G1 + X1*X1-X2*X2))*DXINV
C
        PSX1 = 0.5*G1
        PSX2 = -.5*G2
        PSYY = T1-T2
C
        PDX1 = ((X1+X2)*PSX1 + PSIS - X1*G1 - PSID)*DXINV
        PDX2 = ((X1+X2)*PSX2 + PSIS + X2*G2 + PSID)*DXINV
        PDYY = ((X1+X2)*PSYY - YY*(G1-G2)         )*DXINV
C
        GSUM1 = GAMU(JP,1) + GAMU(JO,1)
        GSUM2 = GAMU(JP,2) + GAMU(JO,2)
        GDIF1 = GAMU(JP,1) - GAMU(JO,1)
        GDIF2 = GAMU(JP,2) - GAMU(JO,2)
C
        GSUM = GAM(JP) + GAM(JO)
        GDIF = GAM(JP) - GAM(JO)
C
        PSI = PSI + QOPI*(PSIS*GSUM + PSID*GDIF)
C
C------ dPsi/dGam
        DZDG(JO) = DZDG(JO) + QOPI*(PSIS-PSID)
        DZDG(JP) = DZDG(JP) + QOPI*(PSIS+PSID)
C
C------ dPsi/dni
        PSNI = PSX1*X1I + PSX2*X2I + PSYY*YYI
        PDNI = PDX1*X1I + PDX2*X2I + PDYY*YYI
        PSI_NI = PSI_NI + QOPI*(GSUM*PSNI + GDIF*PDNI)
C
        QTAN1 = QTAN1 + QOPI*(GSUM1*PSNI + GDIF1*PDNI)
        QTAN2 = QTAN2 + QOPI*(GSUM2*PSNI + GDIF2*PDNI)
C
        DQDG(JO) = DQDG(JO) + QOPI*(PSNI - PDNI)
        DQDG(JP) = DQDG(JP) + QOPI*(PSNI + PDNI)
C
        IF(GEOLIN) THEN
C
C------- dPsi/dn
         DZDN(JO) = DZDN(JO)+ QOPI*GSUM*(PSX1*X1O + PSX2*X2O + PSYY*YYO)
     &                      + QOPI*GDIF*(PDX1*X1O + PDX2*X2O + PDYY*YYO)
         DZDN(JP) = DZDN(JP)+ QOPI*GSUM*(PSX1*X1P + PSX2*X2P + PSYY*YYP)
     &                      + QOPI*GDIF*(PDX1*X1P + PDX2*X2P + PDYY*YYP)
C------- dPsi/dP
         Z_QDOF0 = Z_QDOF0
     &           + QOPI*((PSIS-PSID)*QF0(JO) + (PSIS+PSID)*QF0(JP))
         Z_QDOF1 = Z_QDOF1
     &           + QOPI*((PSIS-PSID)*QF1(JO) + (PSIS+PSID)*QF1(JP))
         Z_QDOF2 = Z_QDOF2
     &           + QOPI*((PSIS-PSID)*QF2(JO) + (PSIS+PSID)*QF2(JP))
         Z_QDOF3 = Z_QDOF3
     &           + QOPI*((PSIS-PSID)*QF3(JO) + (PSIS+PSID)*QF3(JP))
        ENDIF
C
C
   10 CONTINUE
C
   11 CONTINUE
      PSIG = 0.5*YY*(G1-G2) + X2*(T2-APAN) - X1*(T1-APAN)
      PGAM = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
C
      PSIGX1 = -(T1-APAN)
      PSIGX2 =   T2-APAN
      PSIGYY = 0.5*(G1-G2)
      PGAMX1 = 0.5*G1
      PGAMX2 = -.5*G2
      PGAMYY = T1-T2
C
      PSIGNI = PSIGX1*X1I + PSIGX2*X2I + PSIGYY*YYI
      PGAMNI = PGAMX1*X1I + PGAMX2*X2I + PGAMYY*YYI
C
C---- TE panel source and vortex strengths
      SIGTE1 = 0.5*SCS*(GAMU(JP,1) - GAMU(JO,1))
      SIGTE2 = 0.5*SCS*(GAMU(JP,2) - GAMU(JO,2))
      GAMTE1 = -.5*SDS*(GAMU(JP,1) - GAMU(JO,1))
      GAMTE2 = -.5*SDS*(GAMU(JP,2) - GAMU(JO,2))
C
      SIGTE = 0.5*SCS*(GAM(JP) - GAM(JO))
      GAMTE = -.5*SDS*(GAM(JP) - GAM(JO))
C
C---- TE panel contribution to Psi
      PSI = PSI + HOPI*(PSIG*SIGTE + PGAM*GAMTE)
C
C---- dPsi/dGam
      DZDG(JO) = DZDG(JO) - HOPI*PSIG*SCS*0.5
      DZDG(JP) = DZDG(JP) + HOPI*PSIG*SCS*0.5
C
      DZDG(JO) = DZDG(JO) + HOPI*PGAM*SDS*0.5
      DZDG(JP) = DZDG(JP) - HOPI*PGAM*SDS*0.5
C
C---- dPsi/dni
      PSI_NI = PSI_NI + HOPI*(PSIGNI*SIGTE + PGAMNI*GAMTE)
C
      QTAN1 = QTAN1 + HOPI*(PSIGNI*SIGTE1 + PGAMNI*GAMTE1)
      QTAN2 = QTAN2 + HOPI*(PSIGNI*SIGTE2 + PGAMNI*GAMTE2)
C
      DQDG(JO) = DQDG(JO) - HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS)
      DQDG(JP) = DQDG(JP) + HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS)
C
      IF(GEOLIN) THEN
C
C----- dPsi/dn
       DZDN(JO) = DZDN(JO)
     &          + HOPI*(PSIGX1*X1O + PSIGX2*X2O + PSIGYY*YYO)*SIGTE
     &          + HOPI*(PGAMX1*X1O + PGAMX2*X2O + PGAMYY*YYO)*GAMTE
       DZDN(JP) = DZDN(JP)
     &          + HOPI*(PSIGX1*X1P + PSIGX2*X2P + PSIGYY*YYP)*SIGTE
     &          + HOPI*(PGAMX1*X1P + PGAMX2*X2P + PGAMYY*YYP)*GAMTE
C
C----- dPsi/dP
       Z_QDOF0 = Z_QDOF0 + HOPI*PSIG*0.5*(QF0(JP)-QF0(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF0(JP)-QF0(JO))*SDS
       Z_QDOF1 = Z_QDOF1 + HOPI*PSIG*0.5*(QF1(JP)-QF1(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF1(JP)-QF1(JO))*SDS
       Z_QDOF2 = Z_QDOF2 + HOPI*PSIG*0.5*(QF2(JP)-QF2(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF2(JP)-QF2(JO))*SDS
       Z_QDOF3 = Z_QDOF3 + HOPI*PSIG*0.5*(QF3(JP)-QF3(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF3(JP)-QF3(JO))*SDS
C
      ENDIF
C
   12 CONTINUE
C
C**** Freestream terms
      PSI = PSI + QINF*(COSA*YI - SINA*XI)
C
C---- dPsi/dn
      PSI_NI = PSI_NI + QINF*(COSA*NYI - SINA*NXI)
C
      QTAN1 = QTAN1 + QINF*NYI
      QTAN2 = QTAN2 - QINF*NXI
C
C---- dPsi/dQinf
      Z_QINF = Z_QINF + (COSA*YI - SINA*XI)
C
C---- dPsi/dalfa
      Z_ALFA = Z_ALFA - QINF*(SINA*YI + COSA*XI)
C
      IF(.NOT.LIMAGE) RETURN
C
C
C
      DO 20 JO=1, N
        JP = JO+1
C
        JM = JO-1
        JQ = JP+1
C
        IF(JO.EQ.1) THEN
         JM = JO
        ELSE IF(JO.EQ.N-1) THEN
         JQ = JP
        ELSE IF(JO.EQ.N) THEN
         JP = 1
         IF((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2 .LT. SEPS**2) GO TO 22
        ENDIF
C
        DSO = SQRT((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2)
C
C------ skip null panel
        IF(DSO .EQ. 0.0) GO TO 20
C
        DSIO = 1.0 / DSO
C
ccc     APAN = APANEL(JO)
        APAN = PI - APANEL(JO) + 2.0*ALFA
C
        XJO = X(JO) + 2.0*(YIMAGE+Y(JO))*SINA
        YJO = Y(JO) - 2.0*(YIMAGE+Y(JO))*COSA
        XJP = X(JP) + 2.0*(YIMAGE+Y(JP))*SINA
        YJP = Y(JP) - 2.0*(YIMAGE+Y(JP))*COSA
C
        RX1 = XI - XJO
        RY1 = YI - YJO
        RX2 = XI - XJP
        RY2 = YI - YJP
C
        SX = (XJP - XJO) * DSIO
        SY = (YJP - YJO) * DSIO
C
        X1 = SX*RX1 + SY*RY1
        X2 = SX*RX2 + SY*RY2
        YY = SX*RY1 - SY*RX1
C
        RS1 = RX1*RX1 + RY1*RY1
        RS2 = RX2*RX2 + RY2*RY2
C
C------ set reflection flag SGN to avoid branch problems with arctan
        IF(IO.GE.1 .AND. IO.LE.N) THEN
C------- no problem on airfoil surface
         SGN = 1.0
        ELSE
C------- make sure arctan falls between  -/+  Pi/2
         SGN = SIGN(1.0,YY)
        ENDIF
C
C------ set log(r^2) and arctan(x/y), correcting for reflection if any
        G1 = LOG(RS1)
        T1 = ATAN2(SGN*X1,SGN*YY) + (0.5 - 0.5*SGN)*PI
C
        G2 = LOG(RS2)
        T2 = ATAN2(SGN*X2,SGN*YY) + (0.5 - 0.5*SGN)*PI
C
        X1I = SX*NXI + SY*NYI
        X2I = SX*NXI + SY*NYI
        YYI = SX*NYI - SY*NXI
C
        IF(GEOLIN) THEN
         NXO = NX(JO)
         NYO = NY(JO)
         NXP = NX(JP)
         NYP = NY(JP)
C
         X1O =-((RX1-X1*SX)*NXO + (RY1-X1*SY)*NYO)*DSIO-(SX*NXO+SY*NYO)
         X1P = ((RX1-X1*SX)*NXP + (RY1-X1*SY)*NYP)*DSIO
         X2O =-((RX2-X2*SX)*NXO + (RY2-X2*SY)*NYO)*DSIO
         X2P = ((RX2-X2*SX)*NXP + (RY2-X2*SY)*NYP)*DSIO-(SX*NXP+SY*NYP)
         YYO = ((RX1+X1*SY)*NYO - (RY1-X1*SX)*NXO)*DSIO-(SX*NYO-SY*NXO)
         YYP =-((RX1-X1*SY)*NYP - (RY1+X1*SX)*NXP)*DSIO
        ENDIF
C
        IF(JO.EQ.N) GO TO 21
C
        IF(SIGLIN) THEN
C
C------- set up midpoint quantities
         X0 = 0.5*(X1+X2)
         RS0 = X0*X0 + YY*YY
         G0 = LOG(RS0)
         T0 = ATAN2(SGN*X0,SGN*YY) + (0.5 - 0.5*SGN)*PI
C
C------- calculate source contribution to Psi  for  1-0  half-panel
         DXINV = 1.0/(X1-X0)
         PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0)
         PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)
     &        + (X0-X1)*YY) * DXINV
C
         PSX1 =  -(T1-APAN)
         PSX0 =    T0-APAN
         PSYY =  0.5*(G1-G0)
C
         PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV
         PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV
         PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV
C
         DSM = SQRT((X(JP)-X(JM))**2 + (Y(JP)-Y(JM))**2)
         DSIM = 1.0/DSM
C
CCC      SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC      SIG1 = (SIG(JP) - SIG(JM))*DSIM
CCC      SSUM = SIG0 + SIG1
CCC      SDIF = SIG0 - SIG1
C
         SSUM = (SIG(JP) - SIG(JO))*DSIO + (SIG(JP) - SIG(JM))*DSIM
         SDIF = (SIG(JP) - SIG(JO))*DSIO - (SIG(JP) - SIG(JM))*DSIM
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JM) = DZDM(JM) + QOPI*(-PSUM*DSIM + PDIF*DSIM)
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*DSIO - PDIF*DSIO)
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*(DSIO+DSIM)
     &                                          + PDIF*(DSIO-DSIM))
C
C------- dPsi/dni
         PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI
         PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JM) = DQDM(JM) + QOPI*(-PSNI*DSIM + PDNI*DSIM)
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*DSIO - PDNI*DSIO)
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*(DSIO+DSIM)
     &                                          + PDNI*(DSIO-DSIM))
C
C
C------- calculate source contribution to Psi  for  0-2  half-panel
         DXINV = 1.0/(X0-X2)
         PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2)
         PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN)
     &        + (X2-X0)*YY) * DXINV
C
         PSX0 =  -(T0-APAN)
         PSX2 =    T2-APAN
         PSYY =  0.5*(G0-G2)
C
         PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV
         PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV
         PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV
C
         DSP = SQRT((X(JQ)-X(JO))**2 + (Y(JQ)-Y(JO))**2)
         DSIP = 1.0/DSP
C
CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SSUM = SIG2 + SIG0
CCC         SDIF = SIG2 - SIG0
C
         SSUM = (SIG(JQ) - SIG(JO))*DSIP + (SIG(JP) - SIG(JO))*DSIO
         SDIF = (SIG(JQ) - SIG(JO))*DSIP - (SIG(JP) - SIG(JO))*DSIO
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*(DSIP+DSIO)
     &                                          - PDIF*(DSIP-DSIO))
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*DSIO - PDIF*DSIO)
         DZDM(JQ) = DZDM(JQ) + QOPI*( PSUM*DSIP + PDIF*DSIP)
C
C------- dPsi/dni
         PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI
         PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)
     &                                          - PDNI*(DSIP-DSIO))
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO)
         DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP)
C
        ENDIF
C
C------ calculate vortex panel contribution to Psi
        DXINV = 1.0/(X1-X2)
        PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
        PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2-RS1*G1 + X1*X1-X2*X2))*DXINV
C
        PSX1 = 0.5*G1
        PSX2 = -.5*G2
        PSYY = T1-T2
C
        PDX1 = ((X1+X2)*PSX1 + PSIS - X1*G1 - PSID)*DXINV
        PDX2 = ((X1+X2)*PSX2 + PSIS + X2*G2 + PSID)*DXINV
        PDYY = ((X1+X2)*PSYY - YY*(G1-G2)         )*DXINV
C
        GSUM1 = GAMU(JP,1) + GAMU(JO,1)
        GSUM2 = GAMU(JP,2) + GAMU(JO,2)
        GDIF1 = GAMU(JP,1) - GAMU(JO,1)
        GDIF2 = GAMU(JP,2) - GAMU(JO,2)
C
        GSUM = GAM(JP) + GAM(JO)
        GDIF = GAM(JP) - GAM(JO)
C
        PSI = PSI - QOPI*(PSIS*GSUM + PSID*GDIF)
C
C------ dPsi/dGam
        DZDG(JO) = DZDG(JO) - QOPI*(PSIS-PSID)
        DZDG(JP) = DZDG(JP) - QOPI*(PSIS+PSID)
C
C------ dPsi/dni
        PSNI = PSX1*X1I + PSX2*X2I + PSYY*YYI
        PDNI = PDX1*X1I + PDX2*X2I + PDYY*YYI
        PSI_NI = PSI_NI - QOPI*(GSUM*PSNI + GDIF*PDNI)
C
        QTAN1 = QTAN1 - QOPI*(GSUM1*PSNI + GDIF1*PDNI)
        QTAN2 = QTAN2 - QOPI*(GSUM2*PSNI + GDIF2*PDNI)
C
        DQDG(JO) = DQDG(JO) - QOPI*(PSNI - PDNI)
        DQDG(JP) = DQDG(JP) - QOPI*(PSNI + PDNI)
C
        IF(GEOLIN) THEN
C
C------- dPsi/dn
         DZDN(JO) = DZDN(JO)- QOPI*GSUM*(PSX1*X1O + PSX2*X2O + PSYY*YYO)
     &                      - QOPI*GDIF*(PDX1*X1O + PDX2*X2O + PDYY*YYO)
         DZDN(JP) = DZDN(JP)- QOPI*GSUM*(PSX1*X1P + PSX2*X2P + PSYY*YYP)
     &                      - QOPI*GDIF*(PDX1*X1P + PDX2*X2P + PDYY*YYP)
C------- dPsi/dP
         Z_QDOF0 = Z_QDOF0
     &           - QOPI*((PSIS-PSID)*QF0(JO) + (PSIS+PSID)*QF0(JP))
         Z_QDOF1 = Z_QDOF1
     &           - QOPI*((PSIS-PSID)*QF1(JO) + (PSIS+PSID)*QF1(JP))
         Z_QDOF2 = Z_QDOF2
     &           - QOPI*((PSIS-PSID)*QF2(JO) + (PSIS+PSID)*QF2(JP))
         Z_QDOF3 = Z_QDOF3
     &           - QOPI*((PSIS-PSID)*QF3(JO) + (PSIS+PSID)*QF3(JP))
        ENDIF
C
C
   20 CONTINUE
C
   21 CONTINUE
      PSIG = 0.5*YY*(G1-G2) + X2*(T2-APAN) - X1*(T1-APAN)
      PGAM = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
C
      PSIGX1 = -(T1-APAN)
      PSIGX2 =   T2-APAN
      PSIGYY = 0.5*(G1-G2)
      PGAMX1 = 0.5*G1
      PGAMX2 = -.5*G2
      PGAMYY = T1-T2
C
      PSIGNI = PSIGX1*X1I + PSIGX2*X2I + PSIGYY*YYI
      PGAMNI = PGAMX1*X1I + PGAMX2*X2I + PGAMYY*YYI
C
C---- TE panel source and vortex strengths
      SIGTE1 = 0.5*SCS*(GAMU(JP,1) - GAMU(JO,1))
      SIGTE2 = 0.5*SCS*(GAMU(JP,2) - GAMU(JO,2))
      GAMTE1 = -.5*SDS*(GAMU(JP,1) - GAMU(JO,1))
      GAMTE2 = -.5*SDS*(GAMU(JP,2) - GAMU(JO,2))
C
      SIGTE = 0.5*SCS*(GAM(JP) - GAM(JO))
      GAMTE = -.5*SDS*(GAM(JP) - GAM(JO))
C
C---- TE panel contribution to Psi
      PSI = PSI + HOPI*(PSIG*SIGTE - PGAM*GAMTE)
C
C---- dPsi/dGam
      DZDG(JO) = DZDG(JO) - HOPI*PSIG*SCS*0.5
      DZDG(JP) = DZDG(JP) + HOPI*PSIG*SCS*0.5
C
      DZDG(JO) = DZDG(JO) - HOPI*PGAM*SDS*0.5
      DZDG(JP) = DZDG(JP) + HOPI*PGAM*SDS*0.5
C
C---- dPsi/dni
      PSI_NI = PSI_NI + HOPI*(PSIGNI*SIGTE - PGAMNI*GAMTE)
C
      QTAN1 = QTAN1 + HOPI*(PSIGNI*SIGTE1 - PGAMNI*GAMTE1)
      QTAN2 = QTAN2 + HOPI*(PSIGNI*SIGTE2 - PGAMNI*GAMTE2)
C
      DQDG(JO) = DQDG(JO) - HOPI*(PSIGNI*0.5*SCS + PGAMNI*0.5*SDS)
      DQDG(JP) = DQDG(JP) + HOPI*(PSIGNI*0.5*SCS + PGAMNI*0.5*SDS)
C
      IF(GEOLIN) THEN
C
C----- dPsi/dn
       DZDN(JO) = DZDN(JO)
     &          + HOPI*(PSIGX1*X1O + PSIGX2*X2O + PSIGYY*YYO)*SIGTE
     &          - HOPI*(PGAMX1*X1O + PGAMX2*X2O + PGAMYY*YYO)*GAMTE
       DZDN(JP) = DZDN(JP)
     &          + HOPI*(PSIGX1*X1P + PSIGX2*X2P + PSIGYY*YYP)*SIGTE
     &          - HOPI*(PGAMX1*X1P + PGAMX2*X2P + PGAMYY*YYP)*GAMTE
C
C----- dPsi/dP
       Z_QDOF0 = Z_QDOF0 + HOPI*PSIG*0.5*(QF0(JP)-QF0(JO))*SCS
     &                   + HOPI*PGAM*0.5*(QF0(JP)-QF0(JO))*SDS
       Z_QDOF1 = Z_QDOF1 + HOPI*PSIG*0.5*(QF1(JP)-QF1(JO))*SCS
     &                   + HOPI*PGAM*0.5*(QF1(JP)-QF1(JO))*SDS
       Z_QDOF2 = Z_QDOF2 + HOPI*PSIG*0.5*(QF2(JP)-QF2(JO))*SCS
     &                   + HOPI*PGAM*0.5*(QF2(JP)-QF2(JO))*SDS
       Z_QDOF3 = Z_QDOF3 + HOPI*PSIG*0.5*(QF3(JP)-QF3(JO))*SCS
     &                   + HOPI*PGAM*0.5*(QF3(JP)-QF3(JO))*SDS
C
      ENDIF
C
   22 CONTINUE
C
      RETURN
      END


      SUBROUTINE PSWLIN(I,XI,YI,NXI,NYI,PSI,PSI_NI)
C--------------------------------------------------------------------
C     Calculates current streamfunction Psi and tangential velocity
C     Qtan at panel node or wake node I due to freestream and wake
C     sources Sig.  Also calculates sensitivity vectors dPsi/dSig
C     (DZDM) and dQtan/dSig (DQDM).
C
C          Airfoil:  1   < I < N
C          Wake:     N+1 < I < N+NW
C--------------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL NXI, NYI
C
      IO = I
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
      DO 4 JO=N+1, N+NW
        DZDM(JO) = 0.0
        DQDM(JO) = 0.0
    4 CONTINUE
C
      PSI    = 0.
      PSI_NI = 0.
C
      DO 20 JO=N+1, N+NW-1
C
        JP = JO+1
C
        JM = JO-1
        JQ = JP+1
        IF(JO.EQ.N+1) THEN
         JM = JO
        ELSE IF(JO.EQ.N+NW-1) THEN
         JQ = JP
        ENDIF
C
        DSO = SQRT((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2)
        DSIO = 1.0 / DSO
C
        APAN = APANEL(JO)
C
        RX1 = XI - X(JO)
        RY1 = YI - Y(JO)
        RX2 = XI - X(JP)
        RY2 = YI - Y(JP)
C
        SX = (X(JP) - X(JO)) * DSIO
        SY = (Y(JP) - Y(JO)) * DSIO
C
        X1 = SX*RX1 + SY*RY1
        X2 = SX*RX2 + SY*RY2
        YY = SX*RY1 - SY*RX1
C
        RS1 = RX1*RX1 + RY1*RY1
        RS2 = RX2*RX2 + RY2*RY2
C
        IF(IO.GE.N+1 .AND. IO.LE.N+NW) THEN
         SGN = 1.0
        ELSE
         SGN = SIGN(1.0,YY)
        ENDIF
C
        IF(IO.NE.JO .AND. RS1.GT.0.0) THEN
         G1 = LOG(RS1)
         T1 = ATAN2(SGN*X1,SGN*YY) - (0.5 - 0.5*SGN)*PI
        ELSE
         G1 = 0.0
         T1 = 0.0
        ENDIF
C
        IF(IO.NE.JP .AND. RS2.GT.0.0) THEN
         G2 = LOG(RS2)
         T2 = ATAN2(SGN*X2,SGN*YY) - (0.5 - 0.5*SGN)*PI
        ELSE
         G2 = 0.0
         T2 = 0.0
        ENDIF
C
        X1I = SX*NXI + SY*NYI
        X2I = SX*NXI + SY*NYI
        YYI = SX*NYI - SY*NXI
C
C------- set up midpoint quantities
         X0 = 0.5*(X1+X2)
         RS0 = X0*X0 + YY*YY
         G0 = LOG(RS0)
         T0 = ATAN2(SGN*X0,SGN*YY) - (0.5 - 0.5*SGN)*PI
C
C------- calculate source contribution to Psi  for  1-0  half-panel
         DXINV = 1.0/(X1-X0)
         PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0)
         PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)
     &        + (X0-X1)*YY) * DXINV
C
         PSX1 =  -(T1-APAN)
         PSX0 =    T0-APAN
         PSYY =  0.5*(G1-G0)
C
         PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV
         PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV
         PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV
C
         DSM = SQRT((X(JP)-X(JM))**2 + (Y(JP)-Y(JM))**2)
         DSIM = 1.0/DSM
C
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SIG1 = (SIG(JP) - SIG(JM))*DSIM
CCC         SSUM = SIG0 + SIG1
CCC         SDIF = SIG0 - SIG1
C
         SSUM = (SIG(JP) - SIG(JO))*DSIO + (SIG(JP) - SIG(JM))*DSIM
         SDIF = (SIG(JP) - SIG(JO))*DSIO - (SIG(JP) - SIG(JM))*DSIM
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JM) = DZDM(JM) + QOPI*(-PSUM*DSIM + PDIF*DSIM)
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*DSIO - PDIF*DSIO)
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*(DSIO+DSIM)
     &                                          + PDIF*(DSIO-DSIM))
C
C------- dPsi/dni
         PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI
         PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JM) = DQDM(JM) + QOPI*(-PSNI*DSIM + PDNI*DSIM)
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*DSIO - PDNI*DSIO)
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*(DSIO+DSIM)
     &                                          + PDNI*(DSIO-DSIM))
C
C
C------- calculate source contribution to Psi  for  0-2  half-panel
         DXINV = 1.0/(X0-X2)
         PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2)
         PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN)
     &        + (X2-X0)*YY) * DXINV
C
         PSX0 =  -(T0-APAN)
         PSX2 =    T2-APAN
         PSYY =  0.5*(G0-G2)
C
         PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV
         PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV
         PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV
C
         DSP = SQRT((X(JQ)-X(JO))**2 + (Y(JQ)-Y(JO))**2)
         DSIP = 1.0/DSP
C
CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SSUM = SIG2 + SIG0
CCC         SDIF = SIG2 - SIG0
C
         SSUM = (SIG(JQ) - SIG(JO))*DSIP + (SIG(JP) - SIG(JO))*DSIO
         SDIF = (SIG(JQ) - SIG(JO))*DSIP - (SIG(JP) - SIG(JO))*DSIO
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*(DSIP+DSIO)
     &                                          - PDIF*(DSIP-DSIO))
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*DSIO - PDIF*DSIO)
         DZDM(JQ) = DZDM(JQ) + QOPI*( PSUM*DSIP + PDIF*DSIP)
C
C------- dPsi/dni
         PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI
         PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)
     &                                          - PDNI*(DSIP-DSIO))
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO)
         DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP)
C
   20 CONTINUE
C
      RETURN
      END




      SUBROUTINE GGCALC
C--------------------------------------------------------------
C     Calculates two surface vorticity (gamma) distributions
C     for alpha = 0, 90  degrees.  These are superimposed
C     in SPECAL or SPECCL for specified alpha or CL.
C--------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- distance of internal control point ahead of sharp TE
C-    (fraction of smaller panel length adjacent to TE)
      BWT = 0.1
C
      WRITE(*,*) 'Calculating unit vorticity distributions ...'
C
      DO 10 I=1, N
        GAM(I) = 0.
        GAMU(I,1) = 0.
        GAMU(I,2) = 0.
   10 CONTINUE
      PSIO = 0.
C
C---- Set up matrix system for  Psi = Psio  on airfoil surface.
C-    The unknowns are (dGamma)i and dPsio.
      DO 20 I=1, N
C
C------ calculate Psi and dPsi/dGamma array for current node
        CALL PSILIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N,.FALSE.,.TRUE.)
C
        PSIINF = QINF*(COS(ALFA)*Y(I) - SIN(ALFA)*X(I))
C
C------ RES1 = PSI( 0) - PSIO
C------ RES2 = PSI(90) - PSIO
        RES1 =  QINF*Y(I)
        RES2 = -QINF*X(I)
C
C------ dRes/dGamma
        DO 201 J=1, N
          AIJ(I,J) = DZDG(J)
  201   CONTINUE
C
        DO 202 J=1, N
          BIJ(I,J) = -DZDM(J)
  202   CONTINUE
C
C------ dRes/dPsio
        AIJ(I,N+1) = -1.0
C
        GAMU(I,1) = -RES1
        GAMU(I,2) = -RES2
C
   20 CONTINUE
C
C---- set Kutta condition
C-    RES = GAM(1) + GAM(N)
      RES = 0.
C
      DO 30 J=1, N+1
        AIJ(N+1,J) = 0.0
   30 CONTINUE
C
      AIJ(N+1,1) = 1.0
      AIJ(N+1,N) = 1.0
C
      GAMU(N+1,1) = -RES
      GAMU(N+1,2) = -RES
C
C---- set up Kutta condition (no direct source influence)
      DO 32 J=1, N
        BIJ(N+1,J) = 0.
   32 CONTINUE
C
      IF(SHARP) THEN
C----- set zero internal velocity in TE corner 
C
C----- set TE bisector angle
       AG1 = ATAN2(-YP(1),-XP(1)    )
       AG2 = ATANC( YP(N), XP(N),AG1)
       ABIS = 0.5*(AG1+AG2)
       CBIS = COS(ABIS)
       SBIS = SIN(ABIS)
C
C----- minimum panel length adjacent to TE
       DS1 = SQRT( (X(1)-X(2)  )**2 + (Y(1)-Y(2)  )**2 )
       DS2 = SQRT( (X(N)-X(N-1))**2 + (Y(N)-Y(N-1))**2 )
       DSMIN = MIN( DS1 , DS2 )
C
C----- control point on bisector just ahead of TE point
       XBIS = XTE - BWT*DSMIN*CBIS
       YBIS = YTE - BWT*DSMIN*SBIS
ccc       write(*,*) xbis, ybis
C
C----- set velocity component along bisector line
       CALL PSILIN(0,XBIS,YBIS,-SBIS,CBIS,PSI,QBIS,.FALSE.,.TRUE.)
C
CCC--- RES = DQDGj*Gammaj + DQDMj*Massj + QINF*(COSA*CBIS + SINA*SBIS)
       RES = QBIS
C
C----- dRes/dGamma
       DO J=1, N
         AIJ(N,J) = DQDG(J)
       ENDDO
C
C----- -dRes/dMass
       DO J=1, N
         BIJ(N,J) = -DQDM(J)
       ENDDO
C
C----- dRes/dPsio
       AIJ(N,N+1) = 0.
C
C----- -dRes/dUinf
       GAMU(N,1) = -CBIS
C
C----- -dRes/dVinf
       GAMU(N,2) = -SBIS
C
      ENDIF
C
C---- LU-factor coefficient matrix AIJ
      CALL LUDCMP(IQX,N+1,AIJ,AIJPIV)
      LQAIJ = .TRUE.
C
C---- solve system for the two vorticity distributions
      CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,GAMU(1,1))
      CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,GAMU(1,2))
C
C---- set inviscid alpha=0,90 surface speeds for this geometry
      DO 50 I=1, N
        QINVU(I,1) = GAMU(I,1)
        QINVU(I,2) = GAMU(I,2)
   50 CONTINUE
C
      LGAMU = .TRUE.
C
      RETURN
      END



      SUBROUTINE QWCALC
C---------------------------------------------------------------
C     Sets inviscid tangential velocity for alpha = 0, 90
C     on wake due to freestream and airfoil surface vorticity.
C---------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- first wake point (same as TE)
      QINVU(N+1,1) = QINVU(N,1)
      QINVU(N+1,2) = QINVU(N,2)
C
C---- rest of wake
      DO 10 I=N+2, N+NW
        CALL PSILIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_NI,.FALSE.,.FALSE.)
        QINVU(I,1) = QTAN1
        QINVU(I,2) = QTAN2
   10 CONTINUE
C
      RETURN
      END


      SUBROUTINE QDCALC
C-----------------------------------------------------
C     Calculates source panel influence coefficient
C     matrix for current airfoil and wake geometry.
C-----------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      WRITE(*,*) 'Calculating source influence matrix ...'
C
      IF(.NOT.LADIJ) THEN
C
C----- calculate source influence matrix for airfoil surface if it doesn't exist
       DO 10 J=1, N
C
C------- multiply each dPsi/Sig vector by inverse of factored dPsi/dGam matrix
         CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,BIJ(1,J))
C
C------- store resulting dGam/dSig = dQtan/dSig vector
         DO 105 I=1, N
           DIJ(I,J) = BIJ(I,J)
  105    CONTINUE
C
   10  CONTINUE
       LADIJ = .TRUE.
C
      ENDIF
C
C---- set up coefficient matrix of dPsi/dm on airfoil surface
      DO 20 I=1, N
        CALL PSWLIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N)
        DO 202 J=N+1, N+NW
          BIJ(I,J) = -DZDM(J)
  202   CONTINUE
   20 CONTINUE
C
C---- set up Kutta condition (no direct source influence)
      DO 32 J=N+1, N+NW
        BIJ(N+1,J) = 0.
   32 CONTINUE
C
C---- sharp TE gamma extrapolation also has no source influence
      IF(SHARP) THEN
       DO 34 J=N+1, N+NW
         BIJ(N,J) = 0.
   34  CONTINUE
      ENDIF
C
C---- multiply by inverse of factored dPsi/dGam matrix
      DO 40 J=N+1, N+NW
        CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,BIJ(1,J))
   40 CONTINUE
C
C---- set the source influence matrix for the wake sources
      DO 50 I=1, N
        DO 510 J=N+1, N+NW
          DIJ(I,J) = BIJ(I,J)
  510   CONTINUE
   50 CONTINUE
C
C**** Now we need to calculate the influence of sources on the wake velocities
C
C---- calculcate dQtan/dGam and dQtan/dSig at the wake points
      DO 70 I=N+1, N+NW
C
        IW = I-N
C
C------ airfoil contribution at wake panel node
        CALL PSILIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N,.FALSE.,.TRUE.)
C
        DO 710 J=1, N
          CIJ(IW,J) = DQDG(J)
  710   CONTINUE
C  
        DO 720 J=1, N
          DIJ(I,J) = DQDM(J)
  720   CONTINUE
C
C------ wake contribution
        CALL PSWLIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N)
C
        DO 730 J=N+1, N+NW
          DIJ(I,J) = DQDM(J)
  730   CONTINUE
C
   70 CONTINUE
C
C---- add on effect of all sources on airfoil vorticity which effects wake Qtan
      DO 80 I=N+1, N+NW
        IW = I-N
C
C------ airfoil surface source contribution first
        DO 810 J=1, N
          SUM = 0.
          DO 8100 K=1, N
            SUM = SUM + CIJ(IW,K)*DIJ(K,J)
 8100     CONTINUE
          DIJ(I,J) = DIJ(I,J) + SUM
  810   CONTINUE
C
C------ wake source contribution next
        DO 820 J=N+1, N+NW
          SUM = 0.
          DO 8200 K=1, N
            SUM = SUM + CIJ(IW,K)*BIJ(K,J)
 8200     CONTINUE
          DIJ(I,J) = DIJ(I,J) + SUM
  820   CONTINUE
C
   80 CONTINUE
C
C---- make sure first wake point has same velocity as trailing edge
      DO 90 J=1, N+NW
        DIJ(N+1,J) = DIJ(N,J)
   90 CONTINUE
C
      LWDIJ = .TRUE.
C
      RETURN
      END


      SUBROUTINE XYWAKE
C-----------------------------------------------------
C     Sets wake coordinate array for current surface 
C     vorticity and/or mass source distributions.
C-----------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      WRITE(*,*) 'Calculating wake trajectory ...'
C
C---- number of wake points
      NW = N/12 + 10*INT(WAKLEN)
      IF(NW.GT.IWX) THEN
       WRITE(*,*)
     &  'Array size (IWX) too small.  Last wake point index reduced.'
       NW = IWX
      ENDIF
C
      DS1 = 0.5*(S(2) - S(1) + S(N) - S(N-1))
      CALL SETEXP(SNEW(N+1),DS1,WAKLEN*CHORD,NW)
C

c      write(*,*) waklen, chord, waklen*chord
c      write(*,*) ds1
c      do i = n+1, n+nw
c        write(*,*) i-n, snew(i)
c      enddo

      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
C
C---- set first wake point a tiny distance behind TE
      I = N+1
      SX = 0.5*(YP(N) - YP(1))
      SY = 0.5*(XP(1) - XP(N))
      SMOD = SQRT(SX**2 + SY**2)
      NX(I) = SX / SMOD
      NY(I) = SY / SMOD
      X(I) = XTE - 0.0001*NY(I)
      Y(I) = YTE + 0.0001*NX(I)
      S(I) = S(N)
C
C---- calculate streamfunction gradient components at first point
      CALL PSILIN(I,X(I),Y(I),1.0,0.0,PSI,PSI_X,.FALSE.,.FALSE.)
      CALL PSILIN(I,X(I),Y(I),0.0,1.0,PSI,PSI_Y,.FALSE.,.FALSE.)
C
C---- set unit vector normal to wake at first point
      NX(I+1) = -PSI_X / SQRT(PSI_X**2 + PSI_Y**2)
      NY(I+1) = -PSI_Y / SQRT(PSI_X**2 + PSI_Y**2)
C
C---- set angle of wake panel normal
      APANEL(I) = ATAN2( PSI_Y , PSI_X )
C
C---- set rest of wake points
      DO 10 I=N+2, N+NW
        DS = SNEW(I) - SNEW(I-1)
C
C------ set new point DS downstream of last point
        X(I) = X(I-1) - DS*NY(I)
        Y(I) = Y(I-1) + DS*NX(I)
        S(I) = S(I-1) + DS
C
        IF(I.EQ.N+NW) GO TO 10
C
C------- calculate normal vector for next point
         CALL PSILIN(I,X(I),Y(I),1.0,0.0,PSI,PSI_X,.FALSE.,.FALSE.)
         CALL PSILIN(I,X(I),Y(I),0.0,1.0,PSI,PSI_Y,.FALSE.,.FALSE.)
C
         NX(I+1) = -PSI_X / SQRT(PSI_X**2 + PSI_Y**2)
         NY(I+1) = -PSI_Y / SQRT(PSI_X**2 + PSI_Y**2)
C
C------- set angle of wake panel normal
         APANEL(I) = ATAN2( PSI_Y , PSI_X )
C
   10 CONTINUE
C
C---- set wake presence flag and corresponding alpha
      LWAKE = .TRUE.
      AWAKE =  ALFA
C
C---- old source influence matrix is invalid for the new wake geometry
      LWDIJ = .FALSE.
C
      RETURN
      END



      SUBROUTINE STFIND
C-----------------------------------------
C     Locates stagnation point arc length 
C     location SST and panel index IST.
C-----------------------------------------
      INCLUDE 'XFOIL.INC'
C
      DO 10 I=1, N-1
        IF(GAM(I).GE.0.0 .AND. GAM(I+1).LT.0.0) GO TO 11
   10 CONTINUE
C
      WRITE(*,*) 'STFIND: Stagnation point not found. Continuing ...'
      I = N/2
C
   11 CONTINUE
C
      IST = I
      DGAM = GAM(I+1) - GAM(I)
      DS = S(I+1) - S(I)
C
C---- evaluate so as to minimize roundoff for very small GAM(I) or GAM(I+1)
      IF(GAM(I) .LT. -GAM(I+1)) THEN
       SST = S(I)   - DS*(GAM(I)  /DGAM)
      ELSE
       SST = S(I+1) - DS*(GAM(I+1)/DGAM)
      ENDIF
C
C---- tweak stagnation point if it falls right on a node (very unlikely)
      IF(SST .LE. S(I)  ) SST = S(I)   + 1.0E-7
      IF(SST .GE. S(I+1)) SST = S(I+1) - 1.0E-7
C
      SST_GO = (SST  - S(I+1))/DGAM
      SST_GP = (S(I) - SST   )/DGAM
C
      RETURN
      END


      SUBROUTINE IBLPAN
C-------------------------------------------------------------
C     Sets  BL location -> panel location  pointer array IPAN
C-------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- top surface first
      IS = 1
C
      IBL = 1
      DO 10 I=IST, 1, -1
        IBL = IBL+1
        IPAN(IBL,IS) = I
        VTI(IBL,IS) = 1.0
   10 CONTINUE
C
      IBLTE(IS) = IBL
      NBL(IS) = IBL
C
C---- bottom surface next
      IS = 2
C
      IBL = 1
      DO 20 I=IST+1, N
        IBL = IBL+1
        IPAN(IBL,IS) = I
        VTI(IBL,IS) = -1.0
   20 CONTINUE
C
C---- wake
      IBLTE(IS) = IBL
C
      DO 25 IW=1, NW
        I = N+IW
        IBL = IBLTE(IS)+IW
        IPAN(IBL,IS) = I
         VTI(IBL,IS) = -1.0
   25 CONTINUE
C
      NBL(IS) = IBLTE(IS) + NW
C
C---- upper wake pointers (for plotting only)
      DO 35 IW=1, NW
        IPAN(IBLTE(1)+IW,1) = IPAN(IBLTE(2)+IW,2)
         VTI(IBLTE(1)+IW,1) = 1.0
   35 CONTINUE
C
C
      IBLMAX = MAX(IBLTE(1),IBLTE(2)) + NW
      IF(IBLMAX.GT.IVX) THEN
        WRITE(*,*) ' ***  BL array overflow.'
        WRITE(*,*) ' ***  Increase IVX to at least', IBLMAX
        STOP
      ENDIF
C
      LIPAN = .TRUE.
      RETURN
      END


      SUBROUTINE XICALC
C-------------------------------------------------------------
C     Sets BL arc length array on each airfoil side and wake
C-------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
      DATA XFEPS / 1.0E-7 /
C
C---- minimum xi node arc length near stagnation point
      XEPS = XFEPS*(S(N)-S(1))
C
      IS = 1
C
      XSSI(1,IS) = 0.
C
      DO 10 IBL=2, IBLTE(IS)
        I = IPAN(IBL,IS)
        XSSI(IBL,IS) = MAX( SST - S(I) , XEPS )
   10 CONTINUE
C
C
      IS = 2
C
      XSSI(1,IS) = 0.
C
      DO 20 IBL=2, IBLTE(IS)
        I = IPAN(IBL,IS)
        XSSI(IBL,IS) = MAX( S(I) - SST , XEPS )
   20 CONTINUE
C
C
      IS1 = 1
      IS2 = 2
C
      IBL1 = IBLTE(IS1) + 1
      XSSI(IBL1,IS1) = XSSI(IBL1-1,IS1)
C
      IBL2 = IBLTE(IS2) + 1
      XSSI(IBL2,IS2) = XSSI(IBL2-1,IS2)
C
      DO 25 IBL=IBLTE(IS)+2, NBL(IS)
        I = IPAN(IBL,IS)
        DXSSI = SQRT((X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2)
C
        IBL1 = IBLTE(IS1) + IBL - IBLTE(IS)
        IBL2 = IBLTE(IS2) + IBL - IBLTE(IS)
        XSSI(IBL1,IS1) = XSSI(IBL1-1,IS1) + DXSSI
        XSSI(IBL2,IS2) = XSSI(IBL2-1,IS2) + DXSSI
   25 CONTINUE
C
C---- trailing edge flap length to TE gap ratio
      TELRAT = 2.50
C
C---- set up parameters for TE flap cubics
C
ccc   DWDXTE = YP(1)/XP(1) + YP(N)/XP(N)    !!! BUG  2/2/95
C
      CROSP = (XP(1)*YP(N) - YP(1)*XP(N))
     &      / SQRT(  (XP(1)**2 + YP(1)**2)
     &              *(XP(N)**2 + YP(N)**2) )
      DWDXTE = CROSP / SQRT(1.0 - CROSP**2)
C
C---- limit cubic to avoid absurd TE gap widths
      DWDXTE = MAX(DWDXTE,-3.0/TELRAT)
      DWDXTE = MIN(DWDXTE, 3.0/TELRAT)
C
      AA =  3.0 + TELRAT*DWDXTE
      BB = -2.0 - TELRAT*DWDXTE
C
      IF(SHARP) THEN
       DO 30 IW=1, NW
         WGAP(IW) = 0.
   30  CONTINUE
      ELSE
C----- set TE flap (wake gap) array
       IS = 2
       DO 35 IW=1, NW
         IBL = IBLTE(IS) + IW
         ZN = 1.0 - (XSSI(IBL,IS)-XSSI(IBLTE(IS),IS)) / (TELRAT*ANTE)
         WGAP(IW) = 0.
         IF(ZN.GE.0.0) WGAP(IW) = ANTE * (AA + BB*ZN)*ZN**2
   35  CONTINUE
      ENDIF
C
      RETURN
      END


      SUBROUTINE UICALC
C--------------------------------------------------------------
C     Sets inviscid Ue from panel inviscid tangential velocity
C--------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      DO 10 IS=1, 2
        UINV  (1,IS) = 0.
        UINV_A(1,IS) = 0.
        DO 110 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
          UINV  (IBL,IS) = VTI(IBL,IS)*QINV  (I)
          UINV_A(IBL,IS) = VTI(IBL,IS)*QINV_A(I)
  110   CONTINUE
   10 CONTINUE
C
      RETURN
      END


      SUBROUTINE UECALC
C--------------------------------------------------------------
C     Sets viscous Ue from panel viscous tangential velocity
C--------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      DO 10 IS=1, 2
        UEDG(1,IS) = 0.
        DO 110 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
          UEDG(IBL,IS) = VTI(IBL,IS)*QVIS(I)
  110   CONTINUE
   10 CONTINUE
C
      RETURN
      END


      SUBROUTINE QVFUE
C--------------------------------------------------------------
C     Sets panel viscous tangential velocity from viscous Ue
C--------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
          QVIS(I) = VTI(IBL,IS)*UEDG(IBL,IS)
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END


      SUBROUTINE QISET
C-------------------------------------------------------
C     Sets inviscid panel tangential velocity for
C     current alpha.
C-------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
      DO 5 I=1, N+NW
        QINV  (I) =  COSA*QINVU(I,1) + SINA*QINVU(I,2)
        QINV_A(I) = -SINA*QINVU(I,1) + COSA*QINVU(I,2)
    5 CONTINUE
C
      RETURN
      END


      SUBROUTINE GAMQV
      INCLUDE 'XFOIL.INC'
C
      DO 10 I=1, N
        GAM(I)   = QVIS(I)
        GAM_A(I) = QINV_A(I)
   10 CONTINUE
C
      RETURN
      END


      SUBROUTINE STMOVE
C---------------------------------------------------
C     Moves stagnation point location to new panel.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- locate new stagnation point arc length SST from GAM distribution
      ISTOLD = IST
      CALL STFIND
C
      IF(ISTOLD.EQ.IST) THEN
C
C----- recalculate new arc length array
       CALL XICALC
C
      ELSE
C
CCC       WRITE(*,*) 'STMOVE: Resetting stagnation point'
C
C----- set new BL position -> panel position  pointers
       CALL IBLPAN
C
C----- set new inviscid BL edge velocity UINV from QINV
       CALL UICALC
C
C----- recalculate new arc length array
       CALL XICALC
C
C----- set  BL position -> system line  pointers
       CALL IBLSYS
C
       IF(IST.GT.ISTOLD) THEN
C------ increase in number of points on top side (IS=1)
        IDIF = IST-ISTOLD
C
        ITRAN(1) = ITRAN(1) + IDIF
        ITRAN(2) = ITRAN(2) - IDIF
C
C------ move top side BL variables downstream
        DO 110 IBL=NBL(1), IDIF+2, -1
          CTAU(IBL,1) = CTAU(IBL-IDIF,1)
          THET(IBL,1) = THET(IBL-IDIF,1)
          DSTR(IBL,1) = DSTR(IBL-IDIF,1)
          UEDG(IBL,1) = UEDG(IBL-IDIF,1)
  110   CONTINUE            
C
C------ set BL variables between old and new stagnation point
        DUDX = UEDG(IDIF+2,1)/XSSI(IDIF+2,1)
        DO 115 IBL=IDIF+1, 2, -1
          CTAU(IBL,1) = CTAU(IDIF+2,1)
          THET(IBL,1) = THET(IDIF+2,1)
          DSTR(IBL,1) = DSTR(IDIF+2,1)
          UEDG(IBL,1) = DUDX * XSSI(IBL,1)
  115   CONTINUE
C
C------ move bottom side BL variables upstream
        DO 120 IBL=2, NBL(2)
          CTAU(IBL,2) = CTAU(IBL+IDIF,2)
          THET(IBL,2) = THET(IBL+IDIF,2)
          DSTR(IBL,2) = DSTR(IBL+IDIF,2)
          UEDG(IBL,2) = UEDG(IBL+IDIF,2)
  120   CONTINUE            
C
       ELSE
C------ increase in number of points on bottom side (IS=2)
        IDIF = ISTOLD-IST
C
        ITRAN(1) = ITRAN(1) - IDIF
        ITRAN(2) = ITRAN(2) + IDIF
C
C------ move bottom side BL variables downstream
        DO 210 IBL=NBL(2), IDIF+2, -1
          CTAU(IBL,2) = CTAU(IBL-IDIF,2)
          THET(IBL,2) = THET(IBL-IDIF,2)
          DSTR(IBL,2) = DSTR(IBL-IDIF,2)
          UEDG(IBL,2) = UEDG(IBL-IDIF,2)
  210   CONTINUE            
C
C------ set BL variables between old and new stagnation point
        DUDX = UEDG(IDIF+2,2)/XSSI(IDIF+2,2)


c        write(*,*) 'idif Ue xi dudx', 
c     &    idif, UEDG(idif+2,2), xssi(idif+2,2), dudx

        DO 215 IBL=IDIF+1, 2, -1
          CTAU(IBL,2) = CTAU(IDIF+2,2)
          THET(IBL,2) = THET(IDIF+2,2)
          DSTR(IBL,2) = DSTR(IDIF+2,2)
          UEDG(IBL,2) = DUDX * XSSI(IBL,2)
  215   CONTINUE

c        write(*,*) 'Uenew xinew', idif+1, uedg(idif+1,2), xssi(idif+1,2)

C
C------ move top side BL variables upstream
        DO 220 IBL=2, NBL(1)
          CTAU(IBL,1) = CTAU(IBL+IDIF,1)
          THET(IBL,1) = THET(IBL+IDIF,1)
          DSTR(IBL,1) = DSTR(IBL+IDIF,1)
          UEDG(IBL,1) = UEDG(IBL+IDIF,1)
  220   CONTINUE            
       ENDIF
C
C----- tweak Ue so it's not zero, in case stag. point is right on node
       UEPS = 1.0E-7
       DO IS = 1, 2
         DO IBL = 2, NBL(IS)
           I = IPAN(IBL,IS)
           IF(UEDG(IBL,IS).LE.UEPS) THEN
            UEDG(IBL,IS) = UEPS
            QVIS(I) = VTI(IBL,IS)*UEPS
            GAM(I)  = VTI(IBL,IS)*UEPS
           ENDIF
         ENDDO
       ENDDO
C
      ENDIF
C
C---- set new mass array since Ue has been tweaked
      DO 50 IS=1, 2
        DO 510 IBL=2, NBL(IS)
          MASS(IBL,IS) = DSTR(IBL,IS)*UEDG(IBL,IS)
  510   CONTINUE
   50 CONTINUE
C
      RETURN
      END


      SUBROUTINE UESET
C---------------------------------------------------------
C     Sets Ue from inviscid Ue plus all source influence
C---------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
C
          DUI = 0.
          DO 100 JS=1, 2
            DO 1000 JBL=2, NBL(JS)
              J  = IPAN(JBL,JS)
              UE_M = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
              DUI = DUI + UE_M*MASS(JBL,JS)
 1000       CONTINUE
  100     CONTINUE
C
          UEDG(IBL,IS) = UINV(IBL,IS) + DUI
C
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END


      SUBROUTINE DSSET
      INCLUDE 'XFOIL.INC'
C
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          DSTR(IBL,IS) = MASS(IBL,IS) / UEDG(IBL,IS)
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END
      
      
      
C***********************************************************************
C    Module:  xsolve.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************


      SUBROUTINE GAUSS(NSIZ,NN,Z,R,NRHS)
C     *******************************************************
C     *                                                     *
C     *   Solves general NxN system in NN unknowns          *
C     *    with arbitrary number (NRHS) of righthand sides. *
C     *   Assumes system is invertible...                   *
C     *    ...if it isn't, a divide by zero will result.    *
C     *                                                     *
C     *   Z is the coefficient matrix...                    *
C     *     ...destroyed during solution process.           *
C     *   R is the righthand side(s)...                     *
C     *     ...replaced by the solution vector(s).          *
C     *                                                     *
C     *                              Mark Drela  1984       *
C     *******************************************************
C
      DIMENSION Z(NSIZ,NSIZ), R(NSIZ,NRHS)
C
      DO 1 NP=1, NN-1
        NP1 = NP+1
C
C------ find max pivot index NX
        NX = NP
        DO 11 N=NP1, NN
          IF(ABS(Z(N,NP))-ABS(Z(NX,NP))) 11,11,111
  111      NX = N
   11   CONTINUE
C
        PIVOT = 1.0/Z(NX,NP)
C
C------ switch pivots
        Z(NX,NP) = Z(NP,NP)
C
C------ switch rows & normalize pivot row
        DO 12 L=NP1, NN
          TEMP = Z(NX,L)*PIVOT
          Z(NX,L) = Z(NP,L)
          Z(NP,L) = TEMP
   12   CONTINUE
C
        DO 13 L=1, NRHS
          TEMP = R(NX,L)*PIVOT
          R(NX,L) = R(NP,L)
          R(NP,L) = TEMP
   13   CONTINUE
C
C------ forward eliminate everything
        DO 15 K=NP1, NN
          ZTMP = Z(K,NP)
C
C          IF(ZTMP.EQ.0.0) GO TO 15
C
          DO 151 L=NP1, NN
            Z(K,L) = Z(K,L) - ZTMP*Z(NP,L)
  151     CONTINUE
          DO 152 L=1, NRHS
            R(K,L) = R(K,L) - ZTMP*R(NP,L)
  152     CONTINUE
   15   CONTINUE
C
    1 CONTINUE
C
C---- solve for last row
      DO 2 L=1, NRHS
        R(NN,L) = R(NN,L)/Z(NN,NN)
    2 CONTINUE
C
C---- back substitute everything
      DO 3 NP=NN-1, 1, -1
        NP1 = NP+1
        DO 31 L=1, NRHS
          DO 310 K=NP1, NN
            R(NP,L) = R(NP,L) - Z(NP,K)*R(K,L)
  310     CONTINUE
   31   CONTINUE
    3 CONTINUE
C
      RETURN
      END ! GAUSS


      SUBROUTINE CGAUSS(NSIZ,NN,Z,R,NRHS)
C********************************************
C     Solves general complex linear systems.
C********************************************
      COMPLEX Z(NSIZ,NSIZ), R(NSIZ,NRHS)
      COMPLEX PIVOT, TEMP, ZTMP
C
      DO 1 NP=1, NN-1
        NP1 = NP+1
C
C------ find max pivot index NX
        NX = NP
        DO 11 N=NP1, NN
          IF(ABS(Z(N,NP))-ABS(Z(NX,NP))) 11,11,111
  111      NX = N
   11   CONTINUE
C
        PIVOT = (1.0,0.0)/Z(NX,NP)
C
C------ switch pivots
        Z(NX,NP) = Z(NP,NP)
C
C------ switch rows & normalize pivot row
        DO 12 L=NP1, NN
          TEMP = Z(NX,L)*PIVOT
          Z(NX,L) = Z(NP,L)
          Z(NP,L) = TEMP
   12   CONTINUE
C
        DO 13 L=1, NRHS
          TEMP = R(NX,L)*PIVOT
          R(NX,L) = R(NP,L)
          R(NP,L) = TEMP
   13   CONTINUE
C
C------ forward eliminate everything
        DO 15 K=NP1, NN
          ZTMP = Z(K,NP)
C
C          IF(ZTMP.EQ.0.0) GO TO 15
C
          DO 151 L=NP1, NN
            Z(K,L) = Z(K,L) - ZTMP*Z(NP,L)
  151     CONTINUE
          DO 152 L=1, NRHS
            R(K,L) = R(K,L) - ZTMP*R(NP,L)
  152     CONTINUE
   15   CONTINUE
C
    1 CONTINUE
C
C---- solve for last row
      DO 2 L=1, NRHS
        R(NN,L) = R(NN,L)/Z(NN,NN)
    2 CONTINUE
C
C---- back substitute everything
      DO 3 NP=NN-1, 1, -1
        NP1 = NP+1
        DO 31 L=1, NRHS
          DO 310 K=NP1, NN
            R(NP,L) = R(NP,L) - Z(NP,K)*R(K,L)
  310     CONTINUE
   31   CONTINUE
    3 CONTINUE
C
      RETURN
      END ! CGAUSS



      SUBROUTINE LUDCMP(NSIZ,N,A,INDX)
C     *******************************************************
C     *                                                     *
C     *   Factors a full NxN matrix into an LU form.        *
C     *   Subr. BAKSUB can back-substitute it with some RHS.*
C     *   Assumes matrix is non-singular...                 *
C     *    ...if it isn't, a divide by zero will result.    *
C     *                                                     *
C     *   A is the matrix...                                *
C     *     ...replaced with its LU factors.                *
C     *                                                     *
C     *                              Mark Drela  1988       *
C     *******************************************************
C
      DIMENSION A(NSIZ,NSIZ), INDX(NSIZ)
C
      PARAMETER (NVX=500)
      DIMENSION VV(NVX)
C
      IF(N.GT.NVX) STOP 'LUDCMP: Array overflow. Increase NVX.'
C
      DO 12 I=1, N
        AAMAX = 0.
        DO 11 J=1, N
          AAMAX = MAX( ABS(A(I,J)) , AAMAX )
   11   CONTINUE
        VV(I) = 1.0/AAMAX
   12 CONTINUE
C
      DO 19 J=1, N
        DO 14 I=1, J-1
          SUM = A(I,J)
          DO 13 K=1, I-1
            SUM = SUM - A(I,K)*A(K,J)
   13     CONTINUE
          A(I,J) = SUM
   14   CONTINUE
C
        AAMAX = 0.
        DO 16 I=J, N
          SUM = A(I,J)
          DO 15 K=1, J-1
            SUM = SUM - A(I,K)*A(K,J)
   15     CONTINUE
          A(I,J) = SUM
C
          DUM = VV(I)*ABS(SUM)
          IF(DUM.GE.AAMAX) THEN
           IMAX = I
           AAMAX = DUM
          ENDIF
   16   CONTINUE
C
        IF(J.NE.IMAX) THEN
         DO 17 K=1, N
           DUM = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = DUM
   17    CONTINUE
         VV(IMAX) = VV(J)
        ENDIF
C
        INDX(J) = IMAX
        IF(J.NE.N) THEN
         DUM = 1.0/A(J,J)
         DO 18 I=J+1, N
           A(I,J) = A(I,J)*DUM
   18    CONTINUE
        ENDIF
C
   19 CONTINUE
C
      RETURN
      END ! LUDCMP


      SUBROUTINE BAKSUB(NSIZ,N,A,INDX,B)
      DIMENSION A(NSIZ,NSIZ), B(NSIZ), INDX(NSIZ)
C
      II = 0
      DO 12 I=1, N
        LL = INDX(I)
        SUM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0) THEN
         DO 11 J=II, I-1
           SUM = SUM - A(I,J)*B(J)
   11    CONTINUE
        ELSE IF(SUM.NE.0.0) THEN
         II = I
        ENDIF
        B(I) = SUM
   12 CONTINUE
C
      DO 14 I=N, 1, -1
        SUM = B(I)
        IF(I.LT.N) THEN
         DO 13 J=I+1, N
           SUM = SUM - A(I,J)*B(J)
   13    CONTINUE
        ENDIF
        B(I) = SUM/A(I,I)
   14 CONTINUE
C
      RETURN
      END ! BAKSUB



      SUBROUTINE BLSOLV
C-----------------------------------------------------------------
C      Custom solver for coupled viscous-inviscid Newton system:
C
C        A  |  |  .  |  |  .  |    d       R       S
C        B  A  |  .  |  |  .  |    d       R       S
C        |  B  A  .  |  |  .  |    d       R       S
C        .  .  .  .  |  |  .  |    d   =   R - dRe S
C        |  |  |  B  A  |  .  |    d       R       S
C        |  Z  |  |  B  A  .  |    d       R       S
C        .  .  .  .  .  .  .  |    d       R       S
C        |  |  |  |  |  |  B  A    d       R       S
C
C       A, B, Z  3x3  blocks containing linearized BL equation coefficients
C       |        3x1  vectors containing mass defect influence 
C                     coefficients on Ue
C       d        3x1  unknown vectors (Newton deltas for Ctau, Theta, m)
C       R        3x1  residual vectors
C       S        3x1  Re influence vectors
C-----------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      IVTE1 = ISYS(IBLTE(1),1)
C
      VACC1 = VACCEL
      VACC2 = VACCEL * 2.0 / (S(N) - S(1))
      VACC3 = VACCEL * 2.0 / (S(N) - S(1))
C
      DO 1000 IV=1, NSYS
C
        IVP = IV + 1
C
C====== Invert VA(IV) block
C
C------ normalize first row
        PIVOT = 1.0 / VA(1,1,IV)
        VA(1,2,IV) = VA(1,2,IV) * PIVOT
        DO 10 L=IV, NSYS
          VM(1,L,IV) = VM(1,L,IV)*PIVOT
   10   CONTINUE
        VDEL(1,1,IV) = VDEL(1,1,IV)*PIVOT
        VDEL(1,2,IV) = VDEL(1,2,IV)*PIVOT
C
C------ eliminate lower first column in VA block
        DO 15 K=2, 3
          VTMP = VA(K,1,IV)
          VA(K,2,IV) = VA(K,2,IV) - VTMP*VA(1,2,IV)
          DO 150 L=IV, NSYS
            VM(K,L,IV) = VM(K,L,IV) - VTMP*VM(1,L,IV)
  150     CONTINUE
          VDEL(K,1,IV) = VDEL(K,1,IV) - VTMP*VDEL(1,1,IV)
          VDEL(K,2,IV) = VDEL(K,2,IV) - VTMP*VDEL(1,2,IV)
   15   CONTINUE
C
C
C------ normalize second row
        PIVOT = 1.0 / VA(2,2,IV)
        DO 20 L=IV, NSYS
          VM(2,L,IV) = VM(2,L,IV)*PIVOT
   20   CONTINUE
        VDEL(2,1,IV) = VDEL(2,1,IV)*PIVOT
        VDEL(2,2,IV) = VDEL(2,2,IV)*PIVOT
C
C------ eliminate lower second column in VA block
        K = 3
        VTMP = VA(K,2,IV)
        DO 250 L=IV, NSYS
          VM(K,L,IV) = VM(K,L,IV) - VTMP*VM(2,L,IV)
  250   CONTINUE
        VDEL(K,1,IV) = VDEL(K,1,IV) - VTMP*VDEL(2,1,IV)
        VDEL(K,2,IV) = VDEL(K,2,IV) - VTMP*VDEL(2,2,IV)
C
C
C------ normalize third row
        PIVOT = 1.0/VM(3,IV,IV)
        DO 350 L=IVP, NSYS
          VM(3,L,IV) = VM(3,L,IV)*PIVOT
  350   CONTINUE
        VDEL(3,1,IV) = VDEL(3,1,IV)*PIVOT
        VDEL(3,2,IV) = VDEL(3,2,IV)*PIVOT
C
C
C------ eliminate upper third column in VA block
        VTMP1 = VM(1,IV,IV)
        VTMP2 = VM(2,IV,IV)
        DO 450 L=IVP, NSYS
          VM(1,L,IV) = VM(1,L,IV) - VTMP1*VM(3,L,IV)
          VM(2,L,IV) = VM(2,L,IV) - VTMP2*VM(3,L,IV)
  450   CONTINUE
        VDEL(1,1,IV) = VDEL(1,1,IV) - VTMP1*VDEL(3,1,IV)
        VDEL(2,1,IV) = VDEL(2,1,IV) - VTMP2*VDEL(3,1,IV)
        VDEL(1,2,IV) = VDEL(1,2,IV) - VTMP1*VDEL(3,2,IV)
        VDEL(2,2,IV) = VDEL(2,2,IV) - VTMP2*VDEL(3,2,IV)
C
C------ eliminate upper second column in VA block
        VTMP = VA(1,2,IV)
        DO 460 L=IVP, NSYS
          VM(1,L,IV) = VM(1,L,IV) - VTMP*VM(2,L,IV)
  460   CONTINUE
        VDEL(1,1,IV) = VDEL(1,1,IV) - VTMP*VDEL(2,1,IV)
        VDEL(1,2,IV) = VDEL(1,2,IV) - VTMP*VDEL(2,2,IV)
C
C
        IF(IV.EQ.NSYS) GO TO 1000
C
C====== Eliminate VB(IV+1) block, rows  1 -> 3
        DO 50 K=1, 3
          VTMP1 = VB(K, 1,IVP)
          VTMP2 = VB(K, 2,IVP)
          VTMP3 = VM(K,IV,IVP)
          DO 510 L=IVP, NSYS
            VM(K,L,IVP) = VM(K,L,IVP)
     &        - (  VTMP1*VM(1,L,IV)
     &           + VTMP2*VM(2,L,IV)
     &           + VTMP3*VM(3,L,IV) )
  510     CONTINUE
          VDEL(K,1,IVP) = VDEL(K,1,IVP)
     &        - (  VTMP1*VDEL(1,1,IV)
     &           + VTMP2*VDEL(2,1,IV)
     &           + VTMP3*VDEL(3,1,IV) )
          VDEL(K,2,IVP) = VDEL(K,2,IVP)
     &        - (  VTMP1*VDEL(1,2,IV)
     &           + VTMP2*VDEL(2,2,IV)
     &           + VTMP3*VDEL(3,2,IV) )
   50   CONTINUE
C
        IF(IV.EQ.IVTE1) THEN
C------- eliminate VZ block
         IVZ = ISYS(IBLTE(2)+1,2)
C
         DO 55 K=1, 3
           VTMP1 = VZ(K,1)
           VTMP2 = VZ(K,2)
           DO 515 L=IVP, NSYS
             VM(K,L,IVZ) = VM(K,L,IVZ)
     &         - (  VTMP1*VM(1,L,IV)
     &            + VTMP2*VM(2,L,IV) )
  515      CONTINUE
           VDEL(K,1,IVZ) = VDEL(K,1,IVZ)
     &         - (  VTMP1*VDEL(1,1,IV)
     &            + VTMP2*VDEL(2,1,IV) )
           VDEL(K,2,IVZ) = VDEL(K,2,IVZ)
     &         - (  VTMP1*VDEL(1,2,IV)
     &            + VTMP2*VDEL(2,2,IV) )
   55    CONTINUE
        ENDIF
C
        IF(IVP.EQ.NSYS) GO TO 1000
C
C====== Eliminate lower VM column
        DO 60 KV=IV+2, NSYS
          VTMP1 = VM(1,IV,KV)
          VTMP2 = VM(2,IV,KV)
          VTMP3 = VM(3,IV,KV)
C
          IF(ABS(VTMP1).GT.VACC1) THEN
          DO 610 L=IVP, NSYS
            VM(1,L,KV) = VM(1,L,KV) - VTMP1*VM(3,L,IV)
  610     CONTINUE
          VDEL(1,1,KV) = VDEL(1,1,KV) - VTMP1*VDEL(3,1,IV)
          VDEL(1,2,KV) = VDEL(1,2,KV) - VTMP1*VDEL(3,2,IV)
          ENDIF
C
          IF(ABS(VTMP2).GT.VACC2) THEN
          DO 620 L=IVP, NSYS
            VM(2,L,KV) = VM(2,L,KV) - VTMP2*VM(3,L,IV)
  620     CONTINUE
          VDEL(2,1,KV) = VDEL(2,1,KV) - VTMP2*VDEL(3,1,IV)
          VDEL(2,2,KV) = VDEL(2,2,KV) - VTMP2*VDEL(3,2,IV)
          ENDIF
C
          IF(ABS(VTMP3).GT.VACC3) THEN
          DO 630 L=IVP, NSYS
            VM(3,L,KV) = VM(3,L,KV) - VTMP3*VM(3,L,IV)
  630     CONTINUE
          VDEL(3,1,KV) = VDEL(3,1,KV) - VTMP3*VDEL(3,1,IV)
          VDEL(3,2,KV) = VDEL(3,2,KV) - VTMP3*VDEL(3,2,IV)
          ENDIF
C
   60   CONTINUE
C
 1000 CONTINUE
C
C
C
      DO 2000 IV=NSYS, 2, -1
C
C------ eliminate upper VM columns
        VTMP = VDEL(3,1,IV)
        DO 81 KV=IV-1, 1, -1
          VDEL(1,1,KV) = VDEL(1,1,KV) - VM(1,IV,KV)*VTMP
          VDEL(2,1,KV) = VDEL(2,1,KV) - VM(2,IV,KV)*VTMP
          VDEL(3,1,KV) = VDEL(3,1,KV) - VM(3,IV,KV)*VTMP
   81   CONTINUE
C
        VTMP = VDEL(3,2,IV)
        DO 82 KV=IV-1, 1, -1
          VDEL(1,2,KV) = VDEL(1,2,KV) - VM(1,IV,KV)*VTMP
          VDEL(2,2,KV) = VDEL(2,2,KV) - VM(2,IV,KV)*VTMP
          VDEL(3,2,KV) = VDEL(3,2,KV) - VM(3,IV,KV)*VTMP
   82   CONTINUE
C
 2000 CONTINUE
C
      RETURN
      END
      
      
      



      SUBROUTINE SETEXP(S,DS1,SMAX,NN)
C........................................................
C     Sets geometrically stretched array S:
C
C       S(i+1) - S(i)  =  r * [S(i) - S(i-1)]
C
C       S     (output)  array to be set  
C       DS1   (input)   first S increment:  S(2) - S(1)
C       SMAX  (input)   final S value:      S(NN)
C       NN    (input)   number of points
C........................................................
      REAL S(NN)
C
      SIGMA = SMAX/DS1
      NEX = NN-1
      RNEX = FLOAT(NEX)
      RNI = 1.0/RNEX
C
C---- solve quadratic for initial geometric ratio guess
      AAA = RNEX*(RNEX-1.0)*(RNEX-2.0) / 6.0
      BBB = RNEX*(RNEX-1.0) / 2.0
      CCC = RNEX - SIGMA
C
      DISC = BBB**2 - 4.0*AAA*CCC
      DISC = MAX( 0.0 , DISC )
C
      IF(NEX.LE.1) THEN
       STOP 'SETEXP: Cannot fill array.  N too small.'
      ELSE IF(NEX.EQ.2) THEN
       RATIO = -CCC/BBB  +  1.0
      ELSE
       RATIO = (-BBB + SQRT(DISC))/(2.0*AAA)  +  1.0
      ENDIF
C
      IF(RATIO.EQ.1.0) GO TO 11
C
C---- Newton iteration for actual geometric ratio
      DO 1 ITER=1, 100
        SIGMAN = (RATIO**NEX - 1.0) / (RATIO - 1.0)
        RES = SIGMAN**RNI - SIGMA**RNI
        DRESDR = RNI*SIGMAN**RNI
     &         * (RNEX*RATIO**(NEX-1) - SIGMAN) / (RATIO**NEX - 1.0)
C
        DRATIO = -RES/DRESDR
        RATIO = RATIO + DRATIO
C
        IF(ABS(DRATIO) .LT. 1.0E-5) GO TO 11
C
    1 CONTINUE
      WRITE(*,*) 'SETEXP: Convergence failed.  Continuing anyway ...'
C
C---- set up stretched array using converged geometric ratio
   11 S(1) = 0.0
      DS = DS1
      DO 2 N=2, NN
        S(N) = S(N-1) + DS
        DS = DS*RATIO
    2 CONTINUE
C
      RETURN
      END



      FUNCTION ATANC(Y,X,THOLD)
      IMPLICIT REAL (A-H,M,O-Z)
C---------------------------------------------------------------
C     ATAN2 function with branch cut checking.
C
C     Increments position angle of point X,Y from some previous
C     value THOLD due to a change in position, ensuring that the
C     position change does not cross the ATAN2 branch cut
C     (which is in the -x direction).  For example:
C
C       ATANC( -1.0 , -1.0 , 0.75*pi )  returns  1.25*pi , whereas
C       ATAN2( -1.0 , -1.0 )            returns  -.75*pi .
C
C     Typically, ATANC is used to fill an array of angles:
C
C        THETA(1) = ATAN2( Y(1) , X(1) )
C        DO i=2, N
C          THETA(i) = ATANC( Y(i) , X(i) , THETA(i-1) )
C        END DO
C
C     This will prevent the angle array THETA(i) from jumping by 
C     +/- 2 pi when the path X(i),Y(i) crosses the negative x axis.
C
C     Input:
C       X,Y     point position coordinates
C       THOLD   position angle of nearby point
C
C     Output:
C       ATANC   position angle of X,Y
C---------------------------------------------------------------
      DATA  PI /3.1415926535897932384/
      DATA TPI /6.2831853071795864769/
C
C---- set new position angle, ignoring branch cut in ATAN2 function for now
      THNEW = ATAN2( Y , X )
      DTHET = THNEW - THOLD
C
C---- angle change cannot exceed +/- pi, so get rid of any multiples of 2 pi 
      DTCORR = DTHET - TPI*INT( (DTHET + SIGN(PI,DTHET))/TPI )
C
C---- set correct new angle
      ATANC = THOLD + DTCORR
C
      RETURN
      END ! ATANC   
