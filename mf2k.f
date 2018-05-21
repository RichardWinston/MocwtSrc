C     ******************************************************************
C     MAIN CODE FOR U.S. GEOLOGICAL SURVEY MODULAR MODEL -- MODFLOW
C           BY MICHAEL G. MCDONALD AND ARLEN W. HARBAUGH
C     MODFLOW-88 documented in:
C        McDonald, M.G. and Harbaugh, A.W., 1988, A modular three-
C           dimensional finite-difference ground-water flow model:
C           U.S. Geological Survey Techniques of Water Resources
C           Investigations, Book 6, Chapter A1, 586 p.
C     MODFLOW-96 documented in:
C        Harbaugh, A.W. and McDonald, M.G., 1996, User's documentation
C           for the U.S. Geological Survey modular finite-difference
C           ground-water flow model: U.S. Geological Survey Open-File
C           Report 96-485
C     MODFLOW-2000 documented in:
C        Harbaugh, A.W., Banta, E.R., Hill, M.C., and McDonald, M.G.,
C           2000, MODFLOW-2000, the U.S. Geological Survey modular
C           ground-water model--User guide to modularization concepts
C           and the Ground-Water Flow Process: U.S. Geological Survey
C           Open-File Report 00-92
C        Hill, M.C., Banta, E.R., Harbaugh, A.W., and Anderman, E.R.,
C           2000, MODFLOW-2000, the U.S. Geological Survey modular
C           ground-water model--User guide to the Observation,
C           Sensitivity, and Parameter-Estimation Processes and three
C           post-processing programs: U.S. Geological Survey Open-
C           File Report 00-184
C     MOC3D Version 1.0 11/08/96
C     Documented in:
C           Konikow, L.F., Goode, D.J., and Hornberger, G.Z., 1996,
C           A Three-Dimensional Method-of-Characteristics Solute-
C           Transport Model (MOC3D): U.S. Geological Survey Water-
C           Resources Investigations Report 96-4267, 87 p.
C     MOC3D Version 2.0 11/16/98
C     Documented in:
C           Kipp, K.L., Konikow, L.F., and Hornberger, G.Z., 1998,
C           An Implicit Dispersive Transport Algorithm for
C           the MOC3D Solute-Transport Model: U.S. Geological Survey
C           Water-Resources Investigations Report 98-4234, 54 p.
C     MOC3D VERSION 3.0 1999/03/24
C     Documented in:
C          Goode, D.J., 1999, Age, Double Porosity, and Simple Reaction
C          Modifications for MOC3D Ground-Water Transport Model:
C          U.S. Geological Survey Water-Resources Investigations
C          Report 99-4041, 34 p.
C     MOC3D VERSION 3.5 July 2000 
C     Documented in:
C          Heberton, C.I., Russell, T.F., Konikow, L.F., Hornberger, 
C          G.Z., 2000, A Three-Dimensional Finite Volume Eulerian-
C          Lagrangian Localized Adjoint Method (ELLAM) for Solute-
C          Transport Modeling: U.S. Geological Survey Water-Resources 
C          Investigations Report 00-4087, 63 p.
C     Compatibility with LAK and GAGE Packages December 2000
C     Documented in:
C          Merritt, M.L., and Konikow, L.F., 2000, Documentation of
C          a computer program to simulate lake-aquifer interaction
C          using the MODFLOW ground-water flow model and the MOC3D
C          solute-transport model: U.S. Geological Survey Water-Resources 
C          Investigations Report 00-4167, 146 p. 
C
C	GWT (Ground-Water Transport Process) created January 2001
C		GWT VERSION 1.0 Created from merger of MODFLOW-2000 (Version 1.1)
C			and MOC3D (Version 3.5.1)
C		GWT VERSION 1.0.1 2/22/2001 Fixed bug related to computation of
C             saturated thickness in convertible layers
C		GWT VERSION 1.0.2 4/04/2001 FHB and HFB bug fix (swapped unit 
C             numbers)
C		GWT VERSION 1.0.3 5/01/2001 Update to produce unformatted 
C             unstructured "binary" files (see MODFLOW 1.2 update)
C		GWT VERSION 1.0.4 6/19/2001 Added compatibility with DRT package
C		GWT VERSION 1.0.5 7/16/2001 Compatibility with MF2K 1.4
C		GWT VERSION 1.0.6 8/29/2001 Compatibility with MF2K 1.5
C		GWT VERSION 1.0.7 11/9/2001 Compatibility with MF2K 1.6
C		GWT VERSION 1.0.8 12/18/2001 Compatibility with MF2K 1.7
C         GWT VERSION 1.1 2/7/2002 CHFB Package
C     Documented in:
C          Hornberger,G.Z., Konikow, L.F., Harte, P.T., 2002, Simulating
C          horizontal-flow barriers using the MODFLOW Ground-Water 
C          Transport Process: U.S. Geological Survey Open-File Report 
C          02-52, 28 p.
C         GWT VERSION 1.2 2/27/2002 Compatibility with CHD Package
C         GWT VERSION 1.21 3/1/2002 Saturated thickness of cells in unconfined 
C             layers is recomputed for new steady-state stress periods.  This may
C             cause a 'jump' in solute mass at the transition from one stress 
C             period to the next.
C         GWT VERSION 1.3 5/9/2002 Compatibility with MF2K 1.8
C         GWT VERSION 1.4 8/1/2002 Compatibility with MF2K 1.10  
C             Revised GAGE Package for lakes to allow more comprehensive output.
C             Fixed mass balance errors associated with dispersion terms in 
C             z-direction by accounting for different thickness in layer K+1.
C             Improved roundoff error check in particle moving routine.
C             Improved algorithm for replacing particles at subgrid boundaries.  
C             Added automatic update of saturated thickness for simulation of 
C             unconfined aquifer with multiple steady-state stress periods.
C         GWT VERSION 1.5 11/14/2002 -  For implicit method, fixed mass balance 
C             errors associated with dispersion terms in  z-direction by 
C             accounting for different thickness in layer K+1.  Corrected layer-
C             index calculation for CINFL (C' at subgrid boundary) for certain
C             cases.
C         GWT VERSION 1.6 9/4/2003 -  Compatibility with MF2K 1.11.  Added BFLX  
C             package, including rearrangement and splitting of velocity sub-
C             routines.  Cross-product dispersion terms for explicit and implicit 
C             methods fixed for cases with variable cell thickness.  Source term 
C             fixed for implicit method.  CHFB statements merged into normal dis- 
C             persion subroutines.  Losing lake cells treated as strong sources.
C             Half-life printed to output file.  Changed sign of change in mass 
C             stored budget item (now positive is increase in mass stored).             
C     Documented in:
C              Konikow, L.F., and Hornberger, G.Z., 2003, Use of boundary fluxes
C              when simulating solute transport with the MODFLOW Ground-Water 
C              Transport Process: U.S. Geological Survey Open-File Report 
C              03-303, 17 p.     
C         GWT VERSION 1.7 9/23/2003 -  Compatibility with MF2K 1.12.  Fixed bug
C             related to CHD package.  
C         GWT VERSION 1.8 7/6/2004 -  Compatible with MF2K 1.14, including SFR
C             package.  Fixed bug in dispersion routines for both explicit and 
C             implicit methods.  IGENLK now initialized when LAK package inactive.                    
C             Fixed bug for transient flow conditions with ELLAM method.         
C     Documented in:
C              Konikow, L.F., and Hornberger, G.Z., 2003, Use of boundary fluxes
C              Prudic, D.E., Konikow, L.F., and Banta, E.R., 2004, A new 
C              Streamflow-Routing (SFR1) Package to simulate stream-aquifer 
C              interaction with MODFLOW-2000: U.S. Geological Survey Open-File 
C              Report 2004-1042, 95 p.
C         GWT VERSION 1.8.1 11/29/2004 - Removed extraneous adjustment for 
C              z-direction dispersion calculation.  Reset SMALL (a check for loss
C              of precision) to 1E-4.  NL array redimensioned in ELLAM.  Fixed SFR
C              segments with concentrations associated with external sources.       
C         GWT VERSION 1.9 4/14/2006 Compatible with MF2K 1.16, including GMG 
C              solver.  Compatibility of AGE package with ELLAM solver added.
C              SSTR,CBDY, and PTOB packages added.  See version history for more 
C              details
C              Compatibility with MNW package added (see Konikow and Hornberger, 
C              2006), including MNWO package for concentration output at MNWs.  
C     Documented in:
C              Konikow, L.F., and Hornberger, G.Z., 2006, Use of Multi-Node Well
C              (MNW) Package when simulating solute transport with the MODFLOW
C              Ground-Water Transport Process: U.S. Geological Survey Techniques
C              and Methods 6-A15, 34 p.
C         GWT VERSION 1.9.1 11/9/2006 Compatible with MF2K 1.17.01.  Bug in
C              routing solute in MNWs fixed.  ET option 3 now compatible.  Fixed
C              memory deallocation error when using PTOB package.  Changed algorithm
C              for creating particles at subgrid boundaries.  Velocity output for
C              ELLAM with RF>1 fixed.  CBDY package fixed for no lateral subgrid 
C              boundary.
C         GWT VERSION 1.9.2 11/13/2006 Fixed SRF2 bug to allow ISFROPT=1. 
C              boundary.
C         GWT VERSION 1.9.3 3/28/2007 Compatible with MF2K 1.17.02
C         GWT VERSION 1.9.4 4/20/2007 Fixed error in computing stream (SFR) outflow.
C         GWT VERSION 1.9.5 5/2/2007 Lake concentration updated for case of lake 
C              with inflow tributaries only.  Fixed index for CHFB calculation.
C              Changed CHFB allocation.
C         GWT VERSION 1.9.6 5/16/2007 Fixed calculation of concentration entering
C              lakes from streams.  Corrected budget report of mass entering from
C              subgrid boundary.  Fixed memory leak in MNW simulations.
C         GWT VERSION 1.9.7 7/10/2007 Fixed lake array pointer bug.  Initialized 
C              solute-routing arrays.
C         GWT VERSION 1.9.8 10/22/2008 Fixed SFR array initialization.  BFLX with
C              variable grid spacing bug fix.  Ellam bug with Lakes and/or MNWs fixed.
C         GWT VERSION 1.10 8/2/2017 Added Volume-weighted particle-tracking method
C     Documented in:
C              Winston, R.B., Konikow, L.F., and Hornberger, G.Z., 2017, 
C              Volume-weighted particle-tracking method for solute-transport modeling: 
C              Implementation in MODFLOW-GWT: U.S. Geological Survey Techniques and 
C              Methods 6-A##, xxx p.

C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C-------ASSIGN VERSION NUMBER AND DATE
      CHARACTER*40 VERSION
      PARAMETER (VERSION='1.18.01 06/20/2008')
C
C-----DECLARE ARRAY TYPES
      REAL GX, X, RX, XHS
      DOUBLE PRECISION GZ, VAR, Z, RZ
      INTEGER IG, IX, IR
      CHARACTER*10 EQNAM, NIPRNAM
      CHARACTER*12 NAMES, OBSNAM
      CHARACTER*32 MNWSITE
      CHARACTER*32 MNWOLST
      CHARACTER*32 PTOBLST
      CHARACTER*20 WELLID,MNWIID
C
C *** FOR STATIC MEMORY ALLOCATION, THE FOLLOWING PARAMETER AND
C *** DIMENSION STATEMENTS MUST BE UNCOMMENTED.  TO CHANGE THE SIZE OF
C *** AN ARRAY, CHANGE THE VALUE OF THE CORRESPONDING (FORTRAN)
C *** PARAMETER AND RECOMPILE
C      PARAMETER (LENGX=1000000, LENIG=1000000, LENGZ=1000000,
C     &           LENX=2000000, LENIX=1500000, LENZ=1000000,
C     &           LENRX=1000000, LENIR=1000000, LENXHS=1000000,
C     &           NDD=10000, MPRD=100, IPRD=100)
C      DIMENSION GX(LENGX), IG(LENIG), X(LENX), IX(LENIX), RX(LENRX),
C     &          IR(LENIR), GZ(LENGZ), Z(LENZ), XHS(LENXHS),
C     &          EQNAM(MPRD), NIPRNAM(IPRD), NAMES(NDD+MPRD+IPRD),
C     &          OBSNAM(NDD)
C
C *** FOR STATIC MEMORY ALLOCATION, THE FOLLOWING ALLOCATABLE
C *** STATEMENT MUST BE COMMENTED OUT
      ALLOCATABLE GX(:), IG(:), X(:), IX(:), RX(:), IR(:), GZ(:), Z(:),
     &            RZ(:), XHS(:), NIPRNAM(:), EQNAM(:), NAMES(:), 
     &            OBSNAM(:)
C
      ALLOCATABLE MNWSITE(:)
      ALLOCATABLE MNWOLST(:)
      ALLOCATABLE PTOBLST(:)
cgzh varpt
      ALLOCATABLE NPTLAYA(:,:,:),NPTROWA(:,:,:),NPTCOLA(:,:,:)
cgzh varpt
      ALLOCATABLE SRCFAC(:,:,:)
      ALLOCATABLE WELLID(:),MNWIID(:)
C
      PARAMETER (NIUNIT=100)
      PARAMETER (MXPER=1000)
C
      DIMENSION PERLEN(MXPER),NSTP(MXPER),TSMULT(MXPER),ISSFLG(MXPER)
      CHARACTER*16 VBNM(NIUNIT)
      DIMENSION VBVL(4,NIUNIT),IUNIT(NIUNIT),IREWND(NIUNIT)
CGWT
cgzh debug double sbvl
      DOUBLE PRECISION SBVL
      DIMENSION SBVL(6,NIUNIT),JUNIT(NIUNIT)
      DOUBLE PRECISION DECAY
CGWT
      CHARACTER*80 HEADNG(2)
      DOUBLE PRECISION AP, TOL
C
C  UNCOMMENT "INCLUDE mpif.h" DURING TESTING TO DEVELOP MPI CODE IN THIS
C  ROUTINE OR TO ACTIVATE TIMERS AND DEBUG MODE.
C     INCLUDE 'mpif.h'
      INCLUDE 'parallel.inc'
      INCLUDE 'param.inc'
      INCLUDE 'openspec.inc'
CMOCWT
      INCLUDE 'ptwt.inc'
C-------SPECIFY SIZE OF ARRAY TO HOLD SENSITIVITIES FROM ONE
C-------PARAMETER-ESTIMATION ITERATION TO THE NEXT WHEN PES, SEN, AND 
C-------ANY OBS PACKAGE ARE ACTIVE.  IF IUHEAD IS GREATER THAN ZERO, 
C-------LENXHS MAY EQUAL 1.  IF IUHEAD IS LESS THAN OR EQUAL TO ZERO, 
C-------LENXHS MUST BE AT LEAST:
C-------NLAY*NCOL*NROW*(NUMBER OF PARAMETERS TO BE ESTIMATED).
C
      COMMON /BCFCOM/LAYCON(999)
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
      COMMON /LPFCOM/LAYTYP(999),LAYAVG(999),CHANI(999),LAYVKA(999),
     1               LAYWET(999)
CGWT 
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      COMMON/CAIS/MAR(6,6),MAR1(19,19)
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
C
      INTEGER LAYHDT(999)
C
      CHARACTER*4 PIDTMP
      CHARACTER*20 CHEDFM,CDDNFM,CBOUFM
      CHARACTER*200 FNAME, OUTNAM, COMLIN
      CHARACTER*200 MNWNAME
C
      LOGICAL EXISTS, BEFIRST, SHOWPROG, RESETDD, RESETDDNEXT, OBSALL
      INTEGER NPEVT, NPGHB, NPDRN, NPHFB, NPRIV, NPSTR, NPWEL, NPRCH
CMNW2
      INTEGER WEL1FLAG,BYNDFLAG,QSUMFLAG
      INTEGER IUBE(2), IBDT(8)
      INTEGER IOWELL2(3)  ! FOR MNW1 PACKAGE
      DOUBLE PRECISION SMALL  ! FOR MNW1 PACKAGE
      CHARACTER*4 CUNIT(NIUNIT)
      CHARACTER*10 PARNEG(MXPAR)
      DATA CUNIT/'BCF6', 'WEL ', 'DRN ', 'RIV ', 'EVT ', '    ', 'GHB ',  !  7
     &           'RCH ', 'SIP ', 'DE4 ', 'SOR ', 'OC  ', 'PCG ', 'LMG ',  ! 14
     &           'GWT ', 'FHB ', 'RES ', 'STR ', 'IBS ', 'CHD ', 'HFB6',  ! 21
     &           'LAK ', 'LPF ', 'DIS ', 'SEN ', 'PES ', 'OBS ', 'HOB ',  ! 28
     &           'ADV2', 'COB ', 'ZONE', 'MULT', 'DROB', 'RVOB', 'GBOB',  ! 35
     &           'STOB', 'HUF2', 'CHOB', 'ETS ', 'DRT ', 'DTOB', 'GMG ',  ! 42
     &           'HYD ', 'SFR ', 'sfob', 'GAGE', 'LVDA', '    ', 'LMT6',  ! 49
     &           'MNW1', 'DAF ', 'DAFG', 'KDEP', 'SUB ', 'SWT ', '    ',  ! 56
     &           'MNW2', 'MNWI', '    ', 'unc ', '    ', '    ', '    ',  ! 63
     &           37*'    '/
C     ------------------------------------------------------------------
CGWT
      CHARACTER*4 DUNIT(NIUNIT)
      DATA DUNIT/'CRCH', 'CNCA', 'CNCB', 'PRTA', 'PRTB', 'VELA', 'VELB',  !  7
     1           'OBS ', 'AGE ', 'DP  ', 'DK  ', 'CHFB', 'IPDL', 'IPDA',  ! 14
     2           'BFLX', 'CBDY', 'SSTR', 'MNWO', 'PCT ', 'GAM0', 'GAM1',  ! 21
     3           'MEAN', 'PTOB', 'PRTP', 'CCBD', 'MBRP', 'MBIT', 'VBAL',  ! 28
     4           '    ', '    ', '    ', '    ', '    ', '    ', '    ',  ! 35
     5           '    ', '    ', '    ', '    ', '    ', '    ', '    ',  ! 42
     6           '    ', '    ', '    ', '    ', '    ', '    ', '    ',  ! 49
     &           '    ', 50*'    '/                                       ! 56
CGWT                                                                      ! 63
C     ------------------------------------------------------------------
C RBW begin change
      LCVKA = 1
C RBW end change
      JUNIT=0
      MNWOBS=0
      MOCTYPE=0
      CALL PLL1IN
      IF (MYID.EQ.MPROC) WRITE (*,1) VERSION
    1 FORMAT (/,34X,'MODFLOW-2000',/,
     &4X,'U.S. GEOLOGICAL SURVEY MODULAR FINITE-DIFFERENCE',
     &' GROUND-WATER FLOW MODEL',/,29X,'Version ',A/)
      INUNIT = 99
      IBUNIT = 98
      IBOUTS = 97
      IERRU  = 96
      MAXUNIT= INUNIT
C     DEFINE RANGE OF RESERVED FILE UNITS
      MINRSV = 96
      MAXRSV = 99
      IBATCH = 0
CLAK
      NSOL = 1
      DUM=0.0D0
      IDUM=0
cgzh srcfix
      IALFLAG=0
C
      INQUIRE (FILE='modflow.bf',EXIST=EXISTS)
      IF (EXISTS) THEN
        IBATCH = 1
        IF (MYID.EQ.MPROC) THEN
          OPEN (UNIT=IBUNIT,FILE='modflow.bf',STATUS='OLD')
          OPEN (UNIT=IBOUTS,FILE='modbatch.rpt')
          WRITE (IBOUTS,*) ' USGS MODFLOW MODEL BATCH-MODE REPORT'
        ENDIF
      ENDIF
C2------OPEN FILE OF FILE NAMES.
   10 CONTINUE
      IF (MYID.EQ.MPROC) THEN
        IF (IBATCH.GT.0) THEN
          READ (IBUNIT,'(A)',END=11) FNAME
          GO TO 12
   11     CLOSE(IBUNIT)
          CLOSE(IBOUTS)
          FNAME=' '
          GO TO 15
   12     IF (FNAME.EQ.' ') GOTO 10
          WRITE (IBOUTS,'(1X,/1X,A)') FNAME
        ELSE
          FNAME=' '
          COMLIN=' '
C *** Subroutines GETARG and GETCL are extensions to Fortran 90/95 that
C *** allow a program to retrieve command-line arguments.  To enable
C *** Modflow-2000 to read the name of a Name file from the command
C *** line, either GETARG or GETCL must be called, but not both.  As
C *** distributed, the call to GETARG is uncommented.  For compilers
C *** that support GETCL but not GETARG, comment out the call to GETARG
C *** and uncomment the call to GETCL.  The calls to both GETARG and
C *** GETCL may be commented out for compilers that do not support
C *** either extension.
          CALL GETARG(1,COMLIN)
C          CALL GETCL(COMLIN)
          ICOL = 1
          IF(COMLIN.NE.' ') THEN
            FNAME=COMLIN
          ELSE
            WRITE (*,*) ' Enter the name of the NAME FILE: '
            READ (*,'(A)') FNAME
            CALL URWORD(FNAME,ICOL,ISTART,ISTOP,0,N,R,0,0)
            FNAME=FNAME(ISTART:ISTOP)
          ENDIF
          IF (FNAME.EQ.' ') GOTO 15
          INQUIRE (FILE=FNAME,EXIST=EXISTS)
          IF(.NOT.EXISTS) THEN
            NC=INDEX(FNAME,' ')
            FNAME(NC:NC+3)='.nam'
            INQUIRE (FILE=FNAME,EXIST=EXISTS)
            IF(.NOT.EXISTS) THEN
              WRITE (*,480) FNAME(1:NC-1),FNAME(1:NC+3)
              FNAME=' '
            ENDIF
          ENDIF
        ENDIF
        IF (FNAME.EQ.' ') GOTO 15
        INQUIRE (FILE=FNAME,EXIST=EXISTS)
        IF (.NOT.EXISTS) THEN
          IF (IBATCH.GT.0) THEN
            WRITE (IBOUTS,*) ' Specified name file does not exist.'
            WRITE (IBOUTS,*) ' Processing will continue with the next ',
     &                       'name file in modflow.bf.'
          ENDIF
          GOTO 10
        ENDIF
      ENDIF
   15 CONTINUE
  480 FORMAT(1X,'Can''t find name file ',A,' or ',A)
C
C     BROADCAST FNAME AND OPEN FILE FOR WARNINGS AND ERROR MESSAGES
      CALL PLL1FN(FNAME)
      IF (FNAME.EQ.' ') GOTO 120
      CALL PLL1OP(IERRU,IERR)
      OPEN (UNIT=INUNIT,FILE=FNAME,STATUS='OLD',ACTION=ACTION(1))
      IF (MYID.EQ.MPROC) WRITE(*,490)' Using NAME file: ',FNAME
  490 FORMAT(A,A)
C
C  DEFINE (DF) PROCEDURE
      CALL GLO1BAS6DF(INUNIT,IUNIT,CUNIT,IREWND,NIUNIT,IOUTG,IOUT,
     &                VERSION,NCOL,NROW,NLAY,NPER,ITMUNI,ISUMGX,
     &                MXPER,ISUMIG,ISUMGZ,INBAS,LENUNI,ISUMX,ISUMZ,
     &                ISUMIX,LAYHDT,24,IFREFM,INAMLOC,IPRTIM,IBDT,
     &                SHOWPROG,NOTICECOUNT)
      CALL OBS1BAS6DF(IOBS,IOSTAR,IOWTQ,IOWTQDR,IOWTQGB,
     &                IOWTQRV,IOWTQST,IQ1,IUNIT(27),JT,LCCOFF,LCHFB,
     &                LCIPLO,LCIPLP,LCIQOB,LCNDER,LCNQOB,LCOBADV,
     &                LCOBDRN,LCOBGHB,LCOBBAS,LCOBRIV,LCOBSE,LCOBSTR,
     &                LCQCLS,LCROFF,LCSSAD,LCSSCH,LCSSDR,LCSSGB,LCSSGF,
     &                LCSSPI,LCSSRV,LCSSST,LCSSTO,LCWT,LCWTQ,MOBS,NC,ND,
     &                NDMH,NDMHAR,NH,NOBADV,NQ,NQC,NQT,NQT1,NQTDR,
     &                NQTGB,NQTRV,NQTST,NQTCH,NT,NTT2,IOBSUM,LCX,
     &                LCBUF2,NDAR,LCOBDRT,LCSSDT,NQTDT,IOWTQDT,LCSSSF,
     &                NQTSF,LCOBSFR,IOWTQSF,NHT,LCRSQA,LCRSPA,LCBUF1,
     &                LCH,LCHOBS,LCWTQS,LCHANI,LCXND,LCOTIM,OBSALL)
      CALL SEN1BAS6DF(ISENALL,ISEN,IPRINTS,IUNIT(25),LCB1,LCLN,LCSV,NPE,
     &                NPLIST,RCLOSE,IUHEAD,MXSEN,LCSNEW,IOUTG,LCBSCA,
     &                LCISEN)
      CALL PES1BAS6DF(IBEALE,IBEFLG,IFO,IOUB,IPES,
     &                IPR,IPRAR,IPRINT,ITERPF,
     &                ITERPK,ITMXP,IUNIT(26),IYCFLG,JMAX,LASTX,LCDMXA,
     &                LCNIPR,LCNPAR,LCPRM,LCWP,LCWTP,LCWTPS,LCW3,LCW4,
     &                MPR,MPRAR,NPNGAR,SOSC,SOSR,BEFIRST,LCBPRI,LCPARE,
     &                LCAMPA,LCAMCA,LCAAP)
      CALL GWF1HUF2DF(IOHUFHDS,IOHUFFLWS)
Cdep replaced SFR1 with SFR2
      CALL GWF1SFR2DF(NLAKES,NLAKESAR,LKACC7,LCSTAG,LSLAKE,LSTGLD,
     &                LCRNF,ISTRIN,IDSTRT,ISTROT,ISTGNW,LSCOUT,LSCONQ,
     &                LSCNRN,LSCPPT,LSCNTB,LSSLIN,LSCQIN,LSCGW,LSSIN,
     &                LSSOUT,LSCOTO,ISUZN,NUZST,NSTOTRL,NSTRMSQD,
     &                NUMCELL,NUZROW,NUZCOL,NSSLK)
Cdep end of change
      CALL GWF1LAK3DF(NSS,IDSTRT,ISTRIN,ISTROT,LSLAKE,LSAUG,LSPPT,
     &                LSRNF,LSCGWL,LSSLAK,LSSWIN,LSSWOT,LSSPPT,LSCDRW,
     &                LSSRUN,LSGWIN,LSGWOT,LSOVOL,LSKLK,LSDONE,LSLKSM,
     &                LSFLOB,LSRTCO,LSCLKO,LSALKI,LSALKO,NSSAR,LCSEG,
     &                NSSLK,ISLKOTFLW,IDLKOTFLW,IDLKSTAGE)
Cdep  revised GAG5 to include SFR2
      CALL GWF1GAG5DF(NUMGAGE,LSGAGE,NSTRM,ICSTRM,NLAKES,LKACC7,
     &                LCSTAG,LSLAKE,NLAKESAR,NSTRMAR,NSS,NSSAR,LCIVAR,
     &                NUZST,NUMAVE)
Cdep  end of change
      CALL GWF1MNW1DF(LCHANI,LCHK,LCHKCC,LCHUFTHK,LCHY,LCSSHMN,LCTRPY,
     &                NHUFAR)
C
CGWT----DEFINE GWT (SOLUTE TRANSPORT) PROBLEM
      IF(IUNIT(15).GT.0) 
cgzh varpt  compatibility checks within this call
     *   CALL GWT1BAS6DF(NCOL,NROW,NLAY,IOUTG,IOUTS,IUNIT(15),
     *				INMOC,JUNIT,DUNIT,
     *				NSCOL,NSROW,NSLAY,NODESS,NPMAX,NLIMBO,
     *				NEWPTS,NUMOBS,LSOBSW,
     *				ICSTRT,ICONLY,NODISP,DECAY,DIFFUS,NCINFL,
     *				IABOVE,IBELOW,
     *				IDIM,NPTPND,MOCTYPE,
     *				IDKTIM,IDKRF,IDKZO,IDKFO,IDKZS,IDKFS,
     *				AGER8,IDPZO,IDPFO,IDPTIM,IDPPS,NIUNIT,
     *                JUNIT(13),JUNIT(14))
C
CGWT----READ BOUNDARY FLUX FLAGS
      IF(IUNIT(15).GT.0)
     *   CALL GWT1BFLX5DF(JUNIT(15),IUNIT(8),IUNIT(5),IUNIT(39),
     *     IRCHTP,IEVTTP,NCHNDS,IOUTS,LSBFCH)
C
CGWT----READ SSTR FLAG (FIRST PERIOD FOR TRANSPORT = IPERGWT)
      IF(IUNIT(15).GT.0)
     *   CALL GWT1SSTR5DF(JUNIT(17),IOUTS,IPERGWT)
C
CGWT----READ NUMBER OF GWT OBSERVATION LOCATIONS
C  THIS IS NOT A PART OF THE PARAMETER ESTIMATION OBSERVATIONS
      IF(IUNIT(15).GT.0.AND.JUNIT(8).GT.0)
     *   CALL GWT1OBS5DF(NUMOBS,IOUTS,JUNIT(8),IOBSFL)
C
CGWT----READ NUMBER OF PARTICLE OBSERVATION LOCATIONS
C  THIS IS NOT A PART OF THE PARAMETER ESTIMATION OBSERVATIONS
      IF(IUNIT(15).GT.0.AND.JUNIT(23).GT.0)
     *   CALL GWT1PTOB5DF(NUMPTOB,NUMPTOB_MNW,IOUTS,JUNIT(23))
C
CMOCWT
CGWT----ALLOCATE NPT ARRAYS FOR IPDL AND IPDA OPTIONS
cgzh need this here so npmax can be determined before particle allocations
cgzh ***put npt arrays in GX?  may solve this problem***
      IF(IUNIT(15).GT.0) THEN
	  IF(JUNIT(13).GT.0.OR.JUNIT(14).GT.0) THEN
         ALLOCATE(NPTLAYA(NSCOL,NSROW,NSLAY),NPTROWA(NSCOL,NSROW,NSLAY),
     *           NPTCOLA(NSCOL,NSROW,NSLAY))
        ELSE
	   ALLOCATE(NPTLAYA(1,1,1),NPTROWA(1,1,1),
     *           NPTCOLA(1,1,1))
        END IF
      END IF
CGWT----READ INITIAL PARTICLE DENSITY FILE INFO (LIST-BASED) IPDL
cgzh varpt
      IF(IUNIT(15).GT.0.AND.JUNIT(13).GT.0) 
     *  CALL GWT1IPDL1DF(JUNIT(13),NPTLAY,NPTROW,NPTCOL,NPTLIST,
     *                   NPMAX,IOUTS,NPTLAYA,NPTROWA,NPTCOLA,
     *                   NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NLIMBO,IDIM)
CMOCWT
CGWT----READ INITIAL PARTICLE DENSITY FILE INFO (ARRAY-BASED) IPDA
cgzh varpt
      IF(IUNIT(15).GT.0.AND.JUNIT(14).GT.0)
     *   CALL GWT1IPDA1DF(JUNIT(14),NPMAX,IOUTS,NPTLAYA,NPTROWA,NPTCOLA,
     *       NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NODESS,NLIMBO)
C
C  GLOBAL ALLOCATE (AL) PROCEDURE
      CALL GLO1BAS6AL(IUNIT(24),NCNFBD,NBOTM,NCOL,NROW,NLAY,LCBOTM,
     &                LCDELR,LCDELC,ISUMGX,IOUTG,LCHNEW,LCIBOU,LCCR,
     &                LCCC,LCCV,LCRHS,LCHCOF,LCHOLD,LCBUFF,LCSTRT,
     &                ISUMGZ,ISUMIG,ISEN,IOBS,IPES,ISENALL,ITMXP,IPAR,
     &                IUNIT(31),IUNIT(32),NMLTAR,NZONAR,NML,NZN,LCRMLT,
     &                LCIZON,IUNIT(15))
C
C  DYNAMICALLY ALLOCATE GLOBAL ARRAYS GX, GZ, AND IG.  FOR STATIC
C  MEMORY ALLOCATION, THE FOLLOWING THREE ASSIGNMENT AND ONE ALLOCATE
C  STATEMENTS MUST BE COMMENTED OUT
      LENGX = ISUMGX - 1
      LENGZ = ISUMGZ - 1
      LENIG = ISUMIG - 1
      ALLOCATE (GX(LENGX),GZ(LENGZ),IG(LENIG),STAT=ISTAT)
      IF (ISTAT.NE.0) THEN
        WRITE(*,700)ISTAT
  700   FORMAT(1X,'ALLOCATION OF ARRAYS GX, GZ, AND IG FAILED,',
     &  ' RETURNED ERROR MESSAGE NUMBER: ',I6)
        CALL USTOP(' ')
      ENDIF
C
      CALL MEMCHKG(ISUMGX,ISUMIG,ISUMGZ,LENGX,LENIG,LENGZ,IOUTG,IERR,
     &             IERRU)
      IF (IERR.GT.0) CALL PLL1SD(IERR,IERRU,IOUT,IOUTG)
C
C  GLOBAL READ AND PREPARE (RP) PROCEDURE
      CALL GLO1BAS6RP(IUNIT(24),NCOL,NROW,NLAY,GX(LCBOTM),NBOTM,IOUTG,
     1                GX(LCDELR),GX(LCDELC),NPER,PERLEN,NSTP,TSMULT,
     2                ISSFLG,ITRSS,IUNIT(31),IUNIT(32),NMLTAR,NZONAR,
     3                GX(LCRMLT),IG(LCIZON),NML,NZN)
C
C-----NO rewind AL and RP for Ground-Water Flow Process
      IF(IUNIT(23).GT.0)
     1    CALL GWF1LPF1ALG(ISUMX,LCHK,LCVKA,LCSC1,LCSC2,LCHANI,LCVKCB,
     2                     IUNIT(23),NCOL,NROW,NLAY,IOUTG,ILPFCB,LCWETD,
     3                     HDRY,NPLPF,NCNFBD,LCLAYF,IREWND(23),ISUMIX,
     4                     LAYHDT,ITRSS,LCSV,ISEN)
      IF(IUNIT(37).GT.0) THEN
        CALL GWF1HUF2ALG(ISUMX,LCHK,LCVKA,LCSC1,IUNIT(37),ITRSS,NCOL,
     &                   NROW,NLAY,IOUTG,IHUFCB,LCWETD,HDRY,NPER,
     &                   ISSFLG,LCHGUF,IREWND(37),
     &                   NHUF,NPHUF,LCHUFTHK,LCHKCC,ISUMIX,IOHUFHDS,
     &                   IOHUFFLWS,LAYHDT,LCHUFTMP)
        CALL GWF1HUF2LVDA1ALG(ISUMX,ISUMIX,IUNIT(47),IOUTG,NCOL,
     &                        NROW,NLAY,LCVDHD,LCDVDH,LCVDHT,NPLVDA,
     &                        LCA9)
        CALL GWF1HUF2KDEP1ALG(ISUMX,IUNIT(53),IOUTG,NCOL,NROW,
     &                        LCGS,NPKDEP,IFKDEP)
      ENDIF

      IF(IUNIT(9).GT.0)
     1    CALL SIP5ALG(ISUMX,ISUMIX,LCEL,LCFL,LCGL,LCV,LCHDCG,LCLRCH,
     2                 LCW,MXITER,NPARM,NCOL,NROW,NLAY,IUNIT(9),IOUTG,
     3                 IFREFM,IREWND(9))
      IF(IUNIT(10).GT.0)
     1    CALL DE45ALG(ISUMX,ISUMIX,LCAU,LCAL,LCIUPP,LCIEQP,LCD4B,
     2                 LCLRCH,LCHDCG,MXUP,MXLOW,MXEQ,MXBW,IUNIT(10),
     3                 ITMX,ID4DIR,NCOL,NROW,NLAY,IOUTG,ID4DIM,
     4                 IREWND(10))
      IF(IUNIT(11).GT.0)
     1    CALL SOR5ALG(ISUMX,ISUMIX,LCA,LCRES,LCHDCG,LCLRCH,LCIEQP,
     2                 MXITER,NCOL,NLAY,NSLICE,MBW,IUNIT(11),IOUTG,
     3                 IFREFM,IREWND(11))
      IF(IUNIT(13).GT.0)
     1    CALL PCG2ALG(ISUMX,ISUMIX,LCV,LCSS,LCP,LCCD,LCHCHG,LCLHCH,
     2                 LCRCHG,LCLRCH,MXITER,ITER1,NCOL,NROW,NLAY,
     3                IUNIT(13),IOUTG,NPCOND,LCIT1,LCHCSV,IFREFM,
     4                IREWND(13),ISUMZ,LCHPCG)
      IF(IUNIT(14).GT.0)
     1    CALL LMG1ALG(ISUMZ,ISUMIX,LCA,LCIA,LCJA,LCU1,LCFRHS,
     2                 LCIG,ISIZ1,ISIZ2,ISIZ3,ISIZ4,ICG,NCOL,NROW,NLAY,
     3                 IUNIT(14),IOUTG,1)
c      IF(IUNIT(42).GT.0)
c     1    CALL GMG1ALG(NCOL,NROW,NLAY,MXITER,IITER,RCLOSE,HCLOSE,DAMP,
c     2                 IADAMP,IOUTGMG,IUNIT(42),IOUTG)
C
C-----ALLOCATE SPACE FOR SENSITIVITY CALCULATIONS
      IF (ISEN.GT.0)
     &    CALL SEN1BAS6AL(ISUMX,ISUMIX,NCOL,NROW,NLAY,IOUTG,IUHEAD,
     &                    NPLIST,IUNIT(25),IPAR,LCHCLO,LCRCLO,LCLN,
     &                    IPRINTS,LCISEN,LCBU,LCBL,LCB1,ISENALL,
     &                    IREWND(25),LCSNEW,LCSOLD,ISUMZ,ISEN,ISENSU,
     &                    ISENPU,ISENFM,IPES,MXSEN,LCBSCA,ITMXP,MAXUNIT,
     &                    MINRSV,MAXRSV,NSTP,NPER,NTIMES,LCSEND,LCSNDT)
C-----ALLOCATE SPACE FOR PARAMETER-ESTIMATION PROCESS
      IF (IPES.GT.0)
     &    CALL PES1BAS6AL(ISUMX,ISUMZ,ISUMIX,IOUTG,NPLIST,LCC,LCSCLE,
     &                    LCG,LCDD,LCWP,MPR,LCPRM,LCR,LCU,LCGD,
     &                    LCS,NOPT,IPR,LCWTP,LCWTPS,LCW3,LCW4,LCNIPR,
     &                    LCEIGL,LCEIGV,LCEIGW,LCIPNG,IUNIT(26),
     &                    NPNG,MPRAR,IPRAR,NPNGAR,IREWND(26),
     &                    LCPRNT,LCPARE,ITMXP,LCSSPI,LCSSTO,DMAX,TOL,
     &                    SOSC,IOSTAR,NFIT,SOSR,IPRC,IPRINT,LPRINT,CSA,
     &                    FCONV,LASTX,ISEN,IPES,IPAR,IBEFLG,IYCFLG,
     &                    LCDMXA,LCNPAR,LCBPRI,RMARM,IAP,LCAAP,
     &                    LCAMCA,LCAMPA,RMAR)
C-----READ INPUT RELATED TO ALL OBSERVATIONS AND OPEN
C     PARAMETER-VALUE FILE ON IOUB
      IF (IOBS.GT.0)
     &    CALL OBS1BAS6AL(IOUB,IOUTG,ISCALS,ISEN,IUNIT(27),OUTNAM,
     &                    ISOLDX,ISOLDZ,ISOLDI,ISUMX,ISUMZ,ISUMIX,
     &                    OBSALL)
C-----ALLOCATE SPACE FOR HEAD OBSERVATIONS
      IF (IUNIT(28).GT.0)
     &    CALL OBS1BAS6HAL(IUNIT(28),NH,MOBS,MAXM,ISUMX,ISUMIX,LCNDER,
     &                     LCCOFF,LCROFF,LCIOFF,LCJOFF,LCRINT,LCMLAY,
     &                     LCPR,ND,IOUTG,IOBSUM,LCOBBAS,ITMXP,LCSSGF,
     &                     IOBS,NHT)
C-----ALLOCATE SPACE FOR FLOW OBSERVATIONS
      IF (IUNIT(33).GT.0)
     &    CALL OBS1DRN6AL(IUNIT(33),NQ,NQC,NQT,IOUTG,NQDR,NQTDR,IOBSUM,
     &                    LCOBDRN,ITMXP,LCSSDR,ISUMX,IOBS)
      IF (IUNIT(34).GT.0)
     &    CALL OBS1RIV6AL(IUNIT(34),NQ,NQC,NQT,IOUTG,NQRV,NQTRV,IOBSUM,
     &                    LCOBRIV,ITMXP,LCSSRV,ISUMX,IOBS)
      IF (IUNIT(35).GT.0)
     &    CALL OBS1GHB6AL(IUNIT(35),NQ,NQC,NQT,IOUTG,NQGB,NQTGB,IOBSUM,
     &                    LCOBGHB,ITMXP,LCSSGB,ISUMX,IOBS)
      IF (IUNIT(36).GT.0)
     &    CALL OBS1STR6AL(IUNIT(36),NQ,NQC,NQT,IOUTG,NQST,NQTST,IOBSUM,
     &                    LCOBSTR,ITMXP,LCSSST,ISUMX,IOBS)
      IF (IUNIT(38).GT.0)
     &    CALL OBS1BAS6FAL(IUNIT(38),NQ,NQC,NQT,IOUTG,NQCH,NQTCH,IOBSUM,
     &                     LCOBCHD,ITMXP,LCSSCH,ISUMX,IOBS)
      IF (IUNIT(41).GT.0)
     &    CALL OBS1DRT1AL(IUNIT(41),NQ,NQC,NQT,IOUTG,NQDT,NQTDT,IOBSUM,
     &                    LCOBDRT,ITMXP,LCSSDT,ISUMX,IOBS)
C-----ALLOCATE SPACE FOR ADVECTIVE TRAVEL OBSERVATIONS (ADV PACKAGE)
      IF (IUNIT(29).GT.0)
     &    CALL OBS1ADV2AL(IUNIT(29),NPTH,NTT2,IOUTT2,KTDIM,KTFLG,KTREV,
     &                    ADVSTP,IOUTG,LCICLS,LCPRST,NPRST,LCTT2,LCPOFF,
     &                    LCNPNT,ND,ISUMX,ISUMIX,NROW,NCOL,NLAY,
     &                    IOBSUM,LCOBADV,NOBADV,ITMXP,LCSSAD,IOBS,
     &                    FSNK,NBOTM,IUNIT,NIUNIT,LCDRAI,MXDRN,
     &                    NDRAIN,LCRIVR,MXRIVR,LCBNDS,MXBND,NBOUND,
     &                    LCIRCH,LCRECH,ICSTRM_,LCSTRM_,MXSTRM,NSTREM,
     &                    NDRNVL,NGHBVL,NRIVVL,NRIVER,LCHANI,LCHKCC,
     &                    LCHUFTHK,NHUF,LCGS,LCVDHT,LCDVDH,
     &                    LCWELL,NWELVL,MXWELL,NWELLS,ISEN,IADVHUF,
     &                    IMPATHOUT,IMPATHOFF)
C-----ALLOCATE SPACE FOR ALL OBSERVATIONS AND FOR RESIDUALS RELATED TO
C     OBSERVATIONS AND PRIOR INFORMATION. ALSO INITIALIZE SOME ARRAYS
      IF (IOBS.GT.0)
     &    CALL OBS1BAS6AC(EV,ISUMX,ISUMZ,ISUMIX,LCTOFF,NH,LCH,ND,
     &                    LCHOBS,LCWT,NDMH,NDMHAR,LCWTQ,LCWTQS,LCW1,
     &                    LCW2,LCX,NPLIST,LCXD,IPAR,IOUTG,IDRY,
     &                    JDRY,NQ,NQAR,NQC,NQCAR,NQT,NQTAR,NHAR,MOBS,
     &                    MOBSAR,LCIBT,LCNQOB,LCNQCL,LCIQOB,LCQCLS,
     &                    LCIPLO,LCIPLP,IPR,MPR,IPRAR,LCBUF1,LCSSTO,
     &                    ITMXP,LBUFF,LCOBSE,ISOLDX,ISOLDZ,ISOLDI,MXSEN,
     &                    LCBUF2,NDAR,NHT,LCRSQA,LCRSPA,LCXND,LCOTIM)
C
CGWT----ALLOCATE SPACE FOR GWT (SOLUTE TRANSPORT)
      IF(IUNIT(15).GT.0) 
     *  CALL GWT1BAS6AL(ISUMX,ISUMIX,ISUMZ,LSCB,LSLB,
     *					LSIDB,LSAS,LSIAS,LSJAS,LSA,LSIA,LSJA,LSIW,LSRW,
cea
     *                    LSRW1,
     *					LSRHSE,LSRHSO,LSDO,LSNONU,LSSAV,LSXFOR,LSXBAC,
     *                    LSYFOR,LSYBAC,LSCFOR,LSRFOR,LSTFOR,LSCONL,
     *					LSVOL,LSIACT,LSNZIN,LSDIST,LSDIDS,NTFACE,LENW,
     *					LENIW,LSDCC,LSDCR,LSDCL,LSDRR,LSDRC,LSDRL,
     *					LSDLL,LSDLC,LSDLR,LSPC,LSPR,LSPL,LSPCON,
     *					LSPTID,LSLMBO,LSIGNP,LSRF,LSCONC,LSCINT,
     *					LSCHBC,LSCINF,LSSRCS,LSSRCF,LSSNKF,LSTHCK,
     *					LSNPCL,LSNPLD,LSALNG,LSATH,LSCTCF,LSCHDF,
     *					LSLBDY,LSEVTF,LSATV,LSVC,LSVR,LSVL,LSPOR,
     *					LSSUMC,LSCNCN,LSCOLD,LSCAVG,LSPNWC,LSPNWR,
     *					LSPNWL,LSCI,LSCIN,LSCIR,LSCIRL,LSCIRH,LSIND,
     *					LSMRNO,LSMRNZ,LSRHS,LSVA,LSVAD,LSIBOU,
     *					LSAP,LSPP,LSRA,LSRR,LSSS,LSXX,LSWW,LSWORK,
     *					LSRS,NCOL,NROW,NLAY,
     *					NSCOL,NSROW,NSLAY,NPMAX,NLIMBO,NEWPTS,
     *					IOUTG,IOUTS,ICSTRT,ICONLY,NODISP,DIFFUS,
     *					NCINFL,NFACES,MOCTYPE,
     *					JUNIT(11),IDKTIM,LSDKZO,IDKZO,LSDKFO,IDKFO,
     *					LSDKZS,IDKZS,LSDKFS,IDKFS,JUNIT(10),LSDPCON,
     *					LSDPRAT,LSDPPOR,LSDPZO,IDPZO,LSDPFO,IDPFO,
     *                    LSIGLK,
     *                    LSXDMA,SRCDCY,LSWTFC,
cgzh cbdy
     &                    JUNIT(16),LSCINA,LSCINB,LSCINXY,
cgzh mnw
     &                    IUNIT(57),
     *                    LSSRCM,LSSNKM,LSSOLM,ISSIZH)
CMOCWT
CGWT----ALLOCATE SPACE FOR WEIGHTED PARTICLES OPTION
C
      IF(IUNIT(15).GT.0) THEN
c         WRITE(*,*) 'CALL GWT1PTWT1AL'
         CALL GWT1PTWT1AL(ISUMX,ISUMIX,ISUMZ,LSCELV,LSPTWT,
     *        LSSUMW,LSSRCC,LSSRCV,LSBSRC,LSBSNK,
     *        LSBSOL,LSSGMS,LSSGWT,
     *        LSBFMS,LSBFWT,
     *        LSRESW,LSRESC,LSISRC,LSNPOR,
cgzh srcfix2
     *        LS_SRC,LS_SNK,LS_SOL,LSSOL,
     &        LSTOIN,LSTOOT,LSCOIN,LSCOOT,  
     *        NSROW,NSCOL,NSLAY,NPMAX,IOUTG,IOUTS)
      ENDIF
C
CGWT----ALLOCATE SPACE FOR IPDL OR IPDA PACKAGE
cgzh varpt
        IF(JUNIT(13).GT.0.OR.JUNIT(14).GT.0) 
     *    CALL GWT1IPD1AL(ISUMX,NPMAX,IOUTS,IOUTG,
     *      LSPCOR,LSPROR,LSPLOR)
C
CGWT----ALLOCATE SPACE FOR GWT OBSERVATION LOCATIONS
C  THIS IS NOT A PART OF THE PARAMETER ESTIMATION OBSERVATIONS
      IF(JUNIT(8).GT.0.AND.IUNIT(15).GT.0) THEN
         CALL GWT1OBS5AL(ISUMIX,LSOBSW,NUMOBS,
     *					IOUTG,IOUTS)
      ENDIF
C
CGWT----ALLOCATE SPACE FOR GWT PARTICLE OBSERVATION LOCATIONS
C  THIS IS NOT A PART OF THE PARAMETER ESTIMATION OBSERVATIONS
      IF(IUNIT(15).GT.0) THEN
         CALL GWT1PTOB5AL(ISUMIX,LSPTOB,LSPTMN,LSPTUN,NUMPTOB,
     *                    NUMPTOB_MNW,NSCOL,NSROW,NSLAY,
     *					IOUTG,IOUTS,JUNIT(23))
           IF(NUMPTOB_MNW.GT.0) THEN
            ALLOCATE (PTOBLST(NUMPTOB_MNW),STAT=ISTAT)
            IF (ISTAT.NE.0) THEN
              WRITE(*,1713) ISTAT
 1713         FORMAT(1X,'ALLOCATION OF ARRAY PTOBLST FAILED,',/,
     &        ' RETURNED ERROR MESSAGE NUMBER: ',I6)
              CALL USTOP(' ')
            ENDIF
           ELSE
            ALLOCATE (PTOBLST(1))
           ENDIF
      ENDIF
CGWT BFLX
      IF(IUNIT(15).GT.0) THEN
         CALL GWT1BFLX5AL(ISUMIX,LSBFCH,JUNIT(15),NCHNDS,IOUTS)      
      END IF
CGWT----ALLOCATE FOR PERCENT CHANGE AND GAMMA OUTPUT
      IF(IUNIT(15).GT.0) THEN
         CALL GWT1PCT5AL(ISUMX,LSCINI,LSCTNM,JUNIT,NIUNIT,IOUTS,NODESS)      
      END IF
C------DYNAMICALLY ALLOCATE X, Z, IX, XHS, NIPRNAM, EQNAM, NAMES, AND
C      OBSNAM ARRAYS FOR OBS, SEN, AND PES PROCESSES; SOLVERS; AND
C      PACKAGES THAT DO ALLOCATION ONCE ONLY.  FOR STATIC MEMORY
C      ALLOCATION, THE FOLLOWING LINES, THROUGH THE ALLOCATE STATEMENTS,
C      MUST BE COMMENTED OUT
      LENX = ISUMX - 1
      IF(LENX.LT.1) LENX=1
      LENZ = ISUMZ - 1
      IF(LENZ.LT.1) LENZ=1
      LENIX = ISUMIX - 1
cgzh debug increasing IX for debugging with lahey
c      LENIX = ISUMIX - 1 + 1000
      IF(LENIX.LT.1) LENIX=1
      IF (ISEN.NE.0 .AND. IUHEAD.LE.0 .AND. MXSEN.GT.0) THEN
        LENXHS = NCOL*NROW*NLAY*MXSEN
      ELSE
        LENXHS = 1
      ENDIF
      NDD = NDAR
      MPRD = MPRAR
      IPRD = IPRAR
      ALLOCATE (X(LENX),Z(LENZ),IX(LENIX),XHS(LENXHS),STAT=ISTAT)
      IF (ISTAT.NE.0) THEN
        WRITE(*,701)ISTAT
  701   FORMAT(1X,'ALLOCATION OF ARRAYS X, Z, IX, AND XHS FAILED,',/,
     &  ' RETURNED ERROR MESSAGE NUMBER: ',I6)
        CALL USTOP(' ')
      ENDIF
      ALLOCATE (NIPRNAM(IPRAR),EQNAM(MPRAR),NAMES(ND+IPRAR+MPRAR),
     &          OBSNAM(NDAR),STAT=ISTAT)
      IF (ISTAT.NE.0) THEN
        WRITE(*,702)ISTAT
  702   FORMAT(1X,'ALLOCATION OF ARRAYS NIPRNAM, EQNAM, NAMES, AND',
     &  ' OBSNAM FAILED,',/,
     &  ' RETURNED ERROR MESSAGE NUMBER: ',I6)
        CALL USTOP(' ')
      ENDIF
C
C------IF THE ARRAYS ARE NOT BIG ENOUGH THEN STOP.
      CALL MEMCHK(ISUMX,ISUMIX,ISUMZ,LENX,LENIX,LENZ,IOUTG,ISEN,IUHEAD,
     &            LENXHS,NCOL,NROW,NLAY,MXSEN,IERR,IERRU,NDD,NDAR,MPRD,
     &            MPRAR,IPRD,IPRAR)
      IF (IERR.GT.0) CALL PLL1SD(IERR,IERRU,IOUT,IOUTG)
C
      IF (ISEN.GT.0 .OR. ISENALL.LT.0 .OR. IBEFLG.EQ.2)
     &    CALL SEN1BAS6RP(X(LCBL),X(LCBU),FAC,IX(LCISEN),IOUTG,
     &                    IUNIT(25),IX(LCLN),NPE,NPLIST,DETWTP,ISENALL,
     &                    X(LCBSCA),MXSEN)
      IF (IPES.GT.0 .OR. IBEFLG.EQ.2)
     &    CALL PES1BAS6RP(IUNIT(26),IOUTG,NPE,X(LCWP),IX(LCLN),DMAX,
     &                    Z(LCDD),FCONV,EV,MPR,X(LCPRM),IX(LCISEN),
     &                    NPLIST,X(LCWTP),X(LCWTPS),Z(LCW3),Z(LCW4),IPR,
     &                    IX(LCNIPR),DETWTP,ND,ADMX,AP,DMX,NIPRNAM,
     &                    EQNAM,MPRAR,IPRAR,IX(LCIPNG),NPNG,NPNGAR,
     &                    IX(LCIPLO),NAMES,PARNEG,MXPAR,LBUFF,FSTAT,
     &                    X(LCBPRI),IERR,IYCFLG,IX(LCNPAR),ITMXP,IBEFLG)
C
C-----INITIALIZE ARRAYS USED FOR OBSERVATION PROCESS
      IF (IOBS.GT.0) CALL OBS1BAS6RP(ND,NDAR,NDMH,NDMHAR,NQCAR,
     &                               X(LCQCLS),RSQO,RSQOO,RSQP,X(LCWT),
     &                               X(LCWTQ),X(LCWTQS),X(LCOTIM))
C
C-----READ AND PREPARE INFORMATION FOR OBSERVATIONS
C
C-----READ HEAD OBSERVATION DATA
      IF (IUNIT(28).GT.0)
     &    CALL OBS1BAS6HRP(NCOL,NROW,NLAY,NPER,IUNIT(28),IOUTG,OBSNAM,
     &                    NH,IX(LCNDER),JT,IX(LCJOFF),IX(LCIOFF),
     &                    X(LCHOBS),X(LCWT),GX(LCDELR),GX(LCDELC),
     &                    X(LCRINT),X(LCCOFF),X(LCROFF),IX(LCMLAY),
     &                    X(LCPR),MOBS,IERR,X(LCTOFF),EV,EVH,MAXM,NSTP,
     &                    PERLEN,TSMULT,ISSFLG,ITRSS,NHAR,MOBSAR,
     &                    IX(LCIPLO),NAMES,ND,IPR,MPR,X(LCOTIM))
C-----READ HEAD-DEPENDENT-BOUNDARY FLOW-OBSERVATION DATA
      IF (IUNIT(33).GT.0)
     &    CALL OBS1DRN6RP(NCOL,NROW,NPER,IUNIT(33),IOUTG,OBSNAM,NHT,JT,
     &                    IX(LCIBT),IX(LCNQOB),IX(LCNQCL),
     &                    IX(LCIQOB),X(LCQCLS),IERR,X(LCHOBS),X(LCTOFF),
     &                    X(LCWTQ),IOWTQ,IPRN,NDMH,NSTP,PERLEN,
     &                    TSMULT,ISSFLG,ITRSS,NQAR,NQCAR,
     &                    NQTAR,IQ1,NQT1,NDD,IUNIT(3),NQDR,NQTDR,NT,
     &                    NC,IX(LCIPLO),NAMES,ND,IPR,MPR,IOWTQDR,
     &                    X(LCOTIM))
      IF (IUNIT(34).GT.0)
     &    CALL OBS1RIV6RP(NCOL,NROW,NPER,IUNIT(34),IOUTG,OBSNAM,
     &                    NH,JT,IX(LCIBT),IX(LCNQOB),
     &                    IX(LCNQCL),IX(LCIQOB),X(LCQCLS),IERR,
     &                    X(LCHOBS),X(LCTOFF),X(LCWTQ),IOWTQ,IPRN,
     &                    NDMH,NSTP,PERLEN,TSMULT,
     &                    ISSFLG,ITRSS,NQAR,NQCAR,NQTAR,IQ1,NQT1,
     &                    NDD,IUNIT(4),NQRV,NQTRV,NT,NC,IX(LCIPLO),
     &                    NAMES,ND,IPR,MPR,IOWTQRV,X(LCOTIM))
      IF (IUNIT(35).GT.0)
     &    CALL OBS1GHB6RP(NCOL,NROW,NPER,IUNIT(35),IOUTG,OBSNAM,
     &                    NHT,JT,IX(LCIBT),IX(LCNQOB),
     &                    IX(LCNQCL),IX(LCIQOB),X(LCQCLS),IERR,
     &                    X(LCHOBS),X(LCTOFF),X(LCWTQ),IOWTQ,IPRN,
     &                    NDMH,NSTP,PERLEN,TSMULT,
     &                    ISSFLG,ITRSS,NQAR,NQCAR,NQTAR,IQ1,NQT1,
     &                    NDD,IUNIT(7),NQGB,NQTGB,NT,NC,IX(LCIPLO),
     &                    NAMES,ND,IPR,MPR,IOWTQGB,X(LCOTIM))
      IF (IUNIT(36).GT.0)
     &    CALL OBS1STR6RP(NPER,IUNIT(36),IOUTG,OBSNAM,NHT,JT,
     &                    IX(LCIBT),IX(LCNQOB),IX(LCNQCL),IX(LCIQOB),
     &                    X(LCQCLS),IERR,X(LCHOBS),X(LCTOFF),X(LCWTQ),
     &                    IOWTQ,IPRN,NDMH,NSTP,PERLEN,TSMULT,ISSFLG,
     &                    ITRSS,NQAR,NQCAR,NQTAR,IQ1,NQT1,IUNIT(18),
     &                    NQST,NQTST,NT,NC,IX(LCIPLO),NAMES,ND,IPR,
     &                    MPR,IOWTQST,X(LCOTIM))
      IF (IUNIT(38).GT.0)
     &    CALL OBS1BAS6FRP(NCOL,NROW,NPER,IUNIT(38),IOUTG,OBSNAM,
     &                     NHT,JT,IX(LCIBT),IX(LCNQOB),
     &                     IX(LCNQCL),IX(LCIQOB),X(LCQCLS),IERR,
     &                     X(LCHOBS),X(LCTOFF),X(LCWTQ),IOWTQ,IPRN,
     &                     NDMH,NSTP,PERLEN,TSMULT,ISSFLG,ITRSS,NQAR,
     &                     NQCAR,NQTAR,IQ1,NQT1,NDD,NQCH,NQTCH,NT,NC,
     &                     IX(LCIPLO),NAMES,ND,IPR,MPR,IOWTQCH,NLAY,
     &                     X(LCOTIM))
      IF (IUNIT(41).GT.0)
     &    CALL OBS1DRT1RP(NCOL,NROW,NPER,IUNIT(41),IOUTG,OBSNAM,NHT,JT,
     &                    IX(LCIBT),IX(LCNQOB),IX(LCNQCL),
     &                    IX(LCIQOB),X(LCQCLS),IERR,X(LCHOBS),X(LCTOFF),
     &                    X(LCWTQ),IOWTQ,IPRN,NDMH,NSTP,PERLEN,
     &                    TSMULT,ISSFLG,ITRSS,NQAR,NQCAR,
     &                    NQTAR,IQ1,NQT1,NDD,IUNIT(40),NQDT,NQTDT,NT,
     &                    NC,IX(LCIPLO),NAMES,ND,IPR,MPR,IOWTQDT,
     &                    X(LCOTIM))
C
C-----READ ADVECTIVE-TRANSPORT DATA
      IF (IUNIT(29).GT.0)
     &    CALL OBS1ADV2RP(IOUTG,NROW,NCOL,NLAY,
     &                    X(LCPRST),NPRST,NPTH,IX(LCNPNT),NTT2,NH,NQT,
     &                    OBSNAM,IX(LCICLS),X(LCPOFF),X(LCTT2),
     &                    X(LCHOBS),GX(LCDELR),GX(LCDELC),X(LCWTQ),ND,
     &                    KTDIM,IUNIT(29),NDMH,IOWTQ,GX(LCBOTM),
     &                    NBOTM,IX(LCIPLO),NAMES,IPR,MPR,JT,NPADV,
     &                    INAMLOC,IPFLG,IADVHUF,NHUF,X(LCOTIM),
     &                    PERLEN,NPER,NSTP,ISSFLG,IADVPER,
     &                    TDELC,IMPATHOFF)
C-----CHECK OBSERVATION DATA AGAINST ALLOCATED STORAGE
      IF (IOBS.GT.0) CALL OBS1BAS6CK(NC,ND,NQC,NT,NQT,IOUTG,OBSNAM)
C-----CHECK FOR ERRORS, CALCULATE THE WEIGHT MATRIX AND ITS SQUARE-ROOT
      IF (IPAR.GE.-1)
     &    CALL OBS1BAS6QM(NDMH,X(LCWTQ),X(LCWTQS),DTLWTQ,Z(LCW1),
     &                    Z(LCW2),EV,IOWTQ,IPRN,IOUTG,NDMHAR,OBSALL,
     &                    OUTNAM,ND,NH,X(LCWT))
C
C---------SOLVER PACKAGE
      IF(IUNIT(9).GT.0)
     1    CALL SIP5RPG(NPARM,MXITER,ACCL,HCLOSE,X(LCW),IUNIT(9),IPCALC,
     2                 IPRSIP,IOUTG,IFREFM)
      IF(IUNIT(10).GT.0)
     1    CALL DE45RPG(IUNIT(10),MXITER,NITER,ITMX,ACCL,HCLOSE,IFREQ,
     2                IPRD4,IOUTG,MUTD4)
      IF(IUNIT(11).GT.0)
     1    CALL SOR5RPG(MXITER,ACCL,HCLOSE,IUNIT(11),IPRSOR,IOUTG,IFREFM)
      IF(IUNIT(13).GT.0)
     1    CALL PCG2RPG(MXITER,ITER1,HCLOSE,RCLOSE,NPCOND,NBPOL,RELAX,
     2                IPRPCG,IUNIT(13),IOUTG,MUTPCG,NITER,DAMP,IFREFM)
      IF(IUNIT(14).GT.0)
     1    CALL LMG1RPG(IUNIT(14),MXITER,MXCYC,BCLOSE,DAMP,IOUTAMG,IOUTG,
     2                1,ICG,IADAMP,DUP,DLOW,HCLOSE)
C-----CHECK DATA AND CALCULATE CONVERGENCE CRITERIA FOR SENSITIVITIES
      IF (ISEN.GT.0)
     &    CALL SEN1BAS6CM(JT,IOUTG,IX(LCLN),X(LCB1),IERR,NPER,X(LCHCLO),
     &                    X(LCRCLO),HCLOSE,RCLOSE,IPAR,NPE,NPLIST,
     &                    IX(LCISEN),NSTP,PERLEN,TSMULT,IUNIT(10),
     &                    NOTICECOUNT)
C
C-----READ AND PREPARE FOR PACKAGES WITH NO REWIND
      IF(IUNIT(23).GT.0)
     1    CALL GWF1LPF1RPGD(X(LCHK),X(LCVKA),X(LCVKCB),X(LCHANI),
     2                      X(LCSC1),X(LCSC2),IUNIT(23),ITRSS,NCOL,NROW,
     3                      NLAY,IOUTG,X(LCWETD),NPLPF,WETFCT,IWETIT,
     4                      IHDWET,IX(LCLAYF),GX(LCBOTM),NBOTM,
     5                      GX(LCDELR),GX(LCDELC),1,INAMLOC,
     6                      IX(LCISEN),ISEN,NPLIST,
cgzh mocwt add iwdflg -- this is a change to standard mf2k...right?
     7                      IWDFLG)
      IF(IUNIT(37).GT.0)
     &    CALL GWF1HUF2RPGD(IUNIT(37),NCOL,NROW,NLAY,IOUTG,X(LCWETD),
     &                    WETFCT,IWETIT,IHDWET,IX(LCHGUF),
     &                    1,NHUF,NPHUF,X(LCHUFTHK),
     &                    ITRSS)
      IF(IUNIT(47).GT.0)
     &    CALL GWF1HUF2LVDA1RPGD(IUNIT(47),IOUTG,1,NHUF,NPLVDA,NLAY,
     &                           ISEN)
      IF(IUNIT(53).GT.0)
     &    CALL GWF1HUF2KDEP1RPGD(IUNIT(53),IOUTG,1,NPKDEP,IFKDEP,NROW,
     &                           NCOL,X(LCGS),GX(LCBOTM),NHUF)
C
C-------BEGIN ITERATION LOOP FOR PARAMETER ESTIMATION
      DO 105, KITP = 1,ITMXP
        ITERP = KITP
C
C-------SET SENSITIVITY ARRAYS TO ZERO AND STORE ON DISK OR IN MEMORY
        IF (ISEN.GT.0) CALL SEN1BAS6ZS(IUHEAD,LENXHS,NCOL,NPE,NROW,NLAY,
     &                                 Z(LCSNEW),X(LCSOLD),XHS,
     &                                 X(LCSEND),NTIMES)
C-------LOOP TO HERE WHEN CONVERGENCE HAS BEEN ACHIEVED BY TOL CRITERION
 20     CONTINUE
        ITERPK = ITERPK + 1
        ICNVGP = 1
        IF (IPAR.GT.-3) THEN
C-------IF PARAMETER ESTIMATION HAS CONVERGED, SET ITERPF TO
C       CALCULATE HEAD WITH THE NEW PARAMETERS AND THEN STOP
          IF (IFO.GT.0) THEN
            ITERPF = ITERP
          ENDIF
C---------REWIND INPUT FILES
          IF (ITERPK.GT.1)
     1        CALL PES1BAS6RW(INUNIT,FNAME,CUNIT,IREWND,NIUNIT,IOUT,
     2                        IOUTG,VERSION,IX(LCISEN),ITERP,ITERPF,
     3                        LASTX,NPLIST,ITERPK)
        ENDIF
C
C-------INITIALIZE H AND X ARRAYS, AND UNFLAG OMITTED OBSERVATIONS
        IF (IOBS.GT.0) CALL OBS1BAS6FM(X(LCH),ND,NDAR,NDMH,NDMHAR,
     &                                 X(LCWT),X(LCWTQ))
        IF (ISEN.GT.0 .AND. (ITERPF.EQ.0 .OR. LASTX.GT.0))
     &      CALL OBS1BAS6DR(ND,NPE,X(LCX))
C4------ALLOCATE SPACE IN RX AND IR ARRAYS.
        CALL GWF1BAS6ALP(HEADNG,NPER,TOTIM,NCOL,NROW,NLAY,NODES,INBAS,
     1                   IOUT,IXSEC,ICHFLG,IFREFM,ISUMRX,ISUMIR,ISUMRZ,
     2                   LCIOFL,ISTRT,IAPART)
        IF(IUNIT(1).GT.0)
     1      CALL GWF1BCF6ALP(ISUMRX,LCSC1,LCHY,LCSC2,LCTRPY,ITRSS,ISS,
     2                       IUNIT(1),NCOL,NROW,NLAY,IOUT,IBCFCB,LCWETD,
     3                       IWDFLG,LCCVWD,WETFCT,IWETIT,IHDWET,HDRY,
     4                       IAPART,IFREFM,LAYHDT)
        IF(IUNIT(2).GT.0)
     1      CALL GWF1WEL6ALP(ISUMRX,LCWELL,MXWELL,NWELLS,IUNIT(2),IOUT,
     2                       IWELCB,NWELVL,IWELAL,IFREFM,NPWEL,IPWBEG,
     3                       NNPWEL,NOPRWL)
        IF(IUNIT(3).GT.0)
     1      CALL GWF1DRN6ALP(ISUMRX,LCDRAI,MXDRN,NDRAIN,IUNIT(3),IOUT,
     2                       IDRNCB,NDRNVL,IDRNAL,IFREFM,NPDRN,IDRNPB,
     3                       NNPDRN,NOPRDR)
        IF(IUNIT(4).GT.0)
     1      CALL GWF1RIV6ALP(ISUMRX,LCRIVR,MXRIVR,NRIVER,IUNIT(4),IOUT,
     2                       IRIVCB,NRIVVL,IRIVAL,IFREFM,NPRIV,IRIVPB,
     3                       NNPRIV,NOPRRV)
        IF(IUNIT(5).GT.0)
     1      CALL GWF1EVT6ALP(ISUMRX,ISUMIR,LCIEVT,LCEVTR,LCEXDP,LCSURF,
     2                       NCOL,NROW,NEVTOP,IUNIT(5),IOUT,IEVTCB,
     3                       IFREFM,NPEVT,IEVTPF)
        IF(IUNIT(7).GT.0)
     1      CALL GWF1GHB6ALP(ISUMRX,LCBNDS,MXBND,NBOUND,IUNIT(7),IOUT,
     2                       IGHBCB,NGHBVL,IGHBAL,IFREFM,NPGHB,IGHBPB,
     3                       NNPGHB,NOPRGB)
        IF(IUNIT(8).GT.0)
     1      CALL GWF1RCH6ALP(ISUMRX,ISUMIR,LCIRCH,LCRECH,NRCHOP,NCOL,
     2                       NROW,IUNIT(8),IOUT,IRCHCB,IFREFM,NPRCH,
     3                       IRCHPF)
        IF(IUNIT(16).GT.0)
     1      CALL GWF1FHB1ALP(ISUMRX,ISUMIR,LCFLLC,LCBDTM,LCFLRT,LCBDFV,
     2                  LCBDHV,LCHDLC,LCSBHD,NBDTIM,NFLW,NHED,IUNIT(16),
     3                  IOUT,IFHBCB,NFHBX1,NFHBX2,IFHBD3,IFHBD4,IFHBD5,
     4                  IFHBSS,ITRSS,NHEDDIM,NFLWDIM,NBDHVDIM)
        IF(IUNIT(18).GT.0)
     1      CALL GWF1STR6ALP(ISUMRX,ISUMIR,LCSTRM_,ICSTRM_,MXSTRM,
     2                  NSTREM,IUNIT(18),IOUT,ISTCB1STR6,ISTCB2STR6,
     3                  NSSSTR6,NTRIB,NDIV,ICALC,CONSTSTR6,LCTBAR,
     4                  LCTRIB,LCIVAR_,LCFGAR,NPSTR,ISTRPB)
        IF(IUNIT(19).GT.0)
     1      CALL GWF1IBS6ALP(ISUMRX,LCHC,LCSCE,LCSCV,LCSUB,NCOL,
     2                  NROW,NLAY,IIBSCB,IIBSOC,IUNIT(19),IOUT,IBSDIM,
     &                  IUNIT(54))
        IF(IUNIT(54).GT.0)
     1      CALL GWF1SUB1ALP(NROW,NCOL,NLAY,ITERP,ISUBCB,ISUBOC,AC1,AC2,
     2                  ITMIN,NNDB,NDB,NPZ,NN,NND1,ND1,ND2,IDSAVE,
     3                  IDREST,ISSFLG,NPER,NSTP,NSTPT,IUNIT(54),IOUT,
     4                  IUNIT(9),LCV,ISEN)
        IF(IUNIT(55).GT.0)
     1      CALL GWF1SWT1AL(IUNIT(55),IOUT)
        IF(IUNIT(20).GT.0)
     1      CALL GWF1CHD6ALP(ISUMRX,LCCHDS,NCHDS,MXCHD,IUNIT(20),IOUT,
     2                       NCHDVL,IFREFM,NPCHD,IPCBEG,NNPCHD,NOPRCH)
        IF (IUNIT(17).GT.0)
     &      CALL GWF1RES1ALP(ISUMRX,LCIRES,LCIRSL,LCBRES,LCCRES,LCBBRE,
     &                  LCHRES,LCHRSE,IUNIT(17),IOUT,NRES,IRESCB,NRESOP,
     &                  IRESPT,NPTS,NCOL,NROW,ISUMIR)
        IF (IUNIT(21).GT.0)
     &      CALL GWF1HFB6ALP(IUNIT(21),IOUT,ISUMRX,LCHFB,MXACTFB,NHFBNP,
     &                       NPHFB,MXHFB,IHFB,NOPRHB)
Cdep  Changed SFR1 call to SFR2 call
Cdep   ADDED 3 NEW ARGUMENTS TO END OF SFR2ALP CALL STATEMENT FOR
Cdep   REVISIONS TO THE CALCULATION OF LAKE OUTFLOW 
Cdep   June 6, 2006
        IF(IUNIT(44).GT.0) THEN
            CALL GWF1SFR2ALP(ISUMRX,ISUMIR,ISUMRZ,LCSTRM,ICSTRM,
     2             NSTRM,IUNIT(44),IOUT,ISTCB1,ISTCB2,NSS,CONST,MAXPTS,
     3             DLEAK,LCSEG,ICSEG,LCOTSG,LCXSEC,LCIVAR,LCQSTG,
     4             IUNIT(22),ISTRIN,ISTROT,LCOTFLW,LCDVFLW,IUNIT(15),
     5             NSOL,LSCOUT,LSCONQ,LSCNRN,LSCNTB,LSSLIN,LSCQIN,
     6             LSCGW,ISTGLD,IDSTRT,LSSIN,LSSOUT,LSCOTO,LKACC7,
     7             LCSTAG,LSLAKE,ISTGNW,LSCPPT,LCNSEG,NSFRPAR,
     8             NSEGDIM,LCSFRQ,NSTRMAR,NSSAR,LCSUZDPIT,LCSUZDPST,
     9             LCSUZTHIT,LCSUZTHST,LCSUZSPIT,LCSUZSPST,LCSUZFLIT,
     &             LCSUZFLST,LCSUZFLWT,LCSUZSTOR,LCSDELSTR,LCSUZWDTH,
     &             ICSLTRLIT,ICSLTRLST,ICSITRLIT,ICSITRLST,ICSTRLHLD,
     &             ICSNWAVST,ISFROPT,IUZT,ISUZN,NSTRAIL,NSTOTRL,NUZST,
     &             LCSOLSFLX,ICSLOOP,LCSTHTS,LCSTHTR,LCSTHTI,LCSEPS,
     &             NCOL,NROW,ICELEV,LCDPTH,LCWETP,NSFRSETS,LCSUZSEEP,
     &             LCOLDFLBT,NUMCELL,NUZROW,NUZCOL,LCAVWAT,LCWAT1,
     &             NUMAVE,LCAVDPT,LCUHC,LCSFRUZBD,NSSLK,ISLKOTFLW,
     &             IDLKOTFLW,IDLKSTAGE)
        ELSE
          NSEGDIM=1
        ENDIF
CGWT GWT NOT COMPATIBLE WITH UNSATURATED FLOW IN SFR2
        IF(IUNIT(15).GT.0.AND.IUNIT(44).GT.0) THEN
          IF (ISFROPT.GT.1) then
            WRITE(IOUTS,*) ' GWT AND UNSATURATED FLOW IN SFR2
     * ARE NOT COMPATIBLE.'
            STOP 'GWT-SFR2 ERROR'
          END IF
	  END IF        
Cdep   ADDED 4 NEW ARGUMENTS TO END OF LAK3ALP CALL STATEMENT FOR
Cdep   REVISIONS TO THE CALCULATION OF  LAKE STAGE AND OUTFLOW
Cdep   June 6, 2006
Cdep   Added SURFDEPTh to end of  LAK3ALP CALL STATEMENT for
Cdep   revisions to calculations of ground-water discharge from 
Cdep   a dry lake cell  March 23, 2009
CLAK
        IF(IUNIT(22).GT.0) THEN
                     CALL GWF1LAK3ALP(ISUMRX,ISUMIR,LCCOND,ICLAKE,
     2     MXLKND,LKNODE,LCSTAG,IUNIT(22),IOUT,ILKCB,NLAKES,INTRB,
     3     INDV,LCCNDF,LCLKPR,LCLKEV,ISTGLD,ISTGNW,IICS,IISUB,ISILL,
     4     LCWTDR,IFREFM,NROW,NCOL,NLAY,IBNLK,ILKBL,LKACC1,LKACC2,
     5     LKACC3,LKACC4,LKACC5,LKACC6,LKACC7,LKACC8,LKACC9,LKACC10,
     6     LKACC11,LKDRY,IBTMS,LKNCNT,LKKSUB,LKSADJ,LKFLXI,LKNCNS,LKSVT,
     7     LKJCLS,THETA,LCRNF,ITRSS,NSSITR,SSCNCR,LKSSMN,LKSSMX,LKNCN,
     8     LKDSR,LKCNN,LKCHN,IAREN,IUNIT(44),LSOVOL,NSS,IUNIT(15),
     9     LSLAKE,LSPPT,LSRNF,LSAUG,NSOL,IMSUB,IMSUB1,LSCGWL,LSSLAK,
     *     LSSWIN,LSSWOT,LSSPPT,LSCDRW,LSSRUN,LSGWIN,LSGWOT,LSLKSM,
     *     LSKLK,LSDONE,LSFLOB,LSRTCO,LSCLKO,LSALKI,LSALKO,ISTRIN,
     *     ISTROT,LKLMRR,IDSTRT,LKVI,ISTGLD2,LKCLKI,LKCUM1,LKCUM2,
     *     LKCUM3,LKCUM4,LKCUM5,LKCUM6,LKCUM7,LKCUM8,LKCUM9,NSSAR,
     *     NLAKESAR,IUNIT(46),ISTGITR,LKSEP3,IAREATAB,IDPTHTAB,LCSEG,
     *     ISUMRZ,NSSLK,ISLKOTFLW,IDLKOTFLW,IDLKSTAGE,LCEVAPO,
     *     LCFLWIN,LCFLWIT,LCWITDW,LCGWRAT,SURFDEPTH)
CLAK
        ELSE
           MXLKND=1
        END IF
CLAK
        IF(IUNIT(46).GT.0)
     &      CALL GWF1GAG5ALP(IUNIT(46),ISUMIR,ISUMRX,LSGAGE,NUMGAGE,
     &                       IOUT,IUNIT(44),IUNIT(22),LKACC7,LCSTAG,
     &                       LSLAKE,ICSTRM,LCIVAR)
        IF(IUNIT(39).GT.0)
     &      CALL GWF1ETS1ALP(ISUMRX,ISUMIR,LCIETS,LCETSR,LCETSX,LCETSS,
     &                       NCOL,NROW,NETSOP,IUNIT(39),IOUT,IETSCB,
     &                       IFREFM,NPETS,IETSPF,NETSEG,LCPXDP,LCPETM,
     &                       NSEGAR)
        IF(IUNIT(40).GT.0)
     &      CALL GWF1DRT1ALP(ISUMRX,LCDRTF,MXDRT,NDRTCL,IUNIT(40),IOUT,
     &                       IDRTCB,NDRTVL,IDRTAL,IFREFM,NPDRT,IDRTPB,
     &                       NDRTNP,IDRTFL,NOPRDT)
        IF (IUNIT(43).GT.0)
     &      CALL GWF1HYD1ALP(ISUMRX,LCHYDM,NHYDM,IHYDMUN,HYDNOH,
     &                       IUNIT(43),IOUT)
        IF(IUNIT(51).GT.0)
     1      CALL GWF1DAF1ALP(IERR,IUNIT(52)+1,IUNIT(52),IUNIT(51),IOUT,
     2                       IDAFCB,IDAFBK)
        IF(IUNIT(50).GT.0) THEN
          CALL GWF1MNW1AL(ISUMRZ,LCWEL2,MXWEL2,NWELL2,LCHREF,NODES,
     &                    KSPREF,IUNIT(50),IOUT,IMNWCB,IOWELL2,
     &                    NOMOITER,MNWNAME,FNAME)
C
C         Allocate array for MNW1 site IDs
          IF (ITERPK.EQ.1) THEN
            ALLOCATE (MNWSITE(MXWEL2+1),STAT=ISTAT)
            IF (ISTAT.NE.0) THEN
              WRITE(*,703)ISTAT
  703         FORMAT(1X,'ALLOCATION OF ARRAY MNWSITE FAILED,',/,
     &        ' RETURNED ERROR MESSAGE NUMBER: ',I6)
              CALL USTOP(' ')
            ENDIF
          ENDIF
        ELSE
             mxwel2=1
	       ALLOCATE (MNWSITE(1))
        ENDIF
CGWT  FOR TRANSPORT, ALLOCATE SPACE FOR MNWid array
        IF(IUNIT(15).GT.0) THEN
	    IF(IUNIT(50).GT.0) THEN
	      CALL GWT1MNW1AL(ISUMIR,LSMNWI,mxwel2,IOUTS)
          ELSE
            LSMNWI=1
          END IF
cmnw2
	    IF(IUNIT(57).GT.0) THEN
	      CALL GWT1MNW2AL(ISUMIR,LSMNWI,mnwmax,IOUTS)
          ELSE
            LSMNWI=1
          END IF
        END IF
C
CGWT  FOR MNW OBS WITH TRANSPORT (MNWO PACKAGE), ALLOCATE SPACE 
        IF(IUNIT(15).GT.0) THEN
          CALL GWT1MNWO5AL(ISUMRX,ISUMIR,LSMNWU,LSMNWO,LSIQZE,
     *          MNWOBS,mnwmax,
     *          IOUT,IOUTS,JUNIT(18),IUNIT(57))    
           IF(MNWOBS.GT.0) THEN
            ALLOCATE (MNWOLST(MNWOBS),STAT=ISTAT)
            IF (ISTAT.NE.0) THEN
              WRITE(*,1703) ISTAT
 1703         FORMAT(1X,'ALLOCATION OF ARRAY MNWOLST FAILED,',/,
     &        ' RETURNED ERROR MESSAGE NUMBER: ',I6)
              CALL USTOP(' ')
            ENDIF
           ELSE
            ALLOCATE (MNWOLST(1))
           ENDIF
         END IF
C
CGWT----ALLOCATE SPACE FOR RECHARGE CONCENTRATION (GWT)
C
      IF(IUNIT(8).GT.0.AND.IUNIT(15).GT.0) THEN
         CALL GWT1CRCH5AL(ISUMRX,LSCRCH,NSROW,NSCOL,
     *					JUNIT(1),IOUTS,IOUTG,NRCHOP,IRCHTP)
      ENDIF
C
CGWT----ALLOCATE SPACE FOR CHFB PACKAGE
C
      IF(IUNIT(15).GT.0) THEN
         CALL GWT1CHFB5AL(ISUMRX,ISUMIR,NSROW,NSCOL,NSLAY,
     *                    JUNIT(12),LSHFBL,LSHFBD,MXHFB,LSHBCK,
     *					IOUTS,IOUTG)
      ENDIF
CMOCWT
CGWT----ALLOCATE SPACE FOR CONSTANT-CONCENTRATION BOUNDARIES
cgzh ccbd
      IF(IUNIT(15).GT.0.AND.JUNIT(25).GT.0)
     *   CALL GWT1CCBD1AL(ISUMRX,LSCCBD,JUNIT(25),IOUTS,PTWTON,
     *     NSCOL,NSROW,NSLAY)
C--MNW2
        IF(IUNIT(57).GT.0) THEN
          CALL GWF1MNW2AL(isumrz,lcmnw2,lcmnwn,lcmnwi,lcmnwc,mnwmax,
     &                    IUNIT(57),IOUT,IMNWCB,MNWPRNT,NODTOT,
     &                    IFREFM,NLAY,NMNWVL)
C         Allocate array for MNW2 site IDs
          IF (ITERPK.EQ.1) THEN
            ALLOCATE (WELLID(mnwmax+1),STAT=ISTAT)
            IF (ISTAT.NE.0) THEN
              WRITE(*,713)ISTAT
  713         FORMAT(1X,'ALLOCATION OF ARRAY WELLID FAILED,',/,
     &        ' RETURNED ERROR MESSAGE NUMBER: ',I6)
              CALL USTOP(' ')
            ENDIF
          ENDIF
cgzh 6/4/09 bug fix  allocate WELLID of 1 when MNW2 off
        ELSE
          mnwmax=1
          ALLOCATE (WELLID(mnwmax))     
        ENDIF
C--MNWI
        CALL GWF1MNWIAL(IUNIT(58),IUNIT(57),
     & IOUT,LCMNIO,Wel1flag,QSUMflag,BYNDflag,ISUMRZ,MNWOBS,
     & LSMNWO,MNWMAX)    
C
        IF(IUNIT(58).GT.0) THEN
C         Allocate array for MNW2 site IDs in MNWI routine
          IF (ITERPK.EQ.1) THEN
            ALLOCATE (MNWIID(MNWOBS+1),STAT=ISTAT)
            IF (ISTAT.NE.0) THEN
              WRITE(*,773)ISTAT
  773         FORMAT(1X,'ALLOCATION OF ARRAY MNWIID FAILED,',/,
     &        ' RETURNED ERROR MESSAGE NUMBER: ',I6)
              CALL USTOP(' ')
            ENDIF
          ENDIF
cgzh 6/4/09 bug fix  allocate WELLID of 1 when MNWI off
        ELSE
          MNWOBS=1
          ALLOCATE (MNWIID(MNWOBS))     
        ENDIF
C
C------DYNAMICALLY ALLOCATE RX AND IR ARRAYS FOR PACKAGES THAT DO
C      ALLOCATION EVERY PARAMETER-ESTIMATION ITERATION.  FOR STATIC
C      MEMORY ALLOCATION, THE FOLLOWING IF...THEN BLOCK MUST BE
C      COMMENTED OUT
        IF (ITERPK.EQ.1) THEN
          LENRX = ISUMRX - 1
          IF(LENRX.LE.0) LENRX=1
          LENIR = ISUMIR - 1
          IF(LENIR.LE.0) LENIR=1
          LENRZ = ISUMRZ - 1          
          IF(LENRZ.LE.0) LENRZ=1
          ALLOCATE (RX(LENRX),IR(LENIR),RZ(LENRZ),STAT=ISTAT)
          IF (ISTAT.NE.0) THEN
            WRITE(*,704)ISTAT
  704       FORMAT(1X,'ALLOCATION OF ARRAYS RX, IR, RZ FAILED,',/,
     &      ' RETURNED ERROR MESSAGE NUMBER: ',I6)
            CALL USTOP(' ')
          ENDIF
        ENDIF
C
C5------IF THE ARRAYS ARE NOT BIG ENOUGH THEN STOP.
        CALL MEMCHKR(ISUMRX,ISUMRZ,ISUMIR,LENRX,LENRZ,LENIR,IOUT,IERR,
     &               IERRU)
        IF (IERR.GT.0) CALL PLL1SD(IERR,IERRU,IOUT,IOUTG)
C Initialize SFR GWT arrays
        IF(IUNIT(15).GT.0) THEN
          CALL SMOC5Z(RX(LSCONQ),NSEGDIM,NSOL,1,0.0)
          CALL SMOC5Z(RX(LSCPPT),NSEGDIM,NSOL,1,0.0)
          CALL SMOC5Z(RX(LSCNRN),NSS,NSOL,1,0.0)
        END IF 
C
C6------READ AND PREPARE INFORMATION FOR ENTIRE SIMULATION.
C---------BASIC PACKAGE
        CALL GWF1BAS6RPP(IG(LCIBOU),GZ(LCHNEW),GX(LCSTRT),INBAS,HEADNG,
     1                   NCOL,NROW,NLAY,VBVL,IR(LCIOFL),IUNIT(12),
     2                   IHEDFM,IDDNFM,IHEDUN,IDDNUN,IOUT,IPEROC,ITSOC,
     3                   CHEDFM,CDDNFM,IBDOPT,IXSEC,LBHDSV,LBDDSV,
     4                   IFREFM,IBOUUN,LBBOSV,CBOUFM,HNOFLO,NIUNIT,ITS,
     5                   IAUXSV,RESETDD,RESETDDNEXT)
        IF(IUNIT(1).GT.0)
     1      CALL GWF1BCF6RPP(IG(LCIBOU),GZ(LCHNEW),RX(LCSC1),RX(LCHY),
     2                       GX(LCCR),GX(LCCC),GX(LCCV),GX(LCDELR),
     3                       GX(LCDELC),RX(LCSC2),RX(LCTRPY),IUNIT(1),
     4                       ISS,NCOL,NROW,NLAY,IOUT,RX(LCWETD),IWDFLG,
     5                       RX(LCCVWD),HNOFLO)
C-------SUBSTITUTE AND PREPARE FOR PACKAGES WITH NO REWIND
        IF(IUNIT(23).GT.0)
     1      CALL GWF1LPF1SP(IG(LCIBOU),GZ(LCHNEW),GX(LCCR),GX(LCCC),
     2                      GX(LCCV),GX(LCDELR),GX(LCDELC),GX(LCBOTM),
     3                      X(LCHK),X(LCVKA),X(LCVKCB),X(LCHANI),
     4                      X(LCSC1),X(LCSC2),ITRSS,NCOL,NROW,NLAY,IOUT,
     5                      X(LCWETD),NPLPF,NBOTM,GX(LCRMLT),IG(LCIZON),
     6                      NMLTAR,NZONAR,IX(LCLAYF),GX(LCBUFF),ITERPK,
     7                      HNOFLO)
        IF(IUNIT(37).GT.0)
     1     CALL GWF1HUF2SP(IG(LCIBOU),GZ(LCHNEW),GX(LCCR),GX(LCCC),
     2                      GX(LCCV),GX(LCDELR),GX(LCDELC),GX(LCBOTM),
     3                      X(LCHK),X(LCVKA),X(LCSC1),ITRSS,NCOL,NROW,
     4                      NLAY,IOUT,X(LCWETD),NHUF,NBOTM,GX(LCRMLT),
     5                      IG(LCIZON),NMLTAR,NZONAR,X(LCHUFTHK),
     6                      X(LCHKCC),HDRY,0,0,0,IX(LCHGUF),
     7                      X(LCHUFTMP),IUNIT(47),
     8                      X(LCVDHD),X(LCVDHT),IWETIT,
     9                      IHDWET,WETFCT,X(LCGS),X(LCA9),HNOFLO)
C---------FLOW-SIMULATION OPTIONS
        IF(IUNIT(2).GT.0)
     1      CALL GWF1WEL6RPPD(IUNIT(2),IOUTG,NWELVL,IWELAL,NCOL,NROW,
     2                        NLAY,NPWEL,RX(LCWELL),IPWBEG,MXWELL,
     3                        IFREFM,ITERPK,INAMLOC,NOPRWL)
        IF(IUNIT(3).GT.0)
     1      CALL GWF1DRN6RPPD(IUNIT(3),IOUTG,NDRNVL,IDRNAL,NCOL,NROW,
     2                        NLAY,NPDRN,RX(LCDRAI),IDRNPB,MXDRN,IFREFM,
     &                        ITERPK,INAMLOC,NOPRDR)
        IF(IUNIT(4).GT.0)
     1      CALL GWF1RIV6RPPD(IUNIT(4),IOUTG,NRIVVL,IRIVAL,NCOL,NROW,
     2                        NLAY,NPRIV,RX(LCRIVR),IRIVPB,MXRIVR,
     3                        IFREFM,ITERPK,INAMLOC,NOPRRV)
        IF(IUNIT(5).GT.0)
     &      CALL GWF1EVT6RPPD(IUNIT(5),IOUTG,NPEVT,ITERPK,INAMLOC)
        IF(IUNIT(7).GT.0)
     1      CALL GWF1GHB6RPPD(IUNIT(7),IOUTG,NGHBVL,IGHBAL,NCOL,NROW,
     2                        NLAY,NPGHB,RX(LCBNDS),IGHBPB,MXBND,IFREFM,
     &                        ITERPK,INAMLOC,NOPRGB)
        IF(IUNIT(8).GT.0)
     &      CALL GWF1RCH6RPPD(IUNIT(8),IOUTG,NPRCH,ITERPK,INAMLOC)
        IF(IUNIT(16).GT.0)
     &      CALL GWF1FHB1RPP(IG(LCIBOU),NROW,NCOL,NLAY,IR(LCFLLC),
     &                   RX(LCBDTM),NBDTIM,RX(LCFLRT),NFLW,NHED,
     &                   IR(LCHDLC),RX(LCSBHD),IUNIT(16),IOUT, NFHBX1,
     &                   NFHBX2,IFHBD3,IFHBD5,NHEDDIM,NFLWDIM)
Cdep  Replaced SFR1 with SFR2-- Note either BCF or LPF used for SFR2
Cdep  Added three new variables for computing lake outflow in Lake Package
Crgn  changed calls to sfr2 to support HUF package. 3/26/07
        IF(IUNIT(44).GT.0)THEN
          IF(IUNIT(23).GT.0)THEN
              CALL GWF1SFR2RPP(ITRSS,RX(LCSTRM),IR(ICSTRM),
     2         NSTRM,IUNIT(44),IOUTG,RX(LCSEG),IR(ICSEG),NSS,
     3         IR(LCIVAR),IR(LCOTSG),RX(LCOTFLW),RX(LCDVFLW),MAXPTS,
     4         RX(LCXSEC),RX(LCQSTG),IUNIT(15),RX(LSCONQ),
     5         RX(LSCNRN),RX(LSCPPT),NSOL,IOUTS,NSFRPAR,NSEGDIM,ITERPK,
     &         INAMLOC,IG(LCIBOU),NCOL,NROW,NLAY,RZ(LCSUZDPIT),
     &         RZ(LCSUZDPST),RZ(LCSUZTHIT),RZ(LCSUZTHST),RZ(LCSUZSPIT),
     &         RZ(LCSUZSPST),RZ(LCSUZFLIT),RZ(LCSUZFLST),IR(ICSLTRLIT),
     &         IR(ICSLTRLST),IR(ICSITRLIT),IR(ICSITRLST),IR(ICSTRLHLD),
     &         RZ(LCSUZFLWT),RZ(LCSUZSTOR),RZ(LCSDELSTR),RZ(LCSUZWDTH),
     &         IR(ICSNWAVST),NSTOTRL,ISFROPT,IUZT,ISUZN,NUZST,
     &         RZ(LCSOLSFLX),X(LCSC2),RZ(LCSTHTS),RZ(LCSTHTR),
     &         RZ(LCSTHTI),RZ(LCSEPS),GX(LCDELR),GX(LCDELC),
     &         RZ(LCSUZSEEP),RZ(LCOLDFLBT),NUZROW,NUZCOL,RX(LCUHC),
     &         IUNIT(1),IUNIT(37),X(LCSC1),GX(LCBOTM),NBOTM,GX(LCSTRT),
     &         RX(LCSFRUZBD),ISSFLG(1),ITMP,IRDFLG,IPTFLG,NP,IR(LCNSEG),
     &         IUNIT(22),NSSLK,RZ(ISLKOTFLW),RZ(IDLKOTFLW),
     &         RZ(IDLKSTAGE))
            IF(ISFROPT.EQ.2.OR.ISFROPT.EQ.4) THEN
              CALL GWF1SFR2UHC(IR(ICSTRM),NSTRM,RX(LCUHC),X(LCHK),
     1              X(LCVKA),IG(LCIBOU),NROW,NCOL,NLAY,NUZST,IOUTG)
            END IF
          ELSEIF(IUNIT(1).GT.0)THEN
Cdep  Added three new variables for computing lake outflow in Lake Package
            CALL GWF1SFR2RPP(ITRSS, RX(LCSTRM),IR(ICSTRM),
     2       NSTRM,IUNIT(44),IOUTG,RX(LCSEG),IR(ICSEG),NSS,
     3       IR(LCIVAR),IR(LCOTSG),RX(LCOTFLW),RX(LCDVFLW),MAXPTS,
     4       RX(LCXSEC),RX(LCQSTG),IUNIT(15),RX(LSCONQ),
     5       RX(LSCNRN),RX(LSCPPT),NSOL,IOUTS,NSFRPAR,NSEGDIM,ITERPK,
     &       INAMLOC,IG(LCIBOU),NCOL,NROW,NLAY,RZ(LCSUZDPIT),
     &       RZ(LCSUZDPST),RZ(LCSUZTHIT),RZ(LCSUZTHST),RZ(LCSUZSPIT),
     &       RZ(LCSUZSPST),RZ(LCSUZFLIT),RZ(LCSUZFLST),IR(ICSLTRLIT),
     &       IR(ICSLTRLST),IR(ICSITRLIT),IR(ICSITRLST),IR(ICSTRLHLD),
     &       RZ(LCSUZFLWT),RZ(LCSUZSTOR),RZ(LCSDELSTR),RZ(LCSUZWDTH),
     &       IR(ICSNWAVST),NSTOTRL,ISFROPT,IUZT,ISUZN,NUZST,
     &       RZ(LCSOLSFLX),RX(LCSC2),RZ(LCSTHTS),RZ(LCSTHTR),
     &       RZ(LCSTHTI),RZ(LCSEPS),GX(LCDELR),GX(LCDELC),
     &       RZ(LCSUZSEEP),RZ(LCOLDFLBT),NUZROW,NUZCOL,RX(LCUHC),
     &       IUNIT(1),IUNIT(37),RX(LCSC1),GX(LCBOTM),NBOTM,GX(LCSTRT),
     &       RX(LCSFRUZBD),ISSFLG(1),ITMP,IRDFLG,IPTFLG,NP,IR(LCNSEG),
     &       IUNIT(22),NSSLK,RZ(ISLKOTFLW),RZ(IDLKOTFLW),RZ(IDLKSTAGE))
          ELSEIF(IUNIT(37).GT.0)THEN
            CALL GWF1SFR2RPP(ITRSS, RX(LCSTRM),IR(ICSTRM),
     2         NSTRM,IUNIT(44),IOUTG,RX(LCSEG),IR(ICSEG),NSS,
     3         IR(LCIVAR),IR(LCOTSG),RX(LCOTFLW),RX(LCDVFLW),MAXPTS,
     4         RX(LCXSEC),RX(LCQSTG),IUNIT(15),RX(LSCONQ),
     5         RX(LSCNRN),RX(LSCPPT),NSOL,IOUTS,NSFRPAR,NSEGDIM,ITERPK,
     &         INAMLOC,IG(LCIBOU),NCOL,NROW,NLAY,RZ(LCSUZDPIT),
     &         RZ(LCSUZDPST),RZ(LCSUZTHIT),RZ(LCSUZTHST),RZ(LCSUZSPIT),
     &         RZ(LCSUZSPST),RZ(LCSUZFLIT),RZ(LCSUZFLST),IR(ICSLTRLIT),
     &         IR(ICSLTRLST),IR(ICSITRLIT),IR(ICSITRLST),IR(ICSTRLHLD),
     &         RZ(LCSUZFLWT),RZ(LCSUZSTOR),RZ(LCSDELSTR),RZ(LCSUZWDTH),
     &         IR(ICSNWAVST),NSTOTRL,ISFROPT,IUZT,ISUZN,NUZST,
     &         RZ(LCSOLSFLX),X(LCSC1),RZ(LCSTHTS),RZ(LCSTHTR),
     &         RZ(LCSTHTI),RZ(LCSEPS),GX(LCDELR),GX(LCDELC),
     &         RZ(LCSUZSEEP),RZ(LCOLDFLBT),NUZROW,NUZCOL,RX(LCUHC),
     &         IUNIT(1),IUNIT(37),X(LCSC1),GX(LCBOTM),NBOTM,GX(LCSTRT),
     &         RX(LCSFRUZBD),ISSFLG(1),ITMP,IRDFLG,IPTFLG,NP,IR(LCNSEG),
     &         IUNIT(22),NSSLK,RZ(ISLKOTFLW),RZ(IDLKOTFLW),
     &         RZ(IDLKSTAGE))
          ELSE 
            WRITE(*,202)
  202       FORMAT(1X,'ERROR--SFR2 ONLY CAN BE USED WITH BCF, LPF,',
     &             ' OR HUF FLOW PACKAGES. PROGRAM STOPPING')
            CALL USTOP(' ')
          END IF       
        END IF
Cdep  End change from SFR1 to SFR2
        IF(IUNIT(18).GT.0)
     1      CALL GWF1STR6RPPD(IUNIT(18),IOUTG,NCOL,NROW,NLAY,NPSTR,
     2                  RX(LCSTRM_),IR(ICSTRM_),ISTRPB,MXSTRM,ITERPK,
     &                  INAMLOC)
        IF(IUNIT(19).GT.0)
     1      CALL GWF1IBS6RPP(GX(LCDELR),GX(LCDELC),GZ(LCHNEW),RX(LCHC),
     2                  RX(LCSCE),RX(LCSCV),RX(LCSUB),NCOL,NROW,
     3                  NLAY,NODES,IIBSOC,ISUBFM,ICOMFM,IHCFM,
     4                  ISUBUN,ICOMUN,IHCUN,IUNIT(19),IOUT,IBSDIM)
        IF(IUNIT(54).GT.0)
     1      CALL GWF1SUB1RPP(GX(LCDELR),GX(LCDELC),GZ(LCHNEW),
     2                  GX(LCBUFF),NCOL,NROW,NLAY,NODES,NPER,NSTP,
     3                  ISUBOC,NND1,ND1,ND2,NDB,NNDB,NPZ,NN,IDSAVE,
     4                  IDREST,NSTPT,IUNIT(54),IOUT)
        IF(IUNIT(55).GT.0)
     1      CALL GWF1SWT1AR(NCOL,NROW,NLAY,NPER,ISSFLG,NSTP,ITERP,
     2                  ISEN,LAYCBD,IG(LCIBOU),GZ(LCHNEW),GX(LCBOTM),
     3                  GX(LCBUFF),GX(LCDELR),GX(LCDELC),IUNIT(55),IOUT)
        IF(IUNIT(20).GT.0)
     1      CALL GWF1CHD6RPPD(IUNIT(20),IOUTG,NCHDVL,NCOL,NROW,NLAY,
     2                        NPCHD,RX(LCCHDS),IPCBEG,MXCHD,IFREFM,
     &                        ITERPK,INAMLOC,NOPRCH)
C
        IF (IUNIT(21).GT.0)
     &      CALL GWF1HFB6RPPA(GX(LCBOTM),GX(LCCR),GX(LCCC),GX(LCDELR),
     &                        GX(LCDELC),RX(LCHFB),IUNIT(21),MXACTFB,
     &                        NBOTM,NCOL,NROW,NLAY,NODES,NHFBNP,NHFB,
     &                        NPHFB,IOUT,IOUTG,ITERPK,MXHFB,IHFB,LAYHDT,
     &                        INAMLOC,NOPRHB)
c chfb
CGWT----READ HFB DISPERSIVITIES AND DIFFUSION COEFFS FOR GW TRANSPORT
CGWT----DEFINE HFB LOCATION ARRAY FOR GWT DISPERSION CALCULATION
        IF (IUNIT(15).GT.0.AND.IUNIT(21).GT.0) 
     *	  CALL GWT1HFB6RP(IR(LSHFBL),IG(LCIBOU),RX(LSHFBD),
     *					  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS,    
     *	 		  RX(LCHFB),NHFB,MXHFB,JUNIT(12),IR(LSHBCK),MOCTYPE)
C
        IF(IUNIT(39).GT.0)
     &      CALL GWF1ETS1RPPD(IUNIT(39),IOUTG,NPETS,ITERPK,INAMLOC)
        IF(IUNIT(40).GT.0)
     &      CALL GWF1DRT1RPPD(IUNIT(40),IOUTG,NDRTVL,IDRTAL,NCOL,NROW,
     &                        NLAY,NPDRT,RX(LCDRTF),IDRTPB,MXDRT,IFREFM,
     &                        ITERPK,IDRTFL,INAMLOC,NOPRDT)
CLAK
C  REVISED IF STATEMENT
        IF(IUNIT(46).GT.0.AND.NUMGAGE.GT.0)
     &      CALL GWF1GAG5RPP(IR(LSGAGE),NUMGAGE,IOUT,IUNIT(46))
C
C-------CHECK THAT PARAMETER DEFINITIONS ARE COMPLETE
        IF (ITERPK.EQ.1) CALL GLO1BAS6CK(IOUTG,ISEN,NPLIST)
        IF ((ISEN.GT.0 .OR. IBEFLG.EQ.2) .AND. ITERPK.EQ.1)
     &      CALL SEN1BAS6CP(IOUTG,NPLIST,ISENSU,CHEDFM)
        IF (IPES.GT.0)
     &      CALL PES1BAS6CK(X(LCBL),X(LCBU),IX(LCISEN),IOUB,IOUTG,
     &                      IX(LCIPNG),IX(LCLN),NPNG,NPLIST,NPNGAR,
     &                      ITERPK,FAC,FCONV,AP,ADMX,TOL,LAYHDT,NLAY,
     &                      X(LCBSCA),X(LCPARE),ITMXP)
        IF (IUNIT(43).GT.0)
     &      CALL GWF1HYD1RPP(RX(LCHYDM),GX(LCSTRT),NHYDM,NUMH,
     &                       GX(LCDELR),GX(LCDELC),NCOL,NROW,NLAY,
     &                       LCHNEW,LCIBOU,IUNIT(43),IOUT)
        IF(IUNIT(43).GT.0 .AND. IUNIT(19).GT.0)
     &      CALL GWF1HYD1IBS2RPP(RX(LCHYDM),NHYDM,NUMH,GX(LCDELR),
     &                       GX(LCDELC),NCOL,NROW,NLAY,LCIBOU,LCSUB,
     &                       LCHC,IUNIT(43),IOUT)
CGWT----READ AND PREPARE INFORMATION FOR GWT (SOLUTE TRANSPORT)
C       CHECK FOR CONSISTENCY BETWEEN FLOW AND TRANSPORT INPUT DATA
C
C this call depends on BCF or LPF
      IF(IUNIT(15).GT.0) THEN
C for BCF, SC1 SC2 and CHY are in RX
	  IF(IUNIT(1).GT.0) CALL GWT1CKRP6(IG(LCIBOU),GX(LCDELR),
     *			GX(LCDELC),X(LSCB),NTFACE,IX(LSLB),
     *			X(LSXFOR),X(LSXBAC),X(LSYFOR),X(LSYBAC),
C RX calls
     *			RX(LCSC1),RX(LCSC2),X(LSRF),GX(LCBOTM),NBOTM,
     *			X(LSCONC),X(LSCINT),X(LSCHBC),X(LSCINF),
     *			X(LSTHCK),X(LSALNG),X(LSATH),X(LSATV),
     *			X(LSPOR),X(LSPNWC),X(LSPNWR),X(LSPNWL),IX(LSIGNP),
     *			FDTMTH,EPSSLV,IDIREC,NCXIT,MAXIT,
     *			SBVL,NEWPTS,NSCOL,NSROW,NSLAY,NODESS,
     *			NCOL,NROW,NLAY,IOUTS,INMOC,
     *			iunit,lcwell,nwelvl,mxwell,lcdrai,ndrnvl,mxdrn,
     *			lcrivr,nrivvl,mxrivr,lcbnds,nghbvl,mxbnd,
     *            lcdrtf,ndrtvl,mxdrt,ICLAKE,
     *			lcievt,lsevtf,lscrch,lcirch,lcrech,
     *            LCHFB,lshfbl,lshfbd,LSCCBD,
     *			LCFLLC,LCBDFV,LCHDLC,LCBDHV,
     *			ICSTRT,JRF,NODISP,NPNTVL,NPNTCL,NPNTDL,NPNTPL,
     *			IVELFM,ICONFM,IDSPFM,JUNIT(24),cnoflw,
     *			ICONLY,IFXHED,ISS,
     *            NCINFL,INTRPL,IDIM,NPTPND,MOCTYPE,
C RX call
     *			GX(LCCC),RX(LCHY),
     *			IX(LSOBSW),NUMOBS,JUNIT,IOBSFL,
     *			SRCDCY,IDUM,IDUM,IDUM,NPER,DUM,TOTIM,
     *			GZ(LCHNEW),0.0,GX(LCBUFF),
     *			X(LSDKZO),X(LSDKFO),X(LSDKZS),
     *			X(LSDKFS),IDKRF,IDKZO,IDKFO,IDKZS,IDKFS,
     *			x(LSDPCON),x(LSDPRAT),x(LSDPPOR),
     *			X(LSDPZO),X(LSDPFO),IDPZO,IDPFO,IDPPS,IWDFLG,LCSEG,
     *            LCSTRM,ICSTRM,NSTRM,NSS,LCOTFLW,LCIVAR,ICSEG,LCOTSG,
     *            LWRT,LSLAKE,LSPPT,LSAUG,LSRNF,LSCGWL,LSSLAK,LSSWIN,
     *            LSCOUT,LSCONQ,LSCNRN,LSCNTB,LSCQIN,LSCGW,LSSLIN,
     *            LSSWOT,LSSPPT,LSCDRW,LSSRUN,LSGWIN,LSGWOT,LKACC7,
     *            LSOVOL,LSDONE,LSKLK,LSLKSM,LSFLOB,LSRTCO,LSALKI,
     *            LSALKO,LKNODE,LCCNDF,IMSUB,IMSUB1,
     *            LCWTDR,LKACC2,LCRNF,LSCLKO,LSCPPT,
     *            INTRB,INDV,ISTRIN,ISTROT,LCSFRQ,
     *            ISTGLD,ISTGNW,LAYHDT,NIUNIT,
cgage
     *            LKACC1,LKACC8,LKACC9,LKACC5,LKACC6,LKFLXI,LKCNN,
     *            LCSTAG,LKVI,LKCLKI,ISTGLD2,
cgzh varpt
     *            LSPCOR,LSPROR,LSPLOR,REMCRIT,GENCRIT,IRAND,ISRCFIX,
     *            ISEEDPT,ISEEDBD,
cgzh ssfix
     *            ISSFLG,MULTSS,
     *            IX(LSBFCH),NCHNDS,
     *            NZERO,NSKIP,
cgzh cbdy
     *            X(LSCINA),X(LSCINB),X(LSCINXY),IABOVE,IBELOW,
cgzh mnw
     *            LCWEL2,NEVTOP,LCMNW2,LCMNWN,
cgzh zinn%
     *            X(LSCINI),X(LSCTNM),
cea
     *            SRCAGE,
csfr
     1            RX(LSCOUT),RX(LSCONQ),RX(LSCNRN),RX(LSCNTB),
     2            RX(LSCQIN),RX(LSCGW),RX(LSCPPT),
	3            NSOL,NSEGDIM) 
c
C for LPF, SC1 SC2 and CHY are in X
        IF(IUNIT(23).GT.0.OR.IUNIT(37).GT.0) THEN
C          for HUF, set SC2=SC1
	     IF (IUNIT(37).GT.0) LCSC2=LCSC1
           CALL GWT1CKRP6(IG(LCIBOU),GX(LCDELR),
     *			GX(LCDELC),X(LSCB),NTFACE,IX(LSLB),
     *			X(LSXFOR),X(LSXBAC),X(LSYFOR),X(LSYBAC),
C X calls
     *			X(LCSC1),X(LCSC2),X(LSRF),GX(LCBOTM),NBOTM,
     *			X(LSCONC),X(LSCINT),X(LSCHBC),X(LSCINF),
     *			X(LSTHCK),X(LSALNG),X(LSATH),X(LSATV),
     *			X(LSPOR),X(LSPNWC),X(LSPNWR),X(LSPNWL),IX(LSIGNP),
     *			FDTMTH,EPSSLV,IDIREC,NCXIT,MAXIT,
     *			SBVL,NEWPTS,NSCOL,NSROW,NSLAY,NODESS,
     *			NCOL,NROW,NLAY,IOUTS,INMOC,
     *			iunit,lcwell,nwelvl,mxwell,lcdrai,ndrnvl,mxdrn,
     *			lcrivr,nrivvl,mxrivr,lcbnds,nghbvl,mxbnd,
     *            lcdrtf,ndrtvl,mxdrt,ICLAKE,
     *			lcievt,lsevtf,lscrch,lcirch,lcrech,
     *            LCHFB,lshfbl,lshfbd,LSCCBD,
     *			LCFLLC,LCBDFV,LCHDLC,LCBDHV,
     *			ICSTRT,JRF,NODISP,NPNTVL,NPNTCL,NPNTDL,NPNTPL,
     *			IVELFM,ICONFM,IDSPFM,JUNIT(24),cnoflw,
     *			ICONLY,IFXHED,ISS,
     *            NCINFL,INTRPL,IDIM,NPTPND,MOCTYPE,
C X call
     *			GX(LCCC),X(LCHK), 
     *			IX(LSOBSW),NUMOBS,JUNIT,IOBSFL,
     *			SRCDCY,IDUM,IDUM,IDUM,NPER,DUM,TOTIM,
     *			GZ(LCHNEW),0.0,GX(LCBUFF),
     *			X(LSDKZO),X(LSDKFO),X(LSDKZS),
     *			X(LSDKFS),IDKRF,IDKZO,IDKFO,IDKZS,IDKFS,
     *			x(LSDPCON),x(LSDPRAT),x(LSDPPOR),
     *			X(LSDPZO),X(LSDPFO),IDPZO,IDPFO,IDPPS,IWDFLG,LCSEG,
     *            LCSTRM,ICSTRM,NSTRM,NSS,LCOTFLW,LCIVAR,ICSEG,LCOTSG,
     *            LWRT,LSLAKE,LSPPT,LSAUG,LSRNF,LSCGWL,LSSLAK,LSSWIN,
     *            LSCOUT,LSCONQ,LSCNRN,LSCNTB,LSCQIN,LSCGW,LSSLIN,
     *            LSSWOT,LSSPPT,LSCDRW,LSSRUN,LSGWIN,LSGWOT,LKACC7,
     *            LSOVOL,LSDONE,LSKLK,LSLKSM,LSFLOB,LSRTCO,LSALKI,
     *            LSALKO,LKNODE,LCCNDF,IMSUB,IMSUB1,
     *            LCWTDR,LKACC2,LCRNF,LSCLKO,LSCPPT,
     *            INTRB,INDV,ISTRIN,ISTROT,LCSFRQ,
     *            ISTGLD,ISTGNW,LAYHDT,NIUNIT,
cgage
     *            LKACC1,LKACC8,LKACC9,LKACC5,LKACC6,LKFLXI,LKCNN,
     *            LCSTAG,LKVI,LKCLKI,ISTGLD2,
cgzh varpt
     *            LSPCOR,LSPROR,LSPLOR,REMCRIT,GENCRIT,IRAND,ISRCFIX,
     *            ISEEDPT,ISEEDBD,
cgzh ssfix
     *            ISSFLG,MULTSS,
     *            IX(LSBFCH),NCHNDS,
     *            NZERO,NSKIP,
cgzh cbdy
     *            X(LSCINA),X(LSCINB),X(LSCINXY),IABOVE,IBELOW,
cgzh mnw
     *            LCWEL2,NEVTOP,LCMNW2,LCMNWN,
cgzh zinn%
     *            X(LSCINI),X(LSCTNM),
cea
     *            SRCAGE,
csfr 7/9/07 initialize routing variables
     1            RX(LSCOUT),RX(LSCONQ),RX(LSCNRN),RX(LSCNTB),
     2            RX(LSCQIN),RX(LSCGW),RX(LSCPPT),
	3            NSOL,NSEGDIM) 
        END IF
	END IF
cgzh debug output
c      IF(PTWTON.EQ.1) THEN
c        WRITE(*,*) 
c        WRITE(*,*) 'USING MOCWT'
c        WRITE(*,*) 
c      END IF
C
CGWT
cgzh SSTR
c flag for GWT on, used in SSTR package runs
        IGWTON=0
C7------SIMULATE EACH STRESS PERIOD.
        DO 100 KPER = 1, NPER
          KKPER = KPER
cgzh SSTR 
cgzh check to see if GWT starts this stress period, activate flag 
cgzh if SSTR off, IPERGWT set to 1 in GWT1SSTR5DF
          IF(IUNIT(15).GT.0.AND.IPERGWT.EQ.KKPER) THEN
            IGWTON=1
          END IF
C
          CALL GWF1BAS6ST(NSTP(KKPER),DELT,TSMULT(KKPER),PERTIM,KKPER,
     &                    IOUT,PERLEN(KKPER))
          IF(IUNIT(19).GT.0)
     1        CALL GWF1IBS6ST(ISSFLG,KKPER,GZ(LCHNEW),RX(LCHC),NCOL,
     2                        NROW,NLAY,IBSDIM,IOUT)
          IF(IUNIT(54).GT.0)
     1        CALL GWF1SUB1ST(GZ(LCHNEW),NNDB,NDB,ISSFLG,NROW,NCOL,
     1                        NODES,NPER,KPER,NN)
        IF(IUNIT(55).GT.0)
     1      CALL GWF1SWT1ST(IG(LCIBOU),GZ(LCHNEW),GX(LCBOTM),
     2                      GX(LCBUFF),GX(LCDELR),GX(LCDELC),ISSFLG,
     3                      NROW,NCOL,NLAY,NPER,KKPER,IOUT)
C
C7B-----READ AND PREPARE INFORMATION FOR STRESS PERIOD.
C----------READ USING PACKAGE READ AND PREPARE MODULES.
          IF(IUNIT(2).GT.0)
     &        CALL GWF1WEL6RPSS(RX(LCWELL),NWELLS,MXWELL,IUNIT(2),IOUT,
     1                          NWELVL,IWELAL,IFREFM,NCOL,NROW,NLAY,
     2                          NNPWEL,NPWEL,IPWBEG,NOPRWL)
          IF(IUNIT(3).GT.0)
     &        CALL GWF1DRN6RPSS(RX(LCDRAI),NDRAIN,MXDRN,IUNIT(3),IOUT,
     1                          NDRNVL,IDRNAL,IFREFM,NCOL,NROW,NLAY,
     2                          NNPDRN,NPDRN,IDRNPB,NOPRDR)
          IF(IUNIT(4).GT.0)
     &        CALL GWF1RIV6RPSS(RX(LCRIVR),NRIVER,MXRIVR,IUNIT(4),IOUT,
     1                          NRIVVL,IRIVAL,IFREFM,NCOL,NROW,NLAY,
     2                          NNPRIV,NPRIV,IRIVPB,NOPRRV)
          IF(IUNIT(5).GT.0)
     &        CALL GWF1EVT6RPSS(NEVTOP,IR(LCIEVT),RX(LCEVTR),RX(LCEXDP),
     1                          RX(LCSURF),GX(LCDELR),GX(LCDELC),NCOL,
     2                          NROW,IUNIT(5),IOUT,IFREFM,NPEVT,
     3                          GX(LCRMLT),IG(LCIZON),NMLTAR,NZONAR,
     &                          IEVTPF)
          IF(IUNIT(7).GT.0)
     &        CALL GWF1GHB6RPSS(RX(LCBNDS),NBOUND,MXBND,IUNIT(7),IOUT,
     1                          NGHBVL,IGHBAL,IFREFM,NCOL,NROW,NLAY,
     2                          NNPGHB,NPGHB,IGHBPB,NOPRGB)
          IF(IUNIT(8).GT.0)
     &        CALL GWF1RCH6RPSS(NRCHOP,IR(LCIRCH),RX(LCRECH),GX(LCDELR),
     1                          GX(LCDELC),NROW,NCOL,IUNIT(8),IOUT,
     2                          IFREFM,NPRCH,GX(LCRMLT),IG(LCIZON),
     &                          NMLTAR,NZONAR,IRCHPF)
          IF (IUNIT(17).GT.0)
     &        CALL GWF1RES1RPS(IR(LCIRES),IR(LCIRSL),RX(LCBRES),
     &                     RX(LCCRES),RX(LCBBRE),RX(LCHRSE),IG(LCIBOU),
     &                     GX(LCDELR),GX(LCDELC),NRES,NRESOP,NPTS,NCOL,
     &                     NROW,NLAY,PERLEN(KKPER),DELT,NSTP(KKPER),
     &                     TSMULT(KKPER),IUNIT(17),IOUT,KKPER)
          IF (IUNIT(18).GT.0)
     &        CALL GWF1STR6RPSS(RX(LCSTRM_),IR(ICSTRM_),NSTREM,MXSTRM,
     &                    IUNIT(18),IOUT,IR(LCTBAR),NDIV,NSSSTR6,NTRIB,
     &                    IR(LCIVAR_),ICALC,IPTFLG,NCOL,NROW,NLAY,
     &                    NPSTR,ISTRPB)
          IF(IUNIT(20).GT.0)
     &        CALL GWF1CHD6RPSS(RX(LCCHDS),NCHDS,MXCHD,IG(LCIBOU),NCOL,
     &                          NROW,NLAY,IUNIT(20),IOUT,NCHDVL,IFREFM,
     &                          NNPCHD,NPCHD,IPCBEG,NOPRCH)
CGWT----SET CONCENTRATION AT CHD CELLS
          IF(IUNIT(20).GT.0.AND.IUNIT(15).GT.0)
     &        CALL GWT1CHD6RP(RX(LCCHDS),NCHDS,MXCHD,IG(LCIBOU),NCOL,
     &                        NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS,NCHDVL,
     &                        NNPCHD,X(LSCHBC),KKPER,IPERGWT,
     &                        ICONLY,NOPRCH)
Cdep  Changed SFR1 to SFR2
Cdep  Added three new variables for computing lake outflow in Lake Package
          IF(IUNIT(44).GT.0)
     1        CALL GWF1SFR2RPS(RX(LCSTRM),IR(ICSTRM),KKPER,NSTRM,
     2            IUNIT(44),IOUT,RX(LCSEG),IR(ICSEG),IR(LCNSEG),NSS,
     3            IR(LCIVAR),IR(LCOTSG),MAXPTS,IPTFLG,RX(LCXSEC),
     4            RX(LCQSTG),IUNIT(15),RX(LSCONQ),RX(LSCNRN),
     5            RX(LSCPPT),NSOL,IOUTS,NSFRPAR,NSEGDIM,RZ(LCSUZTHST),
     6            RZ(LCSUZFLST),RZ(LCSUZDPST),RZ(LCSUZSPST),
     7            RZ(LCSOLSFLX),RZ(LCSTHTI),RZ(LCSTHTR),RZ(LCSTHTS),
     3            RZ(LCSEPS),ISFROPT,IUZT,ISUZN,NCOL,NROW,NLAY,NUZST,
     &            NSTOTRL,RZ(LCSUZSEEP),RZ(LCSUZSTOR),RX(LCUHC),
     &            RX(LCDPTH),RZ(LCWETP),NUZROW,NUZCOL,IR(ICSNWAVST),
     &            GZ(LCHNEW),IG(LCIBOU),ISSFLG,NPER,ITMP,IRDFLG,NP,
     &            GX(LCBOTM),NBOTM,CONST,IUNIT(22),NLAKESAR,NSSLK,
     &            RZ(ISLKOTFLW),RZ(IDLKOTFLW),RZ(IDLKSTAGE))
Cdep  End of change
          IF (IUNIT(44).GT.0 .AND. ISEN.GT.0)
     &        CALL GWF1SFR2SEN(IOUTG,NPLIST,NSEGDIM,IR(ICSEG),
     &                           IX(LCISEN))
Cdep  Added two new arguments to end of LAK3RPS call statement.
CLAK
          IF(IUNIT(22).GT.0) THEN
            CALL GWF1LAK3RPS(IR(ICLAKE),LKNODE,MXLKND,
     1        IUNIT(22),IOUT,NLAKES,RX(LCSTAG),RX(LCLKPR),RX(LCLKEV),
     2        RX(LCCOND),NTRB,NDV,IR(INTRB),IR(INDV),KKPER,
     3        GX(LCDELR),GX(LCDELC),
     4        NCOL,NROW,NLAY,IR(IICS),RX(LKACC7),GX(LCBOTM),NBOTM,
     5        IR(IISUB),RX(ISILL),ICMX,NCLS,RX(LCWTDR),LWRT,IFREFM,
     6        IR(IBNLK),RX(ILKBL),IR(IBNLK),RX(ILKBL),NODES,
     7        RX(IBTMS),RX(LCRNF),RX(IAREN),IUNIT(44),NSS,
     8        IUNIT(15),RX(LSLAKE),RX(LSAUG),RX(LSPPT),RX(LSRNF),
     9        NSOL,IOUTS,RX(LKSSMN),RX(LKSSMX),ISSFLG(KKPER),RX(LKVI),
     *        RX(LKCLKI),RX(LKCUM1),RX(LKCUM2),RX(LKCUM3),RX(LKCUM4),
     &        RX(LKCUM5),RX(LKCUM6),RX(LKCUM7),RX(LKCUM8),
     &        RX(LKCUM9),IG(LCIBOU),NSSLK,RX(IAREATAB),RX(IDPTHTAB))
            IF (IUNIT(1).GT.0) THEN
              CALL GWF1LAK3BCF6RPS(IOUT,RX(LCCOND),IR(IBNLK),
     1             IR(ICLAKE),RX(LCCNDF),GX(LCDELR),GX(LCDELC),
     2             RX(LCHY),RX(LCTRPY),LAYHDT,MXLKND,NCOL,NROW,NLAY,
     3             LKNODE,IWDFLG,RX(LCCVWD))
            ELSE IF (IUNIT(23).GT.0) THEN
              CALL GWF1LAK3LPF1RPS(IOUT,RX(LCCOND),IR(IBNLK),
     1             IR(ICLAKE),RX(LCCNDF),GX(LCDELR),GX(LCDELC),
     2             X(LCHK),X(LCHANI),LAYHDT,MXLKND,NCOL,NROW,NLAY,
     3             LKNODE,X(LCVKA),X(LCVKCB),GX(LCBOTM),NBOTM)
            ELSE IF(IUNIT(37).GT.0) THEN
              CALL GWF1LAK3HUF1RPS(IOUT,RX(LCCOND),IR(IBNLK),
     1             IR(ICLAKE),RX(LCCNDF),GX(LCDELR),GX(LCDELC),
     2             X(LCHK),X(LCHKCC),LAYHDT,MXLKND,NCOL,NROW,NLAY,
     3             LKNODE,X(LCVKA),GX(LCBOTM),NBOTM)
            ELSE
              WRITE(IOUT,*) 'LAK Package requires BCF, LPF, or HUF'
              CALL USTOP(' ')
            END IF
            IF (IUNIT(44).GT.0)
     &          CALL GWF1LAK3SFR2RPS(NTRB,NDV,NLAKES,IR(INTRB),IR(INDV),
     &                  NSS,NSSLK,IR(LCIVAR),IR(LCOTSG),RX(LCSEG),
     &                  IR(ICSEG),IOUT,NODES,GX(LCBUFF))
          END IF
CLAK
          IF (IUNIT(46).GT.0.AND.(NUMGAGE.GT.0).AND.KKPER.EQ.1)
     &        CALL GWF1GAG5I(IR(LSGAGE),NUMGAGE,IOUT,IUNIT(15),
     &                       RX(LCSTAG),RX(LSLAKE),NLAKES,IR(ICSTRM),
     &                       NSTRM,IR(LCIVAR),DUM,NSOL,RX(LKACC7),
     &                       NLAKESAR,NSTRMAR,NSSAR)
          IF(IUNIT(39).GT.0)
     &        CALL GWF1ETS1RPSS(NETSOP,IR(LCIETS),RX(LCETSR),RX(LCETSX),
     &                          RX(LCETSS),GX(LCDELR),GX(LCDELC),NCOL,
     &                          NROW,IUNIT(39),IOUT,IFREFM,NPETS,
     &                          GX(LCRMLT),IG(LCIZON),NMLTAR,NZONAR,
     &                          IETSPF,NETSEG,RX(LCPXDP),RX(LCPETM),
     &                          NSEGAR)
          IF(IUNIT(40).GT.0)
     &        CALL GWF1DRT1RPSS(RX(LCDRTF),NDRTCL,MXDRT,IUNIT(40),IOUT,
     &                          NDRTVL,IDRTAL,IFREFM,NCOL,NROW,NLAY,
     &                          NDRTNP,NPDRT,IDRTPB,IDRTFL,NRFLOW,
     &                          NOPRDT)
          IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0 .AND. KPER.EQ.1)
     &        CALL GWF1HYD1STR6RPS(IR(ICSTRM_),RX(LCHYDM),NHYDM,NUMH,
     &                         GX(LCDELR),GX(LCDELC),NCOL,NROW,NLAY,
     &                         LCIBOU,LCSTRM_,NSTREM,IUNIT(43),IOUT,
     &                         MXSTRM)
          IF(IUNIT(43).GT.0 .AND. KPER.EQ.1)
     &        CALL GWF1HYD1OT(GZ,LENGZ,RX,LENRX,IG,LENIG,RX(LCHYDM),
     &                        NUMH,IHYDMUN,0.0,HYDNOH,NROW,NCOL,
     &                        ITMUNI,IOUT)
          IF(IUNIT(50).GT.0)
     &        CALL GWF1MNW1RP(MNWSITE,RZ(LCWEL2),NWELL2,MXWEL2,
     &                         GX(LCHOLD),RZ(LCHREF),IG(LCIBOU),
     &                         GX(LCDELR),GX(LCDELC),GX(LCCR),GX(LCCC),
     &                         RX(LCHY),GZ(LCHNEW),HCLOSE,SMALL,HDRY,
     &                         NODES,NROW,NCOL,KPER,KSPREF,IUNIT(50),
     &                         IOUT,IOWELL2,TOTIM,LAYHDT,GX(LCBOTM),
     &                         NBOTM,X(LCHK),IUNIT(1),IUNIT(23),
     &                         IUNIT(37),NLAY,RX(LCTRPY),
     &                         X(LCHKCC),X(LCHANI))
C-----READ MNWO LOCATION (MNW OBSERVATIONS FOR TRANSPORT)
cgzh debug (only read once, not each SP.  MNW Sites checked each SP)
          IF(JUNIT(18).GT.0.AND.KKPER.EQ.IPERGWT)
     &        CALL GWT1MNWO5RP(MNWOLST,IR(LSMNWU),MNWOBS,
     *                  IOUTS,JUNIT(18),
     *                  RZ(lcmnw2),WELLID,mxwel2,nwell2,
     &                  RX(LSMNWO), mnw2, NMNWVL)
CGWT----READ PARTICLE OBSERVATION LOCATIONS
      IF(JUNIT(23).GT.0.AND.KKPER.EQ.IPERGWT)
     *   CALL GWT1PTOB5RP(IX(LSPTOB),IX(LSPTMN),NUMPTOB,NUMPTOB_MNW,
     *                  MNWsite,mxwel2,nwell2,IX(LSPTUN),RZ(LCWEL2),
     *                  PTOBLST,ncol,nrow,nlay,nscol,nsrow,nslay,
     *                  IOUTS,JUNIT(23))
CGWT----DEFINE PTOB MNW LOCATIONS FOR THIS STRESS PERIOD
!      IF(JUNIT(23).GT.0)
!     *   CALL GWT1PTOB5RPS(IX(LSPTMN),NUMPTOB_MNW,
!     *              MNWsite,WELLID,mxwel2,nwell2,IX(LSPTUN),RZ(LCWEL2),
!     *              PTOBLST,ncol,nrow,nlay,nscol,nsrow,nslay,
!     *              IOUTS,JUNIT(23))
CMOCWT
CGWT----READ CONSTANT-CONCENTRATION BOUNDARY CONDITION PACKAGE
cgzh ccbd
      IF(IUNIT(15).GT.0.AND.JUNIT(25).GT.0)
     *   CALL GWT1CCBD1RPS(JUNIT(25),RX(LSCCBD),X(LSCONC),IOUTS,
     *     NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,KKPER,IPERGWT,
     *     X(LSDPCON),X(LSDPZO),JUNIT(10),
     *     X(LSDKZO),X(LSDKZS),JUNIT(11))
C
C--MNW2
          IF(IUNIT(57).GT.0)
     &        CALL GWF1MNW2RP(WELLID,RZ(LCMNW2),RZ(LCMNWN),RZ(LCMNWI),
     +                 GX(LCBOTM),NBOTM,IG(LCIBOU),
     &                 mnwmax,NODTOT,nrow,ncol,nlay,INTTOT,
     &                 iunit(57),iout,iunit(15),nmnw2,small,hclose,
     &                 KKPER,MNWPRNT,ntotnod,RZ(LCMNWC),NMNWVL)
CGWT----DEFINE PTOB MNW LOCATIONS FOR THIS STRESS PERIOD
      IF(JUNIT(23).GT.0)
     *   CALL GWT1PTOB5RPS(IX(LSPTMN),NUMPTOB_MNW,
     *              MNWsite,WELLID,mxwel2,nwell2,IX(LSPTUN),RZ(LCWEL2),
     *              PTOBLST,ncol,nrow,nlay,nscol,nsrow,nslay,
     *              IOUTS,JUNIT(23),RZ(LCMNW2), MNWMAX, IX(LSPTOB), 
     *              RZ(LCMNWN), NODTOT, NMNWVL)
C--MNW2
          IF(IUNIT(58).GT.0.AND.KKPER.EQ.1)
     &        CALL GWF1MNWIRP(MNWOBS,IOUT,RZ(LCMNIO),IUNIT(58),MNWIID,
     &                 MNWMAX,RZ(LCMNW2),WELLID,IUNIT(15),RZ(LSMNWO),
     &                 NMNWVL)   
C-----INITIALIZE SV ARRAY
          IF (ISEN.GT.0 .AND. IUNIT(23).GT.0)
     &        CALL SEN1LPF1SV(IG(LCIZON),KKPER,NCOL,NLAY,NMLTAR,NPLIST,
     &                        NROW,NZONAR,GX(LCRMLT),X(LCSV))
C
CGWT----INITIALIZE IBOUND INFO FOR GWT
          IF(IGWTON.GT.0)
     &        CALL GWT1BAS6CH(IG(LCIBOU),IFXHED,ICONLY,
     *                        NSCOL,NSROW,NSLAY,
     *                        NCOL,NROW,NLAY,IOUTS)
C
cgzh this was moved out of CKRP due to CCBD potentially changing the init conc array
CGWT----PRINT INITIAL CONCENTRATION
          KKSTP = 0
          IF(IGWTON.GT.0.AND.KKPER.EQ.IPERGWT) THEN
            CALL SMOC6C(X(LSCONC),X(LSTHCK),SBVL,SRCDCY,
     *        KKSTP,NSTP,KKPER,NPER,0,0,0.0,
     *        GX(LCDELR),GX(LCDELC),NCOL,NROW,NLAY,
     *        NSCOL,NSROW,NSLAY,IOUTS,JUNIT,PERTIM,TOTIM,0.0,
     *        NPNTCL,ICONFM,ICONLY,
     *        IDKZO,IDKFO,IDKZS,IDKFS,
     *        IDPZO,IDPFO,IUNIT(44),IUNIT(22),IUNIT(40),NIUNIT,
C MOCWT
     *        0.0,IUNIT,MULTSS,
cea
     *        SRCAGE,IUNIT(57))

CDP----PRINT INITIAL DOUBLE POROSITY CONCENTRATION DATA
            IF(JUNIT(10).GT.0.AND.IDPPS.NE.0) CALL SDP6C(X(LSDPCON),
     *        KKSTP,NSTP,KKPER,NPER,0,0,
     *        NSCOL,NSROW,NSLAY,IOUTS,JUNIT,PERTIM,TOTIM,0.0,
     *        NPNTCL,ICONFM,IDPPS,NIUNIT)
C
CGWT----WRITE INITIAL CONDITION TO OBSERVATION WELLS
cgzh debug zinn (had commented this out to get rid of HNEW print at init cond
cgzh send in SUMTCH=0.0
            SUMTCH=0.0
            IF (JUNIT(8).GT.0)
     *        CALL SOBS5O(GZ(LCHNEW),X(LSCONC),SUMTCH,IX(LSOBSW),0,
     *          NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,NUMOBS,IOBSFL,
     *          GZ(LCHNEW),GZ(LCHNEW))
C
          END IF
C         END INIT GWT PRINTS
C
C7C-----SIMULATE EACH TIME STEP.
          DO 90 KSTP = 1, NSTP(KKPER)
            KKSTP = KSTP
C
CGWT----COPY BEGINNING OF FLOW STEP IBOUND FOR ELLAM
c           IF(MOCTYPE.EQ.3.AND.ISSFLG(KKPER).EQ.0.AND.IGWTON.GT.0)
c     *          CALL GWT1IBOU6(IG(LCIBOU),IX(LSIBOU),NODES)
C
CGWT SAVE OLD IBOUND FOR REWETTING 
CMOCWT  MOCWT only
cgzh rewet
           IF(IGWTON.GT.0.AND.PTWTON.GT.0) THEN
             IF(IWDFLG.GT.0) THEN
c do this for all time steps after first one [same as: .NOT.(first ts)]
			 IF(.NOT.(KKPER.EQ.IPERGWT.AND.KKSTP.EQ.1)) THEN
                 CALL GWT1IBOU6(IG(LCIBOU),IX(LSIBOU),NODES)
               END IF
             END IF
           END IF
C7C1----CALCULATE TIME STEP LENGTH. SET HOLD=HNEW.
            CALL GWF1BAS6AD(DELT,TSMULT(KKPER),TOTIM,PERTIM,GZ(LCHNEW),
     1                      GX(LCHOLD),KKSTP,NCOL,NROW,NLAY,ITS)
            IF (IUNIT(20).GT.0)
     &          CALL GWF1CHD6AD(NCHDS,MXCHD,RX(LCCHDS),GZ(LCHNEW),
     &                          GX(LCHOLD),PERLEN(KKPER),PERTIM,NCOL,
     &                          NROW,NLAY,NCHDVL,IOUT)
            IF (IUNIT(1).GT.0)
     &          CALL GWF1BCF6AD(IG(LCIBOU),GX(LCHOLD),GX(LCBOTM),NBOTM,
     &                          RX(LCWETD),IWDFLG,ISSFLG(KKPER),NCOL,
     &                          NROW,NLAY)
            IF (IUNIT(23).GT.0)
     &          CALL GWF1LPF1AD(IG(LCIBOU),GX(LCHOLD),GX(LCBOTM),
     &                          X(LCWETD),ISSFLG(KKPER),NCOL,NROW,NLAY,
     &                          NBOTM)
            IF (IUNIT(37).GT.0)
     &          CALL GWF1HUF2AD(IG(LCIBOU),GX(LCHOLD),GX(LCBOTM),
     &                          X(LCWETD),ISSFLG(KKPER),NCOL,NROW,NLAY,
     &                          NBOTM)
            IF(IUNIT(16).GT.0)
     &          CALL GWF1FHB1AD(GZ(LCHNEW),GX(LCHOLD),NCOL,NROW,NLAY,
     &                      ITRSS,TOTIM,DELT,RX(LCBDTM),NBDTIM,
     &                      RX(LCFLRT),RX(LCBDFV),RX(LCBDHV),NFLW,
     &                      RX(LCSBHD),IR(LCHDLC),NHED,NFHBX1,NFHBX2,
     &                      IFHBD3,IFHBD4,IFHBD5,IFHBSS,NHEDDIM,NFLWDIM,
     &                      NBDHVDIM)
            IF (IUNIT(17).GT.0)
     &          CALL GWF1RES1AD(RX(LCHRES),RX(LCHRSE),IR(LCIRES),
     &                      RX(LCBRES),GX(LCDELR),GX(LCDELC),NRES,
     &                      IRESPT,NCOL,NROW,PERLEN(KKPER),PERTIM,TOTIM,
     &                      KKSTP,KKPER,IOUT)
CLAK
            IF (IUNIT(22).GT.0)
     1          CALL GWF1LAK3AD(KKPER,KKSTP,NLAKES,RX(ISTGLD),
     2                      RX(ISTGNW),RX(LCSTAG),NROW,NCOL,NLAY,
     3                      RX(LSFLOB),IUNIT(15),LKNODE,RX(ISTGLD2))
            IF(IUNIT(51).GT.0)
     1          CALL GWF1DAF1AD(DELT,IERR,ITMUNI,IUNIT(51),IOUT)
            IF(IUNIT(50).GT.0)
     &          CALL GWF1MNW1AD(NWELL2,MXWEL2,RZ(LCWEL2),IG(LCIBOU),
     &                          GX(LCDELR),GX(LCDELC),GX(LCCR),GX(LCCC),
     &                          RX(LCHY),SMALL,HDRY,GZ(LCHNEW),NCOL,
     &                          NROW,NODES,LAYHDT,GX(LCBOTM),NBOTM,
     &                          X(LCHK),IUNIT(1),IUNIT(23),IUNIT(37),
     &                          NLAY,RX(LCTRPY),X(LCHKCC),
     &                          X(LCHANI))
C--MNW2
            IF(IUNIT(57).GT.0) THEN
cgzh mnw2 calculate transmissivity terms needed for MNW conductance calculation
             IF (IUNIT(1).GT.0) THEN
              CALL GWF1MNW2BCF(GX(LCDELR),GX(LCDELC),GX(LCCR),GX(LCCC),
     +                 RX(LCHY),GZ(LCHNEW),ncol,nrow,nlay,Hdry,small,
     &                 LAYHDT,GX(LCBOTM),NBOTM,RX(LCTRPY),
     &                 RZ(LCMNWN),NODTOT,nmnw2,MNWMAX,RZ(LCMNW2),
     &                 WELLID,MNWPRNT,iout,NMNWVL)
             ELSE IF (IUNIT(23).GT.0) THEN
              CALL GWF1MNW2LPF(GX(LCDELR),GX(LCDELC),GZ(LCHNEW),
     +                 ncol,nrow,nlay,Hdry,
     &                 LAYHDT,GX(LCBOTM),NBOTM,X(LCHK),X(LCHANI),
     &                 RZ(LCMNWN),NODTOT,nmnw2,MNWMAX,RZ(LCMNW2),
     &                 WELLID,MNWPRNT,iout,NMNWVL)
             ELSE IF(IUNIT(37).GT.0) THEN
              CALL GWF1MNW2HUF(GX(LCDELR),GX(LCDELC),GZ(LCHNEW),
     +                 ncol,nrow,nlay,Hdry,
     &                 LAYHDT,GX(LCBOTM),NBOTM,X(LCHK),X(LCHKCC),
     &                 RZ(LCMNWN),NODTOT,nmnw2,MNWMAX,RZ(LCMNW2),
     &                 WELLID,MNWPRNT,iout,NMNWVL)
	       ELSE
              write(iout,1000)
 1000 FORMAT(/1X,
     &'***ERROR: MNW2 PACKAGE DOES NOT SUPPORT',/,
     &' SELECTED FLOW PACKAGE',/,
     &' (MNW2 DOES FULLY SUPPORT BCF, LPF, AND HUF PACKAGES)',/,
     &' -- STOP EXECUTION')
              CALL USTOP('MNW2 error-flow package')
             END IF
cgzh mnw2  if BCF (SC1 in RX)
             IF(IUNIT(1).GT.0) THEN
               CALL GWF1MNW2AD(nmnw2,MNWMAX,RZ(LCMNW2),NODTOT,
     +                      RZ(LCMNWN),IG(LCIBOU),RX(LCSC1),
     +                      GX(LCDELR),GX(LCDELC),GX(LCBOTM),NBOTM,
     +                      GZ(LCHNEW),GX(LCHOLD),GX(LCSTRT),
     +                      ncol,nrow,nlay,small,kkstp,
     +                      GX(LCCV),RZ(LCMNWI),INTTOT,IOUT,WELLID,
     +                      IUNIT(1),ISSFLG(KKPER),LAYHDT,
     +                      X(LCHK),X(LCVKA),KKPER,MNWPRNT,NMNWVL)
C           if LPF or HUF (SC1 in X)
             ELSE
               CALL GWF1MNW2AD(nmnw2,MNWMAX,RZ(LCMNW2),NODTOT,
     +                      RZ(LCMNWN),IG(LCIBOU),X(LCSC1),
     +                      GX(LCDELR),GX(LCDELC),GX(LCBOTM),NBOTM,
     +                      GZ(LCHNEW),GX(LCHOLD),GX(LCSTRT),
     +                      ncol,nrow,nlay,small,kkstp,
     +                      GX(LCCV),RZ(LCMNWI),INTTOT,IOUT,WELLID,
     +                      IUNIT(1),ISSFLG(KKPER),LAYHDT,
     +                      X(LCHK),X(LCVKA),KKPER,MNWPRNT,NMNWVL)
             END IF
            END IF
C
C---------INDICATE IN PRINTOUT THAT SOLUTION IS FOR HEADS
            CALL UMESPR('SOLVING FOR HEAD',' ',IOUT)
C-----------SHOW PROGRESS IF REQUESTED
            IF(SHOWPROG)THEN
              IF (ITERPK.GT.1 .AND. NPER.EQ.1) THEN
                WRITE (*,'(A)') ' '
              ENDIF
              WRITE(*,25)KPER,KSTP
            ENDIF
   25 FORMAT('+Solving:  Stress period: ',i5,4x,
     &       'Time step: ',i5,4x,'Ground-Water Flow Eqn.')
C
C7C2----ITERATIVELY FORMULATE AND SOLVE THE FLOW EQUATIONS.
C
            DO 30 KITER = 1, MXITER
              KKITER = KITER
C
C7C2A---FORMULATE THE FINITE DIFFERENCE EQUATIONS.
              CALL GWF1BAS6FM(GX(LCHCOF),GX(LCRHS),NODES)
              IF (IUNIT(1).GT.0)
     &            CALL GWF1BCF6FM(GX(LCHCOF),GX(LCRHS),GX(LCHOLD),
     &                            RX(LCSC1),GZ(LCHNEW),IG(LCIBOU),
     &                            GX(LCCR),GX(LCCC),GX(LCCV),RX(LCHY),
     &                            RX(LCTRPY),GX(LCBOTM),NBOTM,RX(LCSC2),
     &                            GX(LCDELR),GX(LCDELC),DELT,
     &                            ISSFLG(KKPER),KKITER,KKSTP,KKPER,NCOL,
     &                            NROW,NLAY,IOUT,RX(LCWETD),IWDFLG,
     &                            RX(LCCVWD),WETFCT,IWETIT,IHDWET,HDRY,
     &                            GX(LCBUFF))
              IF (IUNIT(23).GT.0)
     &            CALL GWF1LPF1FM(GX(LCHCOF),GX(LCRHS),GX(LCHOLD),
     &                            X(LCSC1),GZ(LCHNEW),IG(LCIBOU),
     &                            GX(LCCR),GX(LCCC),GX(LCCV),X(LCHK),
     &                            X(LCHANI),X(LCVKA),GX(LCBOTM),
     &                            X(LCSC2),GX(LCDELR),GX(LCDELC),DELT,
     &                            ISSFLG(KKPER),KKITER,KKSTP,KKPER,NCOL,
     &                            NROW,NLAY,IOUT,X(LCWETD),WETFCT,
     &                            IWETIT,IHDWET,HDRY,NBOTM,X(LCVKCB))
              IF (IUNIT(37).GT.0)
     &            CALL GWF1HUF2FM(GX(LCHCOF),GX(LCRHS),GX(LCHOLD),
     &                            X(LCSC1),GZ(LCHNEW),IG(LCIBOU),
     &                            GX(LCCR),GX(LCCC),GX(LCCV),X(LCHK),
     &                            X(LCVKA),GX(LCBOTM),GX(LCDELR),
     &                            GX(LCDELC),DELT,ITRSS,ISSFLG(KKPER),
     &                            NCOL,NROW,NLAY,IOUT,X(LCWETD),NBOTM,
     &                            NHUF,GX(LCRMLT),IG(LCIZON),NMLTAR,
     &                            NZONAR,X(LCHUFTHK),X(LCHKCC),HDRY,
     &                            KKITER,KSTP,KPER,X(LCHUFTMP),
     &                            IX(LCHGUF),
     &                            IUNIT(47),X(LCVDHD),X(LCVDHT),IWETIT,
     &                            IHDWET,WETFCT,X(LCGS),X(LCA9),HNOFLO)
              IF (IUNIT(21).GT.0)
     &            CALL GWF1HFB6FM(GX(LCBOTM),GX(LCCC),GX(LCCR),
     &                            GX(LCDELC),GX(LCDELR),RX(LCHFB),
     &                            GZ(LCHNEW),MXACTFB,NBOTM,NCOL,NHFB,
     &                            NLAY,NROW,LAYHDT)
              IF (IUNIT(2).GT.0)
     &            CALL GWF1WEL6FM(NWELLS,MXWELL,GX(LCRHS),RX(LCWELL),
     &                            IG(LCIBOU),NCOL,NROW,NLAY,NWELVL)
              IF (IUNIT(3).GT.0)
     &            CALL GWF1DRN6FM(NDRAIN,MXDRN,RX(LCDRAI),GZ(LCHNEW),
     &                            GX(LCHCOF),GX(LCRHS),IG(LCIBOU),NCOL,
     &                            NROW,NLAY,NDRNVL)
              IF (IUNIT(4).GT.0)
     &            CALL GWF1RIV6FM(NRIVER,MXRIVR,RX(LCRIVR),GZ(LCHNEW),
     &                            GX(LCHCOF),GX(LCRHS),IG(LCIBOU),NCOL,
     &                            NROW,NLAY,NRIVVL)
CLAK
              IF (IUNIT(5).GT.0) THEN
                IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3)
     1                CALL GWF1LAK3ST(0,NCOL,NROW,NLAY,IG(LCIBOU),
     2                            LKNODE,IR(ICLAKE),MXLKND,NLAKES,
     3                            RX(ISTGLD),GX(LCBOTM),NBOTM)
                CALL GWF1EVT6FM(NEVTOP,IR(LCIEVT),RX(LCEVTR),
     &                          RX(LCEXDP),RX(LCSURF),GX(LCRHS),
     &                          GX(LCHCOF),IG(LCIBOU),GZ(LCHNEW),NCOL,
     &                          NROW,NLAY)
                IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3)
     1                CALL GWF1LAK3ST(1,NCOL,NROW,NLAY,IG(LCIBOU),
     2                            LKNODE,IR(ICLAKE),MXLKND,NLAKES,
     3                            RX(ISTGLD),GX(LCBOTM),NBOTM)
              END IF
              IF (IUNIT(7).GT.0)
     &            CALL GWF1GHB6FM(NBOUND,MXBND,RX(LCBNDS),GX(LCHCOF),
     &                            GX(LCRHS),IG(LCIBOU),NCOL,NROW,NLAY,
     &                            NGHBVL)
CLAK
              IF (IUNIT(8).GT.0) THEN
                IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3)
     1                CALL GWF1LAK3ST(0,NCOL,NROW,NLAY,IG(LCIBOU),
     2                            LKNODE,IR(ICLAKE),MXLKND,NLAKES,
     3                            RX(ISTGLD),GX(LCBOTM),NBOTM)
                CALL GWF1RCH6FM(NRCHOP,IR(LCIRCH),RX(LCRECH),
     &                            GX(LCRHS),IG(LCIBOU),NCOL,NROW,NLAY)
                IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3)
     1                CALL GWF1LAK3ST(1,NCOL,NROW,NLAY,IG(LCIBOU),
     2                            LKNODE,IR(ICLAKE),MXLKND,NLAKES,
     3                            RX(ISTGLD),GX(LCBOTM),NBOTM)
              END IF
              IF(IUNIT(16).GT.0)
     &            CALL GWF1FHB1FM(GX(LCRHS),IG(LCIBOU),IR(LCFLLC),
     &                        RX(LCBDFV),NFLW,NCOL,NROW,NLAY,IFHBD4,
     &                        NFLWDIM)
              IF (IUNIT(17).GT.0)
     &            CALL GWF1RES1FM(IR(LCIRES),IR(LCIRSL),RX(LCBRES),
     &                        RX(LCCRES),RX(LCBBRE),RX(LCHRES),
     &                        IG(LCIBOU),GZ(LCHNEW),GX(LCHCOF),
     &                        GX(LCRHS),NRES,NRESOP,NCOL,NROW,NLAY)
              IF(IUNIT(18).GT.0)
     &            CALL GWF1STR6FM(NSTREM,RX(LCSTRM_),IR(ICSTRM_),
     &                        GZ(LCHNEW),GX(LCHCOF),GX(LCRHS),
     &                        IG(LCIBOU),MXSTRM,NCOL,NROW,NLAY,
     &                        NSSSTR6,IR(LCTBAR),NTRIB,RX(LCTRIB),
     &                        IR(LCIVAR_),IR(LCFGAR),ICALC,
     &                        CONSTSTR6)
              IF (IUNIT(19).GT.0)
     1            CALL GWF1IBS6FM(GX(LCRHS),GX(LCHCOF),GZ(LCHNEW),
     2                            GX(LCHOLD),RX(LCHC),RX(LCSCE),
     3                            RX(LCSCV),IG(LCIBOU),NCOL,NROW,NLAY,
     4                            DELT,ISSFLG(KKPER),IBSDIM)
              IF(IUNIT(54).GT.0)
     1            CALL GWF1SUB1FM(GX(LCRHS),GX(LCHCOF),GZ(LCHNEW),
     2                            GX(LCHOLD),IG(LCIBOU),X(LCV),
     3                            GX(LCDELR),GX(LCDELC),NCOL,NROW,
     4                            NODES,DELT,AC1,AC2,HCLOSE,KKITER,
     5                            ITMIN,NN,NND1,ND1,ND2,NDB,NNDB,NPZ,
     6                            ISSFLG(KKPER),IUNIT(9))
              IF(IUNIT(55).GT.0)
     1            CALL GWF1SWT1FM(GX(LCRHS),GX(LCHCOF),IG(LCIBOU),
     2                            GZ(LCHNEW),GX(LCHOLD),GX(LCBOTM),
     3                            NCOL,NROW,NLAY,DELT,ISSFLG(KKPER),
     4                            kkper,iout)
Cdep  Changed SFR1 to SFR2
              IF (IUNIT(44).GT.0)
     1               CALL GWF1SFR2FM(RX(LCSTRM),IR(ICSTRM),
     2        GZ(LCHNEW),GX(LCHOLD),GX(LCHCOF),GX(LCRHS),IG(LCIBOU),
     3        NSTRM,NCOL,NROW,NLAY,IOUT,NSS,NSEGDIM,RX(LCSEG),IR(ICSEG),
     4        IR(LCOTSG),RX(LCXSEC),IR(LCIVAR),RX(LCQSTG),CONST,MAXPTS,
     5        DLEAK,RX(LCOTFLW),RX(LCDVFLW),NLAKESAR,RX(ISTGLD),
     6        RX(ISTRIN),RX(ISTROT),RX(ISTGNW),THETA,RX(LKACC7),
     7        ISSFLG(KKPER),RX(IDSTRT),RX(LCSFRQ),IUNIT(22),KKITER,
     8        RZ(LCSUZDPIT),RZ(LCSUZTHIT),RZ(LCSUZSPIT),RZ(LCSUZFLIT),
     9        RZ(LCSUZDPST),RZ(LCSUZTHST),RZ(LCSUZSPST),RZ(LCSUZFLST),
     &        RZ(LCSEPS),RZ(LCSTHTS),RZ(LCSTHTR),IR(ICSLTRLIT),
     &        IR(ICSITRLIT),IR(ICSLTRLST),IR(ICSITRLST),RZ(LCSUZFLWT),
     &        IR(ICSNWAVST),NSTRAIL,NSTOTRL,IUZT,ISUZN,NUZST,DELT,
     &        RZ(LCSUZWDTH),RZ(LCSOLSFLX),IR(ICSLOOP),RX(LCUHC),
     &        TOTIM,IR(ICELEV),RX(LCDPTH),RZ(LCWETP),NSFRSETS,
     &        RZ(LCSUZSEEP),RZ(LCOLDFLBT),KKPER,KKSTP,NUMCELL,
     &        NUZROW,NUZCOL)
Cdep  End of change
Cdep  Added 
CLAK
Cdep Revised Lake Package call statement  June 4, 2006
Cdep Added SURFDEPTH to end of call statement  March 23, 2009
              IF (IUNIT(22).GT.0)
     *               CALL GWF1LAK3FM(LKNODE,MXLKND,IR(ICLAKE),
     1                      GZ(LCHNEW),GX(LCHCOF),GX(LCRHS),
     2                      IG(LCIBOU),NCOL,NROW,NLAY,NLAKES,
     3                      RX(ISTGLD),RX(LCCNDF),GX(LCBOTM),NBOTM,
     4                      IOUT,DELT,NSS,NTRB,NDV,IR(INTRB),IR(INDV),
     5                      RX(ISTRIN),RX(ISTROT),RX(ISTGNW),
     6                      RX(LCWTDR),RX(LCLKPR),RX(LCLKEV),
     7                      GX(LCDELR),GX(LCDELC),RZ(LKACC1),
     8                      RZ(LKACC2),RZ(LKACC3),RX(LKACC4),
     9                      RX(LKACC5),RX(LKACC6),THETA,RX(LCRNF),
     *                      KKSTP,KKITER,ISSFLG(KKPER),NSSITR,SSCNCR,
     *                      RX(LKSSMN),RX(LKSSMX),RX(IDSTRT),IR(LKNCN),
     *                      RX(LKDSR),RX(LKCNN),RX(LKCHN),RX(IAREN),
     *                      IR(LKLMRR),NSSAR,IUNIT(44),RX(ISTGITR),
     *                      RZ(LKSEP3),RX(IAREATAB),RX(IDPTHTAB),
     *                      NSSLK,RZ(ISLKOTFLW),RZ(IDLKOTFLW),
     *                      RZ(IDLKSTAGE),RX(IBTMS),RX(LCEVAPO),
     *                      RX(LCFLWIN),RX(LCFLWIT),RX(LCWITDW),
     *                      RX(LCGWRAT),SURFDEPTH)
              IF (IUNIT(39).GT.0)
     &            CALL GWF1ETS1FM(NETSOP,IR(LCIETS),RX(LCETSR),
     &                            RX(LCETSX),RX(LCETSS),GX(LCRHS),
     &                            GX(LCHCOF),IG(LCIBOU),GZ(LCHNEW),NCOL,
     &                            NROW,NLAY,NETSEG,RX(LCPXDP),
     &                            RX(LCPETM),NSEGAR)
              IF (IUNIT(40).GT.0)
     &            CALL GWF1DRT1FM(NDRTCL,MXDRT,RX(LCDRTF),GZ(LCHNEW),
     &                            GX(LCHCOF),GX(LCRHS),IG(LCIBOU),NCOL,
     &                            NROW,NLAY,NDRTVL,IDRTFL)
              IF(IUNIT(51).GT.0)
     1            CALL GWF1DAF1FM(IERR,ITMUNI,GZ(LCHNEW),GX(LCHOLD),
     2                            IOUT,IG(LCIBOU),GX(LCHCOF),
     3                            GX(LCRHS),NCOL,NROW,NLAY,KITER,
     4                            IDAFBK)
              IF (IUNIT(50).GT.0)
     &            CALL GWF1MNW1FM(NWELL2,MXWEL2,RZ(LCWEL2),IG(LCIBOU),
     &                            GX(LCDELR),GX(LCDELC),GX(LCCR),
     &                            GX(LCCC),RX(LCHY),SMALL,HDRY,
     &                            GX(LCHCOF),GX(LCRHS),GZ(LCHNEW),NCOL,
     &                            NROW,NODES,KITER,NOMOITER,LAYHDT,
     &                            GX(LCBOTM),NBOTM,X(LCHK),IUNIT(1),
     &                            IUNIT(23),IUNIT(37),NLAY,
     &                       RX(LCTRPY),X(LCHKCC),X(LCHANI))
C--MNW2
            IF(IUNIT(57).GT.0) THEN
cgzh mnw2 calculate transmissivity terms needed for MNW conductance calculation
             IF (IUNIT(1).GT.0) THEN
              CALL GWF1MNW2BCF(GX(LCDELR),GX(LCDELC),GX(LCCR),GX(LCCC),
     +                 RX(LCHY),GZ(LCHNEW),ncol,nrow,nlay,Hdry,small,
     &                 LAYHDT,GX(LCBOTM),NBOTM,RX(LCTRPY),
     &                 RZ(LCMNWN),NODTOT,nmnw2,MNWMAX,RZ(LCMNW2),
     &                 WELLID,MNWPRNT,iout,NMNWVL)
             ELSE IF (IUNIT(23).GT.0) THEN
              CALL GWF1MNW2LPF(GX(LCDELR),GX(LCDELC),GZ(LCHNEW),
     +                 ncol,nrow,nlay,Hdry,
     &                 LAYHDT,GX(LCBOTM),NBOTM,X(LCHK),X(LCHANI),
     &                 RZ(LCMNWN),NODTOT,nmnw2,MNWMAX,RZ(LCMNW2),
     &                 WELLID,MNWPRNT,iout,NMNWVL)
             ELSE IF(IUNIT(37).GT.0) THEN
              CALL GWF1MNW2HUF(GX(LCDELR),GX(LCDELC),GZ(LCHNEW),
     +                 ncol,nrow,nlay,Hdry,
     &                 LAYHDT,GX(LCBOTM),NBOTM,X(LCHK),X(LCHKCC),
     &                 RZ(LCMNWN),NODTOT,nmnw2,MNWMAX,RZ(LCMNW2),
     &                 WELLID,MNWPRNT,iout,NMNWVL)
             END IF
cgzh mnw2  if bcf
             IF(IUNIT(1).GT.0) THEN
               CALL GWF1MNW2FM(nmnw2,ncol,nrow,nlay,GX(LCDELR),
     &             GX(LCDELC),
     &             RX(LCSC1),GX(LCBOTM),NBOTM,IG(LCIBOU),RZ(LCMNWN),
     &             NODTOT,MNWMAX,RZ(LCMNW2),KITER,GZ(LCHNEW),GX(LCHOLD),
     &             GX(LCSTRT),GX(LCHCOF),GX(LCRHS),small,GX(LCCV),
     &             RZ(LCMNWI),INTTOT,IOUT,kkstp,IUNIT(1),ISSFLG(KKPER),
     &             LAYHDT,wellid,X(LCHK),X(LCVKA),KKPER,HDRY,RZ(LCMNWC),
     &             MNWPRNT,hclose,NMNWVL)
             ELSE
cgzh mnw2  if lpf or huf
               CALL GWF1MNW2FM(nmnw2,ncol,nrow,nlay,GX(LCDELR),
     &             GX(LCDELC),
     &             X(LCSC1),GX(LCBOTM),NBOTM,IG(LCIBOU),RZ(LCMNWN),
     &             NODTOT,MNWMAX,RZ(LCMNW2),KITER,GZ(LCHNEW),GX(LCHOLD),
     &             GX(LCSTRT),GX(LCHCOF),GX(LCRHS),small,GX(LCCV),
     &             RZ(LCMNWI),INTTOT,IOUT,kkstp,IUNIT(1),ISSFLG(KKPER),
     &             LAYHDT,wellid,X(LCHK),X(LCVKA),KKPER,HDRY,RZ(LCMNWC),
     &             MNWPRNT,hclose,NMNWVL)
             END IF
            END IF
C
C
C-------IF HNEW=HOLD=0 AND RHS=0, NO NEED TO SOLVE.
              CALL UNOITER(GX(LCRHS),GZ(LCHNEW),NODES,ISA)
              IF (ISA.EQ.0) THEN
                ICNVG = 1
                GOTO 33
              ENDIF
C
C7C2B---MAKE ONE CUT AT AN APPROXIMATE SOLUTION.
              IF (IUNIT(9).GT.0)
     &            CALL SIP5AP(GZ(LCHNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                        GX(LCCV),GX(LCHCOF),GX(LCRHS),X(LCEL),
     &                        X(LCFL),X(LCGL),X(LCV),X(LCW),X(LCHDCG),
     &                        IX(LCLRCH),NPARM,KKITER,HCLOSE,ACCL,
     &                        ICNVG,KKSTP,KKPER,IPCALC,IPRSIP,MXITER,
     &                        NSTP(KKPER),NCOL,NROW,NLAY,NODES,IOUT,0,
     &                        IERR,IERRU)
              IF (IUNIT(10).GT.0)
     &            CALL DE45AP(GZ(LCHNEW),IG(LCIBOU),X(LCAU),X(LCAL),
     &                        IX(LCIUPP),IX(LCIEQP),X(LCD4B),MXUP,
     &                        MXLOW,MXEQ,MXBW,GX(LCCR),GX(LCCC),
     &                        GX(LCCV),GX(LCHCOF),GX(LCRHS),ACCL,
     &                        KKITER,ITMX,MXITER,NITER,HCLOSE,IPRD4,
     &                        ICNVG,NCOL,NROW,NLAY,IOUT,IX(LCLRCH),
     &                        X(LCHDCG),IFREQ,KKSTP,KKPER,DELT,
     &                        NSTP(KKPER),ID4DIR,ID4DIM,MUTD4,IERR,
     &                        IERRU)
              IF (IUNIT(11).GT.0)
     &            CALL SOR5AP(GZ(LCHNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                        GX(LCCV),GX(LCHCOF),GX(LCRHS),X(LCA),
     &                        X(LCRES),IX(LCIEQP),X(LCHDCG),IX(LCLRCH),
     &                        KKITER,HCLOSE,ACCL,ICNVG,KKSTP,KKPER,
     &                        IPRSOR,MXITER,NSTP(KKPER),NCOL,NROW,NLAY,
     &                        NSLICE,MBW,IOUT,0)
              IF (IUNIT(13).GT.0)
     &            CALL PCG2AP(GZ(LCHNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                        GX(LCCV),GX(LCHCOF),GX(LCRHS),Z(LCV),
     &                        Z(LCSS),Z(LCP),X(LCCD),X(LCHCHG),
     &                        IX(LCLHCH),X(LCRCHG),IX(LCLRCH),KKITER,
     &                        NITER,HCLOSE,RCLOSE,ICNVG,KKSTP,KKPER,
     &                        IPRPCG,MXITER,ITER1,NPCOND,NBPOL,
     &                        NSTP(KKPER),NCOL,NROW,NLAY,NODES,RELAX,
     &                        IOUT,MUTPCG,IX(LCIT1),DAMP,GX(LCBUFF),
     &                        X(LCHCSV),IERR,IERRU,Z(LCHPCG))
              IF (IUNIT(14).GT.0)
     &            CALL LMG1AP(GZ(LCHNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                        GX(LCCV),GX(LCHCOF),GX(LCRHS),Z(LCA),
     &                        IX(LCIA),IX(LCJA),Z(LCU1),Z(LCFRHS),
     &                        IX(LCIG),ISIZ1,ISIZ2,ISIZ3,ISIZ4,KKITER,
     &                        BCLOSE,DAMP,ICNVG,KKSTP,KKPER,MXITER,
     &                        MXCYC,NCOL,NROW,NLAY,NODES,HNOFLO,IOUT,
     &                        IOUTAMG,ICG,IADAMP,DUP,DLOW)
c              IF (IUNIT(42).GT.0)
c     &            CALL GMG1AP(GZ(LCHNEW),GX(LCRHS),GX(LCCR),GX(LCCC),
c     &                        GX(LCCV),GX(LCHCOF),HNOFLO,IG(LCIBOU),
c     &                        IITER,MXITER,RCLOSE,HCLOSE,KKITER,KKSTP,
c     &                        KKPER,ICNVG,DAMP,IADAMP,IOUTGMG,IOUT)
              IF (IERR.GT.0) THEN
C               WRITE MESSAGE RELATED TO BEALE'S MEASURE, IF
C               APPROPRIATE, THEN STOP EXECUTION.
                IF (IBEFLG.EQ.2) CALL PES1BAS6ER(IOUT,ITERPK,NPLIST)
                CALL PLL1SD(IERR,IERRU,IOUT,IOUTG)
              ENDIF
C
C7C2C---IF CONVERGENCE CRITERION HAS BEEN MET STOP ITERATING.
              IF (ICNVG.EQ.1) GOTO 33
  30        CONTINUE
            KITER = MXITER
C
C7C2C-----IF CONVERGENCE CRITERION HAS NOT BEEN MET . . .
C---------IF ESTIMATING PARAMETERS OR CALCULATING BEALE'S MEASURE, USE
C         THE AVAILABLE VALUES AND KEEP GOING
            IF (IPES.GT.0 .OR. IBEFLG.EQ.2) THEN
              ICNVG = 1
              ICNVGP = 0
            ENDIF
C---------PRINT THE DATA TABLE AND WARNING MESSAGES AND STOP EXCEPT
C         AS NOTED ABOVE
            IF (IPAR.GT.-3) THEN
              CALL UNOCONV(X(LCBUF1+IPRAR),OBSNAM,X(LCH),
     &                     X(LCHOBS),IOUT,0,IPAR,IPR,KKPER,KKSTP,
     &                     IX(LCLN),MPR,ND,NDMH,NH,IX(LCNIPR),X(LCPRM),
     &                     X(LCBUF1+IPRAR+ND+MPR+IPR),RSQ,
     &                     RSQP,X(LCWP),X(LCWTPS),X(LCWT),X(LCWTQ),
     &                     X(LCWTQS),NPLIST,MPRAR,IPRAR,OUTNAM,
     &                     IX(LCIPLO),EQNAM,NAMES,IX(LCIPLP),NDMHAR,
     &                     NQTDR,NQTRV,NQTGB,NQTST,NQTCH,IOWTQCH,
     &                     IOWTQDR,IOWTQRV,IOWTQGB,IOWTQST,LCOBBAS,
     &                     LCOBDRN,LCOBRIV,LCOBGHB,LCOBSTR,LCOBCHD,
     &                     LCOBADV,X(LCSSGF),X(LCSSDR),X(LCSSRV),
     &                     X(LCSSGB),X(LCSSST),X(LCSSAD),X(LCSSCH),
     &                     X(LCSSPI),X(LCSSTO),ITMXP,IPES,X(LCBPRI),
     &                     ITERP,IERR,IERRU,NTT2,LCOBDRT,X(LCSSDT),
     &                     NQTDT,IOWTQDT,NRSO,NPOST,NNEGT,NRUNS,NQTSF,
     &                     IOWTQSF,LCOBSFR,X(LCSSSF),KTDIM,NHT,
     &                     X(LCOTIM))
              IF (IPAR.EQ.1 .OR. IBEFLG.EQ.2) THEN
C               CONTINUE EXECUTION, BUT WRITE MESSAGE(S) REGARDING
C               NONCONVERGENCE
                CALL PES1BAS6NC(GX(LCCC),GX(LCCR),GX(LCCV),GX(LCHCOF),
     &                          GZ(LCHNEW),IBEFLG,IG(LCIBOU),IOUTG,
     &                          ITERPK,KKPER,KKSTP,NCOL,NLAY,NROW,
     &                          GX(LCRHS))
              ELSE
C               STOP EXECUTION
                CALL PLL1SD(IERR,IERRU,IOUT,IOUTG)
              ENDIF
            ENDIF
C
   33       CONTINUE
CGWT    SET WTFAC FOR ADJUSTING PARTICLES IN WATER TABLE CELLS 
CGWT    AND ADJUST PARTICLE POSITIONS
C  update only if transient or 1st ts of ss stress period
c  same as: (not (SS and KSTP>1))
           IF(IGWTON.GT.0) THEN
            IF (.NOT.(ISSFLG(KKPER).EQ.1.AND.KKSTP.GT.1)) THEN
              CALL GWT1BAS6WTF(IG(LCIBOU),
     *               NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *               GZ(LCHNEW),GX(LCHOLD),GX(LCBOTM),NBOTM,
     *               NP,NPMAX,X(LSPC),X(LSPR),X(LSPL),
     *               IOUTS,X(LSWTFC),ISSFLG(KKPER),KKPER,
     *               IWDFLG,IPERGWT,KKSTP)
            ELSE
C  NO ADJUSTMENT UNLESS WT CHANGED
              CALL SMOC5Z(X(LSWTFC),NSCOL,NSROW,NSLAY,0.0)
		  ENDIF
           END IF
C
CGWT----INITIALIZE TRANSPORT TERMS BASED ON FLOW SOLUTION
C orig    IF(IUNIT(15).GT.0.AND.KKPER.EQ.1.AND.KKSTP.EQ.1) THEN
c orig code would not call again for subsequent SS periods 
c set thck at beginning of simulation 
c   and at beginning of all subsequent SS stress periods 
cgzh SSTR      IF(IUNIT(15).GT.0.AND.KKSTP.EQ.1) THEN
      IF(IGWTON.GT.0.AND.KKSTP.EQ.1) THEN
CMOCWT
cgzh debug update pt weights here (before thck update) for TR period
cgzh first calculate celvol
CMOCWT
cgzh SSTR        IF(KKPER.EQ.1.AND.
        IF(KKPER.EQ.IPERGWT.AND.
     *    ISSFLG(KKPER).EQ.0.AND.PTWTON.EQ.1) THEN
C   CALCULATE PORE VOLUME OF CELL (THCK*AREA*POROSITY)
            CALL GWT1BAS6CV(Z(LSCELV),X(LSTHCK),X(LSPOR),IG(LCIBOU),
     *        NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY)
        END IF
C
cgzh SSTR        IF((KKPER.EQ.1).OR.(KKPER.GT.1.AND.ISSFLG(KKPER).EQ.1)) THEN
        IF((KKPER.EQ.IPERGWT).OR.
     *    (KKPER.GT.IPERGWT.AND.ISSFLG(KKPER).EQ.1)) THEN
          IF(JUNIT(13).GT.0.OR.JUNIT(14).GT.0) THEN
            CALL GWT1BAS6INIT(IUNIT,
     *		GZ(LCHNEW),IG(LCIBOU),GX(LCHOLD),GX(LCBOTM),NBOTM,
     *		X(LSTHCK),NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,
     *		IOUTS,ISSFLG(KKPER),X(LSCONC),X(LSRF),X(LSPOR),SBVL,JRF,
     *		X(LSPC),X(LSPR),X(LSPL),X(LSPCON),IX(LSPTID),IX(LSNPCL),
     *		X(LSPNWC),X(LSPNWR),X(LSPNWL),IX(LSLMBO),
     *		NEWPTS,NPMAX,NLIMBO,NP,MOCTYPE,LAYHDT,NIUNIT,
     *		JUNIT,KKPER,NPER,KKSTP,NSTP(KKPER),
     *        IDUM,NMOV,TIMV,NPNTPL,SUMTCH,
cgzh debug ptwt to pt files
     *        Z(LSPTWT),
cgzh varpt
     *       X(LSPCOR),X(LSPROR),X(LSPLOR),NPTLAYA,NPTROWA,NPTCOLA,IDIM,
cgzh SSTR
     *       IPERGWT,X(LSWTFC),JUNIT(24))
          ELSE
            CALL GWT1BAS6INIT(IUNIT,
     *		GZ(LCHNEW),IG(LCIBOU),GX(LCHOLD),GX(LCBOTM),NBOTM,
     *		X(LSTHCK),NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,
     *		IOUTS,ISSFLG(KKPER),X(LSCONC),X(LSRF),X(LSPOR),SBVL,JRF,
     *		X(LSPC),X(LSPR),X(LSPL),X(LSPCON),IX(LSPTID),IX(LSNPCL),
     *		X(LSPNWC),X(LSPNWR),X(LSPNWL),IX(LSLMBO),
     *		NEWPTS,NPMAX,NLIMBO,NP,MOCTYPE,LAYHDT,NIUNIT,
     *		JUNIT,KKPER,NPER,KKSTP,NSTP(KKPER),
     *        IDUM,NMOV,TIMV,NPNTPL,SUMTCH,
cgzh debug ptwt to pt files
     *        Z(LSPTWT),
cgzh varpt
     *       X(LSPCOR),X(LSPROR),X(LSPLOR),
cgzh these three are dummies for nptlaya etc
     *        IX(LSIGNP),IX(LSIGNP),IX(LSIGNP),
c     *     NPTLAYA,NPTROWA,NPTCOLA,
cgzh SSTR
     *       IDIM,IPERGWT,X(LSWTFC),JUNIT(24))
          END IF
C
CGWT----DETERMINE INITIAL PARTICLE WEIGHTS AT BEGINNING OF SIMULATION
c              IF(KKPER.EQ.1.AND.KKSTP.EQ.1) THEN
cgzh debug update pt weights here (before thck update) for TR period
CMOCWT
cgzh SSTR        IF(KKPER.EQ.1.AND.
        IF(KKPER.EQ.IPERGWT.AND.
     *    ISSFLG(KKPER).EQ.0.AND.PTWTON.EQ.1) THEN
c	        WRITE(*,*) 'CALL PTWT1INITWT'
            CALL PTWT1INITWT(IG(LCIBOU),X(LSPC),X(LSPR),X(LSPL),
     *        IX(LSNPCL),Z(LSCELV),Z(LSPTWT),Z(LSSUMW),
     *        NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NPMAX,NP,IOUTS)
        END IF
C
c only at very beginning (finding init mass)
c orig           IF(JUNIT(10).GT.0) 
cgzh SSTR            IF(JUNIT(10).GT.0.AND.KKPER.EQ.1) 
            IF(JUNIT(10).GT.0.AND.KKPER.EQ.IPERGWT) 
     *        CALL GWT1DP6INIT(SBVL,IG(LCIBOU),X(LSTHCK),
     *          X(LSDPCON),X(LSDPRAT),X(LSDPPOR),X(LSDPZO),X(LSDPFO),
     *          IDPZO,IDPFO,NSCOL,NSROW,NSLAY,
     *          NCOL,NROW,NLAY,JUNIT(10),KKPER,IOUT,IOUTS,NIUNIT)
        END IF
      END IF
C
C7C3----DETERMINE WHICH OUTPUT IS NEEDED.
            CALL GWF1BAS6OC(NSTP(KKPER),KKSTP,ICNVG,IR(LCIOFL),NLAY,
     1                      IBUDFL,ICBCFL,IHDDFL,IUNIT(12),IOUT,KKPER,
     2                      IPEROC,ITSOC,IBDOPT,IXSEC,IFREFM,RESETDD,
     3                      RESETDDNEXT)
C
C7C4----CALCULATE BUDGET TERMS. SAVE CELL-BY-CELL FLOW TERMS.
            MSUM = 1
C7C4A---THE ORIGINAL BCF BUDGET MODULE HAS BEEN REPLACED BY THREE
C7C4A---SUBMODULES: SGWF1BCF6S, SGWF1BCF6F, AND SGWF1BCF6B .
            IF (IUNIT(1).GT.0) THEN
              CALL SGWF1BCF6S(VBNM,VBVL,MSUM,GZ(LCHNEW),IG(LCIBOU),
     &                        GX(LCHOLD),RX(LCSC1),GX(LCBOTM),NBOTM,
     &                        RX(LCSC2),DELT,ISSFLG(KKPER),NCOL,NROW,
     &                        NLAY,KKSTP,KKPER,IBCFCB,ICBCFL,GX(LCBUFF),
     &                        IOUT,PERTIM,TOTIM)
CGWT----UPDATE FLUID STORAGE TERMS FOR TRANSPORT (MOC AND MOCIMP)
cea			IF(MOCTYPE.NE.3.AND.ISSFLG(KKPER).EQ.0.AND.IGWTON.GT.0)
			IF(ISSFLG(KKPER).EQ.0.AND.IGWTON.GT.0)
     *		   CALL SGWT1BAS6UP(GZ(LCHNEW),IG(LCIBOU),GX(LCHOLD),
     *           GX(LCBOTM),NBOTM,
     *		   RX(LCSC1),RX(LCSC2),X(LSTHCK),X(LSPOR),GX(LCBUFF),DELT,
     *		   NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,IOUTS,LAYHDT,
     *		   MOCTYPE,GX(LCDELR),GX(LCDELC),
     *           X(LSRW),X(LSRHSE),X(LSCONC))
C
              CALL SGWF1BCF6F(VBNM,VBVL,MSUM,GZ(LCHNEW),IG(LCIBOU),
     &                        GX(LCCR),GX(LCCC),GX(LCCV),DELT,NCOL,NROW,
     &                        NLAY,KKSTP,KKPER,IBCFCB,GX(LCBUFF),IOUT,
     &                        ICBCFL,PERTIM,TOTIM,GX(LCBOTM),NBOTM,
     &                        ICHFLG)
CGWT----FOR TRANSPORT, INITIALIZE SINK/SOURCE ARRAYS
			IF(IGWTON.GT.0)
     *		   CALL GWT1SRC6FM(X(LSEVTF),X(LSCHDF),
     *	       GX(LCBUFF),IG(LCIBOU),X(LSSRCF),X(LSSRCS),X(LSSNKF),
     *           X(LSCHBC),SRCDCY,NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *           IOUTS,IFXHED,ICONLY,MOCTYPE)
CGWT----FOR TRANSPORT AND FHB, COMPUTE TERMS AT FIXED HEADS
			IF(IGWTON.GT.0.AND.IUNIT(16).GT.0)
     *		   CALL GWT1FHBH1FM(IR(LCHDLC),RX(LCBDHV),NHED,NFHBX2,
     *           IFHBHC,X(LSCHBC),NCOL,NROW,NLAY,KKSTP,KKPER,IPERGWT,
     *		   NSCOL,NSROW,NSLAY,IOUTS,ICONLY)
C
              IBDRET=0
              IC1=1
              IC2=NCOL
              IR1=1
              IR2=NROW
              IL1=1
              IL2=NLAY
CGWT----CALCULATE FACE FLUXES FOR TRANSPORT ONLY ON SUBGRID
			IF(IGWTON.GT.0) THEN
				IBDRET=1
				IC1=ISCOL1
				IC2=ISCOL2
				IR1=ISROW1
				IR2=ISROW2
				IL1=ISLAY1
				IL2=ISLAY2
			END IF
C
              DO 37 IDIR = 1, 3
                CALL SGWF1BCF6B(GZ(LCHNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                          GX(LCCV),NCOL,NROW,NLAY,KKSTP,KKPER,
     &                          IBCFCB,GX(LCBUFF),IOUT,ICBCFL,DELT,
     &                          PERTIM,TOTIM,IDIR,IBDRET,ICHFLG,IC1,IC2,
     &                          IR1,IR2,IL1,IL2,GX(LCBOTM),NBOTM)
CGWT----CALCULATE FLOWS ACROSS BOUNDARIES OF SUBGRID
C       CALCULATE VELOCITIES
			  IF(IGWTON.GT.0) 
     *			CALL GWT1VEL6(X(LSCTCF),GX(LCBUFF),IX(LSLBDY),
     *			IDIR,ID,NCOL,NROW,NLAY,IC1,IC2,IR1,IR2,IL1,
     *			IL2,NFACES,NCINFL,
     *			X(LSTHCK),X(LSPOR),X(LSRF),X(LSVC),X(LSVR),X(LSVL),
     *			IG(LCIBOU),NSCOL,NSROW,NSLAY,IOUTS,NTIMV,
     *			VCMAX,VRMAX,VLMAX,TLMIN,
     *			MAXVCJ,MAXVCI,MAXVCK,MAXVRJ,MAXVRI,MAXVRK,
     *			MAXVLJ,MAXVLI,MAXVLK,
     *			DELT,KKSTP,KKPER,NSTP(KKPER),NPER,TIMV,TOTIM,ITCD,
     *			NPNTVL,IVELFM,JUNIT,SUMTCH,
     *			X(LSDKZO),X(LSDKFO),X(LSDKZS),
     *			X(LSDKFS),IDKZO,IDKFO,IDKZS,IDKFS,
     *			MOCTYPE,TCMIN,TRMIN,GX(LCDELR),GX(LCDELC),NIUNIT,
cgzh vlwt
     *            LAYHDT,GZ(LCHNEW),GX(LCHOLD),GX(LCBOTM),NBOTM,
     *            ISSFLG(KKPER))
C
   37         CONTINUE
            ENDIF
            IF(IUNIT(23).GT.0) THEN                                     
              CALL SGWF1LPF1S(VBNM,VBVL,MSUM,GZ(LCHNEW),IG(LCIBOU),     
     &                        GX(LCHOLD),X(LCSC1),GX(LCBOTM),X(LCSC2),  
     &                        DELT,ISSFLG(KKPER),NCOL,NROW,NLAY,KKSTP,
     &                        KKPER,ILPFCB,ICBCFL,GX(LCBUFF),IOUT,
     &                        PERTIM,TOTIM,NBOTM)
CGWT----UPDATE FLUID STORAGE TERMS FOR TRANSPORT (MOC AND MOCIMP)
cea			IF(MOCTYPE.NE.3.AND.ISSFLG(KKPER).EQ.0.AND.IGWTON.GT.0) 
			IF(ISSFLG(KKPER).EQ.0.AND.IGWTON.GT.0) 
     *		  CALL SGWT1BAS6UP(GZ(LCHNEW),IG(LCIBOU),GX(LCHOLD),
     *            GX(LCBOTM),NBOTM,
     *		    X(LCSC1),X(LCSC2),X(LSTHCK),X(LSPOR),GX(LCBUFF),DELT,
     *		    NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,IOUTS,LAYHDT,MOCTYPE,
     *		    GX(LCDELR),GX(LCDELC),X(LSRW),X(LSRHSE),X(LSCONC))
C
              CALL SGWF1LPF1F(VBNM,VBVL,MSUM,GZ(LCHNEW),IG(LCIBOU),
     &                        GX(LCCR),GX(LCCC),GX(LCCV),GX(LCBOTM),
     &                        DELT,NCOL,NROW,NLAY,KKSTP,KKPER,ILPFCB,
     &                        GX(LCBUFF),IOUT,ICBCFL,PERTIM,TOTIM,
     &                        NBOTM,ICHFLG)
CGWT----FOR TRANSPORT, INITIALIZE SINK/SOURCE ARRAYS
			IF(IGWTON.GT.0)
     *		   CALL GWT1SRC6FM(X(LSEVTF),X(LSCHDF),
     *	       GX(LCBUFF),IG(LCIBOU),X(LSSRCF),X(LSSRCS),X(LSSNKF),
     *           X(LSCHBC),SRCDCY,NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *           IOUTS,IFXHED,ICONLY,MOCTYPE)
CGWT----FOR TRANSPORT AND FHB, COMPUTE TERMS AT FIXED HEADS
			IF(IGWTON.GT.0.AND.IUNIT(16).GT.0)
     *		   CALL GWT1FHBH1FM(IR(LCHDLC),RX(LCBDHV),NHED,NFHBX2,
     *           IFHBHC,X(LSCHBC),NCOL,NROW,NLAY,KKSTP,KKPER,IPERGWT,
     *		   NSCOL,NSROW,NSLAY,IOUTS,ICONLY)
C
              IBDRET=0
              IC1=1
              IC2=NCOL
              IR1=1
              IR2=NROW
              IL1=1
              IL2=NLAY
CGWT----CALCULATE FACE FLUXES FOR TRANSPORT ONLY ON SUBGRID
			IF(IGWTON.GT.0) THEN
				IBDRET=1
				IC1=ISCOL1
				IC2=ISCOL2
				IR1=ISROW1
				IR2=ISROW2
				IL1=ISLAY1
				IL2=ISLAY2
			END IF
              DO 157 IDIR=1,3
                CALL SGWF1LPF1B(GZ(LCHNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                          GX(LCCV),GX(LCBOTM),NCOL,NROW,NLAY,
     &                          KKSTP,KKPER,ILPFCB,GX(LCBUFF),IOUT,
     &                          ICBCFL,DELT,PERTIM,TOTIM,IDIR,IBDRET,
     &                          ICHFLG,IC1,IC2,IR1,IR2,IL1,IL2,NBOTM)
CGWT----CALCULATE FLOWS ACROSS BOUNDARIES OF SUBGRID
C       CALCULATE VELOCITIES
			  IF(IGWTON.GT.0) 
     *			CALL GWT1VEL6(X(LSCTCF),GX(LCBUFF),IX(LSLBDY),
     *			IDIR,ID,NCOL,NROW,NLAY,IC1,IC2,IR1,IR2,IL1,
     *			IL2,NFACES,NCINFL,
     *			X(LSTHCK),X(LSPOR),X(LSRF),X(LSVC),X(LSVR),X(LSVL),
     *			IG(LCIBOU),NSCOL,NSROW,NSLAY,IOUTS,NTIMV,
     *			VCMAX,VRMAX,VLMAX,TLMIN,
     *			MAXVCJ,MAXVCI,MAXVCK,MAXVRJ,MAXVRI,MAXVRK,
     *			MAXVLJ,MAXVLI,MAXVLK,
     *			DELT,KKSTP,KKPER,NSTP(KKPER),NPER,TIMV,TOTIM,ITCD,
     *			NPNTVL,IVELFM,JUNIT,SUMTCH,
     *			X(LSDKZO),X(LSDKFO),X(LSDKZS),
     *			X(LSDKFS),IDKZO,IDKFO,IDKZS,IDKFS,
     *			MOCTYPE,TCMIN,TRMIN,GX(LCDELR),GX(LCDELC),NIUNIT,
cgzh vlwt
     *            LAYHDT,GZ(LCHNEW),GX(LCHOLD),GX(LCBOTM),NBOTM,
     *            ISSFLG(KKPER))
157           CONTINUE
            ENDIF
            IF(IUNIT(37).GT.0) THEN
              CALL SGWF1HUF2S(VBNM,VBVL,MSUM,GZ(LCHNEW),IG(LCIBOU),
     &                        GX(LCHOLD),X(LCSC1),GX(LCBOTM),DELT,
     &                        ISSFLG(KKPER),NCOL,NROW,NLAY,KKSTP,KKPER,
     &                        IHUFCB,ICBCFL,GX(LCBUFF),IOUT,PERTIM,
     &                        TOTIM,NBOTM,X(LCHUFTHK),NHUF,IG(LCIZON),
     &                        NZONAR,GX(LCRMLT),NMLTAR,GX(LCDELR),
     &                        GX(LCDELC))
CGWT----UPDATE FLUID STORAGE TERMS FOR TRANSPORT (MOC AND MOCIMP)
cgzh SSTR			IF(MOCTYPE.NE.3.AND.ISSFLG(KKPER).EQ.0.AND.IUNIT(15).GT.0)
			IF(ISSFLG(KKPER).EQ.0.AND.IGWTON.GT.0)
     *		   CALL SGWT1BAS6UP(GZ(LCHNEW),IG(LCIBOU),GX(LCHOLD),
     *           GX(LCBOTM),NBOTM,
     *		   X(LCSC1),X(LCSC2),X(LSTHCK),X(LSPOR),GX(LCBUFF),DELT,
     *		   NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,IOUTS,LAYHDT,
     *		   MOCTYPE,GX(LCDELR),GX(LCDELC),
     *           X(LSRW),X(LSRHSE),X(LSCONC))
C
              CALL SGWF1HUF2F(VBNM,VBVL,MSUM,GZ(LCHNEW),IG(LCIBOU),
     &                        GX(LCCR),GX(LCCC),GX(LCCV),GX(LCBOTM),
     &                        DELT,NCOL,NROW,NLAY,KKSTP,KKPER,IHUFCB,
     &                        GX(LCBUFF),IOUT,ICBCFL,PERTIM,TOTIM,NBOTM,
     &                        ICHFLG,IUNIT(47),X(LCVDHT))
CGWT----FOR TRANSPORT, INITIALIZE SINK/SOURCE ARRAYS
			IF(IGWTON.GT.0)
     *		   CALL GWT1SRC6FM(X(LSEVTF),X(LSCHDF),
     *	       GX(LCBUFF),IG(LCIBOU),X(LSSRCF),X(LSSRCS),X(LSSNKF),
     *           X(LSCHBC),SRCDCY,NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *           IOUTS,IFXHED,ICONLY,MOCTYPE)
CGWT----FOR TRANSPORT AND FHB, COMPUTE TERMS AT FIXED HEADS
			IF(IGWTON.GT.0.AND.IUNIT(16).GT.0)
     *		   CALL GWT1FHBH1FM(IR(LCHDLC),RX(LCBDHV),NHED,NFHBX2,
     *           IFHBHC,X(LSCHBC),NCOL,NROW,NLAY,KKSTP,KKPER,IPERGWT,
     *		   NSCOL,NSROW,NSLAY,IOUTS,ICONLY)
C
              IBDRET=0
              IC1=1
              IC2=NCOL
              IR1=1
              IR2=NROW
              IL1=1
              IL2=NLAY
CGWT----CALCULATE FACE FLUXES FOR TRANSPORT ONLY ON SUBGRID
			IF(IGWTON.GT.0) THEN
				IBDRET=1
				IC1=ISCOL1
				IC2=ISCOL2
				IR1=ISROW1
				IR2=ISROW2
				IL1=ISLAY1
				IL2=ISLAY2
			END IF
              DO 159 IDIR=1,3
                CALL SGWF1HUF2B(GZ(LCHNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                          GX(LCCV),GX(LCBOTM),NCOL,NROW,NLAY,
     &                          KKSTP,KKPER,IHUFCB,GX(LCBUFF),IOUT,
     &                          ICBCFL,DELT,PERTIM,TOTIM,IDIR,IBDRET,
     &                          ICHFLG,IC1,IC2,IR1,IR2,IL1,IL2,NBOTM,
     &                          IUNIT(47),X(LCVDHT))
CGWT----CALCULATE FLOWS ACROSS BOUNDARIES OF SUBGRID
C       CALCULATE VELOCITIES
			  IF(IGWTON.GT.0) 
     *			CALL GWT1VEL6(X(LSCTCF),GX(LCBUFF),IX(LSLBDY),
     *			IDIR,ID,NCOL,NROW,NLAY,IC1,IC2,IR1,IR2,IL1,
     *			IL2,NFACES,NCINFL,
     *			X(LSTHCK),X(LSPOR),X(LSRF),X(LSVC),X(LSVR),X(LSVL),
     *			IG(LCIBOU),NSCOL,NSROW,NSLAY,IOUTS,NTIMV,
     *			VCMAX,VRMAX,VLMAX,TLMIN,
     *			MAXVCJ,MAXVCI,MAXVCK,MAXVRJ,MAXVRI,MAXVRK,
     *			MAXVLJ,MAXVLI,MAXVLK,
     *			DELT,KKSTP,KKPER,NSTP(KKPER),NPER,TIMV,TOTIM,ITCD,
     *			NPNTVL,IVELFM,JUNIT,SUMTCH,
     *			X(LSDKZO),X(LSDKFO),X(LSDKZS),
     *			X(LSDKFS),IDKZO,IDKFO,IDKZS,IDKFS,
     *			MOCTYPE,TCMIN,TRMIN,GX(LCDELR),GX(LCDELC),NIUNIT,
cgzh vlwt
     *            LAYHDT,GZ(LCHNEW),GX(LCHOLD),GX(LCBOTM),NBOTM,
     *            ISSFLG(KKPER))
C
159           CONTINUE
            ENDIF
CGWT----CALCULATE AND PRINT DISPERSION COEFFICIENTS
            IF(IGWTON.GT.0) THEN
               IF(NODISP.NE.1.OR.DIFFUS.NE.0.0) 
     *			CALL GWT1DSP6(X(LSDCC),X(LSDCR),X(LSDCL),
     *            X(LSDRR),X(LSDRC),X(LSDRL),X(LSDLL),X(LSDLC),X(LSDLR),
     *            X(LSTHCK),X(LSALNG),X(LSATH),X(LSATV),
     *			X(LSPOR),X(LSRF),X(LSVC),X(LSVR),X(LSVL),IG(LCIBOU),
     *			MOCTYPE,DELT,NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *			IOUTS,NTIMD,TIMDC,NODISP,DIFFUS,TIMV,
     *			KKSTP,NSTP(KKPER),KKPER,NPER,IDSPFM,NPNTDL,
     *			GX(LCDELR),GX(LCDELC),
     *            IUNIT(21),RX(LCHFB),NHFB,MXHFB,RX(LSHFBD),JUNIT(12))
            END IF
cea
            IF(MOCTYPE.EQ.3.AND.(KKSTP.NE.1.OR.KKPER.NE.1)) THEN
              IF((KKSTP.EQ.1.AND.
     *           (KKPER.GT.IPERGWT.AND.ISSFLG(KKPER).EQ.1).OR.
     *           ISSFLG(KKPER).EQ.0))
     *         CALL GWT1ELLUP6(X(LSDIST),NODESS,
     *           NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,JRF,
     *           GX(LCDELR),GX(LCDELC),X(LSTHCK),IG(LCIBOU),
     *           X(LSCOLD),X(LSCB),X(LSRHSE),X(LSRHSO),
     *           NTFACE,X(LSVC),X(LSVR),X(LSVL),X(LSRF),
     *           X(LSPOR),TIMV,VCMAX,VRMAX,VLMAX,IOUTS,
     *           IX(LSLB),IX(LSNZIN),NFACES,IX(LSNONU),X(LSSAV),
     *           X(LSXFOR),X(LSXBAC),X(LSYFOR),X(LSYBAC),
     *           X(LSCINF),NCINFL,IABOVE,IBELOW,DECAY,IDIM,NRFLG,
     *           IX(LSIACT),X(LSCFOR),X(LSRFOR),X(LSTFOR),X(LSCONL),
cea     *           X(LSVOL),NODES,X(LSRW),IX(LSJAS),IX(LSIAS),
     *           X(LSVOL),NODES,X(LSRW),X(LSRW1),IX(LSJAS),IX(LSIAS),
     *           X(LSAS),NELTS,LENW,LENIW,KKSTP,KKPER,SBVL,
     *           X(LSDIDS),X(LSCONC),IX(LSIW),NMOV,
     *           X(LSDO),IX(LSIDB),NIUNIT,DAGE,JUNIT(9),AGER8,DELT,
     *           IPERGWT)
            ENDIF
CGWT----APPLY CONSTANT HEAD BOUNDARY FLUXES, BFLX PACKAGE
            IF(IGWTON.GT.0.AND.JUNIT(15).GT.0)
     *		   CALL GWT1BFLX5FM(X(LSCHDF),IG(LCIBOU),
     *           NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *           IOUTS,IFXHED,ICONLY,
     *           IX(LSBFCH),NCHNDS,X(LSVC),X(LSVR),X(LSVL),X(LSTHCK))
C--WELLS
            IF (IUNIT(2).GT.0)
     &          CALL GWF1WEL6BD(NWELLS,MXWELL,VBNM,VBVL,MSUM,RX(LCWELL),
     &                          IG(LCIBOU),DELT,NCOL,NROW,NLAY,KKSTP,
     &                          KKPER,IWELCB,ICBCFL,GX(LCBUFF),IOUT,
     &                          PERTIM,TOTIM,NWELVL,IWELAL,IAUXSV)
CGWT----SAVE WELL FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(2).GT.0.AND.IGWTON.GT.0)
     *		  CALL GWT1WEL5BD(NWELLS,MXWELL,RX(LCWELL),IG(LCIBOU),
     *                          NCOL,NROW,NLAY,KKSTP,KKPER,IPERGWT,
     *						  X(LSSRCF),X(LSSRCS),X(LSSNKF),
     *						  NSCOL,NSROW,NSLAY,IOUTS,NWELVL,
     *                          IWELAL,IWELLC,ICONLY,NOPRWL,
     *                          X(LSVC),X(LSVR),X(LSVL),IFACWEL)
C--DRAINS
            IF (IUNIT(3).GT.0)
     &          CALL GWF1DRN6BD(NDRAIN,MXDRN,VBNM,VBVL,MSUM,RX(LCDRAI),
     &                          DELT,GZ(LCHNEW),NCOL,NROW,NLAY,
     &                          IG(LCIBOU),KKSTP,KKPER,IDRNCB,ICBCFL,
     &                          GX(LCBUFF),IOUT,PERTIM,TOTIM,NDRNVL,
     &                          IDRNAL,IAUXSV)
CGWT----SAVE DRAIN FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(3).GT.0.AND.IGWTON.GT.0)
     *		  CALL GWT1DRN5BD(NDRAIN,MXDRN,RX(LCDRAI),IG(LCIBOU),
     *						  NCOL,NROW,NLAY,X(LSSNKF),
     *                          NSCOL,NSROW,NSLAY,
     *						  IOUTS,NDRNVL,IDRNAL,ICONLY,KKSTP,KKPER,
     *                          NOPRDR,
     *                          X(LSVC),X(LSVR),X(LSVL),IPERGWT,IFACDRN)
C--RIVERS
            IF (IUNIT(4).GT.0)
     &          CALL GWF1RIV6BD(NRIVER,MXRIVR,RX(LCRIVR),IG(LCIBOU),
     &                          GZ(LCHNEW),NCOL,NROW,NLAY,DELT,VBVL,
     &                          VBNM,MSUM,KKSTP,KKPER,IRIVCB,ICBCFL,
     &                          GX(LCBUFF),IOUT,PERTIM,TOTIM,NRIVVL,
     &                          IRIVAL,IAUXSV)
CGWT----SAVE RIVER FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(4).GT.0.AND.IGWTON.GT.0)
     *		  CALL GWT1RIV5BD(NRIVER,MXRIVR,RX(LCRIVR),IG(LCIBOU),
     *						  NCOL,NROW,NLAY,KKSTP,KKPER,IPERGWT,
     *						  X(LSSRCF),X(LSSRCS),X(LSSNKF),
     *						  NSCOL,NSROW,NSLAY,IOUTS,
     *						  NRIVVL,IRIVAL,IRIVRC,ICONLY,NOPRRV,
     *                          X(LSVC),X(LSVR),X(LSVL),IFACRIV)
C--EVAPOTRANSPIRATION
CLAK
            IF (IUNIT(5).GT.0) THEN
              IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3)
     1                CALL GWF1LAK3ST(0,NCOL,NROW,NLAY,IG(LCIBOU),
     2                            LKNODE,IR(ICLAKE),MXLKND,NLAKES,
     3                            RX(ISTGLD),GX(LCBOTM),NBOTM)
              CALL GWF1EVT6BD(NEVTOP,IR(LCIEVT),RX(LCEVTR),RX(LCEXDP),
     &                          RX(LCSURF),IG(LCIBOU),GZ(LCHNEW),NCOL,
     &                          NROW,NLAY,DELT,VBVL,VBNM,MSUM,KKSTP,
     &                          KKPER,IEVTCB,ICBCFL,GX(LCBUFF),IOUT,
     &                          PERTIM,TOTIM)
              IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3)
     1                CALL GWF1LAK3ST(1,NCOL,NROW,NLAY,IG(LCIBOU),
     2                            LKNODE,IR(ICLAKE),MXLKND,NLAKES,
     3                            RX(ISTGLD),GX(LCBOTM),NBOTM)
            END IF
CGWT----SAVE ET FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(5).GT.0.AND.IGWTON.GT.0)
     *		  CALL GWT1EVT6BD(NEVTOP,GX(LCBUFF),X(LSEVTF),IR(LCIEVT),
     *						  IG(LCIBOU),NCOL,NROW,NLAY,
     *						  NSCOL,NSROW,NSLAY,IOUTS,ICONLY,
     *                          X(LSVL),IEVTTP)
C--GENERAL-HEAD BOUNDARIES
            IF (IUNIT(7).GT.0)
     &          CALL GWF1GHB6BD(NBOUND,MXBND,VBNM,VBVL,MSUM,RX(LCBNDS),
     &                          DELT,GZ(LCHNEW),NCOL,NROW,NLAY,
     &                          IG(LCIBOU),KKSTP,KKPER,IGHBCB,ICBCFL,
     &                          GX(LCBUFF),IOUT,PERTIM,TOTIM,NGHBVL,
     &                          IGHBAL,IAUXSV)
CGWT----SAVE GHB FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(7).GT.0.AND.IGWTON.GT.0)
     *		  CALL GWT1GHB5BD(NBOUND,MXBND,RX(LCBNDS),IG(LCIBOU),
     *						  NCOL,NROW,NLAY,KKSTP,KKPER,IPERGWT,
     *						  X(LSSRCF),X(LSSRCS),X(LSSNKF),
     *						  NSCOL,NSROW,NSLAY,IOUTS,
     *						  NGHBVL,IGHBAL,IBNDSC,ICONLY,NOPRGB,
     *                          X(LSVC),X(LSVR),X(LSVL),IFACGHB)
C--RECHARGE
CLAK
            IF (IUNIT(8).GT.0) THEN
              IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3)
     1                CALL GWF1LAK3ST(0,NCOL,NROW,NLAY,IG(LCIBOU),
     2                            LKNODE,IR(ICLAKE),MXLKND,NLAKES,
     3                            RX(ISTGLD),GX(LCBOTM),NBOTM)
              CALL GWF1RCH6BD(NRCHOP,IR(LCIRCH),RX(LCRECH),IG(LCIBOU),
     &                          NROW,NCOL,NLAY,DELT,VBVL,VBNM,MSUM,
     &                          KKSTP,KKPER,IRCHCB,ICBCFL,GX(LCBUFF),
     &                          IOUT,PERTIM,TOTIM)
              IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3)
     1                CALL GWF1LAK3ST(1,NCOL,NROW,NLAY,IG(LCIBOU),
     2                            LKNODE,IR(ICLAKE),MXLKND,NLAKES,
     3                            RX(ISTGLD),GX(LCBOTM),NBOTM)
            END IF
CGWT----SAVE RECHARGE FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(8).GT.0.AND.IGWTON.GT.0.AND.JUNIT(1).GT.0)
     *		  CALL GWT1RCH5BD(NRCHOP,IR(LCIRCH),RX(LSCRCH),GX(LCBUFF),
     *					IG(LCIBOU),X(LSSRCF),
     *                    X(LSSRCS),X(LSSNKF),
     *   		        	NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *					IOUTS,JUNIT(1),KKSTP,ICONLY,
     *                    X(LSVL),IRCHTP)
C--Specified-Flow and Specified-Head Boundary
            IF(IUNIT(16).GT.0)
     &          CALL GWF1FHB1BD(IR(LCFLLC),RX(LCBDFV),NFLW,VBNM,VBVL,
     &                      MSUM,IG(LCIBOU),DELT,NCOL,NROW,NLAY,KKSTP,
     &                      KKPER,IFHBCB,ICBCFL,PERTIM,TOTIM,
     &                      GX(LCBUFF),IOUT,IFHBD4,NFLWDIM)
CGWT----SAVE FLOW-HEAD BOUNDARIES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(16).GT.0.AND.IGWTON.GT.0)
     *		  CALL GWT1FHBF1BD(IR(LCFLLC),RX(LCBDFV),NFLW,IFHBD4,
     *                      NFHBX1,IFHBFC,IG(LCIBOU),
     *                      NCOL,NROW,NLAY,KKSTP,KKPER,IPERGWT,
     *					  X(LSSRCF),X(LSSRCS),X(LSSNKF),
     *					  NSCOL,NSROW,NSLAY,IOUTS,ICONLY,
     *                       X(LSVC),X(LSVR),X(LSVL),IFACFHB)
C--RESERVOIRS
            IF (IUNIT(17).GT.0)
     &          CALL GWF1RES1BD(IR(LCIRES),IR(LCIRSL),RX(LCBRES),
     &                      RX(LCCRES),RX(LCBBRE),RX(LCHRES),IG(LCIBOU),
     &                      GZ(LCHNEW),GX(LCBUFF),VBVL,VBNM,MSUM,KKSTP,
     &                      KKPER,NRES,NRESOP,NCOL,NROW,NLAY,DELT,
     &                      IRESCB,ICBCFL,IOUT)
C--STREAM-AQUIFER RELATIONS (STR6 PACKAGE)
            IF(IUNIT(18).GT.0)
     &          CALL GWF1STR6BD(NSTREM,RX(LCSTRM_),IR(ICSTRM_),
     &                      IG(LCIBOU),MXSTRM,GZ(LCHNEW),NCOL,NROW,
     &                      NLAY,DELT,VBVL,VBNM,MSUM,KKSTP,KKPER,
     &                      ISTCB1STR6,ISTCB2STR6,ICBCFL,GX(LCBUFF),
     &                      IOUT,NTRIB,NSSSTR6,RX(LCTRIB),IR(LCTBAR),
     &                      IR(LCIVAR_),IR(LCFGAR),ICALC,CONSTSTR6,
     &                      IPTFLG)
C--INTERBED STORAGE
            IF (IUNIT(19).GT.0)
     1          CALL GWF1IBS6BD(IG(LCIBOU),GZ(LCHNEW),GX(LCHOLD),
     2                          RX(LCHC),RX(LCSCE),RX(LCSCV),RX(LCSUB),
     3                          GX(LCDELR),GX(LCDELC),NCOL,NROW,NLAY,
     4                          DELT,VBVL,VBNM,MSUM,KSTP,KPER,IIBSCB,
     5                          ICBCFL,GX(LCBUFF),IOUT,ISSFLG(KKPER),
     6                          IBSDIM)
            IF(IUNIT(54).GT.0)
     1          CALL GWF1SUB1BD(IG(LCIBOU),GZ(LCHNEW),GX(LCHOLD),
     2                          GX(LCBUFF),GX(LCDELR),GX(LCDELC),VBVL,
     3                          VBNM,NN,NND1,ND1,ND2,NDB,NNDB,NPZ,NCOL,
     4                          NROW,NLAY,NODES,DELT,MSUM,NIUNIT,KKSTP,
     5                          KKPER,ISUBCB,ICBCFL,ISSFLG(KKPER),IOUT)
            IF(IUNIT(55).GT.0)
     1          CALL GWF1SWT1BD(GX(LCBUFF),IG(LCIBOU),GZ(LCHNEW),
     2                          GX(LCHOLD),GX(LCBOTM),GX(LCDELR),
     3                          GX(LCDELC),NCOL,NROW,NLAY,DELT,VBVL,
     4                          VBNM,MSUM,KSTP,KPER,ICBCFL,
     5                          NIUNIT,ISSFLG(KKPER),IOUT)
C--STREAM-AQUIFER RELATIONS (SFR2 PACKAGE)
C-----ADDED DELEAK MAY 12, 2004
Cdep  Changed SFR1 to SFR2
Cdep  Added new variable to end of list for specified lake outflow
Cdep   June 6, 2006
            IF (IUNIT(44).GT.0)
     1          CALL GWF1SFR2BD(RX(LCSTRM),IR(ICSTRM),GZ(LCHNEW),
     2       IG(LCIBOU),NSTRM,NCOL,NROW,NLAY,NSS,RX(LCSEG),IR(ICSEG),
     3       IR(LCOTSG),RX(LCXSEC),IR(LCIVAR),CONST,MAXPTS,RX(LCQSTG),
     4       RX(LCOTFLW),RX(LCDVFLW),NLAKESAR,RX(ISTGLD),RX(ISTGNW),
     5       RX(LKACC7),RX(ISTRIN),RX(ISTROT),THETA,DELT,KKSTP,KKPER,
     6       VBVL,VBNM,MSUM,ISTCB1,ISTCB2,ICBCFL,GX(LCBUFF),IOUT,IPTFLG,
     7       IUNIT(15),IUNIT(46),IR(LSGAGE),NUMGAGE,PERTIM,TOTIM,
     &       RX(LCSFRQ),NSEGDIM,IUNIT(22),DLEAK,GX(LCHOLD),
     9       RZ(LCSUZDPST),RZ(LCSUZTHST),RZ(LCSUZSPST),RZ(LCSUZFLST),
     &       RZ(LCSEPS),RZ(LCSTHTR),RZ(LCSTHTS),IR(ICSLTRLST),
     &       IR(ICSITRLST),IR(ICSTRLHLD),RZ(LCSUZFLWT),RZ(LCSUZSTOR),
     &       RZ(LCSDELSTR),IR(ICSNWAVST),NSTRAIL,NSTOTRL,IUZT,ISUZN,
     &       NUZST,RZ(LCSOLSFLX),RZ(LCSUZWDTH),IR(ICSLOOP),RX(LCUHC),
     &       IR(ICELEV),RX(LCDPTH),RZ(LCWETP),NSFRSETS,RZ(LCSUZSEEP),
     &       RZ(LCOLDFLBT),NUMCELL,NUZROW,NUZCOL,RX(LCAVWAT),
     &       RX(LCWAT1),NUMAVE,RX(LCAVDPT),RX(LCSFRUZBD),ISSFLG(KKPER),
     &       IBUDFL,RX(IDSTRT))
Cdep  End of change
CGWT----SFR
      IF(IUNIT(44).GT.0.AND.ICNVG.EQ.1.AND.IGWTON.GT.0)
     * CALL GWT1SFR2BD(NSTRM,RX(LCSTRM),IR(ICSTRM),IUNIT(22),RX(LCSEG),
     1  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NSOL,ISTCB1,ICBCFL,
     2  IOUTS,NSS,NSEGDIM,RX(LCOTFLW),IR(LCIVAR),IPTFLG,IR(LCOTSG),
     3  X(LSSRCF),X(LSSRCS),X(LSSNKF),X(LSCONC),RX(LSCPPT),
     4  IR(ICSEG),RX(LSCOUT),RX(LSCONQ),RX(LSCNRN),
     5  RX(LSCNTB),GX(LCBUFF),RX(LSCQIN),RX(LSCGW),RX(LSLAKE),NLAKES,
     *  NLAKESAR)
C--LAKES
Cdep  Added SURFDEPTH to end of call statement  March 23, 2009
            IF (IUNIT(22).GT.0)
     &          CALL GWF1LAK3BD(LKNODE,MXLKND,NODES,IR(ICLAKE),
     &              GZ(LCHNEW),IG(LCIBOU),NCOL,NROW,NLAY,NLAKES,
     &              DELT,NSSAR,NTRB,NDV,IR(INTRB),IR(INDV),RX(ISTRIN),
     &              RX(ISTROT),RX(ISTGLD),RX(ISTGNW),RX(LCCNDF),
     &              RX(LCLKPR),RX(LCLKEV),GX(LCDELR),GX(LCDELC),
     &              GX(LCBOTM),NBOTM,VBVL,VBNM,MSUM,KSTP,KPER,
     &              ILKCB,ICBCFL,GX(LCBUFF),IOUT,RX(LCSTAG),
     &              PERTIM,TOTIM,IR(IICS),IR(IISUB),RX(ISILL),
     &              ICMX,NCLS,RX(LCWTDR),LWRT,RZ(LKACC1),
     &              RZ(LKACC2),RZ(LKACC3),RX(LKACC4),RX(LKACC5),
     &              RX(LKACC6),RX(LKACC7),RX(LKACC8),RX(LKACC9),
     &              RX(LKACC10),RX(LKACC11),IR(LKDRY),RX(IBTMS),
     &              IR(LKNCNT),IR(LKKSUB),RX(LKSADJ),RX(LKFLXI),
     &              IR(LKNCNS),RX(LKSVT),IR(LKJCLS),RX(LCRNF),
     &              THETA,RX(LKCNN),RX(LKCHN),RX(IAREN),
     &              IUNIT(15),KCNT,IR(IMSUB),IR(IMSUB1),
     &              IUNIT(46),NUMGAGE,IR(LSGAGE),RX(LSOVOL),
     &              RX(LSFLOB),ISSFLG(KPER),LAYHDT,
     &              IAUXSV,RX(LKVI),RX(ISTGLD2),RX(LKCUM1),
     &              RX(LKCUM2),RX(LKCUM3),RX(LKCUM4),RX(LKCUM5),
     &              RX(LKCUM6),RX(LKCUM7),RX(LKCUM8),RX(LKCUM9),
     &              HNOFLO,RX(IAREATAB),RX(IDPTHTAB),RX(IDSTRT),
     *              RX(LCEVAPO),RX(LCFLWIN),RX(LCWITDW),RX(LCGWRAT),
     *              NSSLK,RZ(IDLKSTAGE),RZ(ISLKOTFLW),SURFDEPTH)
CGWT----SAVE LAKE FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
cgzh orig            IF(IUNIT(22).GT.0.AND.IUNIT(15).GT.0) 
cgzh orig     *          CALL GWT1LAK3BD(NLAKES,LKNODE,IR(ICLAKE),MXLKND,
            IF(IGWTON.GT.0) 
     *        CALL GWT1LAK3BD(NLAKES,LKNODE,IR(ICLAKE),MXLKND,IUNIT(22),
     1                      NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NSOL,IOUTS,
     2                      X(LSSRCF),X(LSSRCS),X(LSSNKF),RX(LSLAKE),
     3                      RX(LSRTCO),RX(LSFLOB),IX(LSIGLK),GX(LCBUFF),
     3                      JUNIT(15),IG(LCIBOU),X(LSVC),X(LSVR),
     4                      X(LSVL),NLAKESAR)
C--EVAPOTRANSPIRATION WITH SEGMENTED RATE FUNCTION
            IF (IUNIT(39).GT.0)
     &          CALL GWF1ETS1BD(NETSOP,IR(LCIETS),RX(LCETSR),RX(LCETSX),
     &                          RX(LCETSS),IG(LCIBOU),GZ(LCHNEW),NCOL,
     &                          NROW,NLAY,DELT,VBVL,VBNM,MSUM,KKSTP,
     &                          KKPER,IETSCB,ICBCFL,GX(LCBUFF),IOUT,
     &                          PERTIM,TOTIM,NETSEG,RX(LCPXDP),
     &                          RX(LCPETM),NSEGAR)
CGWT----SAVE ETS FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(39).GT.0.AND.IGWTON.GT.0)
     *		  CALL GWT1EVT6BD(NETSOP,GX(LCBUFF),X(LSEVTF),IR(LCIETS),
     *						  IG(LCIBOU),NCOL,NROW,NLAY,
     *						  NSCOL,NSROW,NSLAY,IOUTS,ICONLY,
     *                          X(LSVL),IEVTTP)
C--DRAINS WITH RETURN FLOW
            IF (IUNIT(40).GT.0)
     &          CALL GWF1DRT1BD(NDRTCL,MXDRT,VBNM,VBVL,MSUM,RX(LCDRTF),
     &                          DELT,GZ(LCHNEW),NCOL,NROW,NLAY,
     &                          IG(LCIBOU),KKSTP,KKPER,IDRTCB,ICBCFL,
     &                          GX(LCBUFF),IOUT,PERTIM,TOTIM,NDRTVL,
     &                          IDRTAL,IDRTFL,NRFLOW,IAUXSV)
CGWT----SAVE DRT FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(40).GT.0.AND.IGWTON.GT.0)
     &		  CALL GWT1DRT1BD(NDRTCL,MXDRT,RX(LCDRTF),NCOL,NROW,NLAY,
     &                IG(LCIBOU),IOUTS,NDRTVL,IDRTAL,IDRTFL,
     &				X(LSSRCF),X(LSSRCS),X(LSSNKF),X(LSCONC),
     &                NSCOL,NSROW,NSLAY,ICONLY,KKSTP,KKPER,NOPRDT,
     *                X(LSVC),X(LSVR),X(LSVL),IPERGWT,IFACDRT)
C
            IF(IUNIT(51).GT.0)
     1          CALL GWF1DAF1BD(IUNIT(52)+1,IOUT,ITMUNI,DELT,VBVL,
     2                          VBNM,MSUM,KSTP,KPER,IDAFCB,ICBCFL,
     3                          GX(LCBUFF),PERTIM,TOTIM,NCOL,NROW,
     4                          NLAY,IG(LCIBOU))
C--MULTINODE WELLS
            IF(IUNIT(50).GT.0)
     &          CALL GWF1MNW1BD(MNWSITE,NWELL2,MXWEL2,VBNM,VBVL,MSUM,
     &                          DELT,RZ(LCWEL2),IG(LCIBOU),GZ(LCHNEW),
     &                          NCOL,NROW,NODES,NSTP(KKPER),KKSTP,KKPER,
     &                          IMNWCB,ICBCFL,GX(LCBUFF),IOUT,IOWELL2,
     &                          TOTIM,HDRY,PERTIM)
C--MNW2
            IF(IUNIT(57).GT.0) 
     &          CALL GWF1MNW2BD(IMNWCB,ICBCFL,NAUX,KKSTP,KKPER,
     &                          NCOL,NROW,NLAY,nmnw2,RZ(LCMNW2),iout,
     &                          DELT,PERTIM,TOTIM,IG(LCIBOU),
     &                          MNWMAX,NODTOT,msum,HDRY,VBNM,VBVL,
     &                          GX(LCBUFF),RZ(LCMNWN),GZ(LCHNEW),
     &                          WELLID,MNWPRNT,IGWTON,NMNWVL,
     &                          MNWOBS,IUNIT(15))
C
CGWT----SAVE MNW FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(50).GT.0.AND.IGWTON.GT.0)
     &		  CALL GWT1MNW1BD(NWELL2,mxwel2,RZ(LCWEL2),IG(LCIBOU),
     &             ncol,nrow,nodes,
     +             iouts,X(LSCONC),X(LSSRCF),X(LSSRCS),X(LSSNKF),
     *             nscol,nsrow,nslay,IR(LSMNWI),
     *             X(LSSRCM),X(LSSNKM),X(LSSOLM),icomplex,KKSTP,MNWSITE,
     *             KKPER,IR(LSIQZE),JUNIT(18),IPERGWT)
CGWT----SAVE MNW2 FLUXES FOR TRANSPORT MASS-BALANCE CALCULATION
		  IF(IUNIT(57).GT.0.AND.IGWTON.GT.0)
     &		  CALL GWT1MNW2BD(nmnw2,MNWMAX,RZ(LCMNW2),NODTOT,RZ(LCMNWN),
     &             IG(LCIBOU),
     +             ncol,nrow,nlay,WELLID,MNWPRNT,
     +             iouts,X(LSCONC),X(LSSRCF),X(LSSRCS),X(LSSNKF),
     &             nscol,nsrow,nslay,
     +             X(LSSRCM),X(LSSNKM),X(LSSOLM),KKSTP,
     +             kkper,IPERGWT,NMNWVL)
            IF (IUNIT(43).GT.0)
     &          CALL GWF1HYD1OT(GZ,LENGZ,RX,LENRX,IG,LENIG,RX(LCHYDM),
     &                        NUMH,IHYDMUN,TOTIM,HYDNOH,NROW,NCOL,
     &                        ITMUNI,IOUT)
C--MULTINODE WELLS OUTPUT
            IF(IUNIT(50).GT.0.AND.IUNIT(15).EQ.0)
     &          CALL GWF1MNW1BO(MNWSITE,NWELL2,MXWEL2,
     &                          RZ(LCWEL2),IG(LCIBOU),GZ(LCHNEW),
     &                          NCOL,NROW,NODES,NSTP(KKPER),KKSTP,KKPER,
     &                          IMNWCB,ICBCFL,IOUT,IOWELL2,
     &                          TOTIM,HDRY,IUNIT(15))
C--MULTINODE WELLS OUTPUT (MNW2)
            IF(IUNIT(58).GT.0.AND.IUNIT(15).EQ.0)
     &          CALL GWF1MNW2BO(IMNWCB,ICBCFL,NAUX,KKSTP,KKPER,
     &                          NCOL,NROW,NLAY,nmnw2,RZ(LCMNW2),iout,
     &                          DELT,PERTIM,TOTIM,IG(LCIBOU),
     &                          MNWMAX,NODTOT,msum,HDRY,VBNM,VBVL,
     &                          GX(LCBUFF),RZ(LCMNWN),GZ(LCHNEW),
     &                          WELLID,MNWPRNT,IGWTON,MNWOBS,IUNIT(15)
     &                          ,NMNWVL)
CLMT
CLMT----CALL LINK-MT3DMS SUBROUTINES TO SAVE FLOW-TRANSPORT LINK FILE
CLMT----FOR USE BY MT3DMS FOR TRANSPORT SIMULATION
CLMT
            INCLUDE 'lmt6.inc'
CLMT
C
C7C5---PRINT AND/OR SAVE HEAD AND DRAWDOWN MATRICES.
C------PRINT OVERALL WATER BUDGET.
            CALL GWF1BAS6OT(GZ(LCHNEW),GX(LCSTRT),ISTRT,GX(LCBUFF),
     &                      IR(LCIOFL),MSUM,IG(LCIBOU),VBNM,VBVL,KKSTP,
     &                      KKPER,DELT,PERTIM,TOTIM,ITMUNI,NCOL,NROW,
     &                      NLAY,ICNVG,IHDDFL,IBUDFL,IHEDFM,IHEDUN,
     &                      IDDNFM,IDDNUN,IOUT,CHEDFM,CDDNFM,IXSEC,
     &                      LBHDSV,LBDDSV,IBOUUN,LBBOSV,CBOUFM,ISA,
     &                      RESETDD)
            IF (IUNIT(19).GT.0)
     1          CALL GWF1IBS6OT(NCOL,NROW,NLAY,PERTIM,TOTIM,KKSTP,
     2                          KKPER,NSTP(KKPER),GX(LCBUFF),RX(LCSUB),
     3                          RX(LCHC),IIBSOC,ISUBFM,ICOMFM,IHCFM,
     4                          ISUBUN,ICOMUN,IHCUN,IUNIT(19),IOUT,
     5                          ISSFLG(KKPER),IBSDIM)
            IF(IUNIT(54).GT.0)
     1          CALL GWF1SUB1OT(NCOL,NROW,NLAY,PERTIM,TOTIM,KKSTP,
     2                          KKPER,NSTP(KKPER),GX(LCBUFF),NODES,NN,
     3                          ND1,ND2,NND1,NNDB,NDB,ISSFLG(KKPER),
     4                          IUNIT(54),IOUT)
            IF(IUNIT(55).GT.0)
     1          CALL GWF1SWT1OT(NCOL,NROW,NLAY,PERTIM,TOTIM,KKSTP,KKPER,
     1                      NSTP(KKPER),GX(LCBUFF),IOUT)
C------PRINT AND/OR SAVE HEADS INTERPOLATED TO HYDROGEOLOGIC UNITS
            IF(IUNIT(37).GT.0.AND.(IOHUFHDS.NE.0.OR.IOHUFFLWS.NE.0))
     &          CALL GWF1HUF2OT(IOHUFHDS,IOHUFFLWS,GZ(LCHNEW),IHEDFM,
     &                      IG(LCIBOU),NHUF,NCOL,NROW,NLAY,X(LCHUFTHK),
     &                      GX(LCBOTM),NBOTM,GX(LCCV),GX(LCDELR),
     &                      GX(LCDELC),GX(LCRMLT),NMLTAR,IG(LCIZON),
     &                      NZONAR,KKSTP,KKPER,ISA,ICNVG,IOUT,HNOFLO,
     &                      CHEDFM,DELT, PERTIM,TOTIM,X(LCHUFTMP),
     &                      X(LCGS),ICBCFL,ICHFLG)
C
CGWT----SKIP FOLLOWING SECTION IF SOLUTE-TRANSPORT OPTION NOT SELECTED (GWT)
		  IF(IGWTON.LE.0) GO TO 190
CGWT----ALSO SKIP IF FLOW SOLUTION DID NOT CONVERGE
            IF(ICNVG.EQ.0) THEN
              WRITE(IOUTS,*) 
              WRITE(IOUTS,*) '*** FLOW MODEL FAILED TO CONVERGE; ',
     * 'TRANSPORT MODEL NOT CALLED ***'
              GO TO 190
            END IF
C
C       COMPUTE NEXT TIME STEP AND PRINT VELOCITIES

			  IF(IGWTON.GT.0) 
     *			CALL SMOC6V(X(LSTHCK),GX(LCBUFF),X(LSPOR),X(LSRF),
     *			X(LSVC),X(LSVR),X(LSVL),IG(LCIBOU),DELT,
     *            KKSTP,KKPER,NSTP(KKPER),NPER,
     *            NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS,NTIMV,
     *            TIMV,VCMAX,VRMAX,VLMAX,TLMIN,IDIR,TOTIM,
     *            MAXVCJ,MAXVCI,MAXVCK,MAXVRJ,MAXVRI,MAXVRK,
     *            MAXVLJ,MAXVLI,MAXVLK,ITCD,
     *            NPNTVL,IVELFM,JUNIT,SUMTCH,
     *            DKZO,DKFO,DKZS,
     *            DKFS,IDKZO,IDKFO,IDKZS,IDKFS,
cellam
     *			MOCTYPE,TCMIN,TRMIN,GX(LCDELR),GX(LCDELC),NIUNIT)
cellam  
CGWT DP  COMPUTE TIME STEP STABILITY FOR DOUBLE POROSITY (EXPLICIT ONLY)
		  IF(JUNIT(10).GT.0) THEN
		    IF(MOCTYPE.EQ.1) THEN
                 CALL GWT1DP6ST(IG(LCIBOU),X(LSRF),X(LSPOR),
     1			TIMDP,JMAXDP,IMAXDP,KMAXDP,
     2			X(LSDPRAT),X(LSDPPOR),
     3			NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS)
			END IF
		  END IF
CGWT----COMPUTE TRANSPORT TIME STEP (TIMV) AND NUMBER OF MOVES
CGWT----REQUIRED FOR THIS FLOW TIME STEP (NMOV), INITIALIZE TRANSPORT
CGWT----TIME (SUMTCH), ZERO CNCNC BEFORE MOVE LOOP
CGWT----CALCULATE INITIAL MASS
		  CALL GWT1BAS6ST(X(LSSRCF),X(LSSNKF),
     *	    X(LSRF),X(LSPOR),X(LSTHCK),
     *		VCMAX,VRMAX,VLMAX,TLMIN,TIMV,NTIMV,ITCD,
     *		MAXVCJ,MAXVCI,MAXVCK,MAXVRJ,MAXVRI,MAXVRK,MAXVLJ,MAXVLI,
     *        MAXVLK,TIMDC,NTIMD,
     *		NMOV,SUMTCH,NSCOL,NSROW,NSLAY,IOUTS,KKSTP,KKPER,
     *		TOTIM,DELT,NODISP,DIFFUS,DECAY,ICONLY,MOCTYPE,
     *		JUNIT(10),TIMDP,JMAXDP,IMAXDP,KMAXDP,
     *        KKPER,IG(LCIBOU),X(LSCONC),SBVL,NCOL,NROW,NLAY,JRF,NIUNIT)
C
CGWT INITIALIZE, LOAD AND ASSEMBLE (MOCIMP)
	      IF(MOCTYPE.EQ.2) CALL GWT1LDAS6(IG(LCIBOU),NCOL,NROW,NLAY,
     *		X(LSDCC),X(LSDCR),X(LSDCL),X(LSDRR),X(LSDRC),X(LSDRL),
     *		X(LSDLL),X(LSDLC),X(LSDLR),
     *		NSCOL,NSROW,NSLAY,X(LSTHCK),X(LSPOR),X(LSRF),TIMV,
     *		X(LSVA),X(LSVAD),X(LSRHS),IX(LSCIN),IX(LSMRNO),
     *		FDTMTH,NRN,NODESS,IX(LSNPCL),
     *        IX(LSCI),IX(LSCIR),IX(LSCIRH),IX(LSCIRL),
     *        IX(LSIND),EPSSLV,MAXIT,IDIREC,NBN,IX(LSMRNZ),IX(LSXDMA),
     *        ISSIZH)
C
cgzh srcfix
c if srcfac was allocated, and needs to be again, deallocate
              IF(IALFLAG.EQ.1.AND.
     &         (.NOT.(ISSFLG(KKPER).EQ.1.AND.KKSTP.GT.1))) THEN
                 IF(ALLOCATED(SRCFAC)) DEALLOCATE (SRCFAC)
                 IALFLAG=0
              ELSE
              END IF			
CMOCWT  
            IF(PTWTON.EQ.1) THEN
C   
cgzh debug  update only if transient or 1st ts of ss stress period
cgzh debug  same as: (not (SS and KSTP>1))
              IF (.NOT.(ISSFLG(KKPER).EQ.1.AND.KKSTP.GT.1)) THEN
CMOCWT   SUM FLUX ACROSS SUBGRID BOUNDARIES (ONLY WHEN USING WEIGHTED PARTICLES)
c  	          WRITE(*,*) 'CALL SMOC5BYSRC'
                CALL SMOC5BYSRC(X(LSCTCF),IX(LSLBDY),
     *           X(LSBSRC),X(LSBSOL),X(LSBSNK),X(LSCINF),NCINFL,
cgzh cbdy
     *           X(LSCINA),X(LSCINB),X(LSCINXY),
     *           NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NFACES,IABOVE,IBELOW)
CGWT----UPDATE PARAMETERS NEEDED FOR PARTICLE WEIGHTING SCHEME
c	          WRITE(*,*) 'CALL PTWT1UP'
               IF(JUNIT(13).GT.0.OR.JUNIT(14).GT.0) THEN
		      CALL PTWT1UP(IG(LCIBOU),X(LSTHCK),X(LSPOR),
     *           Z(LSCELV),IX(LSIGNP),Z(LSPTWT),
     *           X(LSVC),X(LSVR),X(LSVL),X(LSRF),
     *           X(LSSRCF),X(LSSNKF),X(LSBSRC),
     *           NPTPND,NPMAX,
     *           IX(LSPTID),X(LSPNWC),X(LSPNWR),X(LSPNWL),NEWPTS,
     *           NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS,
     *           KKPER,ISSFLG(KKPER),
     *           Z(LSSUMW),IX(LSNPCL),
     *           X(LSPC),X(LSPR),X(LSPL),X(LSPCON),LAYHDT,NP,
cgzh varpt
     *           JUNIT(13),JUNIT(14),X(LSPCOR),X(LSPROR),X(LSPLOR),
     *           NPTLAYA,NPTROWA,NPTCOLA,
cgzh ssfix
     *           IDIM,X(LSCONC),IX(LSLMBO),NLIMBO,MULTSS,SBVL,NIUNIT,
     *           REMCRIT,GENCRIT,SRCDCY,
cgzh mnw
     &           X(LSSRCM),X(LSSNKM),IUNIT(50),
cgzh srcfix
     &           IX(LSISRC),MAXSRC,NPCELLMAX,IX(LSNPOR),IPERGWT,
     &           X(LSWTFC),
cgzh et
     *           IUNIT(5),NEVTOP,IR(LCIEVT),X(LSEVTF),ISRCFIX,
c RBW sorption     
     *           JRF)
               ELSE
		      CALL PTWT1UP(IG(LCIBOU),X(LSTHCK),X(LSPOR),
     *           Z(LSCELV),IX(LSIGNP),Z(LSPTWT),
     *           X(LSVC),X(LSVR),X(LSVL),X(LSRF),
     *           X(LSSRCF),X(LSSNKF),X(LSBSRC),
     *           NPTPND,NPMAX,
     *           IX(LSPTID),X(LSPNWC),X(LSPNWR),X(LSPNWL),NEWPTS,
     *           NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS,
     *           KKPER,ISSFLG(KKPER),
     *           Z(LSSUMW),IX(LSNPCL),
     *           X(LSPC),X(LSPR),X(LSPL),X(LSPCON),LAYHDT,NP,
cgzh varpt
     *           JUNIT(13),JUNIT(14),X(LSPCOR),X(LSPROR),X(LSPLOR),
cgzh these three are dummies for nptlaya etc
     *           IX(LSIGNP),IX(LSIGNP),IX(LSIGNP),
c     *           NPTLAYA,NPTROWA,NPTCOLA,
cgzh ssfix
     *           IDIM,X(LSCONC),IX(LSLMBO),NLIMBO,MULTSS,SBVL,NIUNIT,
     *           REMCRIT,GENCRIT,SRCDCY,
cgzh mnw
     &           X(LSSRCM),X(LSSNKM),IUNIT(50),
cgzh srcfix
     &           IX(LSISRC),MAXSRC,NPCELLMAX,IX(LSNPOR),IPERGWT,
     &           X(LSWTFC),
cgzh et
     *           IUNIT(5),NEVTOP,IR(LCIEVT),X(LSEVTF),ISRCFIX,
c RBW sorption     
     *           JRF)
               END IF
C
C
cgzh srcfix
                IALFLAG=1
                IF(ALLOCATED(SRCFAC)) DEALLOCATE (SRCFAC)
			  ALLOCATE(SRCFAC(MAXSRC,6,3))
C
cgzh srcfix
C DEFINE SRCFAC: PROBABILITY FUNCTION FOR CREATING PARTICLES FOR SRCFAC
cgzh srcfix2
C DEFINE SS_SRC and SS_SNK: FLOWS BETWEEN STRONG SOURCES
                CALL PTWT1SRCF(IX(LSISRC),SRCFAC,
     &            X(LSTHCK),X(LSPOR),X(LSRF),IX(LSIGNP),IG(LCIBOU),
     &            X(LSVC),X(LSVR),X(LSVL),
     &            NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,MAXSRC,
     &            RDEL,CDEL,TIMV,
cgzh srcfix2
     &            X(LS_SRC),X(LS_SNK),ISRCFIX,
     &            Z(LSTOIN),Z(LSTOOT),IX(LSCOIN),IX(LSCOOT))    
c              ELSE
c                IALFLAG=1
c                ALLOCATE(SRCFAC(1,6,2))
c                MAXSRC=1
              END IF
CMOCWT
CGWT----DETERMINE INITIAL PARTICLE WEIGHTS AT BEGINNING OF SIMULATION
cgzh debug update pt weights here (after flow) for SS period
c              IF(KKPER.EQ.1.AND.KKSTP.EQ.1) THEN
cgzh SSTR              IF(KKPER.EQ.1.AND.KKSTP.EQ.1.AND.ISSFLG(KKPER).EQ.1) THEN
              IF(KKPER.EQ.IPERGWT.AND.KKSTP.EQ.1) THEN
c	          WRITE(*,*) 'CALL PTWT1INITWT'
                IF(ISSFLG(KKPER).EQ.1)
     *		    CALL PTWT1INITWT(IG(LCIBOU),X(LSPC),X(LSPR),X(LSPL),
     *            IX(LSNPCL),Z(LSCELV),Z(LSPTWT),Z(LSSUMW),
     *            NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NPMAX,NP,IOUTS)
CMOCWT
CGWT----INITIALIZE ARRAYS THAT TRACK RESIDUAL WEIGHTS AND CONCENTRATIONS
cgzh debug new subr for dp array init
                CALL SMOC5ZD(Z(LSRESW),NSCOL,NSROW,NSLAY,0.D0)
                CALL SMOC5Z(X(LSRESC),NSCOL,NSROW,NSLAY,0.0)
              END IF
CMOCWT
CGWT----POPULATE REWETTED CELLS WITH PARTICLES
cgzh rewet
c do this for all time steps after first one [same as: .NOT.(first ts)]
              IF(IWDFLG.GT.0) THEN
			  IF(.NOT.(KKPER.EQ.IPERGWT.AND.KKSTP.EQ.1)) THEN
c for bcf, wetdry is in rx
                  IF (IUNIT(1).GT.0) THEN
	                  CALL PTWT1REWET(IG(LCIBOU),IX(LSIBOU),
     *              X(LSPC),X(LSPR),X(LSPL),X(LSPCON),IX(LSPTID),
     *              X(LSPNWC),X(LSPNWR),X(LSPNWL),NEWPTS,
     *              X(LSPCOR),X(LSPROR),X(LSPLOR),
     *              NPTCOLA,NPTROWA,NPTLAYA,IDIM,
     *              IX(LSNPCL),Z(LSCELV),Z(LSPTWT),Z(LSSUMW),
c
     *              RX(LCWETD),Z(LSSUMC),X(LSVC),X(LSVR),X(LSVL),
     *              X(LSCONC),JUNIT(13),JUNIT(14),X(LSWTFC),
     *              NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NPMAX,NP,IOUTS)
c for lpf and huf, wetdry is in x
                  ELSEIF (IUNIT(23).GT.0.OR.IUNIT(37).GT.0) THEN
	                  CALL PTWT1REWET(IG(LCIBOU),IX(LSIBOU),
     *              X(LSPC),X(LSPR),X(LSPL),X(LSPCON),IX(LSPTID),
     *              X(LSPNWC),X(LSPNWR),X(LSPNWL),NEWPTS,
     *              X(LSPCOR),X(LSPROR),X(LSPLOR),
     *              NPTCOLA,NPTROWA,NPTLAYA,IDIM,
     *              IX(LSNPCL),Z(LSCELV),Z(LSPTWT),Z(LSSUMW),
c
     *              X(LCWETD),Z(LSSUMC),X(LSVC),X(LSVR),X(LSVL),
     *              X(LSCONC),JUNIT(13),JUNIT(14),X(LSWTFC),
     *              NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NPMAX,NP,IOUTS)
                  END IF
			  END IF
              END IF
            ELSE
cgzh always allocate to default when MOCWT is off
             IF(.NOT.(ISSFLG(KKPER).EQ.1.AND.KKSTP.GT.1)) THEN
              IALFLAG=1
              IF(ALLOCATED(SRCFAC)) DEALLOCATE (SRCFAC)
              ALLOCATE(SRCFAC(1,6,2))
              MAXSRC=1
             END IF
CMOCWT  END IF FOR PTWTON=1
            END IF
C
CGWT----PERFORM EACH MOVE AND PRINT OUTPUT
		  DO 188 IMOV=1,NMOV
		    IIMOV=IMOV
	        IF(MOCTYPE.EQ.3) NRFLG=0
CMOCWT  
cgzh debug
      if(kkper.eq.1) then
	continue
	endif
              IF(PTWTON.EQ.1) 
C	CALCULATE AVERAGE SOURCE CONCENTRATION FOR WEIGHTED PARTICLE ROUTINES
C     CALCULATE MASS DECAYED IN SOURCES
c  	          WRITE(*,*) 'CALL PTWT1SRCC'
     *          CALL PTWT1SRCC(X(LSSRCF),X(LSSRCS),X(LSBSRC),X(LSBSOL),
     *           X(LSSRCV),
     *           X(LSSRCC),X(LSRF),NSCOL,NSROW,NSLAY,DECAY,TIMV,SRCDCY,
     *           X(LSDKFO),X(LSDKFS),JUNIT(11),IDKFO,IDKFS,X(LSSOL))
C
             IF(JUNIT(13).GT.0.OR.JUNIT(14).GT.0) THEN
			CALL GWT1MVOT6(X(LSSRCF),X(LSSNKF),X(LSRF),X(LSPOR),
     *          X(LSTHCK),VCMAX,VRMAX,VLMAX,TLMIN,TIMV,
     *          NTIMV,ITCD,ISSFLG(KKPER),
     *		  MAXVCJ,MAXVCI,MAXVCK,MAXVRJ,MAXVRI,MAXVRK,MAXVLJ,MAXVLI,
     *          MAXVLK,TIMDC,NTIMD,NMOV,SUMTCH,NSCOL,NSROW,NSLAY,IOUTS,
     *		  TOTIM,DELT,NODISP,DIFFUS,DECAY,ICONLY,MOCTYPE,
     *		  IIMOV,X(LSCONC),X(LSCOLD),IX(LSNPCL),IX(LSNPLD),
     *          IR(LSHFBL),IX(LSIGLK),
     *		  INTRPL,X(LSPC),X(LSPR),X(LSPL),X(LSPCON),IX(LSPTID),
     *		  X(LSVC),X(LSVR),X(LSVL),X(LSCINF),IG(LCIBOU),GX(LCDELR),
     *          GX(LCDELC),Z(LSSUMC),IX(LSLMBO),X(LSPNWC),X(LSPNWR),
     *          X(LSPNWL),IX(LSIGNP),NEWPTS,NCOL,NROW,NLAY,NPMAX,NLIMBO,
     *          NP,NCINFL,IABOVE,IBELOW,X(LSCNCN),X(LSSRCS),X(LSEVTF),
     *          IR(LCIEVT),SRCDCY,NEVTOP,KKSTP,NSTP(KKPER),KKPER,NPER,
     *          X(LSRS),
     *          FDTMTH,X(LSCAVG),X(LSDCC),X(LSDCR),X(LSDCL),X(LSDRR),
     *          X(LSDRC),X(LSDRL),X(LSDLL),X(LSDLC),X(LSDLR),X(LSVA),
     *          X(LSVAD),X(LSRHS),IX(LSCIN),IX(LSMRNO),X(LSAP),X(LSPP),
     *          X(LSRA),X(LSRR),X(LSSS),X(LSXX),X(LSWW),X(LSWORK),
     *          IX(LSCI),IX(LSCIR),IX(LSCIRH),IX(LSCIRL),NCXIT,EPSSLV,
     *          MAXIT,NBN,NRN,X(LSCHBC),X(LSCHDF),X(LSCTCF),IX(LSLBDY),
     *		  RX(LCRECH),RX(LSCRCH),IR(LCIRCH),NRCHOP,
     *		  RX(LCWELL),MXWELL,NWELLS,NWELVL,IWELLC,
     *		  RX(LCRIVR),MXRIVR,NRIVER,NRIVVL,IRIVRC,
     *		  RX(LCDRAI),MXDRN,NDRAIN,NDRNVL,
     *		  RX(LCBNDS),MXBND,NBOUND,NGHBVL,IBNDSC,
     *          RX(LCDRTF),MXDRT,NDRTCL,NDRTVL,IDRTFL,
     *		  IR(LCFLLC),RX(LCBDFV),NFLW,IFHBD4,IFHBFC,
     *		  IR(LCHDLC),RX(LCBDHV),NHED,NFHBX2,IFHBHC,
     *		  IUNIT,JRF,SBVL,NFACES,JUNIT,PERTIM,
     *		  NPNTCL,ICONFM,GZ(LCHNEW),IX(LSOBSW),
     *          NUMOBS,IOBSFL,NPNTPL,
     *          NODESS,X(LSDKZO),X(LSDKFO),X(LSDKZS),
     *		  X(LSDKFS),IDKZO,IDKFO,IDKZS,IDKFS,
     *		  AGER8,X(LSDPCON),X(LSDPRAT),X(LSDPPOR),
     *		  X(LSDPZO),X(LSDPFO),IDPZO,IDPFO,IDPPS,
     *		  X(LSRHSE),X(LSRHSO),NTFACE,X(LSDO),X(LSDIST),X(LSDIDS),
     *		  IX(LSIDB),IX(LSLB),X(LSXFOR),X(LSXBAC),X(LSYFOR),
     *          X(LSYBAC),X(LSCB),IX(LSNZIN),IX(LSNONU),X(LSSAV),
     *          IX(LSIACT),X(LSCFOR),X(LSRFOR),X(LSTFOR),X(LSCONL),
cea     *          X(LSVOL),X(LSRW),LENW,IX(LSJA),IX(LSIA),
     *          X(LSVOL),X(LSRW),X(LSRW1),LENW,IX(LSJA),IX(LSIA),
     *          X(LSA),IDIM,
     *          NELT,NODES,IX(LSJAS),IX(LSIAS),X(LSAS),NELTS,LENIW,
     *          IX(LSIW),NRFLG,NSTRMAR,RX(LCSTRM),IR(ICSTRM),NSOL,
     *          ISTCB1,
     *          ICBCFL,IR(LSGAGE),NSSAR,RX(LCOTFLW),IR(LCIVAR),IPTFLG,
     *          IR(LCOTSG),IR(ICSEG),RX(LSCOUT),RX(LSCONQ),RX(LSCNRN),
     *          RX(LSCNTB),RX(LSSLIN),RX(LSCQIN),RX(LSCGW),RX(LSLAKE),
     *          NLAKESAR,LKNODE,MXLKND,IR(ICLAKE),NTRB,NDV,IR(INTRB),
     *          IR(INDV),RX(ISTRIN),RX(ISTROT),RX(ISTGLD2),RX(ISTGNW),
     *          RX(LCCNDF),GX(LCBOTM),NBOTM,KCNT,IR(IMSUB),IR(IMSUB1),
     *          RX(LSSPPT),RX(LSSWIN),RX(LSSWOT),RX(LSCGWL),RX(LCWTDR),
     *          RX(LSCDRW),RX(LSAUG),RZ(LKACC2),RX(LKACC7),RX(LSOVOL),
     *          RX(LSPPT),RX(LSGWIN),RX(LSGWOT),RX(LSSRUN),RX(LSRNF),
     *          RX(LSSLAK),RX(LCRNF),THETA,RX(LSLKSM),IR(LSKLK),
     *          IR(LSDONE),RX(LSRTCO),RX(LSFLOB),RX(LSCLKO),RX(LSALKI),
     *          RX(LSALKO),LWRT,NUMGAGE,RX(LSCPPT),RX(LCSEG),
cgage
     *          RZ(LKACC1),RX(LKACC8),RX(LKACC9),RX(LKACC5),
     *          RX(LKACC6),RX(LKFLXI),RX(LKCNN),RX(LCSTAG),RX(LKVI),
     *          RX(LKCLKI),NIUNIT,RX(LCSFRQ),NSEGDIM,
CMOCWT
     *          Z(LSPTWT),Z(LSSUMW),X(LSSRCC),X(LSBSRC),X(LSBSNK),
     *          X(LSBSOL),Z(LSSGMS),Z(LSSGWT),
     *          Z(LSBFMS),Z(LSBFWT),
     *          Z(LSRESW),X(LSRESC),Z(LSCELV),
     *          REMCRIT,IRAND,ISEEDPT,ISEEDBD,
     *          IX(LSXDMA),NZERO,NSKIP,
cgzh varpt
     *          X(LSPCOR),X(LSPROR),X(LSPLOR),
     *          NPTCOLA,NPTROWA,NPTLAYA,
cgzh cbdy
     *          MULTSS,X(LSCINA),X(LSCINB),X(LSCINXY),
cgzh mnw
     *          nwell2,mxwel2,RZ(LCWEL2),IR(LSMNWI),
     *          X(LSSRCM),X(LSSNKM),X(LSSOLM),MNWSITE,
     *          MNWOLST,IR(LSMNWU),MNWOBS,
cgzh mnw2
     *         nmnw2,MNWMAX,RZ(LCMNW2),NODTOT,RZ(LCMNWN),WELLID,MNWPRNT,
     *         Wel1flag,QSUMflag,BYNDflag,MNWIID,RZ(LCMNIO),ntotnod,
     *         iout,RZ(LSMNWO),NMNWVL,
cgzh srcfix
     *          IX(LSISRC),MAXSRC,NPCELLMAX,SRCFAC,IX(LSNPOR),
cgzh srcfix2
     &          X(LS_SRC),X(LS_SNK),X(LS_SOL),ISRCFIX,X(LSSRCV),
     &          X(LSSOL),
     &          Z(LSTOIN),Z(LSTOOT),IX(LSCOIN),IX(LSCOOT),     
cgzh zinn%
     *          X(LSCINI),X(LSCTNM),
cgzh wtfac
     *          X(LSWTFC),JUNIT(24),			
cgzh ptob
     *          IX(LSPTOB),NUMPTOB,IX(LSPTUN),NUMPTOB_MNW,
cea
     *          SAGE,SRCAGE,
cgzh ccbd
     *          JUNIT(25),RX(LSCCBD),IPERGWT,
cgzh mbrpt
     *          JUNIT(26),JUNIT(27))
             ELSE
              CALL GWT1MVOT6(X(LSSRCF),X(LSSNKF),X(LSRF),X(LSPOR),
     *          X(LSTHCK),VCMAX,VRMAX,VLMAX,TLMIN,TIMV,
     *          NTIMV,ITCD,ISSFLG(KKPER),
     *		  MAXVCJ,MAXVCI,MAXVCK,MAXVRJ,MAXVRI,MAXVRK,MAXVLJ,MAXVLI,
     *          MAXVLK,TIMDC,NTIMD,NMOV,SUMTCH,NSCOL,NSROW,NSLAY,IOUTS,
     *		  TOTIM,DELT,NODISP,DIFFUS,DECAY,ICONLY,MOCTYPE,
     *		  IIMOV,X(LSCONC),X(LSCOLD),IX(LSNPCL),IX(LSNPLD),
     *          IR(LSHFBL),IX(LSIGLK),
     *		  INTRPL,X(LSPC),X(LSPR),X(LSPL),X(LSPCON),IX(LSPTID),
     *		  X(LSVC),X(LSVR),X(LSVL),X(LSCINF),IG(LCIBOU),GX(LCDELR),
     *          GX(LCDELC),Z(LSSUMC),IX(LSLMBO),X(LSPNWC),X(LSPNWR),
     *          X(LSPNWL),IX(LSIGNP),NEWPTS,NCOL,NROW,NLAY,NPMAX,NLIMBO,
     *          NP,NCINFL,IABOVE,IBELOW,X(LSCNCN),X(LSSRCS),X(LSEVTF),
     *          IR(LCIEVT),SRCDCY,NEVTOP,KKSTP,NSTP(KKPER),KKPER,NPER,
     *          X(LSRS),
     *          FDTMTH,X(LSCAVG),X(LSDCC),X(LSDCR),X(LSDCL),X(LSDRR),
     *          X(LSDRC),X(LSDRL),X(LSDLL),X(LSDLC),X(LSDLR),X(LSVA),
     *          X(LSVAD),X(LSRHS),IX(LSCIN),IX(LSMRNO),X(LSAP),X(LSPP),
     *          X(LSRA),X(LSRR),X(LSSS),X(LSXX),X(LSWW),X(LSWORK),
     *          IX(LSCI),IX(LSCIR),IX(LSCIRH),IX(LSCIRL),NCXIT,EPSSLV,
     *          MAXIT,NBN,NRN,X(LSCHBC),X(LSCHDF),X(LSCTCF),IX(LSLBDY),
     *		  RX(LCRECH),RX(LSCRCH),IR(LCIRCH),NRCHOP,
     *		  RX(LCWELL),MXWELL,NWELLS,NWELVL,IWELLC,
     *		  RX(LCRIVR),MXRIVR,NRIVER,NRIVVL,IRIVRC,
     *		  RX(LCDRAI),MXDRN,NDRAIN,NDRNVL,
     *		  RX(LCBNDS),MXBND,NBOUND,NGHBVL,IBNDSC,
     *          RX(LCDRTF),MXDRT,NDRTCL,NDRTVL,IDRTFL,
     *		  IR(LCFLLC),RX(LCBDFV),NFLW,IFHBD4,IFHBFC,
     *		  IR(LCHDLC),RX(LCBDHV),NHED,NFHBX2,IFHBHC,
     *		  IUNIT,JRF,SBVL,NFACES,JUNIT,PERTIM,
     *		  NPNTCL,ICONFM,GZ(LCHNEW),IX(LSOBSW),
     *          NUMOBS,IOBSFL,NPNTPL,
     *          NODESS,X(LSDKZO),X(LSDKFO),X(LSDKZS),
     *		  X(LSDKFS),IDKZO,IDKFO,IDKZS,IDKFS,
     *		  AGER8,X(LSDPCON),X(LSDPRAT),X(LSDPPOR),
     *		  X(LSDPZO),X(LSDPFO),IDPZO,IDPFO,IDPPS,
     *		  X(LSRHSE),X(LSRHSO),NTFACE,X(LSDO),X(LSDIST),X(LSDIDS),
     *		  IX(LSIDB),IX(LSLB),X(LSXFOR),X(LSXBAC),X(LSYFOR),
     *          X(LSYBAC),X(LSCB),IX(LSNZIN),IX(LSNONU),X(LSSAV),
     *          IX(LSIACT),X(LSCFOR),X(LSRFOR),X(LSTFOR),X(LSCONL),
cea     *          X(LSVOL),X(LSRW),LENW,IX(LSJA),IX(LSIA),
     *          X(LSVOL),X(LSRW),X(LSRW1),LENW,IX(LSJA),IX(LSIA),
     *          X(LSA),IDIM,
     *          NELT,NODES,IX(LSJAS),IX(LSIAS),X(LSAS),NELTS,LENIW,
     *          IX(LSIW),NRFLG,NSTRMAR,RX(LCSTRM),IR(ICSTRM),NSOL,
     *          ISTCB1,
     *          ICBCFL,IR(LSGAGE),NSSAR,RX(LCOTFLW),IR(LCIVAR),IPTFLG,
     *          IR(LCOTSG),IR(ICSEG),RX(LSCOUT),RX(LSCONQ),RX(LSCNRN),
     *          RX(LSCNTB),RX(LSSLIN),RX(LSCQIN),RX(LSCGW),RX(LSLAKE),
     *          NLAKESAR,LKNODE,MXLKND,IR(ICLAKE),NTRB,NDV,IR(INTRB),
     *          IR(INDV),RX(ISTRIN),RX(ISTROT),RX(ISTGLD2),RX(ISTGNW),
     *          RX(LCCNDF),GX(LCBOTM),NBOTM,KCNT,IR(IMSUB),IR(IMSUB1),
     *          RX(LSSPPT),RX(LSSWIN),RX(LSSWOT),RX(LSCGWL),RX(LCWTDR),
     *          RX(LSCDRW),RX(LSAUG),RZ(LKACC2),RX(LKACC7),RX(LSOVOL),
     *          RX(LSPPT),RX(LSGWIN),RX(LSGWOT),RX(LSSRUN),RX(LSRNF),
     *          RX(LSSLAK),RX(LCRNF),THETA,RX(LSLKSM),IR(LSKLK),
     *          IR(LSDONE),RX(LSRTCO),RX(LSFLOB),RX(LSCLKO),RX(LSALKI),
     *          RX(LSALKO),LWRT,NUMGAGE,RX(LSCPPT),RX(LCSEG),
cgage
     *          RZ(LKACC1),RX(LKACC8),RX(LKACC9),RX(LKACC5),
     *          RX(LKACC6),RX(LKFLXI),RX(LKCNN),RX(LCSTAG),RX(LKVI),
     *          RX(LKCLKI),NIUNIT,RX(LCSFRQ),NSEGDIM,
CMOCWT
     *          Z(LSPTWT),Z(LSSUMW),X(LSSRCC),X(LSBSRC),X(LSBSNK),
     *          X(LSBSOL),Z(LSSGMS),Z(LSSGWT),
     *          Z(LSBFMS),Z(LSBFWT),
     *          Z(LSRESW),X(LSRESC),Z(LSCELV),
     *          REMCRIT,IRAND,ISEEDPT,ISEEDBD,
     *          IX(LSXDMA),NZERO,NSKIP,
cgzh varpt
     *          X(LSPCOR),X(LSPROR),X(LSPLOR),
cgzh these three are dummies for nptlaya etc
     *          IX(LSIGNP),IX(LSIGNP),IX(LSIGNP),
c     *          NPTCOLA,NPTROWA,NPTLAYA,
cgzh cbdy
     *          MULTSS,X(LSCINA),X(LSCINB),X(LSCINXY),
cgzh mnw
     *          nwell2,mxwel2,RZ(LCWEL2),IR(LSMNWI),
     *          X(LSSRCM),X(LSSNKM),X(LSSOLM),MNWSITE,
     *          MNWOLST,IR(LSMNWU),MNWOBS,
cgzh mnw2
     *         nmnw2,MNWMAX,RZ(LCMNW2),NODTOT,RZ(LCMNWN),WELLID,MNWPRNT,
     *         Wel1flag,QSUMflag,BYNDflag,MNWIID,RZ(LCMNIO),ntotnod,
     *         iout,RZ(LSMNWO),NMNWVL,
cgzh srcfix
     *          IX(LSISRC),MAXSRC,NPCELLMAX,SRCFAC,IX(LSNPOR),
cgzh srcfix2
     &          X(LS_SRC),X(LS_SNK),X(LS_SOL),ISRCFIX,X(LSSRCV),
     &          X(LSSOL),
     &          Z(LSTOIN),Z(LSTOOT),IX(LSCOIN),IX(LSCOOT),     
cgzh zinn%
     *          X(LSCINI),X(LSCTNM),		
cgzh wtfac
     *          X(LSWTFC),JUNIT(24),			
cgzh ptob
     *          IX(LSPTOB),NUMPTOB,IX(LSPTUN),NUMPTOB_MNW,
cea
     *          SAGE,SRCAGE,
cgzh ccbd
     *          JUNIT(25),RX(LSCCBD),IPERGWT,
cgzh mbrpt
     *          JUNIT(26),JUNIT(27))
             END IF
C
CGWT----END OF MOVE LOOP
  188       CONTINUE
CGWT ellam update section for after move
cea	      IF(MOCTYPE.EQ.3.AND.ISSFLG(KKPER).EQ.0) THEN
cea			 CALL GWT1ELFLW6(GZ(LCHNEW),IG(LCIBOU),
cea     *             GX(LCHOLD),RX(LCSC1),GX(LCBOTM),NBOTM,RX(LCSC2),
cea     *             DELT,ISSFLG(KKPER),NCOL,NROW,NLAY,GX(LCBUFF),LAYHDT)
c  send old copy of IBOUND into update routine
cea			 CALL SGWT1BAS6UP(GZ(LCHNEW),
cea     *             IX(LSIBOU),GX(LCHOLD),GX(LCBOTM),NBOTM,
cea     *		     X(LCSC1),X(LCSC2),X(LSTHCK),X(LSPOR),GX(LCBUFF),
cea     *             DELT,NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,IOUTS,LAYHDT,
cea     *		     MOCTYPE,GX(LCDELR),GX(LCDELC),
cea     *             X(LSRW),X(LSRHSE),X(LSCONC))
C
cea               CALL GWT1ELLUP6(X(LSDIST),NODESS,
cea     *           NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,JRF,
cgzh debug bugfix: delr and delc in GX
cea     *           GX(LCDELR),GX(LCDELC),X(LSTHCK),IG(LCIBOU),
cea     *           X(LSCOLD),X(LSCB),X(LSRHSE),X(LSRHSO),
cea     *           NTFACE,X(LSVC),X(LSVR),X(LSVL),X(LSRF),
cea     *           X(LSPOR),TIMV,VCMAX,VRMAX,VLMAX,IOUTS,
cea     *           IX(LSLB),IX(LSNZIN),NFACES,IX(LSNONU),X(LSSAV),
cea     *           X(LSXFOR),X(LSXBAC),X(LSYFOR),X(LSYBAC),
cea     *           X(LSCINF),NCINFL,IABOVE,IBELOW,DECAY,IDIM,NRFLG,
cea     *           IX(LSIACT),X(LSCFOR),X(LSRFOR),X(LSTFOR),X(LSCONL),
cea     *           X(LSVOL),NODES,X(LSRW),IX(LSJAS),IX(LSIAS),
cea     *           X(LSAS),NELTS,LENW,LENIW,KKSTP,KKPER,SBVL,
cea     *           X(LSDIDS),X(LSCONC),IX(LSIW),NMOV,
cea     *           X(LSDO),IX(LSIDB),NIUNIT,JUNIT(9),AGER8,DELT)
cgzh debug this is not needed--lcibou never changed
c reset flow ibound pointer
c               lcibou = lsave 
cea	      END IF
C
CGWT----SKIP TO HERE IF NO SOLUTE TRANSPORT
  190 CONTINUE
C--MULTINODE WELLS OUTPUT
            IF(IUNIT(50).GT.0.AND.IGWTON.GT.0)
     &          CALL GWF1MNW1BO(MNWSITE,NWELL2,MXWEL2,
     &                          RZ(LCWEL2),IG(LCIBOU),GZ(LCHNEW),
     &                          NCOL,NROW,NODES,NSTP(KKPER),KKSTP,KKPER,
     &                          IMNWCB,ICBCFL,IOUT,IOWELL2,
     &                          TOTIM,HDRY,IUNIT(15))
C--MULTINODE WELLS OUTPUT (MNW2)
            IF(IUNIT(58).GT.0.AND.IUNIT(15).GT.0)
     &          CALL GWF1MNW2BO(IMNWCB,ICBCFL,NAUX,KKSTP,KKPER,
     &                          NCOL,NROW,NLAY,nmnw2,RZ(LCMNW2),iout,
     &                          DELT,PERTIM,TOTIM,IG(LCIBOU),
     &                          MNWMAX,NODTOT,msum,HDRY,VBNM,VBVL,
     &                          GX(LCBUFF),RZ(LCMNWN),GZ(LCHNEW),
     &                          WELLID,MNWPRNT,IGWTON,MNWOBS,IUNIT(15),
     &                          NMNWVL)
C-------Post-processing of MNW2 output --- MNWI package
        IF (IUNIT(58).GT.0) 
     &          CALL GWF1MNWIOT(Wel1flag,QSUMflag,BYNDflag,
     &                    NSTP(KKPER),kkstp,nmnw2,MNWMAX,RZ(LCMNW2),
     &                    RZ(LCMNWN),NODTOT,WELLID,TOTIM,GZ(LCHNEW),
     &                    ncol,nrow,nlay,MNWOBS,DELT,MNWIID,RZ(LCMNIO),
     &                    KKPER,ntotnod,iout,HDRY,IUNIT(15),0,1,
     &                    RZ(LSMNWO), NMNWVL)
C
C-------OBSERVATION CALCULATIONS
            IF (IPAR.NE.-3 .AND. ND.GT.0) THEN
C
C7C6----IF ITERATION FAILED TO CONVERGE THEN STOP.
              IF (ICNVG.EQ.0) CALL PLL1SD(IERR,IERRU,IOUT,IOUTG)
C-------------SHOW PROGRESS IF REQUESTED
              IF(SHOWPROG)THEN
                IF (KPER.EQ.1 .AND. KSTP.EQ.1 .AND. IPES.GT.0) THEN
                  WRITE(*,26)
                ENDIF
              ENDIF
C
C-------WRITE INITIAL PARAMETER VALUES ON FILE IOUB
              IF (IPES.GT.0)
     &          CALL PES1BAS6WB(X(LCBL),X(LCBU),IX(LCISEN),IOUB,ITERP,
     &                          ITS,IX(LCLN),NPLIST,X(LCBSCA),IOSTAR,
     &                          NPE,X(LCPARE),ITMXP,ITERPF,ITERPK)
C-------------SHOW PROGRESS IF REQUESTED
              IF(SHOWPROG)THEN
                IF (KPER.EQ.1 .AND. KSTP.EQ.1 .AND. IPES.GT.0) THEN
                  WRITE(*,'(A)')' '
                  IF (ITERPK.EQ.1) WRITE(*,'(A)')' '
                ELSEIF (KPER.EQ.NPER .AND. KSTP.EQ.NSTP(KKPER)) THEN
                  WRITE(*,26)
                ENDIF
              ENDIF
C
C-------INTERPOLATE, SAVE AND PRINT DATA FOR OBSERVATIONS.
              IF (IUNIT(28).NE.0 .AND. NH.GT.0)
     &            CALL OBS1BAS6HFM(NH,IX(LCNDER),IX(LCIOFF),IX(LCJOFF),
     &                             IX(LCMLAY),IG(LCIBOU),X(LCRINT),
     &                             OBSNAM,X(LCCOFF),X(LCROFF),
     &                             GX(LCDELR),GX(LCDELC),NCOL,NROW,NLAY,
     &                             X(LCPR),X(LCH),X(LCWT),GZ(LCHNEW),
     &                             IDRY,NPE,X(LCTOFF),MAXM,JDRY,
     &                             IPAR,IOUT,ITS,NHAR,MOBSAR,ND,IPES,
     &                             IYCFLG,GX(LCSTRT))
              IF (NQ.GT.0) THEN
                IF (IUNIT(3).NE.0)
     &              CALL OBS1DRN6FM(NQ,IX(LCNQOB),IX(LCNQCL),
     &                              IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                              GZ(LCHNEW),NCOL,NROW,NLAY,IOUT,
     &                              IG(LCIBOU),NHT,OBSNAM,X(LCH),
     &                              X(LCTOFF),MXDRN,NDRAIN,RX(LCDRAI),
     &                              X(LCWTQ),NDMH,ITS,NQAR,NQCAR,NQTAR,
     &                              NDRNVL,ND)
                IF (IUNIT(4).NE.0)
     &              CALL OBS1RIV6FM(NQ,IX(LCNQOB),IX(LCNQCL),
     &                              IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                              MXRIVR,NRIVER,RX(LCRIVR),GZ(LCHNEW),
     &                              NCOL,NROW,NLAY,IOUT,IG(LCIBOU),NH,
     &                              OBSNAM,X(LCH),X(LCTOFF),X(LCWTQ),
     &                              NDMH,ITS,NQAR,NQCAR,NQTAR,NRIVVL,ND)
                IF (IUNIT(7).NE.0)
     &              CALL OBS1GHB6FM(NQ,IX(LCNQOB),IX(LCNQCL),
     &                              IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                              MXBND,NBOUND,RX(LCBNDS),GZ(LCHNEW),
     &                              NCOL,NROW,NLAY,IOUT,IG(LCIBOU),NHT,
     &                              OBSNAM,X(LCH),X(LCTOFF),ITS,NQAR,
     &                              NQCAR,NQTAR,NGHBVL,ND,X(LCWTQ),NDMH)
                IF (IUNIT(18).NE.0)
     &              CALL OBS1STR6FM(NQ,IX(LCNQOB),IX(LCNQCL),
     &                              IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                              GZ(LCHNEW),NCOL,NROW,NLAY,IOUT,
     &                              IG(LCIBOU),NHT,OBSNAM,X(LCH),
     &                              X(LCTOFF),MXSTRM,NSTREM,RX(LCSTRM_),
     &                              IR(ICSTRM_),X(LCWTQ),NDMH,ITS,NQAR,
     &                              NQCAR,NQTAR,ND)
                IF (IUNIT(38).NE.0)
     &              CALL OBS1BAS6FFM(NQ,IX(LCNQOB),IX(LCNQCL),
     &                               IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                               GZ(LCHNEW),NCOL,NROW,NLAY,
     &                               IG(LCIBOU),NHT,X(LCH),
     &                               X(LCTOFF),ITS,NQAR,NQCAR,NQTAR,
     &                               ICHFLG,GX(LCCR),GX(LCCC),GX(LCCV),
     &                               GX(LCBOTM),NBOTM,LAYHDT,ND,IOUT,
     &                               KKPER)
                IF (IUNIT(40).NE.0)
     &              CALL OBS1DRT1FM(NQ,IX(LCNQOB),IX(LCNQCL),
     &                              IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                              GZ(LCHNEW),NCOL,NROW,NLAY,IOUT,
     &                              IG(LCIBOU),NHT,OBSNAM,X(LCH),
     &                              X(LCTOFF),MXDRT,NDRTCL,RX(LCDRTF),
     &                              X(LCWTQ),NDMH,ITS,NQAR,NQCAR,NQTAR,
     &                              NDRTVL,ND)
              ENDIF
              IF (IUNIT(29).NE.0 .AND. ISSFLG(KKPER).EQ.1)
     &            CALL OBS1ADV2P(NROW,NCOL,NLAY,GX(LCDELC),GX(LCDELR),
     &                           IOUT,GX(LCCR),GX(LCCC),GX(LCCV),
     &                           GZ(LCHNEW),IG(LCIBOU),OBSNAM,X(LCPOFF),
     &                           NHT,NQT,NTT2,NPTH,IX(LCNPNT),KTDIM,
     &                           KTFLG,KTREV,ADVSTP,
     &                           IX(LCICLS),X(LCPRST),NPRST,0,
     &                           GX(LCRMLT),X(LCHK),IG(LCIZON),X(LCH),
     &                           X(LCX),NPE,ND,X(LCTT2),IPRINT,ITERP,
     &                           IOUTT2,MXBND,NBOUND,RX(LCBNDS),NRCHOP,
     &                           IR(LCIRCH),RX(LCRECH),MXSTRM,NSTREM,
     &                           IR(ICSTRM_),RX(LCSTRM_),MXRIVR,NRIVER,
     &                           RX(LCRIVR),MXDRN,NDRAIN,RX(LCDRAI),
     &                           X(LCSV),NMLTAR,NZONAR,GX(LCBOTM),NBOTM,
     &                           RX(LCWELL),NWELVL,MXWELL,NWELLS,
     &                           Z(LCSNEW),X(LCVKA),
     &                           IUNIT(21),RX(LCHFB),MXACTFB,NHFB,
     &                           X(LCHANI),NGHBVL,NRIVVL,NDRNVL,LAYHDT,
     &                           IX(LCLN),NPLIST,ISCALS,FSNK,X(LCWTQ),
     &                           NDMH,X(LCBSCA),X(LCHKCC),X(LCHUFTHK),
     &                           NHUF,IUNIT(23),IUNIT(37),X(LCGS),
     &                           X(LCVDHT),IUNIT(47),X(LCDVDH),NPADV,
     &                           IPFLG,IADVHUF,IADVPER,KKPER,
     &                           IMPATHOUT,TDELC)
              CALL OBS1BAS6SS(IOUT,NPE,NH,OBSNAM,KKPER,KKSTP,X(LCBUF1),
     &                        X(LCX),X(LCH),X(LCWT),X(LCHOBS),IPRINT,
     &                        IFO,ITERP,IPAR,NPER,IX(LCLN),LASTX,ISCALS,
     &                        X(LCWP),MPR,X(LCPRM),RSQ,
     &                        RSQP,IPR,IX(LCNIPR),X(LCWTPS),ND,
     &                        X(LCWTQ),X(LCWTQS),IOWTQ,NDMH,NTT2,KTDIM,
     &                        IOSTAR,NPLIST,NSTP,MPRAR,IPRAR,OUTNAM,
     &                        IX(LCIPLO),EQNAM,NAMES,IX(LCIPLP),NDMHAR,
     &                        NQTDR,NQTRV,NQTGB,NQTST,NQTCH,IOWTQCH,
     &                        IOWTQDR,IOWTQRV,IOWTQGB,IOWTQST,LCOBBAS,
     &                        LCOBDRN,LCOBRIV,LCOBGHB,LCOBSTR,LCOBCHD,
     &                        LCOBADV,X(LCSSGF),X(LCSSDR),X(LCSSRV),
     &                        X(LCSSGB),X(LCSSST),X(LCSSAD),X(LCSSCH),
     &                        X(LCSSPI),X(LCSSTO),ITMXP,IOUTG,X(LCBUF2),
     &                        IPES,X(LCBPRI),X(LCBSCA),X(LCRSQA),
     &                        X(LCRSPA),LCOBDRT,X(LCSSDT),NQTDT,IOWTQDT,
     &                        NRSO,NPOST,NNEGT,NRUNS,NQTSF,IOWTQSF,
     &                        LCOBSFR,X(LCSSSF),NHT,X(LCOTIM),OBSALL)
            ENDIF
C-----------SHOW PROGRESS IF REQUESTED
            IF(SHOWPROG)THEN
              IF (KPER.EQ.NPER .AND. KSTP.EQ.NSTP(KKPER)) THEN
                WRITE(*,26)
              ENDIF
   26         FORMAT('+',77(' '))
            ENDIF
C-----------SKIP OVER SENSITIVITY LOOP?
            IF (IPAR.NE.-3) THEN
              IF (IPAR.EQ.-1) GOTO 90
              IF (IFO.EQ.1 .AND. LASTX.EQ.0) GOTO 90
            ELSE
C7C6------IF ITERATION FAILED TO CONVERGE THEN STOP.
              IF (ICNVG.EQ.0) CALL PLL1SD(IERR,IERRU,IOUT,IOUTG)
              GOTO 90
            ENDIF
C
            CALL PLL1AS(NPE)
C
C-----------LOOP THROUGH THE PARAMETERS THAT ARE TO BE ESTIMATED TO
C           CALCULATE SENSITIVITIES
C
            DO 80 KKIP = 1,NPE
              IP = KKIP
C           TO SUBROUTINE?
C             DO THOSE IP ITEMS ASSIGNED TO THIS PROCESSOR
              IF (IPDO(IP).EQ.MYID) THEN
              ELSE
                GOTO 80
              ENDIF
C-----------ASSIGN PARAMETER-APPROPRIATE CONVERGENCE CRITERIA AND OTHER
C           PARAMETER-SPECIFIC SETTINGS
              CALL SEN1BAS6CC(X(LCHCLO),X(LCRCLO),FAC,HCLOSES,IP,NPLIST,
     &                        RCLOSES,IIPP,PIDTMP,NCOL,NROW,NLAY,IUHEAD,
     &                        Z(LCSNEW),X(LCSOLD),XHS,LENXHS)
C-----------PRINT PARAMETER NAME
              CALL UMESPR('SOLVING PARAMETER SENSITIVITY FOR ',
     &                    PARNAM(IIPP),IOUT)
C-------------SHOW PROGRESS IF REQUESTED
              IF(SHOWPROG)THEN
                WRITE(*,58)KPER,KSTP,PARNAM(IIPP)
   58           FORMAT('+Solving:  Stress period: ',i5,4x,
     &                 'Time step: ',i5,4x,a,' Sensitivity')
              ENDIF
C
C7C2----ITERATIVELY FORMULATE AND SOLVE THE SENSITIVITY EQUATIONS.
C
              DO 60 KITER = 1, MXITER
                KKITER = KITER
C       PREPARE TO CALCULATE SENSITIVITY-EQUATION RHS FOR ONE PARAMETER
                CALL SEN1BAS6FM(NCOL,NLAY,NROW,GX(LCRHS))
C-------------CALCULATE MATRIX AND VECTOR DERIVATIVES, MULTIPLY BY
C-------------HEADS, AND ADD COMPONENTS TO RHS
                IF (PIDTMP.EQ.'GHB ')
     &              CALL SEN1GHB6FM(MXBND,RX(LCBNDS),GZ(LCHNEW),NCOL,
     &                              NROW,NLAY,IG(LCIBOU),GX(LCRHS),
     &                              IIPP,NGHBVL)
                IF (PIDTMP.EQ.'DRN ')
     &              CALL SEN1DRN6FM(MXDRN,RX(LCDRAI),GZ(LCHNEW),NCOL,
     &                              NROW,NLAY,IG(LCIBOU),GX(LCRHS),
     &                              IIPP,NDRNVL)
                IF (PIDTMP.EQ.'RIV ')
     &              CALL SEN1RIV6FM(MXRIVR,RX(LCRIVR),GZ(LCHNEW),NCOL,
     &                              NROW,NLAY,IG(LCIBOU),GX(LCRHS),
     &                              IIPP,NRIVVL)
                IF (PIDTMP.EQ.'STR ')
     &              CALL SEN1STR6FM(NSTREM,MXSTRM,RX(LCSTRM_),
     &                              GZ(LCHNEW),NCOL,NROW,NLAY,
     &                              IG(LCIBOU),GX(LCRHS),IR(ICSTRM_),
     &                              IIPP)
                IF (PIDTMP.EQ.'Q   ')
     &              CALL SEN1WEL6FM(NWELLS,MXWELL,RX(LCWELL),NCOL,NROW,
     &                              NLAY,IG(LCIBOU),GX(LCRHS),IIPP,
     &                              NWELVL)
                IF (PIDTMP.EQ.'HK  ' .OR. PIDTMP.EQ.'VK  ' .OR.
     &              PIDTMP.EQ.'VANI' .OR. PIDTMP.EQ.'SS  ' .OR.
     &              PIDTMP.EQ.'SY  ' .OR. PIDTMP.EQ.'VKCB' .OR.
     &              PIDTMP.EQ.'HANI' .OR. PIDTMP.EQ.'LVDA' .OR.
     &              PIDTMP.EQ.'KDEP' .OR. PIDTMP.EQ.'SYTP') THEN
                  IF (IUNIT(23).NE.0)
     &                CALL SEN1LPF1FM(GX(LCRMLT),GZ(LCHNEW),NCOL,NROW,
     &                                NLAY,ISSFLG(KKPER),PIDTMP,X(LCHK),
     &                                GX(LCDELR),GX(LCDELC),IG(LCIBOU),
     &                                DELT,GX(LCRHS),GX(LCHOLD),
     &                                IG(LCIZON),GX(LCCV),X(LCSV),
     &                                NMLTAR,NZONAR,IIPP,GX(LCBOTM),
     &                                NBOTM,X(LCVKA),IUNIT(21),
     &                                RX(LCHFB),MXACTFB,NHFB,X(LCHANI))
                  IF (IUNIT(37).NE.0) THEN
                    IF (IUNIT(47).EQ.0) THEN
                      CALL SEN1HUF2FM(GZ(LCHNEW),NCOL,NROW,NLAY,PIDTMP,
     &                                X(LCHK),X(LCHKCC),GX(LCDELR),
     &                                GX(LCDELC),IG(LCIBOU),GX(LCRHS),
     &                                GX(LCCV),GX(LCBOTM),NBOTM,
     &                                X(LCHUFTHK),NHUF,IIPP,IG(LCIZON),
     &                                NZONAR,GX(LCRMLT),NMLTAR,
     &                                IUNIT(21),RX(LCHFB),MXACTFB,NHFB,
     &                                GX(LCHOLD),DELT,ISSFLG(KKPER),
     &                                IOUT,X(LCGS))
                    ELSE
                      CALL SEN1HUF2VDFM(GZ(LCHNEW),Z(LCSNEW),IG(LCIBOU),
     &                                  X(LCHK),X(LCHKCC),GX(LCCR),
     &                                  GX(LCCC),GX(LCCV),X(LCVDHT),
     &                                  X(LCVDHD),X(LCDVDH),GX(LCRHS),
     &                                  NCOL,NROW,NLAY,GX(LCDELR),
     &                                  GX(LCDELC),X(LCHUFTHK),NHUF,
     &                                  GX(LCBOTM),NBOTM,IIPP,PIDTMP,
     &                                  IG(LCIZON),NZONAR,GX(LCRMLT),
     &                                  NMLTAR,X(LCGS))
                    ENDIF
                  ENDIF
                ENDIF
                IF (IUNIT(23).NE.0)
     &              CALL SEN1LPF1UN(ISSFLG(KKPER),DELT,NCOL,NROW,NLAY,
     &                              X(LCSOLD),GZ(LCHNEW),Z(LCSNEW),
     &                              GX(LCDELR),GX(LCDELC),IG(LCIBOU),
     &                              GX(LCRHS),X(LCSC1),GX(LCCR),
     &                              GX(LCCC),KKITER,X(LCSC2),X(LCHK),
     &                              GX(LCBOTM),NBOTM,GX(LCHOLD),
     &                              GX(LCCV),X(LCHANI),X(LCVKA))
                IF (IUNIT(37).NE.0 .AND. IUNIT(47).EQ.0)
     &              CALL SEN1HUF2UN(ISSFLG(KKPER),DELT,NCOL,NROW,NLAY,
     &                              X(LCSOLD),GZ(LCHNEW),Z(LCSNEW),
     &                              GX(LCDELR),GX(LCDELC),IG(LCIBOU),
     &                              GX(LCRHS),X(LCSC1),GX(LCCR),
     &                              GX(LCCC),KKITER,X(LCHK),X(LCHKCC),
     &                              GX(LCBOTM),NBOTM,GX(LCHOLD),
     &                              GX(LCCV),X(LCHUFTHK),NHUF,
     &                              IG(LCIZON),NZONAR,GX(LCRMLT),
     &                              NMLTAR,X(LCGS),IOUT)
                IF (PIDTMP.EQ.'HFB ')
     &              CALL SEN1HFB6FM(GX(LCBOTM),GX(LCDELC),GX(LCDELR),
     &                              GZ(LCHNEW),RX(LCHFB),IIPP,MXACTFB,
     &                              MXHFB,NBOTM,NCOL,NLAY,NROW,
     &                              GX(LCRHS),LAYHDT,NHFBNP,IG(LCIBOU))
                IF (PIDTMP.EQ.'RCH ')
     &              CALL SEN1RCH6FM(NCOL,NROW,NLAY,GX(LCDELR),
     &                              GX(LCDELC),GX(LCRMLT),NRCHOP,
     &                              IR(LCIRCH),IG(LCIBOU),GX(LCRHS),
     &                              IG(LCIZON),NMLTAR,NZONAR,IIPP)
                IF (PIDTMP.EQ.'EVT ')
     &              CALL SEN1EVT6FM(NCOL,NROW,NLAY,GX(LCDELR),
     &                              GX(LCDELC),GX(LCRMLT),NEVTOP,
     &                              IR(LCIEVT),IG(LCIBOU),GX(LCRHS),
     &                              RX(LCSURF),RX(LCEXDP),GZ(LCHNEW),
     &                              IG(LCIZON),NMLTAR,NZONAR,IIPP)
                IF (PIDTMP.EQ.'CHD ')
     &              CALL SEN1CHD6FM(MXCHD,RX(LCCHDS),Z(LCSNEW),
     &                              PERLEN(KKPER),PERTIM,NCOL,NROW,NLAY,
     &                              NCHDVL,IOUT,IIPP,IERR,IERRU)
                IF (PIDTMP.EQ.'ETS ')
     &              CALL SEN1ETS1FM(NCOL,NROW,NLAY,GX(LCDELR),
     &                              GX(LCDELC),GX(LCRMLT),NETSOP,
     &                              IR(LCIETS),IG(LCIBOU),GX(LCRHS),
     &                              RX(LCETSS),RX(LCETSX),GZ(LCHNEW),
     &                              IG(LCIZON),NMLTAR,NZONAR,IIPP,
     &                              NETSEG,RX(LCPXDP),RX(LCPETM),NSEGAR)
                IF (PIDTMP.EQ.'DRT ')
     &              CALL SEN1DRT1FM(MXDRT,RX(LCDRTF),GZ(LCHNEW),NCOL,
     &                              NROW,NLAY,IG(LCIBOU),GX(LCRHS),
     &                              IIPP,NDRTVL,IDRTFL)
                IF (IERR.GT.0) GOTO 85
C
C-------IF SNEW=SOLD=0 AND RHS=0, NO NEED TO SOLVE.
                CALL UNOITER(GX(LCRHS),Z(LCSNEW),NODES,ISA)
                IF (ISA.EQ.0) THEN
                  ICNVG = 1
                  GOTO 70
                ENDIF
C
C7C2B---MAKE ONE CUT AT AN APPROXIMATE SOLUTION.
                IF (IUNIT(9).GT.0)
     &              CALL SIP5AP(Z(LCSNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                          GX(LCCV),GX(LCHCOF),GX(LCRHS),X(LCEL),
     &                          X(LCFL),X(LCGL),X(LCV),X(LCW),X(LCHDCG),
     &                          IX(LCLRCH),NPARM,KKITER,HCLOSES,ACCL,
     &                          ICNVG,KKSTP,KKPER,IPCALC,IPRSIP,MXITER,
     &                          NSTP(KKPER),NCOL,NROW,NLAY,NODES,IOUT,3,
     &                          IERR,IERRU)
                IF (IUNIT(10).GT.0)
     &              CALL DE45AP(Z(LCSNEW),IG(LCIBOU),X(LCAU),X(LCAL),
     &                          IX(LCIUPP),IX(LCIEQP),X(LCD4B),MXUP,
     &                          MXLOW,MXEQ,MXBW,GX(LCCR),GX(LCCC),
     &                          GX(LCCV),GX(LCHCOF),GX(LCRHS),ACCL,
     &                          KKITER,ITMX,MXITER,NITER,HCLOSES,IPRD4,
     &                          ICNVG,NCOL,NROW,NLAY,IOUT,IX(LCLRCH),
     &                          X(LCHDCG),0,KKSTP,KKPER,DELT,
     &                          NSTP(KKPER),ID4DIR,ID4DIM,3,IERR,IERRU)
                IF (IUNIT(11).GT.0)
     &              CALL SOR5AP(Z(LCSNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                          GX(LCCV),GX(LCHCOF),GX(LCRHS),X(LCA),
     &                          X(LCRES),IX(LCIEQP),X(LCHDCG),
     &                          IX(LCLRCH),KKITER,HCLOSES,ACCL,ICNVG,
     &                          KKSTP,KKPER,IPRSOR,MXITER,NSTP(KKPER),
     &                          NCOL,NROW,NLAY,NSLICE,MBW,IOUT,3)
                IF (IUNIT(13).GT.0)
     &              CALL PCG2AP(Z(LCSNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                          GX(LCCV),GX(LCHCOF),GX(LCRHS),Z(LCV),
     &                          Z(LCSS),Z(LCP),X(LCCD),X(LCHCHG),
     &                          IX(LCLHCH),X(LCRCHG),IX(LCLRCH),KKITER,
     &                          NITER,HCLOSES,RCLOSES,ICNVG,KKSTP,KKPER,
     &                          IPRPCG,MXITER,ITER1,NPCOND,NBPOL,
     &                          NSTP(KKPER),NCOL,NROW,NLAY,NODES,RELAX,
     &                          IOUT,3,IX(LCIT1),DAMP,GX(LCBUFF),
     &                          X(LCHCSV),IERR,IERRU,Z(LCHPCG))
           IF (IUNIT(14).GT.0)
     &            CALL LMG1AP(Z(LCSNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                        GX(LCCV),GX(LCHCOF),GX(LCRHS),Z(LCA),
     &                        IX(LCIA),IX(LCJA),Z(LCU1),Z(LCFRHS),
     &                        IX(LCIG),ISIZ1,ISIZ2,ISIZ3,ISIZ4,KKITER,
     &                        BCLOSE,DAMP,ICNVG,KKSTP,KKPER,MXITER,
     &                        MXCYC,NCOL,NROW,NLAY,NODES,HNOFLO,IOUT,
     &                        10,ICG,IADAMP,DUP,DLOW)
c              IF (IUNIT(42).GT.0)
c     &            CALL GMG1AP(Z(LCSNEW),GX(LCRHS),GX(LCCR),GX(LCCC),
c     &                        GX(LCCV),GX(LCHCOF),HNOFLO,IG(LCIBOU),
c     &                        IITER,MXITER,RCLOSES,HCLOSES,KKITER,
c     &                        KKSTP,KKPER,ICNVG,DAMP,IADAMP,IOUTGMG,
c     &                        IOUT)
                IF (IERR.GT.0) GOTO 85
C
C7C2C---IF CONVERGENCE CRITERION HAS BEEN MET STOP ITERATING.
                IF (ICNVG.EQ.1) GOTO 70
   60         CONTINUE
              KITER = MXITER
C
C7C2C-----IF CONVERGENCE CRITERION HAS NOT BEEN MET . . .
C---------IF ESTIMATING PARAMETERS OR CALCULATING BEALE'S MEASURE, USE
C         THE AVAILABLE VALUES AND KEEP GOING
              IF (IPES.GT.0 .OR. IBEFLG.EQ.2) THEN
                ICNVG = 1
                ICNVGP = 0
              ENDIF
C---------PRINT THE DATA TABLE AND WARNING MESSAGES AND STOP EXCEPT
C         AS NOTED ABOVE
              CALL UNOCONV(X(LCBUF1+IPRAR),OBSNAM,X(LCH),
     &                     X(LCHOBS),IOUTG,IP,IPAR,IPR,KKPER,KKSTP,
     &                     IX(LCLN),MPR,ND,NDMH,NH,IX(LCNIPR),X(LCPRM),
     &                     X(LCBUF1+IPRAR+ND+MPR+IPR),RSQ,
     &                     RSQP,X(LCWP),X(LCWTPS),X(LCWT),X(LCWTQ),
     &                     X(LCWTQS),NPLIST,MPRAR,IPRAR,OUTNAM,
     &                     IX(LCIPLO),EQNAM,NAMES,IX(LCIPLP),NDMHAR,
     &                     NQTDR,NQTRV,NQTGB,NQTST,NQTCH,IOWTQCH,
     &                     IOWTQDR,IOWTQRV,IOWTQGB,IOWTQST,LCOBBAS,
     &                     LCOBDRN,LCOBRIV,LCOBGHB,LCOBSTR,LCOBCHD,
     &                     LCOBADV,X(LCSSGF),X(LCSSDR),X(LCSSRV),
     &                     X(LCSSGB),X(LCSSST),X(LCSSAD),X(LCSSCH),
     &                     X(LCSSPI),X(LCSSTO),ITMXP,IPES,X(LCBPRI),
     &                     ITERP,IERR,IERRU,NTT2,LCOBDRT,X(LCSSDT),
     &                     NQTDT,IOWTQDT,NRSO,NPOST,NNEGT,NRUNS,NQTSF,
     &                     IOWTQSF,LCOBSFR,X(LCSSSF),KTDIM,NHT,
     &                     X(LCOTIM))
              IF (IPAR.NE.1 .AND. IBEFLG.NE.2) THEN
                IERR = 1
                GOTO 85
              ENDIF
C
   70         CONTINUE
C-------------CHECK ACCURACY OF SENSITIVITY CALCULATIONS
              CALL SEN1BAS6CS(Z(LCSNEW),IG(LCIBOU),GX(LCCR),GX(LCCC),
     &                        GX(LCCV),GX(LCHCOF),GX(LCRHS),NCOL,NROW,
     &                        NLAY,IOUT,X(LCSEND),NPE,NTIMES,IP,ITS)
C
C7C5---PRINT AND/OR SAVE SENSITIVITY MATRICES.
              IF (IPAR.EQ.0 .OR. IPAR.EQ.-2) THEN
                CALL SEN1BAS6OT(IHDDFL,IOUT,ISA,KKSTP,IIPP,PIDTMP,
     &                          Z(LCSNEW),GX(LCBUFF),IR(LCIOFL),
     &                          IG(LCIBOU),KKPER,DELT,PERTIM,TOTIM,
     &                          ITMUNI,NCOL,NROW,NLAY,ICNVG,ISENFM,
     &                          ISENPU,ISENSU,CHEDFM,IXSEC,LBHDSV,
     &                          HNOFLO,IP,NPE,IPRINTS,IERR,X(LCBSCA),
     &                          NPLIST,IX(LCLN))
              ENDIF
C-------OBSERVATION CALCULATIONS
              IF (ND.GT.0) THEN
C
C7C6----IF ITERATION FAILED TO CONVERGE THEN STOP.
                IF (ICNVG.EQ.0) THEN
                  IERR = 1
                  GOTO 85
                ENDIF
C-------INTERPOLATE, SAVE AND PRINT SENSITIVITIES FOR OBSERVATIONS.
                IF (IUNIT(28).NE.0 .AND. NH.GT.0)
     &              CALL OBS1BAS6HDR(NH,IX(LCNDER),IX(LCIOFF),
     &                               IX(LCJOFF),IX(LCMLAY),X(LCRINT),
     &                               NCOL,NROW,NLAY,X(LCPR),X(LCWT),
     &                               Z(LCSNEW),X(LCX),IP,NPE,IX(LCLN),
     &                               X(LCTOFF),MAXM,IPAR,NPLIST,ITS,
     &                               NHAR,MOBSAR,ND)
                IF (NQ.GT.0) THEN
                  CALL OBS1BAS6IQ(X(LCQCLS),NQCAR)
                  IF (IUNIT(3).NE.0)
     &                CALL OBS1DRN6DR(NQ,IX(LCNQOB),IX(LCNQCL),
     &                                IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                                GZ(LCHNEW),IP,Z(LCSNEW),NCOL,NROW,
     &                                NLAY,IOUTG,IG(LCIBOU),NHT,X(LCX),
     &                                OBSNAM,NPE,IX(LCLN),X(LCTOFF),
     &                                MXDRN,NDRAIN,RX(LCDRAI),NPLIST,
     &                                ITS,NQAR,NQCAR,NQTAR,NDRNVL,IERR,
     &                                IERRU,ND)
                  IF (IUNIT(4).NE.0)
     &                CALL OBS1RIV6DR(NQ,IX(LCNQOB),IX(LCNQCL),
     &                                IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                                MXRIVR,NRIVER,RX(LCRIVR),
     &                                GZ(LCHNEW),IP,Z(LCSNEW),NCOL,NROW,
     &                                NLAY,IOUTG,IG(LCIBOU),NH,X(LCX),
     &                                OBSNAM,NPE,IX(LCLN),X(LCTOFF),
     &                                NPLIST,ITS,NQAR,NQCAR,NQTAR,
     &                                NRIVVL,IERR,IERRU,ND)
                  IF (IUNIT(7).NE.0)
     &                CALL OBS1GHB6DR(NQ,IX(LCNQOB),IX(LCNQCL),
     &                                IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                                MXBND,NBOUND,RX(LCBNDS),
     &                                GZ(LCHNEW),IP,Z(LCSNEW),NCOL,NROW,
     &                                NLAY,IOUTG,IG(LCIBOU),NHT,X(LCX),
     &                                OBSNAM,NPE,IX(LCLN),X(LCTOFF),
     &                                NPLIST,ITS,NQAR,NQCAR,NQTAR,
     &                                NGHBVL,IERR,IERRU,ND)
                  IF (IUNIT(18).NE.0)
     &                CALL OBS1STR6DR(NQ,IX(LCNQOB),IX(LCNQCL),
     &                                IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                                GZ(LCHNEW),IP,Z(LCSNEW),NCOL,NROW,
     &                                NLAY,IOUTG,IG(LCIBOU),NHT,X(LCX),
     &                                OBSNAM,NPE,IX(LCLN),X(LCTOFF),
     &                                MXSTRM,NSTREM,RX(LCSTRM_),
     &                                IR(ICSTRM_),NPLIST,ITS,NQAR,NQCAR,
     &                                NQTAR,IERR,IERRU,ND)
                  IF (IUNIT(38).NE.0)
     &                CALL OBS1BAS6FDR(NQ,IX(LCNQOB),IX(LCNQCL),
     &                                 IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                                 GZ(LCHNEW),IP,Z(LCSNEW),NCOL,
     &                                 NROW,NLAY,IG(LCIBOU),NHT,
     &                                 X(LCX),X(LCTOFF),ITS,NQAR,NQCAR,
     &                                 NQTAR,ICHFLG,GX(LCBOTM),NBOTM,
     &                                 PIDTMP,LAYHDT,GX(LCRMLT),NMLTAR,
     &                                 IG(LCIZON),NZONAR,GX(LCDELC),
     &                                 GX(LCDELR),RX(LCHFB),NHFB,
     &                                 IUNIT(21),MXACTFB,X(LCSV),
     &                                 X(LCVKA),X(LCHK),X(LCHANI),
     &                                 GX(LCCR),GX(LCCC),GX(LCCV),NPE,
     &                                 IERR,IERRU,IOUTG,IUNIT(23),
     &                                 IX(LCLN),NPLIST,ND,IUNIT(37),
     &                                 X(LCHKCC),X(LCHUFTHK),NHUF,
     &                                 X(LCGS))
                  IF (IUNIT(40).NE.0)
     &                CALL OBS1DRT1DR(NQ,IX(LCNQOB),IX(LCNQCL),
     &                                IX(LCIQOB),X(LCQCLS),IX(LCIBT),
     &                                GZ(LCHNEW),IP,Z(LCSNEW),NCOL,NROW,
     &                                NLAY,IOUTG,IG(LCIBOU),NHT,X(LCX),
     &                                OBSNAM,NPE,IX(LCLN),X(LCTOFF),
     &                                MXDRT,NDRTCL,RX(LCDRTF),NPLIST,
     &                                ITS,NQAR,NQCAR,NQTAR,NDRTVL,IERR,
     &                                IERRU,ND)
                  IF (IERR.GT.0) GOTO 85
                ENDIF
                IF (IUNIT(29).NE.0 .AND. ISSFLG(KKPER).EQ.1)
     &            CALL OBS1ADV2P(NROW,NCOL,NLAY,GX(LCDELC),GX(LCDELR),
     &                           IOUT,GX(LCCR),GX(LCCC),GX(LCCV),
     &                           GZ(LCHNEW),IG(LCIBOU),OBSNAM,X(LCPOFF),
     &                           NHT,NQT,NTT2,NPTH,IX(LCNPNT),KTDIM,
     &                           KTFLG,KTREV,ADVSTP,
     &                           IX(LCICLS),X(LCPRST),NPRST,IP,
     &                           GX(LCRMLT),X(LCHK),IG(LCIZON),X(LCH),
     &                           X(LCX),NPE,ND,X(LCTT2),IPRINT,ITERP,
     &                           IOUTT2,MXBND,NBOUND,RX(LCBNDS),NRCHOP,
     &                           IR(LCIRCH),RX(LCRECH),MXSTRM,NSTREM,
     &                           IR(ICSTRM_),RX(LCSTRM_),MXRIVR,NRIVER,
     &                           RX(LCRIVR),MXDRN,NDRAIN,RX(LCDRAI),
     &                           X(LCSV),NMLTAR,NZONAR,GX(LCBOTM),NBOTM,
     &                           RX(LCWELL),NWELVL,MXWELL,NWELLS,
     &                           Z(LCSNEW),X(LCVKA),
     &                           IUNIT(21),RX(LCHFB),MXACTFB,NHFB,
     &                           X(LCHANI),NGHBVL,NRIVVL,NDRNVL,LAYHDT,
     &                           IX(LCLN),NPLIST,ISCALS,FSNK,X(LCWTQ),
     &                           NDMH,X(LCBSCA),X(LCHKCC),X(LCHUFTHK),
     &                           NHUF,IUNIT(23),IUNIT(37),X(LCGS),
     &                           X(LCVDHT),IUNIT(47),X(LCDVDH),NPADV,
     &                           IPFLG,IADVHUF,IADVPER,KKPER,
     &                           IMPATHOUT,TDELC)
              ENDIF
C-------IF CONVERGENCE ACHIEVED BY SUM OF SQUARES CRITERIA (SOSC)
              IF (IFO.EQ.2) THEN
C-----------GO TO KSTP LOOP WHEN DONE WITH PARAMETERS
                IF (IP.EQ.NPE) THEN
                  ITERPF = ITERP
                  GOTO 90
                ENDIF
              ENDIF
C-------SAVE CURRENT SENSITIVITY ARRAY
              IF (IFO.EQ.0 .OR. LASTX.NE.0)
     &            CALL SEN1BAS6TM(NCOL,NROW,NLAY,IUHEAD,IP,GX(LCBUFF),
     &                            XHS,LENXHS,Z(LCSNEW))
C
C-------END OF SENSITIVITY LOOP
   80       CONTINUE
C
   85       CONTINUE
            CALL PLL1BR()
            CALL PLL1EH(IERR,IERRU,IOUT,IOUTG,MINERR)
C-----END OF TIME STEP (KSTP) AND STRESS PERIOD (KPER) LOOPS
   90     CONTINUE
  100   CONTINUE
C
C-------SHOW PROGRESS IF REQUESTED
        IF(SHOWPROG)THEN
          WRITE(*,57) '+'
   57     FORMAT(A,77(' '))
        ENDIF
C
        CALL PLL1BR()
        IF (ISEN.GT.0) THEN
          CALL PLL1MX(X(LCX),X(LCXND),NPE,ND)
          IF (IFO.NE.1 .OR. LASTX.GT.0)
     &        CALL SEN1BAS6PD(IOUT,NPE,NPER,NSTP,NTIMES,X(LCSEND),
     &                        X(LCSNDT))
        ENDIF
C
C       Post-processing of MNW list output --- Parses to time series for
C       individual wells
        IF (IUNIT(50).GT.0) CALL GWF1MNW1OT(MNWSITE,RZ(LCWEL2),NWELL2,
     &                                      MXWEL2,IOWELL2,MNWNAME)
C
C       PRINT DATA FOR OBSERVED HEADS AND FLOWS.
        IF (ND.GT.0 .AND. (IFO.NE.1 .OR. LASTX.NE.0))
     &      CALL OBS1BAS6OT(IOUT,IOUTG,NPE,NH,OBSNAM,X(LCBUF1),X(LCX),
     &                      X(LCH),X(LCWT),X(LCHOBS),IPRINT,IFO,ITERP,
     &                      IPAR,IX(LCLN),ISCALS,X(LCWP),MPR,X(LCPRM),
     &                      RSQ,RSQP,RSQO,RSQOO,SOSC,SOSR,IPR,
     &                      IX(LCNIPR),X(LCWTPS),ND,X(LCWTQ),
     &                      X(LCWTQS),IOWTQ,NDMH,NTT2,KTDIM,NPLIST,
     &                      MPRAR,IPRAR,OUTNAM,IX(LCIPLO),EQNAM,NAMES,
     &                      IX(LCIPLP),NDMHAR,NQTDR,NQTRV,NQTGB,NQTST,
     &                      NQTCH,IOWTQCH,IOWTQDR,IOWTQRV,IOWTQGB,
     &                      IOWTQST,LCOBBAS,LCOBDRN,LCOBRIV,LCOBGHB,
     &                      LCOBSTR,LCOBCHD,LCOBADV,X(LCSSGF),X(LCSSDR),
     &                      X(LCSSRV),X(LCSSGB),X(LCSSST),X(LCSSAD),
     &                      X(LCSSCH),X(LCSSPI),X(LCSSTO),ITMXP,
     &                      X(LCBUF2),IPES,X(LCBPRI),X(LCBSCA),LCOBDRT,
     &                      X(LCSSDT),NQTDT,IOWTQDT,NRSO,NPOST,
     &                      NNEGT,NRUNS,NQTSF,IOWTQSF,LCOBSFR,X(LCSSSF),
     &                      NHT,X(LCOTIM),OBSALL)
C       PARALLEL CONVERGENCE TEST
        CALL PLL1CV(IFO)
C-------IF CONVERGENCE ACHIEVED BY SUM OF SQUARES CRITERIA (SOSC)
        IF (IFO.EQ.2) THEN
          ITERPF = ITERP
        ENDIF
C
C-----NONLINEAR REGRESSION BY MODIFIED GAUSS-NEWTON
        IF (IPES.GT.0) THEN
          IF (IYCFLG.LT.1) THEN
C---------EXECUTE ONE GAUSS-NEWTON ITERATION
            IF (MYID.EQ.MPROC) THEN
              CALL PES1GAU1AP(X(LCX),ND,NPE,X(LCHOBS),X(LCWT),X(LCWP),
     &                        Z(LCC),Z(LCSCLE),Z(LCG),X(LCH),Z(LCDD),
     &                        DMAX,CSA,TOL,IND,IFO,AMP,AP,DMX,IOUTG,
     &                        X(LCB1),ITERP,IPRINT,IX(LCLN),MPR,
     &                        X(LCPRM),JMAX,NFIT,Z(LCR),Z(LCGD),
     &                        Z(LCU),NOPT,X(LCXD),Z(LCS),SOSR,
     &                        IX(LCNIPR),IPR,GX(LCBUFF),X(LCWTP),NHT,
     &                        X(LCWTQ),IOWTQ,NDMH,IOSTAR,NPLIST,MPRAR,
     &                        IPRAR,NDMHAR,X(LCBPRI),RMARM,IAP,
     &                        Z(LCDMXA),IX(LCNPAR),X(LCAMPA),X(LCAMCA),
     &                        X(LCAAP),ITMXP,RMAR,IX(LCIPNG),NPNG,
     &                        NPNGAR)
C---------FINAL OUTPUT:
C-----------PRINT SIMULATED EQUIVALENTS AND RESIDUALS IF LEAST-SQUARES
C           COEFFICIENT MATRIX IS SINGULAR OR IF PARAMETER ESTIMATION
C           DOES NOT CONVERGE
              IF (IND.GT.0 .OR. (IFO.EQ.0 .AND. KITP.EQ.ITMXP))
     &            CALL OBS1BAS6OH(X(LCWP),IOUT,NH,X(LCH),X(LCHOBS),
     &                            X(LCWT),OBSNAM,ND,MPR,X(LCPRM),RSQ,
     &                            RSQP,2,IX(LCLN),IPR,IX(LCNIPR),
     &                            X(LCWTPS),X(LCBUF1+IPRAR),
     &                            X(LCBUF1+IPRAR+ND+MPR+IPR),X(LCWTQ),
     &                            X(LCWTQS),NDMH,NTT2,KTDIM,NPLIST,
     &                            MPRAR,IPRAR,OUTNAM,IX(LCIPLO),EQNAM,
     &                            NAMES,IX(LCIPLP),NDMHAR,NQTDR,NQTRV,
     &                            NQTGB,NQTST,NQTCH,IOWTQCH,IOWTQDR,
     &                            IOWTQRV,IOWTQGB,IOWTQST,LCOBBAS,
     &                            LCOBDRN,LCOBRIV,LCOBGHB,LCOBSTR,
     &                            LCOBCHD,LCOBADV,0,X(LCSSGF),X(LCSSDR),
     &                            X(LCSSRV),X(LCSSGB),X(LCSSST),
     &                            X(LCSSAD),X(LCSSCH),X(LCSSPI),
     &                            X(LCSSTO),ITMXP,IPES,X(LCBPRI),
     &                            LCOBDRT,X(LCSSDT),NQTDT,IOWTQDT,
     &                            NRSO,NPOST,NNEGT,NRUNS,NQTSF,IOWTQSF,
     &                            LCOBSFR,X(LCSSSF),NHT,X(LCOTIM))
C-------------SHOW PROGRESS IF REQUESTED
              IF(SHOWPROG)THEN
                WRITE(*,'(A)') ' '
              ENDIF
            ENDIF
          ENDIF
          CALL PLL1BR()
          IF (NUMPROCS.GT.1) THEN
            CALL PLL1CV(IFO)
            CALL PLL1CV(ITERPF)
            CALL PLL1CV(IND)
            CALL PLL1BA(B,MXPAR)
          ENDIF
          IF (IYCFLG.LT.1) THEN
            IF (IFO.GT.0 .AND. IND.EQ.0 .AND. ITERPF.EQ.0) GOTO 20
C
C     IF PARAMETER ESTIMATION DOES NOT CONVERGE, PRINT
C     OBSERVATION-SENSITIVITY TABLE(S)
            IF (ND.GT.0 .AND. ITERP.EQ.ITMXP .AND. IFO.EQ.0)
     &          CALL OBS1BAS6NC(X(LCBUF1),X(LCBUF2),IOUTG,IOWTQ,
     &                          IX(LCIPLO),IPR,ISCALS,ITERP,IX(LCLN),
     &                          MPR,ND,NDMH,NDMHAR,NHT,NPE,NPLIST,
     &                          OBSNAM,OUTNAM,X(LCWT),X(LCWTQ),
     &                          X(LCWTQS),X(LCX),X(LCBSCA),OBSALL)
C
C-----------PRINT FINAL PARAMETER-ESTIMATION OUTPUT
C
C           WRITE CONTRIBUTIONS TO SSWR OF EACH OBSERVATION TYPE AND
C           PRIOR INFORMATION FOR EACH PARAMETER-ESTIMATION ITERATION
C           TO _ss FILE
            IF ((IFO.GT.0 .OR. ITERP.EQ.ITMXP) .AND. OUTNAM.NE.'NONE'
     &          .AND. MYID.EQ.MPROC) THEN
              CALL OBS1BAS6PR1(IFO,IOUTG,ITERPK,ITERSS,ITMXP,IUSS,
     &                         IX(LCNPAR),OUTNAM)
              IF (NH.GT.0) CALL OBS1BAS6HPR(ITERSS,ITMXP,IUSS,X(LCSSGF))
              IF (NQTCH.GT.0) CALL OBS1BAS6FPR(ITERSS,ITMXP,IUSS,
     &                                         X(LCSSCH))
              IF (NQTDR.GT.0) CALL OBS1DRN6PR(ITERSS,ITMXP,IUSS,
     &                                        X(LCSSDR))
              IF (NQTDT.GT.0) CALL OBS1DRT1PR(ITERSS,ITMXP,IUSS,
     &                                        X(LCSSDT))
              IF (NQTRV.GT.0) CALL OBS1RIV6PR(ITERSS,ITMXP,IUSS,
     &                                        X(LCSSRV))
              IF (NQTGB.GT.0) CALL OBS1GHB6PR(ITERSS,ITMXP,IUSS,
     &                                        X(LCSSGB))
              IF (NQTST.GT.0) CALL OBS1STR6PR(ITERSS,ITMXP,IUSS,
     &                                        X(LCSSST))
              IF (NOBADV.GT.0) CALL OBS1ADV2PR(ITERSS,ITMXP,IUSS,
     &                                         X(LCSSAD))
              IF (MPR.GT.0 .OR. IPR.GT.0)
     &            CALL PES1BAS6PR(ITERSS,ITMXP,IUSS,X(LCSSPI))
              CALL OBS1BAS6PR2(IPR,ITERSS,ITMXP,IUSS,MPR,X(LCSSTO))
            ENDIF
C
C           WRITE PARAMETER-ESTIMATION OUTPUT TO GLOBAL FILE
            CALL PES1BAS6OT(Z(LCC),X(LCWT),NPE,RSQ,IOUTG,GX(LCBUFF),ND,
     &                      IPRC,IFO,IND,Z(LCSCLE),X(LCHOBS),X(LCH),
     &                      X(LCB1),X(LCWP),ITERPF,IX(LCLN),MPR,
     &                      X(LCPRM),LPRINT,IDRY,EV,RSQP,VAR,IPR,
     &                      IX(LCNIPR),X(LCWTPS),DETWTP,X(LCBL),X(LCBU),
     &                      Z(LCEIGL),Z(LCEIGV),Z(LCEIGW),NHT,X(LCWTQ),
     &                      X(LCWTQS),DTLWTQ,IOWTQ,NDMH,NPLIST,MPRAR,
     &                      IPRAR,IOUB,IX(LCISEN),IBEALE,ITERP,ITMXP,
     &                      NDMHAR,X(LCPRNT),OUTNAM,X(LCPARE),X(LCSSPI),
     &                      X(LCSSTO),IX(LCNPAR),Z(LCDMXA),X(LCBPRI),
     &                      X(LCBSCA),IPRINT,X(LCAAP),X(LCAMCA),
     &                      X(LCRSQA),X(LCRSPA),X(LCAMPA),ITERPK,OBSALL,
     &                      IUSS)
            IF (IFO.EQ.0 .AND. ITERP.EQ.ITMXP) GOTO 110
C
          ENDIF
        ENDIF
C-------GENERATE INPUT FILE(S) FOR RESAN-2000, BEALE-2000 AND YCINT-2000
        IF (MYID.EQ.MPROC) THEN
          IF (IYCFLG.LT.1 .AND. IPES.GT.0)
     &        CALL PES1BAS6RS(NPE,ND,NDMH,VAR,Z(LCC),X(LCWT),NHT,
     &                        X(LCWTQS),X(LCX),MPR,X(LCPRM),X(LCWP),
     &                        NPLIST,MPRAR,NDMHAR,OUTNAM,X(LCWTPS),
     &                        IPR,IPRAR,IX(LCNIPR),RSQP,IDRY)
          IF (IBEFLG.GT.0 .AND. (IPES.LE.0 .OR. (IPES.GT.0 .AND.
     &        IFO.GT.0)))
     &        CALL PES1BAS6BE(NPE,ND,MPR,VAR,X(LCH),X(LCWT),X(LCX),
     &                        X(LCWP),IX(LCLN),X(LCPRM),X(LCHOBS),
     &                        Z(LCC),IBEALE,ITERPK,IOUT,OBSNAM,
     &                        GX(LCBUFF),NHT,NDMH,X(LCWTQ),NPLIST,MPRAR,
     &                        IBEFLG,OUTNAM,IUBE,BEFIRST,FSTAT,IERR,
     &                        IERRU,NDMHAR,X(LCWTP),IPR,IPRAR,X(LCBPRI),
     &                        IX(LCNIPR))
          IF (IERR.GT.0) GOTO 103
          IF (IYCFLG.GT.-1)
     &        CALL PES1BAS6YC(NPE,ND,MPR,X(LCH),X(LCWT),X(LCX),Z(LCC),
     &                        IOUT,OBSNAM,NHT,NDMH,X(LCWTQ),OUTNAM,
     &                        IYCFLG,IPR,IX(LCIPLO),IERR,IERRU,NDMHAR)
          IF (IERR.GT.0) GOTO 103
        ENDIF
  103   CONTINUE
        CALL PLL1BR()
        CALL PLL1EH(IERR,IERRU,IOUT,IOUTG,MINERR)
        IF (IBEFLG.EQ.2 .AND. IBEALE.NE.0) GOTO 20
        IF (ITERPF.GT.0) GOTO  107
C
C     END OF PARAMETER-ESTIMATION LOOP
  105 CONTINUE
C
  107 CONTINUE
C-------RESIDUAL ANALYSIS
C        OBS1BAS6RE CHANGES H AND MAY CHANGE HOBS
      IF (MYID.EQ.MPROC) THEN
        IF (ND.GT.0)
     &      CALL OBS1BAS6RE(X(LCWP),IOUTG,IOUT,NHT,X(LCH),X(LCHOBS),
     &                      X(LCWT),NDMH,ND,IPAR,MPR,X(LCPRM),IPR,
     &                      IX(LCNIPR),X(LCWTPS),X(LCBUF1),LBUFF,
     &                      X(LCWTQ),X(LCWTQS),NPLIST,MPRAR,IPRAR,
     &                      NDMHAR,NAMES,IX(LCOBSE),X(LCBPRI),RSQP,
     &                      NRSO,NPOST,NNEGT,NRUNS)
C
C       PRINT FINAL PARAMETER-ESTIMATION OUTPUT
        IF (IPES.GT.0 .AND. IYCFLG.LT.1)
     &      CALL PES1BAS6FO(ICNVGP,IFO,IOUTG)
      ENDIF
C
  110 CONTINUE
C     WRITE ANY RECORDS TO BE USED IN RESTARTING FUTURE SIMULATIONS
C       SAVE RESTART RECORDS FOR SUB PACKAGE
      IF(IUNIT(54).GT.0) CALL GWF1SUB1SV(ND2,IDSAVE)
C     DEALLOCATE VBAL
      IF (JUNIT(28).GT.0) CALL GWT1VBAL1DA()      
C8------END OF SIMULATION
CGWT----FOR GWT, ADD EXTRA SPACE FOR SCREEN OUTPUT
	IF(IGWTON.GT.0) WRITE (*,'(//)') 
C
      CALL GLO1BAS6ET(IBDT,IOUTG,IPRTIM,NOTICECOUNT)
      CALL CLOSEFILES(INUNIT,FNAME)
      IF (IBATCH.GT.0) THEN
C       TO USE STATIC MEMORY ALLOCATION, COMMENT OUT THE FOLLOWING
C       DEALLOCATE STATEMENTS
        DEALLOCATE (GX,GZ,IG,X,Z,IX,XHS,RX,IR,NIPRNAM,EQNAM,NAMES,
     &              OBSNAM,RZ)
        IF(IUNIT(50).GT.0) DEALLOCATE(MNWSITE)
C       DEALLOCATE (DA) PROCEDURE
        IF(IUNIT(54).GT.0) CALL GWF1SUB1DA()
c        IF(IUNIT(42).GT.0) CALL GMG1DA()
CGWT       
        IF(IGWTON.GT.0) THEN
          IF(ALLOCATED(SRCFAC)) DEALLOCATE (SRCFAC)
          DEALLOCATE(MNWOLST)
          DEALLOCATE(PTOBLST)
cgzh moved allocation from below to fix batch-mode bug
          IF(JUNIT(13).GT.0.OR.JUNIT(14).GT.0)
     *      DEALLOCATE(NPTLAYA,NPTROWA,NPTCOLA)
	  END IF
        GOTO 10
      ENDIF
C
C     HANDLE WARNINGS AND ERRORS
      CALL PLL1BR()
      CALL PLL1EH(IERR,IERRU,IOUT,IOUTG,MINERR)
      IF (MINERR.LT.0) CALL PLL1SD(IERR,IERRU,IOUT,IOUTG)
      CALL PLL1DE(IERRU,IOUT,IOUTG)
C
  120 CONTINUE
C
      CALL PLL1CL()
      WRITE(*,121)
121   FORMAT(1X,'Normal termination of MODFLOW-2000')
      CALL USTOP(' ')
cgzh debug
      read(*,*)
cgzh varpt this moved above to fix batch-mode bug
c      IF(IUNIT(15).GT.0.AND.(JUNIT(13).GT.0.OR.JUNIT(14).GT.0)) 
c     *  DEALLOCATE(NPTLAYA,NPTROWA,NPTCOLA)
      STOP 
C
      END

