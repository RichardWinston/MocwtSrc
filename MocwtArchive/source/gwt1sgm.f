*DECK SGMRES
      SUBROUTINE SGMRES(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, 
     $     RGWK, LRGW, IGWK, LIGW, RWORK, IWORK )
C***BEGIN PROLOGUE  SGMRES
C***DATE WRITTEN   871001   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SGMRES-S),
C             Non-Symmetric Linear system, Sparse,
C             Iterative Precondition, Generalized Minimum Residual
C***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov
C           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C***PURPOSE  Preconditioned GMRES iterative sparse Ax=b solver.
C            This routine uses the generalized minimum residual 
C            (GMRES) method with preconditioning to solve
C            non-symmetric linear systems of the form: A*x = b.
C***DESCRIPTION
C *Usage:
C      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C      INTEGER   IERR, IUNIT, LRGW, LIGW, IGWK(LIGW)
C      INTEGER   IWORK(USER DEFINED)
C      REAL      B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N), 
C      REAL      RGWK(LRGW), RWORK(USER DEFINED)
C      EXTERNAL  MATVEC, MSOLVE
C
C      CALL SGMRES(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
C     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, 
C     $     RGWK, LRGW, IGWK, LIGW, RWORK, IWORK)
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand side vector.
C X      :INOUT    Real X(N).
C         On input X is your initial guess for the solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Integer A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See "Description", below
C         for more late breaking details...
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MATVEC :EXT      External.
C         Name of a routine which performs the matrix vector multiply
C         Y = A*X given A and X.  The name of the MATVEC routine must 
C         be declared external in the calling program.  The calling 
C         sequence to MATVEC is:
C             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )
C         where N is the number of unknowns, Y is the product A*X
C         upon return, X is an input vector, and NELT is the number of
C         non-zeros in the SLAP IA, JA, A storage for the matrix A.  
C         ISYM is a flag which, if non-zero, denotes that A is 
C         symmetric and only the lower or upper triangle is stored.
C MSOLVE :EXT      External.
C         Name of the routine which solves a linear system Mz = r for
C         z given r with the preconditioning matrix M (M is supplied via
C         RWORK and IWORK arrays.  The name of the MSOLVE routine must 
C         be declared external in the calling program.  The calling 
C         sequence to MSLOVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         Where N is the number of unknowns, R is the right-hand side 
C         vector, and z is the solution upon return.  RWORK is a real 
C         array that can be used to pass necessary preconditioning 
C         information and/or workspace to MSOLVE.  IWORK is an integer 
C         work array for the same purpose as RWORK.
C ITOL   :IN       Integer.
C         Flag to indicate the type of convergence criterion used.
C         ITOL=0  Means the  iteration stops when the test described
C                 below on  the  residual RL  is satisfied.  This is
C                 the  "Natural Stopping Criteria" for this routine.
C                 Other values  of   ITOL  cause  extra,   otherwise
C                 unnecessary, computation per iteration and     are
C                 therefore  much less  efficient.  See  ISSGMR (the
C                 stop test routine) for more information.
C         ITOL=1  Means   the  iteration stops   when the first test
C                 described below on  the residual RL  is satisfied,
C                 and there  is either right  or  no preconditioning
C                 being used.
C         ITOL=2  Implies     that   the  user    is   using    left
C                 preconditioning, and the second stopping criterion
C                 below is used.
C         ITOL=3  Means the  iteration stops   when  the  third test
C                 described below on Minv*Residual is satisfied, and
C                 there is either left  or no  preconditioning begin
C                 used.
C         ITOL=11 is    often  useful  for   checking  and comparing
C                 different routines.  For this case, the  user must
C                 supply  the  "exact" solution or  a  very accurate
C                 approximation (one with  an  error much less  than
C                 TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C                 if ITOL=11, iteration stops when the 2-norm of the
C                 difference between the iterative approximation and
C                 the user-supplied solution  divided by the  2-norm
C                 of the  user-supplied solution  is  less than TOL.
C                 Note that this requires  the  user to  set up  the
C                 "COMMON     /SOLBLK/ SOLN(LENGTH)"  in the calling
C                 routine.  The routine with this declaration should
C                 be loaded before the stop test so that the correct
C                 length is used by  the loader.  This procedure  is
C                 not standard Fortran and may not work correctly on
C                 your   system (although  it  has  worked  on every
C                 system the authors have tried).  If ITOL is not 11
C                 then this common block is indeed standard Fortran.
C TOL    :INOUT    Real.
C         Convergence criterion, as described below.  If TOL is set
C         to zero on input, then a default value of 500*(the smallest
C         positive magnitude, machine epsilon) is used.
C ITMAX  :DUMMY    Integer.
C         Maximum number of iterations in most SLAP routines.  In
C         this routine this does not make sense.  The maximum number
C         of iterations here is given by ITMAX = MAXL*(NRMAX+1).
C         See IGWK for definitions of MAXL and NRMAX.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.  Letting norm() denote the Euclidean
C         norm, ERR is defined as follows..
C
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               for right or no preconditioning, and
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               for left preconditioning.
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               since right or no preconditioning
C                               being used.
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               since left preconditioning is being
C                               used.
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
C                               i=1,n
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
C IERR   :OUT      Integer.
C         Return error flag.
C               IERR = 0 => All went well.
C               IERR = 1 => Insufficient storage allocated for 
C                           RGWK or IGWK.
C               IERR = 2 => Routine SGMREs failed to reduce the norm
C                           of the current residual on its last call,
C                           and so the iteration has stalled.  In
C                           this case, X equals the last computed
C                           approximation.  The user must either
C                           increase MAXL, or choose a different
C                           initial guess.
C               IERR =-1 => Insufficient length for RGWK array.
C                           IGWK(6) contains the required minimum
C                           length of the RGWK array.
C               IERR =-2 => Inconsistent ITOL and JPRE values.
C         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
C         left-hand-side of the relevant stopping test defined
C         below associated with the residual for the current
C         approximation X(L).
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C SB     :IN       Real SB(N).
C         Array of length N containing scale factors for the right 
C         hand side vector B.  If JSCAL.eq.0 (see below), SB need 
C         not be supplied.
C SX     :IN       Real SX(N).
C         Array of length N containing scale factors for the solution
C         vector X.  If JSCAL.eq.0 (see below), SX need not be
C         supplied.  SB and SX can be the same array in the calling
C         program if desired.
C RGWK   :INOUT    Real RGWK(LRGW).
C         Real array of size at least 1 + N*(MAXL+6) + MAXL*(MAXL+3)
C         used for work space by SGMRES.  See below for definition of
C         MAXL.
C         On return, RGWK(1) = RHOL.  See IERR for definition of RHOL.
C LRGW   :IN       Integer.
C         Length of the real workspace, RGWK.  LRGW > 1 + N*(MAXL+6) 
C         + MAXL*(MAXL+3). For the default values, RGWK has size at 
C         least 131 + 16*N.
C IGWK   :INOUT    Integer IGWK(LIGW).  
C         The following IGWK parameters should be set by the user
C         before calling this routine.
C         IGWK(1) = MAXL.  Maximum dimension of Krylov subspace in 
C            which X - X0 is to be found (where, X0 is the initial 
C            guess).  The default value of MAXL is 10.
C         IGWK(2) = KMP.  Maximum number of previous Krylov basis 
C            vectors to which each new basis vector is made orthogonal.
C            The default value of KMP is MAXL.
C         IGWK(3) = JSCAL.  Flag indicating whether the scaling 
C            arrays SB and SX are to be used.
C            JSCAL = 0 => SB and SX are not used and the algorithm 
C               will perform as if all SB(I) = 1 and SX(I) = 1.
C            JSCAL = 1 =>  Only SX is used, and the algorithm
C               performs as if all SB(I) = 1.
C            JSCAL = 2 =>  Only SB is used, and the algorithm
C               performs as if all SX(I) = 1.
C            JSCAL = 3 =>  Both SB and SX are used.
C         IGWK(4) = JPRE.  Flag indicating whether preconditioning 
C            is being used.
C            JPRE = 0  =>  There is no preconditioning.
C            JPRE > 0  =>  There is preconditioning on the right
C               only, and the solver will call routine MSOLVE.
C            JPRE < 0  =>  There is preconditioning on the left
C               only, and the solver will call routine MSOLVE.
C         IGWK(5) = NRMAX.  Maximum number of restarts of the 
C            Krylov iteration.  The default value of NRMAX = 10.
C            if IWORK(5) = -1,  then no restarts are performed (in 
C            this case, NRMAX is set to zero internally).
C         The following IWORK parameters are diagnostic information
C         made available to the user after this routine completes.
C         IGWK(6) = MLWK.  Required minimum length of RGWK array.
C         IGWK(7) = NMS.  The total number of calls to MSOLVE. 
C LIGW   :IN       Integer.
C         Length of the integer workspace, IGWK.  LIGW >= 20. 
C
C *Description:
C       SGMRES solves a linear system A*X = B rewritten in the form:
C
C        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,
C
C       with right preconditioning, or
C
C        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,
C
C       with left preconditioning, where A is an N-by-N real matrix,
C       X  and  B are N-vectors,   SB and SX   are  diagonal scaling
C       matrices,   and M is  a preconditioning    matrix.   It uses
C       preconditioned  Krylov   subpace  methods  based     on  the
C       generalized minimum residual  method (GMRES).   This routine
C       optionally performs  either  the  full     orthogonalization
C       version of the  GMRES  algorithm or an incomplete variant of
C       it.  Both versions use restarting of the linear iteration by
C       default, although the user can disable this feature.
C       
C       The GMRES  algorithm generates a sequence  of approximations
C       X(L) to the  true solution of the above  linear system.  The
C       convergence criteria for stopping the  iteration is based on
C       the size  of the  scaled norm of  the residual  R(L)  =  B -
C       A*X(L).  The actual stopping test is either:
C
C               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B),
C
C       for right preconditioning, or
C
C               norm(SB*(M-inverse)*(B-A*X(L))) .le. 
C                       TOL*norm(SB*(M-inverse)*B),
C
C       for left preconditioning, where norm() denotes the euclidean
C       norm, and TOL is  a positive scalar less  than one  input by
C       the user.  If TOL equals zero  when SGMRES is called, then a
C       default  value  of 500*(the   smallest  positive  magnitude,
C       machine epsilon) is used.  If the  scaling arrays SB  and SX
C       are used, then  ideally they  should be chosen  so  that the
C       vectors SX*X(or SX*M*X) and  SB*B have all their  components
C       approximately equal  to  one in  magnitude.  If one wants to
C       use the same scaling in X  and B, then  SB and SX can be the
C       same array in the calling program.
C
C       The following is a list of the other routines and their
C       functions used by SGMRES:
C       SPIGMR  Contains the main iteration loop for GMRES.
C       SORTH   Orthogonalizes a new vector against older basis vects.
C       SHEQR   Computes a QR decomposition of a Hessenberg matrix.
C       SHELS   Solves a Hessenberg least-squares system, using QR 
C               factors.
C       SRLCAL  Computes the scaled residual RL.
C       SXLCAL  Computes the solution XL.
C       ISSGMR  User-replaceable stopping routine.
C
C       This routine does  not care  what matrix data   structure is
C       used for  A and M.  It simply   calls  the MATVEC and MSOLVE
C       routines, with  the arguments as  described above.  The user
C       could write any type of structure and the appropriate MATVEC
C       and MSOLVE routines.  It is assumed  that A is stored in the
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
C       stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP
C       routines SSDCG and SSICCG are examples of this procedure.
C       
C       Two  examples  of  matrix  data structures  are the: 1) SLAP
C       Triad  format and 2) SLAP Column format.
C       
C       =================== S L A P Triad format ===================
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C *Precision:           Single Precision
C***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh,
C                 "Reduced Storage Matrix Methods In Stiff ODE 
C                 Systems," LLNL report UCRL-95088, Rev. 1, 
C                 June 1987.
C***ROUTINES CALLED  SPIGMR, SORTHO, SHEQR, SHELS, SRCAL, SXLCAL, 
C                    ISSGMR, SNRM2, SDOT, SAXPY, SSCAL, ISAMAX, R1MACH.
C***END PROLOGUE  SGMRES
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
      INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      INTEGER  IERR, IUNIT, LRGW, LIGW, IGWK(LIGW)
      INTEGER  IWORK(*)
      REAL     B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N), RGWK(LRGW)
      REAL     RWORK(*)
      EXTERNAL MATVEC, MSOLVE, R1MACH
      INTEGER JPRE, KMP, MAXL, NMS, MAXLP1, NMSL, NRSTS, NRMAX
      INTEGER I, IFLAG, LR, LDL, LHES, LGMR, LQ, LV, LW
      REAL BNRM, RHOL, SUM
C
C***FIRST EXECUTABLE STATEMENT  SGMRES
      IERR = 0
C   ------------------------------------------------------------------
C         Load method parameters with user values or defaults.
C   ------------------------------------------------------------------
      MAXL = IGWK(1)
      IF (MAXL .EQ. 0) MAXL = 10
      IF (MAXL .GT. N) MAXL = N
      KMP = IGWK(2)
      IF (KMP .EQ. 0) KMP = MAXL
      IF (KMP .GT. MAXL) KMP = MAXL
      JSCAL = IGWK(3)
      JPRE = IGWK(4)
C         Check for consistent values of ITOL and JPRE.
      IF( ITOL.EQ.1 .AND. JPRE.LT.0 ) GOTO 650
      IF( ITOL.EQ.2 .AND. JPRE.GE.0 ) GOTO 650
      NRMAX = IGWK(5)
      IF( NRMAX.EQ.0 ) NRMAX = 10
C         If NRMAX .eq. -1, then set NRMAX = 0 to turn off restarting.
      IF( NRMAX.EQ.-1 ) NRMAX = 0
C         If input value of TOL is zero, set it to its default value.
      IF( TOL.EQ.0.0E0 ) TOL = 500.0*R1MACH(3)
C
C         Initialize counters.
      ITER = 0
      NMS = 0
      NRSTS = 0
C   ------------------------------------------------------------------
C         Form work array segment pointers.
C   ------------------------------------------------------------------
      MAXLP1 = MAXL + 1
      LV = 1
      LR = LV + N*MAXLP1
      LHES = LR + N + 1
      LQ = LHES + MAXL*MAXLP1
      LDL = LQ + 2*MAXL
      LW = LDL + N
      LXL = LW + N
      LZ = LXL + N
C
C         Load igwk(6) with required minimum length of the rgwk array.
      IGWK(6) = LZ + N - 1
      IF( LZ+N-1.GT.LRGW ) GOTO 640
C   ------------------------------------------------------------------
C         Calculate scaled-preconditioned norm of RHS vector b.
C   ------------------------------------------------------------------
      IF (JPRE .LT. 0) THEN
         CALL MSOLVE(N, B, RGWK(LR), NELT, IA, JA, A, ISYM,
     $        RWORK, IWORK)
         NMS = NMS + 1
      ELSE
         CALL SCOPY(N, B, 1, RGWK(LR), 1)
      ENDIF
      IF( JSCAL.EQ.2 .OR. JSCAL.EQ.3 ) THEN
         SUM = 0.E0
         DO 10 I = 1,N
            SUM = SUM + (RGWK(LR-1+I)*SB(I))**2
 10      CONTINUE
         BNRM = SQRT(SUM)
      ELSE
         BNRM = SNRM2(N,RGWK(LR),1)
      ENDIF
      if (bnrm.eq.0.) then
        do 45 i=1,n
           x(i)=0.
45      continue
        goto 600
      endif
C   ------------------------------------------------------------------
C         Calculate initial residual.
C   ------------------------------------------------------------------
      CALL MATVEC(N, X, RGWK(LR), NELT, IA, JA, A, ISYM)
      DO 50 I = 1,N
         RGWK(LR-1+I) = B(I) - RGWK(LR-1+I)
 50   CONTINUE
C   ------------------------------------------------------------------
C         If performing restarting, then load the residual into the 
C         correct location in the Rgwk array.
C   ------------------------------------------------------------------
 100  CONTINUE
      IF( NRSTS.GT.NRMAX ) GOTO 610
      IF( NRSTS.GT.0 ) THEN
C         Copy the curr residual to different loc in the Rgwk array.
         CALL SCOPY(N, RGWK(LDL), 1, RGWK(LR), 1)
      ENDIF
C   ------------------------------------------------------------------
C         Use the SPIGMR algorithm to solve the linear system A*Z = R.
C   ------------------------------------------------------------------
      CALL SPIGMR(N, RGWK(LR), SB, SX, JSCAL, MAXL, MAXLP1, KMP,
     $       NRSTS, JPRE, MATVEC, MSOLVE, NMSL, RGWK(LZ), RGWK(LV),
     $       RGWK(LHES), RGWK(LQ), LGMR, RWORK, IWORK, RGWK(LW),
     $       RGWK(LDL), RHOL, NRMAX, B, BNRM, X, RGWK(LXL), ITOL,
     $       TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
      ITER = ITER + LGMR
      NMS = NMS + NMSL
C
C         Increment X by the current approximate solution Z of A*Z = R.
C
      LZM1 = LZ - 1
      DO 110 I = 1,N
         X(I) = X(I) + RGWK(LZM1+I)
 110  CONTINUE
      IF( IFLAG.EQ.0 ) GOTO 600
      IF( IFLAG.EQ.1 ) THEN
         NRSTS = NRSTS + 1
         GOTO 100
      ENDIF
      IF( IFLAG.EQ.2 ) GOTO 620
C   ------------------------------------------------------------------
C         All returns are made through this section.
C   ------------------------------------------------------------------
C         The iteration has converged.
C
 600  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 0
      RETURN
C
C         Max number((NRMAX+1)*MAXL) of linear iterations performed.
 610  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 1
      RETURN
C
C         GMRES failed to reduce last residual in MAXL iterations.
C         The iteration has stalled.
 620  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 2
      RETURN
C         Error return.  Insufficient length for Rgwk array.
 640  CONTINUE
      ERR = TOL
      IERR = -1
      RETURN
C         Error return.  Inconsistent ITOL and JPRE values.
 650  CONTINUE
      ERR = TOL
      IERR = -2
      RETURN
C------------- LAST LINE OF SGMRES FOLLOWS ----------------------------
      END
*DECK SSDGMR
      SUBROUTINE SSDGMR(N, B, X, NELT, IA, JA, A, ISYM, NSAVE,
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $     IWORK, LENIW )
C***BEGIN PROLOGUE  SSDGMR
C***DATE WRITTEN   880615   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSDGMR-S),
C             Non-Symmetric Linear system, Sparse,
C             Iterative Precondition, Generalized Minimum Residual
C***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov
C           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C***PURPOSE  Diagonally scaled GMRES iterative sparse Ax=b solver.
C            This routine uses the generalized minimum residual 
C            (GMRES) method with diagonal scaling to solve possibly 
C            non-symmetric linear systems of the form: A*x = b.
C***DESCRIPTION
C *Usage:
C      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE
C      INTEGER   ITOL, ITMAX, IERR, IUNIT, LENW, IWORK(LENIW), LENIW
C      REAL      B(N), X(N), A(NELT), TOL, ERR
C      REAL      RWORK(LENW)
C      EXTERNAL  MATVEC, MSOLVE
C
C      CALL SSDGMR(N, B, X, NELT, IA, JA, A, ISYM, NSAVE,
C     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
C     $     RWORK, LENW, IWORK, LENIW)
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand side vector.
C X      :INOUT    Real X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Integer A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "Description", 
C         below.  If the SLAP Triad format is chosen it is changed 
C         internally to the SLAP Column format.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C NSAVE  :IN       Integer.
C         Number of direction vectors to save and orthogonalize against.
C         Must be greater than 1.
C ITOL   :IN       Integer.
C         Flag to indicate the type of convergence criterion used.
C         ITOL=0  Means the  iteration stops when the test described
C                 below on  the  residual RL  is satisfied.  This is
C                 the  "Natural Stopping Criteria" for this routine.
C                 Other values  of   ITOL  cause  extra,   otherwise
C                 unnecessary, computation per iteration and     are
C                 therefore  much less  efficient.  See  ISSGMR (the
C                 stop test routine) for more information.
C         ITOL=1  Means   the  iteration stops   when the first test
C                 described below on  the residual RL  is satisfied,
C                 and there  is either right  or  no preconditioning
C                 being used.
C         ITOL=2  Implies     that   the  user    is   using    left
C                 preconditioning, and the second stopping criterion
C                 below is used.
C         ITOL=3  Means the  iteration stops   when  the  third test
C                 described below on Minv*Residual is satisfied, and
C                 there is either left  or no  preconditioning begin
C                 used.
C         ITOL=11 is    often  useful  for   checking  and comparing
C                 different routines.  For this case, the  user must
C                 supply  the  "exact" solution or  a  very accurate
C                 approximation (one with  an  error much less  than
C                 TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C                 if ITOL=11, iteration stops when the 2-norm of the
C                 difference between the iterative approximation and
C                 the user-supplied solution  divided by the  2-norm
C                 of the  user-supplied solution  is  less than TOL.
C                 Note that this requires  the  user to  set up  the
C                 "COMMON     /SOLBLK/ SOLN(LENGTH)"  in the calling
C                 routine.  The routine with this declaration should
C                 be loaded before the stop test so that the correct
C                 length is used by  the loader.  This procedure  is
C                 not standard Fortran and may not work correctly on
C                 your   system (although  it  has  worked  on every
C                 system the authors have tried).  If ITOL is not 11
C                 then this common block is indeed standard Fortran.
C TOL    :INOUT    Real.
C         Convergence criterion, as described below.  If TOL is set
C         to zero on input, then a default value of 500*(the smallest
C         positive magnitude, machine epsilon) is used.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.  This routine uses the default
C         of NRMAX = ITMAX/NSAVE to determine the when each restart 
C         oshould ccur.  See the description of NRMAX and MAXL in 
C         SGMRES for a full and frightfully interesting discussion of 
C         this topic.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.  Letting norm() denote the Euclidean
C         norm, ERR is defined as follows...
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               for right or no preconditioning, and
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               for left preconditioning.
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               since right or no preconditioning
C                               being used.
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               since left preconditioning is being
C                               used.
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
C                               i=1,n
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
C IERR   :OUT      Integer.
C         Return error flag.
C               IERR = 0 => All went well.
C               IERR = 1 => Insufficient storage allocated for 
C                           RGWK or IGWK.
C               IERR = 2 => Routine SPIGMR failed to reduce the norm
C                           of the current residual on its last call,
C                           and so the iteration has stalled.  In
C                           this case, X equals the last computed
C                           approximation.  The user must either
C                           increase MAXL, or choose a different
C                           initial guess.
C               IERR =-1 => Insufficient length for RGWK array.
C                           IGWK(6) contains the required minimum
C                           length of the RGWK array.
C               IERR =-2 => Inconsistent ITOL and JPRE values.
C         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
C         left-hand-side of the relevant stopping test defined
C         below associated with the residual for the current
C         approximation X(L).
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C RWORK  :WORK    Real RWORK(LENW).
C         Real array of size LENW.
C LENW   :IN       Integer.
C         Length of the real workspace, RWORK.  LENW >= 1 + N*(NSAVE+7) 
C         + NSAVE*(NSAVE+3). For the recommended values of NSAVE (10), 
C         RWORK has size at least 131 + 17*N.
C IWORK  :INOUT    Integer IWORK(USER DEFINED >= 30).
C         Used to hold pointers into the RWORK array.
C         Upon return the following locations of IWORK hold information
C         which may be of use to the user:
C         IWORK(9)  Amount of Integer workspace actually used.
C         IWORK(10) Amount of Real    workspace actually used.
C LENIW  :IN       Integer.
C         Length of the integer workspace IWORK.  LENIW >= 30.
C
C *Description:
C       SSDGMR solves a linear system A*X = B rewritten in the form:
C
C        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,
C
C       with right preconditioning, or
C
C        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,
C
C       with left preconditioning, where a is an n-by-n real matrix,
C       X and  B  are N-vectors,  SB and  SX  are  diagonal  scaling
C       matrices, and  M   is   the  diagonal  of   A.     It   uses
C       preconditioned   Krylov  subpace   methods  based    on  the
C       generalized  minimum residual method (GMRES).   This routine
C       is  a  driver routine  which   assumes a  SLAP matrix   data
C       structure  and   sets  up the  necessary information   to do
C       diagonal preconditioning and  calls  the main GMRES  routine
C       SGMRES   for  the  solution  of the   linear system.  SGMRES
C       optionally   performs   either the   full  orthogonalization
C       version of the GMRES algorithm or an  incomplete  variant of
C       it.  Both versions use restarting of the linear iteration by
C       default, although the user can disable this feature.
C       
C       The GMRES  algorithm generates a sequence  of approximations
C       X(L) to the  true solution of the above  linear system.  The
C       convergence criteria for stopping the  iteration is based on
C       the size  of the  scaled norm of  the residual  R(L)  =  B -
C       A*X(L).  The actual stopping test is either:
C
C               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B),
C
C       for right preconditioning, or
C
C               norm(SB*(M-inverse)*(B-A*X(L))) .le. 
C                       TOL*norm(SB*(M-inverse)*B),
C
C       for left preconditioning, where norm() denotes the euclidean
C       norm, and TOL is  a positive scalar less  than one  input by
C       the user.  If TOL equals zero  when SSDGMR is called, then a
C       default  value  of 500*(the   smallest  positive  magnitude,
C       machine epsilon) is used.  If the  scaling arrays SB  and SX
C       are used, then  ideally they  should be chosen  so  that the
C       vectors SX*X(or SX*M*X) and  SB*B have all their  components
C       approximately equal  to  one in  magnitude.  If one wants to
C       use the same scaling in X  and B, then  SB and SX can be the
C       same array in the calling program.
C
C       The following is a list of the other routines and their
C       functions used by GMRES:
C       SGMRES  Contains the matrix structure independent driver 
C               routine for GMRES.
C       SPIGMR  Contains the main iteration loop for GMRES.
C       SORTH   Orthogonalizes a new vector against older basis vects.
C       SHEQR   Computes a QR decomposition of a Hessenberg matrix.
C       SHELS   Solves a Hessenberg least-squares system, using QR 
C               factors.
C       RLCALC  Computes the scaled residual RL.
C       XLCALC  Computes the solution XL.
C       ISSGMR  User-replaceable stopping routine.
C
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures and SLAP  will figure out  which on
C       is being used and act accordingly.
C       
C       =================== S L A P Triad format ===================
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C *Precision:           Single Precision
C *Side Effects:
C       The SLAP Triad format (IA, JA, A) is modified internally to be
C       the SLAP Column format.  See above.
C***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh,
C                 "Reduced Storage Matrix Methods In Stiff ODE 
C                 Systems," LLNL report UCRL-95088, Rev. 1, 
C                 June 1987.
C***ROUTINES CALLED  SS2Y, SCHKW, SSDS, SGMRES
C***END PROLOGUE  SSDGMR
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
      INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL
      INTEGER  ITMAX, ITER, IERR, IUNIT, LENW, IWORK(LENIW), LENIW
      REAL     B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL SSMV, SSDI
      PARAMETER (LOCRB=1, LOCIB=11)
C
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.
C***FIRST EXECUTABLE STATEMENT  SSDGMR
      IERR = 0
      ERR  = 0.0
      IF( NSAVE.LE.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL SS2Y( N, NELT, IA, JA, A, ISYM )
C         
C         Set up the workspace.  We assume MAXL=KMP=NSAVE.
C         Compute the inverse of the diagonal of the matrix.
      LOCIGW = LOCIB
      LOCIW = LOCIGW + 20
C
      LOCDIN = LOCRB
      LOCRGW = LOCDIN + N
      LOCW = LOCRGW + 1+N*(NSAVE+6)+NSAVE*(NSAVE+3)
C
      IWORK(4) = LOCDIN
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
C         Check the workspace allocations.
      CALL SCHKW( 'SSDGMR', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      CALL SSDS(N, NELT, IA, JA, A, ISYM, RWORK(LOCDIN))
C         
C         Perform the Diagonaly Scaled Generalized Minimum
C         Residual iteration algorithm.  The following SGMRES 
C         defaults are used MAXL = KMP = NSAVE, JSCAL = 0,
C         JPRE = -1, NRMAX = ITMAX/NSAVE
      IWORK(LOCIGW  ) = NSAVE
      IWORK(LOCIGW+1) = NSAVE
      IWORK(LOCIGW+2) = 0
      IWORK(LOCIGW+3) = -1
      IWORK(LOCIGW+4) = ITMAX/NSAVE
      MYITOL = 0
C      
      CALL SGMRES( N, B, X, NELT, IA, JA, A, ISYM, SSMV, SSDI,
     $     MYITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, RWORK,
     $     RWORK(LOCRGW), LENW-LOCRGW, IWORK(LOCIGW), 20,
     $     RWORK, IWORK )
C
      IF( ITER.GT.ITMAX ) IERR = 2
      RETURN
C------------- LAST LINE OF SSDGMR FOLLOWS ----------------------------
      END
*DECK SSLUGM
      SUBROUTINE SSLUGM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE,
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $     IWORK, LENIW )
C***BEGIN PROLOGUE  SSLUGM
C***DATE WRITTEN   880615   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSLUGM-S),
C             Non-Symmetric Linear system, Sparse,
C             Iterative Precondition, Generalized Minimum Residual
C***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov
C           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C***PURPOSE  Incomplete LU GMRES iterative sparse Ax=b solver.
C            This routine uses the generalized minimum residual 
C            (GMRES) method with incomplete LU factorization for 
C            preconditioning to solve possibly non-symmetric linear 
C            systems of the form: Ax = b.
C***DESCRIPTION
C *Usage:
C      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE
C      INTEGER   ITOL, ITMAX, IERR, IUNIT, LENW, IWORK(LENIW), LENIW
C      REAL      B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N), 
C      REAL      RWORK(LENW)
C      EXTERNAL  MATVEC, MSOLVE
C
C      CALL SSLUGM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE,
C     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
C     $     RWORK, LENW, IWORK, LENIW)
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand side vector.
C X      :INOUT    Real X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Integer A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "Description", 
C         below.  If the SLAP Triad format is chosen it is changed 
C         internally to the SLAP Column format.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C NSAVE  :IN       Integer.
C         Number of direction vectors to save and orthogonalize against.
C         Must be greater than 1.
C ITOL   :IN       Integer.
C         Flag to indicate the type of convergence criterion used.
C         ITOL=0  Means the  iteration stops when the test described
C                 below on  the  residual RL  is satisfied.  This is
C                 the  "Natural Stopping Criteria" for this routine.
C                 Other values  of   ITOL  cause  extra,   otherwise
C                 unnecessary, computation per iteration and     are
C                 therefore  much less  efficient.  See  ISSGMR (the
C                 stop test routine) for more information.
C         ITOL=1  Means   the  iteration stops   when the first test
C                 described below on  the residual RL  is satisfied,
C                 and there  is either right  or  no preconditioning
C                 being used.
C         ITOL=2  Implies     that   the  user    is   using    left
C                 preconditioning, and the second stopping criterion
C                 below is used.
C         ITOL=3  Means the  iteration stops   when  the  third test
C                 described below on Minv*Residual is satisfied, and
C                 there is either left  or no  preconditioning begin
C                 used.
C         ITOL=11 is    often  useful  for   checking  and comparing
C                 different routines.  For this case, the  user must
C                 supply  the  "exact" solution or  a  very accurate
C                 approximation (one with  an  error much less  than
C                 TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C                 if ITOL=11, iteration stops when the 2-norm of the
C                 difference between the iterative approximation and
C                 the user-supplied solution  divided by the  2-norm
C                 of the  user-supplied solution  is  less than TOL.
C                 Note that this requires  the  user to  set up  the
C                 "COMMON     /SOLBLK/ SOLN(LENGTH)"  in the calling
C                 routine.  The routine with this declaration should
C                 be loaded before the stop test so that the correct
C                 length is used by  the loader.  This procedure  is
C                 not standard Fortran and may not work correctly on
C                 your   system (although  it  has  worked  on every
C                 system the authors have tried).  If ITOL is not 11
C                 then this common block is indeed standard Fortran.
C TOL    :INOUT    Real.
C         Convergence criterion, as described below.  If TOL is set
C         to zero on input, then a default value of 500*(the smallest
C         positive magnitude, machine epsilon) is used.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.  This routine uses the default
C         of NRMAX = ITMAX/NSAVE to determine the when each restart 
C         should occur.  See the description of NRMAX and MAXL in 
C         SGMRES for a full and frightfully interesting discussion of 
C         this topic.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.  Letting norm() denote the Euclidean
C         norm, ERR is defined as follows...
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               for right or no preconditioning, and
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               for left preconditioning.
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               since right or no preconditioning
C                               being used.
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               since left preconditioning is being
C                               used.
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
C                               i=1,n
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
C IERR   :OUT      Integer.
C         Return error flag.
C               IERR = 0 => All went well.
C               IERR = 1 => Insufficient storage allocated for 
C                           RGWK or IGWK.
C               IERR = 2 => Routine SPIGMR failed to reduce the norm
C                           of the current residual on its last call,
C                           and so the iteration has stalled.  In
C                           this case, X equals the last computed
C                           approximation.  The user must either
C                           increase MAXL, or choose a different
C                           initial guess.
C               IERR =-1 => Insufficient length for RGWK array.
C                           IGWK(6) contains the required minimum
C                           length of the RGWK array.
C               IERR =-2 => Inconsistent ITOL and JPRE values.
C         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
C         left-hand-side of the relevant stopping test defined
C         below associated with the residual for the current
C         approximation X(L).
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C RWORK  :WORK    Real RWORK(LENW).
C         Real array of size LENW.
C LENW   :IN       Integer.
C         Length of the real workspace, RWORK. LENW >= 1 + N*(NSAVE+7)
C         +  NSAVE*(NSAVE+3)+NEL+NU. For the recommended values,  RWORK 
C         has size at least 131 + 17*N + NEL + NU.  Where  NEL is  the
C         number of non- zeros  in  the  lower triangle of  the matrix
C         (including the diagonal).  NU is the  number  of nonzeros in
C         the upper triangle of the matrix (including the diagonal).
C IWORK  :INOUT    Integer IWORK(LENIW).
C         Used to hold pointers into the RWORK array.
C         Upon return the following locations of IWORK hold information
C         which may be of use to the user:
C         IWORK(9)  Amount of Integer workspace actually used.
C         IWORK(10) Amount of Real    workspace actually used.
C LENIW  :IN       Integer.
C         Length of the integer workspace, IWORK.
C         LENIW >= NEL+NU+4*N+32.
C
C *Description:
C       SSLUGM solves a linear system A*X = B rewritten in the form:
C
C        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,
C
C       with right preconditioning, or
C
C        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,
C
C       with left preconditioning, where a is an n-by-n real matrix,
C       X and  B are  N-vectors,  SB and  SX  are   diagonal scaling
C       matrices, and M is the Incomplete LU factorization of A.  It
C       uses preconditioned  Krylov subpace   methods  based on  the
C       generalized minimum residual  method (GMRES).   This routine
C       is a  driver  routine  which  assumes a SLAP   matrix   data
C       structure   and  sets  up  the  necessary  information to do
C       diagonal  preconditioning  and calls the main GMRES  routine
C       SGMRES for the   solution   of the linear   system.   SGMRES
C       optionally   performs  either  the full    orthogonalization
C       version of the  GMRES algorithm or  an incomplete variant of
C       it.  Both versions use restarting of the linear iteration by
C       default, although the user can disable this feature.
C       
C       The GMRES  algorithm generates a sequence  of approximations
C       X(L) to the  true solution of the above  linear system.  The
C       convergence criteria for stopping the  iteration is based on
C       the size  of the  scaled norm of  the residual  R(L)  =  B -
C       A*X(L).  The actual stopping test is either:
C
C               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B),
C
C       for right preconditioning, or
C
C               norm(SB*(M-inverse)*(B-A*X(L))) .le. 
C                       TOL*norm(SB*(M-inverse)*B),
C
C       for left preconditioning, where norm() denotes the euclidean
C       norm, and TOL is  a positive scalar less  than one  input by
C       the user.  If TOL equals zero  when SSLUGM is called, then a
C       default  value  of 500*(the   smallest  positive  magnitude,
C       machine epsilon) is used.  If the  scaling arrays SB  and SX
C       are used, then  ideally they  should be chosen  so  that the
C       vectors SX*X(or SX*M*X) and  SB*B have all their  components
C       approximately equal  to  one in  magnitude.  If one wants to
C       use the same scaling in X  and B, then  SB and SX can be the
C       same array in the calling program.
C
C       The following is a list of the other routines and their
C       functions used by GMRES:
C       SGMRES  Contains the matrix structure independent driver 
C               routine for GMRES.
C       SPIGMR  Contains the main iteration loop for GMRES.
C       SORTH   Orthogonalizes a new vector against older basis vects.
C       SHEQR   Computes a QR decomposition of a Hessenberg matrix.
C       SHELS   Solves a Hessenberg least-squares system, using QR 
C               factors.
C       RLCALC  Computes the scaled residual RL.
C       XLCALC  Computes the solution XL.
C       ISSGMR  User-replaceable stopping routine.
C
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures and SLAP  will figure out  which on
C       is being used and act accordingly.
C       
C       =================== S L A P Triad format ===================
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C *Precision:           Single Precision
C *Side Effects:
C       The SLAP Triad format (IA, JA, A) is modified internally to be
C       the SLAP Column format.  See above.
C***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh,
C                 "Reduced Storage Matrix Methods In Stiff ODE 
C                 Systems," LLNL report UCRL-95088, Rev. 1, 
C                 June 1987.
C***ROUTINES CALLED  SS2Y, SCHKW, SSILUS, SGMRES, SSMV, SSLUI
C***END PROLOGUE  SSLUGM
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
       INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL
      INTEGER  ITMAX, ITER, IERR, IUNIT, LENW, IWORK(LENIW), LENIW
      REAL     B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL SSMV, SSLUI
      PARAMETER (LOCRB=1, LOCIB=11)
C
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.
C***FIRST EXECUTABLE STATEMENT  SSLUGM
      IERR = 0
      ERR  = 0.0
      IF( NSAVE.LE.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL SS2Y( N, NELT, IA, JA, A, ISYM )
C
C         Count number of Non-Zero elements preconditioner ILU matrix.
C         Then set up the work arrays.  We assume MAXL=KMP=NSAVE.
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
C         Don't count diagonal.
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
                  IF( ISYM.NE.0 ) NU = NU + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
C         
      LOCIGW = LOCIB
      LOCIL = LOCIGW + 20
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
C
      LOCL = LOCRB
      LOCDIN = LOCL + NL
      LOCU = LOCDIN + N
      LOCRGW = LOCU + NU
      LOCW = LOCRGW + 1+N*(NSAVE+6)+NSAVE*(NSAVE+3)
C
C         Check the workspace allocations.
      CALL SCHKW( 'SSLUGM', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
C         Compute the Incomplete LU decomposition.
      CALL SSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )
C         
C         Perform the Incomplet LU Preconditioned Generalized Minimum
C         Residual iteration algorithm.  The following SGMRES 
C         defaults are used MAXL = KMP = NSAVE, JSCAL = 0,
C         JPRE = -1, NRMAX = ITMAX/NSAVE
      IWORK(LOCIGW  ) = NSAVE
      IWORK(LOCIGW+1) = NSAVE
      IWORK(LOCIGW+2) = 0
      IWORK(LOCIGW+3) = -1
      IWORK(LOCIGW+4) = ITMAX/NSAVE
      MYITOL = 0
C      
      CALL SGMRES( N, B, X, NELT, IA, JA, A, ISYM, SSMV, SSLUI,
     $     MYITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, RWORK,
     $     RWORK(LOCRGW), LENW-LOCRGW, IWORK(LOCIGW), 20,
     $     RWORK, IWORK )
C
      IF( ITER.GT.ITMAX ) IERR = 2
      RETURN
C------------- LAST LINE OF SSLUGM FOLLOWS ----------------------------
      END
*DECK SHELS
      SUBROUTINE SHELS(A, LDA, N, Q, B)
C***BEGIN PROLOGUE  SHEQR
C***DATE WRITTEN   871001   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SHEQR-S),
C             Non-Symmetric Linear system, Sparse,
C             Iterative Precondition, Generalized Minimum Residual
C***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov
C           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C***PURPOSE  Internal routine for SGMRES.
C***DESCRIPTION
C        This routine is extraced from the LINPACK routine SGESL with
C        changes  due to  the fact  that  A is an  upper   Hessenberg
C        matrix.
C        
C        SHELS solves the least squares problem:
C        
C                   MIN(B-A*X,B-A*X)
C        
C        using the factors computed by SHEQR.
C
C *Usage:
C      INTEGER LDA, N
C      REAL A(LDA,1), B(1), Q(1)
C
C      CALL SHELS(A, LDA, N, Q, B)
C
C *Arguments:
C A       :IN       Real A(LDA,N)
C          The output from SHEQR which contains the upper
C          triangular factor R in the QR decomposition of A.
C LDA     :IN       Integer
C          The leading dimension of the array A.
C N       :IN       Integer
C          A is originally an (N+1) by N matrix.
C Q       :IN       Real Q(2*N)
C          The coefficients of the N givens rotations
C          used in the QR factorization of A.
C B       :INOUT    Real B(N+1)
C          On input, B is the right hand side vector.
C          On output, B is the solution vector X.
C *See Also:
C         SGMRES
C
C***ROUTINES CALLED  SAXPY
C***END PROLOGUE  SHEQR
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
      INTEGER LDA, N
      REAL A(LDA,1), B(1), Q(1)
C
C         Local Variables.
C
      INTEGER IQ, K, KB, KP1
      REAL C, S, T, T1, T2
C
C         minimize(B-A*X,B-A*X).  First form Q*B.
C
      DO 20 K = 1, N
         KP1 = K + 1
         IQ = 2*(K-1) + 1
         C = Q(IQ)
         S = Q(IQ+1)
         T1 = B(K)
         T2 = B(KP1)
         B(K) = C*T1 - S*T2
         B(KP1) = S*T1 + C*T2
 20   CONTINUE
C
C         Now solve  R*X = Q*B.
C
      DO 40 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL SAXPY(K-1, T, A(1,K), 1, B(1), 1)
 40   CONTINUE
      RETURN
C------------- LAST LINE OF SHELS FOLLOWS ----------------------------
      END
*DECK SHEQR
      SUBROUTINE SHEQR(A, LDA, N, Q, INFO, IJOB)
C***BEGIN PROLOGUE  SHEQR
C***DATE WRITTEN   871001   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SHEQR-S),
C             Non-Symmetric Linear system, Sparse,
C             Iterative Precondition, Generalized Minimum Residual
C***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov
C           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C***PURPOSE  Internal routine for SGMRES.
C***DESCRIPTION
C        This   routine  performs  a QR   decomposition  of an  upper
C        Hessenberg matrix A using Givens  rotations.  There  are two
C        options  available: 1)  Performing  a fresh decomposition 2)
C        updating the QR factors by adding a row and  a column to the
C        matrix A.
C
C *Usage:
C      INTEGER LDA, N, INFO, IJOB
C      REAL A(LDA,1), Q(1)
C
C      CALL SHEQR(A, LDA, N, Q, INFO, IJOB)
C
C *Arguments:
C A      :INOUT    Real A(LDA,N)
C         On input, the matrix to be decomposed.
C         On output, the upper triangular matrix R.
C         The factorization can be written Q*A = R, where
C         Q is a product of Givens rotations and R is upper
C         triangular.
C LDA    :IN       Integer
C         The leading dimension of the array A.
C N      :IN       Integer
C         A is an (N+1) by N Hessenberg matrix.
C IJOB   :IN       Integer
C         = 1     means that a fresh decomposition of the
C                 matrix A is desired.
C         .ge. 2  means that the current decomposition of A
C                 will be updated by the addition of a row
C                 and a column.
C Q      :OUT      Real Q(2*N)
C         The factors c and s of each Givens rotation used
C         in decomposing A.
C INFO   :OUT      Integer
C         = 0  normal value.
C         = K  if  A(K,K) .eq. 0.0 .  This is not an error
C           condition for this subroutine, but it does
C           indicate that SHELS will divide by zero
C           if called.
C
C *See Also:
C         SGMRES
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SHEQR
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
      INTEGER LDA, N, INFO, IJOB
      REAL A(LDA,1), Q(1)
C
C         Local Variables.
C
      INTEGER I, IQ, J, K, KM1, KP1, NM1
      REAL C, S, T, T1, T2
C
C***FIRST EXECUTABLE STATEMENT  SHEQR
      IF (IJOB .GT. 1) GO TO 70
C   -------------------------------------------------------------------
C         A new facorization is desired.
C   -------------------------------------------------------------------
C         QR decomposition without pivoting.
C
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
C
C           Compute K-th column of R.
C           First, multiply the K-th column of a by the previous
C           K-1 Givens rotations.
C
         IF (KM1 .LT. 1) GO TO 20
         DO 10 J = 1, KM1
            I = 2*(J-1) + 1
            T1 = A(J,K)
            T2 = A(J+1,K)
            C = Q(I)
            S = Q(I+1)
            A(J,K) = C*T1 - S*T2
            A(J+1,K) = S*T1 + C*T2
 10      CONTINUE
C
C         Compute Givens components C and S.
C
 20      CONTINUE
         IQ = 2*KM1 + 1
         T1 = A(K,K)
         T2 = A(KP1,K)
         IF( T2.EQ.0.0E0 ) THEN
            C = 1.0E0
            S = 0.0E0
         ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
            T = T1/T2
            S = -1.0E0/SQRT(1.0E0+T*T)
            C = -S*T
         ELSE
            T = T2/T1
            C = 1.0E0/SQRT(1.0E0+T*T)
            S = -C*T
         ENDIF
         Q(IQ) = C
         Q(IQ+1) = S
         A(K,K) = C*T1 - S*T2
         IF( A(K,K).EQ.0.0E0 ) INFO = K
 60   CONTINUE
      RETURN
C   -------------------------------------------------------------------
C         The old factorization of a will be updated.  A row and a 
C         column has been added to the matrix A.  N by N-1 is now 
C         the old size of the matrix.
C   -------------------------------------------------------------------
 70   CONTINUE
      NM1 = N - 1
C   -------------------------------------------------------------------
C         Multiply the new column by the N previous Givens rotations.
C   -------------------------------------------------------------------
      DO 100 K = 1,NM1
         I = 2*(K-1) + 1
         T1 = A(K,N)
         T2 = A(K+1,N)
         C = Q(I)
         S = Q(I+1)
         A(K,N) = C*T1 - S*T2
         A(K+1,N) = S*T1 + C*T2
 100  CONTINUE
C   -------------------------------------------------------------------
C         Complete update of decomposition by forming last Givens 
C         rotation, and multiplying it times the column 
C         vector(A(N,N),A(NP1,N)).
C   -------------------------------------------------------------------
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF ( T2.EQ.0.0E0 ) THEN
         C = 1.0E0
         S = 0.0E0
      ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
         T = T1/T2
         S = -1.0E0/SQRT(1.0E0+T*T)
         C = -S*T
      ELSE
         T = T2/T1
         C = 1.0E0/SQRT(1.0E0+T*T)
         S = -C*T
      ENDIF
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
C------------- LAST LINE OF SHEQR FOLLOWS ----------------------------
      END
*DECK SORTH
      SUBROUTINE SORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
C***BEGIN PROLOGUE  SORTH
C***DATE WRITTEN   871001   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SORTH-S),
C             Non-Symmetric Linear system, Sparse,
C             Iterative Precondition, Generalized Minimum Residual
C***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov
C           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C***PURPOSE  Internal routine for SGMRES.
C***DESCRIPTION
C        This routine  orthogonalizes  the  vector  VNEW  against the
C        previous KMP  vectors in the   V array.  It uses  a modified
C        gram-schmidt   orthogonalization procedure with  conditional
C        reorthogonalization.
C
C *Usage:
C      INTEGER N, LL, LDHES, KMP
C      REAL VNEW, V, HES, SNORMW
C      DIMENSION VNEW(1), V(N,1), HES(LDHES,1)
C
C      CALL SORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
C
C *Arguments:
C VNEW   :INOUT    Real VNEW(N)
C         On input, the vector of length n containing a scaled
C         product of the jacobian and the vector v(*,ll).
C         On output, the new vector orthogonal to v(*,i0) to v(*,ll),
C         where i0 = max(1, ll-kmp+1).
C V      :IN       Real V(N,1)
C         The n x ll array containing the previous ll
C         orthogonal vectors v(*,1) to v(*,ll).
C HES    :INOUT    Real HES(LDHES,1)
C         On input, an LL x LL upper hessenberg matrix containing,
C         in HES(I,K), K.lt.LL, the scaled inner products of
C         A*V(*,K) and V(*,i).
C         On return, column LL of HES is filled in with
C         the scaled inner products of A*V(*,LL) and V(*,i).
C LDHES  :IN       Integer
C         The leading dimension of the HES array.
C N      :IN       Integer
C         The order of the matrix A, and the length of VNEW.
C LL     :IN       Integer
C         The current order of the matrix HES.
C KMP    :IN       Integer
C         The number of previous vectors the new vector VNEW
C         must be made orthogonal to (KMP .le. MAXL).
C SNORMW :OUT      REAL
C         Scalar containing the l-2 norm of VNEW.
C
C *See Also:
C         SGMRES
C
C***ROUTINES CALLED  SAXPY
C***END PROLOGUE  SORTH
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
      INTEGER N, LL, LDHES, KMP
      REAL VNEW, V, HES, SNORMW
      DIMENSION VNEW(1), V(N,1), HES(LDHES,1)
C
C         Internal variables.
C
      INTEGER I, I0
      REAL ARG, SUMDSQ, TEM, VNRM
C
C         Get norm of unaltered VNEW for later use.
C***FIRST EXECUTABLE STATEMENT  SORTH
      VNRM = SNRM2(N, VNEW, 1)
C   -------------------------------------------------------------------
C         Perform the modified gram-schmidt procedure on VNEW =A*V(LL).
C         Scaled inner products give new column of HES.
C         Projections of earlier vectors are subtracted from VNEW.
C   -------------------------------------------------------------------
      I0 = MAX0(1,LL-KMP+1)
      DO 10 I = I0,LL
         HES(I,LL) = SDOT(N, V(1,I), 1, VNEW, 1)
         TEM = -HES(I,LL)
         CALL SAXPY(N, TEM, V(1,I), 1, VNEW, 1)
 10   CONTINUE
C   -------------------------------------------------------------------
C         Compute SNORMW = norm of VNEW.  If VNEW is small compared 
C         to its input value (in norm), then reorthogonalize VNEW to 
C         V(*,1) through V(*,LL).  Correct if relative correction 
C         exceeds 1000*(unit roundoff).  Finally, correct SNORMW using 
C         the dot products involved.
C   -------------------------------------------------------------------
      SNORMW = SNRM2(N, VNEW, 1)
      IF (VNRM + 0.001E0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0.0E0
      DO 30 I = I0,LL
         TEM = -SDOT(N, V(1,I), 1, VNEW, 1)
         IF (HES(I,LL) + 0.001E0*TEM .EQ. HES(I,LL)) GO TO 30
         HES(I,LL) = HES(I,LL) - TEM
         CALL SAXPY(N, TEM, V(1,I), 1, VNEW, 1)
         SUMDSQ = SUMDSQ + TEM**2
 30   CONTINUE
      IF (SUMDSQ .EQ. 0.0E0) RETURN
      ARG = AMAX1(0.0E0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
C
      RETURN
C------------- LAST LINE OF SORTH FOLLOWS ----------------------------
      END
*DECK SPIGMR
      SUBROUTINE SPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, 
     $     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR,
     $     RPAR, IPAR, WK, DL, RHOL, NRMAX, B, BNRM, X, XL,
     $     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
C***BEGIN PROLOGUE  SPIGMR
C***DATE WRITTEN   871001   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SPIGMR-S),
C             Non-Symmetric Linear system, Sparse,
C             Iterative Precondition, Generalized Minimum Residual
C***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov
C           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C***PURPOSE  Internal routine for SGMRES.
C***DESCRIPTION
C         This routine solves the linear system A * Z = R0 using a 
C         scaled preconditioned version of the generalized minimum 
C         residual method.  An initial guess of Z = 0 is assumed.
C
C *Usage:
C      EXTERNAL MATVEC, MSOLVE
C      INTEGER N,MAXL,MAXLP1,KMP,JPRE,NMSL,LGMR,IPAR,IFLAG,JSCAL,NRSTS
C      INTEGER NRMAX,ITOL,NELT,ISYM
C      REAL R0,SR,SZ,Z,V,HES,Q,RPAR,WK,DL,RHOL,BNRM,TOL,A,B,X
C      DIMENSION R0(1), SR(1), SZ(1), Z(1), V(N,1),
C     $     HES(MAXLP1,1), Q(1), RPAR(1), IPAR(1), WK(1), DL(1),
C     $     IA(NELT), JA(NELT), A(NELT), B(1), X(1), XL(1)
C
C      CALL SPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, 
C     $     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR,
C     $     RPAR, IPAR, WK, DL, RHOL, NRMAX, B, BNRM, X, XL,
C     $     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
C
C *Arguments:
C R0     :IN       Real R0(N)
C         R0 = the right hand side of the system A*Z = R0.
C         R0 is also used as work space when computing
C         the final approximation.
C         (R0 is the same as V(*,MAXL+1) in the call to SPIGMR.)
C SR     :IN       Real SR(N)
C         SR is a vector of length N containing the nonzero
C         elements of the diagonal scaling matrix for R0.
C SZ     :IN       Real SZ(N)
C         SZ is a vector of length N containing the nonzero
C         elements of the diagonal scaling matrix for Z.
C JSCAL  :IN       Integer
C         A flag indicating whether arrays SR and SZ are used.
C         JSCAL=0 means SR and SZ are not used and the
C                 algorithm will perform as if all
C                 SR(i) = 1 and SZ(i) = 1.
C         JSCAL=1 means only SZ is used, and the algorithm
C                 performs as if all SR(i) = 1.
C         JSCAL=2 means only SR is used, and the algorithm
C                 performs as if all SZ(i) = 1.
C         JSCAL=3 means both SR and SZ are used.
C N      :IN       Integer
C         The order of the matrix A, and the lengths
C         of the vectors SR, SZ, R0 and Z.
C MAXL   :IN       Integer
C         The maximum allowable order of the matrix H.
C MAXLP1 :IN       Integer
C         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
C KMP    :IN       Integer
C         The number of previous vectors the new vector VNEW
C         must be made orthogonal to.  (KMP .le. MAXL)
C NRSTS  :IN       Integer
C         Counter for the number of restarts on the current
C         call to SGMRES.  If NRSTS .gt. 0, then the residual
C         R0 is already scaled, and so scaling of it is
C         not necessary.
C JPRE   :IN       Integer
C         Preconditioner type flag.
C WK     :IN       Real WK(N)
C         A real work array of length N used by routine MATVEC
C         and MSOLVE.
C DL     :INOUT    Real DL(N)
C         On input, a real work array of length N used for calculation
C         of the residual norm RHO when the method is incomplete
C         (KMP.lt.MAXL), and/or when using restarting.
C         On output, the scaled residual vector RL.  It is only loaded
C         when performing restarts of the Krylov iteration.
C NRMAX  :IN       Integer
C         The maximum number of restarts of the Krylov iteration.
C         NRMAX .gt. 0 means restarting is active, while
C         NRMAX = 0 means restarting is not being used.
C B      :IN       Real B(N)
C         The right hand side of the linear system A*X = B.
C BNRM   :IN       Real
C         The scaled norm of b.
C X      :IN       Real X(N)
C         The current approximate solution as of the last
C         restart.
C XL     :IN       Real XL(N)
C         An array of length N used to hold the approximate
C         solution X(L) when ITOL=11.
C ITOL   :IN       Integer
C         A flag to indicate the type of convergence criterion
C         used.  see the driver for its description.
C TOL    :IN       Real
C         The tolerance on residuals R0-A*Z in scaled norm.
C NELT   :IN       Integer
C         The length of arrays IA, JA and A.
C IA     :IN       Integer IA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C JA     :IN       Integer JA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C A      :IN       Real A(NELT)
C         A real array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C ISYM   :IN       Integer
C         A flag to indicate symmetric matrix storage.
C         If ISYM=0, all nonzero entries of the matrix are
C         stored.  If ISYM=1, the matrix is symmetric and
C         only the upper or lower triangular part is stored.
C IUNIT  :IN       Integer
C         The i/o unit number for writing intermediate residual
C         norm values.
C Z      :OUT      Real Z(N)
C         The final computed approximation to the solution
C         of the system A*Z = R0.
C LGMR   :OUT      Integer
C         The number of iterations performed and
C         the current order of the upper hessenberg
C         matrix HES.
C RPAR   :IN       Real RPAR(*)
C         Real work space passed directly to the MSOLVE routine.
C IPAR   :IN       Integer IPAR(*)
C         Integer work space passed directly to the MSOLVE
C         routine.
C NMSL   :OUT      Integer
C         The number of calls to MSOLVE.
C V      :OUT      Real V(N,MAXLP1)
C         The N by (LGMR+1) array containing the LGMR
C         orthogonal vectors V(*,1) to V(*,LGMR).
C HES    :OUT      Real HES(MAXLP1,MAXL)
C         The upper triangular factor of the QR decomposition
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C         entries are the scaled inner-products of A*V(*,I)
C         and V(*,K).
C Q      :OUT      Real Q(2*MAXL)
C         A real array of length 2*MAXL containing the components
C         of the Givens rotations used in the QR decomposition
C         of HES.  It is loaded in SHEQR and used in SHELS.
C RHOL   :OUT      Real
C         A real scalar containing the norm of the final residual.
C IFLAG  :OUT      Integer
C         An integer error flag..
C         0 means convergence in LGMR iterations, LGMR.le.MAXL.
C         1 means the convergence test did not pass in MAXL
C           iterations, but the residual norm is .lt. norm(R0),
C           and so Z is computed.
C         2 means the convergence test did not pass in MAXL
C           iterations, residual .ge. norm(R0), and Z = 0.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.
C
C *See Also:
C         SGMRES
C
C***ROUTINES CALLED  ISSGMR, MATVEC, MSOLVE, SORTH, SRLCAL, SHELS, 
C                    SHEQR, SXLCAL, SAXPY, SCOPY, SSCAL, 
C***END PROLOGUE  SPIGMR
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
      EXTERNAL MATVEC, MSOLVE
      INTEGER N,MAXL,MAXLP1,KMP,JPRE,NMSL,LGMR,IFLAG,JSCAL,NRSTS
      INTEGER NRMAX,ITOL,NELT,ISYM
      REAL RHOL,BNRM,TOL
      REAL R0(*), SR(*), SZ(*), Z(*), V(N,*), HES(MAXLP1,*)
      REAL Q(*), RPAR(*), WK(*), DL(*)
      REAL A(NELT), B(*), X(*), XL(*)
      INTEGER IPAR(*), IA(NELT), JA(NELT)
C
C         Local variables.
C
      INTEGER I, INFO, IP1, I2, J, K, LL, LLP1
      REAL R0NRM,C,DLNRM,PROD,RHO,S,SNORMW,TEM
C
C         Zero out the z array.
C***FIRST EXECUTABLE STATEMENT  SPIGMR
      DO 5 I = 1,N
         Z(I) = 0.0E0
 5    CONTINUE
C
      IFLAG = 0
      LGMR = 0
      NMSL = 0
C         Load ITMAX, the maximum number of iterations.
      ITMAX =(NRMAX+1)*MAXL
C   -------------------------------------------------------------------
C         The initial residual is the vector R0.
C         Apply left precon. if JPRE < 0 and this is not a restart.
C         Apply scaling to R0 if JSCAL = 2 or 3.
C   -------------------------------------------------------------------
      IF ((JPRE .LT. 0) .AND.(NRSTS .EQ. 0)) THEN
         CALL SCOPY(N, R0, 1, WK, 1)
         CALL MSOLVE(N, WK, R0, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      IF (((JSCAL.EQ.2) .OR.(JSCAL.EQ.3)) .AND.(NRSTS.EQ.0)) THEN
         DO 10 I = 1,N
            V(I,1) = R0(I)*SR(I)
 10      CONTINUE
      ELSE
         DO 20 I = 1,N
            V(I,1) = R0(I)
 20      CONTINUE
      ENDIF
      R0NRM = SNRM2(N, V, 1)
      ITER = NRSTS*MAXL
C
C         Call stopping routine ISSGMR.
C
      IF (ISSGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
     $    NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, V(1,1), Z, WK,
     $    RPAR, IPAR, R0NRM, BNRM, SR, SZ, JSCAL,
     $    KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
     $    HES, JPRE) .NE. 0) RETURN
      TEM = 1.0E0/R0NRM
      CALL SSCAL(N, TEM, V(1,1), 1)
C
C         Zero out the HES array.
C
      DO 50 J = 1,MAXL
         DO 40 I = 1,MAXLP1
            HES(I,J) = 0.0E0
 40      CONTINUE
 50   CONTINUE
C   -------------------------------------------------------------------
C         main loop to compute the vectors V(*,2) to V(*,MAXL).
C         The running product PROD is needed for the convergence test.
C   -------------------------------------------------------------------
      PROD = 1.0E0
      DO 90 LL = 1,MAXL
         LGMR = LL
C   -------------------------------------------------------------------
C        Unscale  the  current V(LL)  and store  in WK.  Call routine
C        msolve    to   compute(M-inverse)*WK,   where    M   is  the
C        preconditioner matrix.  Save the answer in Z.   Call routine
C        MATVEC to compute  VNEW  = A*Z,  where  A is  the the system
C        matrix.  save the answer in  V(LL+1).  Scale V(LL+1).   Call
C        routine SORTH  to  orthogonalize the    new vector VNEW   =
C        V(*,LL+1).  Call routine SHEQR to update the factors of HES.
C   -------------------------------------------------------------------
        IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
           DO 60 I = 1,N
              WK(I) = V(I,LL)/SZ(I)
 60        CONTINUE
        ELSE
           CALL SCOPY(N, V(1,LL), 1, WK, 1)
        ENDIF
        IF (JPRE .GT. 0) THEN
           CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
           NMSL = NMSL + 1
           CALL MATVEC(N, Z, V(1,LL+1), NELT, IA, JA, A, ISYM)
        ELSE
           CALL MATVEC(N, WK, V(1,LL+1), NELT, IA, JA, A, ISYM)
        ENDIF
        IF (JPRE .LT. 0) THEN
           CALL SCOPY(N, V(1,LL+1), 1, WK, 1)
           CALL MSOLVE(N,WK,V(1,LL+1),NELT,IA,JA,A,ISYM,RPAR,IPAR)
           NMSL = NMSL + 1
        ENDIF
        IF ((JSCAL .EQ. 2) .OR.(JSCAL .EQ. 3)) THEN
           DO 65 I = 1,N
              V(I,LL+1) = V(I,LL+1)*SR(I)
 65        CONTINUE
        ENDIF
        CALL SORTH(V(1,LL+1), V, HES, N, LL, MAXLP1, KMP, SNORMW)
        HES(LL+1,LL) = SNORMW
        CALL SHEQR(HES, MAXLP1, LL, Q, INFO, LL)
        IF (INFO .EQ. LL) GO TO 120
C   -------------------------------------------------------------------
C         Update RHO, the estimate of the norm of the residual R0-A*ZL.
C         If KMP <  MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
C         necessarily orthogonal for LL > KMP.  The vector DL must then
C         be computed, and its norm used in the calculation of RHO.
C   -------------------------------------------------------------------
        PROD = PROD*Q(2*LL)
        RHO = ABS(PROD*R0NRM)
        IF ((LL.GT.KMP) .AND.(KMP.LT.MAXL)) THEN
           IF (LL .EQ. KMP+1) THEN
              CALL SCOPY(N, V(1,1), 1, DL, 1)
              DO 75 I = 1,KMP
                 IP1 = I + 1
                 I2 = I*2
                 S = Q(I2)
                 C = Q(I2-1)
                 DO 70 K = 1,N
                    DL(K) = S*DL(K) + C*V(K,IP1)
 70              CONTINUE
 75           CONTINUE
           ENDIF
           S = Q(2*LL)
           C = Q(2*LL-1)/SNORMW
           LLP1 = LL + 1
           DO 80 K = 1,N
              DL(K) = S*DL(K) + C*V(K,LLP1)
 80        CONTINUE
           DLNRM = SNRM2(N, DL, 1)
           RHO = RHO*DLNRM
        ENDIF
        RHOL = RHO
C   -------------------------------------------------------------------
C         Test for convergence.  If passed, compute approximation ZL.
C         If failed and LL < MAXL, then continue iterating.
C   -------------------------------------------------------------------
        ITER = NRSTS*MAXL + LGMR
        IF (ISSGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
     $      NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, DL, Z, WK,
     $      RPAR, IPAR, RHOL, BNRM, SR, SZ, JSCAL,
     $      KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
     $      HES, JPRE) .NE. 0) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
C   -------------------------------------------------------------------
C         Rescale so that the norm of V(1,LL+1) is one.
C   -------------------------------------------------------------------
        TEM = 1.0E0/SNORMW
        CALL SSCAL(N, TEM, V(1,LL+1), 1)
 90   CONTINUE
 100  CONTINUE
      IF (RHO .LT. R0NRM) GO TO 150
 120  CONTINUE
      IFLAG = 2
C
C         Load approximate solution with zero.
C
      DO 130 I = 1,N
         Z(I) = 0.E0
 130  CONTINUE
      RETURN
 150  IFLAG = 1
C
C         Tolerance not met, but residual norm reduced.
C
      IF (NRMAX .GT. 0) THEN
C
C        If performing restarting (NRMAX > 0)  calculate the residual
C        vector RL and  store it in the DL  array.  If the incomplete
C        version is being used (KMP < MAXL) then DL has  already been
C        calculated up to a scaling factor.   Use SRLCAL to calculate
C        the scaled residual vector.
C
         CALL SRLCAL(N, KMP, MAXL, MAXL, V, Q, DL, SNORMW, PROD,
     $        R0NRM)
      ENDIF
C   -------------------------------------------------------------------
C         Compute the approximation ZL to the solution.  Since the 
C         vector Z was used as work space, and the initial guess
C         of the linear iteration is zero, Z must be reset to zero.
C   -------------------------------------------------------------------
 200  CONTINUE
      LL = LGMR
      LLP1 = LL + 1
      DO 210 K = 1,LLP1
         R0(K) = 0.0E0
 210  CONTINUE
      R0(1) = R0NRM
      CALL SHELS(HES, MAXLP1, LL, Q, R0)
      DO 220 K = 1,N
         Z(K) = 0.0E0
 220  CONTINUE
      DO 230 I = 1,LL
         CALL SAXPY(N, R0(I), V(1,I), 1, Z, 1)
 230  CONTINUE
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
         DO 240 I = 1,N
            Z(I) = Z(I)/SZ(I)
 240     CONTINUE
      ENDIF
      IF (JPRE .GT. 0) THEN
         CALL SCOPY(N, Z, 1, WK, 1)
         CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      RETURN
C------------- LAST LINE OF SPIGMR FOLLOWS ----------------------------
      END
*DECK SRLCAL
      SUBROUTINE SRLCAL(N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD,
     $     R0NRM)
C***BEGIN PROLOGUE  SRLCAL
C***DATE WRITTEN   871001   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SRLCAL-S),
C             Non-Symmetric Linear system, Sparse,
C             Iterative Precondition, Generalized Minimum Residual
C***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov
C           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C***PURPOSE  Internal routine for SGMRES.
C***DESCRIPTION
C         This routine calculates the scaled residual RL from the 
C         V(I)'s.
C *Usage:
C      INTEGER N, KMP, LL, MAXL
C      REAL V, Q, RL, SNORMW
C      DIMENSION V(N,1), Q(1), RL(N)
C
C      CALL SRLCAL(N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD,
C     $     R0NRM)     
C
C *Arguments:
C N      :IN       Integer
C         The order of the matrix A, and the lengths
C         of the vectors SR, SZ, R0 and Z.
C KMP    :IN       Integer
C         The number of previous V vectors the new vector VNEW
C         must be made orthogonal to. (KMP .le. MAXL)
C LL     :IN       Integer
C         The current dimension of the Krylov subspace.
C MAXL   :IN       Integer
C         The maximum dimension of the Krylov subspace.
C Q      :IN       Real Q(2*MAXL)
C         A real array of length 2*MAXL containing the components
C         of the Givens rotations used in the QR decomposition
C         of HES.  It is loaded in SHEQR and used in SHELS.
C PROD   :IN       Real
C        The product s1*s2*...*sl = the product of the sines of the
C        givens rotations used in the QR factorization of
C        the hessenberg matrix HES.
C R0NRM  :IN       Real
C         The scaled norm of initial residual R0.
C RL     :OUT      Real RL(N)
C         The residual vector RL.  This is either SB*(B-A*XL) if
C         not preconditioning or preconditioning on the right,
C         or SB*(M-inverse)*(B-A*XL) if preconditioning on the
C         left.
C
C *See Also:
C         SGMRES
C
C***ROUTINES CALLED  SCOPY, SSCAL
C***END PROLOGUE  SRLCAL
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
      INTEGER N, KMP, LL, MAXL
      REAL V, Q, RL, SNORMW
      DIMENSION V(N,1), Q(1), RL(N)
C
C         Internal Variables.
C
      INTEGER I, IP1, I2, K
C
C***FIRST EXECUTABLE STATEMENT  SRLCAL
      IF (KMP .EQ. MAXL) THEN
C
C         calculate RL.  Start by copying V(*,1) into RL.
C
         CALL SCOPY(N, V(1,1), 1, RL, 1)
         LLM1 = LL - 1
         DO 20 I = 1,LLM1
            IP1 = I + 1
            I2 = I*2
            S = Q(I2)
            C = Q(I2-1)
            DO 10 K = 1,N
               RL(K) = S*RL(K) + C*V(K,IP1)
 10         CONTINUE
 20      CONTINUE
         S = Q(2*LL)
         C = Q(2*LL-1)/SNORMW
         LLP1 = LL + 1
         DO 30 K = 1,N
            RL(K) = S*RL(K) + C*V(K,LLP1)
 30      CONTINUE
      ENDIF
C
C         When KMP < MAXL, RL vector already partially calculated. 
C         Scale RL by R0NRM*PROD to obtain the residual RL.
C
      TEM = R0NRM*PROD
      CALL SSCAL(N, TEM, RL, 1)
      RETURN
C------------- LAST LINE OF SRLCAL FOLLOWS ----------------------------
      END
*DECK SXLCAL
      SUBROUTINE SXLCAL(N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,
     $     WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR,
     $     NELT, IA, JA, A, ISYM)
C***BEGIN PROLOGUE  SXLCAL
C***DATE WRITTEN   871001   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SXLCAL-S),
C             Non-Symmetric Linear system, Sparse,
C             Iterative Precondition, Generalized Minimum Residual
C***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov
C           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C***PURPOSE  Internal routine for SGMRES.
C***DESCRIPTION
C        This  routine computes the solution  XL,  the current SGMRES
C        iterate, given the  V(I)'s and  the  QR factorization of the
C        Hessenberg  matrix HES.   This routine  is  only called when
C        ITOL=11.
C
C *Usage:
C      EXTERNAL MSOLVE
C      INTEGER N, LGMR, MAXLP1, JSCAL, JPRE, IPAR, NMSL, NELT, IA,
C     $        JA, ISYM
C      REAL X, XL, ZL, HES, Q, V, R0NRM, WK, SZ, RPAR, A
C      DIMENSION X(N), XL(N), ZL(N), HES(MAXLP1,1), Q(1), V(N,1), WK(N),
C     $          SZ(1), RPAR(1), IPAR(1), IA(NELT), JA(NELT), A(NELT)
C      CALL SXLCAL(N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,
C     $     WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR,
C     $     NELT, IA, JA, A, ISYM)
C
C *Arguments:
C N      :IN       Integer
C         The order of the matrix A, and the lengths
C         of the vectors SR, SZ, R0 and Z.
C LGMR   :IN       Integer
C         The number of iterations performed and
C         the current order of the upper Hessenberg
C         matrix HES.
C X      :IN       Real X(N)
C         The current approximate solution as of the last restart.
C ZL     :IN       Real ZL(N)
C         An array of length N used to hold the approximate
C         solution Z(L).
C SZ     :IN       Real SZ(N)
C         A vector of length N containing the nonzero
C         elements of the diagonal scaling matrix for Z.
C JSCAL  :IN       Integer
C         A flag indicating whether arrays SR and SZ are used.
C         JSCAL=0 means SR and SZ are not used and the
C                 algorithm will perform as if all
C                 SR(i) = 1 and SZ(i) = 1.
C         JSCAL=1 means only SZ is used, and the algorithm
C                 performs as if all SR(i) = 1.
C         JSCAL=2 means only SR is used, and the algorithm
C                 performs as if all SZ(i) = 1.
C         JSCAL=3 means both SR and SZ are used.
C MAXLP1 :IN       Integer
C         MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
C         MAXL is the maximum allowable order of the matrix HES.
C JPRE   :IN       Integer
C         The preconditioner type flag.
C WK     :IN       Real WK(N)
C         A real work array of length N.
C NMSL   :IN       Integer
C         The number of calls to MSOLVE.
C V      :IN       Real V(N,MAXLP1)
C         The N by(LGMR+1) array containing the LGMR
C         orthogonal vectors V(*,1) to V(*,LGMR).
C HES    :IN       Real HES(MAXLP1,MAXL)
C         The upper triangular factor of the QR decomposition
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C         entries are the scaled inner-products of A*V(*,i) and V(*,k).
C Q      :IN       Real Q(2*MAXL)
C         A real array of length 2*MAXL containing the components
C         of the givens rotations used in the QR decomposition
C         of HES.  It is loaded in SHEQR.
C R0NRM  :IN       Real
C         The scaled norm of the initial residual for the
C         current call to SPIGMR.
C RPAR   :IN       Real RPAR(*)
C         Real work space passed directly to the MSOLVE routine.
C IPAR   :IN       Integer IPAR(*)
C         Integer work space passed directly to the MSOLVE
C         routine.
C NELT   :IN       Integer
C         The length of arrays IA, JA and A.
C IA     :IN       Integer IA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C JA     :IN       Integer JA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C A      :IN       Real A(NELT)
C         A real array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C ISYM   :IN       Integer
C         A flag to indicate symmetric matrix storage.
C         If ISYM=0, all nonzero entries of the matrix are
C         stored.  If ISYM=1, the matrix is symmetric and
C         only the upper or lower triangular part is stored.
C XL     :OUT      Real XL(N)
C         An array of length N used to hold the approximate
C         solution X(L).
C         Warning: XL and ZL are the same array in the calling routine.
C
C *See Also:
C         SGMRES
C
C***ROUTINES CALLED  MSOLVE, SHELS, SAXPY, SCOPY, SSCAL
C***END PROLOGUE  SXLCAL
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
      EXTERNAL MSOLVE
      INTEGER N, LGMR, MAXLP1, JSCAL, JPRE, IPAR, NMSL, NELT, IA,
     1        JA, ISYM
      REAL X, XL, ZL, HES, Q, V, R0NRM, WK, SZ, RPAR, A
      DIMENSION X(N), XL(N), ZL(N), HES(MAXLP1,*), Q(*), V(N,*), WK(N),
     1          SZ(*), RPAR(*), IPAR(*), IA(NELT), JA(NELT), A(NELT)
C
C         Internal variables.
C
      INTEGER I, K, LL, LLP1
C
C***FIRST EXECUTABLE STATEMENT  SXLCAL
      LL = LGMR
      LLP1 = LL + 1
      DO 10 K = 1,LLP1
         WK(K) = 0.0E0
 10   CONTINUE
      WK(1) = R0NRM
      CALL SHELS(HES, MAXLP1, LL, Q, WK)
      DO 20 K = 1,N
         ZL(K) = 0.0E0
 20   CONTINUE
      DO 30 I = 1,LL
         CALL SAXPY(N, WK(I), V(1,I), 1, ZL, 1)
 30   CONTINUE
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
         DO 40 K = 1,N
            ZL(K) = ZL(K)/SZ(K)
 40      CONTINUE
      ENDIF
      IF (JPRE .GT. 0) THEN
         CALL SCOPY(N, ZL, 1, WK, 1)
         CALL MSOLVE(N, WK, ZL, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
C         calculate XL from X and ZL.
      DO 50 K = 1,N
         XL(K) = X(K) + ZL(K)
 50   CONTINUE
      RETURN
C------------- LAST LINE OF SXLCAL FOLLOWS ----------------------------
      END
*DECK ISSGMR
      FUNCTION ISSGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
     $     NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,
     $     RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL,
     $     KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
     $     HES, JPRE)
C***BEGIN PROLOGUE ISSGMR
C***DATE WRITTEN   871211  (YYMMDD)
C***REVISION DATE  881213  (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=INTEGER(ISSGMR-I)
C             Linear system, Sparse, Stop Test, GMRES
C***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov
C           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C***PURPOSE Generalized Minimum Residual Stop Test.
C         This routine calculates the stop test for the Generalized
C         Minimum RESidual (GMRES) iteration scheme.   It returns a
C         nonzero if  the  error  estimate (the  type  of  which is
C         determined  by   ITOL)  is  less  than the user specified 
C         tolerence TOL.
C***DESCRIPTION
C *Usage:
C      INTEGER KMP, LGMR, MAXL, MAXLP1, JPRE, NMSL
C      REAL DXNRM, RNRM, R0NRM, SNORMW, SOLNRM, PROD
C      DIMENSION B(1), X(1), IA(1), JA(1), A(1), R(1), Z(1), DZ(1)
C      DIMENSION RWORK(1), IWORK(1), SB(1), SX(1), Q(1), V(N,1)
C      DIMENSION HES(MAXLP1,MAXL), XL(1)
C      EXTERNAL MSOLVE
C
C      IF (ISSGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
C     $     NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,
C     $     RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL,
C     $     KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
C     $     HES, JPRE) .NE. 0) THEN ITERATION DONE
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand-side vector.
C X      :IN       Real X(N).
C         Approximate solution vector as of the last restart.
C XL     :OUT      Real XL(N)
C         An array of length N used to hold the approximate
C         solution as of the current iteration.  Only computed by
C         this routine when ITOL=11.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Real A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See "Description", in the SGMRES,
C         SSLUGM and SSDGMR routines for more late breaking details...
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MSOLVE :EXT      External.
C         Name of a routine which solves a linear system Mz = r for  z
C         given r with the preconditioning matrix M (M is supplied via
C         RWORK and IWORK arrays.  The name of the MSOLVE routine must 
C         be declared external in the calling program.  The calling 
C         sequence to MSLOVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         Where N is the number of unknowns, R is the right-hand side 
C         vector, and z is the solution upon return.  RWORK is a real 
C         array that can be used to pass necessary preconditioning 
C         information and/or workspace to MSOLVE.  IWORK is an integer 
C         work array for the same purpose as RWORK.
C NMSL   :INOUT    Integer.
C         A counter for the number of calls to MSOLVE.
C ITOL   :IN       Integer.
C         Flag to indicate the type of convergence criterion used.
C         ITOL=0  Means the  iteration stops when the test described
C                 below on  the  residual RL  is satisfied.  This is
C                 the  "Natural Stopping Criteria" for this routine.
C                 Other values  of   ITOL  cause  extra,   otherwise
C                 unnecessary, computation per iteration and     are
C                 therefore  much less  efficient.  See  ISSGMR (the
C                 stop test routine) for more information.
C         ITOL=1  Means   the  iteration stops   when the first test
C                 described below on  the residual RL  is satisfied,
C                 and there  is either right  or  no preconditioning
C                 being used.
C         ITOL=2  Implies     that   the  user    is   using    left
C                 preconditioning, and the second stopping criterion
C                 below is used.
C         ITOL=3  Means the  iteration stops   when  the  third test
C                 described below on Minv*Residual is satisfied, and
C                 there is either left  or no  preconditioning begin
C                 used.
C         ITOL=11 is    often  useful  for   checking  and comparing
C                 different routines.  For this case, the  user must
C                 supply  the  "exact" solution or  a  very accurate
C                 approximation (one with  an  error much less  than
C                 TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C                 if ITOL=11, iteration stops when the 2-norm of the
C                 difference between the iterative approximation and
C                 the user-supplied solution  divided by the  2-norm
C                 of the  user-supplied solution  is  less than TOL.
C                 Note that this requires  the  user to  set up  the
C                 "COMMON     /SOLBLK/ SOLN(LENGTH)"  in the calling
C                 routine.  The routine with this declaration should
C                 be loaded before the stop test so that the correct
C                 length is used by  the loader.  This procedure  is
C                 not standard Fortran and may not work correctly on
C                 your   system (although  it  has  worked  on every
C                 system the authors have tried).  If ITOL is not 11
C                 then this common block is indeed standard Fortran.
C TOL    :IN       Real.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :IN       Integer.
C         The iteration for which to check for convergence.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.  Letting norm() denote the Euclidean
C         norm, ERR is defined as follows..
C
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               for right or no preconditioning, and
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               for left preconditioning.
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               since right or no preconditioning
C                               being used.
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               since left preconditioning is being
C                               used.
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
C                               i=1,n
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C R      :INOUT    Real R(N).
C         Work array used in calling routine.  It contains
C         information necessary to compute the residual RL = B-A*XL.
C Z      :WORK     Real Z(N).
C         Workspace used to hold the pseudo-residule M z = r.
C DZ     :WORK     Real DZ(N).
C         Workspace used to hold temporary vector(s).
C RWORK  :WORK     Real RWORK(USER DEFINABLE).
C         Real array that can be used by MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINABLE).
C         Integer array that can be used by MSOLVE.
C RNRM   :IN       Real.
C         Norm of the current residual.  Type of norm depends on ITOL.
C BNRM   :IN       Real.
C         Norm of the right hand side.  Type of norm depends on ITOL.
C SB     :IN       Real SB(N).
C         Scaling vector for B.
C SX     :IN       Real SX(N).
C         Scaling vector for X.
C JSCAL  :IN       Integer.
C         Flag indicating if scaling arrays SB and SX are being
C         used in the calling routine SPIGMR.
C         JSCAL=0 means SB and SX are not used and the
C                 algorithm will perform as if all
C                 SB(i) = 1 and SX(i) = 1.
C         JSCAL=1 means only SX is used, and the algorithm
C                 performs as if all SB(i) = 1.
C         JSCAL=2 means only SB is used, and the algorithm
C                 performs as if all SX(i) = 1.
C         JSCAL=3 means both SB and SX are used.
C KMP    :IN       Integer
C         The number of previous vectors the new vector VNEW
C         must be made orthogonal to.  (KMP .le. MAXL)
C LGMR   :IN       Integer
C         The number of GMRES iterations performed on the current call 
C         to SPIGMR (i.e., # iterations since the last restart) and
C         the current order of the upper hessenberg
C         matrix HES.
C MAXL   :IN       Integer
C         The maximum allowable order of the matrix H.
C MAXLP1 :IN       Integer
C         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
C V      :IN       Real V(N,MAXLP1)
C         The N by (LGMR+1) array containing the LGMR
C         orthogonal vectors V(*,1) to V(*,LGMR).
C Q      :IN       Real Q(2*MAXL)
C         A real array of length 2*MAXL containing the components
C         of the Givens rotations used in the QR decomposition
C         of HES.
C SNORMW :IN       Real
C         A scalar containing the scaled norm of VNEW before it
C         is renormalized in SPIGMR.
C PROD   :IN       Real
C        The product s1*s2*...*sl = the product of the sines of the
C        givens rotations used in the QR factorization of
C        the hessenberg matrix HES.
C R0NRM  :IN       Real
C         The scaled norm of initial residual R0.
C HES    :IN       Real HES(MAXLP1,MAXL)
C         The upper triangular factor of the QR decomposition
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C         entries are the scaled inner-products of A*V(*,I)
C         and V(*,K).
C JPRE   :IN       Integer
C         Preconditioner type flag.
C
C *Description
C       When using the GMRES solver,  the preferred value  for ITOL
C       is 0.  This is due to the fact that when ITOL=0 the norm of
C       the residual required in the stopping test is  obtained for
C       free, since this value is already  calculated  in the GMRES
C       algorithm.   The  variable  RNRM contains the   appropriate
C       norm, which is equal to norm(SB*(RL - A*XL))  when right or
C       no   preconditioning is  being  performed,   and equal   to
C       norm(SB*Minv*(RL - A*XL))  when using left preconditioning.
C       Here, norm() is the Euclidean norm.  Nonzero values of ITOL
C       require  additional work  to  calculate the  actual  scaled
C       residual  or its scaled/preconditioned  form,  and/or   the
C       approximate solution XL.  Hence, these values of  ITOL will
C       not be as efficient as ITOL=0.
C
C***ROUTINES CALLED     MSOLVE, SNRM2, SCOPY, 
C***END PROLOG ISSGMR
      INTEGER KMP, LGMR, MAXL, MAXLP1, JPRE, NMSL
      REAL DXNRM, RNRM, R0NRM, SNORMW, SOLNRM, PROD
      DIMENSION B(*), X(*), IA(*), JA(*), A(*), R(*), Z(*), DZ(*)
      DIMENSION RWORK(*), IWORK(*), SB(*), SX(*), Q(*), V(N,*)
      DIMENSION HES(MAXLP1,MAXL), XL(*)
      EXTERNAL MSOLVE
      COMMON /SOLBLK/ SOLN(1)
      SAVE SOLNRM
C         
C***FIRST EXECUTABLE STATEMENT ISSGMR
      ISSGMR = 0
      IF ( ITOL.EQ.0 ) THEN
C
C       Use input from SPIGMR to determine if stop conditions are met.
C
         ERR = RNRM/BNRM
      ENDIF
      IF ( (ITOL.GT.0) .AND. (ITOL.LE.3) ) THEN
C
C       Use SRLCAL to calculate the scaled residual vector. 
C       Store answer in R.
C
         IF ( LGMR.NE.0 ) CALL SRLCAL(N, KMP, LGMR, MAXL, V, Q, R,
     $                                SNORMW, PROD, R0NRM)
         IF ( ITOL.LE.2 ) THEN
C         err = ||Residual||/||RightHandSide||(2-Norms).
            ERR = SNRM2(N, R, 1)/BNRM
C
C         Unscale R by R0NRM*PROD when KMP < MAXL.
C
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN
               TEM = 1.0E0/(R0NRM*PROD)
               CALL SSCAL(N, TEM, R, 1)
            ENDIF
         ELSEIF ( ITOL.EQ.3 ) THEN
C         err = Max |(Minv*Residual)(i)/x(i)|
C         When jpre .lt. 0, r already contains Minv*Residual.
            IF ( JPRE.GT.0 ) THEN
               CALL MSOLVE(N, R, DZ, NELT, IA, JA, A, ISYM, RWORK,
     $              IWORK)
               NMSL = NMSL + 1
            ENDIF
C
C         Unscale R by R0NRM*PROD when KMP < MAXL.
C
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN
               TEM = 1.0E0/(R0NRM*PROD)
               CALL SSCAL(N, TEM, R, 1)
            ENDIF
C
            FUZZ = R1MACH(1)
            IELMAX = 1
            RATMAX = ABS(DZ(1))/AMAX1(ABS(X(1)),FUZZ)
            DO 25 I = 2, N
               RAT = ABS(DZ(I))/AMAX1(ABS(X(I)),FUZZ)
               IF( RAT.GT.RATMAX ) THEN
                  IELMAX = I
                  RATMAX = RAT
               ENDIF
 25         CONTINUE
            ERR = RATMAX
            IF( RATMAX.LE.TOL ) ISSGMR = 1
            IF( IUNIT.GT.0 ) WRITE(IUNIT,1020) ITER, IELMAX, RATMAX
            RETURN
         ENDIF
      ENDIF
      IF ( ITOL.EQ.11 ) THEN
C
C       Use SXLCAL to calculate the approximate solution XL.
C
         IF ( (LGMR.NE.0) .AND. (ITER.GT.0) ) THEN
            CALL SXLCAL(N, LGMR, X, XL, XL, HES, MAXLP1, Q, V, R0NRM,
     $           DZ, SX, JSCAL, JPRE, MSOLVE, NMSL, RWORK, IWORK,
     $           NELT, IA, JA, A, ISYM)
         ELSEIF ( ITER.EQ.0 ) THEN
C         Copy X to XL to check if initial guess is good enough.
            CALL SCOPY(N, X, 1, XL, 1)
         ELSE
C         Return since this is the first call to SPIGMR on a restart.
            RETURN
         ENDIF
C
         IF ((JSCAL .EQ. 0) .OR.(JSCAL .EQ. 2)) THEN
C         err = ||x-TrueSolution||/||TrueSolution||(2-Norms).
            IF ( ITER.EQ.0 ) SOLNRM = SNRM2(N, SOLN, 1)
            DO 30 I = 1, N
               DZ(I) = XL(I) - SOLN(I)
 30         CONTINUE
            ERR = SNRM2(N, DZ, 1)/SOLNRM
         ELSE
            IF (ITER .EQ. 0) THEN
               SOLNRM = 0.E0
               DO 40 I = 1,N
                  SOLNRM = SOLNRM + (SX(I)*SOLN(I))**2
 40            CONTINUE
               SOLNRM = SQRT(SOLNRM)
            ENDIF
            DXNRM = 0.E0
            DO 50 I = 1,N
               DXNRM = DXNRM + (SX(I)*(XL(I)-SOLN(I)))**2
 50         CONTINUE
            DXNRM = SQRT(DXNRM)
C         err = ||SX*(x-TrueSolution)||/||SX*TrueSolution|| (2-Norms).
            ERR = DXNRM/SOLNRM
         ENDIF
      ENDIF
C         
      IF( IUNIT.NE.0 ) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL, MAXL, KMP
         ENDIF
         WRITE(IUNIT,1010) ITER, RNRM/BNRM, ERR
      ENDIF
      IF ( ERR.LE.TOL ) ISSGMR = 1
C         
      RETURN
 1000 FORMAT(' Generalized Minimum Residual(',i3,i3,') for ',
     $     'N, ITOL = ',I5, I5,
     $     /' ITER','   Natral Err Est','   Error Estimate')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7)
 1020 FORMAT(1x,' ITER = ',I5, ' IELMAX = ',I5,
     $     ' |R(IELMAX)/X(IELMAX)| = ',E12.5)
C------------- LAST LINE OF ISSGMR FOLLOWS ----------------------------
      END
*DECK SAXPY
      SUBROUTINE SAXPY (N, SA, SX, INCX, SY, INCY)
C***BEGIN PROLOGUE  SAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A7
C***TYPE      SINGLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SA  single precision scalar multiplier
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C       SY  single precision result (unchanged if N .LE. 0)
C
C     Overwrite single precision SY with single precision SA*SX +SY.
C     For I = 0 to N-1, replace  SY(LY+I*INCY) with SA*SX(LX+I*INCX) +
C       SY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SAXPY
      REAL SX(*), SY(*), SA
C***FIRST EXECUTABLE STATEMENT  SAXPY
      IF (N.LE.0 .OR. SA.EQ.0.0E0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I+1) = SY(I+1) + SA*SX(I+1)
        SY(I+2) = SY(I+2) + SA*SX(I+2)
        SY(I+3) = SY(I+3) + SA*SX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        SY(I) = SA*SX(I) + SY(I)
   70 CONTINUE
      RETURN
      END
*DECK SCOPY
      SUBROUTINE SCOPY (N, SX, INCX, SY, INCY)
C***BEGIN PROLOGUE  SCOPY
C***PURPOSE  Copy a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      SINGLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C       SY  copy of vector SX (unchanged if N .LE. 0)
C
C     Copy single precision SX to single precision SY.
C     For I = 0 to N-1, copy  SX(LX+I*INCX) to SY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SCOPY
      REAL SX(*), SY(*)
C***FIRST EXECUTABLE STATEMENT  SCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        SY(I) = SX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I+1) = SX(I+1)
        SY(I+2) = SX(I+2)
        SY(I+3) = SX(I+3)
        SY(I+4) = SX(I+4)
        SY(I+5) = SX(I+5)
        SY(I+6) = SX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        SY(I) = SX(I)
   70 CONTINUE
      RETURN
      END
*DECK SDOT
      REAL FUNCTION SDOT (N, SX, INCX, SY, INCY)
C***BEGIN PROLOGUE  SDOT
C***PURPOSE  Compute the inner product of two vectors.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A4
C***TYPE      SINGLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C     SDOT  single precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of single precision SX and SY.
C     SDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SDOT
      REAL SX(*), SY(*)
C***FIRST EXECUTABLE STATEMENT  SDOT
      SDOT = 0.0E0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SDOT = SDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        SDOT = SDOT + SX(I)*SY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      SDOT = SDOT + SX(I)*SY(I) + SX(I+1)*SY(I+1) + SX(I+2)*SY(I+2) +
     1              SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        SDOT = SDOT + SX(I)*SY(I)
   70 CONTINUE
      RETURN
      END
*DECK SBHIN
      SUBROUTINE SBHIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS,
     $     IUNIT, JOB )
C***BEGIN PROLOGUE  SBHIN
C***DATE WRITTEN   881107   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SBHIN-S),
C             Linear system, SLAP Sparse, Diagnostics
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Read a Sparse Linear System in the Boeing/Harwell Format.
C            The matrix is read in and if the right hand side is also
C            present in the input file then it too is read in.
C            The matrix is then modified to be in the SLAP Column
C            format.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB
C     REAL    A(NELT), SOLN(N), RHS(N)
C
C     CALL SBHIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB )
C
C *Arguments:
C N      :OUT      Integer
C         Order of the Matrix.
C NELT   :INOUT    Integer.
C         On input NELT is the maximum number of non-zeros that
C         can be stored in the IA, JA, A arrays.
C         On output NELT is the number of non-zeros stored in A.
C IA     :OUT      Integer IA(NELT).
C JA     :OUT      Integer JA(NELT).
C A      :OUT      Real A(NELT).
C         On output these arrays hold the matrix A in the SLAP
C         Triad format.  See "LONG DESCRIPTION", below.
C ISYM   :OUT      Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C SOLN   :OUT      Real SOLN(N).
C         The solution to the linear system, if present.  This array
C         is accessed if and only if JOB to read it in, see below.
C         If the user requests that SOLN be read in, but it is not in
C         the file, then it is simply zeroed out.
C RHS    :OUT      Real RHS(N).
C         The right hand side vector.  This array is accessed if and
C         only if JOB is set to read it in, see below.
C         If the user requests that RHS be read in, but it is not in 
C         the file, then it is simply zeroed out.
C IUNIT  :IN       Integer.
C         Fortran logical I/O device unit number to write the matrix
C         to.  This unit must be connected in a system dependent fashion
C         to a file or the console or you will get a nasty message
C         from the Fortran I/O libraries.
C JOB    :INOUT    Integer.
C         Flag indicating what I/O operations to perform.
C         On input JOB indicates what Input operations to try to 
C         perform.
C         JOB = 0 => Read only the matrix.
C             = 1 => Read matrix and RHS (if present).
C             = 2 => Read matrix and SOLN (if present).
C             = 3 => Read matrix, RHS and SOLN (if present).
C         On output JOB indicates what operations were actually
C         performed.
C               -3 => Unable to parse matrix "CODE" from input file
C                     to determine if only the lower triangle of matrix
C                     is stored.
C               -2 => Number of non-zeros (NELT) too large.
C               -1 => System size (N) too large.
C         JOB =  0 => Read in only the matrix.
C             =  1 => Read in the matrix and RHS.
C             =  2 => Read in the matrix and SOLN.
C             =  3 => Read in the matrix, RHS and SOLN.
C             = 10 => Read in only the matrix *STRUCTURE*, but no 
C                     non-zero entries.  Hence, A(*) is not referenced
C                     and has the return values the same as the input.
C             = 11 => Read in the matrix *STRUCTURE* and RHS.
C             = 12 => Read in the matrix *STRUCTURE* and SOLN.
C             = 13 => Read in the matrix *STRUCTURE*, RHS and SOLN.
C  
C *Precision:           Single Precision
C *Portability:
C         You must make sure that IUNIT is a valid Fortran logical
C         I/O device unit number and that the unit number has been
C         associated with a file or the console.  This is a system
C         dependent function.
C
C***LONG DESCRIPTION
C       The format for the output is as follows.  On  the first line
C       are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT
C       and ISYM are described above.  IRHS is  a flag indicating if
C       the RHS was  written out (1 is  yes, 0 is  no).  ISOLN  is a
C       flag indicating if the SOLN was written out  (1 is yes, 0 is
C       no).  The format for the fist line is: 5i10.  Then comes the
C       NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format
C       for  these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then  comes
C       RHS(I), I = 1, N, if IRHS = 1.  Then  comes SOLN(I), I  = 1,
C       N, if ISOLN = 1.  The format for these lines is: 1X,E16.7.
C
C       =================== S L A P Triad format ===================
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero  matrix  element is  placed   in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SBHIN
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB
      REAL    A(NELT), RHS(N), SOLN(N)
C
C         Local Variables
C
      CHARACTER*80  TITLE
      CHARACTER*3   CODE
      CHARACTER*16  PNTFMT, RINFMT
      CHARACTER*20  NVLFMT, RHSFMT
C
      INTEGER NLINE, NPLS, NRILS, NNVLS, NRHSLS, NROW, NCOL, NIND, NELE
C
C         Read Matrices In BOEING-HARWELL format.
C
C NLINE  Number of Data (after the header) lines in the file.
C NPLS   Number of lines for the Column Pointer data in the file.
C NRILS  Number of lines for the Row indicies in the data file.
C NNVLS  Number of lines for the Matrix elements in the data file.
C NRHSLS Number of lines for the RHS in the data file.
C
C***FIRST EXECUTABLE STATEMENT  SBHIN
      READ(IUNIT,9000) TITLE
      READ(IUNIT,9010) NLINE, NPLS, NRILS, NNVLS, NRHSLS
      READ(IUNIT,9020) CODE, NROW, NCOL, NIND, NELE
      READ(IUNIT,9030) PNTFMT, RINFMT, NVLFMT, RHSFMT
C
      IF( NROW.GT.N ) THEN
         N = NROW
         JOBRET = -1
         GOTO 999
      ENDIF
      IF( NIND.GT.NELT ) THEN
         NELT = NIND
         JOBRET = -2
         GOTO 999
      ENDIF
C
C         Set the parameters.
C
      N    = NROW
      NELT = NIND
      IF( CODE.EQ.'RUA' ) THEN
         ISYM = 0
      ELSE IF( CODE.EQ.'RSA' ) THEN
         ISYM = 1
      ELSE
         JOBRET = -3
         GOTO 999
      ENDIF
      READ(IUNIT,PNTFMT) (JA(I), I = 1, N+1)
      READ(IUNIT,RINFMT) (IA(I), I = 1, NELT)
      JOBRET = 10
      IF( NNVLS.GT.0 ) THEN
         READ(IUNIT,NVLFMT) (A(I),  I = 1, NELT)
         JOBRET = 0
      ENDIF
      IF( NRHSLS.GT.0 .AND. MOD(JOB,2).EQ.1 ) THEN
         READ(5,RHSFMT) (RHS(I), I = 1, N)
         JOBRET = JOBRET + 1
      ENDIF
C
C         Now loop thru the IA(i) array making sure that the Diagonal
C         matrix element appears first in the column.  Then sort the
C         rest of the column in ascending order.
C
      DO 70 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 30 I = IBGN, IEND
            IF( IA(I).EQ.ICOL ) THEN
C         Swap the diag element with the first element in the column.
               ITEMP = IA(I)
               IA(I) = IA(IBGN)
               IA(IBGN) = ITEMP
               TEMP = A(I)
               A(I) = A(IBGN)
               A(IBGN) = TEMP
               GOTO 40
            ENDIF
 30      CONTINUE
 40      IBGN = IBGN + 1
         IF( IBGN.LT.IEND ) THEN
            DO 60 I = IBGN, IEND
               DO 50 J = I+1, IEND
                  IF( IA(I).GT.IA(J) ) THEN
                     ITEMP = IA(I)
                     IA(I) = IA(J)
                     IA(J) = ITEMP
                     TEMP = A(I)
                     A(I) = A(J)
                     A(J) = TEMP
                  ENDIF
 50            CONTINUE
 60         CONTINUE
         ENDIF
 70   CONTINUE
C
C         Set return flag.
 999  JOB = JOBRET
      RETURN
 9000 FORMAT( A80 )
 9010 FORMAT( 5I14 )
 9020 FORMAT( A3, 11X, 4I14 )
 9030 FORMAT( 2A16, 2A20 )
C------------- LAST LINE OF SBHIN FOLLOWS ------------------------------
      END
*DECK SCHKW
      SUBROUTINE SCHKW( NAME, LOCIW, LENIW, LOCW, LENW,
     $     IERR, ITER, ERR )
C***BEGIN PROLOGUE  SCHKW
C***DATE WRITTEN   880225   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  R2
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SCHKW-S),
C             SLAP, Error Checking, Workspace Checking
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP WORK/IWORK Array Bounds Checker.
C            This routine checks the work array lengths  and  inter-
C            faces to the SLATEC  error  handler  if  a  problem  is 
C            found.
C***DESCRIPTION
C *Usage:
C     CHARACTER*(*) NAME
C     INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
C     REAL    ERR
C
C     CALL SCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
C
C *Arguments:
C NAME   :IN       Character*(*).
C         Name of the calling routine.  This is used in the output
C         message, if an error is detected.
C LOCIW  :IN       Integer.
C         Location of the first free element in the integer workspace
C         array.
C LENIW  :IN       Integer.
C         Length of the integer workspace array.
C LOCW   :IN       Integer.
C         Location of the first free element in the real workspace
C         array.
C LENRW  :IN       Integer.
C         Length of the real workspace array.
C IERR   :OUT      Integer.
C         Return error flag.
C               IERR = 0 => All went well.
C               IERR = 1 => Insufficient storage allocated for 
C                           WORK or IWORK.
C ITER   :OUT      Integer.
C         Set to 0 if an error is detected.
C ERR    :OUT      Real.
C         Set to a very large number if an error is detected.
C
C *Precision:           Single Precision
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH, XERRWV
C***END PROLOGUE  SCHKW
      CHARACTER*(*) NAME
      CHARACTER*72 MESG
      INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
      REAL    ERR, R1MACH
c      EXTERNAL R1MACH, XERRWV
      EXTERNAL R1MACH
C
C         Check the Integer workspace situation.
C***FIRST EXECUTABLE STATEMENT  SCHKW
      IERR = 0
      IF( LOCIW.GT.LENIW ) THEN
         IERR = 1
         ITER = 0
         ERR = R1MACH(2)
         MESG = NAME // ': INTEGER work array too short. '//
     $        ' IWORK needs i1: have allocated i2.'
c         CALL XERRWV( MESG, LEN(MESG), 1, 1, 2, LOCIW, LENIW,
c     $        0, 0.0, 0.0 )
      ENDIF
C
C         Check the Real workspace situation.
      IF( LOCW.GT.LENW ) THEN
         IERR = 1
         ITER = 0
         ERR = R1MACH(2)
         MESG = NAME // ': REAL work array too short. '//
     $        ' RWORK needs i1: have allocated i2.'
c         CALL XERRWV( MESG, LEN(MESG), 1, 1, 2, LOCW, LENW,
c     $        0, 0.0, 0.0 )
      ENDIF
      RETURN
C------------- LAST LINE OF SCHKW FOLLOWS ----------------------------
      END
*DECK QS2I1R
      SUBROUTINE QS2I1R( IA, JA, A, N, KFLAG )
C***BEGIN PROLOGUE  QS2I1R
C***DATE WRITTEN   761118   (YYMMDD)
C***REVISION DATE  890125   (YYMMDD)
C***CATEGORY NO.  N6A2A
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=INTEGER(QS2I1R-I),
C             QUICKSORT,SINGLETON QUICKSORT,SORT,SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Kahaner, D. K., (NBS)
C           Seager, M. K., (LLNL) seager@lll-crg.llnl.gov
C           Wisniewski, J. A., (SNLA)
C***PURPOSE  Sort an integer array also moving an integer and real array
C            This routine sorts the integer array IA and makes the same
C            interchanges in the integer array JA and the real array A.
C            The array IA may be sorted in increasing order or decreas-
C            ing order. A slightly modified QUICKSORT algorithm is used.
C
C***DESCRIPTION
C     Written by Rondall E Jones
C     Modified by John A. Wisniewski to use the Singleton QUICKSORT
C     algorithm. date 18 November 1976.
C
C     Further modified by David K. Kahaner
C     National Bureau of Standards
C     August, 1981
C
C     Even further modification made to bring the code up to the 
C     Fortran 77 level and make it more readable and to carry
C     along one integer array and one real array during the sort by
C     Mark K. Seager
C     Lawrence Livermore National Laboratory
C     November, 1987
C     This routine was adapted from the ISORT routine.
C
C     ABSTRACT
C         This routine sorts an integer array IA and makes the same
C         interchanges in the integer array JA and the real array A.  
C         The array a may be sorted in increasing order or decreasing 
C         order.  A slightly modified quicksort algorithm is used.
C
C     DESCRIPTION OF PARAMETERS
C        IA - Integer array of values to be sorted.
C        JA - Integer array to be carried along.
C         A - Real array to be carried along.
C         N - Number of values in integer array IA to be sorted.
C     KFLAG - Control parameter
C           = 1 means sort IA in INCREASING order.
C           =-1 means sort IA in DECREASING order.
C
C***REFERENCES
C     Singleton, R. C., Algorithm 347, "An Efficient Algorithm for 
C     Sorting with Minimal Storage", cacm, Vol. 12, No. 3, 1969, 
C     Pp. 185-187.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  QS2I1R
      DIMENSION IL(21),IU(21)
      INTEGER   IA(N),JA(N),IT,IIT,JT,JJT
      REAL      A(N), TA, TTA
C
C***FIRST EXECUTABLE STATEMENT  QS2I1R
      NN = N
      IF (NN.LT.1) THEN
c         CALL XERROR ( 'QS2I1R- the number of values to be sorted was no
c     $t POSITIVE.',59,1,1)
         RETURN
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = IABS(KFLAG)
      IF ( KK.NE.1 ) THEN
c         CALL XERROR ( 'QS2I1R- the sort control parameter, k, was not 1
c     $ OR -1.',55,2,1)
         RETURN
      ENDIF
C
C     Alter array IA to get decreasing order if needed.
C
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
C
C     Sort IA and carry JA and A along.
C     And now...Just a little black magic...
      M = 1
      I = 1
      J = NN
      R = .375
 210  IF( R.LE.0.5898437 ) THEN
         R = R + 3.90625E-2
      ELSE
         R = R-.21875
      ENDIF
 225  K = I
C
C     Select a central element of the array and save it in location 
C     it, jt, at.
C
      IJ = I + IFIX (FLOAT (J-I) *R)
      IT = IA(IJ)
      JT = JA(IJ)
      TA = A(IJ)
C
C     If first element of array is greater than it, interchange with it.
C
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(I)
         JA(I)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(I)
         A(I)   = TA
         TA     = A(IJ)
      ENDIF
      L=J
C                           
C     If last element of array is less than it, swap with it.
C
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(J)
         JA(J)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(J)
         A(J)   = TA
         TA     = A(IJ)
C
C     If first element of array is greater than it, swap with it.
C
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
            JA(IJ) = JA(I)
            JA(I)  = JT
            JT     = JA(IJ)
            A(IJ)  = A(I)
            A(I)   = TA
            TA     = A(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is 
C     smaller than it.
C
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
C
C     Find an element in the first half of the array which is 
C     greater than it.
C
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
C
C     Interchange these elements.
C
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         JJT   = JA(L)
         JA(L) = JA(K)
         JA(K) = JJT
         TTA   = A(L)
         A(L)  = A(K)
         A(K)  = TTA
         GOTO 240
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted.
C
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
C
C     Begin again on another portion of the unsorted array.
C                                  
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      JT = JA(I+1)
      TA =  A(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      JA(K+1) = JA(K)
      A(K+1)  =  A(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      JA(K+1) = JT
      A(K+1)  = TA
      GO TO 265
C
C     Clean up, if necessary.
C
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
C------------- LAST LINE OF QS2I1R FOLLOWS ----------------------------
      END
*DECK SS2Y
      SUBROUTINE SS2Y(N, NELT, IA, JA, A, ISYM )
C***BEGIN PROLOGUE  SS2Y
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SS2Y-S),
C             Linear system, SLAP Sparse
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP Triad to SLAP Column Format Converter.
C            Routine to convert from the SLAP Triad to SLAP Column
C            format.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     REAL    A(NELT)
C
C     CALL SS2Y( N, NELT, IA, JA, A, ISYM )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of non-zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Real A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "LONG 
C         DESCRIPTION", below.  If the SLAP Triad format is used
C         this format is translated to the SLAP Column format by
C         this routine.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C
C *Precision:           Single Precision
C
C***LONG DESCRIPTION
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures.  If the SLAP Triad format is give
C       as input then this routine transforms it into SLAP Column
C       format.  The way this routine tells which format is given as
C       input is to look at JA(N+1).  If JA(N+1) = NELT+1 then we
C       have the SLAP Column format.  If that equality does not hold
C       then it is assumed that the IA, JA, A arrays contain the 
C       SLAP Triad format.
C       
C       =================== S L A P Triad format ===================
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  QS2I1R
C***END PROLOGUE  SS2Y
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      REAL    A(NELT)
C
C         Check to see if the (IA,JA,A) arrays are in SLAP Column 
C         format.  If it's not then transform from SLAP Triad.
C***FIRST EXECUTABLE STATEMENT  SS2LT
      IF( JA(N+1).EQ.NELT+1 ) RETURN
C
C         Sort into ascending order by COLUMN (on the ja array).
C         This will line up the columns.
C
      CALL QS2I1R( JA, IA, A, NELT, 1 )
C         
C         Loop over each column to see where the column indicies change 
C         in the column index array ja.  This marks the beginning of the 
C         next column.
C         
      JA(1) = 1
      DO 20 ICOL = 1, N-1
         DO 10 J = JA(ICOL)+1, NELT
            IF( JA(J).NE.ICOL ) THEN
               JA(ICOL+1) = J
               GOTO 20
            ENDIF
 10      CONTINUE
 20   CONTINUE
      JA(N+1) = NELT+1
C         
C         Mark the n+2 element so that future calls to a SLAP routine 
C         utilizing the YSMP-Column storage format will be able to tell.
C         
      JA(N+2) = 0
C
C         Now loop thru the ia(i) array making sure that the Diagonal
C         matrix element appears first in the column.  Then sort the
C         rest of the column in ascending order.
C
      DO 70 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 30 I = IBGN, IEND
            IF( IA(I).EQ.ICOL ) THEN
C         Swap the diag element with the first element in the column.
               ITEMP = IA(I)
               IA(I) = IA(IBGN)
               IA(IBGN) = ITEMP
               TEMP = A(I)
               A(I) = A(IBGN)
               A(IBGN) = TEMP
               GOTO 40
            ENDIF
 30      CONTINUE
 40      IBGN = IBGN + 1
         IF( IBGN.LT.IEND ) THEN
            DO 60 I = IBGN, IEND
               DO 50 J = I+1, IEND
                  IF( IA(I).GT.IA(J) ) THEN
                     ITEMP = IA(I)
                     IA(I) = IA(J)
                     IA(J) = ITEMP
                     TEMP = A(I)
                     A(I) = A(J)
                     A(J) = TEMP
                  ENDIF
 50            CONTINUE
 60         CONTINUE
         ENDIF
 70   CONTINUE
      RETURN
C------------- LAST LINE OF SS2Y FOLLOWS ----------------------------
      END
*DECK SCPPLT
      SUBROUTINE SCPPLT( N, NELT, IA, JA, A, ISYM, IUNIT )
C***BEGIN PROLOGUE  SCPPLT
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SCPPLT-S),
C             Linear system, SLAP Sparse, Diagnostics
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Printer Plot of SLAP Column Format Matrix.
C            Routine to print out a SLAP Column format matrix in
C            a "printer plot" graphical representation.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(N+1), ISYM, IUNIT
C     REAL    A(NELT)
C
C     CALL SCPPLT( N, NELT, IA, JA, A, ISYM, IUNIT )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of non-zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(N+1).
C A      :INOUT    Real A(NELT).
C         These arrays should hold the matrix A in the SLAP
C         Column format.  See "LONG DESCRIPTION", below.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C IUNIT  :IN       Integer.
C         Fortran logical I/O device unit number to write the matrix
C         to.  This unit must be connected in a system dependent fashion
C         to a file or the console or you will get a nasty message
C         from the Fortran I/O libraries.
C
C *Precision:           Single Precision
C *Portability:
C         You must make sure that IUNIT is a valid Fortran logical
C         I/O device unit number and that the unit number has been
C         associated with a file or the console.  This is a system
C         dependent function.
C
C***LONG DESCRIPTION
C       This routine prints out a SLAP  Column format matrix  to the
C       Fortran logical I/O unit   number  IUNIT.  The  numbers them
C       selves  are not printed  out, but   rather  a one  character
C       representation of the numbers.   Elements of the matrix that
C       are not represented in the (IA,JA,A)  arrays are  denoted by
C       ' ' character (a blank).  Elements of A that are *ZERO* (and
C       hence  should  really not be  stored) are  denoted  by a '0'
C       character.  Elements of A that are *POSITIVE* are denoted by
C       'D' if they are Diagonal elements  and '#' if  they are off
C       Diagonal  elements.  Elements of  A that are *NEGATIVE* are
C       denoted by 'N'  if they  are Diagonal  elements and  '*' if
C       they are off Diagonal elements.
C       
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SCPPLT
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      REAL    A(NELT)
      CHARACTER*225 CHMAT(225)
C
C         Set up the character matrix...
C***FIRST EXECUTABLE STATEMENT  SCPPLT
      NMAX = MIN( 225, N)
      DO 10 I = 1, NMAX
         CHMAT(I)(1:NMAX) = ' '
 10   CONTINUE
      DO 30 ICOL = 1, NMAX
         JBGN = JA(ICOL)
         JEND = JA(ICOL+1)-1
         DO 20 J = JBGN, JEND
            IROW = IA(J)
            IF( IROW.LE.NMAX ) THEN
               IF( ISYM.NE.0 ) THEN
C         Put in non-sym part as well...
                  IF( A(J).EQ.0.0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = '0'
                  ELSEIF( A(J).GT.0.0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = '#'
                  ELSE
                     CHMAT(IROW)(ICOL:ICOL) = '*'
                  ENDIF
               ENDIF
               IF( IROW.EQ.ICOL ) THEN
C         Diagonal entry.
                  IF( A(J).EQ.0.0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = '0'
                  ELSEIF( A(J).GT.0.0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = 'D'
                  ELSE
                     CHMAT(IROW)(ICOL:ICOL) = 'N'
                  ENDIF
               ELSE
C         Off-Diagonal entry
                  IF( A(J).EQ.0.0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = '0'
                  ELSEIF( A(J).GT.0.0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = '#'
                  ELSE
                     CHMAT(IROW)(ICOL:ICOL) = '*'
                  ENDIF
               ENDIF
            ENDIF
 20      CONTINUE
 30   CONTINUE
C
C         Write out the heading.
      WRITE(IUNIT,1000) N, NELT, FLOAT(NELT)/FLOAT(N*N)
      WRITE(IUNIT,1010) (MOD(I,10),I=1,NMAX)
C
C         Write out the character representations matrix elements.
      DO 40 IROW = 1, NMAX
         WRITE(IUNIT,1020) IROW, CHMAT(IROW)(1:NMAX)
 40   CONTINUE
      RETURN
 1000 FORMAT(/'**** Picture of Column SLAP matrix follows ****'/
     $     ' N, NELT and Density = ',2I10,E16.7)
 1010 FORMAT(4X,255(I1))
 1020 FORMAT(1X,I3,A)
C------------- LAST LINE OF SCPPLT FOLLOWS ----------------------------
      END
*DECK STOUT
      SUBROUTINE STOUT( N, NELT, IA, JA, A, ISYM, SOLN, RHS,
     $     IUNIT, JOB )
C***BEGIN PROLOGUE  STOUT
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(STOUT-S),
C             Linear system, SLAP Sparse, Diagnostics
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Write out SLAP Triad Format Linear System.
C            Routine to write out a SLAP Triad format matrix and
C            right hand side and solution to the system, if known.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB
C     REAL    A(NELT), SOLN(N), RHS(N)
C
C     CALL STOUT( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of non-zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Real A(NELT).
C         These arrays should hold the matrix A in the SLAP
C         Triad format.  See "LONG DESCRIPTION", below.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C SOLN   :IN       Real SOLN(N).
C         The solution to the linear system, if known.  This array
C         is accessed if and only if JOB is set to print it out, 
C         see below.
C RHS    :IN       Real RHS(N).
C         The right hand side vector.  This array is accessed if and
C         only if JOB is set to print it out, see below.
C IUNIT  :IN       Integer.
C         Fortran logical I/O device unit number to write the matrix
C         to.  This unit must be connected in a system dependent fashion
C         to a file or the console or you will get a nasty message
C         from the Fortran I/O libraries.
C JOB    :IN       Integer.
C         Flag indicating what I/O operations to perform.
C         JOB = 0 => Print only the matrix.
C             = 1 => Print matrix and RHS.
C             = 2 => Print matrix and SOLN.
C             = 3 => Print matrix, RHS and SOLN.
C
C *Precision:           Single Precision
C *Portability:
C         You must make sure that IUNIT is a valid Fortran logical
C         I/O device unit number and that the unit number has been
C         associated with a file or the console.  This is a system
C         dependent function.
C
C***LONG DESCRIPTION
C       The format for the output is as follows.  On  the first line
C       are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT
C       and ISYM are described above.  IRHS is  a flag indicating if
C       the RHS was  written out (1 is  yes, 0 is  no).  ISOLN  is a
C       flag indicating if the SOLN was written out  (1 is yes, 0 is
C       no).  The format for the fist line is: 5i10.  Then comes the
C       NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format
C       for  these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then  comes
C       RHS(I), I = 1, N, if IRHS = 1.  Then  comes SOLN(I), I  = 1,
C       N, if ISOLN = 1.  The format for these lines is: 1X,E16.7.
C
C       =================== S L A P Triad format ===================
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero  matrix  element is  placed   in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  STOUT
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB
      REAL    A(NELT), RHS(N), SOLN(N)
C
C         Local variables.
C
      INTEGER IRHS, ISOLN, I
C
C         If RHS and SOLN are to be printed also.
C         Write out the information heading.
C***FIRST EXECUTABLE STATEMENT  STOUT
      IRHS = 0
      ISOLN = 0
      IF( JOB.EQ.1 .OR. JOB.EQ.3 ) IRHS = 1
      IF( JOB.GT.1 ) ISOLN = 1
      WRITE(IUNIT,1000) N, NELT, ISYM, IRHS, ISOLN
C
C         Write out the matrix non-zeros in Triad format.
      DO 10 I = 1, NELT
         WRITE(IUNIT,1010) IA(I), JA(I), A(I)
 10   CONTINUE
C
C         If requested, write out the rhs.
      IF( IRHS.EQ.1 ) THEN
         WRITE(IUNIT,1020) (RHS(I),I=1,N)
      ENDIF
C
C         If requested, write out the soln.
      IF( ISOLN.EQ.1 ) THEN
         WRITE(IUNIT,1020) (SOLN(I),I=1,N)
      ENDIF
      RETURN
 1000 FORMAT(5I10)
 1010 FORMAT(1X,I5,1X,I5,1X,E16.7)
 1020 FORMAT(1X,E16.7)
C------------- LAST LINE OF STOUT FOLLOWS ----------------------------
      END
*DECK STIN
      SUBROUTINE STIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS,
     $     IUNIT, JOB )
C***BEGIN PROLOGUE  STIN
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(STIN-S),
C             Linear system, SLAP Sparse, Diagnostics
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Read in SLAP Triad Format Linear System.
C            Routine to read in a SLAP Triad format matrix and
C            right hand side and solution to the system, if known.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB
C     REAL    A(NELT), SOLN(N), RHS(N)
C
C     CALL STIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB )
C
C *Arguments:
C N      :OUT      Integer
C         Order of the Matrix.
C NELT   :INOUT    Integer.
C         On input NELT is the maximum number of non-zeros that
C         can be stored in the IA, JA, A arrays.
C         On output NELT is the number of non-zeros stored in A.
C IA     :OUT      Integer IA(NELT).
C JA     :OUT      Integer JA(NELT).
C A      :OUT      Real A(NELT).
C         On output these arrays hold the matrix A in the SLAP
C         Triad format.  See "LONG DESCRIPTION", below.
C ISYM   :OUT      Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C SOLN   :OUT      Real SOLN(N).
C         The solution to the linear system, if present.  This array
C         is accessed if and only if JOB to read it in, see below.
C         If the user requests that SOLN be read in, but it is not in
C         the file, then it is simply zeroed out.
C RHS    :OUT      Real RHS(N).
C         The right hand side vector.  This array is accessed if and
C         only if JOB is set to read it in, see below.
C         If the user requests that RHS be read in, but it is not in 
C         the file, then it is simply zeroed out.
C IUNIT  :IN       Integer.
C         Fortran logical I/O device unit number to write the matrix
C         to.  This unit must be connected in a system dependent fashion
C         to a file or the console or you will get a nasty message
C         from the Fortran I/O libraries.
C JOB    :INOUT    Integer.
C         Flag indicating what I/O operations to perform.
C         On input JOB indicates what Input operations to try to 
C         perform.
C         JOB = 0 => Read only the matrix.
C             = 1 => Read matrix and RHS (if present).
C             = 2 => Read matrix and SOLN (if present).
C             = 3 => Read matrix, RHS and SOLN (if present).
C         On output JOB indicates what operations were actually
C         performed.
C         JOB = 0 => Read in only the matrix.
C             = 1 => Read in the matrix and RHS.
C             = 2 => Read in the matrix and SOLN.
C             = 3 => Read in the matrix, RHS and SOLN.
C
C *Precision:           Single Precision
C *Portability:
C         You must make sure that IUNIT is a valid Fortran logical
C         I/O device unit number and that the unit number has been
C         associated with a file or the console.  This is a system
C         dependent function.
C
C***LONG DESCRIPTION
C       The format for the output is as follows.  On  the first line
C       are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT
C       and ISYM are described above.  IRHS is  a flag indicating if
C       the RHS was  written out (1 is  yes, 0 is  no).  ISOLN  is a
C       flag indicating if the SOLN was written out  (1 is yes, 0 is
C       no).  The format for the fist line is: 5i10.  Then comes the
C       NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format
C       for  these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then  comes
C       RHS(I), I = 1, N, if IRHS = 1.  Then  comes SOLN(I), I  = 1,
C       N, if ISOLN = 1.  The format for these lines is: 1X,E16.7.
C
C       =================== S L A P Triad format ===================
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero  matrix  element is  placed   in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  STIN
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB
      REAL    A(NELT), RHS(N), SOLN(N)
C
C         Local variables.
C
      INTEGER IRHS, ISOLN, I, NELTMAX
C
C         Read in the information heading.
C***FIRST EXECUTABLE STATEMENT  STIN
      NELTMAX = NELT
      READ(IUNIT,1000) N, NELT, ISYM, IRHS, ISOLN
      NELT = MIN( NELT, NELTMAX )
C
C         Read in the matrix non-zeros in Triad format.
      DO 10 I = 1, NELT
         READ(IUNIT,1010) IA(I), JA(I), A(I)
 10   CONTINUE
C
C         If requested, read in the rhs.
      JOBRET = 0
      IF( JOB.EQ.1 .OR. JOB.EQ.3 ) THEN
C
C         Check to see if rhs is in the file.
         IF( IRHS.EQ.1 ) THEN
            JOBRET = 1
            READ(IUNIT,1020) (RHS(I),I=1,N)
         ELSE
            DO 20 I = 1, N
               RHS(I) = 0.0
 20         CONTINUE
         ENDIF
      ENDIF
C
C         If requested, read in the soln.
      IF( JOB.GT.1 ) THEN
C
C         Check to see if soln is in the file.
         IF( ISOLN.EQ.1 ) THEN
            JOBRET = JOBRET + 2
            READ(IUNIT,1020) (SOLN(I),I=1,N)
         ELSE
            DO 30 I = 1, N
               SOLN(I) = 0.0
 30         CONTINUE
         ENDIF
      ENDIF
C
      JOB = JOBRET
      RETURN
 1000 FORMAT(5I10)
 1010 FORMAT(1X,I5,1X,I5,1X,E16.7)
 1020 FORMAT(1X,E16.7)
C------------- LAST LINE OF STIN FOLLOWS ----------------------------
      END
*DECK SSDS
      SUBROUTINE SSDS(N, NELT, IA, JA, A, ISYM, DINV)
C***BEGIN PROLOGUE  SSDS
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSDS-S),
C             SLAP Sparse, Diagonal
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Diagonal Scaling Preconditioner SLAP Set Up.
C            Routine to compute the inverse of the diagonal of a matrix
C            stored in the SLAP Column format.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     REAL    A(NELT), DINV(N)
C
C     CALL SSDS( N, NELT, IA, JA, A, ISYM, DINV )
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Real A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C DINV   :OUT      Real DINV(N).
C         Upon return this array holds 1./DIAG(A).
C
C *Description
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       With the SLAP  format  all  of  the   "inner  loops" of this
C       routine should vectorize  on  machines with hardware support
C       for vector   gather/scatter  operations.  Your compiler  may
C       require a compiler directive to  convince it that  there are
C       no  implicit  vector  dependencies.  Compiler directives for
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
C       supplied with the standard SLAP distribution.
C
C *Precision:           Single Precision
C
C *Cautions:
C       This routine assumes that the diagonal of A is all  non-zero
C       and that the operation DINV = 1.0/DIAG(A) will not underflow
C       or overflow.    This  is done so that the  loop  vectorizes.
C       Matricies with zero or near zero or very  large entries will
C       have numerical difficulties  and  must  be fixed before this 
C       routine is called.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSDS
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      REAL    A(NELT), DINV(N)
C
C         Assume the Diagonal elements are the first in each column.
C         This loop should *VECTORIZE*.  If it does not you may have
C         to add a compiler directive.  We do not check for a zero
C         (or near zero) diagonal element since this would interfere 
C         with vectorization.  If this makes you nervous put a check
C         in!  It will run much slower.
C***FIRST EXECUTABLE STATEMENT  SSDS
 1    CONTINUE
      DO 10 ICOL = 1, N
         DINV(ICOL) = 1.0/A(JA(ICOL))
 10   CONTINUE
C         
      RETURN
C------------- LAST LINE OF SSDS FOLLOWS ----------------------------
      END
*DECK SSDSCL
      SUBROUTINE SSDSCL( N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB,
     $     ITOL )
C***BEGIN PROLOGUE  SSDSCL
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSDSCL-S),
C             SLAP Sparse, Diagonal
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Diagonal Scaling of system Ax = b.
C            This routine scales (and unscales) the system Ax = b 
C            by symmetric diagonal scaling.  The new system is:
C             -1/2  -1/2  1/2      -1/2
C            D    AD    (D   x) = D    b
C            when scaling is selected with the JOB parameter.  When
C            unscaling is selected this process is reversed.
C            The true solution is also scaled or unscaled if ITOL is set
C            appropriately, see below.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB, ITOL
C     REAL    A(NELT), DINV(N)
C
C     CALL SSDSCL( N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB, ITOL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Real A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C X      :INOUT    Real X(N).
C         Initial guess that will be later used in the iterative 
C         solution.
C         of the scaled system.
C B      :INOUT    Real B(N).
C         Right hand side vector.
C DINV   :OUT      Real DINV(N).
C         Upon return this array holds 1./DIAG(A).
C JOB    :IN       Integer.
C         Flag indicating weather to scale or not.  JOB nonzero means
C         do scaling.  JOB = 0 means do unscaling.
C ITOL   :IN       Integer.
C         Flag indicating what type of error estimation to do in the 
C         iterative method.  When ITOL = 11 the exact solution from
C         common block solblk will be used.  When the system is scaled
C         then the true solution must also be scaled.  If ITOL is not
C         11 then this vector is not referenced.
C
C *Common Blocks:
C SOLN    :INOUT   Real SOLN(N).  COMMON BLOCK /SOLBLK/
C         The true solution, SOLN, is scaled (or unscaled) if ITOL is
C         set to 11, see above.
C
C *Description
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       With the SLAP  format  all  of  the   "inner  loops" of this
C       routine should vectorize  on  machines with hardware support
C       for vector   gather/scatter  operations.  Your compiler  may
C       require a compiler directive to  convince it that  there are
C       no  implicit  vector  dependencies.  Compiler directives for
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
C       supplied with the standard SLAP distribution.
C
C *Precision:           Single Precision
C
C *Cautions:
C       This routine assumes that the diagonal of A is all  non-zero
C       and that the operation DINV = 1.0/DIAG(A)  will  not  under-
C       flow or overflow. This is done so that the loop  vectorizes.
C       Matricies with zero or near zero or very  large entries will
C       have numerical difficulties  and  must  be fixed before this 
C       routine is called.
C
C *See Also:
C       SSDCG
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSDSCL
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB, ITOL
      REAL    A(NELT), X(N), B(N), DINV(N)
      COMMON /SOLBLK/ SOLN(1)
C
C         SCALING...
C
      IF( JOB.NE.0 ) THEN
         DO 10 ICOL = 1, N
            DINV(ICOL) = 1.0/SQRT( A(JA(ICOL)) )
 10      CONTINUE
      ELSE
C
C         UNSCALING...
C
         DO 15 ICOL = 1, N
            DINV(ICOL) = 1.0/DINV(ICOL)
 15      CONTINUE
      ENDIF
C
      DO 30 ICOL = 1, N
         JBGN = JA(ICOL)
         JEND = JA(ICOL+1)-1
         DI = DINV(ICOL)
         DO 20 J = JBGN, JEND
            A(J) = DINV(IA(J))*A(J)*DI
 20      CONTINUE
 30   CONTINUE
C
      DO 40 ICOL = 1, N
         B(ICOL) = B(ICOL)*DINV(ICOL)
         X(ICOL) = X(ICOL)/DINV(ICOL)
 40   CONTINUE
C
C         Check to see if we need to scale the "true solution" as well.
C
      IF( ITOL.EQ.11 ) THEN
         DO 50 ICOL = 1, N
            SOLN(ICOL) = SOLN(ICOL)/DINV(ICOL)
 50      CONTINUE
      ENDIF
C
      RETURN
      END
*DECK SSD2S
      SUBROUTINE SSD2S(N, NELT, IA, JA, A, ISYM, DINV)
C***BEGIN PROLOGUE  SSD2S
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSD2S-S),
C             SLAP Sparse, Diagonal
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Diagonal Scaling Preconditioner SLAP Normal Eqns Set Up.
C            Routine to compute the inverse of the diagonal of the 
C            matrix A*A'.  Where A is stored in SLAP-Column format.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     REAL    A(NELT), DINV(N)
C
C     CALL SSD2S( N, NELT, IA, JA, A, ISYM, DINV )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Real A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C DINV   :OUT      Real DINV(N).
C         Upon return this array holds 1./DIAG(A*A').
C
C *Description
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       With the SLAP  format  all  of  the   "inner  loops" of this
C       routine should vectorize  on  machines with hardware support
C       for vector   gather/scatter  operations.  Your compiler  may
C       require a compiler directive to  convince it that  there are
C       no  implicit  vector  dependencies.  Compiler directives for
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
C       supplied with the standard SLAP distribution.
C
C *Precision:           Single Precision
C
C *Cautions:
C       This routine assumes that the diagonal of A is all  non-zero
C       and that the operation DINV = 1.0/DIAG(A*A') will not under-
C       flow or overflow. This is done so that the loop  vectorizes.
C       Matricies with zero or near zero or very  large entries will
C       have numerical difficulties  and  must  be fixed before this 
C       routine is called.
C
C *See Also:
C       SSDCGN
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSD2S
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      REAL    A(NELT), DINV(N)
C
C***FIRST EXECUTABLE STATEMENT  SSD2S
      DO 10 I = 1, N
         DINV(I) = 0.
 10   CONTINUE
C
C         Loop over each column.
      DO 40 I = 1, N
         KBGN = JA(I)
         KEND = JA(I+1) - 1
C
C         Add in the contributions for each row that has a non-zero
C         in this column.
         DO 20 K = KBGN, KEND
            DINV(IA(K)) = DINV(IA(K)) + A(K)**2
 20      CONTINUE
         IF( ISYM.EQ.1 ) THEN
C
C         Lower triangle stored by columns => upper triangle stored by
C         rows with Diagonal being the first entry.  Loop across the
C         rest of the row.
            KBGN = KBGN + 1
            IF( KBGN.LE.KEND ) THEN
               DO 30 K = KBGN, KEND
                  DINV(I) = DINV(I) + A(K)**2
 30            CONTINUE
            ENDIF
         ENDIF
 40   CONTINUE
      DO 50 I=1,N
         DINV(I) = 1./DINV(I)
 50   CONTINUE
C         
      RETURN
C------------- LAST LINE OF SSD2S FOLLOWS ----------------------------
      END
*DECK SS2LT
      SUBROUTINE SS2LT( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL )
C***BEGIN PROLOGUE  SS2LT
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SS2LT-S),
C             Linear system, SLAP Sparse, Lower Triangle
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Lower Triangle Preconditioner SLAP Set Up.
C            Routine to store the lower triangle of a matrix stored
C            in the Slap Column format.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     INTEGER NEL, IEL(N+1), JEL(NEL), NROW(N)
C     REAL    A(NELT), EL(NEL)
C
C     CALL SS2LT( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of non-zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Real A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C NEL    :OUT      Integer.
C         Number of non-zeros in the lower triangle of A.   Also 
C         coresponds to the length of the JEL, EL arrays.
C IEL    :OUT      Integer IEL(N+1).
C JEL    :OUT      Integer JEL(NEL).
C EL     :OUT      Real     EL(NEL).
C         IEL, JEL, EL contain the lower triangle of the A matrix
C         stored in SLAP Column format.  See "Description", below
C         for more details bout the SLAP Column format.
C
C *Description
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C *Precision:           Single Precision
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SS2LT
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM 
      INTEGER NEL, IEL(NEL), JEL(NEL)
      REAL    A(NELT), EL(NELT)
C***FIRST EXECUTABLE STATEMENT  SS2LT
      IF( ISYM.EQ.0 ) THEN
C
C         The matrix is stored non-symmetricly.  Pick out the lower
C         triangle.
C
         NEL = 0
         DO 20 ICOL = 1, N
            JEL(ICOL) = NEL+1
            JBGN = JA(ICOL)
            JEND = JA(ICOL+1)-1
            DO 10 J = JBGN, JEND
               IF( IA(J).GE.ICOL ) THEN
                  NEL = NEL + 1
                  IEL(NEL) = IA(J)
                  EL(NEL)  = A(J)
               ENDIF
 10         CONTINUE
 20      CONTINUE
         JEL(N+1) = NEL+1
      ELSE
C
C         The matrix is symmetric and only the lower triangle is 
C         stored.  Copy it to IEL, JEL, EL.
C
         NEL = NELT
         DO 30 I = 1, NELT
            IEL(I) = IA(I)
            EL(I) = A(I)
 30      CONTINUE
         DO 40 I = 1, N+1
            JEL(I) = JA(I)
 40      CONTINUE
      ENDIF
      RETURN
C------------- LAST LINE OF SS2LT FOLLOWS ----------------------------
      END
*DECK SSICS
      SUBROUTINE SSICS(N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL,
     $     EL, D, R, IWARN )
C***BEGIN PROLOGUE  SSICS
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSICS-S),
C             Linear system, SLAP Sparse, Iterative Precondition
C             Incomplete Cholesky Factorization.
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Incompl Cholesky Decomposition Preconditioner SLAP Set Up.
C            Routine to generate the Incomplete Cholesky decomposition,
C            L*D*L-trans, of  a symmetric positive definite  matrix, A,
C            which  is stored  in  SLAP Column format.  The  unit lower
C            triangular matrix L is  stored by rows, and the inverse of
C            the diagonal matrix D is stored.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     INTEGER NEL, IEL(NEL), JEL(N+1), IWARN
C     REAL    A(NELT), EL(NEL), D(N), R(N)
C
C     CALL SSICS( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, R,
C    $    IWARN )
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Real A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower 
C         triangle of the matrix is stored.
C NEL    :OUT      Integer.
C         Number of non-zeros in the lower triangle of A.   Also 
C         coresponds to the length of the JEL, EL arrays.
C IEL    :OUT      Integer IEL(N+1).
C JEL    :OUT      Integer JEL(NEL).
C EL     :OUT      Real     EL(NEL).
C         IEL, JEL, EL contain the unit lower triangular factor  of the
C         incomplete decomposition   of the A  matrix  stored  in  SLAP
C         Row format.   The Diagonal of   ones   *IS*   stored.     See 
C         "Description", below for more details about the SLAP Row fmt.  
C D      :OUT      Real D(N)
C         Upon return this array holds D(I) = 1./DIAG(A).
C R      :WORK     Real R(N).
C         Temporary real workspace needed for the factorization.
C IWARN  :OUT      Integer.
C         This is a warning variable and is zero if the IC factoriza-
C         tion goes well.  It is set to the row index corresponding to
C         the last zero pivot found.  See "Description", below.
C
C *Description
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       ==================== S L A P Row format ====================
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the real
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C       
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|  
C       |51  0 53  0 55|  
C
C       With the SLAP  format some  of  the   "inner  loops" of this
C       routine should vectorize  on  machines with hardware support
C       for vector   gather/scatter  operations.  Your compiler  may
C       require a compiler directive to  convince it that  there are
C       no  implicit  vector  dependencies.  Compiler directives for
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
C       supplied with the standard SLAP distribution.
C
C       The IC  factorization is not  alway exist for SPD matricies.
C       In the event that a zero pivot is found it is set  to be 1.0
C       and the factorization procedes.   The integer variable IWARN
C       is set to the last row where the Diagonal was fudged.  This 
C       eventuality hardly ever occurs in practice
C       
C *Precision:           Single Precision
C
C *See Also:
C       SCG, SSICCG
C***REFERENCES  1. Gene Golub & Charles Van Loan, "Matrix Computations",
C                 John Hopkins University Press; 3 (1983) IBSN 
C                 0-8018-3010-9.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  SSICS
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      INTEGER NEL, IEL(NEL), JEL(NEL)
      REAL    A(NELT), EL(NEL), D(N), R(N)
C         
C         Set the lower triangle in IEL, JEL, EL
C***FIRST EXECUTABLE STATEMENT  SSICS
      IWARN = 0
C         
C         All matrix elements stored in IA, JA, A.  Pick out the lower 
C         triangle (making sure that the Diagonal of EL is one) and
C         store by rows.
C         
      NEL = 1
      IEL(1) = 1
      JEL(1) = 1
      EL(1) = 1.0
      D(1) = A(1)
      DO 30 IROW = 2, N
C         Put in the Diagonal.
         NEL = NEL + 1
         IEL(IROW) = NEL
         JEL(NEL) = IROW
         EL(NEL) = 1.0
         D(IROW) = A(JA(IROW))
C
C         Look in all the lower triangle columns for a matching row.
C         Since the matrix is symmetric, we can look across the 
C         irow-th row by looking down the irow-th column (if it is 
C         stored ISYM=0)...
         IF( ISYM.EQ.0 ) THEN
            ICBGN = JA(IROW)
            ICEND = JA(IROW+1)-1
         ELSE
            ICBGN = 1
            ICEND = IROW-1
         ENDIF
         DO 20 IC = ICBGN, ICEND
            IF( ISYM.EQ.0 ) THEN
               ICOL = IA(IC)
               IF( ICOL.GE.IROW ) GOTO 20
            ELSE
               ICOL = IC
            ENDIF
            JBGN = JA(ICOL)+1
            JEND = JA(ICOL+1)-1
            IF( JBGN.LE.JEND .AND. IA(JEND).GE.IROW ) THEN
               DO 10 J = JBGN, JEND
                  IF( IA(J).EQ.IROW ) THEN
                     NEL = NEL + 1
                     JEL(NEL) = ICOL
                     EL(NEL)  = A(J)
                     GOTO 20
                  ENDIF
 10            CONTINUE
            ENDIF
 20      CONTINUE
 30   CONTINUE
      IEL(N+1) = NEL+1
C
C         Sort ROWS of lower triangle into descending order (count out 
C         along rows out from Diagonal).
C
      DO 60 IROW = 2, N
         IBGN = IEL(IROW)+1
         IEND = IEL(IROW+1)-1
         IF( IBGN.LT.IEND ) THEN
            DO 50 I = IBGN, IEND-1
               DO 40 J = I+1, IEND
                  IF( JEL(I).GT.JEL(J) ) THEN
                     JELTMP = JEL(J)
                     JEL(J) = JEL(I)
                     JEL(I) = JELTMP
                     ELTMP = EL(J)
                     EL(J) = EL(I)
                     EL(I) = ELTMP
                  ENDIF
 40            CONTINUE
 50         CONTINUE
         ENDIF
 60   CONTINUE
C         
C         Perform the Incomplete Cholesky decomposition by looping 
C         over the rows.
C         Scale the first column.  Use the structure of A to pick out
C         the rows with something in column 1.
C
      IRBGN = JA(1)+1
      IREND = JA(2)-1
      DO 65 IRR = IRBGN, IREND
         IR = IA(IRR)
C         Find the index into EL for EL(1,IR).  
C         Hint: it's the second entry.
         I = IEL(IR)+1
         EL(I) = EL(I)/D(1)
 65   CONTINUE
C
      DO 110 IROW = 2, N
C
C         Update the IROW-th diagonal.
C
         DO 66 I = 1, IROW-1
            R(I) = 0.0
 66      CONTINUE
         IBGN = IEL(IROW)+1
         IEND = IEL(IROW+1)-1
         IF( IBGN.LE.IEND ) THEN
            DO 70 I = IBGN, IEND
               R(JEL(I)) = EL(I)*D(JEL(I))
               D(IROW) = D(IROW) - EL(I)*R(JEL(I))
 70         CONTINUE
C
C         Check to see if we gota problem with the diagonal.
C
            IF( D(IROW).LE.0.0 ) THEN
               IF( IWARN.EQ.0 ) IWARN = IROW
               D(IROW) = 1.0
            ENDIF
         ENDIF
C
C         Update each EL(IROW+1:N,IROW), if there are any.
C         Use the structure of A to determine the Non-zero elements
C         of the IROW-th column of EL.
C
         IRBGN = JA(IROW)
         IREND = JA(IROW+1)-1
         DO 100 IRR = IRBGN, IREND
            IR = IA(IRR)
            IF( IR.LE.IROW ) GOTO 100
C         Find the index into EL for EL(IR,IROW)
            IBGN = IEL(IR)+1
            IEND = IEL(IR+1)-1
            IF( JEL(IBGN).GT.IROW ) GOTO 100
            DO 90 I = IBGN, IEND
               IF( JEL(I).EQ.IROW ) THEN
                  ICEND = IEND
 91               IF( JEL(ICEND).GE.IROW ) THEN
                     ICEND = ICEND - 1
                     GOTO 91
                  ENDIF
C         Sum up the EL(IR,1:IROW-1)*R(1:IROW-1) contributions.
                  DO 80 IC = IBGN, ICEND
                     EL(I) = EL(I) - EL(IC)*R(JEL(IC))
 80               CONTINUE
                  EL(I) = EL(I)/D(IROW)
                  GOTO 100
               ENDIF
 90         CONTINUE
C
C         If we get here, we have real problems...
c            CALL XERRWV('SSICS -- A and EL data structure mismatch'//
c     $           ' in row (i1)',53,1,2,1,IROW,0,0,0.0,0.0)
 100     CONTINUE
 110  CONTINUE
C        
C         Replace diagonals by their inverses.
C
      DO 120 I =1, N
         D(I) = 1.0/D(I)
 120  CONTINUE
      RETURN
C------------- LAST LINE OF SSICS FOLLOWS ----------------------------
      END
*DECK SSILUS
      SUBROUTINE SSILUS(N, NELT, IA, JA, A, ISYM, NL, IL, JL,
     $     L, DINV, NU, IU, JU, U, NROW, NCOL)
C***BEGIN PROLOGUE  SSILUS
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSILUS-S),
C             Non-Symmetric Linear system, Sparse, 
C             Iterative Precondition, Incomplete LU Factorization
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up.
C            Routine to generate the incomplete LDU decomposition of a 
C            matrix.  The  unit lower triangular factor L is stored by 
C            rows and the  unit upper triangular factor U is stored by 
C            columns.  The inverse of the diagonal matrix D is stored.
C            No fill in is allowed.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     INTEGER NL, IL(N+1), JL(NL), NU, IU(N+1), JU(NU)
C     INTEGER NROW(N), NCOL(N)
C     REAL    A(NELT), L(NL), U(NU), DINV(N)
C
C     CALL SSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, 
C    $    DINV, NU, IU, JU, U, NROW, NCOL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Real A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower 
C         triangle of the matrix is stored.
C NL     :OUT      Integer.
C         Number of non-zeros in the EL array.
C IL     :OUT      Integer IL(N+1).
C JL     :OUT      Integer JL(NL).
C L      :OUT      Real     L(NL).
C         IL, JL, L  contain the unit ower  triangular factor of  the
C         incomplete decomposition  of some  matrix stored  in   SLAP
C         Row format.     The   Diagonal  of ones  *IS*  stored.  See
C         "DESCRIPTION", below for more details about the SLAP format.
C NU     :OUT      Integer.
C         Number of non-zeros in the U array.     
C IU     :OUT      Integer IU(N+1).
C JU     :OUT      Integer JU(NU).
C U      :OUT      Real     U(NU).
C         IU, JU, U contain   the unit upper triangular factor of the
C         incomplete  decomposition    of some matrix  stored in SLAP
C         Column  format.   The Diagonal of ones   *IS*  stored.  See 
C         "Description", below  for  more  details  about  the   SLAP 
C         format.
C NROW   :WORK     Integer NROW(N).
C         NROW(I) is the number of non-zero elements in the I-th row
C         of L.
C NCOL   :WORK     Integer NCOL(N).
C         NCOL(I) is the number of non-zero elements in the I-th 
C         column of U.
C
C *Description
C       IL, JL, L should contain the unit  lower triangular factor of
C       the incomplete decomposition of the A matrix  stored in SLAP
C       Row format.  IU, JU, U should contain  the unit upper factor
C       of the  incomplete decomposition of  the A matrix  stored in
C       SLAP Column format This ILU factorization can be computed by
C       the SSILUS routine.  The diagonals (which is all one's) are
C       stored.
C
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       ==================== S L A P Row format ====================
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the real
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C       
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|  
C       |51  0 53  0 55|  
C
C *Precision:           Single Precision
C *See Also:
C       SILUR
C***REFERENCES  1. Gene Golub & Charles Van Loan, "Matrix Computations",
C                 John Hopkins University Press; 3 (1983) IBSN 
C                 0-8018-3010-9.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSILUS
cgzh bug fix      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NL, IL(NL), JL(NL)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NL, IL(NL+2), JL(NL)
cgzh bug fix      INTEGER NU, IU(NU), JU(NU), NROW(N), NCOL(N)
      INTEGER NU, IU(NU), JU(NU+2), NROW(N), NCOL(N)
      REAL    A(NELT), L(NL), DINV(N), U(NU)
C         
C         Count number of elements in each row of the lower triangle.
C***FIRST EXECUTABLE STATEMENT  SSILUS
      DO 10 I=1,N
         NROW(I) = 0
         NCOL(I) = 0
 10   CONTINUE
      DO 30 ICOL = 1, N
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               IF( IA(J).LT.ICOL ) THEN
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
                  NROW(IA(J)) = NROW(IA(J)) + 1
                  IF( ISYM.NE.0 ) NCOL(IA(J)) = NCOL(IA(J)) + 1
               ENDIF
 20         CONTINUE
         ENDIF
 30   CONTINUE
      JU(1) = 1
      IL(1) = 1
      DO 40 ICOL = 1, N
         IL(ICOL+1) = IL(ICOL) + NROW(ICOL)
         JU(ICOL+1) = JU(ICOL) + NCOL(ICOL)
         NROW(ICOL) = IL(ICOL)
         NCOL(ICOL) = JU(ICOL)
 40   CONTINUE
C         
C         Copy the matrix A into the L and U structures.
      DO 60 ICOL = 1, N
         DINV(ICOL) = A(JA(ICOL))
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               IROW = IA(J)
               IF( IROW.LT.ICOL ) THEN
C         Part of the upper triangle.
                  IU(NCOL(ICOL)) = IROW
                  U(NCOL(ICOL)) = A(J)
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
C         Part of the lower triangle (stored by row).
                  JL(NROW(IROW)) = ICOL
                  L(NROW(IROW)) = A(J)
                  NROW(IROW) = NROW(IROW) + 1
                  IF( ISYM.NE.0 ) THEN
C         Symmetric...Copy lower triangle into upper triangle as well.
                     IU(NCOL(IROW)) = ICOL
                     U(NCOL(IROW)) = A(J)
                     NCOL(IROW) = NCOL(IROW) + 1
                  ENDIF
               ENDIF
 50         CONTINUE
         ENDIF
 60   CONTINUE
C
C         Sort the rows of L and the columns of U.
      DO 110 K = 2, N
         JBGN = JU(K)
         JEND = JU(K+1)-1
         IF( JBGN.LT.JEND ) THEN
            DO 80 J = JBGN, JEND-1
               DO 70 I = J+1, JEND
                  IF( IU(J).GT.IU(I) ) THEN
                     ITEMP = IU(J)
                     IU(J) = IU(I)
                     IU(I) = ITEMP
                     TEMP = U(J)
                     U(J) = U(I)
                     U(I) = TEMP
                  ENDIF
 70            CONTINUE
 80         CONTINUE
         ENDIF
         IBGN = IL(K)
         IEND = IL(K+1)-1
         IF( IBGN.LT.IEND ) THEN
            DO 100 I = IBGN, IEND-1
               DO 90 J = I+1, IEND
                  IF( JL(I).GT.JL(J) ) THEN
                     JTEMP = JU(I)
                     JU(I) = JU(J)
                     JU(J) = JTEMP
                     TEMP = L(I)
                     L(I) = L(J)
                     L(J) = TEMP
                  ENDIF
 90            CONTINUE
 100        CONTINUE
         ENDIF
 110  CONTINUE
C
C         Perform the incomplete LDU decomposition.
      DO 300 I=2,N
C         
C           I-th row of L
         INDX1 = IL(I)
         INDX2 = IL(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 200
         DO 190 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 180
            INDXR1 = INDX1
            INDXR2 = INDX - 1
            INDXC1 = JU(JL(INDX))
            INDXC2 = JU(JL(INDX)+1) - 1
            IF(INDXC1 .GT. INDXC2) GO TO 180
 160        KR = JL(INDXR1)
 170        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 170
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 160
            ELSEIF(KR .EQ. KC) THEN
               L(INDX) = L(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 160
            ENDIF
 180        L(INDX) = L(INDX)/DINV(JL(INDX))
 190     CONTINUE
C         
C         ith column of u
 200     INDX1 = JU(I)
         INDX2 = JU(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 260
         DO 250 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 240
            INDXC1 = INDX1
            INDXC2 = INDX - 1
            INDXR1 = IL(IU(INDX))
            INDXR2 = IL(IU(INDX)+1) - 1
            IF(INDXR1 .GT. INDXR2) GO TO 240
 210        KR = JL(INDXR1)
 220        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 220
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 210
            ELSEIF(KR .EQ. KC) THEN
               U(INDX) = U(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 210
            ENDIF
 240        U(INDX) = U(INDX)/DINV(IU(INDX))
 250     CONTINUE
C         
C         ith diagonal element
 260     INDXR1 = IL(I)
         INDXR2 = IL(I+1) - 1
         IF(INDXR1 .GT. INDXR2) GO TO 300
         INDXC1 = JU(I)
         INDXC2 = JU(I+1) - 1
         IF(INDXC1 .GT. INDXC2) GO TO 300
 270     KR = JL(INDXR1)
 280     KC = IU(INDXC1)
         IF(KR .GT. KC) THEN
            INDXC1 = INDXC1 + 1
            IF(INDXC1 .LE. INDXC2) GO TO 280
         ELSEIF(KR .LT. KC) THEN
            INDXR1 = INDXR1 + 1
            IF(INDXR1 .LE. INDXR2) GO TO 270
         ELSEIF(KR .EQ. KC) THEN
            DINV(I) = DINV(I) - L(INDXR1)*DINV(KC)*U(INDXC1)
            INDXR1 = INDXR1 + 1
            INDXC1 = INDXC1 + 1
            IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 270
         ENDIF
C         
 300  CONTINUE
C         
C         replace diagonal lts by their inverses.
      DO 430 I=1,N
         DINV(I) = 1./DINV(I)
 430  CONTINUE
C         
      RETURN
C------------- LAST LINE OF SSILUS FOLLOWS ----------------------------
      END
*DECK SSMV
      SUBROUTINE SSMV( N, X, Y, NELT, IA, JA, A, ISYM )
C***BEGIN PROLOGUE  SSMV
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSMV-S),
C             Matrix Vector Multiply, Sparse
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP Column Format Sparse Matrix Vector Product.
C            Routine to calculate the sparse matrix vector product:
C            Y = A*X.
C***DESCRIPTION
C *Usage:
C     INTEGER  N, NELT, IA(NELT), JA(N+1), ISYM
C     REAL     X(N), Y(N), A(NELT)
C
C     CALL SSMV(N, X, Y, NELT, IA, JA, A, ISYM )
C         
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C X      :IN       Real X(N).
C         The vector that should be multiplied by the matrix.
C Y      :OUT      Real Y(N).
C         The product of the matrix and the vector.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(N+1).
C A      :IN       Integer A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C
C *Description
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       With  the SLAP  format  the "inner  loops" of  this  routine
C       should vectorize   on machines with   hardware  support  for
C       vector gather/scatter operations.  Your compiler may require
C       a  compiler directive  to  convince   it that there  are  no
C       implicit vector  dependencies.  Compiler directives  for the
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
C       with the standard SLAP distribution.
C
C *Precision:           Single Precision
C *Cautions:
C     This   routine   assumes  that  the matrix A is stored in SLAP 
C     Column format.  It does not check  for  this (for  speed)  and 
C     evil, ugly, ornery and nasty things  will happen if the matrix 
C     data  structure  is,  in fact, not SLAP Column.  Beware of the 
C     wrong data structure!!!
C
C *See Also:
C       SSMTV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSMV
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      REAL    A(NELT), X(N), Y(N)
C
C         Zero out the result vector.
C***FIRST EXECUTABLE STATEMENT  SSMV
      DO 10 I = 1, N
         Y(I) = 0.0
 10   CONTINUE
C
C         Multiply by A.
C
      DO 30 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 20 I = IBGN, IEND
            Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL)
 20      CONTINUE
 30   CONTINUE
C
      IF( ISYM.EQ.1 ) THEN
C
C         The matrix is non-symmetric.  Need to get the other half in...
C         This loops assumes that the diagonal is the first entry in
C         each column.
C
         DO 50 IROW = 1, N
            JBGN = JA(IROW)+1
            JEND = JA(IROW+1)-1
            IF( JBGN.GT.JEND ) GOTO 50
            DO 40 J = JBGN, JEND
               Y(IROW) = Y(IROW) + A(J)*X(IA(J))
 40         CONTINUE
 50      CONTINUE
      ENDIF
      RETURN
C------------- LAST LINE OF SSMV FOLLOWS ----------------------------
      END
*DECK SSMTV
      SUBROUTINE SSMTV( N, X, Y, NELT, IA, JA, A, ISYM )
C***BEGIN PROLOGUE  SSMTV
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSMTV-S),
C             Matrix transpose Vector Multiply, Sparse
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP Column Format Sparse Matrix (transpose) Vector Prdt.
C            Routine to calculate the sparse matrix vector product:
C            Y = A'*X, where ' denotes transpose.
C***DESCRIPTION
C *Usage:
C     INTEGER  N, NELT, IA(NELT), JA(N+1), ISYM
C     REAL     X(N), Y(N), A(NELT)
C
C     CALL SSMTV(N, X, Y, NELT, IA, JA, A, ISYM )
C         
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C X      :IN       Real X(N).
C         The vector that should be multiplied by the transpose of 
C         the matrix.
C Y      :OUT      Real Y(N).
C         The product of the transpose of the matrix and the vector.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(N+1).
C A      :IN       Integer A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C
C *Description
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       With  the SLAP  format  the "inner  loops" of  this  routine
C       should vectorize   on machines with   hardware  support  for
C       vector gather/scatter operations.  Your compiler may require
C       a  compiler directive  to  convince   it that there  are  no
C       implicit vector  dependencies.  Compiler directives  for the
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
C       with the standard SLAP distribution.
C       
C *Precision:           Single Precision
C *Cautions:
C     This   routine   assumes  that  the matrix A is stored in SLAP 
C     Column format.  It does not check  for  this (for  speed)  and 
C     evil, ugly, ornery and nasty things  will happen if the matrix 
C     data  structure  is,  in fact, not SLAP Column.  Beware of the 
C     wrong data structure!!!
C
C *See Also:
C       SSMV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSMTV
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      REAL    X(N), Y(N), A(NELT)
C
C         Zero out the result vector.
C***FIRST EXECUTABLE STATEMENT  SSMTV
      DO 10 I = 1, N
         Y(I) = 0.0
 10   CONTINUE
C
C         Multiply by A-Transpose.
C         A-Transpose is stored by rows...
      DO 30 IROW = 1, N
         IBGN = JA(IROW)
         IEND = JA(IROW+1)-1
         DO 20 I = IBGN, IEND
            Y(IROW) = Y(IROW) + A(I)*X(IA(I))
 20      CONTINUE
 30   CONTINUE
C
      IF( ISYM.EQ.1 ) THEN
C
C         The matrix is non-symmetric.  Need to get the other half in...
C         This loops assumes that the diagonal is the first entry in
C         each column.
C
         DO 50 ICOL = 1, N
            JBGN = JA(ICOL)+1
            JEND = JA(ICOL+1)-1
            IF( JBGN.GT.JEND ) GOTO 50
            DO 40 J = JBGN, JEND
               Y(IA(J)) = Y(IA(J)) + A(J)*X(ICOL)
 40         CONTINUE
 50      CONTINUE
      ENDIF
      RETURN
C------------- LAST LINE OF SSMTV FOLLOWS ----------------------------
      END
*DECK SSDI
      SUBROUTINE SSDI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C***BEGIN PROLOGUE  SSDI
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213  (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSDI-S),
C             Linear system solve, Sparse, Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Diagonal Matrix Vector Multiply.
C            Routine to calculate the product  X = DIAG*B,  
C            where DIAG is a diagonal matrix.
C***DESCRIPTION
C *Usage:
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Vector to multiply the diagonal by.
C X      :OUT      Real X(N).
C         Result of DIAG*B.
C NELT   :DUMMY    Integer.
C         Retained for compatibility with SLAP MSOLVE calling sequence.
C IA     :DUMMY    Integer IA(NELT).
C         Retained for compatibility with SLAP MSOLVE calling sequence.
C JA     :DUMMY    Integer JA(N+1).
C         Retained for compatibility with SLAP MSOLVE calling sequence.
C  A     :DUMMY    Real A(NELT).
C         Retained for compatibility with SLAP MSOLVE calling sequence.
C ISYM   :DUMMY    Integer.
C         Retained for compatibility with SLAP MSOLVE calling sequence.
C RWORK  :IN       Real RWORK(USER DEFINABLE).
C         Work array holding the diagonal of some matrix to scale
C         B by.  This array must be set by the user or by a call
C         to the slap routine SSDS or SSD2S.  The length of RWORK
C         must be > IWORK(4)+N.
C IWORK  :IN       Integer IWORK(10).
C         IWORK(4) holds the offset into RWORK for the diagonal matrix
C         to scale B by.  This is usually set up by the SLAP pre-
C         conditioner setup routines SSDS or SSD2S.
C
C *Description:
C         This routine is supplied with the SLAP package to perform
C         the  MSOLVE  operation for iterative drivers that require
C         diagonal  Scaling  (e.g., SSDCG, SSDBCG).   It  conforms
C         to the SLAP MSOLVE CALLING CONVENTION  and hence does not
C         require an interface routine as do some of the other pre-
C         conditioners supplied with SLAP.
C
C *Precision:           Single Precision
C *See Also:
C       SSDS, SSD2S
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSDI
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(10)
      REAL    B(N), X(N), A(NELT), RWORK(1)
C
C         Determine where the inverse of the diagonal 
C         is in the work array and then scale by it.
C***FIRST EXECUTABLE STATEMENT  SSDI
      LOCD = IWORK(4) - 1
      DO 10 I = 1, N
         X(I) = RWORK(LOCD+I)*B(I)
 10   CONTINUE
      RETURN
C------------- LAST LINE OF SSDI FOLLOWS ----------------------------
      END
*DECK SSLI
      SUBROUTINE SSLI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK )
C***BEGIN PROLOGUE  SSLI
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSLI-S),
C             Linear system solve, Sparse, Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP MSOLVE for Lower Triangle Matrix.
C            This routine acts as an interface between the SLAP generic 
C            MSOLVE calling convention and the routine that actually 
C                      -1 
C            computes L  B = X.
C
C *Description
C       See the Description of SLLI2 for the gory details.
C***ROUTINES CALLED  SLLI2
C***END PROLOGUE  SSLI
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(10)
      REAL    B(N), X(N), A(NELT), RWORK(1)
C***FIRST EXECUTABLE STATEMENT  SSLI
C         
      NEL = IWORK(1)
      LOCIEL = IWORK(2)
      LOCJEL = IWORK(3)
      LOCEL = IWORK(4)
      CALL SSLI2(N, B, X, NEL, IWORK(LOCIEL), IWORK(LOCJEL),
     $     RWORK(LOCEL))
C         
      RETURN
C------------- LAST LINE OF SSLI FOLLOWS ----------------------------
      END
*DECK SSLI2
      SUBROUTINE SSLI2(N, B, X, NEL, IEL, JEL, EL)
C***BEGIN PROLOGUE  SSLI2
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSLI2-S),
C             Linear system solve, Sparse, Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP for Lower Triangle Matrix Backsolve.
C            Routine to solve a system of the form  Lx = b , where
C            L is a lower triangular matrix.
C***DESCRIPTION
C *Usage:
C     INTEGER N,  NEL, IEL(N+1), JEL(NEL)
C     REAL    B(N), X(N), EL(NEL)
C
C     CALL SSLI2( N, B, X, NEL, IEL, JEL, EL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right hand side vector.
C X      :OUT      Real X(N).
C         Solution to Lx = b.
C NEL    :IN       Integer.
C         Number of non-zeros in the EL array.
C IEL    :IN       Integer IEL(N+1).
C JEL    :IN       Integer JEL(NEL).
C EL     :IN       Real     EL(NEL).
C         IEL, JEL, EL contain the unit lower triangular factor   of
C         the incomplete decomposition   of the A  matrix  stored in 
C         SLAP Row format.  The diagonal of  ones *IS* stored.  This 
C         structure can be set up by the  SS2LT  routine.  See "LONG 
C         DESCRIPTION", below for more details about  the  SLAP  Row 
C         format.
C
C *Description:
C       This routine is supplied with the SLAP package  as a routine
C       to  perform the  MSOLVE operation in  the SIR for the driver
C       routine SSGS.  It must be called via the SLAP MSOLVE calling
C       sequence convention interface routine SSLI.
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
C               **** SLAP MSOLVE CALLING CONVENTION ****
C
C       ==================== S L A P Row format ====================
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the real
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C       
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|  
C       |51  0 53  0 55|  
C
C       With  the SLAP  Row format  the "inner loop" of this routine
C       should vectorize   on machines with   hardware  support  for
C       vector gather/scatter operations.  Your compiler may require
C       a  compiler directive  to  convince   it that there  are  no
C       implicit vector  dependencies.  Compiler directives  for the
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
C       with the standard SLAP distribution.
C
C *Precision: Single Precision
C *See Also:
C         SSLI
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSLI2
      INTEGER N, NEL, IEL(NEL), JEL(NEL)
      REAL    B(N), X(N), EL(NEL)
C
C         Initialize the solution by copying the right hands side
C         into it.
C***FIRST EXECUTABLE STATEMENT  SSLI2
      DO 10 I=1,N
         X(I) = B(I)
 10   CONTINUE
C         
      DO 30 ICOL = 1, N
         X(ICOL) = X(ICOL)/EL(JEL(ICOL))
         JBGN = JEL(ICOL) + 1
         JEND = JEL(ICOL+1) - 1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               X(IEL(J)) = X(IEL(J)) - EL(J)*X(ICOL)
 20         CONTINUE
         ENDIF
 30   CONTINUE
C         
      RETURN
C------------- LAST LINE OF SSLI2 FOLLOWS ----------------------------
      END
*DECK SSLLTI
      SUBROUTINE SSLLTI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C***BEGIN PROLOGUE  SSLLTI
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSLLTI-S),
C             Linear system solve, Sparse, Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP MSOLVE for LDL' (IC) Factorization.
C            This routine acts as an interface between the SLAP generic
C            MSOLVE calling convention and the routine that actually 
C                           -1 
C            computes (LDL')  B = X.
C***DESCRIPTION
C       See the DESCRIPTION of SLLTI2 for the gory details.
C***ROUTINES CALLED SLLTI2
C
C***END PROLOGUE  SSLLTI
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(1)
      REAL    B(1), X(1), A(NELT), RWORK(1)
C         
C***FIRST EXECUTABLE STATEMENT  SSLLTI
      NEL = IWORK(1)
      LOCIEL = IWORK(3)
      LOCJEL = IWORK(2)
      LOCEL  = IWORK(4)
      LOCDIN = IWORK(5)
      CALL SLLTI2(N, B, X, NEL, IWORK(LOCIEL), IWORK(LOCJEL),
     $     RWORK(LOCEL), RWORK(LOCDIN))
C         
      RETURN
C------------- LAST LINE OF SSLLTI FOLLOWS ----------------------------
      END
*DECK SLLTI2
      SUBROUTINE SLLTI2(N, B, X, NEL, IEL, JEL, EL, DINV)
C***BEGIN PROLOGUE  SLLTI2
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SLLTI2-S),
C             Symmetric Linear system solve, Sparse, 
C             Iterative Precondition, Incomplete Factorization
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP back solve routine for LDL' Factorization.
C            Routine to solve a system of the  form  L*D*L' X  =  B,
C            where L is a unit lower triangular  matrix  and  D is a 
C            diagonal matrix and ' means transpose.
C***DESCRIPTION
C *Usage:
C     INTEGER N,  NEL, IEL(N+1), JEL(NEL)
C     REAL    B(N), X(N), EL(NEL), DINV(N)
C
C     CALL SLLTI2( N, B, X, NEL, IEL, JEL, EL, DINV )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right hand side vector.
C X      :OUT      Real X(N).
C         Solution to L*D*L' x = b.
C NEL    :IN       Integer.
C         Number of non-zeros in the EL array.
C IEL    :IN       Integer IEL(N+1).
C JEL    :IN       Integer JEL(NEL).
C EL     :IN       Real     EL(NEL).
C         IEL, JEL, EL contain the unit lower triangular factor   of
C         the incomplete decomposition   of the A  matrix  stored in 
C         SLAP Row format.   The diagonal of ones *IS* stored.  This 
C         structure can be set  up  by  the SS2LT routine. See 
C         "Description", below for more details about the   SLAP Row
C         format.
C DINV   :IN       Real DINV(N).
C         Inverse of the diagonal matrix D.
C
C *Description:
C       This routine is supplied with  the SLAP package as a routine
C       to perform the MSOLVE operation in the SCG iteration routine
C       for  the driver  routine SSICCG.   It must be called via the
C       SLAP  MSOLVE calling sequence  convention  interface routine
C       SSLLI.
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
C               **** SLAP MSOLVE CALLING CONVENTION ****
C
C       IEL, JEL, EL should contain the unit lower triangular factor
C       of  the incomplete decomposition of  the A matrix  stored in
C       SLAP Row format.   This IC factorization  can be computed by
C       the  SSICS routine.  The  diagonal  (which is all one's) is
C       stored.  
C
C       ==================== S L A P Row format ====================
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the real
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C       
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|  
C       |51  0 53  0 55|  
C
C       With  the SLAP  Row format  the "inner loop" of this routine
C       should vectorize   on machines with   hardware  support  for
C       vector gather/scatter operations.  Your compiler may require
C       a  compiler directive  to  convince   it that there  are  no
C       implicit vector  dependencies.  Compiler directives  for the
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
C       with the standard SLAP distribution.
C
C *Precision:           Single Precision
C *See Also:
C       SSICCG, SSICS
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SLLTI2
      INTEGER N, NEL, IEL(NEL), JEL(1)
      REAL    B(N), X(N), EL(NEL), DINV(N)
C         
C         solve  l*y = b,  storing result in x.
C***FIRST EXECUTABLE STATEMENT  SLLTI2
      DO 10 I=1,N
         X(I) = B(I)
 10   CONTINUE
      DO 30 IROW = 1, N
         IBGN = IEL(IROW) + 1
         IEND = IEL(IROW+1) - 1
         IF( IBGN.LE.IEND ) THEN
            DO 20 I = IBGN, IEND
               X(IROW) = X(IROW) - EL(I)*X(JEL(I))
 20         CONTINUE
         ENDIF
 30   CONTINUE
C         
C         Solve  D*Z = Y,  storing result in X.
C         
      DO 40 I=1,N
         X(I) = X(I)*DINV(I)
 40   CONTINUE
C         
C         Solve  L-trans*X = Z.
C
      DO 60 IROW = N, 2, -1
         IBGN = IEL(IROW) + 1
         IEND = IEL(IROW+1) - 1
         IF( IBGN.LE.IEND ) THEN
            DO 50 I = IBGN, IEND
               X(JEL(I)) = X(JEL(I)) - EL(I)*X(IROW)
 50         CONTINUE
         ENDIF
 60   CONTINUE
C         
      RETURN
C------------- LAST LINE OF SLTI2 FOLLOWS ----------------------------
      END
*DECK SSLUI
      SUBROUTINE SSLUI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C***BEGIN PROLOGUE  SSLUI
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSLUI-S),
C             Non-Symmetric Linear system solve, Sparse, 
C             Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP MSOLVE for LDU Factorization.
C            This routine  acts as an  interface between  the   SLAP
C            generic MSLOVE calling convention and the routine  that 
C            actually computes:     -1
C                              (LDU)  B = X.
C***DESCRIPTION
C       See the "DESCRIPTION" of SSLUI2 for the gory details.
C***ROUTINES CALLED  SSLUI2
C***END PROLOGUE  SSLUI
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(*)
      REAL    B(N), X(N), A(NELT), RWORK(*)
C
C         Pull out the locations of the arrays holding the ILU
C         factorization.
C***FIRST EXECUTABLE STATEMENT  SSLUI
      LOCIL = IWORK(1)
      LOCJL = IWORK(2)
      LOCIU = IWORK(3)
      LOCJU = IWORK(4)
      LOCL = IWORK(5)
      LOCDIN = IWORK(6)
      LOCU = IWORK(7)
C
C         Solve the system LUx = b
      CALL SSLUI2(N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL),
     $     RWORK(LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU) )
C         
      RETURN
C------------- LAST LINE OF SSLUI FOLLOWS ----------------------------
      END
*DECK SSLUI2
      SUBROUTINE SSLUI2(N, B, X, IL, JL, L, DINV, IU, JU, U )
C***BEGIN PROLOGUE  SSLUI2
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSLUI2-S),
C             Non-Symmetric Linear system solve, Sparse, 
C             Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP Back solve for LDU Factorization.
C            Routine  to  solve a system of the form  L*D*U X  =  B,
C            where L is a unit  lower  triangular  matrix,  D  is  a 
C            diagonal matrix, and U is a unit upper triangular matrix.
C***DESCRIPTION
C *Usage:
C     INTEGER N, IL(N+1), JL(NL), IU(NU), JU(N+1)
C     REAL    B(N), X(N), L(NL), DINV(N), U(NU)
C
C     CALL SSLUI2( N, B, X, IL, JL, L, DINV, IU, JU, U )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right hand side.
C X      :OUT      Real X(N).
C         Solution of L*D*U x = b.
C NEL    :IN       Integer.
C         Number of non-zeros in the EL array.
C IL     :IN       Integer IL(N+1).
C JL     :IN       Integer JL(NL).
C  L     :IN       Real     L(NL).
C         IL, JL, L contain the unit  lower triangular factor of the
C         incomplete decomposition of some matrix stored in SLAP Row
C         format.  The diagonal of ones *IS* stored.  This structure
C         can   be   set up  by   the  SSILUS routine.   See 
C         "DESCRIPTION", below  for more   details about   the  SLAP
C         format.
C DINV   :IN       Real DINV(N).
C         Inverse of the diagonal matrix D.
C NU     :IN       Integer.
C         Number of non-zeros in the U array.     
C IU     :IN       Integer IU(N+1).
C JU     :IN       Integer JU(NU).
C U      :IN       Real     U(NU).
C         IU, JU, U contain the unit upper triangular factor  of the
C         incomplete decomposition  of  some  matrix stored in  SLAP
C         Column format.   The diagonal of ones  *IS* stored.   This
C         structure can be set up  by the SSILUS routine.  See
C         "DESCRIPTION", below   for  more   details about  the SLAP
C         format.
C
C *Description:
C       This routine is supplied with  the SLAP package as a routine
C       to  perform  the  MSOLVE operation  in   the  SIR and   SBCG
C       iteration routines for  the  drivers SSILUR and SSLUBC.   It
C       must  be called  via   the  SLAP  MSOLVE  calling   sequence
C       convention interface routine SSLUI.
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
C               **** SLAP MSOLVE CALLING CONVENTION ****
C
C       IL, JL, L should contain the unit lower triangular factor of
C       the incomplete decomposition of the A matrix  stored in SLAP
C       Row format.  IU, JU, U should contain  the unit upper factor
C       of the  incomplete decomposition of  the A matrix  stored in
C       SLAP Column format This ILU factorization can be computed by
C       the SSILUS routine.  The diagonals (which is all one's) are
C       stored.
C
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       ==================== S L A P Row format ====================
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the real
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C       
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|  
C       |51  0 53  0 55|  
C
C       With  the SLAP  format  the "inner  loops" of  this  routine
C       should vectorize   on machines with   hardware  support  for
C       vector gather/scatter operations.  Your compiler may require
C       a  compiler directive  to  convince   it that there  are  no
C       implicit vector  dependencies.  Compiler directives  for the
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
C       with the standard SLAP distribution.
C
C *Precision:           Single Precision
C *See Also:
C       SSILUS
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSLUI2
      INTEGER N, IL(1), JL(1), IU(1), JU(1)
      REAL    B(N), X(N), L(1), DINV(N), U(1)
C         
C         Solve  L*Y = B,  storing result in X, L stored by rows.
C***FIRST EXECUTABLE STATEMENT  SSLUI2
      DO 10 I = 1, N
         X(I) = B(I)
 10   CONTINUE
      DO 30 IROW = 2, N
         JBGN = IL(IROW)
         JEND = IL(IROW+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               X(IROW) = X(IROW) - L(J)*X(JL(J))
 20         CONTINUE
         ENDIF
 30   CONTINUE
C         
C         Solve  D*Z = Y,  storing result in X.
      DO 40 I=1,N
         X(I) = X(I)*DINV(I)
 40   CONTINUE
C         
C         Solve  U*X = Z, U stored by columns.
      DO 60 ICOL = N, 2, -1
         JBGN = JU(ICOL)
         JEND = JU(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               X(IU(J)) = X(IU(J)) - U(J)*X(ICOL)
 50         CONTINUE
         ENDIF
 60   CONTINUE
C         
      RETURN
C------------- LAST LINE OF SSLUI2 FOLLOWS ----------------------------
      END
*DECK SSLUTI
      SUBROUTINE SSLUTI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C***BEGIN PROLOGUE  SSLUTI
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSLUTI-S),
C             Linear system solve, Sparse, Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP MTSOLV for LDU Factorization.
C            This routine acts as  an  interface  between  the  SLAP
C            generic MTSOLV calling convention and  the routine that  
C            actually computes:       -T
C                                (LDU)  B = X.
C***DESCRIPTION
C       See the "DESCRIPTION" of SSLUI4 for the gory details.
C***ROUTINES CALLED  SSLUI4
C***END PROLOGUE  SSLUTI
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(10)
      REAL    B(N), X(N), A(N), RWORK(1)
C         
C         Pull out the pointers to the L, D and U matricies and call
C         the workhorse routine.
C***FIRST EXECUTABLE STATEMENT  SSLUTI
      LOCIL = IWORK(1)
      LOCJL = IWORK(2)
      LOCIU = IWORK(3)
      LOCJU = IWORK(4)
      LOCL = IWORK(5)
      LOCDIN = IWORK(6)
      LOCU = IWORK(7)
C
      CALL SSLUI4(N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL),
     $     RWORK(LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU))
C         
      RETURN
C------------- LAST LINE OF SSLUTI FOLLOWS ----------------------------
      END
*DECK SSLUI4
      SUBROUTINE SSLUI4(N, B, X, IL, JL, L, DINV, IU, JU, U )
C***BEGIN PROLOGUE  SSLUI4
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSLUI4-S),
C             Non-Symmetric Linear system solve, Sparse, 
C             Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP back solve for LDU Factorization.
C            Routine to solve a system of the form  (L*D*U)' X =  B,
C            where L is a unit  lower  triangular  matrix,  D  is  a 
C            diagonal matrix, and  U  is  a  unit  upper  triangular 
C            matrix and ' denotes transpose.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NL, IL(N+1), JL(NL), NU, IU(N+1), JU(NU)
C     REAL    B(N), X(N), L(NEL), DINV(N), U(NU)
C
C     CALL SSLUI4( N, B, X, IL, JL, L, DINV, IU, JU, U )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right hand side.
C X      :OUT      Real X(N).
C         Solution of (L*D*U)trans x = b.
C IL     :IN       Integer IL(N+1).
C JL     :IN       Integer JL(NL).
C  L     :IN       Real     L(NL).
C         IL, JL, L contain the unit lower triangular  factor of the
C         incomplete decomposition of some matrix stored in SLAP Row
C         format.  The diagonal of ones *IS* stored.  This structure
C         can    be set  up  by   the  SSILUS routine.   See  
C         "DESCRIPTION",  below for  more  details about  the   SLAP
C         format.
C DINV   :IN       Real DINV(N).
C         Inverse of the diagonal matrix D.
C IU     :IN       Integer IU(N+1).
C JU     :IN       Integer JU(NU).
C U      :IN       Real     U(NU).
C         IU, JU, U contain the  unit upper triangular factor of the
C         incomplete  decomposition of some  matrix stored  in  SLAP
C         Column  format.   The diagonal of  ones *IS* stored.  This
C         structure can be set up by the  SSILUS routine.  See
C         "DESCRIPTION",  below for  more  details  about  the  SLAP
C         format.
C
C *Description:
C       This routine is supplied with the SLAP package as  a routine
C       to  perform  the  MTSOLV  operation  in  the SBCG  iteration
C       routine for the  driver SSLUBC.   It must  be called via the
C       SLAP  MTSOLV calling  sequence convention interface  routine
C       SSLUTI.
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
C               **** SLAP MSOLVE CALLING CONVENTION ****
C
C       IL, JL, L should contain the unit lower triangular factor of
C       the incomplete decomposition of the A matrix  stored in SLAP
C       Row format.  IU, JU, U should contain  the unit upper factor
C       of the  incomplete decomposition of  the A matrix  stored in
C       SLAP Column format This ILU factorization can be computed by
C       the SSILUS routine.  The diagonals (which is all one's) are
C       stored.
C
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       ==================== S L A P Row format ====================
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the real
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C       
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|  
C       |51  0 53  0 55|  
C
C       With  the SLAP  format  the "inner  loops" of  this  routine
C       should vectorize   on machines with   hardware  support  for
C       vector gather/scatter operations.  Your compiler may require
C       a  compiler directive  to  convince   it that there  are  no
C       implicit vector  dependencies.  Compiler directives  for the
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
C       with the standard SLAP distribution.
C
C *Precision:           Single Precision
C *See Also:
C       SSILUS
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSLUI4
      INTEGER N, IL(*), JL(*), IU(*), JU(*)
      REAL    B(N), X(N), L(*), DINV(N), U(*)
C         
C***FIRST EXECUTABLE STATEMENT  SSLUI4
      DO 10 I=1,N
         X(I) = B(I)
 10   CONTINUE
C         
C         Solve  U'*Y = X,  storing result in X, U stored by columns.
      DO 80 IROW = 2, N
         JBGN = JU(IROW)
         JEND = JU(IROW+1) - 1
         IF( JBGN.LE.JEND ) THEN
            DO 70 J = JBGN, JEND
               X(IROW) = X(IROW) - U(J)*X(IU(J))
 70         CONTINUE
         ENDIF
 80   CONTINUE
C         
C         Solve  D*Z = Y,  storing result in X.
      DO 90 I = 1, N
         X(I) = X(I)*DINV(I)
 90   CONTINUE
C         
C         Solve  L'*X = Z, L stored by rows.
      DO 110 ICOL = N, 2, -1
         JBGN = IL(ICOL)
         JEND = IL(ICOL+1) - 1
         IF( JBGN.LE.JEND ) THEN
            DO 100 J = JBGN, JEND
               X(JL(J)) = X(JL(J)) - L(J)*X(ICOL)
 100        CONTINUE
         ENDIF
 110  CONTINUE
      RETURN
C------------- LAST LINE OF SSLUI4 FOLLOWS ----------------------------
      END
*DECK SSMMTI
      SUBROUTINE SSMMTI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK )
C***BEGIN PROLOGUE  SSMMTI
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSMMTI-S),
C             Linear system solve, Sparse, Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP MSOLVE for LDU Factorization of Normal Equations.
C            This routine acts as  an  interface  between  the  SLAP
C            generic MMTSLV calling convention and the routine  that 
C            actually computes:            -1
C                            [(LDU)*(LDU)']  B = X.
C***DESCRIPTION
C       See the "DESCRIPTION" of SSMMI2 for the gory details.
C***ROUTINES CALLED  SSMMI2
C***END PROLOGUE  SSMMTI
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(10)
      REAL    B(N), X(N), A(NELT), RWORK(1)
C
C         Pull out the locations of the arrays holding the ILU
C         factorization.
C***FIRST EXECUTABLE STATEMENT  SSMMTI
      LOCIL = IWORK(1)
      LOCJL = IWORK(2)
      LOCIU = IWORK(3)
      LOCJU = IWORK(4)
      LOCL = IWORK(5)
      LOCDIN = IWORK(6)
      LOCU = IWORK(7)
C         
      CALL SSMMI2(N, B, X, IWORK(LOCIL), IWORK(LOCJL),
     $     RWORK(LOCL), RWORK(LOCDIN), IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCU))
C         
      RETURN
C------------- LAST LINE OF SSMMTI FOLLOWS ----------------------------
      END
*DECK SSMMI2
      SUBROUTINE SSMMI2( N, B, X, IL, JL, L, DINV, IU, JU, U )
C***BEGIN PROLOGUE  SSMMI2
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSMMI2-S),
C             Linear system, Sparse, Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP Back solve for LDU Factorization of Normal Equations.
C            To solve a system of the form (L*D*U)*(L*D*U)' X  =  B,
C            where  L  is a unit lower triangular matrix,  D   is  a 
C            diagonal matrix, and  U  is  a  unit  upper  triangular 
C            matrix and ' denotes transpose.
C***DESCRIPTION
C *Usage:
C     INTEGER N, IL(N+1), JL(NL), IU(N+1), JU(NU)
C     REAL    B(N), X(N), L(NL), DINV(N), U(NU)
C
C     CALL SSMMI2( N, B, X, IL, JE, L, DINV, IU, JU, U )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right hand side.
C X      :OUT      Real X(N).
C         Solution of (L*D*U)(L*D*U)trans x = b.
C IL     :IN       Integer IL(N+1).
C JL     :IN       Integer JL(NL).
C  L     :IN       Real     L(NL).
C         IL, JL, L contain the unit lower  triangular factor of the
C         incomplete decomposition of some matrix stored in SLAP Row
C         format.  The diagonal of ones *IS* stored.  This structure
C         can  be  set up by   the  SSILUS   routine.    See 
C         "DESCRIPTION", below for  more   details   about  the SLAP
C         format.
C DINV   :IN       Real DINV(N).
C         Inverse of the diagonal matrix D.
C IU     :IN       Integer IU(N+1).
C JU     :IN       Integer JU(NU).
C U      :IN       Real     U(NU).
C         IU, JU, U contain the unit upper  triangular factor of the
C         incomplete decomposition  of   some matrix stored in  SLAP
C         Column  format.  The diagonal  of  ones *IS* stored.  This
C         structure can be set up  by the SSILUS routine.  See
C         "DESCRIPTION",  below  for  more  details  about  the SLAP
C         format.
C
C *Description:
C       This routine is supplied with the SLAP package as  a routine
C       to  perform  the  MSOLVE  operation  in  the SBCGN iteration
C       routine for the  driver SSLUCN.   It must  be called via the
C       SLAP  MSOLVE calling  sequence convention interface  routine
C       SSMMTI.
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
C               **** SLAP MSOLVE CALLING CONVENTION ****
C
C       IL, JL, L should contain the unit lower triangular factor of
C       the incomplete decomposition of the A matrix  stored in SLAP
C       Row format.  IU, JU, U should contain  the unit upper factor
C       of the  incomplete decomposition of  the A matrix  stored in
C       SLAP Column format This ILU factorization can be computed by
C       the SSILUS routine.  The diagonals (which is all one's) are
C       stored.
C
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       ==================== S L A P Row format ====================
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the real
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C       
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|  
C       |51  0 53  0 55|  
C
C       With  the SLAP  format  the "inner  loops" of  this  routine
C       should vectorize   on machines with   hardware  support  for
C       vector gather/scatter operations.  Your compiler may require
C       a  compiler directive  to  convince   it that there  are  no
C       implicit vector  dependencies.  Compiler directives  for the
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
C       with the standard SLAP distribution.
C
C *Precision:           Single Precision
C *See Also:
C       SSILUS
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SSMMI2
      INTEGER N, IL(1), JL(1), IU(1), JU(1)
      REAL    B(N), X(N), L(1), DINV(N), U(N)
C         
C         Solve  L*Y = B,  storing result in X, L stored by rows.
C***FIRST EXECUTABLE STATEMENT  SSMMI2
      DO 10 I = 1, N
         X(I) = B(I)
 10   CONTINUE
      DO 30 IROW = 2, N
         JBGN = IL(IROW)
         JEND = IL(IROW+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               X(IROW) = X(IROW) - L(J)*X(JL(J))
 20         CONTINUE
         ENDIF
 30   CONTINUE
C         
C         Solve  D*Z = Y,  storing result in X.
      DO 40 I=1,N
         X(I) = X(I)*DINV(I)
 40   CONTINUE
C         
C         Solve  U*X = Z, U stored by columns.
      DO 60 ICOL = N, 2, -1
         JBGN = JU(ICOL)
         JEND = JU(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               X(IU(J)) = X(IU(J)) - U(J)*X(ICOL)
 50         CONTINUE
         ENDIF
 60   CONTINUE
C         
C         Solve  U'*Y = X,  storing result in X, U stored by columns.
      DO 80 IROW = 2, N
         JBGN = JU(IROW)
         JEND = JU(IROW+1) - 1
         IF( JBGN.LE.JEND ) THEN
            DO 70 J = JBGN, JEND
               X(IROW) = X(IROW) - U(J)*X(IU(J))
 70         CONTINUE
         ENDIF
 80   CONTINUE
C         
C         Solve  D*Z = Y,  storing result in X.
      DO 90 I = 1, N
         X(I) = X(I)*DINV(I)
 90   CONTINUE
C         
C         Solve  L'*X = Z, L stored by rows.
      DO 110 ICOL = N, 2, -1
         JBGN = IL(ICOL)
         JEND = IL(ICOL+1) - 1
         IF( JBGN.LE.JEND ) THEN
            DO 100 J = JBGN, JEND
               X(JL(J)) = X(JL(J)) - L(J)*X(ICOL)
 100        CONTINUE
         ENDIF
 110  CONTINUE
C         
      RETURN
C------------- LAST LINE OF SSMMI2 FOLLOWS ----------------------------
      END
*DECK SNRM2
      REAL FUNCTION SNRM2 (N, SX, INCX)
C***BEGIN PROLOGUE  SNRM2
C***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A3B
C***TYPE      SINGLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
C***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
C             LINEAR ALGEBRA, UNITARY, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C
C     --Output--
C    SNRM2  single precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in SX with storage
C     increment INCX .
C     If N .LE. 0, return with result = 0.
C     If N .GE. 1, then INCX must be .GE. 1
C
C     Four Phase Method using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  SQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief Outline of Algorithm.
C
C     Phase 1 scans zero components.
C     Move to phase 2 when a component is nonzero and .LE. CUTLO
C     Move to phase 3 when a component is .GT. CUTLO
C     Move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI.
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows:
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SNRM2
      INTEGER NEXT
      REAL SX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO, ONE
      SAVE CUTLO, CUTHI, ZERO, ONE
      DATA ZERO, ONE /0.0E0, 1.0E0/
C
      DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
C***FIRST EXECUTABLE STATEMENT  SNRM2
      IF (N .GT. 0) GO TO 10
         SNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C
C                                                 BEGIN MAIN LOOP
C
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF (ABS(SX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF (SX(I) .EQ. ZERO) GO TO 200
      IF (ABS(SX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
C
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / SX(I)) / SX(I)
  105 XMAX = ABS(SX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF (ABS(SX(I)) .GT. CUTLO) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF (ABS(SX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / SX(I))**2
         XMAX = ABS(SX(I))
         GO TO 200
C
  115 SUM = SUM + (SX(I)/XMAX)**2
      GO TO 200
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI / N
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J = I,NN,INCX
      IF (ABS(SX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + SX(J)**2
      SNRM2 = SQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      SNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END
*DECK SSCAL
      SUBROUTINE SSCAL (N, SA, SX, INCX)
C***BEGIN PROLOGUE  SSCAL
C***PURPOSE  Multiply a vector by a constant.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A6
C***TYPE      SINGLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SA  single precision scale factor
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C
C     --Output--
C       SX  single precision result (unchanged if N .LE. 0)
C
C     Replace single precision SX by single precision SA*SX.
C     For I = 0 to N-1, replace SX(IX+I*INCX) with  SA * SX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SSCAL
      REAL SA, SX(*)
      INTEGER I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  SSCAL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        SX(IX) = SA*SX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I+1) = SA*SX(I+1)
        SX(I+2) = SA*SX(I+2)
        SX(I+3) = SA*SX(I+3)
        SX(I+4) = SA*SX(I+4)
   50 CONTINUE
      RETURN
      END
C