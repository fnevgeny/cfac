
*
************************************************************************
*
      SUBROUTINE ACOFZ1(Z,AM,NU,NL,AC,ANC,NDL,IOPT)
*
*  Subroutine to calculate electric dipole radial integrals and
*  transition probabilities for hydrogenic ions.  Results are obtained
*  for all transitions NU,LU - NL,LL for a given NU, NL (NU.ne.NL).
*  Direct recurrence on the radial matrix elements is used.
*  Note that an initialisation call is required.  See IOPT=2 below.
*
*  Reference: P.J. Storey and D.G. Hummer, 1991. CPC.
*
*** Non-standard FORTRAN statement follows
      IMPLICIT NONE
*
*         Subroutine arguments are:
*             Z     =  nuclear charge
*             AM    =  nuclear mass in atomic mass units
*             NU    =  upper state principal quantum number
*             NL    =  lower state principal quantum number
*             AC    =  array of transition probabilities (or radial
*                      matrix elements) for all allowable angular
*                      momentum quantum numbers
*             ANC   =  total transition probability from level NU to NL
*             NDL   =  first dimension of AC array in calling program
*             IOPT  =  option switch:
*                   =  1,  Electric dipole radial matrix elements
*                          in atomic units are stored in AC
*                   =  2,  Initialisation of exponentials and factorials
*                          Input maximum principal quantum number in NU
*
*         Storage of transition probabilities (or radial integrals)
*         is such that
*           AC(L,1)  is the probability for transition NU,L-1 to NL,L
*     and   AC(L,2)  is the probability for transition NU,L   to NL,L-1
*
*  Import
      INTEGER NU,NL,IOPT,NDL
      DOUBLE PRECISION Z,AM
*  Export
      DOUBLE PRECISION AC,ANC
*  Local
      INTEGER I,L,MAX
      DOUBLE PRECISION FAL,ALO,ANU,ANL,
     :                 AI,RP,RM,FP,FM,AL,CU,CL,BU,BL,TWOL,REM,
     :                 ZERO,ONE,TWO,FOUR,TEN,FAL1, ALO1
*
      INTEGER NMAX
      PARAMETER (NMAX=1024)
      DIMENSION FAL(NMAX),ALO(NMAX),AC(NDL,2)
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, FOUR=4.D0, TEN=10.D0)
*
*              Save arrays of factorials and exponentials
*
      SAVE FAL,ALO
      
      IF (IOPT.NE.1 .AND. IOPT.NE.2) THEN
          EXIT
      ENDIF
*
*              Initialise arrays of factorials and exponentials
*
      IF(IOPT.EQ.2) THEN
           MAX=2*NU
           IF(MAX.GT.NMAX) THEN
             WRITE(6,101)
             IOPT=-1
             RETURN
           ENDIF
           FAL(1) = ZERO
           DO I =2,MAX
              AI = REAL(I-1)
              ALO(I-1) = LOG10(AI)
              FAL(I) = ALO(I-1) + FAL(I-1)
           ENDDO
           RETURN
*
      ELSE
*
*              Reduced mass
*
           REM=ONE/(ONE+ONE/(1822.889D0*AM))
*
           DO I=1,NL+1
              AC(I,1)=ZERO
              AC(I,2)=ZERO
           ENDDO
           ANC=ZERO
*
*              Calculation of radial matrix elements, R
*              ****************************************
*
*              Using the notation R(ll,lu) = R(nl,ll,nu,lu)
*
*              Start with direct calculation of R(nl-1,nl)
*              In addition, use                 R(nl,nl-1) = 0
*              Also initialise the factors cu, cl
*
           ANU=REAL(NU)
           ANL=REAL(NL)
*
           IF (NU .GT. NL) THEN
              FAL1 = FAL(NU-NL)
              ALO1 = ALO(NU-NL)
           ELSE
              FAL1 = 0.0
              ALO1 = 0.0
           ENDIF
           FM = FAL(NU+NL+1) - FAL1 - FAL(NL+NL)
           FM = FM/TWO
     :         + (ANL+TWO)*(ALO(NU)+ALO(NL)+0.602059991327962D0)
     :         + (ANU-ANL-TWO)*ALO1 - (ANU+ANL+TWO)*ALO(NU+NL)
           FM = TEN**FM/Z/REM/FOUR
           FP = ZERO
           CU=SQRT((ANU-ANL)*(ANU+ANL))/(ANU*ANL)
           CL=ONE
*
*               Recur on R :
*               Use            R(l,l+1) and R(l+1,l)      (rm, rp)
*               to construct   R(l-1,l) and R(l,l-1)      (fm, fp)
*
           DO L=NL,1,-1
               AL=REAL(L)
               TWOL=TWO*AL
               
               IF(L.NE.NL) THEN
                   CU=SQRT((ANU-AL)*(ANU+AL))/(ANU*AL)
                   CL=SQRT((ANL-AL)*(ANL+AL))/(ANL*AL)
                   FP = ( (TWOL+ONE)*BL*RP + BU*RM) / (TWOL*CU)
                   FM = ( BL*RP + (TWOL+ONE)*BU*RM) / (TWOL*CL)
               ENDIF
*
               AC(L,1)=FP
               AC(L,2)=FM
*
               RP=FP
               RM=FM
               BU=CU
               BL=CL
           ENDDO
*
      ENDIF
*
  101 FORMAT(' INSUFFICIENT WORK SPACE IN ACOFZ1'/
     :       ' INCREASE DIMENSIONS OF FAL AND ALO TO AT LEAST',
     :       ' 2*(MAXIMUM PRINCIPAL QUANTUM NUMBER)'/)
*
      END
*
************************************************************************
*
      SUBROUTINE PIXZ1(Z,AM,NE,NL,PHE,PC,PCP,PNC,NDE,NDL,IOPT,IPCP)
*
*  Subroutine to calculate electric dipole radial integrals and  photo-
*  ionization cross-sections for hydrogenic ions.  Results are obtained
*  for all transitions NL,LL - EU,LU for a given NL and free-electron
*  energy EU.  Direct recurrence on the radial matrix elements is used.
*  Note that an initialisation call is required.   See IOPT=2 below.
*
*  Reference: P.J. Storey and D.G. Hummer, 1991. CPC.
*

********************************************************************
*  Added the option IPCP to indicate whether the partial cross sections
*  are needed. M. F. Gu, 10/25/2001.
*********************************************************************

*** Non-standard FORTRAN statement follows
      IMPLICIT NONE
*
*        Subroutine arguments are:
*            Z     =  nuclear charge
*            AM    =  nuclear mass in atomic mass units
*            NE    =  number of free electron energies
*                     (or photon energies) at which the cross-
*                     -section is to be calculated
*            NL    =  bound state principal quantum number
*            PHE   =  array of free electron energies (in a.u.) or
*                     photon frequencies in Hertz
*            PC    =  array of total photoionization cross-sections
*                     as a function of orbital angular momentum
*                     and frequency
*            PCP   =  array of partial photoionization cross-sections
*                     as a function of orbital angular momentum
*                     and frequency
*            PNC   =  array of photoionization cross-sections
*                     for level NL as a function of frequncy
*            NDE   =  dimension of PHE and PNC arrays in calling program
*            NDL   =  first dimension of PC array in calling program
*            IOPT  =  option switch:
*                  =  0,  array PHE contains free electron energies
*                  =  1,  array PHE contains photon frequencies
*                  =  2,  initialisation of exponentials and factorials,
*                         input maximum principal quantum number in NE
*            IPCP  =  0,  no PCP needed.
*                  =  1,  PCP needed.

*
*      Storage of photoionization cross-sections for states NL,L
*      is such that
*      PC(L+1,IE) is the total cross-section from NL,L summed over
*                    final states at energy IE
*      PCP(L,1,IE) is the cross-section from NL,L   to L-1, at energy IE
*      PCP(L,2,IE) is the cross-section from NL,L-1 to L  , at energy IE
********
* NOTE * If the partial cross-sections PCP are not required, storage can
******** be reduced by removing the array PCP from the calling sequence
*        of PIXZ1,and from the DIMENSION and DOUBLE PRECISION statements
*        Statements using PCP in PIXZ1 must also be deleted; they are
*        marked by **.  The total cross-section for each bound state
*        NL,L will still be returned in PC.
*
*  Import
      INTEGER NL,NE,IOPT,IPCP,NDE,NDL
      DOUBLE PRECISION Z,AM,PHE
*  Export
      DOUBLE PRECISION PC,PCP,PNC
*  Local
      INTEGER I,IE,L,MAX,MUL,POW
      DOUBLE PRECISION FAL,ALO,CON,KAP,ANL,EU,PCP1,PCP2,
     :                 AI,RP,RM,FP,FM,AL,CU,CL,BU,BL,TWOL,REM,
     :                 TM,TP,FMUL,R0,ZERO,ONE,TWO,TEN,D10,DM10
*
      DIMENSION FAL(1000),ALO(1000),PHE(NDE),PC(NDL,NDE),PCP(NDL,2,NDE),
     :          PNC(NDE)
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, TEN=10.D0,
     :           D10=1.D10, DM10=1.D-10)
*
*            Save factorials and exponentials
*
      SAVE FAL,ALO
      
      AI = 0.0
      EU = 0.0
*
*            Initialise arrays of factorials and exponentials
*
      IF(IOPT.EQ.2) THEN
           MAX=2*NE
           IF(MAX.GT.1000) THEN
             WRITE(6,101)
             IOPT=-1
             RETURN
           ENDIF
           FAL(1) = ZERO
           DO 2 I =2,MAX
              AI = REAL(I-1)
              ALO(I-1) = LOG10(AI)
              FAL(I) = ALO(I-1) + FAL(I-1)
    2      CONTINUE
*
      ELSE
*
*
*          Reduced mass
*
         REM=ONE/(ONE+ONE/(1822.889D0*AM))
         DO 3 IE=1,NE
            PNC(IE)=ZERO
            DO 3 I=1,NL+1
               PC(I,IE)=ZERO
**
               IF (IPCP .EQ. 1) THEN 
                 PCP(I,1,IE)=ZERO
                 PCP(I,2,IE)=ZERO
               ENDIF
**
    3    CONTINUE
*
*              Calculation of radial matrix elements, R
*              ****************************************
*
*              Using the notation R(eu,ll,lu) = R(nl,ll,eu,lu)
*              where eu is the free-electron energy
*
*              Start with direct calculation of R(eu,nl-1,nl)
*              In addition, use                 R(eu,nl,nl-1) = 0
*
         ANL=REAL(NL)
* Need bound-free OS in Atomic Units. change the conversion factor.
* M. F. Gu, 10/25/2001.
*        CON=8.5596557D-19*REM*REM*Z*Z/(ANL*ANL)
         CON = 0.6666666667*REM*REM*Z*Z/(ANL*ANL)

*
*              Evaluate log10 ( R(0,nl-1,nl) ) -
*              the value of the radial integral at threshold
*
         R0=(ANL+TWO)*(0.6020599913279624D0+ALO(NL))-FAL(NL+NL)/TWO
     :          -0.8685889638065036D0*ANL-0.504000052812886D0

*
*              Obtain R(eu,nl-1,nl) from R(0,nl-1,nl)
*
         DO 8 IE=1,NE
*
              IF(IOPT.EQ.0) EU=TWO*PHE(IE)/(REM*Z*Z)
              IF(IOPT.EQ.1) THEN
                  EU=3.0396597D-16*PHE(IE)/(REM*Z*Z)-ONE/(ANL*ANL)
              ENDIF
              IF(EU.LT.ZERO) GO TO 8
              KAP=SQRT(EU)
              FM=ONE
              MUL=0
              DO 4 I=1,NL
                  AI=REAL(I)
                  AI=ONE+AI*AI*EU
                  FM=FM*AI
                  IF(FM.GT.D10) THEN
                     POW=INT(LOG10(FM))
                     FM=FM/TEN**POW
                     MUL=MUL+POW
                  ENDIF
    4         CONTINUE
              FP=ZERO
              IF((ANL*KAP).GT.1.D-20) THEN
                  FP=(ANL-ATAN(ANL*KAP)/KAP)
              ENDIF
              FM = R0 + LOG10(FM)/TWO - (ANL+TWO)*LOG10(AI)
     :                + 0.8685889638065036D0*FP + REAL(MUL)/TWO
              FP=ONE
              IF(KAP.GE.0.1D0) THEN
                  FP=6.283185307179586D0/KAP
                  FP=ONE-EXP(-FP)
                  FP=SQRT(ONE/FP)
              ENDIF
              MUL=INT(FM)
              FM=FM-MUL
              FM=FP*TEN**FM/(REM*REM*Z*Z)
              FP = 0.D0
              CU=SQRT(ONE+ANL*ANL*KAP*KAP)/(ANL)
              CL=ONE
*
*               Recur on R :
*               Use            R(eu,l,l+1) and R(eu,l+1,l)      (rm, rp)
*               to construct   R(eu,l-1,l) and R(eu,l,l-1)      (fm, fp)
*
              DO 5 L=NL,1,-1
                  AL=REAL(L)
                  TWOL=TWO*AL
                  IF(L.EQ.NL) GO TO 6
                  CU=SQRT(ONE+AL*AL*KAP*KAP)/(AL)
                  CL=SQRT((ANL-AL)*(ANL+AL))/(ANL*AL)
                  FP = ( (TWOL+ONE)*BL*RP + BU*RM) / (TWOL*CU)
                  FM = ( BL*RP + (TWOL+ONE)*BU*RM) / (TWOL*CL)
*
*               Calculate photoionization cross-sections
*
    6             TP = TEN**MUL
                  TM = FM*TP
                  TP = FP*TP
                  PCP1 = AL/(TWOL+ONE)*TP*AI*TP*CON
                  PCP2 = AL/(TWOL-ONE)*TM*AI*TM*CON
**
                  IF (IPCP .EQ. 1) THEN 
                     PCP(L,1,IE) = PCP1
                     PCP(L,2,IE) = PCP2
                  ENDIF
**
                  IF(L.LT.NL)  PC(L+1,IE) = PC(L+1,IE) + PCP1
                  PC(L,IE)   = PC(L,IE) + PCP2
                  FMUL=ONE
                  IF(FM.GT.D10) THEN
                      FMUL=DM10
                      MUL=MUL+10
                  ENDIF
                  RP=FP*FMUL
                  RM=FM*FMUL
                  BU=CU
                  BL=CL
    5         CONTINUE
*
*               Calculate total cross-section from nl.
*
              PNC(IE) = 0.D0
              DO 7 L=1,NL
                  AL=REAL(L)
                  PNC(IE) = PNC(IE) + (TWO*AL-ONE)*PC(L,IE)
    7         CONTINUE
              PNC(IE) = PNC(IE)/(ANL*ANL)
    8    CONTINUE
         RETURN
      ENDIF
*
  101 FORMAT(' INSUFFICIENT WORK SPACE IN PIXZ1'/
     :       ' INCREASE DIMENSIONS OF FAL AND ALO TO AT LEAST'/
     :       ' 2*(MAXIMUM PRINCIPAL QUANTUM NUMBER)'/)
*
      END
*
************************************************************************
