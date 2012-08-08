module mod_cu_tiedtke_38r2

  use mod_constants

  integer , parameter :: n_vmass = 0       ! Using or not vector mass
  logical , parameter :: lphylin = .false. ! linearized physics is activated ?
  real(dp) , parameter :: rlpal1 = 0.15D0  ! Smoothing coefficient
  real(dp) , parameter :: rlpal2 = 20.0D0  ! Smoothing coefficient
  real(dp) , parameter :: zqmax = 0.5D0

  contains
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
! SUBROUTINE CUADJTQ
!
!          MODIFICATIONS
!          -------------
!          D.SALMOND         CRAY(UK))      12/8/91
!          J.J. MORCRETTE    ECMWF          92-09-18   Update to Cy44
!          J.F. MAHFOUF      ECMWF          96-06-11   Smoothing option
!          D.SALMOND & M.HAMRUD ECMWF       99-06-04   Optimisation
!          J.HAGUE                          03-01-13   MASS Vector Functions
!          J.HAGUE                          03-07-07   More MASS V.F.
!        M.Hamrud              01-Oct-2003 CY28 Cleaning
!        J.Hague & D.Salmond   22-Nov-2005 Optimisations 
!          PURPOSE.
!          --------
!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM SUBROUTINES:
!              *COND*     (T AND Q AT CONDENSATION LEVEL)
!              *CUBASE*   (T AND Q AT CONDENSATION LEVEL)
!              *CUASC*    (T AND Q AT CLOUD LEVELS)
!              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
!              *CUSTRAT*  (T AND Q AT CONDENSATION LEVEL)
!          INPUT ARE UNADJUSTED T AND Q VALUES,
!          IT RETURNS ADJUSTED VALUES OF T AND Q
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           LEVEL
!    *KCALL*        DEFINES CALCULATION AS
!                      KCALL=0  ENV. T AND QS IN*CUINI*
!                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)
!     INPUT PARAMETERS (REAL):
!    *PSP*          PRESSURE                                        PA
!     UPDATED PARAMETERS (REAL):
!    *PT*           TEMPERATURE                                     K
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG
!          EXTERNALS   
!          ---------
!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
!          FOR CONDENSATION CALCULATIONS.
!          THE TABLES ARE INITIALISED IN *SUPHEC*.
!----------------------------------------------------------------------

  subroutine cuadjtq(kidia,kfdia,klon,ktdia,klev,kk,psp,pt,pq,ldflag,kcall)
    implicit none

    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer :: ktdia ! argument not used
    integer , intent(in) :: kk
    real(dp) , dimension(klon) , intent(in) :: psp
    real(dp) , dimension(klon,klev) , intent(inout) :: pt , pq
    logical , dimension(klon) , intent(in) :: ldflag
    integer , intent(in) :: kcall
    integer :: jl , jlen

    real(dp) , dimension(kfdia-kidia+1) :: ztmp0
    real(dp) , dimension(kfdia-kidia+1) :: ztmp1
    real(dp) , dimension(kfdia-kidia+1) :: ztmp2
    real(dp) , dimension(kfdia-kidia+1) :: ztmp3
    real(dp) , dimension(kfdia-kidia+1) :: ztmp4
    real(dp) , dimension(kfdia-kidia+1) :: ztmp5
    real(dp) , dimension(kfdia-kidia+1) :: ztmp6

    real(dp) :: z1s , z2s , zcond , zcond1 , zcor , zfoeewi , zfoeewl , &
                zoealfa , zqsat , ztarg , zqp
    real(dp) :: zl , zi , zf

!dir$ vfunction exphf

    if ( .not. lphylin ) then
      if ( n_vmass > 0 ) then
        jlen = kfdia-kidia+1
      end if
      !
      ! calculate condensation and adjust t and q accordingly
      !
      if ( kcall == 1 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            ! zqsat = foeewmcu(pt(jl,kk))*zqp
            ! foeewmcu(ptare) = c2es*(foealfcu(ptare)* &
            !                   exp(c3les*(ptare-tzero)/(ptare-c4les))+&
            !                 (d_one-foealfcu(ptare))* &
            !                   exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
            zl = d_one/(pt(jl,kk)-c4les)
            zi = d_one/(pt(jl,kk)-c4ies)
            zqsat = c2es *(foealfcu(pt(jl,kk))*exp(c3les*(pt(jl,kk)-tzero)*zl)+&
                    (d_one-foealfcu(pt(jl,kk)))*exp(c3ies*(pt(jl,kk)-tzero)*zi))
            zqsat = zqsat*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = d_one-retv*zqsat
            ! zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
            !         (zcor**d_two+zqsat*foedemcu(pt(jl,kk)))
            ! foedemcu(ptare) = foealfcu(ptare)*c5alvcp*  &
            !                   (d_one/(ptare-c4les)**d_two)+ &
            !                   (d_one-foealfcu(ptare))*  &
            !                   c5alscp*(d_one/(ptare-c4ies)**d_two)
            zf = foealfcu(pt(jl,kk))*c5alvcp*zl**d_two + &
                 (d_one-foealfcu(pt(jl,kk)))*c5alscp*zi**d_two
            zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor)/(zcor**d_two+zqsat*zf)
            ! zcond = max(zcond,0.0D0)
            if ( zcond > 0.0D0 ) then
              pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond
              pq(jl,kk) = pq(jl,kk)-zcond
              ! zqsat = foeewmcu(pt(jl,kk))*zqp
              zl = d_one/(pt(jl,kk)-c4les)
              zi = d_one/(pt(jl,kk)-c4ies)
              zqsat = c2es*(foealfcu(pt(jl,kk)) * &
                      exp(c3les*(pt(jl,kk)-tzero)*zl) + &
                      (d_one-foealfcu(pt(jl,kk))) * &
                      exp(c3ies*(pt(jl,kk)-tzero)*zi))
              zqsat = zqsat*zqp
              zqsat = minj(0.5D0,zqsat)
              zcor = d_one-retv*zqsat
              ! zcond1 = (pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
              !          (zcor**d_two+zqsat*foedemcu(pt(jl,kk)))
              zf = foealfcu(pt(jl,kk))*c5alvcp*zl**d_two + &
                   (d_one-foealfcu(pt(jl,kk)))*c5alscp*zi**d_two
              zcond1 = (pq(jl,kk)*zcor**d_two-zqsat*zcor)/(zcor**d_two+zqsat*zf)
              if ( dabs(zcond) < dlowval ) zcond1 = 0.0D0
              pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
              pq(jl,kk) = pq(jl,kk)-zcond1
            end if
          end if
        end do
      end if

      if ( kcall == 2 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
            zcond = min(zcond,0.0D0)
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
            if ( dabs(zcond) < dlowval ) zcond1 = min(zcond1,0.0D0)
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end if
        end do
      end if

      if ( kcall == 0 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          zqp = d_one/psp(jl)
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(0.5D0,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(0.5D0,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
        end do
      end if

      if ( kcall == 4 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = d_one/(d_one-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end if
        end do
      end if

      if ( kcall == 5 ) then  ! same as 4 but with ldflag all true
!dir$    ivdep
!ocl novrec
        if ( n_vmass <= 0 )  then ! not using vector mass
          do jl = kidia , kfdia
            zqp = d_one/psp(jl)
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = d_one/(d_one-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = d_one/(d_one-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end do
        else ! using vector vmass
          do jl = kidia , kfdia
            ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
            ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
            ztmp3(jl-kidia+1) = pt(jl,kk)-c4les
            ztmp4(jl-kidia+1) = pt(jl,kk)-c4ies
          end do
          call vdiv(ztmp5,ztmp1,ztmp3,jlen)
          call vdiv(ztmp6,ztmp2,ztmp4,jlen)
          call vexp(ztmp1,ztmp5,jlen)
          call vexp(ztmp2,ztmp6,jlen)
          call vrec(ztmp5,ztmp3,jlen)
          call vrec(ztmp6,ztmp4,jlen)
          do jl = kidia , kfdia
            zqp = d_one/psp(jl)
            zqsat = c2es*(foealfa(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                    (d_one-foealfa(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(0.5D0,zqsat)
            zcor = d_one-retv*zqsat
            ! zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
            !         (zcor**d_two+zqsat*foedem(pt(jl,kk)))
            ! foedem(ptare) = foealfa(ptare)*c5alvcp*&
            !                 (d_one/(ptare-c4les)**d_two)+&
            !                 (d_one-foealfa(ptare))*c5alscp*&
            !                 (d_one/(ptare-c4ies)**d_two)
            zf = foealfa(pt(jl,kk))*c5alvcp*(ztmp5(jl-kidia+1)**d_two) + &
                 (d_one-foealfa(pt(jl,kk)))*c5alscp*(ztmp6(jl-kidia+1)**d_two)
            zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor)/(zcor**d_two+zqsat*zf)
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            ztmp0(jl-kidia+1) = zqp
            ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
            ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
            ztmp3(jl-kidia+1) = pt(jl,kk)-c4les
            ztmp4(jl-kidia+1) = pt(jl,kk)-c4ies
          end do
          call vdiv(ztmp5,ztmp1,ztmp3,jlen)
          call vdiv(ztmp6,ztmp2,ztmp4,jlen)
          call vexp(ztmp1,ztmp5,jlen)
          call vexp(ztmp2,ztmp6,jlen)
          call vrec(ztmp5,ztmp3,jlen)
          call vrec(ztmp6,ztmp4,jlen)
          do jl = kidia , kfdia
            zqp = ztmp0(jl-kidia+1)
            zqsat = c2es*(foealfa(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                    (d_one-foealfa(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(0.5D0,zqsat)
            zcor = d_one-retv*zqsat
            ! zcond1=(pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
            !        (zcor**d_two+zqsat*foedem(pt(jl,kk)))
            ! foedem(ptare) = foealfa(ptare)*c5alvcp*&
            !                 (d_one/(ptare-c4les)**d_two)+&
            !                 (d_one-foealfa(ptare))*c5alscp*&
            !                 (d_one/(ptare-c4ies)**d_two)
            zf = foealfa(pt(jl,kk))*c5alvcp*(ztmp5(jl-kidia+1)**d_two) + &
                 (d_one-foealfa(pt(jl,kk)))*c5alscp*(ztmp6(jl-kidia+1)**d_two)
            zcond1 = (pq(jl,kk)*zcor**d_two-zqsat*zcor)/(zcor**d_two+zqsat*zf)
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end do
        end if
      end if

      if ( kcall == 3 ) then
!dir$    ivdep !ocl novrec
        if ( n_vmass <= 0 ) then ! not using vector mass
          do jl = kidia , kfdia
            zqp = d_one/psp(jl)
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = d_one/(d_one-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end do
        else
          do jl = kidia , kfdia
            ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
            ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
            ztmp3(jl-kidia+1) = pt(jl,kk)-c4les
            ztmp4(jl-kidia+1) = pt(jl,kk)-c4ies
          end do
          call vdiv(ztmp5,ztmp1,ztmp3,jlen)
          call vdiv(ztmp6,ztmp2,ztmp4,jlen)
          call vexp(ztmp1,ztmp5,jlen)
          call vexp(ztmp2,ztmp6,jlen)
          do jl = kidia , kfdia
            zqp = d_one/psp(jl)
            zqsat = c2es*(foealfcu(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                    (d_one-foealfcu(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(0.5D0,zqsat)
            zcor = d_one-retv*zqsat
            zcond1 = (pq(jl,kk)*zcor**d_two - &
                     zqsat*zcor)/(zcor**d_two+zqsat*foedemcu(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
            ztmp0(jl-kidia+1) = zqp
            ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
            ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
            ztmp3(jl-kidia+1) = pt(jl,kk)-c4les
            ztmp4(jl-kidia+1) = pt(jl,kk)-c4ies
          end do
          call vdiv(ztmp5,ztmp1,ztmp3,jlen)
          call vdiv(ztmp6,ztmp2,ztmp4,jlen)
          call vexp(ztmp1,ztmp5,jlen)
          call vexp(ztmp2,ztmp6,jlen)
          do jl = kidia,kfdia
            zqp = ztmp0(jl-kidia+1)
            zqsat = c2es*(foealfcu(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                    (d_one-foealfcu(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(0.5D0,zqsat)
            zcor = d_one-retv*zqsat
            zcond1 = (pq(jl,kk)*zcor**d_two - &
                     zqsat*zcor)/(zcor**d_two+zqsat*foedemcu(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end do
        end if
      end if
    else
    !
    !  calculate condensation and adjust t and q accordingly
    !
      if ( kcall == 1 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            ztarg = pt(jl,kk)
            zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = 0.5D0*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                  (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
            zcond = max(zcond,0.0D0)
            if ( dabs(zcond) > dlowval ) then
              pt(jl,kk) = pt(jl,kk) + &
                      ( zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond
              pq(jl,kk) = pq(jl,kk)-zcond
              ztarg = pt(jl,kk)
              zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
              zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
              zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
              zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
              z1s = tanh(rlpal2*(zqsat-zqmax))
              zqsat = 0.5D0*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
              zcor = d_one/(d_one-retv*zqsat)
              zqsat = zqsat*zcor
              z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                    (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
              zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
              pt(jl,kk) = pt(jl,kk) + &
                         (zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond1
              pq(jl,kk) = pq(jl,kk)-zcond1
            end if
          end if
        end do
      end if

      if ( kcall == 2 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            ztarg = pt(jl,kk)
            zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = 0.5D0*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                  (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
            zcond = min(zcond,0.0D0)
            if ( dabs(zcond) > dlowval ) then
              pt(jl,kk) = pt(jl,kk) + &
                         (zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond
              pq(jl,kk) = pq(jl,kk)-zcond
              ztarg = pt(jl,kk)
              zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
              zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
              zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
              zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
              z1s = tanh(rlpal2*(zqsat-zqmax))
              zqsat = 0.5D0*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
              zcor = d_one/(d_one-retv  *zqsat)
              zqsat = zqsat*zcor
              z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                    (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
              zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
              pt(jl,kk) = pt(jl,kk) + &
                          (zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond1
              pq(jl,kk) = pq(jl,kk)-zcond1
            end if
          end if
        end do
      end if

      if ( kcall == 0 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          zqp = d_one/psp(jl)
          ztarg = pt(jl,kk)
          zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
          zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
          zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
          zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
          z1s = tanh(rlpal2*(zqsat-zqmax))
          zqsat = 0.5D0*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
          zcor = d_one/(d_one-retv  *zqsat)
          zqsat = zqsat*zcor
          z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
          zcond1=(pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
          pt(jl,kk) = pt(jl,kk)+(zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
          ztarg = pt(jl,kk)
          zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
          zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
          zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
          zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
          z1s = tanh(rlpal2*(zqsat-zqmax))
          zqsat = 0.5D0*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
          pt(jl,kk) = pt(jl,kk)+(zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
        end do
      end if

      if ( kcall == 4 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            ztarg = pt(jl,kk)
            zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = 0.5D0*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                  (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
            pt(jl,kk) = pt(jl,kk) + &
                       (zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            ztarg = pt(jl,kk)
            zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = 0.5D0*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
            zqsat = min(zqmax,zqsat)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa*c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                  (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
            pt(jl,kk) = pt(jl,kk) + &
                   (zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond1
            pq(jl,kk) = pq(jl,kk) - zcond1
          end if
        end do
      end if
    end if
  end subroutine cuadjtq

  subroutine satur(kidia,kfdia,klon,ktdia,klev,paprsf,pt,pqsat,kflag)  
    implicit none

    integer , intent(in) :: klon 
    integer , intent(in) :: klev 
    integer , intent(in) :: kidia 
    integer , intent(in) :: kfdia 
    integer , intent(in) :: ktdia 
    real(dp) , dimension(klon,klev) , intent(in) :: paprsf , pt
    real(dp) , dimension(klon,klev) , intent(out) :: pqsat
    integer , intent(in) :: kflag 
    integer :: jk , jl , jlen

    real(dp) :: zcor , zew , zfoeew , zqmax , zqs , ztarg
    real(dp) :: zalfa , zfoeewl , zfoeewi
    real(dp) , dimension(kidia:kfdia) :: z_exparg1
    real(dp) , dimension(kidia:kfdia) :: z_exparg2
    real(dp) , dimension(kidia:kfdia) :: z_expout1
    real(dp) , dimension(kidia:kfdia) :: z_expout2

!dir$ vfunction exphf

    !
    ! calculate saturation specific humidity
    !
    if ( lphylin ) then
      do jk = ktdia , klev
        do jl = kidia , kfdia
          ztarg = pt(jl,jk)
          zalfa = foealfa(ztarg)
          zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
          zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
          zfoeew = zalfa*zfoeewl+(d_one-zalfa)*zfoeewi
          zqs = zfoeew/paprsf(jl,jk)
          if ( zqs > zqmax ) then
            zqs = zqmax
          end if
          zcor = d_one/(d_one-retv*zqs)
          pqsat(jl,jk) = zqs*zcor
        end do
      end do
    else
      if ( n_vmass <= 0 ) then ! not using vector mass
        do jk = ktdia , klev
          do jl = kidia , kfdia
            if ( kflag == 1 ) then
              zew = foeewmcu(pt(jl,jk))
            else
              zew = foeewm(pt(jl,jk))
            end if
            zqs = zew/paprsf(jl,jk)
            zqs = min(zqmax,zqs)
            zcor = d_one/(d_one-retv*zqs)
            pqsat(jl,jk) = zqs*zcor
          end do
        end do
      else ! using vector mass
        jlen = kfdia-kidia+1
        do jk = ktdia , klev
          do jl = kidia , kfdia
            z_exparg1(jl) = foeles_v(pt(jl,jk))
            z_exparg2(jl) = foeies_v(pt(jl,jk))
          end do
          call vexp(z_expout1,z_exparg1,jlen)
          call vexp(z_expout2,z_exparg2,jlen)
          do jl = kidia , kfdia
            if ( kflag == 1 ) then
              zew = foeewmcu_v( pt(jl,jk),z_expout1(jl),z_expout2(jl) )
            else
              zew = foeewm_v( pt(jl,jk),z_expout1(jl),z_expout2(jl) )
            end if
            ! zqs = zew/paprsf(jl,jk)
            ! zqs = min(zqmax,zqs)
            !! zcor = _one_/(_one_-retv*zqs)
            ! pqsat(jl,jk) = zqs/(_one_-retv*zqs)
            zqs = min(zqmax*paprsf(jl,jk),zew)
            pqsat(jl,jk) = zqs/(paprsf(jl,jk)-retv*zqs)
          end do
        end do
      end if
    end if
  end subroutine satur

  real(dp) function minj(x,y)
    implicit none
    real(dp) , intent(in) :: x , y
    minj = y - 0.5D0*(dabs(x-y)-(x-y))
  end function minj

  real(dp) function maxj(x,y)
    implicit none
    real(dp) , intent(in) :: x , y
    maxj = y + 0.5D0*(dabs(x-y)+(x-y))
  end function maxj

  real(dp) function foedelta(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedelta = max(0.0D0,sign(d_one,ptare-tzero))
  end function foedelta
  real(dp) function foeew(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeew = c2es*exp((c3les*foedelta(ptare) + &
            c3ies*(d_one-foedelta(ptare)))*(ptare-tzero) / &
            (ptare-(c4les*foedelta(ptare)+c4ies*(d_one-foedelta(ptare)))))
  end function foeew
  real(dp) function foede(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foede = (foedelta(ptare)*c5alvcp+(d_one-foedelta(ptare))*c5alscp) / &
       (ptare-(c4les*foedelta(ptare)+c4ies*(d_one-foedelta(ptare))))**d_two
  end function foede
  real(dp) function foedesu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedesu = (foedelta(ptare)*c5les+(d_one-foedelta(ptare))*c5ies) / &
         (ptare-(c4les*foedelta(ptare)+c4ies*(d_one-foedelta(ptare))))**d_two
  end function foedesu
  real(dp) function foelh(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foelh = foedelta(ptare)*wlhv + (d_one-foedelta(ptare))*wlhs
  end function foelh
  real(dp) function foeldcp(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeldcp = foedelta(ptare)*wlhvocp + (d_one-foedelta(ptare))*wlhsocp
  end function foeldcp
  real(dp) function foealfa(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foealfa = min(d_one,((max(rtice,min(rtwat,ptare))-rtice) * &
                  rtwat_rtice_r)**d_two)
  end function foealfa
  real(dp) function foeewm(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeewm = c2es*(foealfa(ptare)*exp(c3les*(ptare-tzero)/(ptare-c4les))+ &
          (d_one-foealfa(ptare))*exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
  end function foeewm
  real(dp) function foedem(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedem = foealfa(ptare)*c5alvcp*(d_one/(ptare-c4les)**d_two) + &
            (d_one-foealfa(ptare))*c5alscp*(d_one/(ptare-c4ies)**d_two)
  end function foedem
  real(dp) function foeldcpm(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeldcpm = foealfa(ptare)*wlhvocp+(d_one-foealfa(ptare))*wlhsocp
  end function foeldcpm
  real(dp) function foelhm(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foelhm = foealfa(ptare)*wlhv+(d_one-foealfa(ptare))*wlhs
  end function foelhm
  real(dp) function foetb(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foetb = foealfa(ptare)*c3les*(tzero-c4les)*(d_one/(ptare-c4les)**d_two)+ &
      (d_one-foealfa(ptare))*c3ies*(tzero-c4ies)*(d_one/(ptare-c4ies)**d_two)
  end function foetb
  real(dp) function foealfcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foealfcu = min(d_one, &
           ((max(rtice,min(rtwat,ptare))-rtice)*rtwat_rtice_r)**d_two)
  end function foealfcu
  real(dp) function foeewmcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeewmcu = c2es*(foealfcu(ptare)*exp(c3les*(ptare-tzero)/(ptare-c4les))+ &
            (d_one-foealfcu(ptare))*exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
  end function foeewmcu
  real(dp) function foedemcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedemcu = foealfcu(ptare)*c5alvcp*(d_one/(ptare-c4les)**d_two) + &
           (d_one-foealfcu(ptare))*c5alscp*(d_one/(ptare-c4ies)**d_two)
  end function foedemcu
  real(dp) function foeldcpmcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeldcpmcu = foealfcu(ptare)*wlhvocp+(d_one-foealfcu(ptare))*wlhsocp
  end function foeldcpmcu
  real(dp) function foelhmcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foelhmcu = foealfcu(ptare)*wlhv+(d_one-foealfcu(ptare))*wlhs
  end function foelhmcu
  real(dp) function foeewmo(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeewmo = c2es*exp(c3les*(ptare-tzero)/(ptare-c4les))
  end function foeewmo
  real(dp) function foeeliq(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeeliq = c2es*exp(c3les*(ptare-tzero)/(ptare-c4les))
  end function foeeliq
  real(dp) function foeeice(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeeice = c2es*exp(c3ies*(ptare-tzero)/(ptare-c4ies))
  end function foeeice
  real(dp) function foeles_v(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeles_v = c3les*(ptare-tzero)/(ptare-c4les)
  end function foeles_v
  real(dp) function foeies_v(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeies_v = c3ies*(ptare-tzero)/(ptare-c4ies)
  end function foeies_v
  real(dp) function foeewm_v(ptare,exp1,exp2)
    implicit none
    real(dp) , intent(in) :: ptare , exp1 , exp2
    foeewm_v = c2es*(foealfa(ptare)*exp1+(d_one-foealfa(ptare))*exp2)
  end function foeewm_v
  real(dp) function foeewmcu_v(ptare,exp1,exp2)
    implicit none
    real(dp) , intent(in) :: ptare , exp1 , exp2
    foeewmcu_v = c2es*(foealfcu(ptare)*exp1+(d_one-foealfcu(ptare))*exp2)
  end function foeewmcu_v
  subroutine vdiv(z,x,y,n)
    implicit none
    real(dp) , dimension(:) , intent(out) :: z
    real(dp) , dimension(:) , intent(in) :: x , y
    integer , intent(in) :: n
    integer :: i
    do i = 1 , n
      z(i) = x(i)/y(i)
    end do
  end subroutine vdiv
  subroutine vexp(y,x,n)
    implicit none
    real(dp) , dimension(:) , intent(out) :: y
    real(dp) , dimension(:) , intent(in) :: x
    integer , intent(in) :: n
    integer :: i
    do i = 1 , n
      y(i) = dexp(x(i))
    end do
  end subroutine vexp
  subroutine vrec(y,x,n)
    implicit none
    real(dp) , dimension(:) , intent(out) :: y
    real(dp) , dimension(:) , intent(in) :: x
    integer , intent(in) :: n
    integer :: i
    do i = 1 , n
      y(i) = d_one/x(i)
    end do
  end subroutine vrec


end module mod_cu_tiedtke_38r2
