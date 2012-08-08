module mod_cu_tiedte_38r2

  use mod_constants

  integer , parameter :: n_vmass = 0       ! Using or not vector mass
  logical , parameter :: lphylin = .false. ! linearized physics is activated ?
  real(dp) , parameter :: rlpal1 = 0.15D0  ! Smoothing coefficient
  real(dp) , parameter :: rlpal2 = 20.0D0  ! Smoothing coefficient
  real(dp) , parameter :: zqmax = 0.5D0

  contains

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
            zqp = 1.0D0/psp(jl)
            ! zqsat = foeewmcu(pt(jl,kk))*zqp
            ! foeewmcu(ptare) = c2es*(foealfcu(ptare)* &
            !                   exp(c3les*(ptare-tzero)/(ptare-c4les))+&
            !                 (1.0D0-foealfcu(ptare))* &
            !                   exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
            zl = 1.0D0/(pt(jl,kk)-c4les)
            zi = 1.0D0/(pt(jl,kk)-c4ies)
            zqsat = c2es *(foealfcu(pt(jl,kk))*exp(c3les*(pt(jl,kk)-tzero)*zl)+&
                    (1.0D0-foealfcu(pt(jl,kk)))*exp(c3ies*(pt(jl,kk)-tzero)*zi))
            zqsat = zqsat*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = 1.0D0-retv*zqsat
            ! zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
            !         (zcor**d_two+zqsat*foedemcu(pt(jl,kk)))
            ! foedemcu(ptare) = foealfcu(ptare)*c5alvcp*  &
            !                   (1.0D0/(ptare-c4les)**d_two)+ &
            !                   (1.0D0-foealfcu(ptare))*  &
            !                   c5alscp*(1.0D0/(ptare-c4ies)**d_two)
            zf = foealfcu(pt(jl,kk))*c5alvcp*zl**d_two + &
                 (1.0D0-foealfcu(pt(jl,kk)))*c5alscp*zi**d_two
            zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor)/(zcor**d_two+zqsat*zf)
            ! zcond = max(zcond,0.0D0)
            if ( zcond > 0.0D0 ) then
              pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond
              pq(jl,kk) = pq(jl,kk)-zcond
              ! zqsat = foeewmcu(pt(jl,kk))*zqp
              zl = 1.0D0/(pt(jl,kk)-c4les)
              zi = 1.0D0/(pt(jl,kk)-c4ies)
              zqsat = c2es*(foealfcu(pt(jl,kk)) * &
                      exp(c3les*(pt(jl,kk)-tzero)*zl) + &
                      (1.0D0-foealfcu(pt(jl,kk))) * &
                      exp(c3ies*(pt(jl,kk)-tzero)*zi))
              zqsat = zqsat*zqp
              zqsat = minj(0.5D0,zqsat)
              zcor = 1.0D0-retv*zqsat
              ! zcond1 = (pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
              !          (zcor**d_two+zqsat*foedemcu(pt(jl,kk)))
              zf = foealfcu(pt(jl,kk))*c5alvcp*zl**d_two + &
                   (1.0D0-foealfcu(pt(jl,kk)))*c5alscp*zi**d_two
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
            zqp = 1.0D0/psp(jl)
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = 1.0D0/(1.0D0-retv*zqsat)
            zqsat = zqsat*zcor
            zcond = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*foedemcu(pt(jl,kk)))
            zcond = min(zcond,0.0D0)
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = 1.0D0/(1.0D0-retv*zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*foedemcu(pt(jl,kk)))
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
          zqp = 1.0D0/psp(jl)
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(0.5D0,zqsat)
          zcor = 1.0D0/(1.0D0-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*foedem(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(0.5D0,zqsat)
          zcor = 1.0D0/(1.0D0-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*foedem(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
        end do
      end if

      if ( kcall == 4 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = 1.0D0/psp(jl)
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = 1.0D0/(1.0D0-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*foedem(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = 1.0D0/(1.0D0-retv*zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*foedem(pt(jl,kk)))
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
            zqp = 1.0D0/psp(jl)
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = 1.0D0/(1.0D0-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*foedem(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = 1.0D0/(1.0D0-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*foedem(pt(jl,kk)))
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
            zqp = 1.0D0/psp(jl)
            zqsat = c2es*(foealfa(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                    (1.0D0-foealfa(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(0.5D0,zqsat)
            zcor = 1.0D0-retv*zqsat
            ! zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
            !         (zcor**d_two+zqsat*foedem(pt(jl,kk)))
            ! foedem(ptare) = foealfa(ptare)*c5alvcp*&
            !                 (1.0D0/(ptare-c4les)**d_two)+&
            !                 (1.0D0-foealfa(ptare))*c5alscp*&
            !                 (1.0D0/(ptare-c4ies)**d_two)
            zf = foealfa(pt(jl,kk))*c5alvcp*(ztmp5(jl-kidia+1)**d_two) + &
                 (1.0D0-foealfa(pt(jl,kk)))*c5alscp*(ztmp6(jl-kidia+1)**d_two)
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
                    (1.0D0-foealfa(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(0.5D0,zqsat)
            zcor = 1.0D0-retv*zqsat
            ! zcond1=(pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
            !        (zcor**d_two+zqsat*foedem(pt(jl,kk)))
            ! foedem(ptare) = foealfa(ptare)*c5alvcp*&
            !                 (1.0D0/(ptare-c4les)**d_two)+&
            !                 (1.0D0-foealfa(ptare))*c5alscp*&
            !                 (1.0D0/(ptare-c4ies)**d_two)
            zf = foealfa(pt(jl,kk))*c5alvcp*(ztmp5(jl-kidia+1)**d_two) + &
                 (1.0D0-foealfa(pt(jl,kk)))*c5alscp*(ztmp6(jl-kidia+1)**d_two)
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
            zqp = 1.0D0/psp(jl)
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = 1.0D0/(1.0D0-retv*zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*foedemcu(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(0.5D0,zqsat)
            zcor = 1.0D0/(1.0D0-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*foedemcu(pt(jl,kk)))
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
            zqp = 1.0D0/psp(jl)
            zqsat = c2es*(foealfcu(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                    (1.0D0-foealfcu(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(0.5D0,zqsat)
            zcor = 1.0D0-retv*zqsat
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
                    (1.0D0-foealfcu(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(0.5D0,zqsat)
            zcor = 1.0D0-retv*zqsat
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
            zqp = 1.0D0/psp(jl)
            ztarg = pt(jl,kk)
            zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+1.0D0)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(1.0D0-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = 0.5D0*((1.0D0-z1s)*zqsat+(1.0D0+z1s)*zqmax)
            zcor = 1.0D0/(1.0D0-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa *c5alvcp*(1.0D0/(ztarg-c4les)**d_two) + &
                  (1.0D0-zoealfa)*c5alscp*(1.0D0/(ztarg-c4ies)**d_two)
            zcond = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*z2s)
            zcond = max(zcond,0.0D0)
            if ( dabs(zcond) > dlowval ) then
              pt(jl,kk) = pt(jl,kk) + &
                      ( zoealfa*wlhvocp+(1.0D0-zoealfa)*wlhsocp)*zcond
              pq(jl,kk) = pq(jl,kk)-zcond
              ztarg = pt(jl,kk)
              zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+1.0D0)
              zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
              zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
              zqsat = zqp*(zoealfa*zfoeewl+(1.0D0-zoealfa)*zfoeewi)
              z1s = tanh(rlpal2*(zqsat-zqmax))
              zqsat = 0.5D0*((1.0D0-z1s)*zqsat+(1.0D0+z1s)*zqmax)
              zcor = 1.0D0/(1.0D0-retv*zqsat)
              zqsat = zqsat*zcor
              z2s = zoealfa *c5alvcp*(1.0D0/(ztarg-c4les)**d_two) + &
                    (1.0D0-zoealfa)*c5alscp*(1.0D0/(ztarg-c4ies)**d_two)
              zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*z2s)
              pt(jl,kk) = pt(jl,kk) + &
                         (zoealfa*wlhvocp+(1.0D0-zoealfa)*wlhsocp)*zcond1
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
            zqp = 1.0D0/psp(jl)
            ztarg = pt(jl,kk)
            zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+1.0D0)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(1.0D0-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = 0.5D0*((1.0D0-z1s)*zqsat+(1.0D0+z1s)*zqmax)
            zcor = 1.0D0/(1.0D0-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa *c5alvcp*(1.0D0/(ztarg-c4les)**d_two) + &
                  (1.0D0-zoealfa)*c5alscp*(1.0D0/(ztarg-c4ies)**d_two)
            zcond = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*z2s)
            zcond = min(zcond,0.0D0)
            if ( dabs(zcond) > dlowval ) then
              pt(jl,kk) = pt(jl,kk) + &
                         (zoealfa*wlhvocp+(1.0D0-zoealfa)*wlhsocp)*zcond
              pq(jl,kk) = pq(jl,kk)-zcond
              ztarg = pt(jl,kk)
              zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+1.0D0)
              zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
              zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
              zqsat = zqp*(zoealfa*zfoeewl+(1.0D0-zoealfa)*zfoeewi)
              z1s = tanh(rlpal2*(zqsat-zqmax))
              zqsat = 0.5D0*((1.0D0-z1s)*zqsat+(1.0D0+z1s)*zqmax)
              zcor = 1.0D0/(1.0D0-retv  *zqsat)
              zqsat = zqsat*zcor
              z2s = zoealfa *c5alvcp*(1.0D0/(ztarg-c4les)**d_two) + &
                    (1.0D0-zoealfa)*c5alscp*(1.0D0/(ztarg-c4ies)**d_two)
              zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*z2s)
              pt(jl,kk) = pt(jl,kk) + &
                          (zoealfa*wlhvocp+(1.0D0-zoealfa)*wlhsocp)*zcond1
              pq(jl,kk) = pq(jl,kk)-zcond1
            end if
          end if
        end do
      end if

      if ( kcall == 0 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          zqp = 1.0D0/psp(jl)
          ztarg = pt(jl,kk)
          zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+1.0D0)
          zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
          zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
          zqsat = zqp*(zoealfa*zfoeewl+(1.0D0-zoealfa)*zfoeewi)
          z1s = tanh(rlpal2*(zqsat-zqmax))
          zqsat = 0.5D0*((1.0D0-z1s)*zqsat+(1.0D0+z1s)*zqmax)
          zcor = 1.0D0/(1.0D0-retv  *zqsat)
          zqsat = zqsat*zcor
          z2s = zoealfa *c5alvcp*(1.0D0/(ztarg-c4les)**d_two) + &
                (1.0D0-zoealfa)*c5alscp*(1.0D0/(ztarg-c4ies)**d_two)
          zcond1=(pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*z2s)
          pt(jl,kk) = pt(jl,kk)+(zoealfa*wlhvocp+(1.0D0-zoealfa)*wlhsocp)*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
          ztarg = pt(jl,kk)
          zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+1.0D0)
          zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
          zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
          zqsat = zqp*(zoealfa*zfoeewl+(1.0D0-zoealfa)*zfoeewi)
          z1s = tanh(rlpal2*(zqsat-zqmax))
          zqsat = 0.5D0*((1.0D0-z1s)*zqsat+(1.0D0+z1s)*zqmax)
          zcor = 1.0D0/(1.0D0-retv*zqsat)
          zqsat = zqsat*zcor
          z2s = zoealfa *c5alvcp*(1.0D0/(ztarg-c4les)**d_two) + &
                (1.0D0-zoealfa)*c5alscp*(1.0D0/(ztarg-c4ies)**d_two)
          zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*z2s)
          pt(jl,kk) = pt(jl,kk)+(zoealfa*wlhvocp+(1.0D0-zoealfa)*wlhsocp)*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
        end do
      end if

      if ( kcall == 4 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = 1.0D0/psp(jl)
            ztarg = pt(jl,kk)
            zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+1.0D0)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(1.0D0-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = 0.5D0*((1.0D0-z1s)*zqsat+(1.0D0+z1s)*zqmax)
            zcor = 1.0D0/(1.0D0-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa *c5alvcp*(1.0D0/(ztarg-c4les)**d_two) + &
                  (1.0D0-zoealfa)*c5alscp*(1.0D0/(ztarg-c4ies)**d_two)
            zcond = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*z2s)
            pt(jl,kk) = pt(jl,kk) + &
                       (zoealfa*wlhvocp+(1.0D0-zoealfa)*wlhsocp)*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            ztarg = pt(jl,kk)
            zoealfa = 0.5D0*(tanh(rlpal1*(ztarg-mpcrt))+1.0D0)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(1.0D0-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = 0.5D0*((1.0D0-z1s)*zqsat+(1.0D0+z1s)*zqmax)
            zqsat = min(zqmax,zqsat)
            zcor = 1.0D0/(1.0D0-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa*c5alvcp*(1.0D0/(ztarg-c4les)**d_two) + &
                  (1.0D0-zoealfa)*c5alscp*(1.0D0/(ztarg-c4ies)**d_two)
            zcond1 = (pq(jl,kk)-zqsat)/(1.0D0+zqsat*zcor*z2s)
            pt(jl,kk) = pt(jl,kk) + &
                   (zoealfa*wlhvocp+(1.0D0-zoealfa)*wlhsocp)*zcond1
            pq(jl,kk) = pq(jl,kk) - zcond1
          end if
        end do
      end if
    end if
  end subroutine cuadjtq

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
    foedelta = max(0.0D0,sign(1.0D0,ptare-tzero))
  end function foedelta
  real(dp) function foeew(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeew = c2es*exp((c3les*foedelta(ptare) + &
            c3ies*(1.0D0-foedelta(ptare)))*(ptare-tzero) / &
            (ptare-(c4les*foedelta(ptare)+c4ies*(1.0D0-foedelta(ptare)))))
  end function foeew
  real(dp) function foede(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foede = (foedelta(ptare)*c5alvcp+(1.0D0-foedelta(ptare))*c5alscp) / &
       (ptare-(c4les*foedelta(ptare)+c4ies*(1.0D0-foedelta(ptare))))**d_two
  end function foede
  real(dp) function foedesu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedesu = (foedelta(ptare)*c5les+(1.0D0-foedelta(ptare))*c5ies) / &
         (ptare-(c4les*foedelta(ptare)+c4ies*(1.0D0-foedelta(ptare))))**d_two
  end function foedesu
  real(dp) function foelh(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foelh = foedelta(ptare)*wlhv + (1.0D0-foedelta(ptare))*wlhs
  end function foelh
  real(dp) function foeldcp(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeldcp = foedelta(ptare)*wlhvocp + (1.0D0-foedelta(ptare))*wlhsocp
  end function foeldcp
  real(dp) function foealfa(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foealfa = min(1.0D0,((max(rtice,min(rtwat,ptare))-rtice) * &
                  rtwat_rtice_r)**d_two)
  end function foealfa
  real(dp) function foeewm(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeewm = c2es*(foealfa(ptare)*exp(c3les*(ptare-tzero)/(ptare-c4les))+ &
          (1.0D0-foealfa(ptare))*exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
  end function foeewm
  real(dp) function foedem(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedem = foealfa(ptare)*c5alvcp*(1.0D0/(ptare-c4les)**d_two) + &
            (1.0D0-foealfa(ptare))*c5alscp*(1.0D0/(ptare-c4ies)**d_two)
  end function foedem
  real(dp) function foeldcpm(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeldcpm = foealfa(ptare)*wlhvocp+(1.0D0-foealfa(ptare))*wlhsocp
  end function foeldcpm
  real(dp) function foelhm(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foelhm = foealfa(ptare)*wlhv+(1.0D0-foealfa(ptare))*wlhs
  end function foelhm
  real(dp) function foetb(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foetb = foealfa(ptare)*c3les*(tzero-c4les)*(1.0D0/(ptare-c4les)**d_two)+ &
      (1.0D0-foealfa(ptare))*c3ies*(tzero-c4ies)*(1.0D0/(ptare-c4ies)**d_two)
  end function foetb
  real(dp) function foealfcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foealfcu = min(1.0D0, &
           ((max(rtice,min(rtwat,ptare))-rtice)*rtwat_rtice_r)**d_two)
  end function foealfcu
  real(dp) function foeewmcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeewmcu = c2es*(foealfcu(ptare)*exp(c3les*(ptare-tzero)/(ptare-c4les))+ &
            (1.0D0-foealfcu(ptare))*exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
  end function foeewmcu
  real(dp) function foedemcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedemcu = foealfcu(ptare)*c5alvcp*(1.0D0/(ptare-c4les)**d_two) + &
           (1.0D0-foealfcu(ptare))*c5alscp*(1.0D0/(ptare-c4ies)**d_two)
  end function foedemcu
  real(dp) function foeldcpmcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeldcpmcu = foealfcu(ptare)*wlhvocp+(1.0D0-foealfcu(ptare))*wlhsocp
  end function foeldcpmcu
  real(dp) function foelhmcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foelhmcu = foealfcu(ptare)*wlhv+(1.0D0-foealfcu(ptare))*wlhs
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
    foeewm_v = c2es*(foealfa(ptare)*exp1+(1.0D0-foealfa(ptare))*exp2)
  end function foeewm_v
  real(dp) function foeewmcu_v(ptare,exp1,exp2)
    implicit none
    real(dp) , intent(in) :: ptare , exp1 , exp2
    foeewmcu_v = c2es*(foealfcu(ptare)*exp1+(1.0D0-foealfcu(ptare))*exp2)
  end function foeewmcu_v

end module mod_cu_tiedte_38r2
