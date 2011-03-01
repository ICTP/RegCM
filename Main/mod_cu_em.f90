!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      module mod_cu_em
!
! Kerry Emanuel Convective scheme
!
      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_main
      use mod_pmoist
      use mod_cvaria
      use mod_slice
      use mod_rad
      use mod_bats
      use mod_date
!
      private
!
      public :: cupemandrv
      public :: minsig , elcrit , tlcrit , entp , sigd , sigs ,    &
                omtrain , omtsnow , coeffr , coeffs , cu , betae , &
                dtmax , alphae , damp , minorig
!
      real(8) :: alphae , betae , coeffr , coeffs , cu , damp , dtmax , &
               & elcrit , entp , minsig , omtrain , omtsnow , sigd ,    &
               & sigs , tlcrit
!
      integer :: minorig
!
      contains
!
! **********************************************
! **** Driver for Emanuel Convection Scheme ****
! **********************************************
!
      subroutine cupemandrv(j)
! 
      implicit none
!
      integer :: j
      intent (in) j
!
      integer , parameter :: ntra = 0
!
      real(8) :: akclth , aprdiv , cbmf , dtime , fppt , qprime ,       &
               & tprime , uconv , wd
      real(8) , dimension(kz) :: fq , ft , fu , fv , pcup , qcup ,      &
                               & qscup , tcup , ucup , vcup
      real(8) , dimension(kz,1) :: ftra , tra
      integer :: i , iconj , iflag , k , kbase , kclth , kk , ktop
      real(8) , dimension(kzp1) :: phcup
!
      dtime = dt
      uconv = 0.5D0*dt
      aprdiv = 1.0D0/dble(nbatst)
      if ( jyear == jyear0 .and. ktau == 0 ) aprdiv = 1.0D0
      iconj = 0
      do i = 2 , iym2
        if ( icup==99 .or. icup==98 ) then
          if (cumcon%cuscheme(i,j)/=4 ) cycle
        end if
        do k = 1 , kz
          kk = kzp1 - k
          cldlwc(i,k) = 0.0D0       ! Zero out cloud water content
          cldfra(i,k) = 0.0D0       ! Zero out cloud fraction coverage
          tcup(k) = tb3d(i,kk,j)                          ! [k]
          qcup(k) = qvb3d(i,kk,j)/(1.0D0+qvb3d(i,kk,j))   ! [kg/kg]
          qscup(k) = qsb3d(i,kk,j)/(1.0D0+qsb3d(i,kk,j))  ! [kg/kg]
          ucup(k) = ubx3d(i,kk,j)                         ! [m/s]
          vcup(k) = vbx3d(i,kk,j)                         ! [m/s]
          pcup(k) = pb3d(i,kk,j)*10.0D0                   ! [hPa]
        end do
        do k = 1 , kzp1
          kk = kzp1 - k + 1
          phcup(k) = (sigma(kk)*sps2%ps(i,j)+r8pt)*10.0D0 ! [hPa]
        end do
        cbmf = cbmf2d(i,j)                                ! [(kg/m**2)/s]
 
        call cupeman(tcup,qcup,qscup,ucup,vcup,tra,pcup,phcup,kz,kzp1,  &
                   & kzm1,ntra,dtime,iflag,ft,fq,fu,fv,ftra,fppt,wd,    &
                   & tprime,qprime,cbmf,kbase,ktop)
 
        cbmf2d(i,j) = cbmf
 
!       iflag=0: No moist convection; atmosphere stable or surface
!       temperature < 250K or surface humidity is negative.
!       iflag=1: Moist convection occurs.
!       iflag=2: No moist convection: lifted condensation level above
!       200 mb. iflag=3: No moist convection: cloud base higher than
!       level kzm2. iflag=4: Moist convection occurs, but CFL condition
!       on the subsidence warming is violated. (Does not terminate
!       scheme.)
        if ( iflag == 1 .or. iflag == 4 ) then
                                              ! If moist convection
 
!         if (iflag == 4) then                ! If CFL violation
!         print*,'EMAN CFL VIOLATION: ',i,j,cbmf
!         end if
 
!         **** Tendencies
          do k = 1 , kz
            kk = kzp1 - k
            aten%t(i,kk,j) = ft(k)*sps2%ps(i,j) + aten%t(i,kk,j)
            aten%qv(i,kk,j) = fq(k)/(1.0D0-fq(k))* &
                              sps2%ps(i,j)+aten%qv(i,kk,j)
!           There is a bit of an inconsistency here...  The wind
!           tendencies from convection are on cross points, but the
!           model wants them on dot points.
            aten%u(i,kk,j) = fu(k)*sps2%ps(i,j) + aten%u(i,kk,j)
            aten%v(i,kk,j) = fv(k)*sps2%ps(i,j) + aten%v(i,kk,j)
          end do
 
!         **** Cloud fraction and cloud water
          kclth = ktop - kbase + 1
          akclth = 1.0D0/dble(kclth)
          do k = kbase , ktop
            kk = kzp1 - k
            cldlwc(i,kk) = cllwcv
            cldfra(i,kk) = 1.0D0 - (1.0D0-clfrcv)**akclth
          end do
 
!         **** Precipitation
          if ( fppt > 0.0D0 ) then
            sfsta%rainc(i,j) = sfsta%rainc(i,j) + fppt*uconv ! mm
            pptc(i,j)        = pptc(i,j) + fppt*aprdiv       ! mm/s
            iconj = iconj + 1
          end if
 
        end if
 
      end do
 
      icon(j) = iconj
! 
      end subroutine cupemandrv
!
!**************************************************************************
!****                       subroutine cupeman (formerly convect)     *****
!****                          version 4.3C                           *****
!****                          20 May, 2002                           *****
!****                          Kerry Emanuel                          *****
!**************************************************************************
!
!-----------------------------------------------------------------------------
!    *** on input:      ***
!
!     t:   array of absolute temperature (k) of dimension nd, with first
!           index corresponding to lowest model level. note that this array
!           will be altered by the subroutine if dry convective adjustment
!           occurs and if ipbl is not equal to 0.
!
!     q:   array of specific humidity (gm/gm) of dimension nd, with first
!            index corresponding to lowest model level. must be defined
!            at same grid levels as t. note that this array will be altered
!            if dry convective adjustment occurs and if ipbl is not equal to 0.
!
!     qs:  array of saturation specific humidity of dimension nd, with first
!            index corresponding to lowest model level. must be defined
!            at same grid levels as t. note that this array will be altered
!            if dry convective adjustment occurs and if ipbl is not equal to 0.
!
!     u:   array of zonal wind velocity (m/s) of dimension nd, witth first
!            index corresponding with the lowest model level. defined at
!            same levels as t. note that this array will be altered if
!            dry convective adjustment occurs and if ipbl is not equal to 0.
!
!     v:   same as u but for meridional velocity.
!
!     tra: array of passive tracer mixing ratio, of dimensions (nd,ntra),
!            where ntra is the number of different tracers. if no
!            convective tracer transport is needed, define a dummy
!            input array of dimension (nd,1). tracers are defined at
!            same vertical levels as t. note that this array will be altered
!            if dry convective adjustment occurs and if ipbl is not equal to 0.
!
!     p:   array of pressure (mb) of dimension nd, with first
!            index corresponding to lowest model level. must be defined
!            at same grid levels as t.
!
!     ph:  array of pressure (mb) of dimension nd+1, with first index
!            corresponding to lowest level. these pressures are defined at
!            levels intermediate between those of p, t, q and qs. the first
!            value of ph should be greater than (i.e. at a lower level than)
!            the first value of the array p.
!
!     nd:  the dimension of the arrays t,q,qs,p,ph,ft and fq
!
!     nl:  the maximum number of levels to which convection can
!            penetrate, plus 1.
!            nl must be less than or equal to nd-1.
!
!     ntra:the number of different tracers. if no tracer transport
!            is needed, set this equal to 1. (on most compilers, setting
!            ntra to 0 will bypass tracer calculation, saving some cpu.)
!
!     delt: the model time step (sec) between calls to convect
!
!----------------------------------------------------------------------------
!    ***   on output:         ***
!
!     iflag: an output integer whose value denotes the following:
!
!                value                        interpretation
!                -----                        --------------
!                  0               no moist convection; atmosphere is not
!                                  unstable, or surface temperature is less
!                                  than 250 k or surface specific humidity
!                                  is non-positive.
!
!                  1               moist convection occurs.
!
!                  2               no moist convection: lifted condensation
!                                  level is above the 200 mb level.
!
!                  3               no moist convection: cloud base is higher
!                                  then the level nl-1.
!
!                  4               moist convection occurs, but a cfl condition
!                                  on the subsidence warming is violated. this
!                                  does not cause mod_the scheme to terminate.
!
!     ft:   array of temperature tendency (k/s) of dimension nd, defined at same
!             grid levels as t, q, qs and p.
!
!     fq:   array of specific humidity tendencies ((gm/gm)/s) of dimension nd,
!             defined at same grid levels as t, q, qs and p.
!
!     fu:   array of forcing of zonal velocity (m/s^2) of dimension nd,
!             defined at same grid levels as t.
!
!     fv:   same as fu, but for forcing of meridional velocity.
!
!     ftra: array of forcing of tracer content, in tracer mixing ratio per
!             second, defined at same levels as t. dimensioned (nd,ntra).
!
!     precip: scalar convective precipitation rate (mm/s).
!
!     wd:    a convective downdraft velocity scale. for use in surface
!             flux parameterizations. see convect.ps file for details.
!
!     tprime: a convective downdraft temperature perturbation scale (k).
!              for use in surface flux parameterizations. see convect.ps
!              file for details.
!
!     qprime: a convective downdraft specific humidity
!              perturbation scale (gm/gm).
!              for use in surface flux parameterizations. see convect.ps
!              file for details.
!
!     cbmf:   the cloud base mass flux ((kg/m**2)/s). this scalar value must
!              be stored by the calling program and returned to convect at
!              its next call. that is, the value of cbmf must be "remembered"
!              by the calling program between calls to convect.
!
!------------------------------------------------------------------------------
!
!    ***  the parameter na should in general be greater than   ***
!    ***                or equal to  nd + 1                    ***
!
!------------------------------------------------------------------------------
! modifications for regcm:
!   1. units for precipitation were change from mm/day to mm/s
!   2. the thermodynamic constants were made consistent with those
!      of regcm.
!   3. dependence of latent heat of vaporization on tempertature
!      removed because it is neglected in regcm. this is done by
!      setting cpv equal to cl and setting wlhv to the regcm value.
!   4. added cloud base (icb) and cloud top (inb) to the output to
!      compute the cloud fraction and cloud liquid water content.
!   5. each variable is now explicitly declared.  that is, the
!      "implicit none" option was added.
!   6. the value minorig is increased because mod_the thickness of the
!      lowest layer(s) is(are) too small. if the thickness is too
!      small for the given timestep, the mass of the layer is likely
!      to be evacuated.
!   7. a maximum value to the cloud base mass flux has been added.
!
      subroutine cupeman(t,q,qs,u,v,tra,p,ph,nd,na,nl,ntra,delt,iflag,  &
                       & ft,fq,fu,fv,ftra,precip,wd,tprime,qprime,cbmf, &
                       & icb,inb)
!
      implicit none
!
      real(8) :: cbmf , delt , precip , qprime , tprime , wd
      integer :: icb , iflag , inb , na , nd , nl , ntra
      real(8) , dimension(nd) :: fq , ft , fu , fv , p , ph , q , qs ,  &
                               & t , u , v
      real(8) , dimension(nd,1) :: ftra , tra
      intent (in) delt , na , ntra , ph
      intent (out) tprime , wd
      intent (inout) cbmf , fq , ft , ftra , fu , fv , icb , iflag ,    &
                   & inb , precip , q , qprime , qs , t , tra , u , v
!
      real(8) :: a2 , ad , afac , ahm , ahmax , ahmin , alt , altem ,   &
               & alv , alvnew , am , amp1 , anum , asij , asum , awat , &
               & b6 , bf2 , bsum , by , byp , c6 , cape , capem ,       &
               & cbmfold , chi , cl , coeff , cpinv ,                   &
               & cpvmcl , cwat , damps , dbo , dbosum , defrac , dei ,  &
               & delm , delp , delt0 , delti , denom , dhdp , dphinv ,  &
               & dpinv , dtma , dtmin , dtpbl , elacrit , ents , epmax ,&
               & eps , epsi , fac , fqold , frac , ftold , ftraold ,    &
               & fuold , fvold , plcl , qnew , qp1 ,                    &
               & qsm , qstm , qti , rat , rdcp , revap , rh , rm ,      &
               & rowl , scrit , sigt
      real(8) , dimension(na) :: clw , cpn , ep , evap , gz , h , hm ,  &
                               & hp , lv , lvcp , m , mp , qp , sigp ,  &
                               & th , told , tp , tratm , tv , tvp ,    &
                               & up , vp , water , wt
      real(8) , dimension(na,na) :: elij , ment , qent , sij , uent ,   &
                                  & vent
      integer :: i , ihmin , inb1 , ipbl , j , jc , jn , jtt , k , nk
      integer , dimension(na) :: nent
      real(8) :: sjmax , sjmin , smid , smin , stemp , tc , tca ,       &
               & thbar , tnew , traav , tvaplcl , tvpplcl , tvx , tvy , &
               & uav , um , vav , vm , wdtrain , x
      real(8) , dimension(na,na,ntra) :: traent
      real(8) , dimension(na,ntra) :: trap
!
! ----------------------------------------------------------------------
!
!     ***                     specify switches                        
!     ***
!     ***   ipbl: set to zero to bypass dry adiabatic adjustment      
!     *** ***    any other value results in dry adiabatic adjustment   
!     *** ***     (zero value recommended for use in models with       
!     *** ***                   boundary layer schemes)                
!     ***
!     ***   minorig: lowest level from which convection may originate 
!     *** ***     (should be first model level at which t is defined   
!     *** ***      for models using bulk pbl schemes; otherwise, it
!     should *** ***      be the first model level at which t is
!     defined above    *** ***                      the surface layer) 
!     ***
      ipbl = 0
!
!     ***        assign values of thermodynamic constants,        ***
!     ***            gravity, and liquid water density.           ***
!     ***             these should be consistent with             ***
!     ***              those used in calling program              ***
!     ***     note: these are also specified in subroutine tlift  ***
!
      cl = 2500.0D0
      rowl = 1000.0D0
!
      cpvmcl = cl - cpv
      eps = rgas/rwat
      epsi = 1.0D0/eps
      delti = 1.0D0/delt
!
!     ***  initialize output arrays and parameters  ***
!
      do i = 1 , nd
        ft(i) = 0.0D0
        fq(i) = 0.0D0
        fu(i) = 0.0D0
        fv(i) = 0.0D0
        do j = 1 , ntra
          ftra(i,j) = 0.0D0
        end do
      end do
      do i = 1 , nl + 1
        rdcp = (rgas*(1.0D0-q(i))+q(i)*rwat)/(cpd*(1.0D0-q(i))+q(i)*cpv)
        th(i) = t(i)*(1000.0D0/p(i))**rdcp
      end do
      precip = 0.0D0
      wd = 0.0D0
      tprime = 0.0D0
      qprime = 0.0D0
      iflag = 0
!
      if ( ipbl /= 0 ) then
!
!       ***            perform dry adiabatic adjustment            ***
!
        jc = 0
        do i = nl - 1 , 1 , -1
          jn = 0
          asum = th(i)*(1.0D0+q(i)*epsi-q(i))
          do j = i + 1 , nl
            asum = asum + th(j)*(1.0D0+q(j)*epsi-q(j))
            thbar = asum/dble(j+1-i)
            if ( (th(j)*(1.0D0+q(j)*epsi-q(j))) < thbar ) jn = j
          end do
          if ( i == 1 ) jn = max0(jn,2)
          if ( jn /= 0 ) then
            do
              ahm = 0.0D0
              rm = 0.0D0
              um = 0.0D0
              vm = 0.0D0
              do k = 1 , ntra
                tratm(k) = 0.0D0
              end do
              do j = i , jn
                ahm = ahm + (cpd*(1.0D0-q(j))+q(j)*cpv)*t(j)* &
                            (ph(j)-ph(j+1))
                rm = rm + q(j)*(ph(j)-ph(j+1))
                um = um + u(j)*(ph(j)-ph(j+1))
                vm = vm + v(j)*(ph(j)-ph(j+1))
                do k = 1 , ntra
                  tratm(k) = tratm(k) + tra(j,k)*(ph(j)-ph(j+1))
                end do
              end do
              dphinv = 1.0D0/(ph(i)-ph(jn+1))
              rm = rm*dphinv
              um = um*dphinv
              vm = vm*dphinv
              do k = 1 , ntra
                tratm(k) = tratm(k)*dphinv
              end do
              a2 = 0.0D0
              do j = i , jn
                q(j) = rm
                u(j) = um
                v(j) = vm
                do k = 1 , ntra
                  tra(j,k) = tratm(k)
                end do
                rdcp = (rgas*(1.0D0-q(j))+q(j)*rwat) / &
                       (cpd*(1.0D0-q(j))+q(j)*cpv)
                x = (0.001D0*p(j))**rdcp
                told(j) = t(j)
                t(j) = x
                a2 = a2 + (cpd*(1.0D0-q(j))+q(j)*cpv)*x*(ph(j)-ph(j+1))
              end do
              do j = i , jn
                th(j) = ahm/a2
                t(j) = t(j)*th(j)
                tc = told(j) - tzero
                alv = wlhv - cpvmcl*tc
                qs(j) = qs(j) + qs(j)*(1.0D0+qs(j)*(epsi-1.0D0)) * &
                        alv*(t(j)-told(j))/(rwat*told(j)*told(j))
              end do
              if ( ((th(jn+1)*(1.0D0+q(jn+1)*epsi-q(jn+1))) <  &
                    (th(jn)*(1.0D0+q(jn)*epsi-q(jn)))) ) then
                jn = jn + 1
                cycle
              end if
              if ( i == 1 ) jc = jn
              exit
            end do
          end if
        end do
!
!       ***   remove any supersaturation that results from adjustment
!       ***
        if ( jc > 1 ) then
          do j = 1 , jc
            if ( qs(j) < q(j) ) then
              alv = wlhv - cpvmcl*(t(j)-tzero)
              tnew = t(j) + alv*(q(j)-qs(j)) /             &
                      (cpd*(1.0D0-q(j))+cl*q(j)+qs(j) *    &
                      (cpv-cl+alv*alv/(rwat*t(j)*t(j))))
              alvnew = wlhv - cpvmcl*(tnew-tzero)
              qnew = (alv*q(j)-(tnew-t(j)) * &
                     (cpd*(1.0D0-q(j))+cl*q(j)))/alvnew
!rcm          precip=precip+24.*3600.*1.0e5*(ph(j)-ph(j+1))*  ! mm/d
              precip = precip + 1.0D5*(ph(j)-ph(j+1))*(q(j)-qnew)*rgti/ &
                        (delt*rowl)                         ! mm/s
              t(j) = tnew
              q(j) = qnew
              qs(j) = qnew
            end if
          end do
        end if
!
      end if
!
!     *** calculate arrays of geopotential, heat capacity and static
!     energy
      gz(1) = 0.0D0
      cpn(1) = cpd*(1.0D0-q(1)) + q(1)*cpv
      h(1) = t(1)*cpn(1)
      lv(1) = wlhv - cpvmcl*(t(1)-tzero)
      hm(1) = lv(1)*q(1)
      tv(1) = t(1)*(1.0D0+q(1)*epsi-q(1))
      ahmin = 1.0D12
      ihmin = nl
      do i = 2 , nl + 1
        tvx = t(i)*(1.0D0+q(i)*epsi-q(i))
        tvy = t(i-1)*(1.0D0+q(i-1)*epsi-q(i-1))
        gz(i) = gz(i-1) + 0.5D0*rgas*(tvx+tvy)*(p(i-1)-p(i))/ph(i)
        cpn(i) = cpd*(1.0D0-q(i)) + cpv*q(i)
        h(i) = t(i)*cpn(i) + gz(i)
        lv(i) = wlhv - cpvmcl*(t(i)-tzero)
        hm(i) = (cpd*(1.0D0-q(i))+cl*q(i))*(t(i)-t(1))+lv(i)*q(i)+gz(i)
        tv(i) = t(i)*(1.0D0+q(i)*epsi-q(i))
!
!       ***  find level of minimum moist static energy    ***
!
        if ( i >= minorig .and. hm(i) < ahmin .and. &
             hm(i) < hm(i-1) ) then
          ahmin = hm(i)
          ihmin = i
        end if
      end do
      ihmin = min0(ihmin,nl-1)
!
!     ***     find that model level below the level of minimum moist   
!     *** ***  static energy that has the maximum value of moist static
!     energy ***
      ahmax = 0.0D0
      nk = nl
      do i = minorig , ihmin
        if ( hm(i) > ahmax ) then
          nk = i
          ahmax = hm(i)
        end if
      end do
!
!     ***  check whether parcel level temperature and specific humidity
!     *** ***                          are reasonable                  
!     *** ***      skip convection if hm increases monotonically upward
!     ***
      if ( t(nk) < 250.0D0 .or. q(nk) <= 0.0D0 .or. &
           ihmin == (nl-1) ) then
        iflag = 0
        cbmf = 0.0D0
        return
      end if
!
!     ***  calculate lifted condensation level of air at parcel origin
!     level *** ***       (within 0.2% of formula of bolton, mon. wea.
!     rev.,1980)      ***
      rh = q(nk)/qs(nk)
      chi = t(nk)/(1669.0D0-122.0D0*rh-t(nk))
      plcl = p(nk)*(rh**chi)
      if ( plcl < 200.0D0 .or. plcl >= 2000.0D0 ) then
        iflag = 2
        cbmf = 0.0D0
        return
      end if
!
!     ***  calculate first level above lcl (=icb)  ***
!
      icb = nl - 1
      do i = nk + 1 , nl
        if ( p(i) < plcl ) icb = min0(icb,i)
      end do
      if ( icb >= (nl-1) ) then
        iflag = 3
        cbmf = 0.0D0
        return
      end if
!
!     *** find temperature up through icb and test for instability     
!     ***
!     *** subroutine tlift calculates part of the lifted parcel virtual
!     *** ***  temperature, the actual temperature and the adiabatic   
!     *** ***                   liquid water content                   
!     ***
      call tlift(p,t,q,qs,gz,icb,nk,tvp,tp,clw,nd,nl,1)
      do i = nk , icb
        tvp(i) = tvp(i) - tp(i)*q(nk)
      end do
!
!     ***  if there was no convection at last time step and parcel   
!     *** ***       is stable at icb then skip rest of calculation     
!     ***
      if ( dabs(cbmf) < 1.0D-30 .and. tvp(icb) <= (tv(icb)-dtmax) ) then
        iflag = 0
        return
      end if
!
!     ***  if this point is reached, moist convective adjustment is
!     necessary ***
      if ( iflag /= 4 ) iflag = 1
!
!     ***  find the rest of the lifted parcel temperatures          ***
!
      call tlift(p,t,q,qs,gz,icb,nk,tvp,tp,clw,nd,nl,2)
!
!     ***  set the precipitation efficiencies and the fraction of   ***
!     ***          precipitation falling outside of cloud           ***
!     ***      these may be functions of tp(i), p(i) and clw(i)     ***
!
      do i = 1 , nk
        ep(i) = 0.0D0
        sigp(i) = sigs
      end do
      do i = nk + 1 , nl
        tca = tp(i) - tzero
        if ( tca >= 0.0D0 ) then
          elacrit = elcrit
        else
          elacrit = elcrit*(1.0D0-tca/tlcrit)
        end if
        elacrit = dmax1(elacrit,0.0D0)
        epmax = 0.999D0
        ep(i) = epmax*(1.0D0-elacrit/dmax1(clw(i),1.0D-8))
        ep(i) = dmax1(ep(i),0.0D0)
        ep(i) = dmin1(ep(i),epmax)
        sigp(i) = sigs
      end do
!
!     ***       calculate virtual temperature and lifted parcel     ***
!     ***                    virtual temperature                    ***
!
      do i = icb + 1 , nl
        tvp(i) = tvp(i) - tp(i)*q(nk)
      end do
      tvp(nl+1) = tvp(nl) - (gz(nl+1)-gz(nl))*rcpd
!
!     ***        now initialize various arrays used in the computations
!     ***
      do i = 1 , nl + 1
        hp(i) = h(i)
        nent(i) = 0
        water(i) = 0.0D0
        evap(i) = 0.0D0
        wt(i) = omtsnow
        mp(i) = 0.0D0
        m(i) = 0.0D0
        lvcp(i) = lv(i)/cpn(i)
        do j = 1 , nl + 1
          qent(i,j) = q(j)
          elij(i,j) = 0.0D0
          ment(i,j) = 0.0D0
          sij(i,j) = 0.0D0
          uent(i,j) = u(j)
          vent(i,j) = v(j)
          do k = 1 , ntra
            traent(i,j,k) = tra(j,k)
          end do
        end do
      end do
      qp(1) = q(1)
      up(1) = u(1)
      vp(1) = v(1)
      do i = 1 , ntra
        trap(1,i) = tra(1,i)
      end do
      do i = 2 , nl + 1
        qp(i) = q(i-1)
        up(i) = u(i-1)
        vp(i) = v(i-1)
        do j = 1 , ntra
          trap(i,j) = tra(i-1,j)
        end do
      end do
!
!     ***  find the first model level (inb1) above the parcel's      ***
!     ***          highest level of neutral buoyancy                 ***
!     ***     and the highest level of positive cape (inb)           ***
!
      cape = 0.0D0
      capem = 0.0D0
      inb = icb + 1
      inb1 = inb
      byp = 0.0D0
      do i = icb + 1 , nl - 1
        by = (tvp(i)-tv(i))*(ph(i)-ph(i+1))/p(i)
        cape = cape + by
        if ( by >= 0.0D0 ) inb1 = i + 1
        if ( cape > 0.0D0 ) then
          inb = i + 1
          byp = (tvp(i+1)-tv(i+1))*(ph(i+1)-ph(i+2))/p(i+1)
          capem = cape
        end if
      end do
      inb = max0(inb,inb1)
      cape = capem + byp
      defrac = capem - cape
      defrac = dmax1(defrac,0.001D0)
      frac = -cape/defrac
      frac = dmin1(frac,1.0D0)
      frac = dmax1(frac,0.0D0)
!
!     ***   calculate liquid water static energy of lifted parcel   ***
!
      do i = icb , inb
        hp(i) = h(nk) + (lv(i)+(cpd-cpv)*t(i))*ep(i)*clw(i)
      end do
!
!     ***  calculate cloud base mass flux and rates of mixing, m(i), 
!     *** ***                   at each model level                    
!     ***
      dbosum = 0.0D0
!
!     ***     interpolate difference between lifted parcel and      ***
!     ***  environmental temperatures to lifted condensation level  ***
!
      tvpplcl = tvp(icb-1) - rgas*tvp(icb-1)*(p(icb-1)-plcl)            &
              & /(cpn(icb-1)*p(icb-1))
      tvaplcl = tv(icb) + (tvp(icb)-tvp(icb+1))*(plcl-p(icb))           &
              & /(p(icb)-p(icb+1))
      dtpbl = 0.0D0
      do i = nk , icb - 1
        dtpbl = dtpbl + (tvp(i)-tv(i))*(ph(i)-ph(i+1))
      end do
      dtpbl = dtpbl/(ph(nk)-ph(icb))
      dtmin = tvpplcl - tvaplcl + dtmax + dtpbl
      dtma = dtmin
!
!     ***  adjust cloud base mass flux   ***
!
      cbmfold = cbmf
      delt0 = 300.0D0
      damps = damp*delt/delt0
      cbmf = (1.0D0-damps)*cbmf + 0.1D0*alphae*dtma
      cbmf = dmax1(cbmf,0.0D0)
!
!     *** if cloud base mass flux is zero, skip rest of calculation  ***
!
      if ( dabs(cbmf) < 1.0D-30 .and. dabs(cbmfold) < 1.0D-30 ) return
!
!     ***   calculate rates of mixing,  m(i)   ***
!
      m(icb) = 0.0D0
      do i = icb + 1 , inb
        k = min0(i,inb1)
        dbo = dabs(tv(k)-tvp(k)) + entp*0.02D0*(ph(k)-ph(k+1))
        dbosum = dbosum + dbo
        m(i) = cbmf*dbo
      end do
      do i = icb + 1 , inb
        m(i) = m(i)/dbosum
      end do
!
!     ***  calculate entrained air mass flux (ment), total water mixing
!     *** ***     ratio (qent), total condensed water (elij), and
!     mixing     *** ***                        fraction (sij)         
!     ***
      do i = icb + 1 , inb
        qti = q(nk) - ep(i)*clw(i)
        do j = icb , inb
          bf2 = 1.0D0 + lv(j)*lv(j)*qs(j)/(rwat*t(j)*t(j)*cpd)
          anum = h(j) - hp(i) + (cpv-cpd)*t(j)*(qti-q(j))
          denom = h(i) - hp(i) + (cpd-cpv)*(q(i)-qti)*t(j)
          dei = denom
          if ( dabs(dei) < 0.01D0 ) dei = 0.01D0
          sij(i,j) = anum/dei
          sij(i,i) = 1.0D0
          altem = sij(i,j)*q(i) + (1.0D0-sij(i,j))*qti - qs(j)
          altem = altem/bf2
          cwat = clw(j)*(1.0D0-ep(j))
          stemp = sij(i,j)
          if ( (stemp < 0.0D0 .or. stemp > 1.0D0 .or. &
                altem > cwat) .and. j > i ) then
            anum = anum - lv(j)*(qti-qs(j)-cwat*bf2)
            denom = denom + lv(j)*(q(i)-qti)
            if ( dabs(denom) < 0.01D0 ) denom = 0.01D0
            sij(i,j) = anum/denom
            altem = sij(i,j)*q(i) + (1.0D0-sij(i,j))*qti - qs(j)
            altem = altem - (bf2-1.0D0)*cwat
          end if
          if ( sij(i,j) > 0.0D0 .and. sij(i,j) < 0.9D0 ) then
            qent(i,j) = sij(i,j)*q(i) + (1.0D0-sij(i,j))*qti
            uent(i,j) = sij(i,j)*u(i) + (1.0D0-sij(i,j))*u(nk)
            vent(i,j) = sij(i,j)*v(i) + (1.0D0-sij(i,j))*v(nk)
            do k = 1 , ntra
              traent(i,j,k) = sij(i,j)*tra(i,k) + &
                              (1.0D0-sij(i,j))*tra(nk,k)
            end do
            elij(i,j) = altem
            elij(i,j) = dmax1(0.0D0,elij(i,j))
            ment(i,j) = m(i)/(1.0D0-sij(i,j))
            nent(i) = nent(i) + 1
          end if
          sij(i,j) = dmax1(0.0D0,sij(i,j))
          sij(i,j) = dmin1(1.0D0,sij(i,j))
        end do
!
!       ***   if no air can entrain at level i assume that updraft
!       detrains  *** ***   at that level and calculate detrained air
!       flux and properties  ***
        if ( nent(i) == 0 ) then
          ment(i,i) = m(i)
          qent(i,i) = q(nk) - ep(i)*clw(i)
          uent(i,i) = u(nk)
          vent(i,i) = v(nk)
          do j = 1 , ntra
            traent(i,i,j) = tra(nk,j)
          end do
          elij(i,i) = clw(i)
          sij(i,i) = 1.0D0
        end if
      end do
      sij(inb,inb) = 1.0D0
!
!     ***  normalize entrained air mass fluxes to represent equal  ***
!     ***              probabilities of mixing                     ***
!
      do i = icb + 1 , inb
        if ( nent(i) /= 0 ) then
          qp1 = q(nk) - ep(i)*clw(i)
          anum = h(i) - hp(i) - lv(i)*(qp1-qs(i))
          denom = h(i) - hp(i) + lv(i)*(q(i)-qp1)
          if ( dabs(denom) < 0.01D0 ) denom = 0.01D0
          scrit = anum/denom
          alt = qp1 - qs(i) + scrit*(q(i)-qp1)
          if ( alt < 0.0D0 ) scrit = 1.0D0
          scrit = dmax1(scrit,0.0D0)
          asij = 0.0D0
          smin = 1.0D0
          do j = icb , inb
            if ( sij(i,j) > 0.0D0 .and. sij(i,j) < 0.9D0 ) then
              if ( j > i ) then
                smid = dmin1(sij(i,j),scrit)
                sjmax = smid
                sjmin = smid
                if ( smid < smin .and. sij(i,j+1) < smid ) then
                  smin = smid
                  sjmax = dmin1(sij(i,j+1),sij(i,j),scrit)
                  sjmin = dmax1(sij(i,j-1),sij(i,j))
                  sjmin = dmin1(sjmin,scrit)
                end if
              else
                sjmax = dmax1(sij(i,j+1),scrit)
                smid = dmax1(sij(i,j),scrit)
                sjmin = 0.0D0
                if ( j > 1 ) sjmin = sij(i,j-1)
                sjmin = dmax1(sjmin,scrit)
              end if
              delp = dabs(sjmax-smid)
              delm = dabs(sjmin-smid)
              asij = asij + (delp+delm)*(ph(j)-ph(j+1))
              ment(i,j) = ment(i,j)*(delp+delm)*(ph(j)-ph(j+1))
            end if
          end do
          asij = dmax1(1.0D-21,asij)
          asij = 1.0D0/asij
          do j = icb , inb
            ment(i,j) = ment(i,j)*asij
          end do
          bsum = 0.0D0
          do j = icb , inb
            bsum = bsum + ment(i,j)
          end do
          if ( bsum < 1.0D-18 ) then
            nent(i) = 0
            ment(i,i) = m(i)
            qent(i,i) = q(nk) - ep(i)*clw(i)
            uent(i,i) = u(nk)
            vent(i,i) = v(nk)
            do j = 1 , ntra
              traent(i,i,j) = tra(nk,j)
            end do
            elij(i,i) = clw(i)
            sij(i,i) = 1.0D0
          end if
        end if
      end do
!
!     ***  check whether ep(inb)=0, if so, skip precipitating    ***
!     ***             downdraft calculation                      ***
!
      if ( ep(inb) >= 0.0001D0 ) then
!
!       ***  integrate liquid water equation to find condensed water  
!       *** ***                and condensed water flux                
!       ***
        jtt = 2
!
!       ***                    begin downdraft loop                   
!       ***
        do i = inb , 1 , -1
!
!         ***              calculate detrained precipitation           
!         ***
          wdtrain = gti*ep(i)*m(i)*clw(i)
          if ( i > 1 ) then
            do j = 1 , i - 1
              awat = elij(j,i) - (1.0D0-ep(i))*clw(i)
              awat = dmax1(0.0D0,awat)
              wdtrain = wdtrain + gti*awat*ment(j,i)
            end do
          end if
!
!         ***    find rain water and evaporation using provisional   ***
!         ***              estimates of qp(i)and qp(i-1)             ***
!
!
!         ***  value of terminal velocity and coefficient of
!         evaporation for snow   ***
          coeff = coeffs
          wt(i) = omtsnow
!
!         ***  value of terminal velocity and coefficient of
!         evaporation for rain   ***
          if ( t(i) > tzero) then
            coeff = coeffr
            wt(i) = omtrain
          end if
          qsm = 0.50D0*(q(i)+qp(i+1))
          afac = coeff*ph(i)*(qs(i)-qsm)/(1.0D4+2.0D3*ph(i)*qs(i))
          afac = dmax1(afac,0.0D0)
          sigt = sigp(i)
          sigt = dmax1(0.0D0,sigt)
          sigt = dmin1(1.0D0,sigt)
          b6 = 100.0D0*(ph(i)-ph(i+1))*sigt*afac/wt(i)
          c6 = (water(i+1)*wt(i+1)+wdtrain/sigd)/wt(i)
          revap = 0.5D0*(-b6+dsqrt(b6*b6+4.0D0*c6))
          evap(i) = sigt*afac*revap
          water(i) = revap*revap
!
!         ***  calculate precipitating downdraft mass flux under     ***
!         ***              hydrostatic approximation                 ***
!
          if ( i /= 1 ) then
            dhdp = (h(i)-h(i-1))/(p(i-1)-p(i))
            dhdp = dmax1(dhdp,10.0D0)
            mp(i) = 100.0D0*rgti*lv(i)*sigd*evap(i)/dhdp
            mp(i) = dmax1(mp(i),0.0D0)
!
!           ***   add small amount of inertia to downdraft             
!           ***
            fac = 20.0D0/(ph(i-1)-ph(i))
            mp(i) = (fac*mp(i+1)+mp(i))/(1.0D0+fac)
!
!           ***      force mp to decrease linearly to zero             
!           *** ***      between about 950 mb and the surface          
!           ***
            if ( p(i) > (0.949D0*p(1)) ) then
              jtt = max0(jtt,i)
              mp(i) = mp(jtt)*(p(1)-p(i))/(p(1)-p(jtt))
            end if
          end if
!
!         ***       find mixing ratio of precipitating downdraft     ***
!
          if ( i /= inb ) then
            if ( i == 1 ) then
              qstm = qs(1)
            else
              qstm = qs(i-1)
            end if
            if ( mp(i) > mp(i+1) ) then
              rat = mp(i+1)/mp(i)
              qp(i) = qp(i+1)*rat + q(i)*(1.0D0-rat)                    &
                    & + 100.D0*rgti*sigd*(ph(i)-ph(i+1))*(evap(i)/mp(i))
              up(i) = up(i+1)*rat + u(i)*(1.0D0-rat)
              vp(i) = vp(i+1)*rat + v(i)*(1.0D0-rat)
              do j = 1 , ntra
                trap(i,j) = trap(i+1,j)*rat + trap(i,j)*(1.0D0-rat)
              end do
            else if ( mp(i+1) > 0.0D0 ) then
              qp(i) = (gz(i+1)-gz(i)+qp(i+1)*(lv(i+1)+t(i+1)*(cl-cpd))  &
                    & +cpd*(t(i+1)-t(i)))/(lv(i)+t(i)*(cl-cpd))
              up(i) = up(i+1)
              vp(i) = vp(i+1)
              do j = 1 , ntra
                trap(i,j) = trap(i+1,j)
              end do
            end if
            qp(i) = dmin1(qp(i),qstm)
            qp(i) = dmax1(qp(i),0.0D0)
          end if
        end do
!
!       ***  calculate surface precipitation in mm/s     ***
!
!rcm    precip=precip+wt(1)*sigd*water(1)*3600.*24000./(rowl*g)  ! mm/d
        precip = precip + wt(1)*sigd*water(1)*1000.0D0/(rowl*gti)
                                                        ! mm/s
      end if
!
!
!     ***  calculate downdraft velocity scale and surface temperature
!     and  *** ***                    water vapor fluctuations         
!     ***
      wd = betae*dabs(mp(icb))*0.01D0*rgas*t(icb)/(sigd*p(icb))
      qprime = 0.5D0*(qp(1)-q(1))
      tprime = wlhv*qprime*rcpd
!
!     ***  calculate tendencies of lowest level potential temperature 
!     *** ***                      and mixing ratio                    
!     ***
      dpinv = 0.01D0/(ph(1)-ph(2))
      am = 0.0D0
      if ( nk == 1 ) then
        do k = 2 , inb
          am = am + m(k)
        end do
      end if
      if ( (2.0D0*gti*dpinv*am) >= delti ) iflag = 4
      ft(1) = ft(1) + gti*dpinv*am*(t(2)-t(1)+(gz(2)-gz(1))/cpn(1))
      ft(1) = ft(1) - lvcp(1)*sigd*evap(1)
      ft(1) = ft(1) + sigd*wt(2)*(cl-cpd)*water(2)*(t(2)-t(1))          &
            & *dpinv/cpn(1)
      fq(1) = fq(1) + gti*mp(2)*(qp(2)-q(1))*dpinv + sigd*evap(1)
      fq(1) = fq(1) + gti*am*(q(2)-q(1))*dpinv
      fu(1) = fu(1) + gti*dpinv*(mp(2)*(up(2)-u(1))+am*(u(2)-u(1)))
      fv(1) = fv(1) + gti*dpinv*(mp(2)*(vp(2)-v(1))+am*(v(2)-v(1)))
      do j = 1 , ntra
        ftra(1,j) = ftra(1,j) + gti*dpinv*(mp(2)*(trap(2,j)-tra(1,j)) + &
                    am*(tra(2,j)-tra(1,j)))
      end do
      do j = 2 , inb
        fq(1) = fq(1) + gti*dpinv*ment(j,1)*(qent(j,1)-q(1))
        fu(1) = fu(1) + gti*dpinv*ment(j,1)*(uent(j,1)-u(1))
        fv(1) = fv(1) + gti*dpinv*ment(j,1)*(vent(j,1)-v(1))
        do k = 1 , ntra
          ftra(1,k) = ftra(1,k) + gti*dpinv*ment(j,1)                   &
                    & *(traent(j,1,k)-tra(1,k))
        end do
      end do
!
!     ***  calculate tendencies of potential temperature and mixing
!     ratio  *** ***               at levels above the lowest level    
!     ***
!     ***  first find the net saturated updraft and downdraft mass
!     fluxes  *** ***                      through each level          
!     ***
      do i = 2 , inb
        dpinv = 0.01D0/(ph(i)-ph(i+1))
        cpinv = 1.0D0/cpn(i)
        amp1 = 0.0D0
        ad = 0.0D0
        if ( i >= nk ) then
          do k = i + 1 , inb + 1
            amp1 = amp1 + m(k)
          end do
        end if
        do k = 1 , i
          do j = i + 1 , inb + 1
            amp1 = amp1 + ment(k,j)
          end do
        end do
        if ( (2.0D0*gti*dpinv*amp1) >= delti ) iflag = 4
        do k = 1 , i - 1
          do j = i , inb
            ad = ad + ment(j,k)
          end do
        end do
        ft(i) = ft(i) + gti*dpinv*(amp1*(t(i+1)-t(i)+ &
                   (gz(i+1)-gz(i))*cpinv)-ad*(t(i)-t(i-1)+ &
                   (gz(i)-gz(i-1))*cpinv)) - sigd*lvcp(i)*evap(i)
        ft(i) = ft(i) + gti*dpinv*ment(i,i)                             &
              & *(hp(i)-h(i)+t(i)*(cpv-cpd)*(q(i)-qent(i,i)))*cpinv
        ft(i) = ft(i) + sigd*wt(i+1)*(cl-cpd)*water(i+1)*(t(i+1)-t(i))  &
              & *dpinv*cpinv
        fq(i) = fq(i) + gti*dpinv*(amp1*(q(i+1)-q(i))-ad*(q(i)-q(i-1)))
        fu(i) = fu(i) + gti*dpinv*(amp1*(u(i+1)-u(i))-ad*(u(i)-u(i-1)))
        fv(i) = fv(i) + gti*dpinv*(amp1*(v(i+1)-v(i))-ad*(v(i)-v(i-1)))
        do k = 1 , ntra
          ftra(i,k) = ftra(i,k)                                         &
                    & + gti*dpinv*(amp1*(tra(i+1,k)-tra(i,k))-ad*       &
                    & (tra(i,k)-tra(i-1,k)))
        end do
        do k = 1 , i - 1
          awat = elij(k,i) - (1.0D0-ep(i))*clw(i)
          awat = dmax1(awat,0.0D0)
          fq(i) = fq(i) + gti*dpinv*ment(k,i)*(qent(k,i)-awat-q(i))
          fu(i) = fu(i) + gti*dpinv*ment(k,i)*(uent(k,i)-u(i))
          fv(i) = fv(i) + gti*dpinv*ment(k,i)*(vent(k,i)-v(i))
          do j = 1 , ntra
            ftra(i,j) = ftra(i,j) + gti*dpinv*ment(k,i)                 &
                      & *(traent(k,i,j)-tra(i,j))
          end do
        end do
        do k = i , inb
          fq(i) = fq(i) + gti*dpinv*ment(k,i)*(qent(k,i)-q(i))
          fu(i) = fu(i) + gti*dpinv*ment(k,i)*(uent(k,i)-u(i))
          fv(i) = fv(i) + gti*dpinv*ment(k,i)*(vent(k,i)-v(i))
          do j = 1 , ntra
            ftra(i,j) = ftra(i,j) + gti*dpinv*ment(k,i)                 &
                      & *(traent(k,i,j)-tra(i,j))
          end do
        end do
        fq(i) = fq(i) + sigd*evap(i)+ gti*(mp(i+1)                      &
              & *(qp(i+1)-q(i))-mp(i)*(qp(i)-q(i-1)))*dpinv
        fu(i) = fu(i) + gti*(mp(i+1)*(up(i+1)-u(i))-mp(i)*              &
              & (up(i)-u(i-1)))*dpinv
        fv(i) = fv(i) + gti*(mp(i+1)*(vp(i+1)-v(i))-mp(i)*              &
              & (vp(i)-v(i-1)))*dpinv
        do j = 1 , ntra
          ftra(i,j) = ftra(i,j)                                         &
                    & + gti*dpinv*(mp(i+1)*(trap(i+1,j)-tra(i,j))-mp(i) &
                    & *(trap(i,j)-trap(i-1,j)))
        end do
      end do
!
!     *** adjust tendencies at top of convection layer to reflect  ***
!     ***       actual position of the level zero cape             ***
!
      fqold = fq(inb)
      fq(inb) = fq(inb)*(1.0D0-frac)
      fq(inb-1) = fq(inb-1)                                             &
                & + frac*fqold*((ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb)))&
                & *lv(inb)/lv(inb-1)
      ftold = ft(inb)
      ft(inb) = ft(inb)*(1.0D0-frac)
      ft(inb-1) = ft(inb-1)                                             &
                & + frac*ftold*((ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb)))&
                & *cpn(inb)/cpn(inb-1)
      fuold = fu(inb)
      fu(inb) = fu(inb)*(1.0D0-frac)
      fu(inb-1) = fu(inb-1)                                             &
                & + frac*fuold*((ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb)))
      fvold = fv(inb)
      fv(inb) = fv(inb)*(1.0D0-frac)
      fv(inb-1) = fv(inb-1)                                             &
                & + frac*fvold*((ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb)))
      do k = 1 , ntra
        ftraold = ftra(inb,k)
        ftra(inb,k) = ftra(inb,k)*(1.0D0-frac)
        ftra(inb-1,k) = ftra(inb-1,k) + frac*ftraold*(ph(inb)-ph(inb+1))&
                      & /(ph(inb-1)-ph(inb))
      end do
!
!     ***   very slightly adjust tendencies to force exact   ***
!     ***     enthalpy, momentum and tracer conservation     ***
!
      ents = 0.0D0
      uav = 0.0D0
      vav = 0.0D0
      do i = 1 , inb
        ents = ents + (cpn(i)*ft(i)+lv(i)*fq(i))*(ph(i)-ph(i+1))
        uav = uav + fu(i)*(ph(i)-ph(i+1))
        vav = vav + fv(i)*(ph(i)-ph(i+1))
      end do
      ents = ents/(ph(1)-ph(inb+1))
      uav = uav/(ph(1)-ph(inb+1))
      vav = vav/(ph(1)-ph(inb+1))
      do i = 1 , inb
        ft(i) = ft(i) - ents/cpn(i)
        fu(i) = (1.0D0-cu)*(fu(i)-uav)
        fv(i) = (1.0D0-cu)*(fv(i)-vav)
      end do
      do k = 1 , ntra
        traav = 0.0D0
        do i = 1 , inb
          traav = traav + ftra(i,k)*(ph(i)-ph(i+1))
        end do
        traav = traav/(ph(1)-ph(inb+1))
        do i = 1 , inb
          ftra(i,k) = ftra(i,k) - traav
        end do
      end do
!
!     ***           return           ***
!
      end subroutine cupeman
!
! Calculate lifting level temperature
!
      subroutine tlift(p,t,q,qs,gz,icb,nk,tvp,tpk,clw,nd,nl,kk)

      implicit none
!
      integer :: icb , kk , nd , nk , nl
      real(8) , dimension(nd) :: clw , gz , p , q , qs , t , tpk , tvp
      intent (in) gz , icb , kk , nd , nk , nl , p , q , qs , t
      intent (out) tvp
      intent (inout) clw , tpk
!
      real(8) :: ah0 , ahg , alv , cl , cpinv , cpp ,                   &
               & cpvmcl , denom , eps , epsi , es , qg , rg ,           &
               & s , tc , tg
      integer :: i , j , nsb , nst
!
!     ***   assign values of thermodynamic constants     ***
!
      cl = 2500.0D0
!
      cpvmcl = cl - cpv
      eps = rgas/rwat
      epsi = 1.0D0/eps
!
!     ***  calculate certain parcel quantities, including static energy
!     ***
      ah0 = (cpd*(1.0D0-q(nk))+cl*q(nk))*t(nk) + q(nk)   &
          & *(wlhv-cpvmcl*(t(nk)-tzero)) + gz(nk)
      cpp = cpd*(1.0D0-q(nk)) + q(nk)*cpv
      cpinv = 1.0D0/cpp
!
      if ( kk == 1 ) then
!
!       ***   calculate lifted parcel quantities below cloud base   ***
!
        do i = 1 , icb - 1
          clw(i) = 0.0D0
        end do
        do i = nk , icb - 1
          tpk(i) = t(nk) - (gz(i)-gz(nk))*cpinv
          tvp(i) = tpk(i)*(1.0D0+q(nk)*epsi)
        end do
      end if
!
!     ***  find lifted parcel quantities above cloud base    ***
!
      nst = icb
      nsb = icb
      if ( kk == 2 ) then
        nst = nl
        nsb = icb + 1
      end if
      do i = nsb , nst
        tg = t(i)
        qg = qs(i)
        alv = wlhv - cpvmcl*(t(i)-tzero)
        do j = 1 , 2
          s = cpd + alv*alv*qg/(rwat*t(i)*t(i))
          s = 1.0D0/s
          ahg = cpd*tg + (cl-cpd)*q(nk)*t(i) + alv*qg + gz(i)
          tg = tg + s*(ah0-ahg)
          tg = dmax1(tg,35.0D0)
          tc = tg - tzero
          denom = 243.5D0 + tc
          if ( tc >= 0.0D0 ) then
            es = 6.112D0*dexp(17.67D0*tc/denom)
          else
            es = dexp(23.33086D0-6111.72784D0/tg+0.15215D0*dlog(tg))
          end if
          qg = eps*es/(p(i)-es*(1.0D0-eps))
        end do
        alv = wlhv - cpvmcl*(t(i)-tzero)
        tpk(i) = (ah0-(cl-cpd)*q(nk)*t(i)-gz(i)-alv*qg)*rcpd
        clw(i) = q(nk) - qg
        clw(i) = dmax1(0.0D0,clw(i))
        rg = qg/(1.0D0-q(nk))
        tvp(i) = tpk(i)*(1.0D0+rg*epsi)
      end do
!
      end subroutine tlift
!
      end module mod_cu_em
