
!c -------------------------------------------------------------------
      SUBROUTINE  gas_phase(j)

!abt added below for restart purposes
      use mod_chem_interface   , only : nvar, chemten
      use mod_dynparam
      use mod_constants
      use mod_runparams, only : dtchem, dsigma, a,r8pt
      use mod_indices
      use mod_mainchem , only : chia, chib
      use mod_bats     , only : coszrs
      use mod_date
      use mod_main     , only : atm1, sps1, mddom
      use mod_cvaria   , only : chic
      use mod_sun      , only : zenitm
!abt above
      IMPLICIT NONE

      include 'molec_wght.param'


      real*8 taa(2:iym2,1:kz),psaa(2:iym2,1:kz),rhoa(2:iym2,1:kz)
      real*8 chemin(2:iym2,1:kz,1:jxp,nvar),chemox(2:iym2,1:kz,1:jxp,nvar)
      real*8 jphoto(2:iym2,1:kz,56)
      real(8) , dimension(iy,jxp) :: psdot
      real*8 avo,airmw,cfactor,ccfactor,pfact,ro3,rh2o2,kb
      real*8 zena(2:iym2),zenithan(250)
      real*8 vdep(iy-1,ntr)
      real*8 fact1,fact2,srctemp
      real*8 mwa
      integer idatein,tstart
      integer itr,i,j,k,nn
      real*8 TOD   !abt added time of day
      INTEGER, PARAMETER :: jvO2=1,jvO3a=2,jvO3b=3,jvNO2=4,jvNO3a=5     &
            ,jvNO3b=6,jvN2O5a=7,jvN2O5b=8,jvN2O=9,jvHO2=10             &
            ,jvH2O2=11,jvHNO2=12,jvHNO3=13,jvHNO4=14,jvCH2Oa=15        &
            ,jvCH2Ob=16,jvCH3CHOa=17,jvCH3CHOb=18,jvCH3CHOc=19         &
            ,jvC2H5CHO=20,jvCHOCHO=21,jvCH3COCHO=22                    &
            ,jvCH3COCH3=23,jvCH3OOH=24,jvCH3ONO2=25,jvPAN=26


            chemin(:,:,j,:) = 0.0
            chemox(:,:,j,:) = 0.0

          if(j .eq. 1) then 
!          write(*,*)'gasphas',maxval(chib(2:iym2,kz,1:jxp,idms)),idms
          end if
!cah      write(*,*)ntime,j
      avo    = 6.022E23
      mwa    = 28.97
      kb     = 1.380658E-19
      TOD    = xtime/60.0 + gmt/1.0
!c      write(*,*)'start',1e9*chib(20,18,37,io3)

   
      do k  =1,kz
      do i  = 2,iym2
        taa(i,k) = atm1%t(i,k,j)/sps1%ps(i,j) 
        psaa(i,k)= (sps1%ps(i,j)*a(k)+ptop)
!        psaa(i,k) = (sps1%ps(i,j)*a(k)+r8pt)
      end do
      end do 

      DO k  =1,kz
      DO i  = 2,iym2
      cfactor              = psaa(i,k)*10./(kb*taa(i,k))
!      cfactor              = cfactor/psdot(i,j)!sps1%ps(i,j)
      cfactor              = cfactor/sps1%ps(i,j)
      call zenitm(coszrs,iy,j)
      zena(i) = acosd(coszrs(i))
!c non transported species
      chemin(i,k,j,ind_H2O)  = atm1%qv(i,k,j)*cfactor

      
      chemin(i,k,j,ind_O3)   = chib(i,k,j,io3)*cfactor*mwa/W_O3
      chemin(i,k,j,ind_NO2)  = chib(i,k,j,ino2)*cfactor*mwa/W_NO2
      chemin(i,k,j,ind_NO)   = chib(i,k,j,ino)*cfactor*mwa/W_NO
      chemin(i,k,j,ind_CO)   = chib(i,k,j,ico)*cfactor*mwa/W_CO
      chemin(i,k,j,ind_H2O2) = chib(i,k,j,ih2o2)*cfactor*mwa/W_H2O2
      chemin(i,k,j,ind_HNO3) = chib(i,k,j,ihno3)*cfactor*mwa/W_HNO3
      chemin(i,k,j,ind_N2O5) = chib(i,k,j,in2o5)*cfactor*mwa/W_N2O5
      chemin(i,k,j,ind_SO2)  = chib(i,k,j,iso2)*cfactor*mwa/W_SO2
      chemin(i,k,j,ind_SULF) = chib(i,k,j,iso4)*cfactor*mwa/W_SULF
      chemin(i,k,j,ind_DMS)  = chib(i,k,j,idms)*cfactor*mwa/W_DMS
      chemin(i,k,j,ind_HCHO) = chib(i,k,j,ihcho)*cfactor*mwa/W_HCHO
      chemin(i,k,j,ind_ALD2) = chib(i,k,j,iald2)*cfactor*mwa/W_ALD2
      chemin(i,k,j,ind_ISOP) = chib(i,k,j,iisop)*cfactor*mwa/W_ISOP
      chemin(i,k,j,ind_C2H6) = chib(i,k,j,ic2h6)*cfactor*mwa/W_C2H6
      chemin(i,k,j,ind_PAR)  = chib(i,k,j,ipar)*cfactor*mwa/W_C3H8
      chemin(i,k,j,ind_ETHE) = chib(i,k,j,iethe)*cfactor*mwa/W_ETHENE
      chemin(i,k,j,ind_PRPE) = chib(i,k,j,iolt)*cfactor*mwa/W_OLT
      chemin(i,k,j,ind_BUTE) = chib(i,k,j,ioli)*cfactor*mwa/W_OLI
      chemin(i,k,j,ind_TOLU) = chib(i,k,j,itolue)*cfactor*mwa/W_TOLU
      chemin(i,k,j,ind_XYLE) = chib(i,k,j,ixyl)*cfactor*mwa/W_XYLE
      chemin(i,k,j,ind_PAN)  = chib(i,k,j,ipan)*cfactor*mwa/W_PAN
      chemin(i,k,j,ind_CH4)  = chib(i,k,j,ich4)*cfactor*mwa/W_CH4
      chemin(i,k,j,ind_MOH)  = chib(i,k,j,imoh)*cfactor*mwa/W_MOH
      chemin(i,k,j,ind_ACET) = chib(i,k,j,iacet)*cfactor*mwa/W_ACET

      END DO
      END DO
      idatein = (lyear-1900)*10000+lmonth*100+lday

!      write(*,*)minval(psaa(2:ilxm,:)),minval(taa(2:ilxm,:)),j
      if (j .eq. 1)then
!      write(*,*)'Before',maxval(chemin(2:iym2,kz,ind_DMS))
!      write(*,*)'Before',maxval(chib(2:iym2,kz,j,io3))
      end if
!      CALL chemistry(j,chemin,chemox,taa,psaa,zena    &
!                    ,ktau,idatein,TOD)
      CALL chemistry(j,chemin(:,:,j,:),chemox(:,:,j,:),      &
                       taa,psaa,zena,ktau,idatein,TOD)
!      if (j .eq. 5) then
!      write(*,*)'after',maxval(chemox(2:iym2,kz,ind_O3))
!      end if

!cc      CALL chemistry(j,chemin,chemox,taa,psaa,zena
!cc     *               ,ktau,idatein)

      DO k=1,kz !,kx
      DO i=2,iym2

      cfactor              = psaa(i,k)*10./(kb*taa(i,k))
      pfact                = sps1%ps(i,j)/(cfactor*28.97)
!c      pfact                = psb(i,j)/(cfactor)

!      !  pfact to convert chemox to kg/kg

      chia(i,k,j,io3)    = chemox(i,k,j,ind_O3)*pfact*W_O3
      chia(i,k,j,ino2)   = chemox(i,k,j,ind_NO2)*pfact*W_NO2
      chia(i,k,j,ino)    = chemox(i,k,j,ind_NO)*pfact*W_NO
      chia(i,k,j,ico)    = chemox(i,k,j,ind_CO)*pfact*W_CO
      chia(i,k,j,ih2o2)  = chemox(i,k,j,ind_H2O2)*pfact*W_H2O2
      chia(i,k,j,ihno3)  = chemox(i,k,j,ind_HNO3)*pfact*W_HNO3
      chia(i,k,j,in2o5)  = chemox(i,k,j,ind_N2O5)*pfact*W_N2O5
      chia(i,k,j,iso2)   = chemox(i,k,j,ind_SO2)*pfact*W_SO2
      chia(i,k,j,iso4)   = chemox(i,k,j,ind_SULF)*pfact*W_SULF
      chia(i,k,j,idms)   = chemox(i,k,j,ind_DMS)*pfact*W_DMS
      chia(i,k,j,ihcho)  = chemox(i,k,j,ind_HCHO)*pfact*W_HCHO
      chia(i,k,j,iald2)  = chemox(i,k,j,ind_ALD2)*pfact*W_ALD2
      chia(i,k,j,iisop)  = chemox(i,k,j,ind_ISOP)*pfact*W_ISOP
      chia(i,k,j,ic2h6)  = chemox(i,k,j,ind_C2H6)*pfact*W_C2H6
      chia(i,k,j,ipar)   = chemox(i,k,j,ind_PAR)*pfact*W_C3H8
      chia(i,k,j,itolue) = chemox(i,k,j,ind_TOLU)*pfact*W_TOLU
      chia(i,k,j,ixyl)   = chemox(i,k,j,ind_XYLE)*pfact*W_XYLE
      chia(i,k,j,iethe)  = chemox(i,k,j,ind_ETHE)*pfact*W_ETHENE
      chia(i,k,j,ipan)   = chemox(i,k,j,ind_PAN)*pfact*W_PAN
      chia(i,k,j,ich4)   = chemox(i,k,j,ind_CH4)*pfact*W_CH4
      chia(i,k,j,iolt)   = chemox(i,k,j,ind_PRPE)*pfact*W_OLT
      chia(i,k,j,ioli)   = chemox(i,k,j,ind_BUTE)*pfact*W_OLI
      chia(i,k,j,imoh)   = chemox(i,k,j,ind_MOH)*pfact*W_MOH
      chia(i,k,j,iacet)  = chemox(i,k,j,ind_ACET)*pfact*W_ACET

      END DO  ! kx
      END DO  ! ilxm


554     format(I6,5x,5(1pe10.3,2x))      
!cccccccccccccccccccccccccccccccccccccccccccccccc       

!c      CLOSE(88)
        return
        END
