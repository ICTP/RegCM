      subroutine aeroppt(rh,pint,
     &        tauxar_mix, tauasc_mix, gtota_mix, ftota_mix,
     &     tauxar_mix_cs, tauasc_mix_cs, gtota_mix_cs, ftota_mix_cs )

      implicit none

#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/parrad.cb'
      include 'Commons/param2.cb'
      include 'Commons/aeroppt.cb'
#ifdef MPP1
      include 'Common2/trachem.cb'
#else
      include 'Commons/trachem.cb'
#endif
      include 'Commons/crdcon.cb'

      real*8  pint(plond,plevp)   ! Interface pressure
      real*8  rh(plond,plevr),prop(4),frac(4)
c  output arguments
c  radiative properties required by radded

      real*8  tauxar_mix(plond,0:plev,nspi)  ! aerosol extinction optical depth
      real*8  tauasc_mix(plond,0:plev,nspi)
      real*8  gtota_mix(plond,0:plev,nspi)
      real*8  ftota_mix(plond,0:plev,nspi)

      real*8  tauxar_mix_cs(plond,nspi)  ! aerosol extinction optical depth
      real*8  tauasc_mix_cs(plond,nspi)
      real*8  gtota_mix_cs(plond,nspi)
      real*8  ftota_mix_cs(plond,nspi)

      real*8  aermmr        ! Radiation level aerosol mass mixing ratio
      real*8  aermmb        ! background aerosol
      common/aermr2/ aermmr(plond,plevr,ntr),aermmb(plond,plevr)

C
C Aerosol radiative properties (local arrays)
C
      real*8  uaer(plond,0:plev,ntr)
      real*8  tauxar(plond,0:plev,ntr)
      real*8  wa(plond,0:plev,ntr)      ! aerosol single scattering albedo
      real*8  ga(plond,0:plev,ntr)      ! aerosol assymetry parameter
      real*8  fa(plond,0:plev,ntr)      ! aerosol forward scattered fraction

      real*8  utaer(plond,ntr)
      real*8  tauaer(plond,ntr)         ! total column aerosol extinction
      real*8  aertau(plond,ntr)         ! total column aerosol extinction, vis band
      real*8  aerprf(plond,0:plev,ntr)  ! total column aerosol extinction, vis band
      real*8  tauprf                ! total column aerosol extinction, vis band
      real*8  waer(plond,ntr)           ! aerosol single scattering albedo
      real*8  gaer(plond,ntr)           ! aerosol asymmetry parameter
      real*8  faer(plond,ntr)           ! aerosol forward scattering fraction

      real*8  path
      real*8  aervtot(plond,plevr),  aermtot(plond,plevr)

      integer ns,i,j,k,itr,ibin,i1,i2,i3,i4

      real*8  rhoso4,rhobc,rhooc,rhodust
      parameter ( rhoso4 = 1.76)
      parameter ( rhobc = 1.)
      parameter ( rhooc = 1.)
      parameter ( rhodust = 2.5)
c options of feedback

c SPECTRAL LOOP
Ctrapuv
      do itr=1,ntr
        do k=0,plev
          do i=1,plond
            tauxar(i,k,itr) = 0.
            wa(i,k,itr) = 0.
            ga(i,k,itr) = 0.
            fa(i,k,itr) = 0.
          enddo
        enddo
      enddo
Ctrapuv_
      DO ns= 1, nspi

c---------------
        do itr=1,ntr
        do i= 1, plond
c set above top value
          uaer(i,0,itr) = 0.
          tauxar(i,0,itr) = 0.
          wa(i,0,itr)      = 0.
          ga(i,0,itr)      = 0.
          fa(i,0,itr)      = 0.
c initialize
          utaer(i,itr)    = 0.0
          tauaer(i,itr)   = 0.0
          waer(i,itr)     = 0.0
          gaer(i,itr)     = 0.0
          faer(i,itr)     = 0.0
        enddo
        enddo

        do i= 1, plond
          tauxar_mix_cs(i,ns) =0.
          tauasc_mix_cs(i,ns) =0.
          gtota_mix_cs(i,ns) =0.
          ftota_mix_cs(i,ns) =0.
        enddo

        do k=0,plev
        do i= 1, plond
          tauxar_mix(i,k,ns) =0.
          tauasc_mix(i,k,ns) =0.
          gtota_mix(i,k,ns) =0.
          ftota_mix(i,k,ns) =0.
        end do
        end do

        IF ( idirect .ge.1) THEN
c calculate optical properties of each aerosol component

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                         OPTION 1  MELANGE EXTERNE
          if (mixtype==1) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          do k=1,plev
            do i=1,plond
              path = (pint(i,k+1) - pint(i,k))/gravit

              ibin = 0
              do itr=1,ntr
                uaer(i,k,itr)=0.
                if (rh(i,k).lt.0.0.or.rh(i,k).gt.1.0) then
                  print*,i,k,rh(i,k),'  RH WARNING !!!!!'
                endif

                IF( chtrname(itr) .eq. 'DUST') then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path

c gaffe au facteur !!
c gaffe au ntr/bins

                  ibin = ibin+1
                  if (ibin.gt.4) print*,'DUST OP PBLEME !!!!'

                  tauxar(i,k,itr) = 1.e5*uaer(i,k,itr)* ksdust(ns,ibin)
                  wa(i,k,itr) = wsdust(ns,ibin)
                  ga(i,k,itr) = gsdust(ns,ibin)
                  fa(i,k,itr) = gsdust(ns,ibin) * gsdust(ns,ibin)

                ELSE IF (chtrname(itr).eq.'SO4')  THEN 
                  uaer(i,k,itr) = aermmr(i,k,itr)*path

                  tauxar(i,k,itr) = 1.e5*uaer(i,k,itr)*ksbase(ns)* exp(
     $                kscoef(ns,1)
     $                + kscoef(ns,2) / ( rh(i,k) + kscoef(ns,3) )
     $                + kscoef(ns,4) / ( rh(i,k) + kscoef(ns,5) )
     $                )
C

                  wa(i,k,itr) = 1.0 - wsbase(ns) * exp(
     $                wscoef(ns,1)
     $                + wscoef(ns,2) / ( rh(i,k) + wscoef(ns,3) )
     $                + wscoef(ns,4) / ( rh(i,k) + wscoef(ns,5) )
     $                )
C
                  ga(i,k,itr) = gsbase(ns) * exp ( gscoef(ns,1)
     $                + gscoef(ns,2) / ( rh(i,k) + gscoef(ns,3) )
     $                + gscoef(ns,4) / ( rh(i,k) + gscoef(ns,5) )
     $                )
C
                  fa(i,k,itr)     = ga(i,k,itr)*ga(i,k,itr)

                ELSE IF ( chtrname(itr).eq.'OC_HL' ) THEN

                  uaer(i,k,itr) = aermmr(i,k,itr)*path
c Humidity effect !
                  tauxar(i,k,itr) = 1.e5 *uaer(i,k,itr)*ksoc_hl(ns)
     &                              *(1-rh(i,k))**(-0.2)
                  wa(i,k,itr) =  wsoc_hl(ns)
                  ga(i,k,itr) = gsoc_hl(ns)
                  fa(i,k,itr)     = ga(i,k,itr)*ga(i,k,itr)

                ELSE IF(chtrname(itr).eq.'BC_HL') THEN

                  uaer(i,k,itr) = aermmr(i,k,itr)*path
c Humidity effect !
                  tauxar(i,k,itr) = 1.e5 *uaer(i,k,itr)*ksbc_hl(ns)
     &                            *(1-rh (i,k) )**(-0.25)
                  wa(i,k,itr) =  wsbc_hl(ns)
                  ga(i,k,itr) = gsbc_hl(ns)
                  fa(i,k,itr)     = ga(i,k,itr)*ga(i,k,itr)
C
                ELSE IF(chtrname(itr).eq.'OC_HB' ) THEN

                  uaer(i,k,itr) = aermmr(i,k,itr)*path
                  tauxar(i,k,itr) = 1.e5 * uaer(i,k,itr) * ksoc_hb(ns)
                  wa(i,k,itr) = wsoc_hb(ns)               
                  ga(i,k,itr) = gsoc_hb(ns)
                  fa(i,k,itr) = gsoc_hb(ns) * gsoc_hb(ns)

                ELSE IF( chtrname(itr).eq.'BC_HB') THEN

                  uaer(i,k,itr) = aermmr(i,k,itr)*path
c Absorbing aerosols (soot type)
                  tauxar(i,k,itr) = 1.e5 * uaer(i,k,itr) * ksbc_hb(ns)
                  wa(i,k,itr) = wsbc_hb(ns)
                  ga(i,k,itr) = gsbc_hb(ns)
                  fa(i,k,itr) = gsbc_hb(ns) * gsbc_hb(ns)

                END IF
              end do  ! end tracer loop 
            end do
          end do              
      
C optical properties for the clear sky diagnostic 

          do i=1,plond
            do itr=1,ntr  
              do k=1,plev
                utaer(i,itr)  = utaer(i,itr)  + uaer(i,k,itr)
                tauaer(i,itr) = tauaer(i,itr) + tauxar(i,k,itr)
                waer(i,itr)   = waer(i,itr) + wa(i,k,itr)*uaer(i,k,itr)
                gaer(i,itr)   = gaer(i,itr) + ga(i,k,itr)*uaer(i,k,itr)
                faer(i,itr)   = faer(i,itr) + fa(i,k,itr)*uaer(i,k,itr)
              end do
              if( utaer(i,itr).le. 1.e-10 ) utaer(i,itr)= 1.e-10
              waer(i,itr)     = waer(i,itr) / utaer(i,itr)
              gaer(i,itr)     = gaer(i,itr) / utaer(i,itr)
              faer(i,itr)     = faer(i,itr) / utaer(i,itr)
            end do
          end do

cCalculate the EXTERNAL Mixing of aerosols 
c  melange externe               
     
          do i=1,plond
            do itr=1,ntr
c             only for climatic feedback allowed             

              do k=0,plev
                tauxar_mix(i,k,ns)= tauxar_mix(i,k,ns)+ tauxar(i,k,itr)
                tauasc_mix(i,k,ns) = tauasc_mix(i,k,ns) +  
     &                           tauxar(i,k,itr) * wa(i,k,itr) 
                gtota_mix(i,k,ns)  = gtota_mix(i,k,ns) + 
     &                     ga(i,k,itr) * tauxar(i,k,itr) * wa(i,k,itr)
                ftota_mix(i,k,ns)  =  ftota_mix(i,k,ns) +
     &                    fa(i,k,itr) * tauxar(i,k,itr) *  wa(i,k,itr)
              end do

c Clear sky (always calcuated if idirect >=1 for diagnostic radiative forcing)

              tauxar_mix_cs(i,ns) =tauxar_mix_cs(i,ns) + tauaer(i,itr)
              tauasc_mix_cs(i,ns) =  tauasc_mix_cs(i,ns)  +
     &                       tauaer(i,itr) * waer(i,itr)
              gtota_mix_cs(i,ns)  = gtota_mix_cs(i,ns) +
     &                       gaer(i,itr) * tauaer(i,itr)* waer(i,itr)
              ftota_mix_cs(i,ns)  =  ftota_mix_cs(i,ns)  +
     &                       faer(i,itr) * tauaer(i,itr) *waer(i,itr)
            end do
          end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                        FIN OPTION 1  MELANGE EXTERNE
c                            OPTION 2  MELANGE INTERNE
          elseif (mixtype ==2) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            do i=1,plond
              do k=1,plev
                path = (pint(i,k+1) - pint(i,k))/gravit
                if (rh(i,k).lt.0.0.or.rh(i,k).gt.1.0)
     &             print*,i,k,rh(i,k),'  RH WARNING !!!!!'

c  sum of hydrophilic aerosols
                aervtot(i,k)=0.
                aermtot(i,k)=0.

                if(iso4.ne.0) then
                  aervtot(i,k)= aervtot(i,k)+ aermmr(i,k,iso4)/rhoso4
                  aermtot(i,k)= aermtot(i,k)+ aermmr(i,k,iso4)
                end if

                if(ibchl.ne.0)then
                  aervtot(i,k)= aervtot(i,k)+ aermmr(i,k,ibchl)/rhobc
                  aermtot(i,k)= aermtot(i,k)+ aermmr(i,k,ibchl)
                end if

                if(iochl.ne.0)then
                  aervtot(i,k)= aervtot(i,k)+ aermmr(i,k,iochl)/rhooc
                  aermtot(i,k)= aermtot(i,k)+ aermmr(i,k,iochl)
                end if

                if(idust(1).ne.0) then
                  aervtot(i,k)=aervtot(i,k)+aermmr(i,k,idust(1))/rhodust
                  aermtot(i,k)=aermtot(i,k)+aermmr(i,k,idust(1))
                end if

c minimum quantity of total aerosol 

                if ( aermtot(i,k) .gt. 1.E-14) then  
c indexes in the internal mixing table
                  prop(1)=(aermmr(i,k,iso4) / rhoso4) / aervtot(i,k)
                  prop(2)=(aermmr(i,k,ibchl) / rhobc) / aervtot(i,k)
                  prop(3)=(aermmr(i,k,iochl) / rhooc) / aervtot(i,k)
                  prop(4)=(aermmr(i,k,idust(1)) /rhodust) /aervtot(i,k)
                  frac(1)=fraction(prop(1))
                  frac(2)=fraction(prop(2))
                  frac(3)=fraction(prop(3))
                  frac(4)=fraction(prop(4))
C FIND THE GREATEST FRACTIONAL PART	   

                  if (iso4 .ne. 0) then
                    i1 = nint (10 * prop(1) ) +1
                  else
                    i1=0+1
                  end if

                  if (ibchl .ne. 0) then
                    i2 = nint (10 * prop(2) ) +1
                  else
                    i2=0+1
                  end if

                  if (iochl .ne. 0) then
                    i3 = nint (10 * prop(3) ) +1
                  else
                    i3=0+1
                  end if

                  if (idust(1) .ne. 0) then
                    i4 = nint (10 * prop(4) ) +1
                  else
                    i4=0+1
                  end if

c final optical parameters
                  if (i1+i2+i3+i4 .eq. 13 ) i4=i4+1
                  if (i1+i2+i3+i4 .eq. 15) then
                    if (i4.ne.1) i4=i4-1
                  end if   

                  if (i1+i2+i3+i4 .eq. 15) then
                    if (i1.ne.1) i1=i1-1
                  end if   

                  if (i1+i2+i3+i4 .eq. 15) then
                    if (i3.ne.1) i3=i3-1
                  end if   

                  if (i1+i2+i3+i4 .eq. 15) then
                    print*, 'On ne doit jamais arriver ici !!!'
                    STOP
                  end if   

                  if (i1+i2+i3+i4 .ne. 14 ) then
                    print*, i1,i2,i3,i4, i1+i2+i3+i4
                    print*,idust(1),iochl,ibchl,iso4
                    print*,'OC HL', aermmr (i,k,iochl) / rhooc
                    print*,'BC HL', aermmr (i,k,ibchl) / rhobc
                    print*,'SO4', aermmr (i,k,iso4) / rhoso4
                    print*,'DUST', aermmr (i,k,idust(1)) / rhodust
                    print*,'VOL TOT', aervtot(i,k)
                    print*,'OC HL%',10*(aermmr(i,k,iochl)/rhooc)
     &                                /aervtot(i,k)
                    print*,'BC HL%',10*(aermmr(i,k,ibchl)/rhobc)
     &                                /aervtot(i,k)
                    print*,'SO4 %',10*(aermmr(i,k,iso4)/rhoso4)
     &                                /aervtot(i,k)
                    print*,'SO4 %',nint(10*(aermmr(i,k,iso4)/rhoso4)
     &                                    /aervtot(i,k))
                    print*,'DUST %',10*(aermmr (i,k,idust(1))/rhodust)
     &                                /aervtot(i,k)
                  print*,'DUST %',nint(10*(aermmr(i,k,idust(1))/rhodust)
     &                                   /aervtot(i,k))
                    stop
                  end if
                  tauxar_mix(i,k,ns) = dextmix(1,ns,i4,i2,i3,i1)*
     &                                 aermtot(i,k)*path* 1e5
                  tauasc_mix(i,k,ns) = dssamix(1,ns,i4,i2,i3,i1)*
     &                                  tauxar_mix(i,k,ns)
                  gtota_mix(i,k,ns)  = dgmix(1, ns, i4,i2,i3,i1)*
     &                                 tauasc_mix(i,k,ns)*
     &                                 tauxar_mix(i,k,ns)
                  ftota_mix(i,k,ns) =  dgmix(1, ns, i4,i2,i3,i1)*
     &                                 dgmix(1, ns, i4,i2,i3,i1)*
     &                                 tauasc_mix(i,k,ns)*
     &                                 tauxar_mix(i,k,ns)

C clear sky dignostic ( 

                  utaer(i,1)  = utaer(i,1) +  aermtot(i,k)*path
        
                  tauaer(i,1)= tauaer(i,1)+ dextmix( 1,ns,i4,i2,i3,i1)*
     &                          aermtot(i,k)*path *1e5           
                  waer(i,1) = waer(i,1) + 
     &                dssamix(1,ns,i4,i2,i3,i1)*aermtot(i,k)*path
                  gaer(i,1)   = gaer(i,1) +
     &                 dgmix (1, ns, i4,i2,i3,i1)* aermtot(i,k)*path
                  faer(i,1)   =gaer(i,1) +
     &                 dgmix(1, ns, i4,i2,i3,i1)*
     &                 dgmix(1, ns, i4,i2,i3,i1)*
     &                 aermtot(i,k)*path

                end if ! end minimum concentration conditions
              end do ! end k loop

              if( utaer(i,1).gt. 1.e-12 ) then 
                waer(i,1)     = waer(i,1) / utaer(i,1)
                gaer(i,1)     = gaer(i,1) / utaer(i,1)
                faer(i,1)     = faer(i,1) / utaer(i,1)  
              end if
c clear sky final effective optical properties

              tauxar_mix_cs(i,ns)= tauaer(i,1)
              tauasc_mix_cs(i,ns)= waer(i,1) * tauaer(i,1)
              gtota_mix_cs(i,ns)  = gaer(i,1) * waer(i,1) * tauaer(i,1)
              ftota_mix_cs(i,ns)  = faer(i,1) * waer(i,1) * tauaer(i,1)
            end do ! end i loop
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                         FIN OPTION 2  MELANGE INTERNE
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        END IF
c end spectral loop
      END DO

      return
      end
