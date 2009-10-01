      subroutine aerout(jslc,   tauxar_mix,  tauasc_mix, gtota_mix,
     &                 aeradfo, aeradfos) 

      implicit none
      integer jslc,i,k

#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param2.cb'
      include 'Commons/date.cb'
#ifdef MPP1
      include 'Common2/trachem.cb'
#else
      include 'Commons/trachem.cb'
#endif
      include 'Commons/aeroppt.cb'
      include 'Commons/parrad.cb'

      real*8  tauxar_mix(plond,0:plev,nspi)  ! aerosol extinction optical depth
      real*8  tauasc_mix(plond,0:plev,nspi)
      real*8  gtota_mix(plond,0:plev,nspi)
      real*8  aeradfo(plond) ,   aeradfos(plond) 
      integer ntim

c      Viible cas ns=8

      do k=1,kx
        do i=2,ix-1
#ifdef MPP1
          aerext(i-1,k,jslc)= tauxar_mix(i,k,8)        
          aerssa(i-1,k,jslc)= tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc)=  gtota_mix(i,k,8) 
#else
          aerext(i-1,k,jslc-1)= tauxar_mix(i,k,8)        
          aerssa(i-1,k,jslc-1)= tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc-1)=  gtota_mix(i,k,8) 
#endif
        enddo
      enddo
c CARE :Average the radiative forcing between chem output time steps (in
c hour) according to radfrq (in min), aertarf is reset to 0 at each chem
c output (cf output.f)
c
      ntim = 60 * chemfrq / radfrq
c
c aersol radative forcing (care cgs to mks after radiation scheme !)        
      do i=2,ix-1
#ifdef MPP1
        aertarf(i-1,jslc)= aertarf(i-1,jslc)+ aeradfo(i) *1.e-3/ntim
        aersrrf(i-1,jslc)= aersrrf(i-1,jslc)+aeradfos(i)*1.e-3/ntim
#else
        aertarf(i-1,jslc-1)= aertarf(i-1,jslc-1)+ aeradfo(i) *1.e-3/ntim
        aersrrf(i-1,jslc-1)= aersrrf(i-1,jslc-1)+aeradfos(i)*1.e-3/ntim
#endif
      enddo   

      return
      end
