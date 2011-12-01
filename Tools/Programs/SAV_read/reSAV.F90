program resav

  use mod_dynparam

  implicit none
!
  integer , parameter :: iutrst = 14
!
  character(256) :: namelistfile , prgname , savfile
  integer :: ierr
  type (rcm_time_and_date) :: idatex
  real(8) :: xbctime
  integer :: ktau , ntime
  real(8) , allocatable , dimension(:,:,:) :: ub0 , vb0 , qb0 , tb0 , so0
  real(8) , allocatable , dimension(:,:) :: ps0 , ts0
  real(8) , allocatable , dimension(:,:,:) :: ua , va , ta  , qva , qca
  real(8) , allocatable , dimension(:,:,:) :: ub , vb , tb  , qvb , qcb
  real(8) , allocatable , dimension(:,:) :: psa , psb
  real(8) , allocatable , dimension(:,:) :: tga , tgb
  real(8) , allocatable , dimension(:,:) :: rainc , rainnc
  real(8) , allocatable , dimension(:,:,:) :: rsheat , rswat , tbase
  real(8) , allocatable , dimension(:,:) :: cldefi , cbmf2d
  real(8) , allocatable , dimension(:,:) :: hfx , qfx , uvdrag
#ifndef BAND
  real(8) :: tdini , tdadv , tqini , tqadv , tqeva , tqrai
#endif
  real(8) , allocatable , dimension(:,:,:,:) :: absnxt , abstot
  real(8) , allocatable , dimension(:,:,:) :: emstot
  real(8) , allocatable , dimension(:,:,:) :: fcc
#ifdef CLM
  real(8) , allocatable , dimension(:,:) :: sols2d , soll2d , solsd2d , &
        solld2d , aldirs2d , aldirl2d , aldifs2d , aldifl2d , lndcat2d
#endif
  real(8) , allocatable , dimension(:,:) :: sol2d , solvd2d , solvs2d , &
        sabv2d , flw2d , flwd2d , fsw2d , sinc2d , pptc , pptnc , &
        prca2d , prnca2d , ssw2da , sdelqk2d , sdeltk2d , sfracb2d , &
        sfracs2d , sfracv2d , svegfrac2d
  real(8) , allocatable , dimension(:,:,:) :: tlef2d , ssw2d , srw2d , &
        tg2d , tgb2d , scv2d , gwet2d , sag2d , sice2d , dew2d , ircp2d , &
        taf2d , emiss2d
  integer , allocatable , dimension(:,:) :: veg2d , ldmsk
  integer , allocatable , dimension(:,:,:) :: veg2d1 , ocld2d
  real(8) , allocatable , dimension(:,:,:) :: heatrt , o3prof , swt2d
  real(8) , allocatable , dimension(:,:) :: tgbb , zpbl
  real(8) , allocatable , dimension(:,:,:,:) :: chia , chib
  real(8) , allocatable , dimension(:,:,:,:) :: remlsc , remcvc
  real(8) , allocatable , dimension(:,:,:) :: remdrd
  real(8) , allocatable , dimension(:) :: tchiad , tchitb , tchie
  real(8) , allocatable , dimension(:,:,:) :: dstor , hstor
  real(8) , allocatable , dimension(:,:) :: uj1 , uj2 , ujlx , ujl
  real(8) , allocatable , dimension(:,:) :: vj1 , vj2 , vjlx , vjl
  real(8) , allocatable , dimension(:,:) :: ui1 , ui2 , uilx , uil
  real(8) , allocatable , dimension(:,:) :: vi1 , vi2 , vilx , vil
  integer :: ibltyp , iboudy , igcc , ichem , icup , iocnflx , &
             ipptls , ipgf , iemiss , lakemod , idcsst , iseaice , &
             idesseas , iconvlwp , iocnrough
  character(3) :: scenario
!
  namelist /physicsparam/ ibltyp , iboudy , icup , igcc , ipgf ,    &
    iemiss , lakemod , ipptls , iocnflx , iocnrough , ichem,        &
    scenario , idcsst , iseaice , idesseas , iconvlwp
!
! Read input global namelist
!
  ibltyp = 1
  iboudy = 1
  icup = 1
  ipptls = 1
  igcc = 1
  ipgf = 1
  iemiss = 1
  iocnflx = 1
  lakemod = 0
  ichem = 0
  scenario = 'A1B'
  idcsst = 0
  iseaice = 0
  idesseas = 1
  high_nudge = 3.0D0
  medium_nudge = 2.0D0
  low_nudge = 1.0D0
  iconvlwp = 1

  call getarg(0, prgname)
  call getarg(1, namelistfile)
  call getarg(2, savfile)
  call initparam(namelistfile, ierr)
  if ( ierr/=0 ) then
    write ( 6, * ) 'Parameter initialization not completed'
    write ( 6, * ) 'Usage : '
    write ( 6, * ) ' ', trim(prgname), ' regcm.in SAV_file'
    write ( 6, * ) ' '
    write ( 6, * ) 'Check * argument * and * namelist * syntax'
    stop
  end if
  read (ipunit, physicsparam)
  close (ipunit)
!
! Allocate space
!
  allocate(ub0(iy,kz,jx))
  allocate(vb0(iy,kz,jx))
  allocate(qb0(iy,kz,jx))
  allocate(tb0(iy,kz,jx))
  if ( ehso4 ) then
    allocate(so0(iy,kz,jx))
  end if
  allocate(ps0(iy,jx))
  allocate(ts0(iy,jx))
!
  allocate(ua(iy,kz,jx))
  allocate(va(iy,kz,jx))
  allocate(ta(iy,kz,jx))
  allocate(qva(iy,kz,jx))
  allocate(qca(iy,kz,jx))
  allocate(ub(iy,kz,jx))
  allocate(vb(iy,kz,jx))
  allocate(tb(iy,kz,jx))
  allocate(qvb(iy,kz,jx))
  allocate(qcb(iy,kz,jx))
  allocate(psa(iy,jx))
  allocate(psb(iy,jx))
  allocate(tga(iy,jx))
  allocate(tgb(iy,jx))
  allocate(rainc(iy,jx))
  allocate(rainnc(iy,jx))
  if ( icup == 1 ) then
    allocate(rsheat(iy,kz,jx))
    allocate(rswat(iy,kz,jx))
  end if
  if ( icup == 3 ) then
    allocate(tbase(iy,kz,jx))
    allocate(cldefi(iy,jx))
  end if
  if ( icup == 4 .or. icup == 99 .or. icup == 98 ) then
    allocate(cbmf2d(iy,jx))
  end if
  allocate(hfx(iy,jx))
  allocate(qfx(iy,jx))
  allocate(uvdrag(iy,jx))
  allocate(tgbb(iy,jx))
  allocate(zpbl(iy,jx))
#ifdef BAND
  allocate(absnxt(iym1,kz,4,jx))
  allocate(abstot(iym1,kzp1,kz + 1,jx))
  allocate(emstot(iym1,kzp1,jx))
#else
  allocate(absnxt(iym1,kz,4,jxm1))
  allocate(abstot(iym1,kzp1,kzp1,jxm1))
  allocate(emstot(iym1,kzp1,jxm1))
#endif
  if ( ipptls == 1 ) allocate(fcc(iy,kz,jx))
#ifdef CLM
#ifdef BAND
  allocate(sols2d(iym1,jx))
  allocate(soll2d(iym1,jx))
  allocate(solsd2d(iym1,jx))
  allocate(solld2d(iym1,jx))
  allocate(aldirs2d(iym1,jx))
  allocate(aldirl2d(iym1,jx))
  allocate(aldifs2d(iym1,jx))
  allocate(aldifl2d(iym1,jx))
  allocate(lndcat2d(iym1,jx))
#else
  allocate(sols2d(iym1,jxm1))
  allocate(soll2d(iym1,jxm1))
  allocate(solsd2d(iym1,jxm1))
  allocate(solld2d(iym1,jxm1))
  allocate(aldirs2d(iym1,jxm1))
  allocate(aldirl2d(iym1,jxm1))
  allocate(aldifs2d(iym1,jxm1))
  allocate(aldifl2d(iym1,jxm1))
  allocate(lndcat2d(iym1,jxm1))
#endif
#endif
#ifdef BAND
  allocate(sol2d(iym1,jx))
  allocate(solvd2d(iym1,jx))
  allocate(solvs2d(iym1,jx))
  allocate(sabv2d(iym1,jx))
  allocate(tlef2d(nnsg,iym1,jx))
  allocate(ssw2d(nnsg,iym1,jx))
  allocate(srw2d(nnsg,iym1,jx))
  allocate(tg2d(nnsg,iym1,jx))
  allocate(tgb2d(nnsg,iym1,jx))
  allocate(scv2d(nnsg,iym1,jx))
  allocate(gwet2d(nnsg,iym1,jx))
  allocate(sag2d(nnsg,iym1,jx))
  allocate(sice2d(nnsg,iym1,jx))
  allocate(dew2d(nnsg,iym1,jx))
  allocate(ircp2d(nnsg,iym1,jx))
  allocate(veg2d(iym1,jx))
  allocate(ldmsk(iym1,jx))
  allocate(veg2d1(nnsg,iym1,jx))
  allocate(heatrt(ym1,kz,jx))
  allocate(o3prof(iym1,kzp1,jx))
  allocate(flw2d(iym1,jx))
  allocate(flwd2d(iym1,jx))
  allocate(fsw2d(iym1,jx))
  allocate(swt2d(nnsg,iym1,jx))
  allocate(sinc2d(iym1,jx))
  allocate(taf2d(nnsg,iym1,jx))
  allocate(ocld2d(nnsg,iym1,jx))
  allocate(emiss2d(nnsg,iym1,jx))
  allocate(pptnc(iym1,jx))
  allocate(pptc(iym1,jx))
  allocate(prca2d(iym1,jx))
  allocate(prnca2d(iym1,jx))
  if ( ichem == 1 ) then
    allocate(ssw2da(iym1,jx))
    allocate(sdelqk2d(iym1,jx))
    allocate(sdeltk2d(iym1,jx))
    allocate(sfracb2d(iym1,jx))
    allocate(sfracs2d(iym1,jx))
    allocate(sfracv2d(iym1,jx))
    allocate(svegfrac2d(iym1,jx))
  end if
#else
  allocate(sol2d(iym1,jxm1))
  allocate(solvd2d(iym1,jxm1))
  allocate(solvs2d(iym1,jxm1))
  allocate(sabv2d(iym1,jxm1))
  allocate(tlef2d(nnsg,iym1,jxm1))
  allocate(ssw2d(nnsg,iym1,jxm1))
  allocate(srw2d(nnsg,iym1,jxm1))
  allocate(tg2d(nnsg,iym1,jxm1))
  allocate(tgb2d(nnsg,iym1,jxm1))
  allocate(scv2d(nnsg,iym1,jxm1))
  allocate(gwet2d(nnsg,iym1,jxm1))
  allocate(sag2d(nnsg,iym1,jxm1))
  allocate(sice2d(nnsg,iym1,jxm1))
  allocate(dew2d(nnsg,iym1,jxm1))
  allocate(ircp2d(nnsg,iym1,jxm1))
  allocate(veg2d(iym1,jxm1))
  allocate(ldmsk(iym1,jxm1))
  allocate(veg2d1(nnsg,iym1,jxm1))
  allocate(heatrt(iym1,kz,jxm1))
  allocate(o3prof(iym1,kzp1,jxm1))
  allocate(flw2d(iym1,jxm1))
  allocate(flwd2d(iym1,jxm1))
  allocate(fsw2d(iym1,jxm1))
  allocate(swt2d(nnsg,iym1,jxm1))
  allocate(sinc2d(iym1,jxm1))
  allocate(taf2d(nnsg,iym1,jxm1))
  allocate(ocld2d(nnsg,iym1,jxm1))
  allocate(emiss2d(nnsg,iym1,jxm1))
  allocate(pptnc(iym1,jxm1))
  allocate(pptc(iym1,jxm1))
  allocate(prca2d(iym1,jxm1))
  allocate(prnca2d(iym1,jxm1))
  if ( ichem == 1 ) then
    allocate(ssw2da(iym1,jxm1))
    allocate(sdelqk2d(iym1,jxm1))
    allocate(sdeltk2d(iym1,jxm1))
    allocate(sfracb2d(iym1,jxm1))
    allocate(sfracs2d(iym1,jxm1))
    allocate(sfracv2d(iym1,jxm1))
    allocate(svegfrac2d(iym1,jxm1))
  end if
#endif
  if ( ichem == 1 ) then
    allocate(chia(iy,kz,jx,ntr))
    allocate(chib(iy,kz,jx,ntr))
    allocate(remlsc(iy,kz,jx,ntr))
    allocate(remcvc(iy,kz,jx,ntr))
    allocate(remdrd(iy,jx,ntr))
#ifndef BAND
    allocate(tchiad(ntr))
    allocate(tchitb(ntr))
    allocate(tchie(ntr))
#endif
  end if
  allocate(dstor(iy,jx,nsplit))
  allocate(hstor(iy,jx,nsplit))
#ifndef BAND
  allocate(uj1(iy,kz))
  allocate(uj2(iy,kz))
  allocate(ujl(iy,kz))
  allocate(ujlx(iy,kz))
  allocate(vj1(iy,kz))
  allocate(vj2(iy,kz))
  allocate(vjl(iy,kz))
  allocate(vjlx(iy,kz))
#endif
  allocate(ui1(kz,jx))
  allocate(ui2(kz,jx))
  allocate(uil(kz,jx))
  allocate(uilx(kz,jx))
  allocate(vi1(kz,jx))
  allocate(vi2(kz,jx))
  allocate(vil(kz,jx))
  allocate(vilx(kz,jx))
!
  open (iutrst, file=savfile, form='unformatted',status='old')
  read (iutrst) ktau, xbctime, idatex, ntime
  if ( ehso4 ) then
    read (iutrst) ub0, vb0, qb0, tb0, ps0, ts0, so0
  else
    read (iutrst) ub0, vb0, qb0, tb0, ps0, ts0
  end if
  read (iutrst) ua
  read (iutrst) va
  read (iutrst) ta
  read (iutrst) qva
  read (iutrst) qca
  read (iutrst) ub
  read (iutrst) vb
  read (iutrst) tb
  read (iutrst) qvb
  read (iutrst) qcb
  read (iutrst) psa, psb
  read (iutrst) tga, tgb, rainc, rainnc
  if ( icup == 1 ) then
    read (iutrst) rsheat, rswat
  end if
  if ( icup == 3 ) then
    read (iutrst) tbase, cldefi
  end if
  if ( icup == 4 .or. icup == 99 .or. icup == 98 ) then
    read (iutrst) cbmf2d
  end if
  read (iutrst) hfx, qfx, uvdrag
#ifndef BAND
  read (iutrst) tdini , tdadv , tqini , tqadv , tqeva , tqrai
#endif
  read (iutrst) absnxt, abstot, emstot
  if ( ipptls == 1 ) read (iutrst) fcc
#ifdef CLM
  read (iutrst) sols2d
  read (iutrst) soll2d
  read (iutrst) solsd2d
  read (iutrst) solld2d
  read (iutrst) aldirs2d
  read (iutrst) aldirl2d
  read (iutrst) aldifs2d
  read (iutrst) aldifl2d
  read (iutrst) lndcat2d
#endif
  read (iutrst) sol2d
  read (iutrst) solvd2d
  read (iutrst) solvs2d
  read (iutrst) sabv2d
  read (iutrst) tlef2d
  read (iutrst) ssw2d
  read (iutrst) srw2d
  read (iutrst) tg2d
  read (iutrst) tgb2d
  read (iutrst) scv2d
  read (iutrst) gwet2d
  read (iutrst) sag2d
  read (iutrst) sice2d
  read (iutrst) dew2d
  read (iutrst) ircp2d
  read (iutrst) veg2d
  read (iutrst) ldmsk
  read (iutrst) veg2d1
  read (iutrst) heatrt
  read (iutrst) o3prof
  read (iutrst) tgbb
  read (iutrst) flw2d
  read (iutrst) flwd2d
  read (iutrst) fsw2d
  read (iutrst) swt2d
  read (iutrst) sinc2d
  read (iutrst) taf2d
  read (iutrst) ocld2d
  read (iutrst) emiss2d
  read (iutrst) pptnc, pptc, prca2d, prnca2d
  if ( iocnflx == 2 ) read (iutrst) zpbl
  if ( ichem == 1 ) then
    read (iutrst) chia
    read (iutrst) chib
!   cumul removal terms (3d, 2d)
    read (iutrst) remlsc
    read (iutrst) remcvc
    read (iutrst) remdrd
!   cumul ad, dif, emis terms ( scalar)
    read (iutrst) ssw2da
    read (iutrst) sdeltk2d
    read (iutrst) sdelqk2d
    read (iutrst) sfracv2d
    read (iutrst) sfracb2d
    read (iutrst) sfracs2d
    read (iutrst) svegfrac2d
#ifndef BAND
    read (iutrst) tchiad
    read (iutrst) tchitb
    read (iutrst) tchie
  end if
#endif
  read (iutrst) dstor
  read (iutrst) hstor
#ifndef BAND
  read (iutrst) uj1, uj2, ujlx, ujl
#endif
  read (iutrst) ui1, ui2, uilx, uil
#ifndef BAND
  read (iutrst) vj1, vj2, vjlx, vjl
#endif
  read (iutrst) vi1, vi2, vilx, vil
  close(iutrst)

  print *, 'Done'

end program resav
