!
! OASIS parameters
!
 &oasisparam
 write_restart_option = 2, ! 0 => no restart file writing
                           ! 1 => write restart files at the first time
                           ! 2 => write restart files at the last time
                           ! 3 => both 1 & 2
 l_write_grids = .true.,   ! For writing grids.nc, areas.nc, masks.nc.
 oasis_sync_lag = 0,       ! Synchronisation lag with other components (sec)
                           ! > 0 => Regcm starts late
                           ! < 0 => Regcm starts in advance
                           ! Should be a multipe of the coupling period
                           !   in the namcouple.
                           ! In the namcouple, RUNTIME must have
                           !   the run duration + oasis_sync_lag.
                           !------ NAMCOUPLE FIELD ENTRIES ------
                           ! field    | grid
                           !-------------------------------------
 l_cpl_im_sst  = .false.,  ! RCM_SST  | rcim     
 l_cpl_im_wz0  = .false.,  ! RCM_WZ0  | rcim     
 l_cpl_im_wust = .false.,  ! RCM_WUST | rcim     
 l_cpl_ex_u10m = .false.,  ! RCM_U10M | rcin/rcim
 l_cpl_ex_v10m = .false.,  ! RCM_V10M | rcin/rcim
 l_cpl_ex_wspd = .false.,  ! RCM_WSPD | rcin/rcim
 l_cpl_ex_wdir = .false.,  ! RCM_WDIR | rcin/rcim
 l_cpl_ex_t2m  = .false.,  ! RCM_T2M  | rcin/rcim
 l_cpl_ex_q2m  = .false.,  ! RCM_Q2M  | rcin/rcim
 l_cpl_ex_slp  = .false.,  ! RCM_SLP  | rcen/rcem
 l_cpl_ex_taux = .false.,  ! RCM_TAUX | rcin/rcim
 l_cpl_ex_tauy = .false.,  ! RCM_TAUY | rcin/rcim
 l_cpl_ex_z0   = .false.,  ! RCM_Z0   | rcin/rcim
 l_cpl_ex_ustr = .false.,  ! RCM_USTR | rcin/rcim
 l_cpl_ex_evap = .false.,  ! RCM_EVAP | rcin/rcim
 l_cpl_ex_prec = .false.,  ! RCM_PREC | rcin/rcim
 l_cpl_ex_nuwa = .false.,  ! RCM_NUWA | rcin/rcim
 l_cpl_ex_ulhf = .false.,  ! RCM_ULHF | rcin/rcim
 l_cpl_ex_ushf = .false.,  ! RCM_USHF | rcin/rcim
 l_cpl_ex_uwlw = .false.,  ! RCM_UWLW | rcin/rcim
 l_cpl_ex_dwlw = .false.,  ! RCM_DWLW | rcin/rcim
 l_cpl_ex_nulw = .false.,  ! RCM_NULW | rcin/rcim
 l_cpl_ex_uwsw = .false.,  ! RCM_UWSW | rcin/rcim
 l_cpl_ex_dwsw = .false.,  ! RCM_DWSW | rcin/rcim
 l_cpl_ex_ndsw = .false.,  ! RCM_NDSW | rcin/rcim
 l_cpl_ex_rhoa = .false.,  ! RCM_RHOA | rcin/rcim
                           !------ NAMCOUPLE FIELD ENTRIES ------
 /
