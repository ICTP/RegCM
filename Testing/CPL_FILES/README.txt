-----------------------------------------------
Notes about Usage of Coupled Model (RegCM-ROMS)
-----------------------------------------------
------------------------
Limitations
------------------------

The current version of the coupled model does not support,

- Sub-BATS
- Nesting in ROMS model
- RegCM and CLM combination (only support BATS)
- ROMS-ICE version must be used without BULK_FLUX option
- ROMS model grid cannot be outside of the RegCM model grid.
- RegCM BAND option is not supported. There is a plan to add identical mesh
  support for BAND configuration in the future releases.

also,

- Coupled model works with ROMS 3.5 (with patch)
- Use ESMF 5.2.0rp3 beta as ESMF library. ESMF 5.2.0r2 has bug and it will not
  work with coupled model. The ESMF 5.2.0rp3 can be retrieved using following
  commands.
  
  1) cvs -d:pserver:anonymous@esmf.cvs.sourceforge.net:/cvsroot/esmf login
     (when asked for a password simply press enter)

  2) cvs -z3 -d:pserver:anonymous@esmf.cvs.sourceforge.net:/cvsroot/esmf
     co -P -rESMF_5_2_0rp3_beta_snapshot_02 esmf 

- User needs to provide atmospheric forcing files to the ROMS model.
  The dummy.ncl can be used to create dummy forcing file for this
  purposes.

------------------------
Ocean Model Setup (ROMS)
------------------------

1 - User must set following CPP flag in the configuration file
    to activate coupled model support.

    #define REGCM_COUPLING

2 - Following options is used to provide correct forcing information 
    to ocean component.

    #define BULK_FLUXES
    #define SHORTWAVE
    #define SPECIFIC_HUMIDITY

    or user can choose not to use BULK_FLUXES. In this case, ICE 
    enabled ROMS code can not be used. We are still working on to enable
    ICE model without BULK_FLUXES option.

    Be aware that the heat fluxes might not conserved with BULK_FLUXES.
    If you want to used conserved version of the coupled model, compile
    ROMS without BULK_FLUXES option.

3 - SPECIFIC_HUMIDITY option must be used because mixing ratio is in
    kg/kg unit in RegCM side.

4 - Do not specify any option related with the processing of longwave 
    radiation in the configuration file, unless it is really necessary. 
    By default, ROMS uses net longwave flux provided by RegCM. If you 
    set the longwave flux option as,

    #define LONGWAVE_OUT

    Be aware that there might be a problem in the conservation of heat 
    fluxes between model components if this option is enabled.

5 - In grid file, value of the "spherical" variable must be set as "T".
    This is crucial to create grid definition in ESMF side.

6 - If user want to restart coupled model. The restart information must
    be written in daily basis in ROMS side. The ./restart.sh can help to
    restart coupled model.

---
FAQ
---

1 - How can i apply the patch for coupled model?

    cd roms-3.5
    cp ../regcm/Testing/CPL_FILES/roms-3.5.patch .
    patch -p1 -i roms-3.5.patch

2 - How can i download ROMS? (assumed already registered to ROMS site)

    svn co --username [ROMS_USER_NAME] https://www.myroms.org/svn/src/tags/roms-3.5 roms-3.5
     
