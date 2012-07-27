-----------------------------------------------
Notes about Usage of Coupled Model (RegCM-ROMS)
-----------------------------------------------
------------------------
Limitations
------------------------

The current version of the coupled model does not support,

- Nesting in both RegCM (sub BATs) and ROMS side
- RegCM and CLM combination (only support BATs)

also,

- Coupled model works with ROMS 3.5 (with patch)
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
