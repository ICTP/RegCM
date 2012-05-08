#! /bin/tcsh

#set echo verbose
###############################################################################
# This is the master script for pre-processing the CMIP5 emissions
###############################################################################

set start_date=2000
set end_date=2006

##  SPECIFICATIONS ##
set data_dir=/home/netapp/clima-users/users/apozzer/emissions/data/CMIP5
echo $data_dir
@ nscenario=1
set scenario=('RCP26'  'RCP45'  'RCP60'  'RCP85')

#
###############################################################################

@ nspecs = 19 
#species name (output): CBMZ species
set spec = ( 'BC'  'OC'  'NH3' 'NO' 'CO' 'CH4'  'SO2' 'C2H6' 'PAR1' 'PAR2' 'ACET' 'HCHO' 'ETHE' 'OLT' 'OLI' 'ALD2' 'TOL' 'XYL' 'MOH'  )
#species name (input): CMIP emission in data_dir
set ispec = ('BC'  'OC'  'NH3' 'NO' 'CO' 'CH4'  'SO2' 'ethane' 'propane' 'butane' 'ketones' 'formaldehyde' 'ethene' 'propene' 'other_alkenes_and_alkynes' 'other_alkanals' 'toluene' 'xylene' 'alcohols' )


# The molar weights should be in the same order as specs ! the NMVOC are in gr(C)/mol 
#  set spec =('BC'  'OC'  'NH3' 'NO'   'CO'  'CH4'  'SO2'  'HCHO'  'CH3OH' 'C2H6' 'C3H8' 'CH3COCH3' 'CH3CHO'  'CH3COOH' 'HCOOH'   'C2H4'   'C3H6'   'C4H10'    'MEK' )
set mwt =   ( 1.0    1.0   1.0   1.0   1.0    1.0    1.0    1.0  1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 )
set factor_C=1.742
#  set spec =(      'BC'  'OC'  'NH3' 'NO'   'CO'  'CH4'  'SO2'  'HCHO'  'CH3OH' 'C2H6'   'C3H8' 'CH3COCH3' 'CH3CHO'  'CH3COOH' 'HCOOH'   'C2H4'   'C3H6'   'C4H10'    'MEK' )
set scale_bb_tmp = (1.0   1.0    1.0   1.0    1.0    1.0   1.0  1.0 1.0  1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0)
set scale_ff_tmp = (1.0   1.0    1.0   1.0    1.0    1.0   1.0  1.0 1.0  1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0)


@ ntypes = 3
set type = ( 'ff' 'bf' 'bb')

if ( ! -d output   ) mkdir output
if ( ! -d output/nml   ) mkdir output/nml

#link up directory with original files
ln -s $data_dir ./input

# link up tools

ln -sf tools/emcre_v0.5/src/emcre .
ln -sf tools/snr .

# specify location of templates:

set templates=./template

@ nscen =1
while ( $nscen <= $nscenario )

  @ ns = 1
  while ( $ns <= $nspecs )
 
    @ nt = 1
    while ( $nt <= $ntypes )

    echo ${spec[$ns]}_${type[$nt]}
 
       cp ./template/CMIP5_${scenario[$nscen]}_${type[$nt]}.nml  ./CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml
       ./snr CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml output_xyz  ${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nc
       ./snr CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml species_xyz  ${spec[$ns]}
       ./snr CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml input_xyz  ${ispec[$ns]}
       ./snr CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml mm_xyz  ${mwt[$ns]}
       ./snr CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml start_xyz $start_date
       ./snr CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml end_xyz $end_date

       if ( $type[$nt] == 'ff' ) then
          set scale_ff=${scale_ff_tmp[$ns]}
          ./snr CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml scale_xyz $scale_ff
       else if ( $type[$nt] == 'bf' ) then
          set scale_bb=${scale_bb_tmp[$ns]}
          ./snr CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml scale_xyz $scale_bb
       else if ( $type[$nt] == 'bb' ) then
          if (${ispec[$ns]} == "NMVOC") then
            set scale_bb=`echo "scale=4;${scale_bb_tmp[$ns]}/${factor_C}"| bc`
          else
            set scale_bb=${scale_bb_tmp[$ns]}
          endif 
          ./snr CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml scale_xyz $scale_bb
       endif

       ./emcre CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml
       mv CMIP5_${scenario[$nscen]}_${spec[$ns]}_${type[$nt]}.nml ./output

       nextstep:
    
    @ nt++
    end # on nt

   # now merge the files from different types (bf & ff) to one single file for anthropogenic emissions
   # bf is after ff. If bf exist then also ff exist
   ###### LAND ########
   set name_ff = ./output/${scenario[$nscen]}_${spec[$ns]}\_bf\.nc
   set name_bf = ./output/${scenario[$nscen]}_${spec[$ns]}\_ff\.nc
   set name_bb = ./output/${scenario[$nscen]}_${spec[$ns]}\_bb\.nc
   set name_ant = ./output/${scenario[$nscen]}_${spec[$ns]}_ant.nc
   set name_emis = ./output/${scenario[$nscen]}_${spec[$ns]}_emis.nc
   ncbo --op_typ=add $name_bf $name_ff $name_ant
   ncbo --op_typ=add $name_bb $name_ant $name_emis


# remove at the end
   rm -f ./output/${scenario[$nscen]}_${spec[$ns]}_bf.nc ./output/${scenario[$nscen]}_${spec[$ns]}_ff.nc ./output/${scenario[$nscen]}_${spec[$ns]}\_bb\.nc ./output/${scenario[$nscen]}_${spec[$ns]}_ant.nc 
   mv  ./output/*.nml ./output/nml
  
  @ ns++
  end # on ns

# perform some final lumping for CBMZ here 
  set par1 = ./output/${scenario[$nscen]}_PAR1_emis.nc
  set par2 = ./output/${scenario[$nscen]}_PAR2_emis.nc
  set par = ./output/${scenario[$nscen]}_PAR_emis.nc
  ncrename -v PAR1_flux,PAR_flux $par1 toto.nc
  ncrename -v PAR2_flux,PAR_flux $par2 tutu.nc
  ncbo --op_typ=add toto.nc tutu.nc $par
  
rm -f $par1 $par2 tutu.nc toto.nc 

@ nscen++
end # on nscen


             
# clean up links

unlink input
rm emcre
rm snr
             

exit

