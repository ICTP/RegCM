#!/bin/csh -f                                     
set mydir=$PWD                                    
cd ../Main                                        
make clean                                        
./MAKECODE                                        
make                                              
cd $mydir                                         
mv ../Main/regcm .                                
/bin/ln -sf ../Input/DOMAIN.INFO fort.10          
/bin/ln -sf    ../Input/ICBC1990010100 fort.101   
./regcm<./regcm.in                                
