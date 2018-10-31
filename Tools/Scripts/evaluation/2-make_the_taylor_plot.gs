***/home/netapp-clima/users/ctorma/Scripts_scratch/Plot/Alps/Taylor/Pre/011/new)

'reinit'

'set display color white'


*'subplot 1 1 1 -xy 1 -tall 1'

index.1='r95'
index.2='sdii'
index.3='ddf'
index.4='cdd'

title.1='a'
title.2='b'
title.3='c'
title.4='d'

panel.1='subpg l4 '1
panel.2='subpg l4 '3
panel.3='subpg l4 '2
panel.4='subpg l4 '4

  'run 'panel.1
'set strsiz 0.15'
'set font 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.1'_ERAINT  -l -3. -m 3 -z 0.15 -t "A" -c 1 -levs 0.5 1. 1.5 2 2.5 3 -levr -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.1'_m_RCM  -l -3. -m 3 -z 0.15 -t "A" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.1'_h_RCM  -l -3. -m 3 -z 0.15 -t "A" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.1'_ERAINT  -l -3. -m 2 -z 0.15 -t "A" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.1'_m_RCM  -l -3. -m 2 -z 0.15 -t "A" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.1'_h_RCM  -l -3. -m 2 -z 0.15 -t "A" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.1'_m_medRCM  -l -3. -m 3 -z 0.15 -t "A" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.1'_m_medRCM  -l -3. -m 2 -z 0.15 -t "A" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.1'_h_medRCM  -l -3. -m 3 -z 0.15 -t "A" -c 7'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.1'_h_medRCM  -l -3. -m 2 -z 0.15 -t "A" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.1'_ERAINT  -l -3. -m 3 -z 0.15 -t "E" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.1'_m_RCM  -l -3. -m 3 -z 0.15 -t "E" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.1'_h_RCM  -l -3. -m 3 -z 0.15 -t "E" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.1'_ERAINT  -l -3. -m 2 -z 0.15 -t "E" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.1'_m_RCM  -l -3. -m 2 -z 0.15 -t "E" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.1'_h_RCM  -l -3. -m 2 -z 0.15 -t "E" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.1'_m_medRCM  -l -3. -m 3 -z 0.15 -t "E" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.1'_m_medRCM  -l -3. -m 2 -z 0.15 -t "E" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.1'_h_medRCM  -l -3. -m 3 -z 0.15 -t "E" -c 7'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.1'_h_medRCM  -l -3. -m 2 -z 0.15 -t "E" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.1'_ERAINT  -l -3. -m 3 -z 0.15 -t "C" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.1'_m_RCM  -l -3. -m 3 -z 0.15 -t "C" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.1'_h_RCM  -l -3. -m 3 -z 0.15 -t "C" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.1'_ERAINT  -l -3. -m 2 -z 0.15 -t "C" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.1'_m_RCM  -l -3. -m 2 -z 0.15 -t "C" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.1'_h_RCM  -l -3. -m 2 -z 0.15 -t "C" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.1'_m_medRCM  -l -3. -m 3 -z 0.15 -t "C" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.1'_m_medRCM  -l -3. -m 2 -z 0.15 -t "C" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.1'_h_medRCM  -l -3. -m 3 -z 0.15 -t "C" -c 7'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.1'_h_medRCM  -l -3. -m 2 -z 0.15 -t "C" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.1'_ERAINT  -l -3. -m 3 -z 0.15 -t "I" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.1'_m_RCM  -l -3. -m 3 -z 0.15 -t "I" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.1'_h_RCM  -l -3. -m 3 -z 0.15 -t "I" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.1'_ERAINT  -l -3. -m 2 -z 0.15 -t "I" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.1'_m_RCM  -l -3. -m 2 -z 0.15 -t "I" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.1'_h_RCM  -l -3. -m 2 -z 0.15 -t "I" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.1'_m_medRCM  -l -3. -m 3 -z 0.15 -t "I" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.1'_m_medRCM  -l -3. -m 2 -z 0.15 -t "I" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.1'_h_medRCM  -l -3. -m 3 -z 0.15 -t "I" -c 7'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.1'_h_medRCM  -l -3. -m 2 -z 0.15 -t "I" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.1'_ERAINT  -l -3. -m 3 -z 0.15 -t "U" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.1'_m_RCM  -l -3. -m 3 -z 0.15 -t "U" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.1'_h_RCM  -l -3. -m 3 -z 0.15 -t "U" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.1'_ERAINT  -l -3. -m 2 -z 0.15 -t "U" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.1'_m_RCM  -l -3. -m 2 -z 0.15 -t "U" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.1'_h_RCM  -l -3. -m 2 -z 0.15 -t "U" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.1'_m_medRCM  -l -3. -m 3 -z 0.15 -t "U" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.1'_m_medRCM  -l -3. -m 2 -z 0.15 -t "U" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.1'_h_medRCM  -l -3. -m 3 -z 0.15 -t "U" -c 7'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.1'_h_medRCM  -l -3. -m 2 -z 0.15 -t "U" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.1'_ERAINT  -l -3. -m 3 -z 0.15 -t "F" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.1'_m_RCM  -l -3. -m 3 -z 0.15 -t "F" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.1'_h_RCM  -l -3. -m 3 -z 0.15 -t "F" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.1'_ERAINT  -l -3. -m 2 -z 0.15 -t "F" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.1'_m_RCM  -l -3. -m 2 -z 0.15 -t "F" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.1'_h_RCM  -l -3. -m 2 -z 0.15 -t "F" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.1'_m_medRCM  -l -3. -m 3 -z 0.15 -t "F" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.1'_m_medRCM  -l -3. -m 2 -z 0.15 -t "F" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.1'_h_medRCM  -l -3. -m 3 -z 0.15 -t "F" -c 7'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.1'_h_medRCM  -l -3. -m 2 -z 0.15 -t "F" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.1'_ERAINT  -l -3. -m 3 -z 0.15 -t "G" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.1'_m_RCM  -l -3. -m 3 -z 0.15 -t "G" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.1'_h_RCM  -l -3. -m 3 -z 0.15 -t "G" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.1'_ERAINT  -l -3. -m 2 -z 0.15 -t "G" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.1'_m_RCM  -l -3. -m 2 -z 0.15 -t "G" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.1'_h_RCM  -l -3. -m 2 -z 0.15 -t "G" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.1'_m_medRCM  -l -3. -m 3 -z 0.15 -t "G" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.1'_m_medRCM  -l -3. -m 2 -z 0.15 -t "G" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.1'_h_medRCM  -l -3. -m 3 -z 0.15 -t "G" -c 7'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.1'_h_medRCM  -l -3. -m 2 -z 0.15 -t "G" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.1'_ERAINT  -l -3. -m 3 -z 0.15 -t "N" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.1'_m_RCM  -l -3. -m 3 -z 0.15 -t "N" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.1'_h_RCM  -l -3. -m 3 -z 0.15 -t "N" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.1'_ERAINT  -l -3. -m 2 -z 0.15 -t "N" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.1'_m_RCM  -l -3. -m 2 -z 0.15 -t "N" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.1'_h_RCM  -l -3. -m 2 -z 0.15 -t "N" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.1'_m_medRCM  -l -3. -m 3 -z 0.15 -t N" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.1'_m_medRCM  -l -3. -m 2 -z 0.15 -t "N" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.1'_h_medRCM  -l -3. -m 3 -z 0.15 -t "N" -c 7'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.1'_h_medRCM  -l -3. -m 2 -z 0.15 -t "N" -c 1'


'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.1'_ERAINT  -l -3. -m 3 -z 0.15 -t "S" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.1'_m_RCM  -l -3. -m 3 -z 0.15 -t "S" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.1'_h_RCM  -l -3. -m 3 -z 0.15 -t "S" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.1'_ERAINT  -l -3. -m 2 -z 0.15 -t "S" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.1'_m_RCM  -l -3. -m 2 -z 0.15 -t "S" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.1'_h_RCM  -l -3. -m 2 -z 0.15 -t "S" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.1'_m_medRCM  -l -3. -m 3 -z 0.15 -t "S" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.1'_m_medRCM  -l -3. -m 2 -z 0.15 -t "S" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.1'_h_medRCM  -l -3. -m 3 -z 0.15 -t "S" -c 7'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.1'_h_medRCM  -l -3. -m 2 -z 0.15 -t "S" -c 1'

'set strsiz 0.22'
'drawstr -p 6 corr -c 1 -k 6 -t "Normalized Standard Deviation" Correlation'
'set strsiz 0.25'
'drawstr -p 2 -c 1 -t ('title.1')'


k=2
kmax=4
while(k<=kmax)



  'run 'panel.k
'set strsiz 0.15'
'set font 1'

* 0.11
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.k'_ERAINT  -l 1.5 -m 3 -z 0.15 -t "A" -c 1 -levs 0.5 1. 1.5 -levr 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.k'_m_RCM  -l 1.5 -m 3 -z 0.15 -t "A" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.k'_h_RCM  -l 1.5 -m 3 -z 0.15 -t "A" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.k'_ERAINT  -l 1.5 -m 2 -z 0.15 -t "A" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.k'_m_RCM  -l 1.5 -m 2 -z 0.15 -t "A" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.k'_h_RCM  -l 1.5 -m 2 -z 0.15 -t "A" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.k'_m_medRCM  -l 1.5 -m 3 -z 0.15 -t "A" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.k'_m_medRCM  -l 1.5 -m 2 -z 0.15 -t "A" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.k'_h_medRCM  -l 1.5 -m 3 -z 0.15 -t "A" -c 7'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Alps/'index.k'_h_medRCM  -l 1.5 -m 2 -z 0.15 -t "A" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.k'_ERAINT  -l 1.5 -m 3 -z 0.15 -t "E" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.k'_m_RCM  -l 1.5 -m 3 -z 0.15 -t "E" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.k'_h_RCM  -l 1.5 -m 3 -z 0.15 -t "E" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.k'_ERAINT  -l 1.5 -m 2 -z 0.15 -t "E" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.k'_m_RCM  -l 1.5 -m 2 -z 0.15 -t "E" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.k'_h_RCM  -l 1.5 -m 2 -z 0.15 -t "E" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.k'_m_medRCM  -l 1.5 -m 3 -z 0.15 -t "E" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.k'_m_medRCM  -l 1.5 -m 2 -z 0.15 -t "E" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.k'_h_medRCM  -l 1.5 -m 3 -z 0.15 -t "E" -c 7'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Spain/'index.k'_h_medRCM  -l 1.5 -m 2 -z 0.15 -t "E" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.k'_ERAINT  -l 1.5 -m 3 -z 0.15 -t "C" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.k'_m_RCM  -l 1.5 -m 3 -z 0.15 -t "C" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.k'_h_RCM  -l 1.5 -m 3 -z 0.15 -t "C" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.k'_ERAINT  -l 1.5 -m 2 -z 0.15 -t "C" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.k'_m_RCM  -l 1.5 -m 2 -z 0.15 -t "C" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.k'_h_RCM  -l 1.5 -m 2 -z 0.15 -t "C" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.k'_m_medRCM  -l 1.5 -m 3 -z 0.15 -t "C" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.k'_m_medRCM  -l 1.5 -m 2 -z 0.15 -t "C" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.k'_h_medRCM  -l 1.5 -m 3 -z 0.15 -t "C" -c 7'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Carpat/'index.k'_h_medRCM  -l 1.5 -m 2 -z 0.15 -t "C" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.k'_ERAINT  -l 1.5 -m 3 -z 0.15 -t "I" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.k'_m_RCM  -l 1.5 -m 3 -z 0.15 -t "I" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.k'_h_RCM  -l 1.5 -m 3 -z 0.15 -t "I" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.k'_ERAINT  -l 1.5 -m 2 -z 0.15 -t "I" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.k'_m_RCM  -l 1.5 -m 2 -z 0.15 -t "I" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.k'_h_RCM  -l 1.5 -m 2 -z 0.15 -t "I" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.k'_m_medRCM  -l 1.5 -m 3 -z 0.15 -t "I" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.k'_m_medRCM  -l 1.5 -m 2 -z 0.15 -t "I" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.k'_h_medRCM  -l 1.5 -m 3 -z 0.15 -t "I" -c 7'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy/'index.k'_h_medRCM  -l 1.5 -m 2 -z 0.15 -t "I" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.k'_ERAINT  -l 1.5 -m 3 -z 0.15 -t "U" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.k'_m_RCM  -l 1.5 -m 3 -z 0.15 -t "U" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.k'_h_RCM  -l 1.5 -m 3 -z 0.15 -t "U" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.k'_ERAINT  -l 1.5 -m 2 -z 0.15 -t "U" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.k'_m_RCM  -l 1.5 -m 2 -z 0.15 -t "U" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.k'_h_RCM  -l 1.5 -m 2 -z 0.15 -t "U" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.k'_m_medRCM  -l 1.5 -m 3 -z 0.15 -t "U" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.k'_m_medRCM  -l 1.5 -m 2 -z 0.15 -t "U" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.k'_h_medRCM  -l 1.5 -m 3 -z 0.15 -t "U" -c 7'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/England/'index.k'_h_medRCM  -l 1.5 -m 2 -z 0.15 -t "U" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.k'_ERAINT  -l 1.5 -m 3 -z 0.15 -t "F" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.k'_m_RCM  -l 1.5 -m 3 -z 0.15 -t "F" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.k'_h_RCM  -l 1.5 -m 3 -z 0.15 -t "F" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.k'_ERAINT  -l 1.5 -m 2 -z 0.15 -t "F" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.k'_m_RCM  -l 1.5 -m 2 -z 0.15 -t "F" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.k'_h_RCM  -l 1.5 -m 2 -z 0.15 -t "F" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.k'_m_medRCM  -l 1.5 -m 3 -z 0.15 -t "F" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.k'_m_medRCM  -l 1.5 -m 2 -z 0.15 -t "F" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.k'_h_medRCM  -l 1.5 -m 3 -z 0.15 -t "F" -c 7'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/France/'index.k'_h_medRCM  -l 1.5 -m 2 -z 0.15 -t "F" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.k'_ERAINT  -l 1.5 -m 3 -z 0.15 -t "G" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.k'_m_RCM  -l 1.5 -m 3 -z 0.15 -t "G" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.k'_h_RCM  -l 1.5 -m 3 -z 0.15 -t "G" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.k'_ERAINT  -l 1.5 -m 2 -z 0.15 -t "G" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.k'_m_RCM  -l 1.5 -m 2 -z 0.15 -t "G" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.k'_h_RCM  -l 1.5 -m 2 -z 0.15 -t "G" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.k'_m_medRCM  -l 1.5 -m 3 -z 0.15 -t "G" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.k'_m_medRCM  -l 1.5 -m 2 -z 0.15 -t "G" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.k'_h_medRCM  -l 1.5 -m 3 -z 0.15 -t "G" -c 7'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Germany/'index.k'_h_medRCM  -l 1.5 -m 2 -z 0.15 -t "G" -c 1'

'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.k'_ERAINT  -l 1.5 -m 3 -z 0.15 -t "N" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.k'_m_RCM  -l 1.5 -m 3 -z 0.15 -t "N" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.k'_h_RCM  -l 1.5 -m 3 -z 0.15 -t "N" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.k'_ERAINT  -l 1.5 -m 2 -z 0.15 -t "N" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.k'_m_RCM  -l 1.5 -m 2 -z 0.15 -t "N" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.k'_h_RCM  -l 1.5 -m 2 -z 0.15 -t "N" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.k'_m_medRCM  -l 1.5 -m 3 -z 0.15 -t N" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.k'_m_medRCM  -l 1.5 -m 2 -z 0.15 -t "N" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.k'_h_medRCM  -l 1.5 -m 3 -z 0.15 -t "N" -c 7'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Norway/'index.k'_h_medRCM  -l 1.5 -m 2 -z 0.15 -t "N" -c 1'


'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.k'_ERAINT  -l 1.5 -m 3 -z 0.15 -t "S" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.k'_m_RCM  -l 1.5 -m 3 -z 0.15 -t "S" -c 11'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.k'_h_RCM  -l 1.5 -m 3 -z 0.15 -t "S" -c 2'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.k'_ERAINT  -l 1.5 -m 2 -z 0.15 -t "S" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.k'_m_RCM  -l 1.5 -m 2 -z 0.15 -t "S" -c 1'
'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.k'_h_RCM  -l 1.5 -m 2 -z 0.15 -t "S" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.k'_m_medRCM  -l 1.5 -m 3 -z 0.15 -t "S" -c 9'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.k'_m_medRCM  -l 1.5 -m 2 -z 0.15 -t "S" -c 1'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.k'_h_medRCM  -l 1.5 -m 3 -z 0.15 -t "S" -c 7'
*'taylor.gs -i /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Sweden/'index.k'_h_medRCM  -l 1.5 -m 2 -z 0.15 -t "S" -c 1'

'set strsiz 0.22'
'drawstr -p 6 corr -c 1 -k 6 -t "Normalized Standard Deviation" Correlation'
'set strsiz 0.25'
'drawstr -p 2 -c 1 -t ('title.k')'


k=k+1
endwhile


  'set vpage off'
'set strsiz 0.085'
'set string 1 l 3'

*'draw string 0.2 4.5 Resolution'

'draw string 0.19 4.3     ERA-Interim'
'draw string 0.19 4.1     RCM44-All'
'draw string 0.19 3.9     RCM11-All'
'draw string 0.19 3.7     RCM11-MED'

'set line 1'
'draw mark 3 0.1 4.3 0.1'
'set line 11'
'draw mark 3 0.1 4.1 0.1'
'set line 2'
'draw mark 3 0.1 3.9 0.1'
'set line 7'
'draw mark 3 0.1 3.7 0.1'
'set line 1'
'draw mark 2 0.1 4.3 0.1'
'draw mark 2 0.1 4.1 0.1'
'draw mark 2 0.1 3.9 0.1'
'draw mark 2 0.1 3.7 0.1'
*'draw mark 5 0.1 6.7 0.1'
*'draw mark 3 0.1 6.5 0.1'

***************************




'draw string 0.2 3.1 Regions'

'draw string 0.29 2.9     Alps'           
'draw string 0.29 2.7     Carpathians'
'draw string 0.29 2.5     France'
'draw string 0.29 2.3     Germany'
'draw string 0.29 2.1     Italy'
'draw string 0.29 1.9     Norway'
'draw string 0.29 1.7     Spain'
'draw string 0.29 1.5     Sweden'
'draw string 0.29 1.3     UK'

'draw string 0.05 2.9     A'           
'draw string 0.05 2.7     C'
'draw string 0.05 2.5     F'
'draw string 0.05 2.3     G'
'draw string 0.05 2.1     I'
'draw string 0.05 1.9     N'
'draw string 0.05 1.7     E'
'draw string 0.05 1.5     S'
'draw string 0.05 1.3     U'


***************************


*Draw box
'set line 1 1 4'
* -> up
'draw line 0 4.7 1.2 4.7'
'set line 1 1 4'
* down, left
'draw line 0 4.7 0 1.0'
'set line 1 1 4'
* up, right
'draw line 1.2 4.7 1.2 1.0'
'set line 1 1 4'
* -> down
'draw line 0 1.0 1.2 1.0'

* separating 1
'draw line 0.25 3.3 1 3.3'

* separating 2
*'draw line 0.25 4.7 1 4.7'

'enable print Fig_8.gmf '
'print'
'disable'

'!gxps -c -i Fig_8.gmf  -o Fig_8.eps'
'!ps2pdf -DEPSCrop Fig_8.eps Fig_8.pdf'
'!convert -density 600 -units PixelsPerInch -rotate 90 Fig_8.eps Fig_8.jpg'
'!rm *.gmf'

*'legend'
