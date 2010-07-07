#$/bin/bash

TIME="/usr/bin/time -p"
MPI="mpirun"
REGCM="./bin/regcmMPI regcm.in"

YEAR=1990
MSTEP=
RUNS=3

# monthly loop
for i in `seq 1 12`; do
	cp regcm.ok regcm.in

	if [ $i -eq 1 ]; then
		sed -i s/!RST/"false"/g regcm.in
	else
		sed -i s/!RST/"true"/g regcm.in
	fi
	
	sed -i s/!YRS/$YEAR/g regcm.in
	sed -i s/!MS/`printf "%02d" $i`/g regcm.in

	if [ $i -lt 12 ]; then
		sed -i s/!YRE/$YEAR/g regcm.in
		let ii=$((i+1))
		sed -i s/!ME/`printf "%02d" $ii`/g regcm.in
	else
		let YY=$((YEAR+1))
		sed -i s/!YRE/$YY/g regcm.in
		sed -i s/!ME/01/g regcm.in
	fi

	#cp regcm.in test$i

	# runs loop
	for j in `seq 1 $RUNS`; do
		$TIME $MPI $REGCM &> y${YEAR}_m`printf "%02d" $i`_run$j.out
	done

done
