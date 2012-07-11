#!/bin/bash

time=0
njobs=0
max_time=0
for pdir in P*_* ; do
    for dir in $pdir/$@ ; do
	if [[ -e $dir/log.txt ]] ; then
	    x1=`grep 'Time in seconds' $dir/log.txt`
	    itime=${x1##* }
	    time=`expr $time + $itime`
	    if [ $itime -gt $max_time ] ; then
		max_time=$itime
	    fi
	    njobs=`expr $njobs + 1`
	fi
    done
done

echo 'number of jobs found is' $njobs
printf "total time: %02d days %02d hours %02d minutes %02d seconds\n" "$((time/86400))" "$((time/3600%24))" "$((time/60%60))" "$((time%60))"

time=`expr $time/$njobs`
printf "average time per integration channel: %02d days %02d hours %02d minutes %02d seconds\n" "$((time/86400))" "$((time/3600%24))" "$((time/60%60))" "$((time%60))"

printf "the longest job took: %02d days %02d hours %02d minutes %02d seconds\n" "$((max_time/86400))" "$((max_time/3600%24))" "$((max_time/60%60))" "$((max_time%60))"
