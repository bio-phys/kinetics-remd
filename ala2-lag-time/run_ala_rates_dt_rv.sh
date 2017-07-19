#!/bin/bash
dt_dir_l="dt1  dt10  dt100  dt1000  dt2  dt20  dt200  dt25  dt2500  dt5  dt50  dt500  dt5000"

#dt_dir_l=`ls -d dt*`

for d in $dt_dir_l;
    do
    mkdir -p $d
    cd $d
    dt=${d#*t}
    echo $dt
    #pwd
    python ../ala_rates_dt_c.py $dt "_27jan16"
    cd ..
done

# python ../ala_rates_dt_c.py 5000 "test"

