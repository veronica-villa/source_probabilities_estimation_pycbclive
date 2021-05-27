#!/bin/bash

for TRIGGER_FILE in FINAL_COMBINE_STATMAP/*.hdf;
do
    GPS_TIME="$(echo ${TRIGGER_FILE} | cut -d'-' -f3)"
    for bankfile in BANK/*.hdf;
    do
        if [[ $bankfile = *${GPS_TIME}* ]];
        then
            BANK_FILE=${bankfile}
         fi
    done
    for sngltrigger in TRIGGER_MERGE/*.hdf;
    do
        if [[ $sngltrigger = *H1*${GPS_TIME}* ]];
        then
            H1_TRIGGER_FILE=${sngltrigger}
        elif [[ $sngltrigger = *L1*${GPS_TIME}* ]];
        then
       	    L1_TRIGGER_FILE=${sngltrigger}
       	elif [[ $sngltrigger = *V1*${GPS_TIME}* ]];
        then
            V1_TRIGGER_FILE=${sngltrigger}
        fi
    done
    ../scripts/pycbc_source_probability_offline \
        --trigger-file ${TRIGGER_FILE} \
        --bank-file ${BANK_FILE} \
        --single-detector-triggers ${H1_TRIGGER_FILE} ${L1_TRIGGER_FILE} \
                                   ${V1_TRIGGER_FILE} \
        --search-tag PYCBC_AllSky \
        --ifar-threshold 1.36986E-3 \
        --src-class-mass-gap 3.0 3.0 \
        --src-class-eff-to-lum-distance 0.74899 \
        --src-class-lum-distance-to-delta -0.51557 -0.32195 \
        --src-class-lal-cosmology \
        --verbose
done
