#!/bin/bash

source /opt/ilcsoft/muonc/init_ilcsoft.sh 
cd /tmp

python /code/LCIOmacros/study_efficiency.py -i /data/${2}/reco/${2}_reco_${1}.slcio -o /data/${2}/histos/${2}_eff_${1}.root
python /code/LCIOmacros/study_tracks.py -i /data/${2}/reco/${2}_reco_${1}.slcio -o /data/${2}/histos/${2}_ntup_${1}.root
