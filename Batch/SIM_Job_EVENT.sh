#!/bin/bash

source /opt/ilcsoft/muonc/init_ilcsoft.sh 
cd /tmp
mkdir sim_${2}_${1}
cd sim_${2}_${1}
#cp /data/${2}/hepmc_3TeV/tag_1_pythia8_events.hepmc ${2}_${1}_events.hepmc 

cp /code/Sim/sim_steer_DV_CONDOR.py sim_EVENT_${1}.xml
sed 's/EVENTSTOSKIP/'${1}'/g' sim_EVENT_${1}.xml > sim_EVENT_${1}_temp.xml
mv  sim_EVENT_${1}_temp.xml sim_EVENT_${1}.xml
sed 's/OUTFILENAME/"\/data\/'${2}'\/sim\/'${2}'_sim_'${1}'.slcio"/g' sim_EVENT_${1}.xml > sim_EVENT_${1}_temp.xml
mv  sim_EVENT_${1}_temp.xml sim_EVENT_${1}.xml

ddsim --steeringFile sim_EVENT_${1}.xml --inputFiles /data/${2}/hepmc_3TeV/tag_1_pythia8_events.hepmc

cd ..
rm -rf sim_${2}_${1}
#rm sim_EVENT_${1}.xml ${2}_${1}_events.hepmc 