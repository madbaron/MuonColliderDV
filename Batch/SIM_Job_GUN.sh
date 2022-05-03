#!/bin/bash

source /opt/ilcsoft/muonc/init_ilcsoft.sh 
cd /tmp
mkdir sim_${2}_${1}
cd sim_${2}_${1}

cp /code/Sim/sim_steer_${2}_CONDOR.py sim_EVENT_${2}_${1}.py
sed 's/EVENTSTOSKIP/'${1}'/g' sim_EVENT_${2}_${1}.py > sim_EVENT_${2}_${1}_temp.py
mv  sim_EVENT_${2}_${1}_temp.py sim_EVENT_${2}_${1}.py
sed 's/OUTFILENAME/"\/data\/'${2}'\/sim\/'${2}'_sim_'${1}'.slcio"/g' sim_EVENT_${2}_${1}.py > sim_EVENT_${2}_${1}_temp.py
mv  sim_EVENT_${2}_${1}_temp.py sim_EVENT_${2}_${1}.py

ddsim --steeringFile sim_EVENT_${2}_${1}.py

cd ..
rm -rf sim_${2}_${1}
