#!/bin/bash

source /opt/ilcsoft/muonc/init_ilcsoft.sh 
export MARLIN_DLL=$MARLIN_DLL:/code/MySimpleVertex/build/libMySimpleVertex.so
cd /tmp
mkdir reco_${2}_${1}
cd reco_${2}_${1}

cp /code/Reco/acts_CONDOR.xml reco_${2}_${1}.xml
sed 's/INFILENAME/'${1}'/g' reco_${2}_${1}.xml > reco_${2}_${1}_temp.xml
mv  reco_${2}_${1}_temp.xml reco_${2}_${1}.xml
sed 's/TYPEEVENT/'${2}'/g' reco_${2}_${1}.xml > reco_${2}_${1}_temp.xml
mv  reco_${2}_${1}_temp.xml reco_${2}_${1}.xml

Marlin reco_${2}_${1}.xml

cd ..
rm -rf reco_${2}_${1}
#rm reco_EVENT_${1}.xml
#rm lctuple_${2}_actsseededckf_${1}.root
