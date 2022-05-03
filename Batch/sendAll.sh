for entry in `ls ${1} | grep submit `; do
    condor_submit ${1}/$entry
done