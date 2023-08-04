#!/env/bash

#declare -a samples=("Data" "DataB" "DataC" "DataD" "DataE" "DataF" "DataG")
declare -a samples=("Data")
year="2022"

for d in $(ls -d /ceph/cms/store/user/$USER/Run3ScoutingOutput/output*/)
do
    for s in ${samples[@]}
    do
	echo `python3 scripts/addHistosScouting.py --inDir $d --inSamples $s --year ${year} >& /dev/null &`
    done
done

