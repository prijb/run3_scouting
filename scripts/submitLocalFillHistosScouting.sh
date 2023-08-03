#!/env/bash

### Data
#DataB: 291987 -> 1.2 with pace=250000
#DataC: 14797731 -> 59.2 with pace=250000
#DataD: 7707588 -> 30.8 with pace=250000
#DataE: 11616645 -> 46.5 with pace=250000
#DataF: 43122993 -> 172.5 with pace=250000
#DataG: 3669414 -> 14.7 with pace=250000

pace=1000000

if [ $# -lt 1 ]
then
    echo "Please, specify data set"
    return 0
fi

# Data 2022B
if [ $1 == "DataB" ]
then
    nentries=291987
    tj=0
    while [ $(( tj * pace )) -le ${nentries} ]
    do
	echo `python3 fillHistosScouting.py --inDir /ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/ --partialUnblinding --inSample DataB --splitIndex ${tj} --splitPace ${pace} >& /dev/null &`
	tj=$(( ${tj} + 1 ))
    done
# Data 2022C
elif [ $1 == "DataC" ]
then
    nentries=14797731
    tj=0
    while [ $(( tj * pace )) -le ${nentries} ]
    do
	echo `python3 fillHistosScouting.py --inDir /ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/ --partialUnblinding --inSample DataC --splitIndex ${tj} --splitPace ${pace} >& /dev/null &`
	tj=$(( ${tj} + 1 ))
    done
# Data 2022D
elif [ $1 == "DataD" ]
then
    nentries=7707588
    tj=0
    while [ $(( tj * pace )) -le ${nentries} ]
    do
	echo `python3 fillHistosScouting.py --inDir /ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/ --partialUnblinding --inSample DataD --splitIndex ${tj} --splitPace ${pace} >& /dev/null &`
	tj=$(( ${tj} + 1 ))
    done
# Data 2022E
elif [ $1 == "DataE" ]
then
    nentries=11616645
    tj=0
    while [ $(( tj * pace )) -le ${nentries} ]
    do
	echo `python3 fillHistosScouting.py --inDir /ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/ --partialUnblinding --inSample DataE --splitIndex ${tj} --splitPace ${pace} >& /dev/null &`
	tj=$(( ${tj} + 1 ))
    done
# Data 2022F
elif [ $1 == "DataF" ]
then
    nentries=43122993
    tj=0
    while [ $(( tj * pace )) -le ${nentries} ]
    do
	echo `python3 fillHistosScouting.py --inDir /ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/ --partialUnblinding --inSample DataF --splitIndex ${tj} --splitPace ${pace} >& /dev/null &`
	tj=$(( ${tj} + 1 ))
    done
# Data 2022G
elif [ $1 == "DataG" ]
then
    nentries=3669414
    tj=0
    while [ $(( tj * pace )) -le ${nentries} ]
    do
	echo `python3 fillHistosScouting.py --inDir /ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/ --partialUnblinding --inSample DataG --splitIndex ${tj} --splitPace ${pace} >& /dev/null &`
	tj=$(( ${tj} + 1 ))
    done
else
    echo "Please, specify existing data set."
    return 0
fi
