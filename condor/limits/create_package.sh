#!/bin/bash

INDIR="datacards_all_Feb-21-2024"

# Start from ScoutingRun3/.
rm tmp_create_package/ -rf
mkdir -p tmp_create_package
cd tmp_create_package

mkdir -p ScoutingRun3
cp ../HiggsAnalysis/ ScoutingRun3/. -r # Copy relevant folders
tar -cf - ../combineScripts ../${INDIR} | tar -xf - -C ScoutingRun3/. # Copy cpp folder without the plot folders

tar -chJf package.tar.gz ScoutingRun3

mv package.tar.gz ../.
rm -rf tmp_create_package
