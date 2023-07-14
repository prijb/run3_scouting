#!/bin/bash

SCRAMARCH=slc7_amd64_gcc10
CMSSWVERSION=CMSSW_12_6_0_patch1

# Compile code
cd cpp
make
cd ..

# Create package directory
rm tmp_create_package/ -rf
mkdir -p tmp_create_package
cd tmp_create_package

# Prepare package directory
mkdir -p ScoutingRun3
cp ../*.py ../data ScoutingRun3/. -r # Copy relevant folders
tar -cf - --exclude=temp_data* ../cpp | tar -xf - -C ScoutingRun3/. # Copy cpp folder without the plot folders
tar -chJf package.tar.gz ScoutingRun3
mv package.tar.gz ../.
#rm -rf tmp_create_package # Leave the package be to able to easily check what is done
