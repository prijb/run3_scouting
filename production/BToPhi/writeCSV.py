import os
import csv
import argparse

### Script only enabled for ctau configuration



# Init line params
DATASET_ = "BToPhi_MPhi-{}_ctau-{}mm_TuneCP5_13p6TeV_pythia8"
FRAGMENT_ = "genFragments/Generator/Pythia/BToLLPhi/13p6TeV/BToPhi_MPhi-{}_ctau-{}mm_TuneCP5_13p6TeV-pythia8_cfi.py"

# Generate production grid for BTPhi
mphi_grid = [0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.25, 1.5, 2.0, 2.85, 3.35, 4.0, 5.0]
ctau_grid = [0.0, 0.1, 1, 10, 100]

grid = []
for m in mphi_grid:
    for t in ctau_grid:
        mlabel = (str(m)).replace('.', 'p')
        tlabel = (str(t)).replace('.', 'p')
        grid.append([mlabel, tlabel])

print("Grid to produce:")
print(grid)

# Create the CSV file
with open("request_BToPhi.csv", "w") as file_:
    writer = csv.writer(file_)
    writer.writerow(["dataset","fragment","events","generator"])
    for p,point in enumerate(grid):
        mass = point[0]
        ctau = point[1]
        dataset = DATASET_.format(mass, ctau)
        fragment = FRAGMENT_.format(mass, ctau)
        events = "300000"
        generator = "pythia8"
        gridpack = "-"
        writer.writerow([dataset,fragment,events,generator])
        #line="{},{},{},{},{}".format(dataset, fragment, events, generator, gridpack)

#csv.close()
