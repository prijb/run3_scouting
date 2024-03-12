nevents = 300000
njobs = 1000
eventsperjob = int(nevents/njobs)
with open('sub-tmp.sub', 'w') as file:
    for i in range(1, njobs+1):
        file.write(f"$ENV(MCOUTPUTDIR) output.root {i} {eventsperjob} $ENV(CMSSWVER) slc7_amd64_gcc10 $ENV(GENFILE) $ENV(ERA)\n")
