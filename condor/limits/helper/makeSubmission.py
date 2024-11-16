import os
import ROOT

ROOT.gROOT.ProcessLine(".L cpp/helper.C+")

# Select the model
sigModel = "HTo2ZdTo2mu2x"

# Create the file
subfile = open("runLimits_%s_onCondor.sub"%(sigModel),"w")
subfile.write("executable      = $ENV(HOMEDIR)/condor/limits/condor_executable.sh\n")
subfile.write("output          = $ENV(HOMEDIR)/condor/limits/limits_logs/job.$(ClusterId).$(ProcId).out\n")
subfile.write("error           = $ENV(HOMEDIR)/condor/limits/limits_logs/job.$(ClusterId).$(ProcId).err\n")
subfile.write("log             = $ENV(HOMEDIR)/condor/limits/limits_logs/job.$(ClusterId).$(ProcId).log\n")
subfile.write("\n")
subfile.write("getenv = True\n")
subfile.write("""+JobFlavour = "workday"\n""")
subfile.write("""+DESIRED_Sites = "T2_US_UCSD"\n""")
subfile.write("""+SingularityImage = "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64"\n""")
subfile.write("transfer_input_files = $ENV(HOMEDIR)/package.tar.gz\n")
subfile.write("should_transfer_files = YES\n")
subfile.write("when_to_transfer_output = ON_EXIT\n")
subfile.write("x509userproxy=$ENV(X509_USER_PROXY)\n")
subfile.write("use_x509userproxy = True\n")
subfile.write("\n")
subfile.write("queue arguments from (\n")
subfile.write("##\n")


# Add the points
if sigModel=="HTo2ZdTo2mu2x":
    #sigCTaus = [1]
    #lastmass = 0.5
    #while (lastmass < 45.0):
    #    lastmass = 1.04*lastmass
    #    if not ROOT.passMassVeto(lastmass):
    #        continue
    #    m = lastmass
    #    for t in sigCTaus:
    #        if ((m < 1.0 and t > 10) or (m < 30.0 and t > 100)):
    #            continue
    #        line = "$ENV(SCOUTINGSNTINPUTDIRLIM) $ENV(SCOUTINGSNTOUTPUTDIRLIM) HTo2ZdTo2mu2x asymptotic $ENV(PERIOD) %.3f %i\n"%(m, t)
    #        subfile.write(line)
    with open('data/HZdZd_limitgrid.txt', 'r') as f:
        masses = f.readlines()[0].split(',')[:-1]
        sigCTaus = [1, 10, 100]
        for mass in masses:
            m = float(mass)
            for t in sigCTaus:
                if ((m < 1.0 and t > 10) or (m < 30.0 and t > 100)):
                    continue
                line = "$ENV(SCOUTINGSNTINPUTDIRLIM) $ENV(SCOUTINGSNTOUTPUTDIRLIM) HTo2ZdTo2mu2x asymptotic $ENV(PERIOD) %.3f %i\n"%(m, t)
                subfile.write(line)

subfile.write("##\n")
subfile.write(")")
