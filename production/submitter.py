'''
Submitter script to produce samples for
https://gitlab.com/simonknapen/dark_showers_tool/-/blob/master/cards/test_scenarioA.cmnd
'''
#!/usr/bin/env python3

from metis.CondorTask import CondorTask
from metis.Sample import DirectorySample,DummySample
from metis.Path import Path
from metis.StatsParser import StatsParser
import time
import os


if __name__ == '__main__':

    total_summary = {}

    reqname = 'DarkShower_ScenarioA_default'  # Scouting For Scouting

    fragment = "signal/darkshower_fragment.py"  # This defines the physics
    events_per_point = 20
    events_per_job = 10
    njobs = int(events_per_point)//events_per_job

    tag = "v0p0"
    campaign = "Run3Summer22GS"

    task = CondorTask(
        sample = DummySample(
            dataset=f"/{reqname}/{campaign}_{tag}/GENSIM",
            N=njobs,
            nevents=int(events_per_point),
        )
        output_name = "gen.root",
        executable = "executables/test_signal_production.sh",
        tarfile = "package.tar.gz",  # contains precompiled CMSSW_12_4_12 and pythia
        additional_input_files = fragment,
        #scram_arch = "slc7_amd64_gcc630",
        open_dataset = False,
        files_per_output = 1,
        #arguments = gridpack,
        condor_submit_params = {
            "sites":"T2_US_UCSD", #
            "classads": [
                ["param_nevents",events_per_job],
                ["metis_extraargs",""],
                ["JobBatchName",reqname],
                ["SingularityImage", "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/cc7:x86_64-latest"],
                ],
            "requirements_line": 'Requirements = (HAS_SINGULARITY=?=True)'  # && (HAS_CVMFS_cms_cern_ch =?= true) && {extra_requirements})'.format(extra_requirements=extra_requirements),
            },
        tag = tag,
        min_completion_fraction = 0.90,
        )

    task.process()
    total_summary[task.get_sample().get_datasetname()] = task.get_task_summary()

    StatsParser(data=total_summary, webdir="~/public_html/dump/scouting/").do()
