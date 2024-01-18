#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore")

import os
import time
import glob
import numpy as np
import coffea.processor as processor

from coffea.nanoevents import NanoEventsFactory, BaseSchema, NanoAODSchema
from coffea.nanoevents.methods import vector

import hist
import matplotlib.pyplot as plt
import mplhep as hep #matplotlib wrapper for easy plotting in HEP
plt.style.use(hep.style.CMS)

import awkward as ak #just lets you do lists but faster i guess
ak.behavior.update(vector.behavior)

from particle import Particle

from utils import choose, finalizePlotDir, delta_phi, delta_r, delta_r2

import coffea.processor as processor

from coffea.nanoevents import NanoEventsFactory, BaseSchema
from coffea.nanoevents.methods import vector


class MuonProcessor(processor.ProcessorABC):
    def __init__(self):

        #mass_axis = hist.axis.Regular(5000, 0.0, 50, name="mass", label=r"$M(\mu\mu)\ (GeV)$")
        mass_axis       = hist.axis.Regular(50, 0.0, 50, name="mass", label=r"$M(\mu\mu)\ (GeV)$")
        pt_axis         = hist.axis.Regular(20, 0.0, 200, name="pt", label=r"$p_{T}(\mu\mu)\ (GeV)$")
        dataset_axis    = hist.axis.StrCategory([], name="dataset", label="Dataset", growth=True)
        lxy_axis        = hist.axis.Regular(200, 0, 100, name="l", label=r"$L_{xy}$")
        x_axis          = hist.axis.Regular(500, -115, 115, name="x", label=r"$x$")
        y_axis          = hist.axis.Regular(500, -115, 115, name="y", label=r"$y$")
        iso_axis        = hist.axis.Regular(50, 0, 1, name="iso", label=r"$L_{xy}$")

        self.make_output = lambda: {
            # book histograms
            "dimuon": hist.Hist(
                dataset_axis,
                pt_axis,
                mass_axis,
            ),
            "lxy": hist.Hist(
                dataset_axis,
                lxy_axis,
            ),
            "mu1_lxy_iso_mass": hist.Hist(
                lxy_axis,
                iso_axis,
                mass_axis,
            ),
            "mu2_lxy_iso_mass": hist.Hist(
                lxy_axis,
                iso_axis,
                mass_axis,
            ),
            "vertex": hist.Hist(
                x_axis,
                y_axis,
            ),
            "EventCount": processor.value_accumulator(int),
            }

    def process(self, events):

        dataset = events.metadata["dataset"]
        output = self.make_output()

        dimuon_sel = ak.num(events.Muon) == 2
        dimuon     = choose(events.Muon, 2)
        dimuon_tight = (((dimuon['0'].charge*dimuon['1'].charge)<0) & (dimuon['0'].svIdx == dimuon['1'].svIdx) & (dimuon['0'].svIdx>=0))

        # going from 572392 --> 280 events, 0.05% efficiency
        dimuon_sel = ((ak.num(dimuon[dimuon_tight])>0) & (ak.num(events.SV)>0))

        evs = events[((ak.num(events.Muon)>1)&(ak.num(events.SV)>0)&ak.all(events.Muon.svIdx<ak.num(events.SV), axis=1))]
        mu = evs.Muon[((evs.Muon.pt>3)&(evs.Muon.mediumId==True))]

        dimuon       = choose(mu, 2)
        dimuon_tight = (((dimuon['0'].charge*dimuon['1'].charge)<0) & (dimuon['0'].svIdx == dimuon['1'].svIdx) & (dimuon['0'].svIdx>=0))

        SV = evs.SV[dimuon[dimuon_tight]['0'].svIdx]
        PV = evs.PV

        lxy = np.sqrt((SV.x-PV.x)*(SV.x-PV.x) + (SV.y-PV.y)*(SV.y-PV.y))

        output["lxy"].fill(
            dataset = dataset,
            l = ak.flatten(lxy[ak.num(dimuon[dimuon_tight])>0]),
        )

        return output

    def postprocess(self, accumulator):
        return accumulator



if __name__ == '__main__':

    redirector_ucsd = 'root://xcache-redirector.t2.ucsd.edu:2042/'
    redirector_fnal = 'root://cmsxrootd.fnal.gov/'
    plot_dir = os.path.expandvars('/home/users/$USER/public_html/scouting/data/')
    finalizePlotDir(plot_dir)

    #fileset = {"Run2022D": glob.glob("/ceph/cms/store/user/evourlio/ScoutingRun3Output/looperOutput_20230704/output_Data_2022_*.root")}


    # debugging
    if False:
        #/store/data/Run2022D/Muon/NANOAOD/ReRecoNanoAODv11-v1/2550000/5294bc85-66a5-44b0-8857-50498838beb7.root
        events = NanoEventsFactory.from_root(
            'root://cmsxrootd.fnal.gov//store/data/Run2022D/Muon/NANOAOD/ReRecoNanoAODv11-v1/2550000/5294bc85-66a5-44b0-8857-50498838beb7.root',
            #'root://cmsxrootd.fnal.gov//store/data/Run2022D/Muon/NANOAOD/PromptNanoAODv10_v2-v1/60000/5de1482d-bcf2-4007-9bdf-2dc9ff497ed4.root',  # no Muon.svIdx in nano v10
            schemaclass = NanoAODSchema,
            treepath='Events',
        ).events()

        evs = events[((ak.num(events.Muon)>1)&(ak.num(events.SV)>0)&ak.all(events.Muon.svIdx<ak.num(events.SV), axis=1))]
        mu = evs.Muon[((evs.Muon.pt>3)&(evs.Muon.mediumId==True))]

        dimuon       = choose(mu, 2)
        dimuon_tight = (((dimuon['0'].charge*dimuon['1'].charge)<0) & (dimuon['0'].svIdx == dimuon['1'].svIdx) & (dimuon['0'].svIdx>=0))

        evs.SV[dimuon[dimuon_tight]['0'].svIdx]

    else:

        test_files = [
                '/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/bab0de38-e4ee-4739-bea5-32e4c1c870d4.root',
                '/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/3b57e408-5cdc-4bdb-ab79-b86f30470962.root',
                #'/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/b65a68bd-2960-40dc-aa73-3d509bc0dbb1.root',
                #'/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/d71983f9-7bf1-45c9-b1fd-435172ce20f4.root',
                #'/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/310d3e0c-d78d-44e2-9a98-c0f6de09ea0a.root',
                #'/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/255baaa8-d8b7-4d66-84ae-54e0a527f89a.root',
                #'/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/080364cb-b0d7-48b8-8b72-ae4f5641c381.root',
                #'/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/e1ba1b7c-7b45-4124-8ca5-0ab7916d1010.root',
                #'/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/fb024b20-d5e5-488a-ac06-10ad766469f3.root',
                #'/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/cb8843f1-c3b5-482a-8ffd-47b2ab3df386.root',
                #'/store/data/Run2022D/EGamma/NANOAOD/ReRecoNanoAODv11-v1/2540000/5236be10-4dcf-446a-ab72-e6d4a249d2af.root',
        ]
        fileset = {
                'Run2022D': [redirector_fnal + f for f in test_files]
        }

        exe = processor.FuturesExecutor(workers=30)
        runner = processor.Runner(
                        exe,
                        #retries=3,
                        schema=NanoAODSchema,
                        chunksize=100000,
                        maxchunks=None,
                )


        print("Running all events with coffea processor")
        tic = time.time()

        output = runner(
                        fileset,
                        treename="Events",
                        processor_instance=MuonProcessor()
                )

        elapsed = time.time() - tic
        print(f"Finished in {elapsed:.1f}s")
        print(f"Total events {round(output['EventCount']/1e6,2)}M")
        print(f"Events/s: {output['EventCount'] / elapsed:.0f}")
