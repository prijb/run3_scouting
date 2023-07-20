#!/usr/bin/env python3


import warnings
warnings.filterwarnings("ignore")

import os
import time
import glob
import numpy as np
import coffea.processor as processor

from coffea.nanoevents import NanoEventsFactory, BaseSchema
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

def get_muons(events):
    return ak.zip({
        'pt': events['GenMuon_pt'],
        'eta': events['GenMuon_eta'],
        'phi': events['GenMuon_phi'],
        'mass': events['GenMuon_m'],
        'charge': (-1)*abs(events['GenMuon_pdgId'])/events["GenMuon_pdgId"],
        'lxy': events['GenMuon_lxy'],
        'vx': events['GenMuon_vx'],
        'vy': events['GenMuon_vy'],
        'vz': events['GenMuon_vz'],
        'status': events['GenMuon_status'],
        'motherIdx': events['GenMuon_mother_idx'],
        'motherId': events['GenMuon_motherId'],
        'grandmotherIdx': events['GenMuon_grandmother_idx'],
        'grandmotherId': events['GenMuon_grandmotherId'],
        }, with_name="PtEtaPhiMLorentzVector")

def get_particles(events):
    return ak.zip({
        'pt': events['GenPart_pt'],
        'eta': events['GenPart_eta'],
        'phi': events['GenPart_phi'],
        'mass': events['GenPart_m'],
        'charge': (-1)*abs(events['GenPart_Id'])/events["GenPart_Id"],
        'pdgId': events['GenPart_Id'],
        'vx': events['GenPart_vx'],
        'vy': events['GenPart_vy'],
        'vz': events['GenPart_vz'],
        #'lxy': events['GenPart_lxy'],  # not there yet
        'motherIdx': events['GenPart_mother_idx'],
        #'grandmotherIdx': events['GenMuon_grandmother_idx'],  # not there yet
        }, with_name="PtEtaPhiMLorentzVector")


# Look at ProcessorABC to see the expected methods and what they are supposed to do
class MuonProcessor(processor.ProcessorABC):
    def __init__(self):

        mass_axis = hist.axis.Regular(5000, 0.0, 50, name="mass", label=r"$M(\mu\mu)\ (GeV)$")
        pt_axis = hist.axis.Regular(20, 0.0, 200, name="pt", label=r"$p_{T}(\mu\mu)\ (GeV)$")
        dataset_axis = hist.axis.StrCategory([], name="dataset", label="Dataset", growth=True)
        lxy_axis = hist.axis.Regular(50, 0, 100, name="l", label=r"$L_{xy}$")

        self.make_output = lambda: {
            # book histograms
            "dimuon": hist.Hist(
                dataset_axis,
                pt_axis,
                mass_axis,
            ),
            "lxy": hist.Hist(
                lxy_axis,
            ),
            "EventCount": processor.value_accumulator(int),
            }

    def process(self, events):

        dataset = events.metadata["dataset"]
        output = self.make_output()

        dimuon_sel = ak.num(events.Muon_pt[((events.Muon_trackRelIso < 0.1) & (events.Muon_mindr>0.3))]) == 2
        sel_events = events[dimuon_sel]
        #dimuon = choose(get_muons(events), 2)
        #OS_dimuon   = dimuon[((dimuon['0'].charge*dimuon['1'].charge)<0)]

        #output["dimuon"].fill(
        #    dataset=dataset,
        #    mass=ak.flatten(OS_dimuon.mass, axis=1),
        #    pt=ak.flatten(OS_dimuon.pt, axis=1),
        #    )
        output["lxy"].fill(
            l = ak.flatten(sel_events.SV_lxy[sel_events.Muon_bestAssocSVIdx][:,:1]),
        )

        output["EventCount"] = len(events)

        return output

    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':

    plot_dir = os.path.expandvars('/home/users/$USER/public_html/scouting/data/')
    finalizePlotDir(plot_dir)

    fileset = {"Run2022D": glob.glob("/ceph/cms/store/user/evourlio/ScoutingRun3Output/looperOutput_20230704/output_Data_2022_*.root")}

    exe = processor.FuturesExecutor(workers=30)
    runner = processor.Runner(
                    exe,
                    #retries=3,
                    schema=BaseSchema,
                    chunksize=50000,
                    maxchunks=None,
            )


    print("Running all events with coffea processor")
    tic = time.time()

    output = runner(
                    fileset,
                    treename="tout",
                    processor_instance=MuonProcessor()
            )

    elapsed = time.time() - tic
    print(f"Finished in {elapsed:.1f}s")
    print(f"Total events {round(output['EventCount']/1e6,2)}M")
    print(f"Events/s: {output['EventCount'] / elapsed:.0f}")

    lxy_hist = output["lxy"]

    #events = NanoEventsFactory.from_root(
    #    '/ceph/cms/store/user/evourlio/ScoutingRun3Output/looperOutput_20230704/output_Data_2022_0To69.root',
    #    schemaclass = BaseSchema,
    #    treepath='tout',
    #).events()

    #dimuon_sel = ak.num(events.Muon_pt[((events.Muon_trackRelIso < 0.1) & (events.Muon_mindr>0.3))]) == 2

    #sel_events = events[dimuon_sel]

    #lxy_axis = hist.axis.Regular(50, 0, 100, name="l", label=r"$N$")

    #lxy_hist = hist.Hist(lxy_axis)
    #lxy_hist.fill(
    #    ak.flatten(sel_events.SV_lxy[sel_events.Muon_bestAssocSVIdx][:,:1]))

    fig, ax = plt.subplots(figsize=(8, 8))

    lxy_hist.plot1d(
        histtype="step",
        stack=False,
        ax=ax,
    )

    hep.cms.label("Preliminary",data=False,lumi='X',com=13,loc=0,ax=ax,fontsize=15,)
    ax.set_yscale('log')
    ax.set_xlabel(r"$L_{xy}$")
    ax.set_ylabel(r"Events")
    #plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/dimuon_lxy2.png')
