#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore")

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

from utils import choose, finalizePlotDir

def get_muons(events):
    return ak.zip({
        'pt': events['Muon_pt'],
        'eta': events['Muon_eta'],
        'phi': events['Muon_phi'],
        'mass': events['Muon_m'],
        'charge': events["Muon_charge"],
        'trackIso': events['Muon_trackIso'],
        }, with_name="PtEtaPhiMLorentzVector")

# Look at ProcessorABC to see the expected methods and what they are supposed to do
class MuonProcessor(processor.ProcessorABC):
    def __init__(self):

        mass_axis = hist.axis.Regular(5000, 0.0, 50, name="mass", label=r"M(\mu\mu)\ (GeV)")
        pt_axis = hist.axis.Regular(20, 0.0, 200, name="pt", label=r"p_{T}(\mu\mu)\ (GeV)")
        dataset_axis = hist.axis.StrCategory([], name="dataset", label="Dataset", growth=True)

        self.make_output = lambda: {
            # book histograms
            "dimuon": hist.Hist(
                dataset_axis,
                pt_axis,
                mass_axis,
            ),
            "EventCount": processor.value_accumulator(int),
            }

    def process(self, events):

        dataset = events.metadata["dataset"]
        output = self.make_output()

        dimuon = choose(get_muons(events), 2)
        OS_dimuon   = dimuon[((dimuon['0'].charge*dimuon['1'].charge)<0)]

        output["dimuon"].fill(
            dataset=dataset,
            mass=ak.flatten(OS_dimuon.mass, axis=1),
            pt=ak.flatten(OS_dimuon.pt, axis=1),
            )

        output["EventCount"] = len(events)

        return output

    def postprocess(self, accumulator):
        return accumulator


if __name__ == '__main__':

    events = NanoEventsFactory.from_root(
        '/ceph/cms/store/user/namin/babies4mu/ScoutingCaloMuon_Run2017_2017C_RAW_4muv4/output_10.root',
        schemaclass = BaseSchema,
        treepath='Events',
        entry_stop = 1000).events()

    muon        = get_muons(events)
    dimuon      = choose(muon, 2)
    OS_dimuon   = dimuon[((dimuon['0'].charge*dimuon['1'].charge)<0)]

    mass_axis = hist.axis.Regular(20, 0.0, 50, name="mass", label=r"M(\mu\mu)\ (GeV)")
    mass_hist = hist.Hist(mass_axis)
    mass_hist.fill(ak.flatten(OS_dimuon.mass, axis=1))

    fig, ax = plt.subplots(figsize=(8, 8))
    hep.histplot([mass_hist.values()],
            mass_hist.axes[0].edges,
            histtype="step",
            stack=False,
            ax=ax,)

    hep.cms.label("Preliminary",data=False,lumi='X',com=14,loc=0,ax=ax,fontsize=15,)

    plt.legend(loc=0)

    plot_dir = '/home/users/dspitzba/public_html/scouting/'
    finalizePlotDir(plot_dir)

    fig.savefig(f'{plot_dir}/dimuon_mass.png')

    run_all = True
    if run_all:

        fileset = {
            "Run2017C": glob.glob('/ceph/cms/store/user/namin/babies4mu/ScoutingCaloMuon_Run2017_2017C_RAW_4muv4/*.root'),
            "Run2017D": glob.glob('/ceph/cms/store/user/namin/babies4mu/ScoutingCaloMuon_Run2017_2017D_RAW_4muv4/*.root'),
            "Run2018B": glob.glob('/ceph/cms/store/user/namin/babies4mu/ScoutingCaloMuon_Run2018_2018B_RAW_4muv5/*.root'),
            "Run2018C": glob.glob('/ceph/cms/store/user/namin/babies4mu/ScoutingCaloMuon_Run2018_2018C_RAW_4muv4/*.root'),
        }

        exe = processor.FuturesExecutor(workers=8)
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
                    treename="Events",
                    processor_instance=MuonProcessor()
            )

        elapsed = time.time() - tic
        print(f"Finished in {elapsed:.1f}s")
        print(f"Total events {round(output['EventCount']/1e6,2)}M")
        print(f"Events/s: {output['EventCount'] / elapsed:.0f}")

        fig, ax = plt.subplots(figsize=(8, 8))
        hep.histplot(
            [output['dimuon'][:,:,::100j][{'dataset':sum, 'pt':sum}].values()],
            output['dimuon'][:,:,::100j][{'dataset':sum, 'pt':sum}].axes[0].edges,
            histtype="step",
            stack=False,ax=ax,
        )

        hep.cms.label("Preliminary",data=False,lumi='X',com=14,loc=0,ax=ax,fontsize=15,)

        plt.legend(loc=0)

        fig.savefig(f'{plot_dir}/dimuon_mass_all.png')

        fig, ax = plt.subplots(figsize=(8, 8))
        hep.histplot(
            [output['dimuon'][:, :, 2j:4j][{'dataset':sum, 'pt':sum}].values()],
            output['dimuon'][:, :, 2j:4j][{'dataset':sum, 'pt':sum}].axes[0].edges,
            histtype="errorbar",
            stack=False,ax=ax,
        )

        hep.cms.label("Preliminary",data=False,lumi='X',com=14,loc=0,ax=ax,fontsize=15,)

        plt.legend(loc=0)

        fig.savefig(f'{plot_dir}/dimuon_mass_zoomed.png')
