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
        'pt': events['Muon_pt'],
        'eta': events['Muon_eta'],
        'phi': events['Muon_phi'],
        'mass': 0.106*ak.ones_like(events['Muon_pt']),
        'charge': events['Muon_ch'],
        'vertexIdx': events['Muon_bestAssocSVIdx'],
        'lxy': events['SV_lxy'][events.Muon_bestAssocSVIdx],
        'vx': events['SV_x'][events.Muon_bestAssocSVIdx],
        'vy': events['SV_y'][events.Muon_bestAssocSVIdx],
        'vz': events['SV_z'][events.Muon_bestAssocSVIdx],
        'trackRelIso': events['Muon_trackRelIso'],
        'mindr': events['Muon_mindr'],
        'trkLayers': events["Muon_trkLayers"],
        'chi2Ndof': events["Muon_chi2Ndof"],
        }, with_name="PtEtaPhiMLorentzVector")

def get_particles(events):
    return ak.zip({
        'pt': events['GenPart_pt'],
        'eta': events['GenPart_eta'],
        'phi': events['GenPart_phi'],
        'mass': events['GenPart_m'],
        'charge': (-1)*abs(events['GenPart_pdgId'])/events["GenPart_pdgId"],
        'pdgId': events['GenPart_pdgId'],
        'vx': events['GenPart_vx'],
        'vy': events['GenPart_vy'],
        'vz': events['GenPart_vz'],
        'lxy': events['GenPart_lxy'],
        'motherIdx': events['GenPart_motherIndex'],
        'motherPdgId': events['GenPart_motherPdgId'],
        #'grandmotherIdx': events['GenMuon_grandmother_idx'],  # not there yet
        }, with_name="PtEtaPhiMLorentzVector")


# Look at ProcessorABC to see the expected methods and what they are supposed to do
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

        muon = get_muons(events)
        genpart = get_particles(events)

        muon_loose = muon
        muon_tight = muon[((muon.pt>3) & (abs(muon.eta)<2.4) & (muon.trackRelIso < 0.1) & (muon.mindr > 0.3) & (muon.chi2Ndof<3) & (muon.trkLayers>5))]

        #dimuon_sel = ak.num(muon.pt[((muon.trackRelIso < 0.1) & (muon.mindr>0.3))]) == 2
        dimuon_sel = (ak.num(muon_tight) == 2)  # & (ak.num((muon.trackRelIso<0.1)&(muon.mindr<0.3))==2))
        dimuon_sel_loose = (ak.num(muon) == 2)  # & (ak.num((muon.trackRelIso<0.1)&(muon.mindr<0.3))==2))
        sel_events = events[dimuon_sel_loose]
        #sel_muons  = muon[dimuon_sel]
        dimuon     = choose(muon_tight, 2)
        dimuon_tight = (((dimuon['0'].charge*dimuon['1'].charge)<0) & (dimuon['0'].vertexIdx == dimuon['1'].vertexIdx))


        mu1_idx    = ak.singletons(ak.argmax(muon_tight.pt, axis=1))
        mu2_idx    = ak.singletons(ak.argmin(muon_tight.pt, axis=1))
        #dimuon = choose(get_muons(events), 2)
        #OS_dimuon   = dimuon[((dimuon['0'].charge*dimuon['1'].charge)<0)]

        #output["dimuon"].fill(
        #    dataset=dataset,
        #    mass=ak.flatten(OS_dimuon.mass, axis=1),
        #    pt=ak.flatten(OS_dimuon.pt, axis=1),
        #    )
        #print(mu1_idx)
        #print(sel_muons.lxy)
        output["mu1_lxy_iso_mass"].fill(
            l = ak.flatten(muon_tight.lxy[mu1_idx][dimuon_sel]),
            iso = ak.flatten(muon_tight.trackRelIso[mu1_idx][dimuon_sel]),
            mass = ak.flatten(dimuon[dimuon_sel].mass[:,:1])
        )
        output["mu2_lxy_iso_mass"].fill(
            l = ak.flatten(muon_tight.lxy[mu2_idx][dimuon_sel]),
            iso = ak.flatten(muon_tight.trackRelIso[mu2_idx][dimuon_sel]),
            mass = ak.flatten(dimuon[dimuon_sel].mass[:,:1])
        )

        output["lxy"].fill(
            dataset = dataset,
            l = ak.flatten(muon_tight.lxy[:,:1]),
        )

        output["vertex"].fill(
            x = ak.flatten(sel_events.SV_x[sel_events.SV_lxy>2], axis=1),
            y = ak.flatten(sel_events.SV_y[sel_events.SV_lxy>2], axis=1)
            #x = ak.flatten(sel_events.SV_x, axis=1),
            #y = ak.flatten(sel_events.SV_y, axis=1)
        )

        output["EventCount"] = len(events)

        return output

    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':

    plot_dir = os.path.expandvars('/home/users/$USER/public_html/scouting/data/')
    finalizePlotDir(plot_dir)

    #fileset = {"Run2022D": glob.glob("/ceph/cms/store/user/evourlio/ScoutingRun3Output/looperOutput_20230704/output_Data_2022_*.root")}

    fileset = {
        "Run2022B": glob.glob("/ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/output_DataB_*.root"),
        "Run2022C": glob.glob("/ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/output_DataC_*.root"),
        "Run2022D": glob.glob("/ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/output_DataD_*.root"),
        "Run2022E": glob.glob("/ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/output_DataE_*.root"),
        "Run2022F": glob.glob("/ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/output_DataF_*.root"),
        "Run2022G": glob.glob("/ceph/cms/store/user/evourlio/Run3ScoutingOutput/looperOutput_20230722_10Percent/output_DataG_*.root"),
    }

    exe = processor.FuturesExecutor(workers=30)
    runner = processor.Runner(
                    exe,
                    #retries=3,
                    schema=BaseSchema,
                    chunksize=100000,
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

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_yscale('log')
    ax.set_xlabel(r"$L_{xy}$")
    ax.set_ylabel(r"Events")
    #plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/dimuon_lxy3.png')


    from matplotlib.colors import LogNorm
    fig, ax = plt.subplots(figsize=(8, 8))
    output['mu1_lxy_iso_mass'][{'mass':sum}].plot2d(
        ax=ax,
        norm=LogNorm(),
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$L_{xy}$")
    ax.set_ylabel(r"$Isolation$")
    fig.savefig(f'{plot_dir}/mu1_lxy_iso.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['mu1_lxy_iso_mass'][{'l':sum, 'iso':sum}].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$M$")
    ax.set_ylabel(r"$Events$")
    fig.savefig(f'{plot_dir}/dimuon_mass.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['vertex'].plot2d(
        ax=ax,
        norm=LogNorm(),
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    fig.savefig(f'{plot_dir}/vertex.png')


    fig, ax = plt.subplots(figsize=(8, 8))
    output['vertex'][-40j:40j:1j, -40j:40j:1j].plot2d(
        ax=ax,
        norm=LogNorm(),
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    fig.savefig(f'{plot_dir}/vertex_zoom.png')
