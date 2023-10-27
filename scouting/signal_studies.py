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

from utils import choose, finalizePlotDir, delta_phi, delta_r, delta_r2, match, choose

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
        'lxy_recomp': np.sqrt(events['SV_x'][events.Muon_bestAssocSVIdx]**2 + events['SV_y'][events.Muon_bestAssocSVIdx]**2),
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
        'lxy_recomp': np.sqrt(events['GenPart_vx']**2 + events['GenPart_vy']**2),
        'motherIdx': events['GenPart_motherIndex'],
        'motherPdgId': events['GenPart_motherPdgId'],
        #'grandmotherIdx': events['GenMuon_grandmother_idx'],  # not there yet
        }, with_name="PtEtaPhiMLorentzVector")

def get_SV(events):
    return ak.zip({
        'x': events["SV_x"],
        'y': events["SV_y"],
        'z': events["SV_z"],
        'lxy': events["SV_lxy"],
        'r_raw': np.sqrt(events["SV_x"]**2 + events["SV_y"]**2),
        'r': np.sqrt((events["SV_x"]-events["PV_x"])**2 + (events["SV_y"]-events["PV_y"])**2),
        'phi_raw': np.arctan(events.SV_y/events.SV_x),
        'phi': np.arctan((events.SV_y-events.PV_y)/(events.SV_x-events.PV_x)),
        #'eta':
    })

# Look at ProcessorABC to see the expected methods and what they are supposed to do
class MuonProcessor(processor.ProcessorABC):
    def __init__(self):

        #mass_axis = hist.axis.Regular(5000, 0.0, 50, name="mass", label=r"$M(\mu\mu)\ (GeV)$")
        mass_axis       = hist.axis.Regular(50, 0.0, 50, name="mass", label=r"$M(\mu\mu)\ (GeV)$")
        multi_axis      = hist.axis.Integer(0, 10, name="n", label="Multiplicity")
        multi2_axis      = hist.axis.Integer(0, 18, name="n2", label="Multiplicity")
        idx_axis        = hist.axis.Integer(-1, 10, name="idx", label="Index")
        #multi_axis      = hist.axis.Regular(6, -0.5, 5.5, name="n", label="Multiplicity")
        pt_axis         = hist.axis.Regular(20, 0.0, 50, name="pt", label=r"$p_{T}(\mu\mu)\ (GeV)$")
        pt2_axis         = hist.axis.Regular(30, 0.0, 30, name="pt2", label=r"$p_{T}(\mu\mu)\ (GeV)$")
        pt_axis_fine    = hist.axis.Regular(30, 0.0, 30, name="pt", label=r"$p_{T}\ (GeV)$")
        eta_axis        = hist.axis.Regular(50, -2.5, 2.5, name="eta", label=r"$\eta$")
        dataset_axis    = hist.axis.StrCategory([], name="dataset", label="Dataset", growth=True)
        lxy_axis        = hist.axis.Regular(20, 0, 10, name="l", label=r"$L_{xy}$")
        lxy_long_axis   = hist.axis.Regular(20, 0, 100, name="l", label=r"$L_{xy}$")
        delta_lxy_axis  = hist.axis.Regular(20, -100, 100, name='l', label=r"$\Delta L_{xy}")
        delta_axis      = hist.axis.Regular(20, -0.1, 0.1, name='delta', label=r"$\Delta")
        dphi_axis       = hist.axis.Regular(16, 0., 3.2, name='dphi', label=r"$\Delta \varphi")
        x_axis          = hist.axis.Regular(500, -115, 115, name="x", label=r"$x$")
        y_axis          = hist.axis.Regular(500, -115, 115, name="y", label=r"$y$")
        x_zoom_axis     = hist.axis.Regular(50, -30, 30, name="x", label=r"$x$")
        y_zoom_axis     = hist.axis.Regular(50, -30, 30, name="y", label=r"$y$")
        iso_axis        = hist.axis.Regular(50, 0, 1, name="iso", label=r"$L_{xy}$")
        dr_axis         = hist.axis.Regular(40, 0, 4, name="dr", label=r"$Delta R$")
        dr_fine_axis    = hist.axis.Regular(50, 0, 1, name="dr", label=r"$Delta R$")

        self.make_output = lambda: {
            # book histograms
            "point_phi": hist.Hist(
                dataset_axis,
                dphi_axis,
            ),

            "quadmuon": hist.Hist(
                dataset_axis,
                mass_axis,
                #multi_axis,
            ),

            "nmuon": hist.Hist(
                dataset_axis,
                multi_axis,
            ),
            "ndarkphoton": hist.Hist(
                dataset_axis,
                multi_axis,
            ),
            "ngenmuon": hist.Hist(
                dataset_axis,
                multi_axis,
            ),
            "n_genmuon_vs_darkphoton": hist.Hist(
                dataset_axis,
                multi2_axis,
                multi_axis,
            ),
            "ngenmuon_np": hist.Hist(
                dataset_axis,
                multi_axis,
            ),
            "dimuon": hist.Hist(
                dataset_axis,
                pt_axis,
                mass_axis,
            ),
            "dimuon_matched": hist.Hist(
                dataset_axis,
                pt_axis,
                mass_axis,
            ),
            "digenmuon": hist.Hist(
                dataset_axis,
                pt_axis,
                mass_axis,
            ),
            "digenmuon_siblings": hist.Hist(
                dataset_axis,
                pt_axis,
                pt2_axis,
            ),
            "lxy": hist.Hist(
                dataset_axis,
                lxy_axis,
                idx_axis,
                #
                #multi_axis,
            ),
            "lxy_long": hist.Hist(
                dataset_axis,
                lxy_long_axis,
            ),
            "mu1_lxy_iso_mass": hist.Hist(
                dataset_axis,
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
            "SV": hist.Hist(
                dataset_axis,
                lxy_axis,
            ),
            "nSV": hist.Hist(
                dataset_axis,
                multi_axis,
            ),
            "genmuon_dr": hist.Hist(
                dataset_axis,
                dr_axis,
            ),
            "genmuon_dr_min": hist.Hist(
                dataset_axis,
                dr_axis,
            ),
            "genmuon": hist.Hist(
                dataset_axis,
                pt_axis_fine,
                lxy_long_axis,
                eta_axis,
                dr_fine_axis,
            ),
            "genmuon_matched": hist.Hist(
                pt_axis_fine,
                lxy_long_axis,
                eta_axis,
                dr_fine_axis,
            ),
            "sv_res": hist.Hist(
                dataset_axis,
                x_zoom_axis,
                y_zoom_axis,
                delta_lxy_axis,
            ),
            "sv_res_tight": hist.Hist(
                x_zoom_axis,
                y_zoom_axis,
            ),
            "pt_res": hist.Hist(
                dataset_axis,
                delta_axis,
            ),
            "darkphoton": hist.Hist(
                dataset_axis,
                dr_axis,
            ),
            "EventCount": processor.value_accumulator(int),
            }

    def process(self, events):

        dataset = events.metadata["dataset"]
        output = self.make_output()

        genpart = get_particles(events)
        genmuon = genpart[((abs(genpart.pdgId)==13)&(abs(genpart.motherPdgId)==999999))]
        genmuon['mindr'] = ak.min(genmuon.metric_table(genmuon)[genmuon.metric_table(genmuon)>0], axis=2)
        genmuon_np = genpart[((abs(genpart.pdgId)==13)&(abs(genpart.motherPdgId)!=999999))]
        muon = get_muons(events)
        sv = get_SV(events)

        #global_sel = ((ak.count(genmuon.pt, axis=1)==2) & (ak.min(genmuon.mindr, axis=1)<0.2))
        #global_sel = ((ak.count(genmuon.pt, axis=1)>=4))
        #global_sel = ((ak.count(genmuon.pt, axis=1)>=4) & (ak.count(muon.pt, axis=1)>3))
        #global_sel = ((ak.count(genmuon.pt, axis=1)>=4) & (ak.count(muon.pt, axis=1)>3) & (ak.count(sv.x, axis=1)>2))
        #events = events[global_sel]

        sv = get_SV(events)
        muon = get_muons(events)

        genpart = get_particles(events)
        genmuon = genpart[((abs(genpart.pdgId)==13)&(abs(genpart.motherPdgId)==999999))]
        genmuon_np = genpart[((abs(genpart.pdgId)==13)&(abs(genpart.motherPdgId)!=999999))]
        genmuon['mindr'] = ak.min(genmuon.metric_table(genmuon)[genmuon.metric_table(genmuon)>0], axis=2)

        darkphoton = genpart[(abs(genpart.pdgId)==999999)]
        diA = choose(darkphoton, 2)
        diA_dr = delta_r(diA['0'], diA['1'])
        diA_dr_min = ak.flatten(diA_dr, axis=1)

        genmuon_matched = genmuon[match(genmuon, muon, deltaRCut=0.10)]  # muons can be very close
        digenmuon = choose(genmuon, 2)

        digenmuon_os = digenmuon[(digenmuon['0'].charge*digenmuon['1'].charge)<0]
        digenmuon_dr = delta_r(digenmuon['0'], digenmuon['1'])
        digenmuon_dr_min = ak.min(digenmuon_dr, axis=1)

        genmuon_siblings = digenmuon[(digenmuon['0'].motherIdx==digenmuon['1'].motherIdx)]
        genmuon_sibling_max = ak.max(ak.concatenate([ak.singletons(genmuon_siblings['0'].pt), ak.singletons(genmuon_siblings['1'].pt)], axis=2), axis=2)
        genmuon_sibling_min = ak.min(ak.concatenate([ak.singletons(genmuon_siblings['0'].pt), ak.singletons(genmuon_siblings['1'].pt)], axis=2), axis=2)

        muon['genmatch'] = muon.nearest(genmuon, threshold=0.1)
        muon['ismatched'] = match(muon, genmuon, deltaRCut=0.1)
        sv_res_x = muon[muon.ismatched].genmatch.vx - muon[muon.ismatched].vx
        sv_res_y = muon[muon.ismatched].genmatch.vy - muon[muon.ismatched].vy
        sv_res_lxy = muon[muon.ismatched].genmatch.lxy - muon[muon.ismatched].lxy
        pt_res = (muon[muon.ismatched].genmatch.pt - muon[muon.ismatched].pt)/muon[muon.ismatched].genmatch.pt

        #muon_matched = muon[match(muon, genmuon, deltaRCut=0.2)]
        #sv_res_x = muon_matched

        muon_loose = muon
        muon_tight = muon[((muon.pt>3) & (abs(muon.eta)<2.4) & (muon.trackRelIso < 0.1) & (muon.mindr > 0.3) & (muon.chi2Ndof<3) & (muon.trkLayers>5))]
        muon_tight_matched = muon_tight[muon_tight.ismatched]

        sv_res_x_tight = muon_tight[muon_tight.ismatched].genmatch.vx - muon_tight[muon_tight.ismatched].vx
        sv_res_y_tight = muon_tight[muon_tight.ismatched].genmatch.vy - muon_tight[muon_tight.ismatched].vy


        #dimuon_sel = ak.num(muon.pt[((muon.trackRelIso < 0.1) & (muon.mindr>0.3))]) == 2
        dimuon_sel = (ak.num(muon_tight) == 2)  # & (ak.num((muon.trackRelIso<0.1)&(muon.mindr<0.3))==2))
        dimuon_sel_loose = (ak.num(muon) == 2)  # & (ak.num((muon.trackRelIso<0.1)&(muon.mindr<0.3))==2))
        sel_events = events[dimuon_sel_loose]
        #sel_muons  = muon[dimuon_sel]
        dimuon     = choose(muon_tight, 2)
        dimuon_tight = dimuon[(((dimuon['0'].charge*dimuon['1'].charge)<0) & (dimuon['0'].vertexIdx == dimuon['1'].vertexIdx))]
        dimuon_matched = choose(muon_tight_matched, 2)
        dimuon_tight_goodmatch = dimuon[(((dimuon_matched['0'].charge*dimuon_matched['1'].charge)<0) & (dimuon_matched['0'].vertexIdx == dimuon_matched['1'].vertexIdx) & (dimuon_matched['0'].genmatch.motherIdx == dimuon_matched['1'].genmatch.motherIdx))]

        quad_muon = choose(muon, 4)


        point_phi = delta_phi(dimuon_tight, sv[dimuon_tight['0'].vertexIdx])

        mu1_idx    = ak.singletons(ak.argmax(muon_tight.pt, axis=1))
        mu2_idx    = ak.singletons(ak.argmin(muon_tight.pt, axis=1))
        #dimuon = choose(get_muons(events), 2)
        #OS_dimuon   = dimuon[((dimuon['0'].charge*dimuon['1'].charge)<0)]

        output["point_phi"].fill(
            dataset = dataset,
            dphi = ak.flatten(point_phi, axis=1),
        )

        output["quadmuon"].fill(
            dataset = dataset,
            mass = ak.flatten(quad_muon.mass, axis=1),
            #n = ak.count(sv.x, axis=1),
        )

        output["digenmuon_siblings"].fill(
            dataset=dataset,
            pt=ak.flatten(genmuon_sibling_max),
            pt2=ak.flatten(genmuon_sibling_min),
        )

        output["nmuon"].fill(
            dataset=dataset,
            n=ak.count(muon.pt, axis=1),
        )

        output["ndarkphoton"].fill(
            dataset=dataset,
            n=ak.count(darkphoton.pt, axis=1),
        )

        output["SV"].fill(
            dataset=dataset,
            l=ak.flatten(sv.r, axis=1),
        )

        output["nSV"].fill(
            dataset=dataset,
            n=ak.count(sv.r, axis=1),
        )

        output["ngenmuon"].fill(
            dataset=dataset,
            n=ak.count(genmuon.pt, axis=1),
        )

        output["n_genmuon_vs_darkphoton"].fill(
            dataset=dataset,
            n2=ak.count(genmuon.pt, axis=1),
            n=ak.count(darkphoton.pt, axis=1),
        )

        output["ngenmuon_np"].fill(
            dataset=dataset,
            n=ak.count(genmuon_np.pt, axis=1),
        )

        output["darkphoton"].fill(
            dataset=dataset,
            dr=diA_dr_min,
        )

        output["dimuon"].fill(
            dataset=dataset,
            mass=ak.flatten(dimuon_tight.mass, axis=1),
            pt=ak.flatten(dimuon_tight.pt, axis=1),
            )

        output["dimuon_matched"].fill(
            dataset=dataset,
            mass=ak.flatten(dimuon_tight_goodmatch.mass, axis=1),
            pt=ak.flatten(dimuon_tight_goodmatch.pt, axis=1),
            )
        #print(mu1_idx)
        #print(sel_muons.lxy)
        output["sv_res"].fill(
            dataset = dataset,
            x = ak.flatten(sv_res_x),
            y = ak.flatten(sv_res_y),
            l = ak.flatten(sv_res_lxy),
        ),
        output["pt_res"].fill(
            dataset = dataset,
            delta = ak.flatten(pt_res),
        ),
        output["sv_res_tight"].fill(
            x = ak.flatten(sv_res_x_tight),
            y = ak.flatten(sv_res_y_tight),
        ),

        output["digenmuon"].fill(
            dataset=dataset,
            mass=ak.flatten(digenmuon_os.mass),
            pt=ak.flatten(digenmuon_os.pt),
        ),

        output["mu1_lxy_iso_mass"].fill(
            dataset=dataset,
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
            l = ak.flatten(muon_tight_matched.lxy),
            idx = ak.flatten(muon_tight_matched.vertexIdx),
            #l = ak.flatten(muon_tight_matched[ak.count(genmuon.pt,axis=1)==2].lxy),
            #idx = ak.flatten(muon_tight_matched[ak.count(genmuon.pt,axis=1)==2].vertexIdx),
            #n = ak.count(genmuon.pt, axis=1),
        )

        output["lxy_long"].fill(
            dataset = dataset,
            l = ak.flatten(muon_tight_matched.lxy),
        )

        output["vertex"].fill(
            x = ak.flatten(events.SV_x, axis=1),
            y = ak.flatten(events.SV_y, axis=1),  # was sel_events, don't know why?
        )

        output["genmuon_dr"].fill(
            dataset = dataset,
            dr = ak.flatten(digenmuon_dr),
        )

        output["genmuon_dr_min"].fill(
            dataset = dataset,
            dr = digenmuon_dr_min,
        )

        output["genmuon"].fill(
            dataset = dataset,
            pt = ak.flatten(genmuon.pt),
            eta = ak.flatten(genmuon.eta),
            l = ak.flatten(genmuon.lxy),
            dr = ak.flatten(genmuon.mindr),
        ),

        output["genmuon_matched"].fill(
            pt = ak.flatten(genmuon_matched.pt),
            eta = ak.flatten(genmuon_matched.eta),
            l = ak.flatten(genmuon_matched.lxy),
            dr = ak.flatten(genmuon_matched.mindr),
        ),

        output["EventCount"] = len(events)

        return output

    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':

    plot_dir = os.path.expandvars('/home/users/$USER/public_html/scouting/ScenB1_test/')
    #plot_dir = os.path.expandvars('/home/users/$USER/public_html/scouting/ScenA/')
    #plot_dir = os.path.expandvars('/home/users/$USER/public_html/scouting/HTo2Zd/')
    finalizePlotDir(plot_dir)

    #fileset = {"Run2022D": glob.glob("/ceph/cms/store/user/evourlio/ScoutingRun3Output/looperOutput_20230704/output_Data_2022_*.root")}

    fileset = {
        #"ScenB1_mpi5p0_mA1p2_ctau5p0": glob.glob("/ceph/cms/store/user/fernance/Run3ScoutingOutput/Signal_ScenarioB1_mpi5p0_mA1p2_ctau5p0/output_*.root"),
        #"ScenB1_mpi5p0_mA1p6_ctau5p0": glob.glob("/ceph/cms/store/user/fernance/Run3ScoutingOutput/Signal_ScenarioB1_mpi5p0_mA1p6_ctau5p0/output_*.root"),
        #"ScenB1_mpi5p0_mA2p4_ctau5p0": glob.glob("/ceph/cms/store/user/fernance/Run3ScoutingOutput/Signal_ScenarioB1_mpi5p0_mA2p4_ctau5p0/output_*.root"),
        #"HTo2ZdTo2mu2x_MZd10_Epsilon1e-06": glob.glob("/ceph/cms/store/user/fernance/Run3ScoutingOutput/looperOutput_Sep-06-2023/output_Signal_HTo2ZdTo2mu2x_MZd10_Epsilon1e-06_2022*.root"),
        #"HTo2ZdTo2mu2x_MZd10_Epsilon1e-07": glob.glob("/ceph/cms/store/user/fernance/Run3ScoutingOutput/looperOutput_Sep-06-2023/output_Signal_HTo2ZdTo2mu2x_MZd10_Epsilon1e-07_2022*.root"),
        #"HTo2ZdTo2mu2x_MZd10_Epsilon5e-07": glob.glob("/ceph/cms/store/user/fernance/Run3ScoutingOutput/looperOutput_Sep-06-2023/output_Signal_HTo2ZdTo2mu2x_MZd10_Epsilon5e-07_2022*.root"),
        #"HTo2ZdTo2mu2x_MZd10_Epsilon3e-08": glob.glob("/ceph/cms/store/user/fernance/Run3ScoutingOutput/looperOutput_Sep-06-2023/output_Signal_HTo2ZdTo2mu2x_MZd10_Epsilon3e-08_2022*.root"),
        #"ScenA_mpi5p0_mA1p2_ctau5p0": glob.glob("/ceph/cms/store/user/fernance/Run3ScoutingOutput/looperOutput_Sep-22-2023/output_Signal_ScenarioA_mpi5p0_mA1p2_ctau5p0_2022_*.root"),
        #"ScenA_mpi5p0_mA1p6_ctau5p0": glob.glob("/ceph/cms/store/user/fernance/Run3ScoutingOutput/looperOutput_Sep-22-2023/output_Signal_ScenarioA_mpi5p0_mA1p6_ctau5p0_2022_*.root"),
        #"ScenA_mpi5p0_mA2p4_ctau5p0": glob.glob("/ceph/cms/store/user/fernance/Run3ScoutingOutput/looperOutput_Sep-22-2023/output_Signal_ScenarioA_mpi5p0_mA2p4_ctau5p0_2022_*.root"),
        #"ScenB1_30_9p9_4p8_ctau_26": glob.glob("/ceph/cms/store/user/dspitzba/Run3ScoutingOutput/testLooper/output_Signal_ScenB1_30_9p9_4p8_ctau_26_2022_*.root"),
        "ScenB1_30_9p9_4p8_ctau_26": glob.glob("/ceph/cms/store/user/dspitzba/Run3ScoutingOutput/testLooper/output_Signal_ScenB1_30_9p9_4p8_ctau_26_v7_2022_*.root"),
        "ScenA_20_5p0_1p2_ctau_23": glob.glob("/ceph/cms/store/user/dspitzba/Run3ScoutingOutput/testLooper/output_Signal_ScenA_20_5p0_1p2_ctau_23_2022_*.root"),
    }

    local = False
    if local:
        events = NanoEventsFactory.from_root(
                #fileset['ScenA_20_5p0_1p2_ctau_23'][0],
                fileset['ScenB1_30_9p9_4p8_ctau_26'][0],
                #fileset['ScenB1_mpi5p0_mA1p6_ctau5p0'][0],
                #"/ceph/cms/store/user/fernance/Run3ScoutingOutput/looperOutput_Sep-22-2023/output_Signal_ScenarioA_mpi5p0_mA1p2_ctau5p0_2022_0To99.root",
                #"/ceph/cms/store/user/fernance/Run3ScoutingOutput/looperOutput_Sep-22-2023/output_Signal_ScenarioB1_mpi5p0_mA1p2_ctau5p0_2022_0To99.root",
                schemaclass = BaseSchema,
                treepath='tout',
                #entry_stop = 1000,
        ).events()

        muons = get_muons(events)
        genpart = get_particles(events)
        genmuon = genpart[(abs(genpart.pdgId)==13)]
        sv = get_SV(events)

        raise NotImplementedError


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
    output['point_phi'].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_yscale('log')
    ax.set_xlabel(r"$\Delta \varphi$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/point_phi.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['quadmuon'][{"dataset": "ScenB1_30_9p9_4p8_ctau_26"}].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        #overlay='dataset',
        #density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    #ax.set_yscale('log')
    ax.set_xlabel(r"$M(4\mu)\ (GeV)$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/quadmuon_mas.png')



    fig, ax = plt.subplots(figsize=(8, 8))
    output['lxy'][{'idx':sum}].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_yscale('log')
    ax.set_xlabel(r"$L_{xy}$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/muon_lxy.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['lxy'][{'l':sum}].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_yscale('log')
    ax.set_xlabel('SV vertex index')
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/muon_vertexIdx.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['SV'].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_yscale('log')
    ax.set_xlabel('SV Lxy (inclusive)')
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/SV_lxy_inclusive.png')


    from matplotlib.colors import LogNorm
    #fig, ax = plt.subplots(figsize=(8, 8))
    #output['mu1_lxy_iso_mass'][{'mass':sum}].plot2d(
    #    ax=ax,
    #    norm=LogNorm(),
    #)
    #hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    #ax.set_xlabel(r"$L_{xy}$")
    #ax.set_ylabel(r"$Isolation$")
    #fig.savefig(f'{plot_dir}/mu1_lxy_iso.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['lxy_long'].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    #ax.set_yscale('log')
    ax.set_xlabel(r"$L_{xy}$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/muon_lxy_ext.png')


    fig, ax = plt.subplots(figsize=(8, 8))
    output['digenmuon_siblings'][{'pt2':sum}].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    #ax.set_yscale('log')
    ax.set_xlabel(r"$p_{T}$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/genmuon_siblings_pt_max.png')


    fig, ax = plt.subplots(figsize=(8, 8))
    output['digenmuon_siblings'][{'pt':sum}].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    #ax.set_yscale('log')
    ax.set_xlabel(r"$p_{T}$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/genmuon_siblings_pt_min.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['ndarkphoton'].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$N_{A'}$")
    ax.set_ylabel(r"$Events$")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/ndarkphoton.png')


    fig, ax = plt.subplots(figsize=(8, 8))
    output['nmuon'].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$N_{muon}$")
    ax.set_ylabel(r"$Events$")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/nmuon.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['nSV'].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$N_{SV}$")
    ax.set_ylabel(r"$Events$")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/nSV.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['nSV'].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$N_{SV}$")
    ax.set_ylabel(r"$Events$")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/nSV_abs.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['ngenmuon'].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$N_{muon}$")
    ax.set_ylabel(r"$Events$")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/ngenmuon.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['ngenmuon_np'].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$N_{muon}$")
    ax.set_ylabel(r"$Events$")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/ngenmuon_np.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['mu1_lxy_iso_mass'][{'l':sum, 'iso':sum}].plot1d(
        histtype="step",
        stack=False,
        ax=ax,
        overlay='dataset',
        density=True,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$M$")
    ax.set_ylabel(r"$Events$")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/dimuon_mass.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['digenmuon_siblings'][{"dataset":"ScenA_20_5p0_1p2_ctau_23"}].plot2d(
        ax=ax,
        norm=LogNorm(),
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$p_{T}\ (leading)$")
    ax.set_ylabel(r"$p_{T}\ (trailing)$")
    fig.savefig(f'{plot_dir}/genmuon_siblings_pt_vs_pt.png')

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


    fig, ax = plt.subplots(figsize=(8, 8))
    output['genmuon_dr'].plot1d(
        histtype="step",
        stack=False,
        overlay='dataset',
        ax=ax,
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_yscale('log')
    ax.set_xlabel(r"$\Delta R$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/genmuon_dr.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['genmuon_dr_min'].plot1d(
        histtype="step",
        stack=False,
        overlay='dataset',
        ax=ax,
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_yscale('log')
    ax.set_xlabel(r"$min \Delta R(\mu\mu)\ (gen)$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/genmuon_dr_min.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['digenmuon'][{'pt':sum}].plot1d(
        histtype="step",
        stack=False,
        overlay='dataset',
        ax=ax,
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    #ax.set_yscale('log')
    ax.set_xlabel(r"$Delta R$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/digenmuon_mass.png')


    fig, ax = plt.subplots(figsize=(8, 8))
    output['dimuon'][{'pt':sum}].plot1d(
        histtype="step",
        stack=False,
        overlay='dataset',
        ax=ax,
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    #ax.set_yscale('log')
    ax.set_xlabel(r"$Delta R$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/dimuon_mass_tight.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['dimuon_matched'][{'pt':sum}].plot1d(
        histtype="step",
        stack=False,
        overlay='dataset',
        ax=ax,
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    #ax.set_yscale('log')
    ax.set_xlabel(r"$Delta R$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/dimuon_mass_tight_matched.png')


    fig, ax = plt.subplots(figsize=(8, 8))
    output['darkphoton'].plot1d(
        histtype="step",
        stack=False,
        overlay='dataset',
        ax=ax,
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    #ax.set_yscale('log')
    ax.set_xlabel(r"$Delta R$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/darkphoton_dR_min.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['genmuon'][{'l':sum, 'eta':sum, 'dr':sum}].plot1d(
        histtype="step",
        stack=False,
        overlay='dataset',
        ax=ax,
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$p_{T}$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/genmuon_pt.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['genmuon'][{'pt':sum, 'eta':sum, 'dr':sum}].plot1d(
        histtype="step",
        stack=False,
        overlay='dataset',
        ax=ax,
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$L_{xy}$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/genmuon_lxy.png')

    #
    fig = plt.figure(figsize=(10, 8))
    output["genmuon_matched"][{'l':sum, 'eta':sum, 'dr':sum}].plot_ratio(
        output["genmuon"][{'l':sum, 'eta':sum, 'dr':sum, 'dataset':sum}],
        rp_uncert_draw_type="line",
        rp_uncertainty_type="efficiency",
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/muon_efficiency_pt.png')

    fig = plt.figure(figsize=(10, 8))
    output["genmuon_matched"][{'pt':sum, 'eta':sum, 'dr':sum}].plot_ratio(
        output["genmuon"][{'pt':sum, 'eta':sum, 'dr':sum, 'dataset':sum}],
        rp_uncert_draw_type="line",
        rp_uncertainty_type="efficiency",
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    fig.axes[0].set_yscale('log')
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/muon_efficiency_lxy.png')

    fig = plt.figure(figsize=(10, 8))
    output["genmuon_matched"][{'l':sum, 'pt':sum, 'dr':sum}].plot_ratio(
        output["genmuon"][{'l':sum, 'pt':sum, 'dr':sum, 'dataset':sum}],
        rp_uncert_draw_type="line",
        rp_uncertainty_type="efficiency",
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/muon_efficiency_eta.png')

    fig = plt.figure(figsize=(10, 8))
    output["genmuon_matched"][{'l':sum, 'pt':sum, 'eta':sum}].plot_ratio(
        output["genmuon"][{'l':sum, 'pt':sum, 'eta':sum, 'dataset':sum}],
        rp_uncert_draw_type="line",
        rp_uncertainty_type="efficiency",
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/muon_efficiency_dr.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['genmuon'][{'l':sum, 'eta':sum, 'dataset':sum}].plot2d(
        ax=ax,
        norm=LogNorm(),
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$p_{T}$")
    ax.set_ylabel(r"$min \Delta R$")
    fig.savefig(f'{plot_dir}/pt_vs_minDR.png')

    # efficiency for pt vs dR

    fig, ax = plt.subplots(figsize=(8, 8))
    (output['genmuon_matched'][{'l':sum, 'eta':sum}]/output['genmuon'][{'l':sum, 'eta':sum, 'dataset':sum}]).plot2d(
        ax=ax,
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$p_{T}$")
    ax.set_ylabel(r"$min \Delta R$")
    fig.savefig(f'{plot_dir}/efficiency_pt_vs_minDR.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['sv_res'][{"dataset":"ScenA_20_5p0_1p2_ctau_23", 'l':sum}].plot2d(
        ax=ax,
        norm=LogNorm(),
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    fig.savefig(f'{plot_dir}/SV_res_ScenA.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['sv_res'][{"dataset":"ScenB1_30_9p9_4p8_ctau_26", 'l':sum}].plot2d(
        ax=ax,
        norm=LogNorm(),
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    fig.savefig(f'{plot_dir}/SV_res_ScenB.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['sv_res'][{"dataset":sum, 'l':sum}].plot2d(
        ax=ax,
        norm=LogNorm(),
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    fig.savefig(f'{plot_dir}/SV_res.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['sv_res'][{'x':sum, 'y':sum}].plot1d(
        histtype="step",
        stack=False,
        overlay='dataset',
        ax=ax,
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$\Delta L_{xy}$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/SV_res_lxy.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['pt_res'].plot1d(
        histtype="step",
        stack=False,
        overlay='dataset',
        ax=ax,
        density=True,
    )

    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$\Delta p_{T}$")
    ax.set_ylabel(r"Events")
    plt.legend(loc=0)
    fig.savefig(f'{plot_dir}/pt_res.png')

    fig, ax = plt.subplots(figsize=(8, 8))
    output['n_genmuon_vs_darkphoton'][{"dataset":"ScenB1_30_9p9_4p8_ctau_26"}].plot2d(
        ax=ax,
        norm=LogNorm(),
    )
    hep.cms.label("Preliminary",data=True,lumi='X',com=13.6,loc=0,ax=ax,fontsize=15,)
    ax.set_xlabel(r"$N_{\mu}$")
    ax.set_ylabel(r"$N_{A'}$")
    fig.savefig(f'{plot_dir}/n_genmuon_vs_darkphoton.png')


    #eff = output['genmuon_matched'][{'l':sum, 'eta':sum}].values()/output['genmuon'][{'l':sum, 'eta':sum}].values()
    #fig, ax = plt.subplots(1,1,figsize=(10,10))

    #im = ax.matshow(eff)
    #ax.set_ylabel(r'$p_{T}$')
    #ax.set_xlabel(r'$min \Delta R$')
    ##for i in range(eff.shape[0]):
    ##    for j in range(eff.shape[1]):
    ##        c = eff[i,j]
    ##        ax.text(i, j, "%.2f"%c, va='center', ha='center')
    #cbar = ax.figure.colorbar(im)
    #cbar.ax.tick_params(labelsize=12)
    #fig.savefig(f'{plot_dir}/efficiency_pt_vs_minDR.png')
