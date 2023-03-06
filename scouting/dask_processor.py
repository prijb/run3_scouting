#!/usr/bin/env python3
import socket

from distributed import Client, LocalCluster
from dask.distributed import progress
from lpcjobqueue import LPCCondorCluster
import uproot

from coffea.nanoevents.methods import vector
import awkward as ak
ak.behavior.update(vector.behavior)

import hist

from scouting.utils import das_wrapper, redirectors, choose

def get_muons(events, branch="Run3ScoutingMuons_hltScoutingMuonPacker__HLT.obj."):
    from coffea.nanoevents.methods import vector
    import awkward as ak
    ak.behavior.update(vector.behavior)
    return ak.zip({
        'pt': events[f'{branch}pt_'],
        'eta': events[f'{branch}eta_'],
        'phi': events[f'{branch}phi_'],
        'mass': events[f'{branch}m_'],
        'charge': events[f'{branch}charge_']
        #'trackIso': events[f'{branch}trackIso'],
        }, with_name="PtEtaPhiMLorentzVector")

test_files = [
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/e68468b2-1e5a-4f1c-8f68-1dac20389864.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/1d967be3-0308-4000-a6d4-75997faff68c.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/632bc65a-f8c7-40d2-a172-0fbe1de8d0e1.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/4cc9d927-4d88-4d44-8c51-763ddea5f439.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/c009bea1-ff3f-4316-a61d-ef9c69c5bfd9.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/80b42ed3-9a22-419b-a36b-f833dc9dbdf1.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/482/00000/16a0f548-52e3-4638-ba20-9777d853bbcf.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/6202e446-445d-448e-9b25-9dba5356bb5d.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/abf8278d-eaf6-4c11-84b2-7fddf98ef999.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/b905b6d6-f9ac-46ef-979e-620e5a153407.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/57c46af4-c899-4b98-bd13-6f2dcd11bc1d.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/48e8802d-53ce-445d-806c-84821b2fb479.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/cd3932b2-bcfc-4dad-9be6-022daf9a6550.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/710fff40-67bd-4a41-af0f-9e6e3ee47732.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/c30d2157-fe2c-414f-b6a5-51a5925867ea.root',
]

test_files = [redirectors["fnal"] + x for x in test_files]

def get_nevents(f_in):
    with uproot.open(f_in) as f:
        events = f["Events"]
        return len(events["double_hltScoutingPFPacker_pfMetPt_HLT./double_hltScoutingPFPacker_pfMetPt_HLT.obj"].arrays())

def make_simple_hist(f_in):
    #pt_axis = hist.axis.Regular(20, 0.0, 200, name="pt", label=r"$p_{T}^{miss\ (GeV)$")
    mass_axis = hist.axis.Regular(500, 0.0, 50, name="mass", label=r"$M(\mu\mu)\ (GeV)$")
    dataset_axis = hist.axis.StrCategory([], name="dataset", label="Dataset", growth=True)
    h = hist.Hist(
        dataset_axis,
        #pt_axis,
        mass_axis,
    )
    try:
        with uproot.open(f_in, timeout=300) as f:
            events = f["Events"]
            muons = events["Run3ScoutingMuons_hltScoutingMuonPacker__HLT./Run3ScoutingMuons_hltScoutingMuonPacker__HLT.obj"].arrays()
            muons4 = get_muons(muons)
            dimuon = choose(muons4, 2)
            OS_dimuon   = dimuon[((dimuon['0'].charge*dimuon['1'].charge)<0)]
            h.fill(
                dataset=f_in,
                mass=ak.flatten(OS_dimuon.mass, axis=1),
                #pt=ak.flatten(OS_dimuon.pt, axis=1),
                #pt=events["double_hltScoutingPFPacker_pfMetPt_HLT./double_hltScoutingPFPacker_pfMetPt_HLT.obj"].arrays(),
            )
    except OSError:
        print(f"Could not open file {f_in}, skipping.")
        #raise
    return h


    
if __name__ == '__main__':

    import argparse

    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('--cluster', action='store_true', default=False, help="Run on a cluster")
    argParser.add_argument('--workers', action='store', default=4,  help="Set the number of workers for the DASK cluster")
    args = argParser.parse_args()

    local = not args.cluster
    workers = int(args.workers)
    host = socket.gethostname()

    print("Starting cluster")
    if local:
        cluster = LocalCluster(
            n_workers=2,
            threads_per_worker=1,
        )
    else:
        if host.count('ucsd'):
            raise NotImplementedError("Can't yet use a condor cluster on UAF. Please run locally.")
        cluster = LPCCondorCluster(
            transfer_input_files="scouting",
        )

    cluster.adapt(minimum=0, maximum=workers)
    client = Client(cluster)

    files = das_wrapper('/ScoutingPFRun3/Run2022C-v1/RAW', query='file', mask='site=T2_US_Caltech')  # MIT sucks, T2_US_Caltech only has the empty files copied over...
    all_files = [redirectors["fnal"] + x for x in files][15:30]
    #files = das_wrapper('/ScoutingPFRun3/Run2022C-v1/RAW', query='file', mask=' | grep file.name, file.nevents')  # can't filter only files on Caltech...
    #new_files = []
    #for f in files:
    #    new_files += das_wrapper(f, qualifier="file", query="file", mask=" | grep file.name, file.nevents")
    # adding the redirector + weed out useless empty files
    # unfortunately,
    #all_files = [redirectors["fnal"] + x.split()[0] for x in files if int(x.split()[1])>10][:15]  # only the first few are crappy

    #all_files = test_files

    print("Computing the results")
    futures = client.map(make_simple_hist, all_files)  # .reduction does not work
    # NOTE: accumulation below is potentially memory intensive
    # This is WIP
    # https://docs.dask.org/en/stable/bag.html
    # reduce? https://stackoverflow.com/questions/70563132/distributed-chained-computing-with-dask-on-a-high-failure-rate-cluster
    progress(futures)
    print("Gathering")
    print(client.gather(futures))
    results = client.gather(futures)

    print("Accumulating")
    total_hist = sum(results)
    #total_hist[{"dataset":sum, "pt":sum}].show(columns=100)
    total_hist[{"dataset":sum}].show(columns=100)
