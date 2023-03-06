#!/usr/bin/env python3
import json
import pandas as pd
import uproot
import awkward as ak
import numpy as np
import copy

from utils import das_wrapper, redirectors

def get_muons(events, branch="Run3ScoutingMuons_hltScoutingMuonPacker__HLT.obj."):
    return ak.zip({
        'pt': events[f'{branch}pt_'],
        'eta': events[f'{branch}eta_'],
        'phi': events[f'{branch}phi_'],
        'mass': events[f'{branch}m_'],
        #'trackIso': events[f'{branch}trackIso'],
        }, with_name="PtEtaPhiMLorentzVector")


run2022 = [
    '/ScoutingPFRun3/Run2022A-v1/RAW',
#    '/ScoutingPFRun3/Run2022B-v1/RAW',
#    '/ScoutingPFRun3/Run2022C-v1/RAW',
#    '/ScoutingPFRun3/Run2022D-v1/RAW',
#    '/ScoutingPFRun3/Run2022E-v1/RAW',
#    '/ScoutingPFRun3/Run2022F-v1/RAW',
#    '/ScoutingPFRun3/Run2022G-v1/RAW',
]

if __name__ == '__main__':

    samples = {}
    for sample in run2022:
        print(f"Working on sample {sample}")
        res = das_wrapper(sample, query='summary')
        samples[sample] = json.loads(res[0])[0]

    samples_df = pd.DataFrame(samples).transpose()

    total_events = samples_df['num_event'].sum()

    print(f"Run 3 scouting has a total of {total_events/1e9:.2f} billion total events.")


    files = das_wrapper('/ScoutingPFRun3/Run2022A-v1/RAW', query='file')

    # FIXME this is not memory efficient.
    # rather, write a simple processor function that builds muons and reduces to histograms
    #
    first = True
    for f_in in files[:10]:
        print(f_in)
        with uproot.open(redirectors['fnal']+f_in) as f:
            events = f["Events"]
            try:
                muons = events["Run3ScoutingMuons_hltScoutingMuonPacker__HLT./Run3ScoutingMuons_hltScoutingMuonPacker__HLT.obj"].arrays()
                #f_in = uproot.open(redirectors['fnal']+files[0])
                #f_in["Events"].keys()
                #muons = f_in["Events"]["Run3ScoutingMuons_hltScoutingMuonPacker__HLT./Run3ScoutingMuons_hltScoutingMuonPacker__HLT.obj"].arrays()
                if first:
                    muons4 = get_muons(muons)
                    first = False
                else:
                    muons4 = ak.concatenate([muons4, copy.deepcopy(get_muons(muons))], axis=0)
                del muons

                muons4 = muons4[ak.num(muons4, axis=1)>0]
            except OSError:
                print(f"Couldn't open {f_in}")
