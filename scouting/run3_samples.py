#!/usr/bin/env python3
import json
import pandas as pd
import uproot
import awkward as ak
import numpy as np
import copy

from utils import das_wrapper, redirectors

run2022 = [
    '/ScoutingPFRun3/Run2022A-v1/RAW',
    '/ScoutingPFRun3/Run2022B-v1/RAW',
    '/ScoutingPFRun3/Run2022C-v1/RAW',
    '/ScoutingPFRun3/Run2022D-v1/RAW',
    '/ScoutingPFRun3/Run2022E-v1/RAW',
    '/ScoutingPFRun3/Run2022F-v1/RAW',
    '/ScoutingPFRun3/Run2022G-v1/RAW',
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

