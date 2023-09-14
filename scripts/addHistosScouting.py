import ROOT
import os,sys,json
import argparse
from datetime import date    
import numpy as np
import copy
sys.path.append('utils')
import plotUtils

ROOT.gROOT.SetBatch(1)

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--inDir", default=os.environ.get("PWD")+"/outputHistograms_"+today, help="Choose output directory. Default: '"+os.environ.get("PWD")+"/outputHistograms_"+today+"'")
parser.add_argument("--inSamples", default=[], nargs="+", help="Choose sample(s); for all data samples in input directory, choose 'Data'")
parser.add_argument("--generator", default=False, action="store_true", help="Plot GEN-level histograms")
parser.add_argument("--year", default="2022", help="Year to be processes. Default: 2022")
args = parser.parse_args()

samples = args.inSamples
if len(samples)<1:
    print("Please, specify list of samples to be plotted")
    exit()

indir  = args.inDir
outdir = indir

hname = "histograms"
if args.generator:
    hname = "histograms_GEN"
for s in samples:
    if os.path.isfile("%s/%s_%s_%s_all.root"%(indir,hname,s,args.year)):
        os.system("rm -f "+indir+"/"+hname+"_"+s+"_"+args.year+"_all.root")
    if not os.path.isfile("%s/%s_%s_%s_all.root"%(indir,hname,s,args.year)):
        os.system("hadd "+indir+"/"+hname+"_"+s+"_"+args.year+"_all.root $(find "+indir+" -name '"+hname+"_"+s+"*_"+args.year+"_*.root' | awk '!/all/')")
