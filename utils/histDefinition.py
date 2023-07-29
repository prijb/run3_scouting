import ROOT

# Histogram definition
def H1(name,nbins,low,high,xtitle,ytitle,labels=[]):
    hname = ROOT.TH1D(name,"",nbins,low,high) 
    hname.GetXaxis().SetTitle(xtitle)
    hname.GetYaxis().SetTitle(ytitle)
    if len(labels)>0:
        for b in range(1,len(labels)+1):
            hname.GetXaxis().SetBinLabel(b,labels[b-1])
    return hname

def H2(name,nbinsX,lowX,highX,xtitle,nbinsY,lowY,highY,ytitle,ztitle):
    hname = ROOT.TH2D(name,"",nbinsX,lowX,highX,nbinsY,lowY,highY)
    hname.GetXaxis().SetTitle(xtitle)
    hname.GetYaxis().SetTitle(ytitle)
    hname.GetZaxis().SetTitle(ztitle)
    return hname

def hist2dDefinition(nbinsX, lowX, highX, xtitle2d, nbinsY, lowY, highY, ytitle2d, ztitle2d, variablesXY):
    # 2D histograms
    nbinsX     ["yvsx"] = 2000
    lowX       ["yvsx"] = -100 
    highX      ["yvsx"] = 100 
    nbinsY     ["yvsx"] = 2000
    lowY       ["yvsx"] = -100 
    highY      ["yvsx"] = 100 
    xtitle2d   ["yvsx"] = "x (from PV) [cm]"
    ytitle2d   ["yvsx"] = "y (from PV) [cm]"
    ztitle2d   ["yvsx"] = "Number of SVs"
    variablesXY["yvsx"] = ["t.SV_x[v]","t.SV_y[v]"]

def hist1dDefinition(nbins, low, high, xtitle, ytitle, labels, variable):
    # 1D histograms

    #Displaced vertices

    nbins   ["nsv"] = 10
    low     ["nsv"] = 0
    high    ["nsv"] = 10
    xtitle  ["nsv"] = "Number of SVs"
    ytitle  ["nsv"] = "Events"
    variable["nsv"] = "nSVs"

    nbins   ["svchi2"] = 100
    low     ["svchi2"] = 0
    high    ["svchi2"] = 10
    xtitle  ["svchi2"] = "SV #chi^{2}"
    ytitle  ["svchi2"] = "Events / 0.1"
    variable["svchi2"] = "t.SV_chi2Ndof[v]"

    nbins   ["svchi2prob"] = 100
    low     ["svchi2prob"] = 0
    high    ["svchi2prob"] = 1
    xtitle  ["svchi2prob"] = "SV #chi^{2} probability"
    ytitle  ["svchi2prob"] = "Events / 0.01"
    variable["svchi2prob"] = "t.SV_prob[v]"

    nbins   ["svxerr"] = 500
    low     ["svxerr"] = 0
    high    ["svxerr"] = 0.5
    xtitle  ["svxerr"] = "SV x error [cm]"
    ytitle  ["svxerr"] = "Events / 0.001 cm"
    variable["svxerr"] = "t.SV_xe[v]"
    
    nbins   ["svyerr"] = 500
    low     ["svyerr"] = 0
    high    ["svyerr"] = 0.5
    xtitle  ["svyerr"] = "SV y error [cm]"
    ytitle  ["svyerr"] = "Events / 0.001 cm"
    variable["svyerr"] = "t.SV_ye[v]"
    
    nbins   ["svzerr"] = 500
    low     ["svzerr"] = 0
    high    ["svzerr"] = 0.5
    xtitle  ["svzerr"] = "SV z error [cm]"
    ytitle  ["svzerr"] = "Events / 0.001 cm"
    variable["svzerr"] = "t.SV_ze[v]"
    
    nbins   ["lxy"] = 1000
    low     ["lxy"] = 0
    high    ["lxy"] = 100.0
    xtitle  ["lxy"] = "l_{xy} (from PV) [cm]"
    ytitle  ["lxy"] = "Events / 0.1 cm"
    variable["lxy"] = "lxy"
    
    nbins   ["minlxy"] = 1000
    low     ["minlxy"] = 0
    high    ["minlxy"] = 100.0
    xtitle  ["minlxy"] = "min l_{xy} (from PV) [cm]"
    ytitle  ["minlxy"] = "Events / 0.1 cm"
    variable["minlxy"] = "minlxy"
    
    nbins   ["maxlxy"] = 1000
    low     ["maxlxy"] = 0
    high    ["maxlxy"] = 100.0
    xtitle  ["maxlxy"] = "max l_{xy} (from PV) [cm]"
    ytitle  ["maxlxy"] = "Events / 0.1 cm"
    variable["maxlxy"] = "maxlxy"
    
    nbins   ["l3d"] = 1000
    low     ["l3d"] = 0
    high    ["l3d"] = 100.0
    xtitle  ["l3d"] = "l_{3D} (from PV) [cm]"
    ytitle  ["l3d"] = "Events / 0.1 cm"
    variable["l3d"] = "t.SV_l3d[v]"

    nbins   ["mindx"] = 100
    low     ["mindx"] = 0
    high    ["mindx"] = 1
    xtitle  ["mindx"] = "min D_{x} (SV_{i}, SV_{j}) [cm]"
    ytitle  ["mindx"] = "Events / 0.01 cm"    
    variable["mindx"] = "t.SV_mindx[v]"

    nbins   ["mindy"] = 100
    low     ["mindy"] = 0
    high    ["mindy"] = 1
    xtitle  ["mindy"] = "min D_{y} (SV_{i}, SV_{j}) [cm]"
    ytitle  ["mindy"] = "Events / 0.01 cm"    
    variable["mindy"] = "t.SV_mindy[v]"

    nbins   ["mindz"] = 100
    low     ["mindz"] = 0
    high    ["mindz"] = 1
    xtitle  ["mindz"] = "min D_{z} (SV_{i}, SV_{j}) [cm]"
    ytitle  ["mindz"] = "Events / 0.01 cm"    
    variable["mindz"] = "t.SV_mindz[v]"

    nbins   ["maxdx"] = 100
    low     ["maxdx"] = 0
    high    ["maxdx"] = 1
    xtitle  ["maxdx"] = "max D_{x} (SV_{i}, SV_{j}) [cm]"
    ytitle  ["maxdx"] = "Events / 0.01 cm"    
    variable["maxdx"] = "t.SV_maxdx[v]"

    nbins   ["maxdy"] = 100
    low     ["maxdy"] = 0
    high    ["maxdy"] = 1
    xtitle  ["maxdy"] = "max D_{y} (SV_{i}, SV_{j}) [cm]"
    ytitle  ["maxdy"] = "Events / 0.01 cm"    
    variable["maxdy"] = "t.SV_maxdy[v]"

    nbins   ["maxdz"] = 100
    low     ["maxdz"] = 0
    high    ["maxdz"] = 1
    xtitle  ["maxdz"] = "max D_{z} (SV_{i}, SV_{j}) [cm]"
    ytitle  ["maxdz"] = "Events / 0.01 cm"    
    variable["maxdz"] = "t.SV_maxdz[v]"

    #Muons

    nbins   ["nmusel"] = 10
    low     ["nmusel"] = 0
    high    ["nmusel"] = 10
    xtitle  ["nmusel"] = "Number of muons (p_{T}>3 GeV, |#eta|<2.4)"
    ytitle  ["nmusel"] = "Events"
    variable["nmusel"] = "nMuSel"

    nbins   ["nmusv"] = 10
    low     ["nmusv"] = 0
    high    ["nmusv"] = 10
    xtitle  ["nmusv"] = "Number of muons (p_{T}>3 GeV, |#eta|<2.4, from SV)"
    ytitle  ["nmusv"] = "Events"
    variable["nmusv"] = "nMuAss"

    nbins   ["nmuosv"] = 10
    low     ["nmuosv"] = 0
    high    ["nmuosv"] = 10
    xtitle  ["nmuosv"] = "Number of muons (p_{T}>3 GeV, |#eta|<2.4, from overlapping SV)"
    ytitle  ["nmuosv"] = "Events"
    variable["nmuosv"] = "nMuAssOverlap"

    nbins   ["mupt"] = 50
    low     ["mupt"] = 0
    high    ["mupt"] = 150
    xtitle  ["mupt"] = "Muon p_{T} [GeV]"
    ytitle  ["mupt"] = "Events / 3 GeV"
    variable["mupt"] = "t.Muon_pt[m]"

    nbins   ["mueta"] = 50
    low     ["mueta"] = -2.5
    high    ["mueta"] = 2.5
    xtitle  ["mueta"] = "Muon #eta"
    ytitle  ["mueta"] = "Events / 0.1"
    variable["mueta"] = "t.Muon_eta[m]"

    nbins   ["muphi"] = 64
    low     ["muphi"] = -3.2
    high    ["muphi"] = 3.2
    xtitle  ["muphi"] = "Muon #phi"
    ytitle  ["muphi"] = "Events / 0.1"
    variable["muphi"] = "t.Muon_phiCorr[m]"

    nbins   ["much"] = 3
    low     ["much"] = -1.5
    high    ["much"] = 1.5
    xtitle  ["much"] = "Muon charge"
    ytitle  ["much"] = "Events"
    variable["much"] = "t.Muon_ch[m]"

    nbins   ["mumindr"] = 100
    low     ["mumindr"] = 0
    high    ["mumindr"] = 5
    xtitle  ["mumindr"] = "min #DeltaR(#mu_{i}, #mu_{j})"
    ytitle  ["mumindr"] = "Events / 0.05"
    variable["mumindr"] = "t.Muon_mindr[m]"

    nbins   ["mumaxdr"] = 100
    low     ["mumaxdr"] = 0
    high    ["mumaxdr"] = 5
    xtitle  ["mumaxdr"] = "max #DeltaR(#mu_{i}, #mu_{j})"
    ytitle  ["mumaxdr"] = "Events / 0.05"
    variable["mumaxdr"] = "t.Muon_maxdr[m]"

    nbins   ["muchi2ndof"] = 100
    low     ["muchi2ndof"] = 0
    high    ["muchi2ndof"] = 10
    xtitle  ["muchi2ndof"] = "Muon #chi^{2}/ndof"
    ytitle  ["muchi2ndof"] = "Events / 0.1"
    variable["muchi2ndof"] = "t.Muon_chi2Ndof[m]"

    labels  ["mutype"] = ["G & T","G & !T","!G & T","SA & !G & !T","!SA & !G & !T"]
    nbins   ["mutype"] = len(labels["mutype"])
    low     ["mutype"] = 0
    high    ["mutype"] = len(labels["mutype"])
    xtitle  ["mutype"] = "Muon type"
    ytitle  ["mutype"] = "Events"
    variable["mutype"] = "muonType(t.Muon_isGlobal[m], t.Muon_isTracker[m], t.Muon_isStandAlone[m])"

    nbins   ["muecaliso"] = 40
    low     ["muecaliso"] = 0
    high    ["muecaliso"] = 20
    xtitle  ["muecaliso"] = "Muon ECAL isolation [GeV]"
    ytitle  ["muecaliso"] = "Events / 0.5 GeV"
    variable["muecaliso"] = "t.Muon_ecalIso[m]"

    nbins   ["muhcaliso"] = 40
    low     ["muhcaliso"] = 0
    high    ["muhcaliso"] = 20
    xtitle  ["muhcaliso"] = "Muon HCAL isolation [GeV]"
    ytitle  ["muhcaliso"] = "Events / 0.5 GeV"
    variable["muhcaliso"] = "t.Muon_hcalIso[m]"

    nbins   ["mutrkiso"] = 40
    low     ["mutrkiso"] = 0
    high    ["mutrkiso"] = 20
    xtitle  ["mutrkiso"] = "Muon track isolation [GeV]"
    ytitle  ["mutrkiso"] = "Events / 0.5 GeV"
    variable["mutrkiso"] = "t.Muon_trackIso[m]"

    nbins   ["mupfalliso"] = 40
    low     ["mupfalliso"] = 0
    high    ["mupfalliso"] = 20
    xtitle  ["mupfalliso"] = "Muon PF-all isolation (#delta#beta) [GeV]"
    ytitle  ["mupfalliso"] = "Events / 0.5 GeV"
    variable["mupfalliso"] = "t.Muon_PFIsoAll[m]"

    nbins   ["mupfchgiso"] = 40
    low     ["mupfchgiso"] = 0
    high    ["mupfchgiso"] = 20
    xtitle  ["mupfchgiso"] = "Muon PF-charged isolation (#delta#beta) [GeV]"
    ytitle  ["mupfchgiso"] = "Events / 0.5 GeV"
    variable["mupfchgiso"] = "t.Muon_PFIsoChg[m]"

    nbins   ["muecalreliso"] = 40
    low     ["muecalreliso"] = 0
    high    ["muecalreliso"] = 2
    xtitle  ["muecalreliso"] = "Muon ECAL isolation / p_{T}"
    ytitle  ["muecalreliso"] = "Events / 0.05"
    variable["muecalreliso"] = "t.Muon_ecalRelIso[m]"

    nbins   ["muhcalreliso"] = 40
    low     ["muhcalreliso"] = 0
    high    ["muhcalreliso"] = 2
    xtitle  ["muhcalreliso"] = "Muon HCAL isolation / p_{T}"
    ytitle  ["muhcalreliso"] = "Events / 0.05"
    variable["muhcalreliso"] = "t.Muon_hcalRelIso[m]"

    nbins   ["mutrkreliso"] = 40
    low     ["mutrkreliso"] = 0
    high    ["mutrkreliso"] = 2
    xtitle  ["mutrkreliso"] = "Muon track isolation / p_{T}"
    ytitle  ["mutrkreliso"] = "Events / 0.05"
    variable["mutrkreliso"] = "t.Muon_trackRelIso[m]"

    nbins   ["mupfallreliso"] = 40
    low     ["mupfallreliso"] = 0
    high    ["mupfallreliso"] = 2
    xtitle  ["mupfallreliso"] = "Muon PF-all relisolation (#delta#beta) / p_{T}"
    ytitle  ["mupfallreliso"] = "Events / 0.05"
    variable["mupfallreliso"] = "t.Muon_PFRelIsoAll[m]"

    nbins   ["mupfchgreliso"] = 40
    low     ["mupfchgreliso"] = 0
    high    ["mupfchgreliso"] = 2
    xtitle  ["mupfchgreliso"] = "Muon PF-charged isolation (#delta#beta) / p_{T}"
    ytitle  ["mupfchgreliso"] = "Events / 0.05"
    variable["mupfchgreliso"] = "t.Muon_PFRelIsoChg[m]"

    nbins   ["mumindrjet"] = 100
    low     ["mumindrjet"] = 0
    high    ["mumindrjet"] = 5
    xtitle  ["mumindrjet"] = "min #DeltaR(#mu, PF jet)"
    ytitle  ["mumindrjet"] = "Events / 0.05"
    variable["mumindrjet"] = "t.Muon_mindrJet[m]"

    nbins   ["mumindrpfc"] = 100
    low     ["mumindrpfc"] = 0
    high    ["mumindrpfc"] = 0.1
    xtitle  ["mumindrpfc"] = "min #DeltaR(#mu, PF candidate)"
    ytitle  ["mumindrpfc"] = "Events / 0.001"
    variable["mumindrpfc"] = "t.Muon_mindrPF[m]"

    nbins   ["mudxy"] = 100
    low     ["mudxy"] = 0
    high    ["mudxy"] = 5.0
    xtitle  ["mudxy"] = "Muon |d_{xy}| [cm]"
    ytitle  ["mudxy"] = "Events / 0.05 cm"
    variable["mudxy"] = "abs(t.Muon_dxy[m])"

    nbins   ["mudxysig"] = 100
    low     ["mudxysig"] = 0
    high    ["mudxysig"] = 50
    xtitle  ["mudxysig"] = "Muon |d_{xy}|/#sigma_{xy}"
    ytitle  ["mudxysig"] = "Events / 0.5"
    variable["mudxysig"] = "abs(t.Muon_dxysig[m])"

    nbins   ["mudz"] = 300
    low     ["mudz"] = 0
    high    ["mudz"] = 30.0
    xtitle  ["mudz"] = "Muon |d_{z}| [cm]"
    ytitle  ["mudz"] = "Events / 0.1 cm"
    variable["mudz"] = "abs(t.Muon_dz[m])"

    nbins   ["mudzsig"] = 100
    low     ["mudzsig"] = 0
    high    ["mudzsig"] = 50
    xtitle  ["mudzsig"] = "Muon |d_{z}|/#sigma_{z}"
    ytitle  ["mudzsig"] = "Events / 0.5"
    variable["mudzsig"] = "abs(t.Muon_dzsig[m])"

    nbins   ["musahits"] = 50
    low     ["musahits"] = 0
    high    ["musahits"] = 50
    xtitle  ["musahits"] = "Number of valid SA muon hits"
    ytitle  ["musahits"] = "Events"
    variable["musahits"] = "t.Muon_saHits[m]"

    nbins   ["musamatchedstats"] = 10
    low     ["musamatchedstats"] = 0
    high    ["musamatchedstats"] = 10
    xtitle  ["musamatchedstats"] = "Number of SA muon matched stations"
    ytitle  ["musamatchedstats"] = "Events"
    variable["musamatchedstats"] = "t.Muon_saMatchedStats[m]"

    nbins   ["muhits"] = 50
    low     ["muhits"] = 0
    high    ["muhits"] = 50
    xtitle  ["muhits"] = "Number of valid muon hits"
    ytitle  ["muhits"] = "Events"
    variable["muhits"] = "t.Muon_muHits[m]"

    nbins   ["muchambs"] = 25
    low     ["muchambs"] = 0
    high    ["muchambs"] = 25
    xtitle  ["muchambs"] = "Number of muon chambers"
    ytitle  ["muchambs"] = "Events"
    variable["muchambs"] = "t.Muon_muChambs[m]"

    nbins   ["mucscdt"] = 20
    low     ["mucscdt"] = 0
    high    ["mucscdt"] = 20
    xtitle  ["mucscdt"] = "Number of muon chambers (CSC or DT)"
    ytitle  ["mucscdt"] = "Events"
    variable["mucscdt"] = "t.Muon_muCSCDT[m]"

    nbins   ["mumatches"] = 10
    low     ["mumatches"] = 0
    high    ["mumatches"] = 10
    xtitle  ["mumatches"] = "Number of muon matches"
    ytitle  ["mumatches"] = "Events"
    variable["mumatches"] = "t.Muon_muMatch[m]"

    nbins   ["mumatchedstats"] = 10
    low     ["mumatchedstats"] = 0
    high    ["mumatchedstats"] = 10
    xtitle  ["mumatchedstats"] = "Number of muon matched stations"
    ytitle  ["mumatchedstats"] = "Events"
    variable["mumatchedstats"] = "t.Muon_muMatchedStats[m]"

    nbins   ["muexpmatchedstats"] = 10
    low     ["muexpmatchedstats"] = 0
    high    ["muexpmatchedstats"] = 10
    xtitle  ["muexpmatchedstats"] = "Number of muon expected matched stations"
    ytitle  ["muexpmatchedstats"] = "Events"
    variable["muexpmatchedstats"] = "t.Muon_muExpMatchedStats[m]"

    nbins   ["mumatchedstatsmexp"] = 20
    low     ["mumatchedstatsmexp"] = -10
    high    ["mumatchedstatsmexp"] = 10
    xtitle  ["mumatchedstatsmexp"] = "Number of muon matched stations - expected"
    ytitle  ["mumatchedstatsmexp"] = "Events"
    variable["mumatchedstatsmexp"] = "(t.Muon_muMatchedStats[m]-t.Muon_muExpMatchedStats[m])"

    nbins   ["mumatchedrpc"] = 5
    low     ["mumatchedrpc"] = 0
    high    ["mumatchedrpc"] = 5
    xtitle  ["mumatchedrpc"] = "Number of muon matched RPC layers"
    ytitle  ["mumatchedrpc"] = "Events"
    variable["mumatchedrpc"] = "t.Muon_muMatchedRPC[m]"

    nbins   ["mupixelhits"] = 10
    low     ["mupixelhits"] = 0
    high    ["mupixelhits"] = 10
    xtitle  ["mupixelhits"] = "Number of pixel hits"
    ytitle  ["mupixelhits"] = "Events"
    variable["mupixelhits"] = "t.Muon_pixHits[m]"

    nbins   ["mupixellayers"] = 10
    low     ["mupixellayers"] = 0
    high    ["mupixellayers"] = 10
    xtitle  ["mupixellayers"] = "Number of pixel layers"
    ytitle  ["mupixellayers"] = "Events"
    variable["mupixellayers"] = "t.Muon_pixLayers[m]"

    nbins   ["mustriphits"] = 30
    low     ["mustriphits"] = 0
    high    ["mustriphits"] = 30
    xtitle  ["mustriphits"] = "Number of strip"
    ytitle  ["mustriphits"] = "Events"
    variable["mustriphits"] = "t.Muon_stripHits[m]"

    nbins   ["mutrackerlayers"] = 25
    low     ["mutrackerlayers"] = 0
    high    ["mutrackerlayers"] = 25
    xtitle  ["mutrackerlayers"] = "Number of tracker layers"
    ytitle  ["mutrackerlayers"] = "Events"
    variable["mutrackerlayers"] = "t.Muon_trkLayers[m]"

    #Di-muon

    nbins   ["dimumass"] = 1500
    low     ["dimumass"] = 0
    high    ["dimumass"] = 150
    xtitle  ["dimumass"] = "m_{#mu#mu} [GeV]"
    ytitle  ["dimumass"] = "Events / 0.1 GeV"
    variable["dimumass"] = "mass"

    nbins   ["dimupt"] = 200
    low     ["dimupt"] = 0
    high    ["dimupt"] = 600
    xtitle  ["dimupt"] = "p_{T}^{#mu#mu} [GeV]"
    ytitle  ["dimupt"] = "Events / 3 GeV"
    variable["dimupt"] = "pt"

    nbins   ["dimudr"] = 100
    low     ["dimudr"] = 0
    high    ["dimudr"] = 5
    xtitle  ["dimudr"] = "#DeltaR(#mu, #mu)"
    ytitle  ["dimudr"] = "Events / 0.05"
    variable["dimudr"] = "drmm"

    nbins   ["dimudphi"] = 32
    low     ["dimudphi"] = 0
    high    ["dimudphi"] = 3.2
    xtitle  ["dimudphi"] = "|#Delta#phi(#mu, #mu)| [rad]"
    ytitle  ["dimudphi"] = "Events / 0.1 rad"
    variable["dimudphi"] = "dpmm"

    nbins   ["dimudeta"] = 50
    low     ["dimudeta"] = 0
    high    ["dimudeta"] = 5
    xtitle  ["dimudeta"] = "|#Delta#eta(#mu, #mu)|"
    ytitle  ["dimudeta"] = "Events / 0.1"
    variable["dimudeta"] = "demm"

    nbins   ["dimudetadphiratio"] = 20
    low     ["dimudetadphiratio"] = -5
    high    ["dimudetadphiratio"] = 5
    xtitle  ["dimudetadphiratio"] = "log_{10}|#Delta#eta(#mu, #mu) / #Delta#phi(#mu, #mu)|"
    ytitle  ["dimudetadphiratio"] = "Events / 0.1"
    variable["dimudetadphiratio"] = "ROOT.TMath.Log10(dedpmm)"

    nbins   ["dimu3dangle"] = 32
    low     ["dimu3dangle"] = 0
    high    ["dimu3dangle"] = 3.2
    xtitle  ["dimu3dangle"] = "|3D angle(#mu, #mu)|"
    ytitle  ["dimu3dangle"] = "Events / 0.1"
    variable["dimu3dangle"] = "a3dmm"

    nbins   ["dimudphisv"] = 32
    low     ["dimudphisv"] = 0
    high    ["dimudphisv"] = 3.2
    xtitle  ["dimudphisv"] = "|#Delta#phi(#vec{(#mu, #mu)}, #vec{SV})| [rad]"
    ytitle  ["dimudphisv"] = "Events / 0.1 rad"
    variable["dimudphisv"] = "dphisv"

    nbins   ["dimudetasv"] = 50
    low     ["dimudetasv"] = 0
    high    ["dimudetasv"] = 5
    xtitle  ["dimudetasv"] = "|#Delta#eta(#vec{(#mu, #mu)}, #vec{SV})|"
    ytitle  ["dimudetasv"] = "Events / 0.1"
    variable["dimudetasv"] = "detasv"

    nbins   ["dimudetadphisvratio"] = 20
    low     ["dimudetadphisvratio"] = -5
    high    ["dimudetadphisvratio"] = 5
    xtitle  ["dimudetadphisvratio"] = "log_{10}(|#Delta#eta(#vec{(#mu, #mu)}, #vec{SV})|/|#Delta#phi(#vec{(#mu, #mu)}, #vec{SV})|)"
    ytitle  ["dimudetadphisvratio"] = "Events / 0.1"
    variable["dimudetadphisvratio"] = "ROOT.TMath.Log10(detadphisv)"

    nbins   ["dimu3danglesv"] = 32
    low     ["dimu3danglesv"] = 0
    high    ["dimu3danglesv"] = 3.2
    xtitle  ["dimu3danglesv"] = "|3D angle(#vec{(#mu, #mu)}, #vec{SV})|"
    ytitle  ["dimu3danglesv"] = "Events / 0.1"
    variable["dimu3danglesv"] = "a3dsv"

    # Four-muon

    nbins   ["fourmumass"] = 1500
    low     ["fourmumass"] = 0
    high    ["fourmumass"] = 150
    xtitle  ["fourmumass"] = "m_{4#mu} [GeV]"
    ytitle  ["fourmumass"] = "Events / 0.1 GeV"
    variable["fourmumass"] = "mass"

    nbins   ["fourmuminmass"] = 1500
    low     ["fourmuminmass"] = 0
    high    ["fourmuminmass"] = 150
    xtitle  ["fourmuminmass"] = "m_{#mu#mu}^{max} [GeV]"
    ytitle  ["fourmuminmass"] = "Events / 0.1 GeV"
    variable["fourmuminmass"] = "minmass"

    nbins   ["fourmumaxmass"] = 1500
    low     ["fourmumaxmass"] = 0
    high    ["fourmumaxmass"] = 150
    xtitle  ["fourmumaxmass"] = "m_{#mu#mu}^{min} [GeV]"
    ytitle  ["fourmumaxmass"] = "Events / 0.1 GeV"
    variable["fourmumaxmass"] = "maxmass"

    nbins   ["fourmuavgmass"] = 1500
    low     ["fourmuavgmass"] = 0
    high    ["fourmuavgmass"] = 150
    xtitle  ["fourmuavgmass"] = "<m_{#mu#mu}> [GeV]"
    ytitle  ["fourmuavgmass"] = "Events / 0.1 GeV"
    variable["fourmuavgmass"] = "avgmass"

    nbins   ["fourmureldmass"] = 500
    low     ["fourmureldmass"] = 0
    high    ["fourmureldmass"] = 5
    xtitle  ["fourmureldmass"] = "(m_{#mu#mu}^{max} - m_{#mu#mu}^{min})/<m_{#mu#mu}>"
    ytitle  ["fourmureldmass"] = "Events / 0.01"
    variable["fourmureldmass"] = "reldmass"

    nbins   ["fourmupt"] = 200
    low     ["fourmupt"] = 0
    high    ["fourmupt"] = 600
    xtitle  ["fourmupt"] = "p_{T}^{4#mu} [GeV]"
    ytitle  ["fourmupt"] = "Events / 3 GeV"
    variable["fourmupt"] = "pt"

    nbins   ["fourmuminpt"] = 200
    low     ["fourmuminpt"] = 0
    high    ["fourmuminpt"] = 600
    xtitle  ["fourmuminpt"] = "min p_{T}^{#mu#mu} [GeV]"
    ytitle  ["fourmuminpt"] = "Events / 3 GeV"
    variable["fourmuminpt"] = "minpt"

    nbins   ["fourmumaxpt"] = 200
    low     ["fourmumaxpt"] = 0
    high    ["fourmumaxpt"] = 600
    xtitle  ["fourmumaxpt"] = "max p_{T}^{#mu#mu} [GeV]"
    ytitle  ["fourmumaxpt"] = "Events / 3 GeV"
    variable["fourmumaxpt"] = "maxpt"

    nbins   ["fourmumindr"] = 100
    low     ["fourmumindr"] = 0
    high    ["fourmumindr"] = 5
    xtitle  ["fourmumindr"] = "min #DeltaR(#mu, #mu)"
    ytitle  ["fourmumindr"] = "Events / 0.05"
    variable["fourmumindr"] = "mindrmm"

    nbins   ["fourmumaxdr"] = 100
    low     ["fourmumaxdr"] = 0
    high    ["fourmumaxdr"] = 5
    xtitle  ["fourmumaxdr"] = "max #DeltaR(#mu, #mu)"
    ytitle  ["fourmumaxdr"] = "Events / 0.05"
    variable["fourmumaxdr"] = "maxdrmm"

    nbins   ["fourmumindphi"] = 32
    low     ["fourmumindphi"] = 0
    high    ["fourmumindphi"] = 3.2
    xtitle  ["fourmumindphi"] = "min |#Delta#phi(#mu, #mu)| [rad]"
    ytitle  ["fourmumindphi"] = "Events / 0.1 rad"
    variable["fourmumindphi"] = "mindpmm"

    nbins   ["fourmumaxdphi"] = 32
    low     ["fourmumaxdphi"] = 0
    high    ["fourmumaxdphi"] = 3.2
    xtitle  ["fourmumaxdphi"] = "max |#Delta#phi(#mu, #mu)| [rad]"
    ytitle  ["fourmumaxdphi"] = "Events / 0.1 rad"
    variable["fourmumaxdphi"] = "maxdpmm"

    nbins   ["fourmumindeta"] = 50
    low     ["fourmumindeta"] = 0
    high    ["fourmumindeta"] = 5
    xtitle  ["fourmumindeta"] = "min |#Delta#eta(#mu, #mu)|"
    ytitle  ["fourmumindeta"] = "Events / 0.1"
    variable["fourmumindeta"] = "mindemm"

    nbins   ["fourmumaxdeta"] = 50
    low     ["fourmumaxdeta"] = 0
    high    ["fourmumaxdeta"] = 5
    xtitle  ["fourmumaxdeta"] = "max |#Delta#eta(#mu, #mu)|"
    ytitle  ["fourmumaxdeta"] = "Events / 0.1"
    variable["fourmumaxdeta"] = "maxdemm"

    nbins   ["fourmumindetadphiratio"] = 20
    low     ["fourmumindetadphiratio"] = -5
    high    ["fourmumindetadphiratio"] = 5
    xtitle  ["fourmumindetadphiratio"] = "log_{10}(min(|#Delta#eta(#mu, #mu) / #Delta#phi(#mu, #mu)|))"
    ytitle  ["fourmumindetadphiratio"] = "Events / 0.1"
    variable["fourmumindetadphiratio"] = "ROOT.TMath.Log10(mindedpmm)"

    nbins   ["fourmumaxdetadphiratio"] = 20
    low     ["fourmumaxdetadphiratio"] = -5
    high    ["fourmumaxdetadphiratio"] = 5
    xtitle  ["fourmumaxdetadphiratio"] = "log_{10}(max(|#Delta#eta(#mu, #mu) / #Delta#phi(#mu, #mu)|))"
    ytitle  ["fourmumaxdetadphiratio"] = "Events / 0.1"
    variable["fourmumaxdetadphiratio"] = "ROOT.TMath.Log10(maxdedpmm)"

    nbins   ["fourmumin3dangle"] = 32
    low     ["fourmumin3dangle"] = 0
    high    ["fourmumin3dangle"] = 3.2
    xtitle  ["fourmumin3dangle"] = "min |3D angle(#mu, #mu)|"
    ytitle  ["fourmumin3dangle"] = "Events / 0.1"
    variable["fourmumin3dangle"] = "mina3dmm"

    nbins   ["fourmumax3dangle"] = 32
    low     ["fourmumax3dangle"] = 0
    high    ["fourmumax3dangle"] = 3.2
    xtitle  ["fourmumax3dangle"] = "max |3D angle(#mu, #mu)|"
    ytitle  ["fourmumax3dangle"] = "Events / 0.1"
    variable["fourmumax3dangle"] = "maxa3dmm"

    nbins   ["fourmudphisv"] = 32
    low     ["fourmudphisv"] = 0
    high    ["fourmudphisv"] = 3.2
    xtitle  ["fourmudphisv"] = "|#Delta#phi(#vec{4#mu}, #vec{SV})| [rad]"
    ytitle  ["fourmudphisv"] = "Events / 0.1 rad"
    variable["fourmudphisv"] = "dphisv"

    nbins   ["fourmumindphisv"] = 32
    low     ["fourmumindphisv"] = 0
    high    ["fourmumindphisv"] = 3.2
    xtitle  ["fourmumindphisv"] = "min |#Delta#phi(#vec{(#mu, #mu)}, #vec{SV})| [rad]"
    ytitle  ["fourmumindphisv"] = "Events / 0.1 rad"
    variable["fourmumindphisv"] = "mindphisv"

    nbins   ["fourmumaxdphisv"] = 32
    low     ["fourmumaxdphisv"] = 0
    high    ["fourmumaxdphisv"] = 3.2
    xtitle  ["fourmumaxdphisv"] = "max |#Delta#phi(#vec{(#mu, #mu)}, #vec{SV})| [rad]"
    ytitle  ["fourmumaxdphisv"] = "Events / 0.1 rad"
    variable["fourmumaxdphisv"] = "maxdphisv"

    nbins   ["fourmudetasv"] = 50
    low     ["fourmudetasv"] = 0
    high    ["fourmudetasv"] = 5
    xtitle  ["fourmudetasv"] = "|#Delta#eta(#vec{4#mu}, #vec{SV})|"
    ytitle  ["fourmudetasv"] = "Events / 0.1"
    variable["fourmudetasv"] = "detasv"

    nbins   ["fourmumindetasv"] = 50
    low     ["fourmumindetasv"] = 0
    high    ["fourmumindetasv"] = 5
    xtitle  ["fourmumindetasv"] = "min |#Delta#eta(#vec{(#mu, #mu)}, #vec{SV})|"
    ytitle  ["fourmumindetasv"] = "Events / 0.1"
    variable["fourmumindetasv"] = "mindetasv"

    nbins   ["fourmumaxdetasv"] = 50
    low     ["fourmumaxdetasv"] = 0
    high    ["fourmumaxdetasv"] = 5
    xtitle  ["fourmumaxdetasv"] = "max |#Delta#eta(#vec{(#mu, #mu)}, #vec{SV})|"
    ytitle  ["fourmumaxdetasv"] = "Events / 0.1"
    variable["fourmumaxdetasv"] = "maxdetasv"

    nbins   ["fourmudetadphisvratio"] = 20
    low     ["fourmudetadphisvratio"] = -5
    high    ["fourmudetadphisvratio"] = 5
    xtitle  ["fourmudetadphisvratio"] = "log_{10}(|#Delta#eta(#vec{4#mu}, #vec{SV})|/|#Delta#phi(#vec{4#mu}, #vec{SV})|)"
    ytitle  ["fourmudetadphisvratio"] = "Events / 0.1"
    variable["fourmudetadphisvratio"] = "ROOT.TMath.Log10(detadphisv)"

    nbins   ["fourmumindetadphisvratio"] = 20
    low     ["fourmumindetadphisvratio"] = -5
    high    ["fourmumindetadphisvratio"] = 5
    xtitle  ["fourmumindetadphisvratio"] = "log_{10}(min(|#Delta#eta(#vec{(#mu, #mu)}, #vec{SV})|/|#Delta#phi(#vec{(#mu, #mu)}, #vec{SV})|))"
    ytitle  ["fourmumindetadphisvratio"] = "Events / 0.1"
    variable["fourmumindetadphisvratio"] = "ROOT.TMath.Log10(mindetadphisv)"

    nbins   ["fourmumaxdetadphisvratio"] = 20
    low     ["fourmumaxdetadphisvratio"] = -5
    high    ["fourmumaxdetadphisvratio"] = 5
    xtitle  ["fourmumaxdetadphisvratio"] = "log_{10}(max(|#Delta#eta(#vec{(#mu, #mu)}, #vec{SV})|/|#Delta#phi(#vec{(#mu, #mu)}, #vec{SV})|))"
    ytitle  ["fourmumaxdetadphisvratio"] = "Events / 0.1"
    variable["fourmumaxdetadphisvratio"] = "ROOT.TMath.Log10(maxdetadphisv)"

    nbins   ["fourmu3danglesv"] = 32
    low     ["fourmu3danglesv"] = 0
    high    ["fourmu3danglesv"] = 3.2
    xtitle  ["fourmu3danglesv"] = "|3D angle(#vec{4#mu}, #vec{SV})|"
    ytitle  ["fourmu3danglesv"] = "Events / 0.1"
    variable["fourmu3danglesv"] = "a3dsv"

    nbins   ["fourmumin3danglesv"] = 32
    low     ["fourmumin3danglesv"] = 0
    high    ["fourmumin3danglesv"] = 3.2
    xtitle  ["fourmumin3danglesv"] = "min |3D angle(#vec{(#mu, #mu)}, #vec{SV})|"
    ytitle  ["fourmumin3danglesv"] = "Events / 0.1"
    variable["fourmumin3danglesv"] = "mina3dsv"

    nbins   ["fourmumax3danglesv"] = 32
    low     ["fourmumax3danglesv"] = 0
    high    ["fourmumax3danglesv"] = 3.2
    xtitle  ["fourmumax3danglesv"] = "max |3D angle(#vec{(#mu, #mu)}, #vec{SV})|"
    ytitle  ["fourmumax3danglesv"] = "Events / 0.1"
    variable["fourmumax3danglesv"] = "mina3dsv"

# Histogram booking
### Add your histograms here

def histBooking():
    #
    histname = []
    histtype = dict()
    #
    histname2d = []
    histtype2d = dict()
    #
    # Displaced vertices
    #
    histname.append("h_nsvsel")
    histtype[histname[-1]]="nsv"
    #
    histname.append("hsvsel_chi2ndof")
    histtype[histname[-1]]="svchi2"
    #
    histname.append("hsvsel_chi2prob")
    histtype[histname[-1]]="svchi2prob"
    #
    histname.append("hsvsel_xerr")
    histtype[histname[-1]]="svxerr"
    #
    histname.append("hsvsel_yerr")
    histtype[histname[-1]]="svyerr"
    #
    histname.append("hsvsel_zerr")
    histtype[histname[-1]]="svzerr"
    #
    histname.append("hsvsel_lxy")
    histtype[histname[-1]]="lxy"
    #
    histname.append("hsvsel_l3d")
    histtype[histname[-1]]="l3d"
    #
    histname.append("hsvsel_mindx")
    histtype[histname[-1]]="mindx"
    #
    histname.append("hsvsel_mindy")
    histtype[histname[-1]]="mindy"
    #
    histname.append("hsvsel_mindz")
    histtype[histname[-1]]="mindz"
    #
    histname.append("hsvsel_maxdx")
    histtype[histname[-1]]="maxdx"
    #
    histname.append("hsvsel_maxdy")
    histtype[histname[-1]]="maxdy"
    #
    histname.append("hsvsel_maxdz")
    histtype[histname[-1]]="maxdz"
    #
    histname2d.append("hsvsel_yvsx")
    histtype2d[histname2d[-1]]="yvsx"
    ##
    histname.append("h_nsvselass")
    histtype[histname[-1]]="nsv"
    #
    histname.append("hsvselass_chi2ndof")
    histtype[histname[-1]]="svchi2"
    #
    histname.append("hsvselass_chi2prob")
    histtype[histname[-1]]="svchi2prob"
    #
    histname.append("hsvselass_xerr")
    histtype[histname[-1]]="svxerr"
    #
    histname.append("hsvselass_yerr")
    histtype[histname[-1]]="svyerr"
    #
    histname.append("hsvselass_zerr")
    histtype[histname[-1]]="svzerr"
    #
    histname.append("hsvselass_lxy")
    histtype[histname[-1]]="lxy"
    #
    histname.append("hsvselass_l3d")
    histtype[histname[-1]]="l3d"
    #
    histname.append("hsvselass_mindx")
    histtype[histname[-1]]="mindx"
    #
    histname.append("hsvselass_mindy")
    histtype[histname[-1]]="mindy"
    #
    histname.append("hsvselass_mindz")
    histtype[histname[-1]]="mindz"
    #
    histname.append("hsvselass_maxdx")
    histtype[histname[-1]]="maxdx"
    #
    histname.append("hsvselass_maxdy")
    histtype[histname[-1]]="maxdy"
    #
    histname.append("hsvselass_maxdz")
    histtype[histname[-1]]="maxdz"
    #
    histname2d.append("hsvselass_yvsx")
    histtype2d[histname2d[-1]]="yvsx"
    #
    ##
    histname.append("h_nsvselass_osv")
    histtype[histname[-1]]="nsv"
    #
    histname.append("hsvselass_osv_chi2ndof")
    histtype[histname[-1]]="svchi2"
    #
    histname.append("hsvselass_osv_chi2prob")
    histtype[histname[-1]]="svchi2prob"
    #
    histname.append("hsvselass_osv_xerr")
    histtype[histname[-1]]="svxerr"
    #
    histname.append("hsvselass_osv_yerr")
    histtype[histname[-1]]="svyerr"
    #
    histname.append("hsvselass_osv_zerr")
    histtype[histname[-1]]="svzerr"
    #
    histname.append("hsvselass_osv_lxy")
    histtype[histname[-1]]="lxy"
    #
    histname.append("hsvselass_osv_l3d")
    histtype[histname[-1]]="l3d"
    #
    histname.append("hsvselass_osv_mindx")
    histtype[histname[-1]]="mindx"
    #
    histname.append("hsvselass_osv_mindy")
    histtype[histname[-1]]="mindy"
    #
    histname.append("hsvselass_osv_mindz")
    histtype[histname[-1]]="mindz"
    #
    histname.append("hsvselass_osv_maxdx")
    histtype[histname[-1]]="maxdx"
    #
    histname.append("hsvselass_osv_maxdy")
    histtype[histname[-1]]="maxdy"
    #
    histname.append("hsvselass_osv_maxdz")
    histtype[histname[-1]]="maxdz"
    #
    histname2d.append("hsvselass_osv_yvsx")
    histtype2d[histname2d[-1]]="yvsx"
    #
    ##
    # Muons
    #
    histname.append("h_nmuonssel")
    histtype[histname[-1]]="nmusel"
    #
    histname.append("h_nmuonsass")
    histtype[histname[-1]]="nmusv"
    #
    histname.append("h_nmuonsassoverlap")
    histtype[histname[-1]]="nmuosv"
    #
    histname.append("hmuon_pt")
    histtype[histname[-1]]="mupt"
    #
    histname.append("hmuon_eta")
    histtype[histname[-1]]="mueta"
    #
    histname.append("hmuon_phi")
    histtype[histname[-1]]="muphi"
    #
    histname.append("hmuon_ch")
    histtype[histname[-1]]="much"
    #
    histname.append("hmuon_mindr")
    histtype[histname[-1]]="mumindr"
    #
    histname.append("hmuon_maxdr")
    histtype[histname[-1]]="mumaxdr"
    #
    histname.append("hmuon_chi2ndof")
    histtype[histname[-1]]="muchi2ndof"
    #
    histname.append("hmuon_type")
    histtype[histname[-1]]="mutype"
    #
    histname.append("hmuon_ecaliso")
    histtype[histname[-1]]="muecaliso"
    #
    histname.append("hmuon_ecalreliso")
    histtype[histname[-1]]="muecalreliso"
    #
    histname.append("hmuon_hcaliso")
    histtype[histname[-1]]="muhcaliso"
    #
    histname.append("hmuon_hcalreliso")
    histtype[histname[-1]]="muhcalreliso"
    #
    histname.append("hmuon_trackiso")
    histtype[histname[-1]]="mutrkiso"
    #
    histname.append("hmuon_trackreliso")
    histtype[histname[-1]]="mutrkreliso"
    #
    histname.append("hmuon_pfiso0p3all")
    histtype[histname[-1]]="mupfalliso"
    #
    histname.append("hmuon_pfreliso0p3all")
    histtype[histname[-1]]="mupfallreliso"
    #
    histname.append("hmuon_pfiso0p3chg")
    histtype[histname[-1]]="mupfchgiso"
    #
    histname.append("hmuon_pfreliso0p3chg")
    histtype[histname[-1]]="mupfchgreliso"
    #
    histname.append("hmuon_mindrjet")
    histtype[histname[-1]]="mumindrjet"
    #
    histname.append("hmuon_mindrpfc")
    histtype[histname[-1]]="mumindrpfc"
    #
    histname.append("hmuon_dxy")
    histtype[histname[-1]]="mudxy"
    #
    histname.append("hmuon_dxysig")
    histtype[histname[-1]]="mudxysig"
    #
    histname.append("hmuon_dz")
    histtype[histname[-1]]="mudz"
    #
    histname.append("hmuon_dzsig")
    histtype[histname[-1]]="mudzsig"
    #
    histname.append("hmuon_nsahits")
    histtype[histname[-1]]="musahits"
    #
    histname.append("hmuon_nsamatchedstats")
    histtype[histname[-1]]="musamatchedstats"
    #
    histname.append("hmuon_nmuhits")
    histtype[histname[-1]]="muhits"
    #
    histname.append("hmuon_nmuchambs")
    histtype[histname[-1]]="muchambs"
    #
    histname.append("hmuon_nmuCSCorDT")
    histtype[histname[-1]]="mucscdt"
    #
    histname.append("hmuon_nmumatches")
    histtype[histname[-1]]="mumatches"
    #
    histname.append("hmuon_nmumatchedstats")
    histtype[histname[-1]]="mumatchedstats"
    #
    histname.append("hmuon_nmuexpmatchedstats")
    histtype[histname[-1]]="muexpmatchedstats"
    #
    histname.append("hmuon_nmumatchedstatsmexp")
    histtype[histname[-1]]="mumatchedstatsmexp"
    #
    histname.append("hmuon_nmumatchedRPC")
    histtype[histname[-1]]="mumatchedrpc"
    #
    histname.append("hmuon_npixelhits")
    histtype[histname[-1]]="mupixelhits"
    #
    histname.append("hmuon_npixellayers")
    histtype[histname[-1]]="mupixellayers"
    #
    histname.append("hmuon_nstriphits")
    histtype[histname[-1]]="mustriphits"
    #
    histname.append("hmuon_ntrackerlayers")
    histtype[histname[-1]]="mutrackerlayers"
    #
    ##
    #
    histname.append("hselmuon_pt")
    histtype[histname[-1]]="mupt"
    #
    histname.append("hselmuon_eta")
    histtype[histname[-1]]="mueta"
    #
    histname.append("hselmuon_phi")
    histtype[histname[-1]]="muphi"
    #
    histname.append("hselmuon_ch")
    histtype[histname[-1]]="much"
    #
    histname.append("hselmuon_mindr")
    histtype[histname[-1]]="mumindr"
    #
    histname.append("hselmuon_maxdr")
    histtype[histname[-1]]="mumaxdr"
    #
    histname.append("hselmuon_chi2ndof")
    histtype[histname[-1]]="muchi2ndof"
    #
    histname.append("hselmuon_type")
    histtype[histname[-1]]="mutype"
    #
    histname.append("hselmuon_ecaliso")
    histtype[histname[-1]]="muecaliso"
    #
    histname.append("hselmuon_ecalreliso")
    histtype[histname[-1]]="muecalreliso"
    #
    histname.append("hselmuon_hcaliso")
    histtype[histname[-1]]="muhcaliso"
    #
    histname.append("hselmuon_hcalreliso")
    histtype[histname[-1]]="muhcalreliso"
    #
    histname.append("hselmuon_trackiso")
    histtype[histname[-1]]="mutrkiso"
    #
    histname.append("hselmuon_trackreliso")
    histtype[histname[-1]]="mutrkreliso"
    #
    histname.append("hselmuon_pfiso0p3all")
    histtype[histname[-1]]="mupfalliso"
    #
    histname.append("hselmuon_pfreliso0p3all")
    histtype[histname[-1]]="mupfallreliso"
    #
    histname.append("hselmuon_pfiso0p3chg")
    histtype[histname[-1]]="mupfchgiso"
    #
    histname.append("hselmuon_pfreliso0p3chg")
    histtype[histname[-1]]="mupfchgreliso"
    #
    histname.append("hselmuon_mindrjet")
    histtype[histname[-1]]="mumindrjet"
    #
    histname.append("hselmuon_mindrpfc")
    histtype[histname[-1]]="mumindrpfc"
    #
    histname.append("hselmuon_dxy")
    histtype[histname[-1]]="mudxy"
    #
    histname.append("hselmuon_dxysig")
    histtype[histname[-1]]="mudxysig"
    #
    histname.append("hselmuon_dz")
    histtype[histname[-1]]="mudz"
    #
    histname.append("hselmuon_dzsig")
    histtype[histname[-1]]="mudzsig"
    #
    histname.append("hselmuon_nsahits")
    histtype[histname[-1]]="musahits"
    #
    histname.append("hselmuon_nsamatchedstats")
    histtype[histname[-1]]="musamatchedstats"
    #
    histname.append("hselmuon_nmuhits")
    histtype[histname[-1]]="muhits"
    #
    histname.append("hselmuon_nmuchambs")
    histtype[histname[-1]]="muchambs"
    #
    histname.append("hselmuon_nmuCSCorDT")
    histtype[histname[-1]]="mucscdt"
    #
    histname.append("hselmuon_nmumatches")
    histtype[histname[-1]]="mumatches"
    #
    histname.append("hselmuon_nmumatchedstats")
    histtype[histname[-1]]="mumatchedstats"
    #
    histname.append("hselmuon_nmuexpmatchedstats")
    histtype[histname[-1]]="muexpmatchedstats"
    #
    histname.append("hselmuon_nmumatchedstatsmexp")
    histtype[histname[-1]]="mumatchedstatsmexp"
    #
    histname.append("hselmuon_nmumatchedRPC")
    histtype[histname[-1]]="mumatchedrpc"
    #
    histname.append("hselmuon_npixelhits")
    histtype[histname[-1]]="mupixelhits"
    #
    histname.append("hselmuon_npixellayers")
    histtype[histname[-1]]="mupixellayers"
    #
    histname.append("hselmuon_nstriphits")
    histtype[histname[-1]]="mustriphits"
    #
    histname.append("hselmuon_ntrackerlayers")
    histtype[histname[-1]]="mutrackerlayers"
    #
    ##
    #
    histname.append("hselmuon_osv_pt")
    histtype[histname[-1]]="mupt"
    #
    histname.append("hselmuon_osv_eta")
    histtype[histname[-1]]="mueta"
    #
    histname.append("hselmuon_osv_phi")
    histtype[histname[-1]]="muphi"
    #
    histname.append("hselmuon_osv_ch")
    histtype[histname[-1]]="much"
    #
    histname.append("hselmuon_osv_mindr")
    histtype[histname[-1]]="mumindr"
    #
    histname.append("hselmuon_osv_maxdr")
    histtype[histname[-1]]="mumaxdr"
    #
    histname.append("hselmuon_osv_chi2ndof")
    histtype[histname[-1]]="muchi2ndof"
    #
    histname.append("hselmuon_osv_type")
    histtype[histname[-1]]="mutype"
    #
    histname.append("hselmuon_osv_ecaliso")
    histtype[histname[-1]]="muecaliso"
    #
    histname.append("hselmuon_osv_ecalreliso")
    histtype[histname[-1]]="muecalreliso"
    #
    histname.append("hselmuon_osv_hcaliso")
    histtype[histname[-1]]="muhcaliso"
    #
    histname.append("hselmuon_osv_hcalreliso")
    histtype[histname[-1]]="muhcalreliso"
    #
    histname.append("hselmuon_osv_trackiso")
    histtype[histname[-1]]="mutrkiso"
    #
    histname.append("hselmuon_osv_trackreliso")
    histtype[histname[-1]]="mutrkreliso"
    #
    histname.append("hselmuon_osv_pfiso0p3all")
    histtype[histname[-1]]="mupfalliso"
    #
    histname.append("hselmuon_osv_pfreliso0p3all")
    histtype[histname[-1]]="mupfallreliso"
    #
    histname.append("hselmuon_osv_pfiso0p3chg")
    histtype[histname[-1]]="mupfchgiso"
    #
    histname.append("hselmuon_osv_pfreliso0p3chg")
    histtype[histname[-1]]="mupfchgreliso"
    #
    histname.append("hselmuon_osv_mindrjet")
    histtype[histname[-1]]="mumindrjet"
    #
    histname.append("hselmuon_osv_mindrpfc")
    histtype[histname[-1]]="mumindrpfc"
    #
    histname.append("hselmuon_osv_dxy")
    histtype[histname[-1]]="mudxy"
    #
    histname.append("hselmuon_osv_dxysig")
    histtype[histname[-1]]="mudxysig"
    #
    histname.append("hselmuon_osv_dz")
    histtype[histname[-1]]="mudz"
    #
    histname.append("hselmuon_osv_dzsig")
    histtype[histname[-1]]="mudzsig"
    #
    histname.append("hselmuon_osv_nsahits")
    histtype[histname[-1]]="musahits"
    #
    histname.append("hselmuon_osv_nsamatchedstats")
    histtype[histname[-1]]="musamatchedstats"
    #
    histname.append("hselmuon_osv_nmuhits")
    histtype[histname[-1]]="muhits"
    #
    histname.append("hselmuon_osv_nmuchambs")
    histtype[histname[-1]]="muchambs"
    #
    histname.append("hselmuon_osv_nmuCSCorDT")
    histtype[histname[-1]]="mucscdt"
    #
    histname.append("hselmuon_osv_nmumatches")
    histtype[histname[-1]]="mumatches"
    #
    histname.append("hselmuon_osv_nmumatchedstats")
    histtype[histname[-1]]="mumatchedstats"
    #
    histname.append("hselmuon_osv_nmuexpmatchedstats")
    histtype[histname[-1]]="muexpmatchedstats"
    #
    histname.append("hselmuon_osv_nmumatchedstatsmexp")
    histtype[histname[-1]]="mumatchedstatsmexp"
    #
    histname.append("hselmuon_osv_nmumatchedRPC")
    histtype[histname[-1]]="mumatchedrpc"
    #
    histname.append("hselmuon_osv_npixelhits")
    histtype[histname[-1]]="mupixelhits"
    #
    histname.append("hselmuon_osv_npixellayers")
    histtype[histname[-1]]="mupixellayers"
    #
    histname.append("hselmuon_osv_nstriphits")
    histtype[histname[-1]]="mustriphits"
    #
    histname.append("hselmuon_osv_ntrackerlayers")
    histtype[histname[-1]]="mutrackerlayers"
    #
    ##
    #
    histname.append("hselmuon_fourmu_pt")
    histtype[histname[-1]]="mupt"
    #
    histname.append("hselmuon_fourmu_eta")
    histtype[histname[-1]]="mueta"
    #
    histname.append("hselmuon_fourmu_phi")
    histtype[histname[-1]]="muphi"
    #
    histname.append("hselmuon_fourmu_ch")
    histtype[histname[-1]]="much"
    #
    histname.append("hselmuon_fourmu_mindr")
    histtype[histname[-1]]="mumindr"
    #
    histname.append("hselmuon_fourmu_maxdr")
    histtype[histname[-1]]="mumaxdr"
    #
    histname.append("hselmuon_fourmu_chi2ndof")
    histtype[histname[-1]]="muchi2ndof"
    #
    histname.append("hselmuon_fourmu_type")
    histtype[histname[-1]]="mutype"
    #
    histname.append("hselmuon_fourmu_ecaliso")
    histtype[histname[-1]]="muecaliso"
    #
    histname.append("hselmuon_fourmu_ecalreliso")
    histtype[histname[-1]]="muecalreliso"
    #
    histname.append("hselmuon_fourmu_hcaliso")
    histtype[histname[-1]]="muhcaliso"
    #
    histname.append("hselmuon_fourmu_hcalreliso")
    histtype[histname[-1]]="muhcalreliso"
    #
    histname.append("hselmuon_fourmu_trackiso")
    histtype[histname[-1]]="mutrkiso"
    #
    histname.append("hselmuon_fourmu_trackreliso")
    histtype[histname[-1]]="mutrkreliso"
    #
    histname.append("hselmuon_fourmu_pfiso0p3all")
    histtype[histname[-1]]="mupfalliso"
    #
    histname.append("hselmuon_fourmu_pfreliso0p3all")
    histtype[histname[-1]]="mupfallreliso"
    #
    histname.append("hselmuon_fourmu_pfiso0p3chg")
    histtype[histname[-1]]="mupfchgiso"
    #
    histname.append("hselmuon_fourmu_pfreliso0p3chg")
    histtype[histname[-1]]="mupfchgreliso"
    #
    histname.append("hselmuon_fourmu_mindrjet")
    histtype[histname[-1]]="mumindrjet"
    #
    histname.append("hselmuon_fourmu_mindrpfc")
    histtype[histname[-1]]="mumindrpfc"
    #
    histname.append("hselmuon_fourmu_dxy")
    histtype[histname[-1]]="mudxy"
    #
    histname.append("hselmuon_fourmu_dxysig")
    histtype[histname[-1]]="mudxysig"
    #
    histname.append("hselmuon_fourmu_dz")
    histtype[histname[-1]]="mudz"
    #
    histname.append("hselmuon_fourmu_dzsig")
    histtype[histname[-1]]="mudzsig"
    #
    histname.append("hselmuon_fourmu_nsahits")
    histtype[histname[-1]]="musahits"
    #
    histname.append("hselmuon_fourmu_nsamatchedstats")
    histtype[histname[-1]]="musamatchedstats"
    #
    histname.append("hselmuon_fourmu_nmuhits")
    histtype[histname[-1]]="muhits"
    #
    histname.append("hselmuon_fourmu_nmuchambs")
    histtype[histname[-1]]="muchambs"
    #
    histname.append("hselmuon_fourmu_nmuCSCorDT")
    histtype[histname[-1]]="mucscdt"
    #
    histname.append("hselmuon_fourmu_nmumatches")
    histtype[histname[-1]]="mumatches"
    #
    histname.append("hselmuon_fourmu_nmumatchedstats")
    histtype[histname[-1]]="mumatchedstats"
    #
    histname.append("hselmuon_fourmu_nmuexpmatchedstats")
    histtype[histname[-1]]="muexpmatchedstats"
    #
    histname.append("hselmuon_fourmu_nmumatchedstatsmexp")
    histtype[histname[-1]]="mumatchedstatsmexp"
    #
    histname.append("hselmuon_fourmu_nmumatchedRPC")
    histtype[histname[-1]]="mumatchedrpc"
    #
    histname.append("hselmuon_fourmu_npixelhits")
    histtype[histname[-1]]="mupixelhits"
    #
    histname.append("hselmuon_fourmu_npixellayers")
    histtype[histname[-1]]="mupixellayers"
    #
    histname.append("hselmuon_fourmu_nstriphits")
    histtype[histname[-1]]="mustriphits"
    #
    histname.append("hselmuon_fourmu_ntrackerlayers")
    histtype[histname[-1]]="mutrackerlayers"
    #
    ##
    #
    histname.append("hselmuon_fourmu_osv_pt")
    histtype[histname[-1]]="mupt"
    #
    histname.append("hselmuon_fourmu_osv_eta")
    histtype[histname[-1]]="mueta"
    #
    histname.append("hselmuon_fourmu_osv_phi")
    histtype[histname[-1]]="muphi"
    #
    histname.append("hselmuon_fourmu_osv_ch")
    histtype[histname[-1]]="much"
    #
    histname.append("hselmuon_fourmu_osv_mindr")
    histtype[histname[-1]]="mumindr"
    #
    histname.append("hselmuon_fourmu_osv_maxdr")
    histtype[histname[-1]]="mumaxdr"
    #
    histname.append("hselmuon_fourmu_osv_chi2ndof")
    histtype[histname[-1]]="muchi2ndof"
    #
    histname.append("hselmuon_fourmu_osv_type")
    histtype[histname[-1]]="mutype"
    #
    histname.append("hselmuon_fourmu_osv_ecaliso")
    histtype[histname[-1]]="muecaliso"
    #
    histname.append("hselmuon_fourmu_osv_ecalreliso")
    histtype[histname[-1]]="muecalreliso"
    #
    histname.append("hselmuon_fourmu_osv_hcaliso")
    histtype[histname[-1]]="muhcaliso"
    #
    histname.append("hselmuon_fourmu_osv_hcalreliso")
    histtype[histname[-1]]="muhcalreliso"
    #
    histname.append("hselmuon_fourmu_osv_trackiso")
    histtype[histname[-1]]="mutrkiso"
    #
    histname.append("hselmuon_fourmu_osv_trackreliso")
    histtype[histname[-1]]="mutrkreliso"
    #
    histname.append("hselmuon_fourmu_osv_pfiso0p3all")
    histtype[histname[-1]]="mupfalliso"
    #
    histname.append("hselmuon_fourmu_osv_pfreliso0p3all")
    histtype[histname[-1]]="mupfallreliso"
    #
    histname.append("hselmuon_fourmu_osv_pfiso0p3chg")
    histtype[histname[-1]]="mupfchgiso"
    #
    histname.append("hselmuon_fourmu_osv_pfreliso0p3chg")
    histtype[histname[-1]]="mupfchgreliso"
    #
    histname.append("hselmuon_fourmu_osv_mindrjet")
    histtype[histname[-1]]="mumindrjet"
    #
    histname.append("hselmuon_fourmu_osv_mindrpfc")
    histtype[histname[-1]]="mumindrpfc"
    #
    histname.append("hselmuon_fourmu_osv_dxy")
    histtype[histname[-1]]="mudxy"
    #
    histname.append("hselmuon_fourmu_osv_dxysig")
    histtype[histname[-1]]="mudxysig"
    #
    histname.append("hselmuon_fourmu_osv_dz")
    histtype[histname[-1]]="mudz"
    #
    histname.append("hselmuon_fourmu_osv_dzsig")
    histtype[histname[-1]]="mudzsig"
    #
    histname.append("hselmuon_fourmu_osv_nsahits")
    histtype[histname[-1]]="musahits"
    #
    histname.append("hselmuon_fourmu_osv_nsamatchedstats")
    histtype[histname[-1]]="musamatchedstats"
    #
    histname.append("hselmuon_fourmu_osv_nmuhits")
    histtype[histname[-1]]="muhits"
    #
    histname.append("hselmuon_fourmu_osv_nmuchambs")
    histtype[histname[-1]]="muchambs"
    #
    histname.append("hselmuon_fourmu_osv_nmuCSCorDT")
    histtype[histname[-1]]="mucscdt"
    #
    histname.append("hselmuon_fourmu_osv_nmumatches")
    histtype[histname[-1]]="mumatches"
    #
    histname.append("hselmuon_fourmu_osv_nmumatchedstats")
    histtype[histname[-1]]="mumatchedstats"
    #
    histname.append("hselmuon_fourmu_osv_nmuexpmatchedstats")
    histtype[histname[-1]]="muexpmatchedstats"
    #
    histname.append("hselmuon_fourmu_osv_nmumatchedstatsmexp")
    histtype[histname[-1]]="mumatchedstatsmexp"
    #
    histname.append("hselmuon_fourmu_osv_nmumatchedRPC")
    histtype[histname[-1]]="mumatchedrpc"
    #
    histname.append("hselmuon_fourmu_osv_npixelhits")
    histtype[histname[-1]]="mupixelhits"
    #
    histname.append("hselmuon_fourmu_osv_npixellayers")
    histtype[histname[-1]]="mupixellayers"
    #
    histname.append("hselmuon_fourmu_osv_nstriphits")
    histtype[histname[-1]]="mustriphits"
    #
    histname.append("hselmuon_fourmu_osv_ntrackerlayers")
    histtype[histname[-1]]="mutrackerlayers"
    #
    ##
    # Di-muon
    #
    histname.append("hdimuon_mass")
    histtype[histname[-1]]="dimumass"
    #
    histname.append("hdimuon_pt")
    histtype[histname[-1]]="dimupt"
    #
    histname.append("hdimuon_lxy")
    histtype[histname[-1]]="lxy"
    #
    histname.append("hdimuon_dr")
    histtype[histname[-1]]="dimudr"
    #
    histname.append("hdimuon_dphi")
    histtype[histname[-1]]="dimudphi"
    #
    histname.append("hdimuon_deta")
    histtype[histname[-1]]="dimudeta"
    #
    histname.append("hdimuon_detadphiratio")
    histtype[histname[-1]]="dimudetadphiratio"
    #
    histname.append("hdimuon_3dangle")
    histtype[histname[-1]]="dimu3dangle"
    #
    histname.append("hdimuon_dphisv")
    histtype[histname[-1]]="dimudphisv"
    #
    histname.append("hdimuon_detasv")
    histtype[histname[-1]]="dimudetasv"
    #
    histname.append("hdimuon_detadphisvratio")
    histtype[histname[-1]]="dimudetadphisvratio"
    #
    histname.append("hdimuon_3danglesv")
    histtype[histname[-1]]="dimu3danglesv"
    #
    ##
    # Di-muon from overlapping SV
    #
    histname.append("hdimuon_osv_mass")
    histtype[histname[-1]]="dimumass"
    #
    histname.append("hdimuon_osv_pt")
    histtype[histname[-1]]="dimupt"
    #
    histname.append("hdimuon_osv_lxy")
    histtype[histname[-1]]="lxy"
    #
    histname.append("hdimuon_osv_dr")
    histtype[histname[-1]]="dimudr"
    #
    histname.append("hdimuon_osv_dphi")
    histtype[histname[-1]]="dimudphi"
    #
    histname.append("hdimuon_osv_deta")
    histtype[histname[-1]]="dimudeta"
    #
    histname.append("hdimuon_osv_detadphiratio")
    histtype[histname[-1]]="dimudetadphiratio"
    #
    histname.append("hdimuon_osv_3dangle")
    histtype[histname[-1]]="dimu3dangle"
    #
    histname.append("hdimuon_osv_dphisv")
    histtype[histname[-1]]="dimudphisv"
    #
    histname.append("hdimuon_osv_detasv")
    histtype[histname[-1]]="dimudetasv"
    #
    histname.append("hdimuon_osv_detadphisvratio")
    histtype[histname[-1]]="dimudetadphisvratio"
    #
    histname.append("hdimuon_osv_3danglesv")
    histtype[histname[-1]]="dimu3danglesv"
    #
    ##
    # Four-muon from overlapping SV
    #
    histname.append("hfourmuon_osv_mass")
    histtype[histname[-1]]="fourmumass"
    #
    histname.append("hfourmuon_osv_pt")
    histtype[histname[-1]]="fourmupt"
    #
    histname.append("hfourmuon_osv_lxy")
    histtype[histname[-1]]="lxy"
    #
    histname.append("hfourmuon_osv_mindr")
    histtype[histname[-1]]="fourmumindr"
    #
    histname.append("hfourmuon_osv_maxdr")
    histtype[histname[-1]]="fourmumaxdr"
    #
    histname.append("hfourmuon_osv_mindphi")
    histtype[histname[-1]]="fourmumindphi"
    #
    histname.append("hfourmuon_osv_maxdphi")
    histtype[histname[-1]]="fourmumaxdphi"
    #
    histname.append("hfourmuon_osv_mindeta")
    histtype[histname[-1]]="fourmumindeta"
    #
    histname.append("hfourmuon_osv_maxdeta")
    histtype[histname[-1]]="fourmumaxdeta"
    #
    histname.append("hfourmuon_osv_mindetadphiratio")
    histtype[histname[-1]]="fourmumindetadphiratio"
    #
    histname.append("hfourmuon_osv_maxdetadphiratio")
    histtype[histname[-1]]="fourmumaxdetadphiratio"
    #
    histname.append("hfourmuon_osv_min3dangle")
    histtype[histname[-1]]="fourmumin3dangle"
    #
    histname.append("hfourmuon_osv_max3dangle")
    histtype[histname[-1]]="fourmumax3dangle"
    #
    histname.append("hfourmuon_osv_dphisv")
    histtype[histname[-1]]="fourmudphisv"
    #
    histname.append("hfourmuon_osv_detasv")
    histtype[histname[-1]]="dimudetasv"
    #
    histname.append("hfourmuon_osv_detadphisvratio")
    histtype[histname[-1]]="fourmudetadphisvratio"
    #
    histname.append("hfourmuon_osv_3danglesv")
    histtype[histname[-1]]="fourmu3danglesv"
    #
    ##
    # Four-muon
    #
    histname.append("hfourmuon_mass")
    histtype[histname[-1]]="fourmumass"
    #
    histname.append("hfourmuon_maxmass")
    histtype[histname[-1]]="fourmumaxmass"
    #
    histname.append("hfourmuon_minmass")
    histtype[histname[-1]]="fourmuminmass"
    #
    histname.append("hfourmuon_avgmass")
    histtype[histname[-1]]="fourmuavgmass"
    #
    histname.append("hfourmuon_reldmass")
    histtype[histname[-1]]="fourmureldmass"
    #
    histname.append("hfourmuon_pt")
    histtype[histname[-1]]="fourmupt"
    #
    histname.append("hfourmuon_minpt")
    histtype[histname[-1]]="fourmuminpt"
    #
    histname.append("hfourmuon_maxpt")
    histtype[histname[-1]]="fourmumaxpt"
    #
    histname.append("hfourmuon_minlxy")
    histtype[histname[-1]]="minlxy"
    #
    histname.append("hfourmuon_maxlxy")
    histtype[histname[-1]]="maxlxy"
    #
    histname.append("hfourmuon_mindr")
    histtype[histname[-1]]="fourmumindr"
    #
    histname.append("hfourmuon_maxdr")
    histtype[histname[-1]]="fourmumaxdr"
    #
    histname.append("hfourmuon_mindphi")
    histtype[histname[-1]]="fourmumindphi"
    #
    histname.append("hfourmuon_maxdphi")
    histtype[histname[-1]]="fourmumaxdphi"
    #
    histname.append("hfourmuon_mindeta")
    histtype[histname[-1]]="fourmumindeta"
    #
    histname.append("hfourmuon_maxdeta")
    histtype[histname[-1]]="fourmumaxdeta"
    #
    histname.append("hfourmuon_mindetadphiratio")
    histtype[histname[-1]]="fourmumindetadphiratio"
    #
    histname.append("hfourmuon_maxdetadphiratio")
    histtype[histname[-1]]="fourmumaxdetadphiratio"
    #
    histname.append("hfourmuon_min3dangle")
    histtype[histname[-1]]="fourmumin3dangle"
    #
    histname.append("hfourmuon_max3dangle")
    histtype[histname[-1]]="fourmumax3dangle"
    #
    histname.append("hfourmuon_mindphisv")
    histtype[histname[-1]]="fourmumindphisv"
    #
    histname.append("hfourmuon_maxdphisv")
    histtype[histname[-1]]="fourmumaxdphisv"
    #
    histname.append("hfourmuon_mindetasv")
    histtype[histname[-1]]="fourmumindetasv"
    #
    histname.append("hfourmuon_maxdetasv")
    histtype[histname[-1]]="fourmumaxdetasv"
    #
    histname.append("hfourmuon_mindetadphisvratio")
    histtype[histname[-1]]="fourmumindetadphisvratio"
    #
    histname.append("hfourmuon_maxdetadphisvratio")
    histtype[histname[-1]]="fourmumaxdetadphisvratio"
    #
    histname.append("hfourmuon_min3danglesv")
    histtype[histname[-1]]="fourmumin3danglesv"
    #
    histname.append("hfourmuon_max3danglesv")
    histtype[histname[-1]]="fourmumax3danglesv"
    #
    return histname,histtype,histname2d,histtype2d

def histInitialization():
    nbins    = dict()
    low      = dict()
    high     = dict()
    xtitle   = dict()
    ytitle   = dict()
    labels   = dict()
    variable = dict()
    hist1dDefinition(nbins, low, high, xtitle, ytitle, labels, variable)
    
    nbinsX      = dict()
    lowX        = dict()
    highX       = dict()
    xtitle2d    = dict()
    nbinsY      = dict()    
    lowY        = dict()
    highY       = dict()
    ytitle2d    = dict()
    ztitle2d    = dict()
    variablesXY = dict()
    hist2dDefinition(nbinsX, lowX, highX, xtitle2d, nbinsY, lowY, highY, ytitle2d, ztitle2d, variablesXY)

    histname = []
    histtype = dict()
    histname2d = []
    histtype2d = dict()
    histname, histtype, histname2d, histtype2d = histBooking()

    hists1d     = []
    variables1d = dict()
    for hn in histname:
        if histtype[hn] in labels.keys():
            th = H1(hn,nbins[histtype[hn]],low[histtype[hn]],high[histtype[hn]],xtitle[histtype[hn]],ytitle[histtype[hn]],labels[histtype[hn]])
        else:
            th = H1(hn,nbins[histtype[hn]],low[histtype[hn]],high[histtype[hn]],xtitle[histtype[hn]],ytitle[histtype[hn]])
        variables1d[hn] = variable[histtype[hn]]
        hists1d.append(th)
    hists2d     = []
    variables2d = dict()
    for hn in histname2d:
        th = H2(hn,
                nbinsX[histtype2d[hn]],lowX[histtype2d[hn]],highX[histtype2d[hn]],xtitle2d[histtype2d[hn]],
                nbinsY[histtype2d[hn]],lowY[histtype2d[hn]],highY[histtype2d[hn]],ytitle2d[histtype2d[hn]],
                ztitle2d[histtype2d[hn]])
        variables2d[hn] = variablesXY[histtype2d[hn]]
        hists2d.append(th)
    return hists1d, variables1d, hists2d, variables2d
