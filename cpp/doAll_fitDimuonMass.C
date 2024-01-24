{
  gROOT->ProcessLine(".L ./cpp/fit_dimuon.C+");  // Macro that performs the selection

  bool useData = true;
  //bool useData = false;
  bool useSignalMC = true;
  //bool useSignalMC = true;
  bool writeWS = true;
  //asdf
  float mF = 350.0;
  float mL = 2000.0;

  TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Jan-22-2024_noFourMuonMassDiffSel/";
  //std::string inDir_str = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Dec-14-2023_allCuts_v2/";
  vector<TString> dNames = { };
  dNames.push_back("d_FourMu_sep");
  dNames.push_back("d_FourMu_osv");
  dNames.push_back("d_Dimuon_lxy0p0to0p5_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy0p0to0p5_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy0p0to0p5_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy0p0to0p5_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy0p5to2p7_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy0p5to2p7_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy0p5to2p7_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy0p5to2p7_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy2p7to6p5_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy2p7to6p5_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy2p7to6p5_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy2p7to6p5_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy6p5to11p0_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy6p5to11p0_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy6p5to11p0_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy6p5to11p0_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy11p0to16p0_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy11p0to16p0_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy11p0to16p0_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy11p0to16p0_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy16p0to70p0_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy16p0to70p0_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy16p0to70p0_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy16p0to70p0_iso1_pthigh");
  //dNames.push_back("");

  /* (to be uncommented when adding 2023 and splitting in eras)
  vector<TString> years = { };
  years.push_back("2022");
  years.push_back("2022EE");
  years.push_back("2023");
  */

  vector<TString> samples = { };
  vector<TString> sigmodels = { };
  vector<TString> sigsamples = { };
  vector<float> sigmasses_2mu = { };
  vector<float> sigmasses_4mu = { };

  if ( useData ) {
    samples.push_back("Data");
  }

  // Signals (To be modified: Needs to be more general but this is provisional)
  //vector<TString> sigModels = {"Signal_HTo2ZdTo2mu2x_MZd-%s_ctau-%smm"};
  vector<float> sigMass = {2.0, 7.0};
  vector<float> sigCtau = {1, 10, 100};
  for ( unsigned int m=0; m<sigMass.size(); m++ ) {
    TString massString = Form("%.1f",sigMass[m]); 
    massString.ReplaceAll(".", "p");
    for ( unsigned int t=0; t<sigCtau.size(); t++ ) {
      TString ctauString = Form("%.0f",sigCtau[t]); 
      sigsamples.push_back(Form("Signal_HTo2ZdTo2mu2x_MZd-%s_ctau-%smm",massString.Data(),ctauString.Data()));
      sigmasses_2mu.push_back(sigMass[m]);
      sigmasses_4mu.push_back(125.); // hardcoded: to be revisited
      std::cout << Form("Signal_HTo2ZdTo2mu2x_MZd-%s_ctau-%smm",massString.Data(),ctauString.Data()) << std::endl;
    }
  }

  /*
  if ( !(useSignalMC) ) {
    sigMass = {mF};
    float tm = mF;
    while ( tm < mL ) {
      if(tm<400.00) tm+=5.0;
      else if(tm<700) tm+=10.0;
      else if(tm<1000.0) tm+=15.0;
      else if (tm<1500.0) tm+=25.0;
      else tm+=50.0;
      sigMass.push_back(tm);
    }
  }
  */

  for ( unsigned int d=0; d<dNames.size(); d++ ) {
    std::cout << "Loading dataset: " << dNames[d] << std::endl;
    TString year = "2022";
    RooDataSet *mmumu_bkg;
    TSystemDirectory dir(inDir, inDir);
    TList* files = dir.GetListOfFiles(); 
    int idata = 0;
    for (const auto& file : *files) {
      std::cout << "File: " << file->GetName() << std::endl;
      TString filename = file->GetName();
      if (!filename.BeginsWith("histograms_Data") || filename.EndsWith("all.root"))
        continue;
      TString inFile = Form("%s/%s",inDir.Data(),filename.Data());
      TFile fin(inFile, "READ");
      if (idata == 0) {
        mmumu_bkg = (RooDataSet*) fin.Get(dNames[d])->Clone();
      } else {
        RooDataSet *tds_other = (RooDataSet*) fin.Get(dNames[d])->Clone();
        mmumu_bkg->append( *tds_other );
      }
      fin.Close();
      idata++;
    }
   
    mmumu_bkg->SetName(dNames[d]+"_data_"+year);
     
    //RooDataSet mmumu_bkg = getRooDataset(inDir_str, "histograms_Data", "all.root", year, dNames[d]);
    //TString inFile = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Dec-14-2023_allCuts_v2/histograms_DataF_2022_100.root";
    //TString inFile = Form("%s/histograms_Data_%s_%s.root",inDir.Data(),year.Data());
    //TFile fin(inFile, "READ");
    //std::cout << "File loaded: " << inFile << std::endl;
    //fin.ls();
    //mmumu_sig = (RooDataSet*) fin.Get(dNames[d])->Clone();
    //mmumu_bkg = (RooDataSet*) fin.Get(dNames[d])->Clone();
    //fin.Close();

    // Signal
    vector<RooDataSet> mmumu_sig = { }; 
    for ( int isample=0; isample<sigsamples.size(); isample++ ) {
      TString sample = sigsamples[isample];
      cout<<"Sample: "<< sample << endl;
      if ( useSignalMC ) {
        TString inFile = Form("%s/histograms_%s_%s_0.root",inDir.Data(),sample.Data(),year.Data());
        std::cout << "Reading signal file: " <<  inFile << std::endl;
        TFile fin(inFile);
        RooDataSet *tds = (RooDataSet*) fin.Get(dNames[d])->Clone();
        tds->SetName(dNames[d]+"_"+sample+"_"+year);
        mmumu_sig.push_back( *tds );
        fin.Close();
      } else {
	RooDataSet *tds = (RooDataSet*) mmumu_bkg->emptyClone(dNames[d]+"_"+sample+"_"+year, "");
	mmumu_sig.push_back( *tds );
      }

      // Create worksapce, import data and model
      std::cout << "Creating workspace" << std::endl;
      TString outDir = "fitResults";
      RooWorkspace wfit("wfit","workspace"); 

      // Fit invariant mass
      std::cout << "Prepare to fit..." << std::endl;
      if (dNames[d].BeginsWith("d_FourMu_")) {
        fitmass(mmumu_sig[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_4mu[isample], wfit, "dcbfastg", outDir);
        fitmass(*mmumu_bkg, "Background", true, false, false, sigsamples[isample], sigmasses_4mu[isample], wfit, "", outDir); 
      } else {
        fitmass(mmumu_sig[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], wfit, "dcbfastg", outDir);
        fitmass(*mmumu_bkg, "Background", true, false, false, sigsamples[isample], sigmasses_2mu[isample], wfit, "", outDir); 
      }
    
      // Print workspace contents
      std::cout << "Workspace contents: " << std::endl;
      wfit.Print();

      // Save the workspace into a ROOT file
      if ( writeWS ) {
        TString fwsname = Form("%s/%s_workspace.root",outDir.Data(),mmumu_sig[isample].GetName());
        TFile *fws = new TFile(fwsname, "RECREATE");
        fws->cd();
        cout << "Writing workspace..." << endl;
        wfit.Write();
        fws->Close();
      }
      cout<<endl;
    } 
  }

  /*
  if ( !mergeYears ) {
    for ( int iyear=0; iyear<years.size(); iyear++ ) {
	for ( int isample=0; isample<sigsamples.size(); isample++ ) {
	  TString sample = sigsamples[isample];
	  cout<<"Sample: "<<sample<<endl;
	  if ( useSignalMC ) {
	    TString inFile = Form("%s/output_%s_%s.root",inDir.Data(),sample.Data(),year.Data());
	    TFile fin(inFile);
	    RooDataSet *tds = (RooDataSet*) fin.Get(dNames[d])->Clone();
	    tds->SetName(dNames[d]+"_"+sample+"_"+year);
	    mmumu_sig.push_back( *tds );
	    fin.Close();
	  }
	  else {
	    RooDataSet *tds = (RooDataSet*) mmumu_bkg->emptyClone(dNames[d]+"_"+sample+"_"+year, "");
	    mmumu_sig.push_back( *tds );
	  }
	  // Create workspace, import data and model
	  TString outDir = "fitResults";
	  RooWorkspace wfit("wfit","workspace");

	  fitmass(mmumu_sig[isample], sample, false, true, false, sigmodels[isample], sigmasses[isample], wfit, "dcbfastg", outDir);
	  if ( samples.size() > 0 ) {
	    if ( useData )
	      fitmass(*mmumu_bkg, "Background", true, false, false, sigmodels[isample], sigmasses[isample], wfit, "", outDir);
	    else
	      fitmass(*mmumu_bkg, "Background", false, false, false, sigmodels[isample], sigmasses[isample], wfit, "", outDir);
	  }

	  // Print workspace contents
	  wfit.Print();

	  if ( writeWS ) {
	    // Save the workspace into a ROOT file
	    TString fwsname = Form("%s/%s_workspace.root",outDir.Data(),mmumu_sig[isample].GetName());
	    TFile *fws = new TFile(fwsname, "RECREATE");
	    fws->cd();
	    cout << "Writing workspace..." << endl;
	    wfit.Write();
	    fws->Close();
	  }
	  cout<<endl;
	}
      }
    }
  }
  else {
    for ( unsigned int d=0; d<dNames.size(); d++ ) {  
      RooDataSet *mmumu_bkg;
      cout<<endl;
      cout<<"Dataset: "<<d<<endl;
      cout<<"------------"<<endl;
      for ( int isample=0; isample<samples.size(); isample++ ) {
      	TString sample = samples[isample];
      	cout<<"Sample: "<<sample<<endl;
      	for ( int iyear=0; iyear<years.size(); iyear++ ) {
      	  TString year = years[iyear];
      	  cout<<endl;
      	  cout<<"Year: "<<year<<endl;
      	  TString inFile = Form("%s/output_%s_%s.root",inDir.Data(),sample.Data(),year.Data());
      	  TFile fin(inFile);
      	  if ( isample==0 && iyear==0 ) {
      	    mmumu_bkg = (RooDataSet*) fin.Get(dNames[d])->Clone();
      	  }
      	  else {
      	    RooDataSet *tds_other = (RooDataSet*) fin.Get(dNames[d])->Clone();
      	    mmumu_bkg->append( *tds_other );
      	  }
      	  fin.Close();
      	}
      }
      if ( useData )
      	mmumu_bkg->SetName(dNames[d]+"_data_allyears");
      else
      	mmumu_bkg->SetName(dNames[d]+"_BGMC_allyears");
      vector<RooDataSet*> mmumu_sig = { };
      for ( int isample=0; isample<sigsamples.size(); isample++ ) {
	TString sample = sigsamples[isample];
	cout<<"Sample: "<<sample<<endl;
	if ( useSignalMC ) {
	  for ( int iyear=0; iyear<years.size(); iyear++ ) {
	    TString year = years[iyear];
	    cout<<endl;
	    cout<<"Year: "<<year<<endl;
	    TString inFile = Form("%s/output_%s_%s.root",inDir.Data(),sample.Data(),year.Data());
	    TFile fin(inFile);
	    RooDataSet *tds;
	    if ( iyear==0 ) {
	      tds = (RooDataSet*) fin.Get(dNames[d])->Clone();
	      mmumu_sig.push_back(tds);
	    }
	    else {
	      RooDataSet *tds_other = (RooDataSet*) fin.Get(dNames[d])->Clone();
	      mmumu_sig[isample]->append( *tds_other );
	    }
	    fin.Close();
	  }
	  mmumu_sig[isample]->SetName(dNames[d]+"_"+sample+"_allyears");
	}
	else { 
	  mmumu_sig.push_back( (RooDataSet*) mmumu_bkg->emptyClone(dNames[d]+"_"+sample+"_allyears", "") );
	}
	// Create workspace, import data and model
	TString outDir = "fitResults_fullRange";
	RooWorkspace wfit("wfit","workspace");

	fitmass(*mmumu_sig[isample], sample, false, true, false, sigmodels[isample], sigmasses[isample], wfit, "dcbfastg", outDir);
	if ( samples.size() > 0 ) {
	  if ( useData )
	    fitmass(*mmumu_bkg, "Background", true, false, false, sigmodels[isample], sigmasses[isample], wfit, "", outDir);
	  else
	    fitmass(*mmumu_bkg, "Background", false, false, false, sigmodels[isample], sigmasses[isample], wfit, "", outDir);
	}

	// Print workspace contents
	wfit.Print();

	if ( writeWS ) {
	  // Save the workspace into a ROOT file
	  TString fwsname = Form("%s/%s_workspace.root",outDir.Data(),mmumu_sig[isample]->GetName());
	  TFile *fws = new TFile(fwsname, "RECREATE");
	  fws->cd();
	  cout << "Writing workspace..." << endl;
	  wfit.Write();
	  fws->Close();
	}
	cout<<endl;
      }
    }
  }
  */
}
