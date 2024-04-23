{
  gROOT->ProcessLine(".L ./cpp/fit_dimuon.C+");  // Macro that performs the fits
  gROOT->ProcessLine(".L ./cpp/tools/param.C+");  // Macro that handles the spline

  bool useData = true;
  bool useSignalMC = true;
  bool mergeEras = true;
  bool writeWS = false;
  TString model = "HTo2ZdTo2mu2x";
  
  // Dir with the RooDataSets
  TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Apr-03-2024_onlySignal";

  // Names of the search regions we want to parametrize
  vector<TString> dNames = { };
  dNames.push_back("d_Dimuon_full_inclusive");

  // Eras (to be uncommented when adding 2023 and splitting in eras)
  vector<TString> eras;
  eras.push_back("2022");
  eras.push_back("2022postEE");
  //years.push_back("2023"); // Currently not available
  //years.push_back("2023BPix"); // Currently not available

  TString sigTemplate;
  //vector<TString> sigsamples = { };


  // Signals (To be modified: Needs to be more general but this is provisional)
  vector<float> sigMass;
  vector<float> sigCtau;
  vector<float> selectedPoint;
  if ( model=="HTo2ZdTo2mu2x" ) {
    sigMass = {0.5, 0.7, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 22.0, 24.0, 30.0, 34.0, 40.0, 44.0, 50.0};
    //sigMass = {0.5, 0.7, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0};
    //vector<float> sigCtau = {1, 10, 100, 1000};
    sigCtau = {1, 10, 100};
    sigTemplate = "Signal_HTo2ZdTo2mu2x_MZd-%s_ctau-%smm";
    for (unsigned int t=0; t<sigCtau.size(); t++) {
      for (unsigned int m=0; m<sigMass.size(); m++) {
        if ( (sigMass[m] < 1.0 && sigCtau[t] > 10) || (sigMass[m] < 30.0 && sigCtau[t] > 100))
           selectedPoint.push_back(false);
        else
           selectedPoint.push_back(true);
      }
    }
  }

  std::cout << "Control" << std::endl;

  // Loop over signals 
  vector<vector<RooDataSet>> mmumu_sigs {{}};
  for ( unsigned int d=0; d<dNames.size(); d++ ) {
    for ( unsigned int iera=0; iera<eras.size(); iera++ ) {
      TString era = eras[iera];
      TString year = "2022"; // hardcoded
      for (unsigned int t=0; t<sigCtau.size(); t++) {
        for (unsigned int m=0; m<sigMass.size(); m++) {
          unsigned int isample = t*sigMass.size() + m;
          if (iera == 0) {
            vector<RooDataSet> tds_aux{};
            mmumu_sigs.push_back( tds_aux );
          }
          if (!selectedPoint[isample])
            continue;
          TString massString = Form("%.1f",sigMass[m]);
          massString.ReplaceAll(".", "p");
          TString ctauString = Form("%.0f",sigCtau[t]);
          TString sample = Form(sigTemplate.Data(), massString.Data(),ctauString.Data());
          TString inFile = Form("%s/histograms_%s_%s_%s_0.root",inDir.Data(),sample.Data(),era.Data(),year.Data());
          TFile fin(inFile);
          RooDataSet *tds = (RooDataSet*) fin.Get(dNames[d])->Clone();
          tds->SetName(dNames[d]+"_"+sample+"_"+year+"_"+era);
          mmumu_sigs[isample].push_back( *tds );
          std::cout << "For isample: " <<  isample << std::endl;
          std::cout << "Reading signal file: " <<  inFile << ", with dataset with entries: " << tds->sumEntries() << std::endl;
          fin.Close();
          cout << mmumu_sigs[isample][iera].numEntries() << endl;
        }
      }
    }

   
   TFile *finit = new TFile("utils/signalFitParameters.root", "RECREATE");
   finit->Close();

   if (mergeEras) {
     TString outDir = "paramResults_allEras";
     vector<TSpline3> splineList;
     vector<TGraph> graphList;
     vector<RooDataSet> mmumu_sig_merged = {};
     for (unsigned int t=0; t<sigCtau.size(); t++) {

       TString ctauString = Form("%.0f",sigCtau[t]);

       // Graphs
       TGraph *gsigma = new TGraph();
       TGraph *gnL = new TGraph();
       TGraph *gnR = new TGraph();
       TGraph *gaL = new TGraph();
       TGraph *gaR = new TGraph();
       TGraph *gmean = new TGraph();

       gsigma->SetName(Form("gsigma_%smm",ctauString.Data()));
       gnL->SetName(Form("gnL_%smm",ctauString.Data()));
       gnR->SetName(Form("gnR_%smm",ctauString.Data()));
       gaL->SetName(Form("gaL_%smm",ctauString.Data()));
       gaR->SetName(Form("gaR_%smm",ctauString.Data()));
       gmean->SetName(Form("gmean_%smm",ctauString.Data()));

       for (unsigned int m=0; m<sigMass.size(); m++) {

         unsigned int isample = t*sigMass.size() + m;
         TString massString = Form("%.1f",sigMass[m]);
         TString sample = Form(sigTemplate.Data(), massString.Data(),ctauString.Data());
         if (!selectedPoint[isample])
           continue;
         // Create worksapce, import data and model
         std::cout << "Creating container workspace" << std::endl;
         RooWorkspace wfit("wfit","workspace"); 
         int idx = 0;
         for (unsigned int iera=0; iera<eras.size(); iera++ ) {
           if (iera==0) {
             mmumu_sig_merged.push_back(mmumu_sigs[isample][iera]);
           } else {
             idx = mmumu_sig_merged.size() - 1;
             mmumu_sig_merged[idx].append(mmumu_sigs[isample][iera]);
           }
           cout << iera << " ...filling...  "<< mmumu_sigs[isample][iera].GetName() << endl;
         }  
         mmumu_sig_merged[idx].SetName(dNames[d]+"_"+sample+"_allEras");
         cout << "Merged dataset for signal " << sample << " with entries " << mmumu_sig_merged[idx].sumEntries() << endl;
         // Fit invariant mass
         std::cout << "Prepare to fit..." << std::endl;
         if (dNames[d].BeginsWith("d_FourMu_")) {
           fitmass(mmumu_sig_merged[idx], "Signal", false, true, true, sample , 125., wfit, true, "dcbfastg", outDir);
         } else {
           fitmass(mmumu_sig_merged[idx], "Signal", false, true, true, sample, sigMass[m], wfit, false, "dcbfastg", outDir);
         }
    
         // Print workspace contents
         std::cout << "Temporal workspace contents: " << std::endl;
         wfit.Print();

         // Save the workspace into a ROOT file
         //if ( writeWS ) {
         //  TString fwsname = Form("%s/%s_workspace.root",outDir.Data(),mmumu_sig_merged[isample].GetName());
         //  TFile *fws = new TFile(fwsname, "RECREATE");
         //  fws->cd();
         //  cout << "Writing workspace..." << endl;
         //  wfit.Write();
         //  fws->Close();
         //}

         // Get the values
         RooRealVar *sigma = wfit.var("sigma_ch-1");
         RooRealVar *nL = wfit.var("nL_ch-1");
         RooRealVar *nR = wfit.var("nR_ch-1");
         RooRealVar *aL = wfit.var("alphaL_ch-1");
         RooRealVar *aR = wfit.var("alphaR_ch-1");
         RooRealVar *mean = wfit.var("mean_ch-1");

         // Fill graphs
         gsigma->AddPoint(sigMass[m], sigma->getVal());
         gnL->AddPoint(sigMass[m], nL->getVal());
         gnR->AddPoint(sigMass[m], nR->getVal());
         gaL->AddPoint(sigMass[m], aL->getVal());
         gaR->AddPoint(sigMass[m], aR->getVal());
         gmean->AddPoint(sigMass[m], mean->getVal());

       }

       // Create the Tspline's for the parameters
       TSpline3 *splines = new TSpline3("splines",gsigma->GetX(), gsigma->GetY(), gsigma->GetN(), "b2e2", 0, 0); 
       splines->SetName(Form("splines_%smm", ctauString.Data()));
       TSpline3 *splinenL = new TSpline3("splinenL",gnL->GetX(), gnL->GetY(), gnL->GetN(), "b2e2", 0, 0); 
       splinenL->SetName(Form("splinenL_%smm", ctauString.Data()));
       TSpline3 *splinenR = new TSpline3("splinenR",gnR->GetX(), gnR->GetY(), gnR->GetN(), "b2e2", 0, 0); 
       splinenR->SetName(Form("splinenR_%smm", ctauString.Data()));
       TSpline3 *splineaL = new TSpline3("splineaL",gaL->GetX(), gaL->GetY(), gaL->GetN(), "b2e2", 0, 0); 
       splineaL->SetName(Form("splineaL_%smm", ctauString.Data()));
       TSpline3 *splineaR = new TSpline3("splineaR",gaR->GetX(), gaR->GetY(), gaR->GetN(), "b2e2", 0, 0); 
       splineaR->SetName(Form("splineaR_%smm", ctauString.Data()));
       TSpline3 *splinem = new TSpline3("splinem", gmean->GetX(), gmean->GetY(), gmean->GetN(), "b2e2", 0, 0); 
       splinem->SetName(Form("splinem_%smm", ctauString.Data()));

        // Write splines
        TFile *fs = new TFile("utils/signalFitParameters.root", "UPDATE");
        fs->cd();
        splines->Write();
        splinem->Write();
        splinenL->Write();
        splinenR->Write();
        splineaL->Write();
        splineaR->Write();
        gsigma->Write();
        gmean->Write();
        gnL->Write();
        gnR->Write();
        gaL->Write();
        gaR->Write();
        fs->Close();

       //splineList.push_back(*splines); splineList.push_back(*splinem); splineList.push_back(*splinenL); splineList.push_back(*splinenR); splineList.push_back(*splineaL); splineList.push_back(*splineaR);
       //graphList.push_back(*gsigma); graphList.push_back(*gmean); graphList.push_back(*gnL); graphList.push_back(*gnR); graphList.push_back(*gaL); graphList.push_back(*gaR);

     }
   }

   
   //for ( int isample=0; isample<sigMass.size()+sigCtau.size(); isample++ ) {
   //  mmumu_sigs[isample].clear();
   //} 
   //mmumu_sigs.clear();
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
