#include "TFile.h"
#include "TDirectory.h"
#include "RooDataSet.h"
#include "TKey.h"
#include "RooAbsData.h"
#include <iostream>

void printRooDataSetEvents(const char* fileName) {
    // Open ROOT file
    TFile *file = TFile::Open(fileName);

    if (!file || file->IsZombie()) {
        std::cerr << "No file: " << fileName << std::endl;
        return;
    }

    // Access "toys"
    TDirectory* dir = (TDirectory*)file->Get("toys");
    if (!dir) {
        std::cerr << "No 'toys' in file." << std::endl;
        return;
    }

    TIter next(dir->GetListOfKeys());
    TKey *key;
    TH1F *h = new TH1F("CMS_channel", "CMS_channel==29;Number of events in toy RooDataSet; Number of toys", 10, 0, 10);
    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();

        RooDataSet* data = dynamic_cast<RooDataSet*>(obj);
	RooDataSet* filteredData = dynamic_cast<RooDataSet*>(data->reduce("CMS_channel==29 && mfit > 6.4 && mfit < 7.6"));
        if (data) {
            std::cout << "RooDataSet: " << filteredData->GetName() 
                      << " contains " << filteredData->numEntries() << " events." << std::endl;
	    h->Fill(filteredData->numEntries());
        }
    }

    // Plot
    TCanvas c1("c1","");
    h->Draw("HIST");
    c1.SaveAs("channel29.png");

    // Cerrar el archivo ROOT
    file->Close();
}

