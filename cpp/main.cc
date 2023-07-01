#include "run3ScoutingLooper.C"

unsigned char2unsigned(const char *c) {
  char p = *c;
  unsigned res = 0;
  while (p) {
      res = res*10 + (p - '0');
      c++;
      p = *c;
  }
  return res;
}
int char2int(const char *c) {
  return (*c == '-') ? -char2unsigned(c+1) : char2unsigned(c);
}


std::vector<TString> getFiles(std::string inputDir, int startFile, int nFiles) {
  std::vector<TString> files;
  const fs::path dir{inputDir};
  unsigned int iFile=0;
  for (auto const& file : fs::directory_iterator{dir}) {
    if (iFile<startFile) {
      iFile++;
      continue;
    }
    if (!TString(file.path()).Contains(".root"))
      continue;
    else {
      files.push_back(TString(file.path()));
    }
    iFile++;
    if (iFile == startFile+nFiles)
      break;
  }
  return files;
}

int main(int argc, char **argv) {
  // Arguments
  const char* outdir      = ( argc > 1  ? argv[1]            : "temp_data" );
  TString year            = ( argc > 2  ? argv[2]            : "2022" );
  TString sampleArg       = ( argc > 3  ? argv[3]            : "Data" );
  int startFile           = ( argc > 4  ? char2int(argv[4])  : 0 );
  int nFiles              = ( argc > 5  ? char2int(argv[5])  : 1000000 ); // Large number as default to run over all files

  std::vector<TString> files;
  TString process;
  // Sample list: Data
  if ( sampleArg=="Data" ) {
    files = getFiles("/ceph/cms/store/user/legianni/testRAWScouting_0/ScoutingPFRun3/crab_skim_2022D_0/230613_184336/0000/", startFile, nFiles);
    process = "Data";
  }
  // Sample list: Signal
  if ( sampleArg=="Signal" ) {
    files = getFiles("/ceph/cms/store/user/isuarez/ProjectMetis/DarkShower_ScenarioA_default_Run3Summer22GS_v0p30_AODSIM_v0p30/", startFile, nFiles);
    process = "Signal";
  }

  std::cout << "################################## \n";
  std::cout << "Number of files to process: " << files.size() << "\n";
  std::cout << "################################## \n";
  run3ScoutingLooper(files, year, process, outdir, "_"+std::to_string(startFile)+"To"+std::to_string(startFile+nFiles-1));

  return 0;
}

