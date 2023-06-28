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


int main(int argc, char **argv) {
  // Arguments
  const char* outdir      = ( argc > 1  ? argv[1]            : "temp_data" );
  TString year            = ( argc > 2  ? argv[2]            : "2022" );
  int startFile           = ( argc > 3  ? char2int(argv[3])  : 0 );
  int nFiles              = ( argc > 4  ? char2int(argv[4])  : 1000000 ); // Large number as default to run over all files

  // Sample list: Data
  std::vector<TString> dataFiles;
  const fs::path dataDir{"/ceph/cms/store/user/legianni/testRAWScouting_0/ScoutingPFRun3/crab_skim_2022D_0/230613_184336/0000/"};
  unsigned int iFile=0;
  for (auto const& dataFile : fs::directory_iterator{dataDir}) {
    if (iFile<startFile) {
      iFile++;
      continue;
    }
    if (!TString(dataFile.path()).Contains(".root"))
      continue;
    else {
      dataFiles.push_back(TString(dataFile.path()));
    }
    iFile++;
    if (iFile == startFile+nFiles)
      break;
  }

  std::cout << "################################## \n";
  std::cout << "Number of files to process: " << dataFiles.size() << "\n";
  std::cout << "################################## \n";
  run3ScoutingLooper(dataFiles, year, "Data", outdir, "_"+std::to_string(startFile)+"To"+std::to_string(startFile+nFiles-1));

  return 0;
}

