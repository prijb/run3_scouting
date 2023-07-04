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
    files = getFiles("/ceph/cms/store/user/legianni/testRAWScouting_2/ScoutingPFRun3/crab_skim__2022D_2/230703_062612/0000/", startFile, nFiles); // 999 files
    //files = getFiles("/ceph/cms/store/user/legianni/testRAWScouting_2/ScoutingPFRun3/crab_skim__2022D_2/230703_062612/0001/", startFile, nFiles); // 5 files
    process = "Data";
  }
  //
  // Sample list: Signal
  if ( sampleArg=="Signal_ScenA_v0p30" ) {
    files = getFiles("/ceph/cms/store/user/legianni/testRAWScouting_2/DarkShower_ScenarioA_default_Run3Summer22GS_v0p30_AODSIM_v0p30/crab_skim__2022_DarkShower_ScenarioA_default_Run3Summer22GS_v0p30_AODSIM_v0p30/230703_184107/0000/", startFile, nFiles); // 61 files
    process = "Signal_ScenA_v0p30";
  }
  if ( sampleArg=="Signal_ScenA_v1p3" ) {
    files = getFiles("/ceph/cms/store/user/legianni/testRAWScouting_2/DarkShower_ScenarioA_default_Run3Summer22GS_v1p3_AODSIM_v1p3/crab_skim__2022_DarkShower_ScenarioA_default_Run3Summer22GS_v1p3_AODSIM_v1p3/230703_184125/0000/", startFile, nFiles); // 69 files
    process = "Signal_ScenA_v1p3";
  }
  if ( sampleArg=="Signal_ScenA_v1p4" ) {
    files = getFiles("/ceph/cms/store/user/legianni/testRAWScouting_2/DarkShower_ScenarioA_default_Run3Summer22GS_v1p4_AODSIM_v1p4/crab_skim__2022_DarkShower_ScenarioA_default_Run3Summer22GS_v1p4_AODSIM_v1p4/230703_184144/0000/", startFile, nFiles); // 95 files
    process = "Signal_ScenA_v1p4";
  }
  if ( sampleArg=="Signal_ScenA_v1p5" ) {
    files = getFiles("/ceph/cms/store/user/legianni/testRAWScouting_2/DarkShower_ScenarioA_default_Run3Summer22GS_v1p5_AODSIM_v1p5/crab_skim__2022_DarkShower_ScenarioA_default_Run3Summer22GS_v1p5_AODSIM_v1p5/230703_184201/0000/", startFile, nFiles); // 93 files
    process = "Signal_ScenA_v1p5";
  }
  if ( sampleArg=="Signal_ScenB_v0p32" ) {
    files = getFiles("/ceph/cms/store/user/legianni/testRAWScouting_2/DarkShower_ScenarioB_default_Run3Summer22GS_v0p32_AODSIM_v0p32/crab_skim__2022_DarkShower_ScenarioB_default_Run3Summer22GS_v0p32_AODSIM_v0p32/230703_184218/0000/", startFile, nFiles); // 48 files
    process = "Signal_ScenB_v0p32";
  }
  if ( sampleArg=="Signal_ScenC_v0p34" ) {
    files = getFiles("/ceph/cms/store/user/legianni/testRAWScouting_2/DarkShower_ScenarioC_default_Run3Summer22GS_v0p34_AODSIM_v0p34/crab_skim__2022_DarkShower_ScenarioC_default_Run3Summer22GS_v0p34_AODSIM_v0p34/230703_184235/0000/", startFile, nFiles); // 46 files
    process = "Signal_ScenC_v0p34";
  }

  std::cout << "################################## \n";
  std::cout << "Number of files to process: " << files.size() << "\n";
  std::cout << "################################## \n";
  run3ScoutingLooper(files, year, process, outdir, "_"+std::to_string(startFile)+"To"+std::to_string(startFile+nFiles-1));

  return 0;
}

