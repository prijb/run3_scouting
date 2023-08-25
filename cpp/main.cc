#include <cstdlib>
#include <iostream>
#include <fstream>
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


std::vector<TString> getFiles(const std::string inputDir, const int startFile, const int nFiles, const bool isCondor) {
  std::vector<TString> files;
  std::string fullInputDir;
  unsigned int iFile=0;
  if (!isCondor) {
    if (inputDir.rfind("/ceph/cms", 0) == 0)
      fullInputDir = inputDir;
    else
      fullInputDir = "/ceph/cms"+inputDir;
    const fs::path dir{fullInputDir};
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
  }
  else {
    std::string temp_str(inputDir);
    if (inputDir.rfind("/ceph/cms", 0) == 0)
      fullInputDir = temp_str.replace(temp_str.find("/ceph/cms"),sizeof("/ceph/cms")-1,"");
    else
      fullInputDir = inputDir;
    std::string command;
    command = "xrdfs redirector.t2.ucsd.edu:1095 ls ";
    command += fullInputDir;
    command += " > infiles.txt";
    std::system(command.c_str());
    std::ifstream infiles("infiles.txt");
    std::string line;
    while(getline(infiles, line)) {
      if (iFile<startFile) {
	iFile++;
	continue;
      }
      if (!TString(line.c_str()).Contains(".root"))
	continue;
      else {
	files.push_back(TString("davs://redirector.t2.ucsd.edu:1095"+line));
      }
      iFile++;
      if (iFile == startFile+nFiles)
	break;
    }
  }
  std::system("rm -f infiles.txt");
  return files;
}

int main(int argc, char **argv) {
  // Arguments
  const char* outdir      = ( argc > 1  ? argv[1]            : "temp_data" );
  TString year            = ( argc > 2  ? argv[2]            : "2022" );
  TString sampleArg       = ( argc > 3  ? argv[3]            : "DataF" );
  int startFile           = ( argc > 4  ? char2int(argv[4])  : 0 );
  int nFiles              = ( argc > 5  ? char2int(argv[5])  : 1000000 ); // Large number as default to run over all files
  bool isCondor           = ( argc > 6  ? bool(argv[6])      : false );
  std::vector<TString> files;
  TString process;
  // Sample list: Data
  if ( sampleArg=="DataB" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Data/2022B/", startFile, nFiles, isCondor);
    process = "DataB";
  }
  if ( sampleArg=="DataC" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Data/2022C/", startFile, nFiles, isCondor);
    process = "DataC";
  }
  if ( sampleArg=="DataD" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Data/2022D/", startFile, nFiles, isCondor);
    process = "DataD";
  }
  if ( sampleArg=="DataE" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Data/2022E/", startFile, nFiles, isCondor);
    process = "DataE";
  }
  if ( sampleArg=="DataF" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Data/2022F/", startFile, nFiles, isCondor);
    process = "DataF";
  }
  if ( sampleArg=="DataG" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Data/2022G/", startFile, nFiles, isCondor);
    process = "DataG";
  }
  //
  // Sample list: Signal
  if ( sampleArg=="Signal_ScenA_v0p30" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Signal/Signal_ScenA_v0p30/", startFile, nFiles, isCondor); // 61 files
    process = "Signal_ScenA_v0p30";
  }
  if ( sampleArg=="Signal_ScenA_v1p3" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Signal/Signal_ScenA_v1p3/", startFile, nFiles, isCondor); // 69 files
    process = "Signal_ScenA_v1p3";
  }
  if ( sampleArg=="Signal_ScenA_v1p4" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Signal/Signal_ScenA_v1p4/", startFile, nFiles, isCondor); // 95 files
    process = "Signal_ScenA_v1p4";
  }
  if ( sampleArg=="Signal_ScenA_v1p5" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Signal/Signal_ScenA_v1p5", startFile, nFiles, isCondor); // 93 files
    process = "Signal_ScenA_v1p5";
  }
  if ( sampleArg=="Signal_ScenB_v0p32" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Signal/Signal_ScenB_v0p32/", startFile, nFiles, isCondor); // 48 files
    process = "Signal_ScenB_v0p32";
  }
  if ( sampleArg=="Signal_ScenC_v0p34" ) {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jul-29-2023/Signal/Signal_ScenC_v0p34/", startFile, nFiles, isCondor); // 46 files
    process = "Signal_ScenC_v0p34";
  }

  std::cout << "################################## \n";
  std::cout << "Number of files to process: " << files.size() << "\n";
  std::cout << "################################## \n";
  run3ScoutingLooper(files, year, process, outdir, "_"+std::to_string(startFile)+"To"+std::to_string(startFile+nFiles-1));

  return 0;
}
