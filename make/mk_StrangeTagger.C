#include "StrangeTagger.h"
#include <fstream>
R__LOAD_LIBRARY(StrangeTagger_C.so);

void mk_StrangeTagger(){
  TChain *c = new TChain("tree");
  string filename;
  ifstream fin("input_files/mcFiles_MuoRun2.txt");
  //input_files/dataFiles_MuoRun2.txt
  //input_files/mcFiles_MuoRun2.txt
  //input_files/mcFiles_stlocal_emilia.txt
  //input_files/mcFiles_stlocal_emilia_udsample.txt
  //input_files/mcFiles_stlocal_emilia_cssample.txt
  while (fin >> filename) { c->AddFile(filename.c_str()); }
  StrangeTagger s(c);
  s.Loop();
}