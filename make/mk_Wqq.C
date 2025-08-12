#include "../interface/Wqq.h"
#include <fstream>
R__LOAD_LIBRARY(src/Wqq_C.so);

void mk_Wqq(){
  TChain *c = new TChain("tree");
  string filename;
  ifstream fin("input_files/mcFiles_MuoRun2_Mikael.txt");
  //input_files/dataFiles_MuoRun2.txt
  //input_files/mcFiles_MuoRun2.txt
  //input_files/mcFiles_stlocal_emilia.txt
  //input_files/mcFiles_stlocal_emilia_udsample.txt
  //input_files/mcFiles_stlocal_emilia_cssample.txt
  while (fin >> filename) { c->AddFile(filename.c_str()); }
  Wqq s(c);
  s.Loop();
}