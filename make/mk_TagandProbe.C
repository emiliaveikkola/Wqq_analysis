#include "TagandProbe.h"
#include <fstream>
R__LOAD_LIBRARY(TagandProbe_C.so);

void mk_TagandProbe(){
  TChain *c = new TChain("tree");
  string filename;
  ifstream fin("input_files/mcFiles_MuoRun2_Mikael_1718.txt");
  //input_files/mcFiles_MuoRun2_Mikael.txt
  //input_files/dataFiles_MuoRun2_Mikael.txt
  //input_files/dataFiles_MuoRun2.txt
  //input_files/mcFiles_MuoRun2.txt
  while (fin >> filename) { c->AddFile(filename.c_str()); }
  TagandProbe s(c);
  s.Loop();
}