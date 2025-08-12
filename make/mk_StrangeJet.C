#include "../interface/StrangeJet.h"
#include <fstream>
R__LOAD_LIBRARY(src/StrangeJet_C.so);

void mk_StrangeJet(){
  TChain *c = new TChain("Events");
  string filename;
  ifstream fin("input_files/mcFiles_local_emilia.txt");
  while (fin >> filename) { c->AddFile(filename.c_str()); }
  StrangeJet s(c);
  s.Loop();
}
