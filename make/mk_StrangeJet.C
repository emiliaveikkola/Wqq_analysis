#include "StrangeJet.h"
#include <fstream>
R__LOAD_LIBRARY(StrangeJet_C.so);

void mk_StrangeJet(){
  TChain *c = new TChain("Events");
  string filename;
  ifstream fin("input_files/mcFiles_local_emilia.txt");
  while (fin >> filename) { c->AddFile(filename.c_str()); }
  //c->MakeClass("StrangeJetNew");
  StrangeJet s(c);
  s.Loop();
}
