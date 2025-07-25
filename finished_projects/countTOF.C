#include <stdio.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include "TChain.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"

#include "AddTrees.h"

void countTOF() {
  TChain chain("twotauchain");
  const Char_t* base_dir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
  AddTrees(chain, base_dir);

  Long64_t total = chain.GetEntries();
  Long64_t nEntries = std::min(total, static_cast<Long64_t>(1e9));
  printf("Chain has %lld entries (processing %lld)\n", total, nEntries);

  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("fRunNumber", 1);
  chain.SetBranchStatus("fTrkTOFexpMom", 1);
  chain.SetBranchStatus("fTrkTPCinnerParam", 1);

  Int_t   runNumber[2];    
  Float_t TOFexpMom[2];    
  Float_t innerParam[2];   

  chain.SetBranchAddress("fRunNumber",        runNumber);
  chain.SetBranchAddress("fTrkTOFexpMom",     TOFexpMom);
  chain.SetBranchAddress("fTrkTPCinnerParam", innerParam);

  std::map<Int_t, Long64_t> withTOF;
  std::map<Int_t, Long64_t> withoutTOF;
  std::map<Int_t, Long64_t> withP03;

  for (Long64_t i = 0; i < nEntries; ++i) {
    chain.GetEntry(i);

    Int_t rn = runNumber[0]; 
    bool hasTOF = false;
    bool hasP03 = false;

    
    for (int j = 0; j < 2; ++j) {
      if (TOFexpMom[j] >= 0) {
        hasTOF = true;
        break;
      }
    }

    
    for (int j = 0; j < 2; ++j) {
      if (innerParam[j] >= 0.3) { 
        hasP03 = true;
        break;
      }
    }

    if (hasTOF) ++withTOF[rn];
    else        ++withoutTOF[rn];

    if (hasP03) ++withP03[rn];
  }

  std::set<Int_t> allRuns;
  for (auto& p : withTOF)    allRuns.insert(p.first);
  for (auto& p : withoutTOF) allRuns.insert(p.first);
  for (auto& p : withP03)    allRuns.insert(p.first);


  std::ofstream jsonFile("run_counts.json");
  jsonFile << std::fixed << std::setprecision(2);
  jsonFile << "{\n";
  bool first = true;
  for (auto rn : allRuns) {
    if (!first) jsonFile << ",\n";
    first = false;
    Long64_t totalRunEntries = withTOF[rn] + withoutTOF[rn];
    jsonFile << "  \"" << rn << "\": {\n"
            << "    \"withTOF\": " << (100.0 * withTOF[rn] / totalRunEntries) << ",\n"
            << "    \"withoutTOF\": " << (100.0 * withoutTOF[rn] / totalRunEntries) << ",\n"
            << "    \"withP03\": " << (100.0 * withP03[rn] / totalRunEntries) << "\n"
            << "  }";
  }
  jsonFile << "\n}\n";
  jsonFile.close();
  printf("\nSaved run statistics to run_counts.json\n");
}