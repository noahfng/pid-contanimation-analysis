// countTOF.C
#include <stdio.h>
#include <map>
#include <set>
#include <vector>

#include "TChain.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

void countTOF() {
  // ─── 0) Adjust this to your data directory ───────────────────
  const TString base_dir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";

  // ─── 1) Build the TChain ─────────────────────────────────────
  TChain chain("twotauchain");
  {
    TSystemDirectory top("top", base_dir);
    TList *subdirs = top.GetListOfFiles();
    for (int i = 0; i < subdirs->GetEntries(); ++i) {
      TSystemFile *sf = (TSystemFile*)subdirs->At(i);
      TString dname = sf->GetName();
      if (!sf->IsDirectory() || !dname.BeginsWith("hy_")) continue;

      TString subpath = base_dir + "/" + dname;
      TSystemDirectory sd("sd", subpath);
      TList *files = sd.GetListOfFiles();
      for (int j = 0; j < files->GetEntries(); ++j) {
        TSystemFile *f = (TSystemFile*)files->At(j);
        TString fname = f->GetName();
        if (f->IsDirectory() ||
            !fname.BeginsWith("RLAnalysisTree") ||
            !fname.EndsWith(".root")) continue;

        TString full = subpath + "/" + fname;
        TFile tf(full, "READ");
        TList *keys = tf.GetListOfKeys();
        for (int k = 0; k < keys->GetEntries(); ++k) {
          TKey *key = (TKey*)keys->At(k);
          TString kn = key->GetName();
          if (kn.BeginsWith("DF_")) {
            chain.Add(full + "/" + kn + "/O2tautwotrack");
            break;
          }
        }
        tf.Close();
      }
    }
  }

  // ─── 2) Enable only the branches we need ─────────────────────
  chain.SetBranchStatus("*",            0);
  chain.SetBranchStatus("fRunNumber",    1);
  chain.SetBranchStatus("fTrkTOFexpMom", 1);

  // ─── 3) Determine the total entries and reset to tree 0 ─────
  Long64_t nEntries = chain.GetEntries();
  if (nEntries <= 0) {
    printf("No entries in chain!\n");
    return;
  }
  chain.LoadTree(0);  // reset to first tree before TTreeReader

  // ─── 4) Set up TTreeReader for safe access ───────────────────
  TTreeReader              reader(&chain);
  TTreeReaderValue<Int_t>  run   (reader, "fRunNumber");
  TTreeReaderArray<Float_t> tofArr(reader, "fTrkTOFexpMom");

  // ─── 5) Loop over *all* entries ─────────────────────────────
  std::map<Int_t,Long64_t> withTOF, withoutTOF;
  Long64_t count = 0;
  while (reader.Next()) {
    Int_t rn = *run;
    bool has = false;
    for (auto v : tofArr) {
      if (v >= 0) { has = true; break; }
    }
    if (has)       ++withTOF   [rn];
    else           ++withoutTOF[rn];
    ++count;
  }

  // ─── 6) Print summary table ─────────────────────────────────
  std::set<Int_t> allRuns;
  for (auto &p : withTOF)    allRuns.insert(p.first);
  for (auto &p : withoutTOF) allRuns.insert(p.first);

  printf("%10s | %12s | %14s\n",
         "RunNumber", "With TOF data", "Without TOF data");
  printf("---------------------------------------------\n");
  for (auto rn : allRuns) {
    printf("%10d | %12lld | %14lld\n",
           rn,
           withTOF   [rn],
           withoutTOF[rn]);
  }
}


