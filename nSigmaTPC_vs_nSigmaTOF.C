#include <algorithm>

#include "TChain.h"
#include "TH2F.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"  
#include "TROOT.h"
#include "TMarker.h"
#include "TLegend.h"

#include <AddTrees.h>
#include <helpers.h>

void nSigmaTPC_vs_nSigmaTOF() {
    auto help = new helper();
    const Int_t nParts = helper::nParts;
    const Int_t NtrkMax = help->NtrkMax;
    const Bool_t manualPredictPeaks = false; // Set to true if you want to manually enter peak positions
    const Int_t   nBins          = 200;
    const Double_t xMin          = -30, xMax = 30;
    const Double_t yMin          = -30, yMax = 30;
    const Double_t pStart = 0.4, pEnd = 0.7, step = 0.1;
    const Double_t nEntriesLimit  = 1e7;

    gStyle->SetPalette(kRainBow);
    gROOT->SetBatch(!manualPredictPeaks);
    TChain chain("twotauchain");
    AddTrees(chain, help->base_dir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);
    chain.SetBranchStatus("fTrkTOFexpMom", 1);

    for (Int_t i = 0; i < nParts; ++i){
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", help->pNames[i]), 1);
        chain.SetBranchStatus(Form("fTrkTOFnSigma%s", help->pNames[i]), 1);
    }
    std::vector<Float_t> inner(NtrkMax);
    std::vector<Float_t> tofExpMom(NtrkMax);
    std::vector<std::vector<Float_t>> tpcNS(nParts, std::vector<Float_t>(NtrkMax));
    std::vector<std::vector<Float_t>> tofNS(nParts, std::vector<Float_t>(NtrkMax));

    chain.SetBranchAddress("fTrkTPCinnerParam", inner.data());
    chain.SetBranchAddress("fTrkTOFexpMom", tofExpMom.data());

    for (Int_t i = 0; i < nParts; ++i) {
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", help->pNames[i]), tpcNS[i].data());
        chain.SetBranchAddress(Form("fTrkTOFnSigma%s", help->pNames[i]), tofNS[i].data());
    }
    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(nEntriesLimit));

    struct Peak { Double_t mTPC, mTOF, sTPC, sTOF; };
    for (Int_t i = 0; i < nParts; ++i) {
        TCanvas* c = new TCanvas(Form("c_%s", help->pNames[i]), Form("n#sigma %s", help->pNames[i]), 800, 600);
        c->SetLeftMargin(0.125);
        gStyle->SetOptStat(0);
        c->Print(Form("nSigma2D_%s.pdf[", help->pNames[i]));  

        const Int_t nSteps = static_cast<Int_t>(std::floor((pEnd - pStart) / step + 0.5));
        for (int iStep = 0; iStep < nSteps; ++iStep) {
            Double_t pMin = pStart + iStep * step;
            Double_t pMax = std::min(pMin + step, pEnd);
            Double_t pMid = 0.5 * (pMin + pMax);

            TH2F* h2 = new TH2F(
                Form("h2_%s_%.1f-%.1f", help->pNames[i], pMin, pMax),
                Form("n#sigma_{TPC} vs n#sigma_{TOF} (%s), %.2f < p < %.2f;n#sigma_{TPC};n#sigma_{TOF}",
                     help->pNames[i], pMin, pMax),
                nBins, xMin, xMax,
                nBins, yMin, yMax);
            h2->Sumw2();

            for (Long64_t ev = 0; ev < nEntries; ++ev) {
                chain.GetEntry(ev);
                for (Int_t t = 0; t < NtrkMax; ++t) {
                    if (tofExpMom[t] < pMin || tofExpMom[t] >= pMax) continue;
                    Double_t ntpc = tpcNS[i][t];
                    Double_t ntof = tofNS[i][t];
                    if (!TMath::IsNaN(ntpc) && !TMath::IsNaN(ntof)) {
                        h2->Fill(ntpc, ntof);
                    }
                }
            }


            c->Clear();
            h2->Draw("COLZ");
            c->SetLogz();

            std::vector<Peak> peaks;
            if (manualPredictPeaks) {
                c->Modified(); c->Update(); c->RaiseWindow();
                Int_t nManual;
                std::cout << "How many peaks? ";
                std::cin  >> nManual;

                for (Int_t p = 0; p < nManual; ++p) {
                    Peak pk;
                    std::cout << "Enter mean x y for peak " << (p+1) << ": ";
                    std::cin  >> pk.mTPC >> pk.mTOF;
                    pk.sTPC = pk.sTOF = 2.0;
                    peaks.push_back(pk);
                }
            }
            else {
                Double_t dRef = help->getTPCSignal(pMid*1000, help->pMasses[i], 1.0);
                Double_t resoRefTPC = help->getReso(helper::kTPC, help->pNames[i], pMid);
                Double_t resoRefTOF = help->getReso(helper::kTOF, help->pNames[i], pMid);
                Double_t bRef = pMid/TMath::Sqrt(pMid*pMid + help->pMasses[i]*help->pMasses[i]);

                for (Int_t hyp = 0; hyp < nParts; ++hyp) {
                    Double_t dHyp = help->getTPCSignal(pMid*1000, help->pMasses[hyp], 1.0);
                    Double_t muTPC = (dHyp/dRef - 1.0)/(resoRefTPC/dRef);

                    Double_t bHyp = pMid/TMath::Sqrt(pMid*pMid + help->pMasses[hyp]*help->pMasses[hyp]);
                    Double_t muTOF = (bRef - bHyp)/(bHyp*bHyp*help->resoTOF[i]);

                    Double_t resoHypTPC = help->getReso(helper::kTPC, help->pNames[hyp], pMid);
                    Double_t sigma0TPC  = (resoHypTPC/resoRefTPC)*(dHyp/dRef);
                    sigma0TPC = std::clamp(sigma0TPC, 0.5, 15.0);

                    Double_t resoHypTOF = help->getReso(helper::kTOF, help->pNames[hyp], pMid);
                    Double_t sigma0TOF  = (help->resoTOF[hyp]/help->resoTOF[i])*(1.0/(bHyp*bHyp));
                    sigma0TOF = std::clamp(sigma0TOF, 0.5, 15.0);

                    if (muTPC>=xMin && muTPC<=xMax && muTOF>=yMin && muTOF<=yMax)
                        peaks.push_back(Peak{muTPC, muTOF, sigma0TPC, sigma0TOF});
                }
            }

            TLegend* leg = new TLegend(0, 0.10, 0.15, 0.30);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            std::ostringstream expr;
            for(Int_t i=0; i < peaks.size(); ++i) {
            if(i>0) expr << " + ";
            Int_t b = 5*i;
            expr <<
                Form("([%d]/(2*TMath::Pi()*[%d]*[%d]))*exp(-0.5*((x-[%d])*(x-[%d])/[ %d ]/[ %d ] +"
                    "(y-[%d])*(y-[%d])/[ %d ]/[ %d ]))",
                    b+0, b+2, b+4,
                    b+1, b+1, b+2, b+2,
                    b+3, b+3, b+4, b+4
                    );
            }
            std::string formula = expr.str();

            TF2 *fsum = new TF2("fsum",
                    formula.c_str(),
                    xMin, xMax,
                    yMin, yMax);

            for(Int_t i=0; i < peaks.size(); ++i) {
            Int_t b = 5*i;
            const auto &pk = peaks[i];
            TMarker* mk = new TMarker(pk.mTPC, pk.mTOF, 20);
            mk->SetMarkerColor(help->colors[i]);   // same color[i] you use for lines
            mk->SetMarkerSize(1);
            mk->SetMarkerStyle(kMultiply);
            mk->Draw();


            Int_t bin = h2->FindBin(peaks[i].mTPC, peaks[i].mTOF);
            Double_t A0 = h2->GetBinContent(bin);

            fsum->SetParameter(b+0, A0);
            fsum->SetParameter(b+1, peaks[i].mTPC);
            fsum->SetParameter(b+2, peaks[i].sTPC);
            fsum->SetParameter(b+3, peaks[i].mTOF);
            fsum->SetParameter(b+4, peaks[i].sTOF);

            fsum->SetParLimits(b+2, peaks[i].sTPC*0.5, peaks[i].sTPC*1.5);
            fsum->SetParLimits(b+4, peaks[i].sTOF*0.5, peaks[i].sTOF*1.5);

            fsum->SetParLimits(b+1, peaks[i].mTPC - 1, peaks[i].mTPC + 1);
            fsum->SetParLimits(b+3, peaks[i].mTOF - 1, peaks[i].mTOF + 1);
            }

            h2->Fit(fsum, "QR0");  
            fsum->SetLineColor(kRed);
            fsum->SetNpx(200);
            fsum->SetNpy(200);

            for(int i=0; i < peaks.size(); ++i) {
                Int_t b = 5*i;
                Double_t A  = fsum->GetParameter(b+0);
                Double_t muX= fsum->GetParameter(b+1);
                Double_t sX = fsum->GetParameter(b+2);
                Double_t muY= fsum->GetParameter(b+3);
                Double_t sY = fsum->GetParameter(b+4);

                TF2 *f1 = new TF2(Form("f1_%d", i),
                    "[0]/(2*TMath::Pi()*[2]*[4]) *"
                    "exp(-0.5*((x-[1])*(x-[1])/[2]/[2] +"
                                "(y-[3])*(y-[3])/[4]/[4]))",
                    xMin, xMax, yMin, yMax);

                f1->SetParameters(A, muX, sX, muY, sY);
                f1->SetLineColor(help->colors[i]);
                f1->SetLineWidth(1);
                f1->SetNpx(200);
                f1->SetNpy(200);
                f1->SetContour(6);
                f1->Draw("SAME CONT3");
                leg->AddEntry(f1, help->pNames[i], "l");
            }
            leg->Draw();
            c->Print(Form("nSigma2D_%s.pdf", help->pNames[i]));
            delete h2;
        }

        c->Print(Form("nSigma2D_%s.pdf]", help->pNames[i]));
        delete c;
    }
}
