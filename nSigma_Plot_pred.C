#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TF1.h"  
#include "TSpectrum.h"
#include "TSystem.h"
#include "RooFit.h"

#include <AddTrees.h>

Double_t bethe_bloch_aleph(Double_t bg, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5) {
    Double_t beta = bg / TMath::Sqrt(1.0 + bg*bg);
    Double_t aa   = TMath::Power(beta, p4);
    Double_t bb   = TMath::Log(p3 + TMath::Power(1.0/bg, p5));
    return (p2 - aa - bb) * p1 / aa;
}

Double_t get_expected_signal(Double_t p, Double_t mass, Double_t charge) {
    const Double_t mMIP = 50.0;
    const Double_t params[5] = {0.19310481, 4.26696118, 0.00522579, 2.38124907, 0.98055396};
    const Double_t chFact = 2.3;

    Double_t bg = p / mass;
    if (bg < 0.05) return -999.;
    Double_t bethe = mMIP
                   * bethe_bloch_aleph(bg,
                                       params[0], params[1], params[2],
                                       params[3], params[4])
                   * TMath::Power(charge, chFact);
    return bethe >= 0 ? bethe : -999.;
}

std::vector<double> topBinCenters(TH1 *h, int nWanted)
{
   std::vector<std::pair<double,int>> bins;        
   for (int b = 1; b <= h->GetNbinsX(); ++b)
       bins.emplace_back(h->GetBinContent(b), b);

   std::partial_sort(bins.begin(), bins.begin()+nWanted, bins.end(),
                     std::greater<>());              

   std::vector<double> xc;
   for (int i = 0; i < nWanted; ++i)
       xc.push_back(h->GetXaxis()->GetBinCenter(bins[i].second));

   std::sort(xc.begin(), xc.end());                 
   return xc;
}

void nSigma_Plot_pred(double pStart = 0.1, double pEnd = 2.0, double step = 0.1, double muWindow = 1.0)   
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1);

    const char *baseDir = "/home/nfingerle/SMI/UD_LHC23_pass4_SingleGap/0106/B";
    TChain chain("twotauchain");
    AddTrees(chain, baseDir);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("fTrkTPCinnerParam", 1);

    const char *subs[5] = {"El", "Mu", "Pi", "Ka", "Pr"};
    for (int i = 0; i < 5; ++i)
        chain.SetBranchStatus(Form("fTrkTPCnSigma%s", subs[i]), 1);

    Float_t inner[2];
    Float_t tpcNS[5][2];
    chain.SetBranchAddress("fTrkTPCinnerParam", inner);
    for (int i = 0; i < 5; ++i)
        chain.SetBranchAddress(Form("fTrkTPCnSigma%s", subs[i]), tpcNS[i]);

    const Int_t   nBins   = 200;
    const Double_t xMin   = -15.0, xMax = 15.0;
    const Int_t   nParts  = 5;
    const TString names[nParts]   = {"e", "#mu", "#pi", "K", "p"};
    const Int_t   colors[nParts] = {kBlue, kBlueYellow, kGreen + 2,
                                    kOrange + 7, kViolet};
    const Double_t resoTPC[nParts] = {0.085, 0.072, 0.074, 0.09, 0.08}; 
    const Double_t masses[nParts]  = {0.00051099895, 0.1056583755,
                                      0.13957039,    0.493677,
                                      0.93827208816};

    Long64_t nEntries = std::min(chain.GetEntries(), static_cast<Long64_t>(1e6));

    for (int ref = 0; ref < nParts; ++ref) {

        TString pdfName = Form("nSigmaTPC_%s.pdf", names[ref].Data());
        TCanvas *c = new TCanvas("c","n#sigma("+names[ref]+")",950,700);
        c->SetLeftMargin(0.15); c->SetGrid(); c->SetLogy();
        c->Print(pdfName+"[");

        for (double pMin = pStart; pMin < pEnd; pMin += step) {
            double pMax = pMin + step;
            TH1F *h = new TH1F(Form("h_%s_%g",names[ref].Data(),pMin),
                Form("n#sigma_{%s}, %.1f < p < %.1f GeV/c; n#sigma; Counts",
                     names[ref].Data(), pMin, pMax),
                nBins,xMin,xMax);
            h->SetMarkerStyle(kFullCircle); h->SetMarkerSize(0.75);
            h->SetMarkerColor(kBlack); h->SetLineColor(kBlack);

            for(Long64_t ev=0; ev<nEntries; ++ev){
                chain.GetEntry(ev);
                for(int t=0;t<2;++t){
                    float pG=inner[t];
                    if(pG<pMin||pG>=pMax) continue;
                    if(!TMath::IsNaN(tpcNS[ref][t])) h->Fill(tpcNS[ref][t]);
                }
            }

            double pMid = 0.5*(pMin+pMax);
            double dRef = get_expected_signal(pMid*1000., masses[ref]*1000, 1.0);
            std::vector<double> muSeed; std::vector<int> hypId;
            double yMax = 1.05*h->GetMaximum();
            for(int hyp=0; hyp<nParts; ++hyp){
                double dHyp = get_expected_signal(pMid*1000., masses[hyp]*1000, 1.0);
                if(dRef<0||dHyp<0) continue;
                double mu = (dHyp/dRef -1.0)/resoTPC[hyp];
                if(mu<xMin||mu>xMax) continue;

                TLine *l=new TLine(mu,0,mu,yMax); l->SetLineColor(colors[hyp]); l->Draw();
                TLatex tx; tx.SetTextColor(colors[hyp]); tx.SetTextSize(0.03); tx.SetTextAlign(22);
                tx.DrawLatex(mu,0.92*yMax,names[hyp]);
                muSeed.push_back(mu); hypId.push_back(hyp);
            }
            int nG = muSeed.size();
            if(nG==0){ c->Print(pdfName); delete h; continue; }

            std::ostringstream form; for(int i=0;i<nG;++i){ if(i) form<<"+"; form<<"gaus("<<3*i<<")"; }
            TF1 *sum=new TF1("sum",form.str().c_str(),xMin,xMax);
            for(int i=0;i<nG;++i){
                int bin=h->FindBin(muSeed[i]); double amp=h->GetBinContent(bin);
                sum->SetParameters(3*i, amp, muSeed[i], 0.30);
                sum->SetParLimits(3*i+0, 0.0, h->GetMaximum()*1.2);
                sum->SetParLimits(3*i+1, muSeed[i]-muWindow, muSeed[i]+muWindow);
                sum->SetParLimits(3*i+2, 0, 100.0);
            }
            h->Fit(sum,"QR0");
            c->Clear(); h->Draw("E1");
            sum->SetLineColor(kRed); sum->SetLineWidth(2); sum->SetNpx(500); sum->Draw("same");
            for(int i=0;i<nG;++i){ if(sum->GetParameter(3*i)<=0) continue;
                TF1 *g=new TF1(Form("g_%d_%d",ref,i),"gaus",xMin,xMax);
                g->SetParameters(&sum->GetParameters()[3*i]);
                g->SetLineColor(colors[hypId[i]]); g->SetLineStyle(2); g->Draw("same"); }
            for(size_t i=0;i<muSeed.size();++i){ TLine *l=new TLine(muSeed[i],0,muSeed[i],yMax); l->SetLineColor(colors[hypId[i]]); l->Draw(); }
            TPaveText *pt=new TPaveText(0.02,0.90,0.25,0.99,"NDC");
            pt->AddText(Form("#chi^{2}/NDF = %.2f",sum->GetChisquare()/sum->GetNDF()));
            pt->AddText(Form("N_{Gauss} = %d",nG)); pt->SetFillColorAlpha(0,0); pt->Draw("same");

            c->Print(pdfName);
            delete sum; delete h;
        }
        c->Print(pdfName+"]"); delete c;
    }
}

