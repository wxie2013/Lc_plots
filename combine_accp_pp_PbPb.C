#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TCut.h>
#include <TLatex.h>

void combine_accp_pp_PbPb(){
TFile *f_pp = TFile::Open("/depot/cms/users/chand140/Final_results_05172022/combine_accp_pp.root");
TFile *f_PbPb_0_90 = TFile::Open("/depot/cms/users/wxie/xiao147/private/2018_PbPb_Lc_rereco_Dntuple/CMSSW_10_3_3_patch1/src/Bfinder/Bfinder/Dntuple/analysis_steps/accp/ROOT/several_turns/third_turn_PVz/combine_accp/accp_P_cen0_90_thirdturn_PVz.root");
TFile *f_PbPb_0_10 = TFile::Open("/depot/cms/users/wxie/xiao147/private/2018_PbPb_Lc_rereco_Dntuple/CMSSW_10_3_3_patch1/src/Bfinder/Bfinder/Dntuple/analysis_steps/accp/ROOT/several_turns/third_turn_PVz/combine_accp/accp_P_cen0_10_thirdturn_PVz.root");
TFile *f_PbPb_10_30 = TFile::Open("/depot/cms/users/wxie/xiao147/private/2018_PbPb_Lc_rereco_Dntuple/CMSSW_10_3_3_patch1/src/Bfinder/Bfinder/Dntuple/analysis_steps/accp/ROOT/several_turns/third_turn_PVz/combine_accp/accp_P_cen10_30_thirdturn.root");
TFile *f_PbPb_30_50 = TFile::Open("/depot/cms/users/wxie/xiao147/private/2018_PbPb_Lc_rereco_Dntuple/CMSSW_10_3_3_patch1/src/Bfinder/Bfinder/Dntuple/analysis_steps/accp/ROOT/several_turns/third_turn_PVz/combine_accp/accp_P_cen30_50_thirdturn.root");
TFile *f_PbPb_50_90 = TFile::Open("/depot/cms/users/wxie/xiao147/private/2018_PbPb_Lc_rereco_Dntuple/CMSSW_10_3_3_patch1/src/Bfinder/Bfinder/Dntuple/analysis_steps/accp/ROOT/several_turns/third_turn_PVz/combine_accp/accp_P_cen50_90_thirdturn.root");


TH1F *h_PbPb_0_90 = (TH1F*) f_PbPb_0_90->Get("hrecoeff_combine")->Clone("h_PbPb_0_90");
TH1F *h_PbPb_0_10 = (TH1F*) f_PbPb_0_10->Get("hrecoeff_combine")->Clone("h_PbPb_0_10");
TH1F *h_PbPb_10_30 = (TH1F*) f_PbPb_10_30->Get("hrecoeff_combine")->Clone("h_PbPb_10_30");
TH1F *h_PbPb_30_50 = (TH1F*) f_PbPb_30_50->Get("hrecoeff_combine")->Clone("h_PbPb_30_50");
TH1F *h_PbPb_50_90 = (TH1F*) f_PbPb_50_90->Get("hrecoeff_combine")->Clone("h_PbPb_50_90");
TH1F *h_pp = (TH1F*) f_pp->Get("hrecoeff_combine")->Clone("h_pp");

TCanvas *c1 = new TCanvas("c1","multipads",600,600);
c1->SetLeftMargin(0.12);
c1->SetRightMargin(0.05);
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);
const Int_t NBINS = 10;
Double_t edges[NBINS + 1]={3,4,5,6,8,10,12.5,15,20,30,40};
TH1 *h_empty = new TH1D("h_empty","h_empty",NBINS,edges);
h_empty->Sumw2();
h_empty->Draw();
h_pp->Draw("same");
h_PbPb_0_90->Draw("same");
h_PbPb_0_10->Draw("same");
h_PbPb_10_30->Draw("same");
h_PbPb_30_50->Draw("same");
h_PbPb_50_90->Draw("same");
h_pp->SetMarkerColor(1);
h_pp->SetMarkerStyle(20);
h_pp->SetMarkerSize(1.2);
h_pp->SetLineColor(1);

h_PbPb_0_90->SetMarkerStyle(29);
h_PbPb_0_90->SetMarkerColor(2);
h_PbPb_0_90->SetLineColor(2);
h_PbPb_0_90->SetMarkerSize(1.5);

h_PbPb_0_10->SetMarkerStyle(21);
h_PbPb_0_10->SetMarkerColor(kAzure+1);
h_PbPb_0_10->SetMarkerSize(1.2);
h_PbPb_0_10->SetLineColor(kAzure+1);

h_PbPb_10_30->SetMarkerStyle(22);
h_PbPb_10_30->SetMarkerColor(kOrange+0);
h_PbPb_10_30->SetMarkerSize(1.5);
h_PbPb_10_30->SetLineColor(kOrange+0);

h_PbPb_30_50->SetMarkerStyle(23);
h_PbPb_30_50->SetMarkerColor(kGreen+1);
h_PbPb_30_50->SetMarkerSize(1.5);
h_PbPb_30_50->SetLineColor(kGreen+1);

h_PbPb_50_90->SetMarkerStyle(33);
h_PbPb_50_90->SetMarkerColor(kPink+9);
h_PbPb_50_90->SetMarkerSize(1.8);
h_PbPb_50_90->SetLineColor(kPink+9);

h_empty->GetXaxis()->CenterTitle();
h_empty->GetYaxis()->CenterTitle();
h_empty->SetAxisRange(0,0.4,"Y");
h_empty->SetXTitle("p_{T} (GeV/c)");
h_empty->SetYTitle("#alpha #times #epsilon");
h_empty->SetTitleOffset(1.4,"Y");
h_empty->SetTitleOffset(1.0,"X");
h_empty->GetXaxis()->SetTitleSize(0.045);
h_empty->GetYaxis()->SetTitleSize(0.045);
h_empty->GetXaxis()->SetTitleFont(42);
h_empty->GetYaxis()->SetTitleFont(42);
h_empty->GetXaxis()->SetLabelFont(42);
h_empty->GetYaxis()->SetLabelFont(42);
h_empty->GetXaxis()->SetLabelSize(0.04);
h_empty->GetYaxis()->SetLabelSize(0.04);

auto leg0 = new TLegend(0.15,0.75,0.4, 0.88);
leg0->SetTextSize(0.035);
leg0->AddEntry(h_PbPb_0_90,"PbPb: 0-90%","lep");
leg0->AddEntry(h_PbPb_0_10,"PbPb: 0-10%","lep");
leg0->AddEntry(h_PbPb_10_30,"PbPb: 10-30%","lep");
leg0->SetBorderSize(0);
leg0->Draw();

auto leg1 = new TLegend(0.45,0.75,0.88,0.88);
leg1->SetTextSize(0.035);
leg1->AddEntry(h_PbPb_30_50,"PbPb: 30-50%","lep");
leg1->AddEntry(h_PbPb_50_90,"PbPb: 50-90%","lep");
leg1->AddEntry(h_pp,"pp","lep");
leg1->SetBorderSize(0);
leg1->Draw();

TLatex Tl;
Tl.SetNDC();
Tl.SetTextAlign(12);
Tl.SetTextSize(0.05);
Tl.SetTextFont(42);
Tl.DrawLatex(0.12,0.93, "#font[61]{CMS }");
Tl.DrawLatex(0.81,0.93, "#scale[0.8]{5.02 TeV}");

//TLatex* tex;
//tex = new TLatex(0.16,0.79,"PbPb: HYDJET");
//tex->SetNDC();
//tex->SetTextFont(42);
//tex->SetTextSize(0.05*0.75);
//tex->SetLineWidth(2);
//tex->Draw();
//
//tex = new TLatex(0.16,0.74,"pp: PYTHIA8");
//tex->SetNDC();
//tex->SetTextFont(42);
//tex->SetTextSize(0.05*0.75);
//tex->SetLineWidth(2);
//tex->Draw();

TLatex Tl2;
Tl.SetNDC();
Tl.SetTextAlign(12);
Tl.SetTextSize(0.05*0.75);
Tl.SetTextFont(42);
Tl.DrawLatex(0.23,0.93, "#font[52]{Simulations}");
//Tl.DrawLatex(0.16,0.86, "#font[52]{Work in progress}");

c1->SaveAs("plots/accp_combine_2017_2018_pp_PbPb_withPVz_reweight_deposit_updated.pdf");

}
