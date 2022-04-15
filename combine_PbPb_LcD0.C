#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TF1.h>
#include <TNamed.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TGraphErrors.h>
void combine_PbPb_LcD0(){

  TFile *f_pp_crosssection = new TFile("/depot/cms/users/wxie/xiao147/private/2017_pp_Lc_speedup_version_2_toDtuple/CMSSW_9_4_10/src/Bfinder/Bfinder/Dntuple/analysis_firstturn_04072020/all_results/ROOT/manually_weight_chang1251520_055_working077/pp_LcD0_withsys.root");
  TFile *f_PbPb_yield = new TFile("/depot/cms/users/wxie/xiao147/private/2018_PbPb_Lc_rereco_Dntuple/CMSSW_10_3_3_patch1/src/Bfinder/Bfinder/Dntuple/analysis_steps/results/ROOT/manually_pp_reweight/LcD0_090_010_withsys.root");
  TFile *f_greco_PbPb_0_100 = new TFile("/depot/cms/users/wxie/xiao147/private/2018_PbPb_Lc_rereco_Dntuple/CMSSW_10_3_3_patch1/src/Bfinder/Bfinder/Dntuple/analysis_steps/theory_predictions/ROOT/Greco_PbPb_0_100.root");
  TFile *f_greco_PbPb_0_100_minijet = new TFile("/depot/cms/users/wxie/xiao147/private/2018_PbPb_Lc_rereco_Dntuple/CMSSW_10_3_3_patch1/src/Bfinder/Bfinder/Dntuple/analysis_steps/theory_predictions/ROOT/Greco_PbPb_0_100_minijet.root");

  TFile *f_ralf_020_PbPb_high = new TFile("/depot/cms/users/wxie/xiao147/private/2018_PbPb_Lc_rereco_Dntuple/CMSSW_10_3_3_patch1/src/Bfinder/Bfinder/Dntuple/analysis_steps/theory_predictions/ROOT/MinHe_Rapp_PbPb_cen0_20_high.root");
  TFile *f_ralf_020_PbPb_low = new TFile("/depot/cms/users/wxie/xiao147/private/2018_PbPb_Lc_rereco_Dntuple/CMSSW_10_3_3_patch1/src/Bfinder/Bfinder/Dntuple/analysis_steps/theory_predictions/ROOT/MinHe_Rapp_PbPb_cen0_20_low.root");



  TGraphAsymmErrors *pp_crosssection_stat = (TGraphAsymmErrors*) f_pp_crosssection->Get("pp_crosssection_stat");
  //add x-error to be consistent with other plots
  const Int_t NBINS_1 = 9;
  double x_err[NBINS_1] ={0.5,0.5,0.5,1,1,1.25,1.25,2.5,5};
  for(int i = 0; i<pp_crosssection_stat->GetN(); i++) {
    pp_crosssection_stat->SetPointEXhigh(i, x_err[i]);
    pp_crosssection_stat->SetPointEXlow(i, x_err[i]);
  }

  TGraphErrors *pp_crosssection_sys = (TGraphErrors*)f_pp_crosssection->Get("pp_crosssection_sys");
  TGraphAsymmErrors *stat_0_90 = (TGraphAsymmErrors*) f_PbPb_yield->Get("stat_0_90");
  TGraphAsymmErrors *stat_0_10 = (TGraphAsymmErrors*) f_PbPb_yield->Get("stat_0_10");
  TGraphErrors *sys_0_90 = (TGraphErrors*)f_PbPb_yield->Get("sys_0_90");
  TGraphErrors *sys_0_10 = (TGraphErrors*)f_PbPb_yield->Get("sys_0_10");

  TGraphErrors *greco_PbPb = (TGraphErrors*)f_greco_PbPb_0_100->Get("MyGraph");
  TGraphErrors *greco_PbPb_minijets = (TGraphErrors*)f_greco_PbPb_0_100_minijet->Get("MyGraph");

  TGraph *ralf_PbPb_high = (TGraph*) f_ralf_020_PbPb_high->Get("MyGraph");
  TGraph *ralf_PbPb_low = (TGraph*) f_ralf_020_PbPb_low->Get("MyGraph");
  TGraph *ralf_PbPb  = new TGraph(120);

  //try to get a band from Ralf pp predictions
  double x[60]={0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1,10.3,10.5,10.7,10.9,11.1,11.3,11.5,11.7,11.9};
  for (int j=0; j<60; j++)
  {
    ralf_PbPb->SetPoint(j,x[j],ralf_PbPb_high->Eval(x[j]));
    ralf_PbPb->SetPoint(60+j,x[60-j-1],ralf_PbPb_low->Eval(x[60-j-1]));

  }

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  c1->SetLeftMargin(0.155);
  c1->SetRightMargin(0.03);
  c1->SetTopMargin(0.06);
  c1->SetBottomMargin(0.12);

  const Int_t NBINS = 10+2; //+2 means leaving space for plot edge
  Double_t edges[NBINS + 1]={2,3,4,5,6,8,10,12.5,15,20,30,40, 40.5};
  TH1F *h_empty =  new TH1F("h_empty","h_empty",NBINS,edges);
  h_empty->Draw();
  h_empty->GetXaxis()->CenterTitle();
  h_empty->GetYaxis()->CenterTitle();
  h_empty->GetXaxis()->SetTitleOffset(1.0);
  h_empty->GetYaxis()->SetTitleOffset(1.6);
  h_empty->GetXaxis()->SetLabelOffset(0.007);
  h_empty->GetYaxis()->SetLabelOffset(0.007);
  h_empty->GetXaxis()->SetTitleSize(0.04);
  h_empty->GetYaxis()->SetTitleSize(0.04);
  h_empty->GetXaxis()->SetTitleFont(42);
  h_empty->GetYaxis()->SetTitleFont(42);
  h_empty->GetXaxis()->SetLabelFont(42);
  h_empty->GetYaxis()->SetLabelFont(42);
  h_empty->GetXaxis()->SetLabelSize(0.04);
  h_empty->GetYaxis()->SetLabelSize(0.04);
  h_empty->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_empty->GetYaxis()->SetTitle("(#Lambda_{c}^{+} + #Lambda_{c}^{#font[122]{\55}}) / (D^{0}+ #bar{D^{0}})");
  h_empty->SetAxisRange(0,1.2,"Y");
  h_empty->Draw("E");
  h_empty->SetMarkerColor(1);
  h_empty->SetMarkerStyle(28);
  h_empty->SetMarkerSize(2);
  h_empty->SetLineColor(1);
  //gPad->SetLogy();

  ralf_PbPb->SetFillColor(18);
  ralf_PbPb->SetFillStyle(3001);
  ralf_PbPb->SetLineColor(18);
  ralf_PbPb->Draw("Fsame");

  pp_crosssection_sys->SetMarkerColor(1);
  pp_crosssection_sys->SetMarkerStyle(20);
  pp_crosssection_sys->SetMarkerSize(1.8);
  pp_crosssection_sys->SetLineWidth(0);
  pp_crosssection_sys->SetFillStyle(1001);
  pp_crosssection_sys->SetLineColor(1);
  pp_crosssection_sys->SetFillColorAlpha(16,0.6);
  pp_crosssection_sys->GetXaxis()->CenterTitle();
  pp_crosssection_sys->GetYaxis()->CenterTitle();
  pp_crosssection_sys->GetXaxis()->SetTitleOffset(1.0);
  pp_crosssection_sys->GetYaxis()->SetTitleOffset(1.0);
  pp_crosssection_sys->GetXaxis()->SetLabelOffset(0.007);
  pp_crosssection_sys->GetYaxis()->SetLabelOffset(0.007);
  pp_crosssection_sys->GetXaxis()->SetTitleSize(0.045);
  pp_crosssection_sys->GetYaxis()->SetTitleSize(0.045);
  pp_crosssection_sys->GetXaxis()->SetTitleFont(42);
  pp_crosssection_sys->GetYaxis()->SetTitleFont(42);
  pp_crosssection_sys->GetXaxis()->SetLabelFont(42);
  pp_crosssection_sys->GetXaxis()->SetLabelFont(42);
  pp_crosssection_sys->GetXaxis()->SetLabelSize(0.04);
  pp_crosssection_sys->GetYaxis()->SetLabelSize(0.04);
  pp_crosssection_sys->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  pp_crosssection_sys->GetYaxis()->SetTitle("R_{AA}");
  pp_crosssection_sys->Draw("2");

  pp_crosssection_stat->SetMarkerStyle(20);
  pp_crosssection_stat->SetMarkerColor(1);
  pp_crosssection_stat->SetLineColor(1);
  pp_crosssection_stat->Draw("psame");
  pp_crosssection_stat->SetMarkerSize(1.8);



  sys_0_90->SetMarkerColor(2);
  sys_0_90->SetMarkerStyle(29);
  sys_0_90->SetMarkerSize(1.8);
  sys_0_90->SetLineWidth(0);
  sys_0_90->SetFillStyle(1001);
  sys_0_90->SetLineColor(2);
  sys_0_90->SetFillColorAlpha(2,0.2);
  sys_0_90->GetXaxis()->CenterTitle();
  sys_0_90->GetYaxis()->CenterTitle();
  sys_0_90->GetXaxis()->SetTitleOffset(1.0);
  sys_0_90->GetYaxis()->SetTitleOffset(1.0);
  sys_0_90->GetXaxis()->SetLabelOffset(0.007);
  sys_0_90->GetYaxis()->SetLabelOffset(0.007);
  sys_0_90->GetXaxis()->SetTitleSize(0.045);
  sys_0_90->GetYaxis()->SetTitleSize(0.045);
  sys_0_90->GetXaxis()->SetTitleFont(42);
  sys_0_90->GetYaxis()->SetTitleFont(42);
  sys_0_90->GetXaxis()->SetLabelFont(42);
  sys_0_90->GetXaxis()->SetLabelFont(42);
  sys_0_90->GetXaxis()->SetLabelSize(0.04);
  sys_0_90->GetYaxis()->SetLabelSize(0.04);
  sys_0_90->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  sys_0_90->GetYaxis()->SetTitle("R_{AA}");
  sys_0_90->Draw("2same");
  stat_0_90->SetMarkerStyle(29);
  stat_0_90->SetMarkerColor(2);
  stat_0_90->SetLineColor(2);
  stat_0_90->SetMarkerSize(1.8);

  sys_0_10->SetMarkerColor(kAzure+1);
  sys_0_10->SetMarkerStyle(21);
  sys_0_10->SetMarkerSize(1.6);
  sys_0_10->SetLineWidth(0);
  sys_0_10->SetFillStyle(1001);
  sys_0_10->SetLineColor(9);
  sys_0_10->SetFillColorAlpha(kAzure+1,0.4);
  sys_0_10->GetXaxis()->CenterTitle();
  sys_0_10->GetYaxis()->CenterTitle();
  sys_0_10->GetXaxis()->SetTitleOffset(1.0);
  sys_0_10->GetYaxis()->SetTitleOffset(1.0);
  sys_0_10->GetXaxis()->SetLabelOffset(0.007);
  sys_0_10->GetYaxis()->SetLabelOffset(0.007);
  sys_0_10->GetXaxis()->SetTitleSize(0.045);
  sys_0_10->GetYaxis()->SetTitleSize(0.045);
  sys_0_10->GetXaxis()->SetTitleFont(42);
  sys_0_10->GetYaxis()->SetTitleFont(42);
  sys_0_10->GetXaxis()->SetLabelFont(42);
  sys_0_10->GetXaxis()->SetLabelFont(42);
  sys_0_10->GetXaxis()->SetLabelSize(0.06);
  sys_0_10->GetYaxis()->SetLabelSize(0.06);
  sys_0_10->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  sys_0_10->GetYaxis()->SetTitle("R_{AA}");
  sys_0_10->Draw("2same");
  stat_0_10->SetMarkerStyle(21);
  stat_0_10->SetMarkerColor(kAzure+1);
  stat_0_10->SetLineColor(kAzure+1);
  stat_0_10->SetMarkerSize(1.6);
  stat_0_10->Draw("psame");
  stat_0_90->Draw("psame");

  //greco_PbPb->SetLineWidth(3);
  //greco_PbPb->SetLineColor(1);
  //greco_PbPb->SetLineStyle(9);
  //greco_PbPb->Draw("same");
  //greco_PbPb_minijets->SetLineWidth(3);
  //greco_PbPb_minijets->SetLineColor(4);

  TLatex Tl;
  Tl.SetNDC();
  Tl.SetTextAlign(12);
  Tl.SetTextSize(0.05);
  Tl.SetTextFont(42);
  Tl.DrawLatex(0.18,0.91, "#font[61]{CMS}");
  Tl.DrawLatex(0.56,0.97, "#scale[0.8]{PbPb 0.58 nb^{-1} (5.02 TeV)}");

  TLatex Tl2;
  Tl.SetNDC();
  Tl.SetTextAlign(12);
  Tl.SetTextSize(0.05*0.75);
  Tl.SetTextFont(42);
  Tl.DrawLatex(0.29,0.91, "#font[52]{Preliminary}");

  auto leg = new TLegend(0.49,0.70,0.92,0.90);
  leg->AddEntry(pp_crosssection_sys,"pp","fpe");
  leg->AddEntry(sys_0_90,"0-90% PbPb","fp");
  leg->AddEntry(sys_0_10,"0-10% PbPb","fp");
  //leg->AddEntry(greco_PbPb,"EPJC78 (2018) 348","l");
  leg->AddEntry(ralf_PbPb,"PRL124 (2020) 042301","f");
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->Draw();
  TLegendEntry *header1 = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header1->SetTextSize(0.045);

  TLatex* tex4;
  tex4 = new TLatex(0.5,0.65,"Global uncertainty");
  tex4->SetNDC();
  tex4->SetTextFont(42);
  tex4->SetTextSize(0.04);
  tex4->SetLineWidth(2);
  tex4->Draw();

  TLatex* tex5;
  tex5 = new TLatex(0.5,0.61,"pp: 6.8%");
  tex5->SetNDC();
  tex5->SetTextFont(42);
  tex5->SetTextSize(0.04);
  tex5->SetLineWidth(2);
  tex5->Draw();

  TLatex* tex6;
  tex6 = new TLatex(0.5,0.56,"PbPb: 7.4%");
  tex6->SetNDC();
  tex6->SetTextFont(42);
  tex6->SetTextSize(0.04);
  tex6->SetLineWidth(2);
  tex6->Draw();

  TLatex* tex_6;
  tex_6 = new TLatex(0.18,0.85,"|y| < 1");
  tex_6->SetNDC();
  tex_6->SetTextFont(42);
  tex_6->SetTextSize(0.04);
  tex_6->SetLineWidth(2);
  tex_6->Draw();

  c1->SaveAs("plots/LcD0_pp_PbPb.pdf");







}
