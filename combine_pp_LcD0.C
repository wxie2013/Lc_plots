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
void combine_pp_LcD0(){

  TFile *f_pp_crosssection = new TFile("root_file/pp_LcD0_withsys.root");
  TFile *f_CR2 = new TFile("root_file/pp_Lc_D0_ratio_data_withstatistic_CR2predictions.root");
  TFile *f_greco_low = new TFile("root_file/Plumari_pp_low.root");
  TFile *f_greco_high = new TFile("root_file/Plumari_pp_high.root");

  TFile *f_ralf_pp_high =  new TFile("root_file/MinHe_Rapp_pp_RQM_high.root");
  TFile *f_ralf_pp_low = new TFile("root_file/MinHe_Rapp_pp_RQM_low.root");

  TGraphAsymmErrors *pp_crosssection_stat = (TGraphAsymmErrors*) f_pp_crosssection->Get("pp_crosssection_stat");
  //add x-error to be consistent with other plots
  const Int_t NBINS_1 = 9;
  double x_err[NBINS_1] ={0.5,0.5,0.5,1,1,1.25,1.25,2.5,5};
  for(int i = 0; i<pp_crosssection_stat->GetN(); i++) {
    pp_crosssection_stat->SetPointEXhigh(i, x_err[i]);
    pp_crosssection_stat->SetPointEXlow(i, x_err[i]);
  }

  TGraphErrors *pp_crosssection_sys = (TGraphErrors*)f_pp_crosssection->Get("pp_crosssection_sys");
  TH1F *h_CR2 = (TH1F*)f_CR2->Get("h_CR2_ratio")->Clone("h_CR2");
  TGraph *ralf_pp_high = (TGraph*) f_ralf_pp_high->Get("MyGraph");
  TGraph *ralf_pp_low = (TGraph*) f_ralf_pp_low->Get("MyGraph");
  TGraph *ralf_pp = new TGraph(120);

  TGraph *greco_pp_high = (TGraph*) f_greco_high->Get("MyGraph");
  TGraph *greco_pp_low = (TGraph*) f_greco_low->Get("MyGraph");
  TGraph *greco_pp =  new TGraph(198);

  //try to get a band from Ralf pp predictions
  double x[60]={0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1,10.3,10.5,10.7,10.9,11.1,11.3,11.5,11.7,11.9};
  for (int j=0; j<60; j++)
  {
    ralf_pp->SetPoint(j,x[j],ralf_pp_high->Eval(x[j]));
    ralf_pp->SetPoint(60+j,x[60-j-1],ralf_pp_low->Eval(x[60-j-1]));
  }//for

  //try to get a band from Plumari predictions
  double x_greco[99]={0.1,0.20505,0.31010,0.41515,0.52020,0.62525,0.73030,0.83535,0.94040,1.04545,1.15051,1.25556,1.36061,1.46566,1.57071,1.67576,1.78081,1.88586,1.99091,2.09596,2.20101,2.30606,2.41111,2.51616,2.62121,2.72626,2.83131,2.93636,3.04141,3.14646,3.25152,3.35657,3.46162,3.56667,3.67172,3.77677,3.88182,3.98687,4.09192,4.19697,4.30202,4.40707,4.51212,4.61717,4.72222,4.82727,4.93232,5.03737,5.14242,5.24747,5.35253,5.45758,5.56263,5.66768,5.77273,5.87778,5.98283,6.08788,6.19293,6.29798,6.40303,6.50808,6.61313,6.71818,6.82323,6.92828,7.03333,7.13838,7.24343,7.34848,7.45354,7.55859,7.66364,7.76869,7.87374,7.97879,8.08384,8.18889,8.29394,8.39899,8.50404,8.60909,8.71414,8.81919,8.92424,9.02929,9.13434,9.23939,9.34444,9.44949,9.55455,9.6596,9.76465,9.8697,9.97475,10.0798,10.1848,10.2899,10.3949};
  for (int j=0;j<99;j++)
  {
    greco_pp->SetPoint(j,x_greco[j],greco_pp_high->Eval(x_greco[j]));
    greco_pp->SetPoint(99+j,x_greco[99-j-1],greco_pp_low->Eval(x_greco[99-j-1]));
  }

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  c1->SetLeftMargin(0.155);
  c1->SetRightMargin(0.03);
  c1->SetTopMargin(0.06);
  c1->SetBottomMargin(0.12);

  const Int_t NBINS = 9+2; // +2 leave space for plot edge
  Double_t edges[NBINS + 1]={2, 3,4,5,6,8,10,12.5,15,20,30,30.5};
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
  h_empty->SetAxisRange(0,1.0,"Y");
  h_empty->Draw("E");
  h_empty->SetMarkerColor(1);
  h_empty->SetMarkerStyle(28);
  h_empty->SetMarkerSize(2);
  h_empty->SetLineColor(1);

  h_CR2->SetLineColor(kMagenta-4);
  h_CR2->SetMarkerStyle(24);
  h_CR2->SetMarkerColor(kMagenta-4);
  h_CR2->Draw("Esame");

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

  ralf_pp->SetFillColor(kOrange-2);
  ralf_pp->SetFillColorAlpha(kOrange,0.4);
  ralf_pp->SetLineColor(kOrange-2);
  ralf_pp->SetFillStyle(1001);
  ralf_pp->Draw("Fsame");
  greco_pp->SetFillColor(kBlue-7);
  greco_pp->SetFillColorAlpha(kBlue-7,0.7);
  greco_pp->SetLineColor(kBlue-7);
  greco_pp->SetFillStyle(1001);
  greco_pp->Draw("Fsame");

  TLatex Tl;
  Tl.SetNDC();
  Tl.SetTextAlign(12);
  Tl.SetTextSize(0.05);
  Tl.SetTextFont(42);
  Tl.DrawLatex(0.18,0.91, "#font[61]{CMS}");
  Tl.DrawLatex(0.62,0.97, "#scale[0.8]{pp 252 nb^{-1} (5.02 TeV)}");//pp

  TLatex Tl2;
  Tl.SetNDC();
  Tl.SetTextAlign(12);
  Tl.SetTextSize(0.05*0.75);
  Tl.SetTextFont(42);
  Tl.DrawLatex(0.29,0.91, "#font[52]{Preliminary}");

  auto leg = new TLegend(0.49,0.70,0.92,0.90);
  leg->AddEntry(pp_crosssection_sys,"Data","fpe");
  leg->AddEntry(h_CR2,"CR2 prediction","p");
  leg->AddEntry(greco_pp,"PLB821 (2021) 136622","f");
  leg->AddEntry(ralf_pp,"PLB795 (2019) 117","f");
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetHeader("pp");
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->Draw();
  TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header->SetTextSize(0.045);

  TLatex* tex4;
  tex4 = new TLatex(0.5, 0.65, "Global uncertainty: 6.8%");
  tex4->SetNDC();
  tex4->SetTextFont(42);
  tex4->SetTextSize(0.04);
  tex4->SetLineWidth(2);
  tex4->Draw();

  TLatex* tex_6;
  tex_6 = new TLatex(0.18,0.85,"|y| < 1");
  tex_6->SetNDC();
  tex_6->SetTextFont(42);
  tex_6->SetTextSize(0.04);
  tex_6->SetLineWidth(2);
  tex_6->Draw();

  c1->SaveAs("plots/LcD0_pp.pdf");
}
