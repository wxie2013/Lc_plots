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

TGraphAsymmErrors get_graph_from_histo(TH1F*);

void pp_pt_spectra_data_model(){

  TFile *f_pp_crosssection = new TFile("root_file/pp_crosssection_withsys.root");
  TFile *f_PbPb_yield = new TFile("root_file/differential_yield_090_010_withsys.root");
  TFile *f_pythia8 = new TFile("root_file/cross_CR2_soft_nonDiffractive_ON_1B_EvtGen.root");
  TFile *f_pp_GMVFNS_1y1 = new TFile("root_file/GMVFNS_pp_model.root");
  TFile *f_pp_GMVFNS_299 = new TFile("root_file/GMVFNS_pp_model_229.root");

  //
  TH1F *h_pythia8_prompt = (TH1F*)f_pythia8->Get("hLc"); // inc Lc
  TH1F *h_pythia8_nonprompt = (TH1F*)f_pythia8->Get("hB2Lc");
  h_pythia8_prompt->Add(h_pythia8_nonprompt, -1); // prompt Lc
  h_pythia8_prompt->Scale(1e-6);  // from pb to ub
  TGraphAsymmErrors gra_tmp = get_graph_from_histo(h_pythia8_prompt);
  TGraphAsymmErrors *pythia8 = &gra_tmp;
  //
  TGraphAsymmErrors *pp_crosssection_stat = (TGraphAsymmErrors*) f_pp_crosssection->Get("pp_crosssection_stat");
  TGraphAsymmErrors *pp_GMVFNS_1y1 = (TGraphAsymmErrors*)f_pp_GMVFNS_1y1->Get("GM_VFNS_model_pp_Michael_1y1");
  TGraphAsymmErrors *pp_GMVFNS_299 = (TGraphAsymmErrors*)f_pp_GMVFNS_299->Get("GM_VFNS_model_pp_Michael_229");
  TGraphErrors *pp_crosssection_sys = (TGraphErrors*)f_pp_crosssection->Get("pp_crosssection_sys");
  TGraphAsymmErrors *stat_0_90 = (TGraphAsymmErrors*) f_PbPb_yield->Get("stat_0_90");
  TGraphAsymmErrors *stat_0_10 = (TGraphAsymmErrors*) f_PbPb_yield->Get("stat_0_10");
  TGraphAsymmErrors *stat_10_30 = (TGraphAsymmErrors*) f_PbPb_yield->Get("stat_10_30");
  TGraphAsymmErrors *stat_30_50 = (TGraphAsymmErrors*) f_PbPb_yield->Get("stat_30_50");
  TGraphAsymmErrors *stat_50_90 = (TGraphAsymmErrors*) f_PbPb_yield->Get("stat_50_90");
  TGraphErrors *sys_0_90 = (TGraphErrors*)f_PbPb_yield->Get("sys_0_90");
  TGraphErrors *sys_0_10 = (TGraphErrors*)f_PbPb_yield->Get("sys_0_10");
  TGraphErrors *sys_10_30 = (TGraphErrors*)f_PbPb_yield->Get("sys_10_30");
  TGraphErrors *sys_30_50 = (TGraphErrors*)f_PbPb_yield->Get("sys_30_50");
  TGraphErrors *sys_50_90 = (TGraphErrors*)f_PbPb_yield->Get("sys_50_90");

  //TCanvas *c1 = new TCanvas("c1","c1",600,600);
  TCanvas *c1 = new TCanvas("c1","c1",800,900);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  TPad *pad1 = new TPad("pad1","top pad",0.0,0.3,1.0,1.0);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.05);
  pad1->SetLeftMargin(0.155);
  pad1->Draw();

  TPad *pad2 = new TPad("pad2","bottom pad",0.0,0.0,1.0,0.3);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.35);
  pad2->SetRightMargin(0.05);
  pad2->SetLeftMargin(0.155);
  pad2->Draw();
  /*
     c1->SetLeftMargin(0.155);
     c1->SetRightMargin(0.03);
     c1->SetTopMargin(0.06);
     c1->SetBottomMargin(0.12);
     */

  pad1->cd();
  const Int_t NBINS = 9+2; //+2 mean leave space for plot edges
  Double_t edges[NBINS + 1]={2.0, 3, 4,5,6,8,10,12.5,15,20,30, 30.5};
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
  h_empty->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}}(or #frac{1}{T_{AA}} #frac{dN}{dp_{T}}) (#mub GeV^{-1}c)");
  h_empty->SetAxisRange(4E-5,400,"Y");
  h_empty->SetMinimum(1.2e-3);
  h_empty->Draw("E");
  h_empty->SetMarkerColor(1);
  h_empty->SetMarkerStyle(28);
  h_empty->SetMarkerSize(2);
  h_empty->SetLineColor(1);
  gPad->SetLogy();

  pp_crosssection_sys->SetMarkerColor(1);
  pp_crosssection_sys->SetMarkerStyle(20);
  pp_crosssection_sys->SetMarkerSize(1.8);
  pp_crosssection_sys->SetLineWidth(0);
  pp_crosssection_sys->SetFillStyle(1001);
  pp_crosssection_sys->SetLineColor(1);
  pp_crosssection_sys->SetFillColorAlpha(16,0.4);
  pp_crosssection_sys->GetXaxis()->CenterTitle();
  pp_crosssection_sys->GetYaxis()->CenterTitle();
  pp_crosssection_sys->GetXaxis()->SetTitleOffset(1.0);
  pp_crosssection_sys->GetYaxis()->SetTitleOffset(1.0);
  pp_crosssection_sys->GetXaxis()->SetLabelOffset(0.007);
  pp_crosssection_sys->GetYaxis()->SetLabelOffset(0.007);
  pp_crosssection_sys->GetXaxis()->SetTitleSize(0.035);
  pp_crosssection_sys->GetYaxis()->SetTitleSize(0.035);
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
  //sys_0_90->Draw("2same");
  stat_0_90->SetMarkerStyle(29);
  stat_0_90->SetMarkerColor(2);
  stat_0_90->SetLineColor(2);
  //stat_0_90->Draw("psame");
  stat_0_90->SetMarkerSize(1.8);

  sys_0_10->SetMarkerColor(kAzure+1);
  sys_0_10->SetMarkerStyle(21);
  sys_0_10->SetMarkerSize(1.6);
  sys_0_10->SetLineWidth(0);
  sys_0_10->SetFillStyle(1001);
  sys_0_10->SetLineColor(8);
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
  //sys_0_10->Draw("2same");
  stat_0_10->SetMarkerStyle(21);
  stat_0_10->SetMarkerColor(kAzure+1);
  stat_0_10->SetLineColor(kAzure+1);
  stat_0_10->SetMarkerSize(1.6);
  //stat_0_10->Draw("psame");


  sys_10_30->SetMarkerColor(kOrange+0);
  sys_10_30->SetMarkerStyle(22);
  sys_10_30->SetMarkerSize(1.6);
  sys_10_30->SetLineWidth(0);
  sys_10_30->SetFillStyle(1001);
  sys_10_30->SetLineColor(8);
  sys_10_30->SetFillColorAlpha(kOrange+0,0.4);
  sys_10_30->GetXaxis()->CenterTitle();
  sys_10_30->GetYaxis()->CenterTitle();
  sys_10_30->GetXaxis()->SetTitleOffset(1.0);
  sys_10_30->GetYaxis()->SetTitleOffset(1.0);
  sys_10_30->GetXaxis()->SetLabelOffset(0.007);
  sys_10_30->GetYaxis()->SetLabelOffset(0.007);
  sys_10_30->GetXaxis()->SetTitleSize(0.045);
  sys_10_30->GetYaxis()->SetTitleSize(0.045);
  sys_10_30->GetXaxis()->SetTitleFont(42);
  sys_10_30->GetYaxis()->SetTitleFont(42);
  sys_10_30->GetXaxis()->SetLabelFont(42);
  sys_10_30->GetXaxis()->SetLabelFont(42);
  sys_10_30->GetXaxis()->SetLabelSize(0.06);
  sys_10_30->GetYaxis()->SetLabelSize(0.06);
  sys_10_30->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  sys_10_30->GetYaxis()->SetTitle("R_{AA}");
  //sys_10_30->Draw("2same");
  stat_10_30->SetMarkerStyle(22);
  stat_10_30->SetMarkerColor(kOrange+0);
  stat_10_30->SetLineColor(kOrange+0);
  stat_10_30->SetMarkerSize(1.6);
  //stat_10_30->Draw("psame");

  sys_30_50->SetMarkerColor(kGreen+1);
  sys_30_50->SetMarkerStyle(39);
  sys_30_50->SetMarkerSize(1.6);
  sys_30_50->SetLineWidth(0);
  sys_30_50->SetFillStyle(1001);
  sys_30_50->SetLineColor(8);
  sys_30_50->SetFillColorAlpha(kGreen+1,0.4);
  sys_30_50->GetXaxis()->CenterTitle();
  sys_30_50->GetYaxis()->CenterTitle();
  sys_30_50->GetXaxis()->SetTitleOffset(1.0);
  sys_30_50->GetYaxis()->SetTitleOffset(1.0);
  sys_30_50->GetXaxis()->SetLabelOffset(0.007);
  sys_30_50->GetYaxis()->SetLabelOffset(0.007);
  sys_30_50->GetXaxis()->SetTitleSize(0.045);
  sys_30_50->GetYaxis()->SetTitleSize(0.045);
  sys_30_50->GetXaxis()->SetTitleFont(42);
  sys_30_50->GetYaxis()->SetTitleFont(42);
  sys_30_50->GetXaxis()->SetLabelFont(42);
  sys_30_50->GetXaxis()->SetLabelFont(42);
  sys_30_50->GetXaxis()->SetLabelSize(0.06);
  sys_30_50->GetYaxis()->SetLabelSize(0.06);
  sys_30_50->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  sys_30_50->GetYaxis()->SetTitle("R_{AA}");
  //sys_30_50->Draw("2same");
  stat_30_50->SetMarkerStyle(39);
  stat_30_50->SetMarkerColor(kGreen+1);
  stat_30_50->SetLineColor(kGreen+1);
  stat_30_50->SetMarkerSize(1.6);
  //stat_30_50->Draw("psame");


  sys_50_90->SetMarkerColor(kPink+9);
  sys_50_90->SetMarkerStyle(33);
  sys_50_90->SetMarkerSize(1.6);
  sys_50_90->SetLineWidth(0);
  sys_50_90->SetFillStyle(1001);
  sys_50_90->SetLineColor(8);
  sys_50_90->SetFillColorAlpha(kPink+9,0.4);
  sys_50_90->GetXaxis()->CenterTitle();
  sys_50_90->GetYaxis()->CenterTitle();
  sys_50_90->GetXaxis()->SetTitleOffset(1.0);
  sys_50_90->GetYaxis()->SetTitleOffset(1.0);
  sys_50_90->GetXaxis()->SetLabelOffset(0.007);
  sys_50_90->GetYaxis()->SetLabelOffset(0.007);
  sys_50_90->GetXaxis()->SetTitleSize(0.045);
  sys_50_90->GetYaxis()->SetTitleSize(0.045);
  sys_50_90->GetXaxis()->SetTitleFont(42);
  sys_50_90->GetYaxis()->SetTitleFont(42);
  sys_50_90->GetXaxis()->SetLabelFont(42);
  sys_50_90->GetXaxis()->SetLabelFont(42);
  sys_50_90->GetXaxis()->SetLabelSize(0.06);
  sys_50_90->GetYaxis()->SetLabelSize(0.06);
  sys_50_90->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  sys_50_90->GetYaxis()->SetTitle("R_{AA}");
  //sys_50_90->Draw("2same");
  stat_50_90->SetMarkerStyle(33);
  stat_50_90->SetMarkerColor(kPink+9);
  stat_50_90->SetLineColor(kPink+9);
  stat_50_90->SetMarkerSize(1.6);
  //stat_50_90->Draw("psame");

  pythia8->SetMarkerColor(kOrange+5);
  pythia8->SetLineColor(kOrange+5);
  pythia8->SetMarkerStyle(46);
  pythia8->SetMarkerSize(1.6);
  pythia8->Draw("psame");

  pp_GMVFNS_1y1->SetMarkerColor(kCyan-7);
  pp_GMVFNS_1y1->SetMarkerStyle(24);
  pp_GMVFNS_1y1->SetMarkerSize(1.6);
  pp_GMVFNS_1y1->SetLineColor(kCyan-7);
  pp_GMVFNS_1y1->Draw("psame");

  pp_GMVFNS_299->SetMarkerColor(kViolet+1);
  pp_GMVFNS_299->SetMarkerStyle(26);
  pp_GMVFNS_299->SetMarkerSize(1.6);
  pp_GMVFNS_299->SetLineColor(kViolet+1);
  pp_GMVFNS_299->Draw("psame");

  TLatex* tex;
  tex =  new TLatex(0.79,0.84,"#font[61]{#frac{#Lambda_{c}^{+} + #Lambda_{c}^{#font[122]{\55}}}{2}}");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.05);
  tex->SetLineWidth(4);
  tex->Draw();

  TLatex Tl;
  Tl.SetNDC();
  Tl.SetTextAlign(12);
  Tl.SetTextSize(0.05);
  Tl.SetTextFont(42);
  Tl.DrawLatex(0.18,0.92, "#font[61]{CMS}");
  Tl.DrawLatex(0.67,0.97, "#scale[0.8]{pp 252 nb^{-1} (5.02 TeV)}");//pp

  TLatex Tl2;
  Tl.SetNDC();
  Tl.SetTextAlign(12);
  Tl.SetTextSize(0.05*0.75);
  Tl.SetTextFont(42);
  Tl.DrawLatex(0.29,0.92, "#font[52]{Preliminary}");

  auto leg = new TLegend(0.59,0.62,0.81,0.83);
  leg->AddEntry(pp_crosssection_sys,"Data","fpe");
  leg->AddEntry(pp_GMVFNS_1y1,"GMVFNS-1","p");
  leg->AddEntry(pp_GMVFNS_299,"GMVFNS-2","p");
  leg->AddEntry(pythia8,"PYTHIA 8 + CR","p");
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetHeader("pp");
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->Draw();
  TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header->SetTextSize(0.045);
  /*
     auto leg1 = new TLegend(0.59,0.34,0.81,0.6);
     leg1->AddEntry(sys_0_90,"0-90%","fp");
     leg1->AddEntry(sys_0_10,"0-10%","fp");
     leg1->AddEntry(sys_10_30,"10-30%","fp");
     leg1->AddEntry(sys_30_50,"30-50%","fp");
     leg1->AddEntry(sys_50_90,"50-90%","fp");
     leg1->SetTextSize(0.04);
     leg1->SetTextFont(42);
     leg1->SetHeader("PbPb Cent.");
     leg1->SetTextFont(42);
     leg1->SetBorderSize(0);
     leg1->Draw();
     TLegendEntry *header1 = (TLegendEntry*)leg1->GetListOfPrimitives()->First();
     header1->SetTextSize(0.045);
     */

  TLatex* tex4;
  tex4 = new TLatex(0.29,0.79,"Global uncertainty");
  tex4->SetNDC();
  tex4->SetTextFont(42);
  tex4->SetTextSize(0.04);
  tex4->SetLineWidth(2);
  tex4->Draw();

  TLatex* tex5;
  tex5 = new TLatex(0.29,0.75,"pp: 13.3%");
  tex5->SetNDC();
  tex5->SetTextFont(42);
  tex5->SetTextSize(0.04);
  tex5->SetLineWidth(2);
  tex5->Draw();
  /*
     TLatex* tex6;
     tex6 = new TLatex(0.29,0.7,"PbPb: 16.0%");
     tex6->SetNDC();
     tex6->SetTextFont(42);
     tex6->SetTextSize(0.04);
     tex6->SetLineWidth(2);
     tex6->Draw();
     */


  TLatex* tex_6;
  tex_6 = new TLatex(0.18,0.87,"|y| < 1");
  tex_6->SetNDC();
  tex_6->SetTextFont(42);
  tex_6->SetTextSize(0.04);
  tex_6->SetLineWidth(2);
  tex_6->Draw();

  pad2->cd();
  TH1 *h_empty2 = new TH1D("h_empty2","h_empty2",NBINS,edges);

  h_empty2->SetAxisRange(0,7,"Y");
  h_empty2->SetNdivisions(5,"Y");
  h_empty2->GetXaxis()->CenterTitle();
  h_empty2->GetYaxis()->CenterTitle();
  h_empty2->GetXaxis()->SetTitleOffset(0.96);
  h_empty2->GetYaxis()->SetTitleOffset(0.42);
  h_empty2->GetXaxis()->SetLabelOffset(0.004);
  h_empty2->GetYaxis()->SetLabelOffset(0.007);
  h_empty2->GetXaxis()->SetTitleSize(0.1);
  h_empty2->GetYaxis()->SetTitleSize(0.1);
  h_empty2->GetXaxis()->SetTitleFont(42);
  h_empty2->GetYaxis()->SetTitleFont(42);
  h_empty2->GetXaxis()->SetLabelFont(42);
  h_empty2->GetYaxis()->SetLabelFont(42);
  h_empty2->GetXaxis()->SetLabelSize(0.1);
  h_empty2->GetYaxis()->SetLabelSize(0.1);
  h_empty2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_empty2->GetYaxis()->SetTitle("#frac{pp Data}{GM-VFNS}");
  h_empty2->Draw("E");
  h_empty2->SetMarkerColor(2);
  h_empty2->SetMarkerStyle(20);
  h_empty2->SetMarkerSize(2.2);
  h_empty2->SetLineColor(2);
  h_empty2->SetMaximum(8.8);
  h_empty2->Draw();


  const Int_t NBINS_1 = 9; 
  Double_t edges_1[NBINS_1 + 1]={3,4,5,6,8,10,12.5,15,20,30};
  double zero[NBINS_1] ={0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double x_err[NBINS_1] ={0.5,0.5,0.5,1,1,1.25,1.25,2.5,5};
  double x[NBINS_1];
  double py[NBINS_1];
  double err_high[NBINS_1];
  double err_low[NBINS_1];
  double err_sys[NBINS_1];
  double y_GM1[NBINS_1];
  double err_high_GM1[NBINS_1];
  double err_low_GM1[NBINS_1];
  double y_GM299[NBINS_1];
  double err_high_GM299[NBINS_1];
  double err_low_GM299[NBINS_1];
  double y_pythia[NBINS_1];
  double err_high_pythia[NBINS_1];
  double err_low_pythia[NBINS_1];
  double ratio_GM1[NBINS_1];
  double ratio_GM299[NBINS_1];
  double ratio_pythia[NBINS_1];
  double ratio_pythia_high[NBINS_1];
  double ratio_pythia_low[NBINS_1];
  double ratio_data_GM1_high[NBINS_1];
  double ratio_data_GM1_low[NBINS_1];
  double ratio_data_GM299_high[NBINS_1];
  double ratio_data_GM299_low[NBINS_1];
  double pythia_sys_low[NBINS_1];
  double pythia_sys_high[NBINS_1];
  double GM1_sys_low[NBINS_1];
  double GM1_sys_high[NBINS_1];
  double GM299_sys_low[NBINS_1];
  double GM299_sys_high[NBINS_1];
  double pythia_model_low[NBINS_1];
  double pythia_model_high[NBINS_1];
  double GM1_model_low[NBINS_1];
  double GM1_model_high[NBINS_1];
  double GM299_model_low[NBINS_1];
  double GM299_model_high[NBINS_1];

  for (int i=0; i<NBINS_1; i++)
  {
    x[i]=(edges_1[i]+edges_1[i+1])/2;
    py[i]= pp_crosssection_stat->Eval(x[i]);
    err_high[i]=pp_crosssection_stat->GetErrorYhigh(i);
    err_low[i]=pp_crosssection_stat->GetErrorYlow(i);
    err_sys[i] = pp_crosssection_sys->GetErrorY(i);

    y_pythia[i]=pythia8->Eval(x[i]);
    err_high_pythia[i]=pythia8->GetErrorYhigh(i);
    err_low_pythia[i]=pythia8->GetErrorYlow(i);

    y_GM1[i]=pp_GMVFNS_1y1->Eval(x[i]);
    err_high_GM1[i] = pp_GMVFNS_1y1->GetErrorYhigh(i);
    err_low_GM1[i] = pp_GMVFNS_1y1->GetErrorYlow(i);

    y_GM299[i] = pp_GMVFNS_299->Eval(x[i]);
    err_high_GM299[i] = pp_GMVFNS_299->GetErrorYhigh(i);
    err_low_GM299[i] = pp_GMVFNS_299->GetErrorYlow(i);

    // (data +/ stat)/model
    ratio_GM1[i]=py[i]/y_GM1[i];
    ratio_GM299[i]=py[i]/y_GM299[i];
    ratio_pythia[i]=py[i]/y_pythia[i];	
    ratio_data_GM1_high[i] = err_high[i]/y_GM1[i];
    ratio_data_GM1_low[i] = err_low[i]/y_GM1[i];
    ratio_data_GM299_high[i] = err_high[i]/y_GM299[i];
    ratio_data_GM299_low[i] = err_low[i]/y_GM299[i];
    ratio_pythia_high[i]=err_high[i]/y_pythia[i];
    ratio_pythia_low[i]=err_low[i]/y_pythia[i];

    // (data +/ sys)/model
    GM1_sys_low[i] = err_sys[i]/y_GM1[i];
    GM1_sys_high[i] = err_sys[i]/y_GM1[i];
    GM299_sys_high[i] = err_sys[i]/y_GM299[i];
    GM299_sys_low[i] = err_sys[i]/y_GM299[i];
    pythia_sys_high[i]=err_sys[i]/y_pythia[i];
    pythia_sys_low[i]=err_sys[i]/y_pythia[i];

    // data/(model +/ theory error)
    GM1_model_low[i] = abs(py[i]/(y_GM1[i]+err_high_GM1[i]) - py[i]/y_GM1[i]);
    GM1_model_high[i] = abs(py[i]/(y_GM1[i]-err_low_GM1[i]) - py[i]/y_GM1[i]);
    GM299_model_high[i] = abs(py[i]/(y_GM299[i]+err_high_GM299[i]) - py[i]/y_GM299[i]);
    GM299_model_low[i] = abs(py[i]/(y_GM299[i]-err_low_GM299[i]) - py[i]/y_GM299[i]);
    pythia_model_high[i]= abs(py[i]/(y_pythia[i] + err_high_pythia[i]) - py[i]/y_pythia[i]);
    pythia_model_low[i]= abs(py[i]/(y_pythia[i] - err_low_pythia[i]) - py[i]/y_pythia[i]);
  }

  TGraphAsymmErrors *g_cross_stat_GM1 = new TGraphAsymmErrors(9,x,ratio_GM1,x_err,x_err,ratio_data_GM1_low,ratio_data_GM1_high);
  TGraphAsymmErrors *g_cross_stat_GM299 = new TGraphAsymmErrors(9,x,ratio_GM299,x_err,x_err,ratio_data_GM299_low,ratio_data_GM299_high);
  TGraphAsymmErrors *g_cross_stat_pythia = new TGraphAsymmErrors(9,x,ratio_pythia,x_err,x_err,ratio_pythia_low,ratio_pythia_high);

  TGraphAsymmErrors *g_cross_GM1 = new TGraphAsymmErrors(9,x,ratio_GM1,zero,zero,GM1_sys_low,GM1_sys_high);
  TGraphAsymmErrors *g_cross_GM299 = new TGraphAsymmErrors(9,x,ratio_GM299,zero,zero,GM299_sys_low,GM299_sys_high);
  TGraphAsymmErrors *g_cross_pythia = new TGraphAsymmErrors(9,x,ratio_pythia,zero,zero,pythia_sys_low,pythia_sys_high);

  TGraphAsymmErrors *g_cross_model_GM1 = new TGraphAsymmErrors(9,x,ratio_GM1,x_err,x_err,GM1_model_low,GM1_model_high);
  TGraphAsymmErrors *g_cross_model_GM299 = new TGraphAsymmErrors(9,x,ratio_GM299,x_err,x_err,GM299_model_low,GM299_model_high);
  TGraphAsymmErrors *g_cross_model_pythia = new TGraphAsymmErrors(9,x,ratio_pythia,x_err,x_err,pythia_model_low,pythia_model_high);

  g_cross_model_pythia->SetMarkerStyle(46);
  g_cross_model_pythia->SetMarkerSize(1);
  g_cross_model_pythia->SetLineWidth(0);
  g_cross_model_pythia->SetFillStyle(1001);
  g_cross_model_pythia->SetLineColor(kOrange+5);
  g_cross_model_pythia->SetFillColorAlpha(kOrange+5,0.5);
  g_cross_model_pythia->Draw("2same");

  g_cross_model_GM1->SetMarkerStyle(24);
  g_cross_model_GM1->SetMarkerSize(1);
  g_cross_model_GM1->SetLineWidth(0);
  g_cross_model_GM1->SetFillStyle(1001);
  g_cross_model_GM1->SetLineColor(kCyan-7);
  g_cross_model_GM1->SetFillColorAlpha(kCyan-7,0.4);
  g_cross_model_GM1->Draw("2same");


  g_cross_model_GM299->SetMarkerStyle(26);
  g_cross_model_GM299->SetMarkerSize(1);
  g_cross_model_GM299->SetLineWidth(0);
  g_cross_model_GM299->SetFillStyle(1001);
  g_cross_model_GM299->SetLineColor(kViolet+1);
  g_cross_model_GM299->SetFillColorAlpha(kViolet+1,0.3);
  g_cross_model_GM299->Draw("2same");

  g_cross_pythia->SetMarkerStyle(46);
  g_cross_pythia->SetMarkerSize(1);
  g_cross_pythia->SetFillStyle(1001);
  g_cross_pythia->SetLineColor(kOrange+5);
  g_cross_pythia->SetMarkerColor(kOrange+5);
  g_cross_pythia->SetFillColorAlpha(kOrange+5,0.5);
  g_cross_pythia->Draw("[]same");

  g_cross_GM1->SetMarkerStyle(24);
  g_cross_GM1->SetMarkerSize(1);
  g_cross_GM1->SetFillStyle(1001);
  g_cross_GM1->SetLineColor(kCyan-7);
  g_cross_GM1->SetMarkerColor(kCyan-7);
  g_cross_GM1->SetFillColorAlpha(kCyan-7,0.4);
  g_cross_GM1->Draw("[]same");

  g_cross_GM299->SetMarkerStyle(26);
  g_cross_GM299->SetMarkerSize(1);
  g_cross_GM299->SetFillStyle(1001);
  g_cross_GM299->SetLineColor(kViolet+1);
  g_cross_GM299->SetMarkerColor(kViolet+1);
  g_cross_GM299->SetFillColorAlpha(kViolet+1,0.3);
  g_cross_GM299->Draw("[]same");

  g_cross_stat_GM1->SetMarkerStyle(24);
  g_cross_stat_GM1->SetMarkerColor(kCyan-7);
  g_cross_stat_GM1->SetMarkerSize(1);
  g_cross_stat_GM1->SetLineColor(kCyan-7);
  g_cross_stat_GM1->Draw("psame");

  g_cross_stat_GM299->SetMarkerStyle(26);
  g_cross_stat_GM299->SetMarkerColor(kViolet+1);
  g_cross_stat_GM299->SetMarkerSize(1);
  g_cross_stat_GM299->SetLineColor(kViolet+1);
  g_cross_stat_GM299->Draw("psame");

  g_cross_stat_pythia->SetMarkerStyle(46);
  g_cross_stat_pythia->SetMarkerColor(kOrange+5);
  g_cross_stat_pythia->SetMarkerSize(1);
  g_cross_stat_pythia->SetLineColor(kOrange+5);
  g_cross_stat_pythia->Draw("psame");


  //

  TLine *l = new TLine (2,1,30,1);
  l->SetLineStyle(2);
  l->Draw("same");

  c1->SaveAs("plots/crosssection_pp_PbPb_onlypp.pdf");
}

// convert a histogram to a TGraphAsymmErrors
TGraphAsymmErrors get_graph_from_histo(TH1F* h)
{
  const int n = h->GetNbinsX();
  Double_t x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];
  for(int i = 0; i<n; i++) {
    x[i] = h->GetBinCenter(i+1);
    y[i] = h->GetBinContent(i+1);
    exl[i] = 0;
    exh[i] = 0;
    eyl[i] = h->GetBinError(i+1);
    eyh[i] = h->GetBinError(i+1);
  }

  return TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
}
