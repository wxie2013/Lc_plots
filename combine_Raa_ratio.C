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

void combine_Raa_ratio(){
TFile *f_Raa_data = new TFile("root_file/Raa_090_010_withsys.root");

TGraphAsymmErrors *stat_0_90 = (TGraphAsymmErrors*) f_Raa_data->Get("stat_0_90");
TGraphAsymmErrors *stat_0_10 = (TGraphAsymmErrors*) f_Raa_data->Get("stat_0_10");
TGraphAsymmErrors *stat_10_30 = (TGraphAsymmErrors*) f_Raa_data->Get("stat_10_30");
TGraphAsymmErrors *stat_30_50 = (TGraphAsymmErrors*) f_Raa_data->Get("stat_30_50");
TGraphAsymmErrors *stat_50_90 = (TGraphAsymmErrors*) f_Raa_data->Get("stat_50_90");

TGraphErrors *sys_0_90 = (TGraphErrors*)f_Raa_data->Get("sys_0_90");
TGraphErrors *sys_0_10 = (TGraphErrors*)f_Raa_data->Get("sys_0_10");
TGraphErrors *sys_10_30 = (TGraphErrors*)f_Raa_data->Get("sys_10_30");
TGraphErrors *sys_30_50 = (TGraphErrors*)f_Raa_data->Get("sys_30_50");
TGraphErrors *sys_50_90 = (TGraphErrors*)f_Raa_data->Get("sys_50_90");

TCanvas *c1 = new TCanvas("c1","c1",600,600);
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);
c1->SetLeftMargin(0.155);
c1->SetRightMargin(0.03);
c1->SetTopMargin(0.06);
c1->SetBottomMargin(0.12);

const Int_t NBINS = 6+2; //+2 leave space for plot edge
Double_t edges[NBINS + 1]={5, 6,8,10,12.5,15,20,30, 30.5};
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
h_empty->GetYaxis()->SetTitle("R_{AA}");
h_empty->SetAxisRange(0,2.3,"Y");
h_empty->Draw("E");
h_empty->SetMarkerColor(1);
h_empty->SetMarkerStyle(28);
h_empty->SetMarkerSize(2);
h_empty->SetLineColor(1);

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
stat_0_90->Draw("psame");
stat_0_90->SetMarkerSize(1.8);

sys_0_10->SetMarkerColor(kAzure+1);
sys_0_10->SetMarkerStyle(21);
sys_0_10->SetMarkerSize(1.6);
sys_0_10->SetLineWidth(0);
sys_0_10->SetFillStyle(1001);
sys_0_10->SetLineColor(kAzure+1);
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

sys_10_30->SetMarkerColor(kOrange+0);
sys_10_30->SetMarkerStyle(22);
sys_10_30->SetMarkerSize(1.6);
sys_10_30->SetLineWidth(0);
sys_10_30->SetFillStyle(1001);
sys_10_30->SetLineColor(kOrange+0);
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
sys_10_30->Draw("2same");
stat_10_30->SetMarkerStyle(22);
stat_10_30->SetMarkerColor(kOrange+0);
stat_10_30->SetLineColor(kOrange+0);
stat_10_30->SetMarkerSize(1.6);
stat_10_30->Draw("psame");


sys_30_50->SetMarkerColor(kGreen+1);
sys_30_50->SetMarkerStyle(23);
sys_30_50->SetMarkerSize(1.6);
sys_30_50->SetLineWidth(0);
sys_30_50->SetFillStyle(1001);
sys_30_50->SetLineColor(kGreen+1);
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
sys_30_50->Draw("2same");
stat_30_50->SetMarkerStyle(23);
stat_30_50->SetMarkerColor(kGreen+1);
stat_30_50->SetLineColor(kGreen+1);
stat_30_50->SetMarkerSize(1.6);
stat_30_50->Draw("psame");

sys_50_90->SetMarkerColor(kPink+9);
sys_50_90->SetMarkerStyle(33);
sys_50_90->SetMarkerSize(1.6);
sys_50_90->SetLineWidth(0);
sys_50_90->SetFillStyle(1001);
sys_50_90->SetLineColor(kPink+9);
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
sys_50_90->Draw("2same");
stat_50_90->SetMarkerStyle(33);
stat_50_90->SetMarkerColor(kPink+9);
stat_50_90->SetLineColor(kPink+9);
stat_50_90->SetMarkerSize(1.6);
stat_50_90->Draw("psame");

//stat_0_10->Draw("psame");
//stat_0_90->Draw("psame");

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
Tl.DrawLatex(0.18,0.91, "#font[61]{CMS}");
Tl.DrawLatex(0.37,0.97, "#scale[0.8]{PbPb 0.58 nb^{-1}, pp 252 nb^{-1} (5.02 TeV)}");//pp

TLatex Tl2;
Tl.SetNDC();
Tl.SetTextAlign(12);
Tl.SetTextSize(0.05*0.75);
Tl.SetTextFont(42);
Tl.DrawLatex(0.29,0.91, "#font[52]{Preliminary}");

TLatex* tex_6;
tex_6 = new TLatex(0.18,0.85,"|y| < 1");
tex_6->SetNDC();
tex_6->SetTextFont(42);
tex_6->SetTextSize(0.04);
tex_6->SetLineWidth(2);
tex_6->Draw();



auto leg1 = new TLegend(0.18,0.56,0.4,0.8);
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

TLine *l = new TLine(5,1,30.5,1);
l->SetLineStyle(7);
l->Draw("same");

double x_global[1]={29.9};
double y_global[1]={1};
double zero_global[1]={0.6};
double sys_global[1]={0.192};
TGraphErrors *global_uncer = new TGraphErrors(1,x_global,y_global,zero_global,sys_global);
global_uncer->SetFillStyle(1001);
global_uncer->SetLineColor(16);
global_uncer->SetFillColorAlpha(kGray+2,0.4);
global_uncer->SetName("global_uncer");
global_uncer->Draw("5same");

TLatex* tex_global;
tex_global = new TLatex(0.66,0.5,"Global uncertainty");
tex_global->SetNDC();
tex_global->SetTextFont(42);
tex_global->SetTextSize(0.038);
tex_global->SetTextColor(kGray+2);
tex_global->SetLineWidth(2);
tex_global->Draw();

c1->SaveAs("plots/Raa_0_90_010_combination_withpT.pdf");

}
