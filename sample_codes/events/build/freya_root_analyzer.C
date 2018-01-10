#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TROOT.h>
#include <iostream>
#include <cstdio>
#include "TFile.h"
#include "TChain.h"
#include <algorithm>

TFile *vetsex = new TFile("Cf252_sf_10k.dat.root", "READ");
TTree *mytree = (TTree *) gROOT->FindObject("FreyaTree");

//
Double_t pi = 3.14159265359;


//Photon emission multiplicities
TH1D *hframe_p_mult_0;
TH1D *hframe_p_mult_1;
TH1D *hframe_p_mult_2;
TH2D *hframe_p_multi;
TH3F *hframe_p_multi3D;
TH1D *h_p_mult_total;
int nbins_h_p_mult_total  = 20;

TH1D *hframe_ph_E_0;
TH1D *hframe_ph_E_1;
TH1D *hframe_ph_E_2;
TH1D *h_ph_E_total;
// TH1D *h_n_E_Boltzmann;
int nbins_h_ph_E_total  = 500;
int max_h_ph_E  = 7; // in MeV

void create_frames();

/*
List of meanings:

iAf1: fission fragment 1, mass number
iAf2: fission fragment 2, mass number
iAp1: product fission fragment 1, mass number (after initial fragment has decayed ?)
iAp2: product fission fragment 2, mass number (after initial fragment has decayed ?)
Ekin1: kinetic energy of 1 fragment
Ekin2: kinetic energy of second fragment
n0: neutrons from compound nuclei
n1: neutrons from first fission fragment
n2: neutrons from second fission fragment
m0: # pre-fission photons
m1: # photons from fragment 1
m2: # photons from fragment 2
P0: energy of neutrons emmited pre-fission (?)
P1: energy of neutrons emmited by F1
P2: energy of neutrons emitted by F2
Q0: energy of gammas emitted pre-fission
Q1: energy of gammas emitted by F1
Q2: energy of gamma emitted by F2
*/


void freya_root_analyzer() {

create_frames();

Double_t mean;
Double_t norm = 2;
Double_t nu;
Double_t value;

//////////////////////////////////////////////////////
//Photon multiplicities
//////////////////////////////////////////////////////

//Total photon multiplicity
TCanvas *c4 = new TCanvas("c4", "Photon multiplicities",150,10,990,660);

mytree->Draw("m0>>hframe_p_mult_0");
mean = hframe_p_mult_0->GetMean();
cout << "\n" << endl;
cout << "mean_number_of_photons_FF0: " << mean << endl;

mytree->Draw("m1>>hframe_p_mult_1");
mean = hframe_p_mult_1->GetMean();
cout << "mean_number_of_photons_FF1: " << mean << endl;

mytree->Draw("m2>>hframe_p_mult_2");
mean = hframe_p_mult_2->GetMean();
cout << "mean_number_of_photons_FF1: " << mean << endl;

mytree->Draw("m1:m2>>hframe_p_multi","","col");
mytree->Draw("m0:m2:m1>>hframe_p_multi3D","","lego");


Int_t lastbin1 = h_p_mult_total->GetNbinsX();
Double_t n_m0;
Double_t n_m1;
Double_t n_m2;
Double_t n_mtot;
Double_t weight_m;
for (int k=0;k<lastbin1;k++){
n_m0 = k-1;
hframe_p_multi3D->GetZaxis()->SetRange(k,k);
TH2D* pxy_m = (TH2D*) hframe_p_multi3D->Project3D("xy");
pxy_m->Draw("col");
  for (int i=0;i<lastbin1;i++){
    n_m1 = i-1; // bin 1 contains evetns of multiplicity 0
    auto px = pxy_m->ProjectionX("px",i,i);
      for (int j=0;j<lastbin1;j++){
      weight_m = px->GetBinContent(j);
      n_m2 = j-1; //bin 1 contains evetns of multiplicity 0
      if (n_m0<0){n_m0=0;} // bin 0 events (underflow) have to come in underflow
      if (n_m1<0){n_m1=0;} // bin 0 events (underflow) have to come in underflow
      if (n_m2<0){n_m2=0;} // bin 0 events (underflow) have to come in underflow
      n_mtot = n_m0+n_m1+n_m2;
      h_p_mult_total->Fill(n_mtot,weight_m);
      }
  }
}

// normalize the multiplicites
norm = h_p_mult_total->Integral();
h_p_mult_total->Scale(1/norm);

Double_t moments[5];

// factorial moments
std::fill_n(moments, 5, 0.);

// bin zero is underflow bin!
for(int i=1;i<nbins_h_p_mult_total+1;i++){
nu = i-1; // bin 1 is 0 neutrons ...
value = h_p_mult_total->GetBinContent(i);

//moments[0] = 0; // dummy
moments[1] += nu * value;
moments[2] += nu * (nu-1) * value;
moments[3] += nu * (nu-1) * (nu-2) * value;
moments[4] += nu * (nu-1) * (nu-2) * (nu-3) * value;
}

Double_t p_multiplicity = moments[1];

cout << "\n Average photon multiplicity: " << p_multiplicity << endl;

/*
cout << "\n Photons Moments" << endl;
// write to terminal
for(int i=1;i<5;i++){
cout << "Moment" << i << " " << moments[i] << endl;
}
*/

// x axis is such that [0,1] is multiplicity 0, ...
h_p_mult_total->GetXaxis()->SetTitle("Photon multiplicity PM");
h_p_mult_total->GetYaxis()->SetTitle("P(PM)");
h_p_mult_total->SetTitle("Photon multiplicity, total");
// h_n_mult_total->GetXaxis()->SetRangeUser(0.,nbins_h_n_mult_total);
h_p_mult_total->Draw();
c4->Update();


//Photon multiplicity, for different fragments
TCanvas *pMults = new TCanvas("pMults","Photon Multiplicities(2)",150,10,990,660);

norm = hframe_p_mult_2->Integral();
hframe_p_mult_2->SetLineStyle(15);
hframe_p_mult_2->SetLineColor(6);
hframe_p_mult_2->SetTitle("");
hframe_p_mult_2->Scale(1/norm);
hframe_p_mult_2->GetXaxis()->SetTitle("Photon multiplicity PM");
hframe_p_mult_2->GetYaxis()->SetTitle("P(PM)");
hframe_p_mult_2->SetTitle("Photonmultiplicity, per fragment");
hframe_p_mult_2->GetXaxis()->SetRangeUser(0,10);
hframe_p_mult_2->Draw();

norm = hframe_p_mult_1->Integral();
hframe_p_mult_1->Scale(1/norm);
hframe_p_mult_1->SetLineColor(2);
hframe_p_mult_1->SetLineStyle(10);
hframe_p_mult_1->Draw("same");

norm = hframe_p_mult_0->Integral();
hframe_p_mult_0->Scale(1/norm);
hframe_p_mult_0->SetLineColor(1);
hframe_p_mult_0->SetLineStyle(10);
hframe_p_mult_0->Draw("same");

//Plot total n multiplicity in same plot as above
h_p_mult_total->Draw("same");

auto legend = new TLegend(0.7,0.75,0.9,0.9);
legend->SetTextSize(0.03);
// legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
legend->AddEntry(hframe_p_mult_0,     "Prefission","l");
legend->AddEntry(hframe_p_mult_1,     "Fragment 1","l");
legend->AddEntry(hframe_p_mult_2,     "Fragment 2","l");
legend->AddEntry(h_p_mult_total,"Total","l");
legend->Draw();


/////////////////////////////////////////////////////////////////
// Photon energies
/////////////////////////////////////////////////////////////////

//Plot of photon energy spectrum
TCanvas *c6 = new TCanvas("c6","Photon Energy Spectrum ",150,10,990,660);
mytree->Draw("Q0>>hframe_ph_E_0");
mytree->Draw("Q1>>hframe_ph_E_1");
mytree->Draw("Q2>>hframe_ph_E_2");
h_ph_E_total->Add(hframe_ph_E_0,1.0);
h_ph_E_total->Add(hframe_ph_E_1,1.0);
h_ph_E_total->Add(hframe_ph_E_2,1.0);

//Average gamma ray energies
cout << "\n" << endl;
mean = hframe_ph_E_0->GetMean();
cout << "Mean photon energy F0: " << mean << " MeV" << endl;

mean = hframe_ph_E_1->GetMean();
cout << "Mean photon energy F1: " << mean << " MeV" << endl;

mean = hframe_ph_E_2->GetMean();
cout << "Mean photon energy F2: " << mean << " MeV" << endl;

cout << "\n" << endl;


// at the moment if no photon is emitted, this is tabulated in bin 0
cout << "Workaround as long as underflow bin is /WRONG/" << endl;
h_ph_E_total->SetBinContent(1,0);

c6->SetLogy();
h_ph_E_total->GetXaxis()->SetTitle("Photon Energy En [MeV]");
h_ph_E_total->GetYaxis()->SetTitle("Number of Photons");
h_ph_E_total->SetTitle("Photon spectrum");
h_ph_E_total->Draw("E");


//Total gamma ray energy


}


void create_frames() {
gStyle->SetOptStat(0);
int nbins;
int maxbin;

nbins = nbins_h_p_mult_total;
maxbin = nbins;
hframe_p_mult_0 = new TH1D("hframe_p_mult_0","",nbins,-0.5,maxbin-.5);
hframe_p_mult_1 = new TH1D("hframe_p_mult_1","",nbins,-0.5,maxbin-.5);
hframe_p_mult_2 = new TH1D("hframe_p_mult_2","",nbins,-0.5,maxbin-0.5);
hframe_p_multi  = new TH2D("hframe_p_multi","",nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5);
hframe_p_multi3D  = new TH3F("hframe_p_multi3D","",nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5);
h_p_mult_total = new TH1D("h_p_mult_total","",nbins,-0.5,maxbin-.5);

nbins = nbins_h_ph_E_total;
maxbin = max_h_ph_E;
hframe_ph_E_0 = new TH1D("hframe_ph_E_0","",nbins,0,maxbin);
hframe_ph_E_1 = new TH1D("hframe_ph_E_1","",nbins,0,maxbin);
hframe_ph_E_2 = new TH1D("hframe_ph_E_2","",nbins,0,maxbin);
h_ph_E_total = new TH1D("h_ph_E_total","",nbins,0,maxbin);

}

Double_t n_Boltzman_fct(Double_t *x, Double_t *par) {

    Double_t energy =x[0];
    Double_t T = par[0];
    Double_t binWidth = par[1];
    Double_t result;
    // result = T;
    result = binWidth * 2./sqrt(pi) * pow(T,-3./2.) * sqrt(energy) * exp(- energy / T );
    return result; 
    }

Double_t ph_Boltzman_fct(Double_t *x, Double_t *par) {

    Double_t energy =x[0];
    Double_t T = par[0];
    Double_t binWidth = par[1];
    Double_t result;
    // result = T;
    result = binWidth * pow(T,-2.) * energy * exp(- energy / T );
    return result; 
    }