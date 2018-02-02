#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TROOT.h>
#include <iostream>
#include <cstdio>
#include "TFile.h"
#include "TChain.h"
#include <algorithm>

//CHANGED VERSION, THAT I AM GOING TO WORK OUT


TFile *vetsex = new TFile("Cf252.dat.root", "READ");
TTree *mytree = (TTree *) gROOT->FindObject("FreyaTree");

//
Double_t pi = 3.14159265359;

//Run in root with .x fission_script_root6.C

//initialise histograms

//Fragment yield histograms
TH1D *hframe_fragyield;
TH1D *hframe_fragyield2;
TH1D *htotal_fragyield;

//
TH1D *hframe_prodyield;
TH1D *hframe_prodyield2;
TH1D *htotal_prodyield;

//Neutron yields
TH1D *hframe_n_yield;
TH1D *hframe_n_yield2;
TH1D *htotal_n_yield;


//Kinetic energy distribution histograms
TH2D *hframe_Ekin1;
TH2D *hframe_Ekin2;
TH2D *h_Ekin;

TH1D *hframe_kin;
TH1D *hframe_kin2;
TH1D *htotal_kin; 

//Neutron emission multiplicities
TH1D *hframe_n0;
TH1D *hframe_n1;
TH1D *hframe_n2;
TH2D *hframe_n_multi;
TH3F *hframe_n_multi3D;
TH1D *h_n_mult_total;
int nbins_h_n_mult_total  = 12;

//Photon emission multiplicities
TH1D *hframe_p_mult_0;
TH1D *hframe_p_mult_1;
TH1D *hframe_p_mult_2;
TH2D *hframe_p_multi;
TH3F *hframe_p_multi3D;
TH1D *h_p_mult_total;
int nbins_h_p_mult_total  = 20;

//Energy of neutrons?
TH1D *hframe_n_E_0;
TH1D *hframe_n_E_1;
TH1D *hframe_n_E_2;
TH1D *h_n_E_total;
// TH1D *h_n_E_Boltzmann;
int nbins_h_n_E_total  = 500;
int max_h_n_E  = 20; // in MeV

TH1D *hframe_ph_E_0;
TH1D *hframe_ph_E_1;
TH1D *hframe_ph_E_2;
TH1D *h_ph_E_total;
// TH1D *h_n_E_Boltzmann;
int nbins_h_ph_E_total  = 500;
int max_h_ph_E  = 7; // in MeV

Double_t n_Boltzman_fct(Double_t *x, Double_t *par);
Double_t ph_Boltzman_fct(Double_t *x, Double_t *par);
// void n_Boltzman();

//initialise functions
void create_frames();

void fission_script_root6() {

int F = 100000;

create_frames();

/*
List of meanings:
iZ0: compound nucleus, charge
iZ1: fission fragment 1, charge
iZ2: fission fragment 2, charge
n0: neutrons from compound nuclei
n1: neutrons from first fission fragment
n2: neutrosn from second fission fragment
m0: # pre-fission photons
m1: # photons from fragment 1
m2: # photons from fragment 2
iA0: compound nucleus mass number
iAf1: fission fragment 1, mass number
iAf2: fission fragment 2, mass number
iAp1: product fission fragment 1, mass number (after initial fragment has decayed ?)
iAp2: product fission fragment 2, mass number (after initial fragment has decayed ?)
Ekin1: kinetic energy of 1 fragment
Ekin2: kinetic energy of second fragment
E0: exitation energy of compound nucleus
E1: excitation energy of F1
E2: excitation energy of F2 
P0: energy of neutrons emmited pre-fission
P1: energy of neutrons emmited by F1
P2: energy of neutrons emitted by F2
P0(x,y,z): direction of neutrons emmited by compound nucleus
P1(x,y,z): directions of neutrons emmited by F1
P2(x,y,z): directions of neutrons emmited by F2 
Q0: energy of gammas emitted pre-fission
Q1: energy of gammas emitted by F1
Q2: energy of gamma emitted by F2
Q0(x,y,z): directions of photons emmitted by compound nucleus
Q1(x,y,z): directions of photons emmitted by F1
Q2(x,y,z): directions of photons emmited by F2
*/

////////////////////////////////////////////////////////////
//      Fragment yields
////////////////////////////////////////////////////////////

//Fragment yield, A, at moment of fission
TCanvas *c1 = new TCanvas("c1","Fragment Yield",150,10,990,660);
mytree->Draw("iAf1>>hframe_fragyield");
mytree->Draw("iAf2>>hframe_fragyield2");

//Add fragmentyield to histogram, x1 scale
htotal_fragyield->Add(hframe_fragyield,1.0);
htotal_fragyield->Add(hframe_fragyield2,1.0);

// 
htotal_fragyield->Sumw2(); //sum of squared weights, adds uncertanity, square root of number of counts
Double_t norm = 2; //normalized to 2 -> because two distributions
Double_t scale = norm/(htotal_fragyield->Integral());
htotal_fragyield->Scale(scale); //Scale normalized histogram

htotal_fragyield->SetLineColor(23);
htotal_fragyield->SetMarkerColor(23);
htotal_fragyield->SetMarkerStyle(4);
htotal_fragyield->SetMarkerSize(0.2);
htotal_fragyield->Draw("ep");
htotal_fragyield->GetXaxis()->SetTitle("A");
htotal_fragyield->GetYaxis()->SetTitle("Yield of fragment");
htotal_fragyield->SetTitle("Fission fragment yield");
//c1->Print("../plot/FragmentYield.pdf");


///////////////////////////////////////////////////
//Neutron yields
///////////////////////////////////////////////////
 

//Total neutron yield
TCanvas *c01 = new TCanvas("c01", "Neutron Yield", 150, 10, 990, 660);
mytree->Draw("n1>>hframe_n_yield");
mytree->Draw("n2>>hframe_n_yield2");

htotal_n_yield->Add(hframe_n_yield,1.0);
htotal_n_yield->Add(hframe_n_yield2,1.0);

htotal_n_yield->Sumw2();
Double_t norm_n = 1;
Double_t scale_n = norm_n/(htotal_n_yield->Integral());
htotal_n_yield->Scale(scale_n);
htotal_n_yield->SetLineColor(23);
htotal_n_yield->SetMarkerColor(23);
htotal_n_yield->SetMarkerStyle(4);
htotal_n_yield->SetMarkerSize(0.2);
htotal_n_yield->Draw();
htotal_n_yield->GetXaxis()->SetTitle("Neutrons");
htotal_n_yield->GetYaxis()->SetTitle("Yield of neutrons");
htotal_n_yield->SetTitle("Total neutron yield");

//Gets rid of error bars in x dir
//TH1::SetDefaultSumw2();
//gStyle->SetErrorX(0.0001);


//Neutron yield for each fragment
TCanvas *c02 = new TCanvas("c02", "Neutron Yield", 150, 10, 990, 660);
mytree->Draw("n1>>hframe_n_yield");
mytree->Draw("n2>>hframe_n_yield2");

hframe_n_yield->Scale(scale_n);
hframe_n_yield2->Scale(scale_n);

hframe_n_yield->SetMarkerStyle(4);
hframe_n_yield->SetMarkerSize(1.0);
hframe_n_yield->SetMarkerColor(45);
hframe_n_yield->Draw("ep");
hframe_n_yield2->Draw("same,ep");
hframe_n_yield->GetXaxis()->SetTitle("Neutrons");
hframe_n_yield->GetYaxis()->SetTitle("Yield/Fission");
hframe_n_yield->SetTitle("Neutron yield per fission");

   auto legend1 = new TLegend(0.1,0.8,0.3,0.9);
   legend1->SetTextSize(0.02);
   // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend1->AddEntry(hframe_n_yield,"F1","ep");
   legend1->AddEntry(hframe_n_yield2,"F2","ep");
   legend1->Draw();


//Fragment yield, A, after decay (?)
TCanvas *c2 = new TCanvas("c2","Product Yield",150,10,990,660);
mytree->Draw("iAp1>>hframe_prodyield");
mytree->Draw("iAp2>>hframe_prodyield2");
htotal_prodyield->Add(hframe_prodyield,1.0);
htotal_prodyield->Add(hframe_prodyield2,1.0);

htotal_prodyield->Sumw2();
norm = 2;
scale = norm/(htotal_prodyield->Integral());
htotal_prodyield->Scale(scale);

htotal_prodyield->SetMarkerStyle(4);
htotal_prodyield->SetMarkerSize(0.2);
htotal_prodyield->SetMarkerColor(kRed);
htotal_prodyield->SetLineColor(kRed);
htotal_prodyield->Draw("ep");
htotal_fragyield->Draw("same,ep");
htotal_prodyield->GetXaxis()->SetTitle("Mass number");
htotal_prodyield->GetYaxis()->SetTitle("Yield per fission");
htotal_prodyield->GetYaxis()->SetDecimals(2);
htotal_prodyield->GetXaxis()->SetRange(60,180);
htotal_prodyield->SetTitle("Fission yield, fragments and products ");

   auto legend = new TLegend(0.1,0.8,0.3,0.9);
   legend->SetTextSize(0.02);
   // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry(htotal_fragyield,"Fragments","ep");
   legend->AddEntry(htotal_prodyield,"Products","ep");
   legend->Draw();

//c2->Print("../plot/ProdYield.pdf");



/////////////////////////////////////////////////////
// Kinetic Energies of fragments, as a function of fragment mass
/////////////////////////////////////////////////////


TCanvas *cEkin = new TCanvas("cEkin","Ekin",150,10,990,660);

// mytree->Draw("E1kin:iAf1>>hframe_Ekin1");
// hframe_Ekin1->Draw("col");
// mytree->Draw("E2kin:iAf2>>hframe_Ekin2");


mytree->Draw("E1kin:iAf1>>hframe_Ekin1","","col");
mytree->Draw("E2kin:iAf2>>hframe_Ekin2","","col");
h_Ekin->Add(hframe_Ekin1,1.0);
h_Ekin->Add(hframe_Ekin2,1.0);
h_Ekin->Draw("col");

//Find mean fragment energy
const Int_t minMass = 0;
const Int_t maxMass = 200;

Int_t nPoints = 0;
Int_t bin;
int i=0;

Double_t mean_arr[maxMass];
Double_t rms_arr[maxMass];
Double_t Af_arr[maxMass];

Double_t mean;

for(int i=minMass; i<maxMass; i++)
{
    auto temp = h_Ekin->ProjectionY("temp",i,i+1);
    mean = temp->GetMean();

    if(mean<=0){continue;}
    Af_arr[nPoints] = i;
    mean_arr[nPoints]   = temp->GetMean();
    rms_arr[nPoints]= temp->GetRMS();
    nPoints ++;
} 


//Error bars
// TGraphErrors* gr = new TGraphErrors(noMasses+noMasses2,mass_arr,mean,0,rms_stuff);
TGraphErrors* gr = new TGraphErrors(nPoints,Af_arr,mean_arr,0,rms_arr);
gr->GetXaxis()->SetTitle("Mass number");
gr->GetYaxis()->SetTitle("Mean fragment energy [MeV]");
gr->SetTitle("");
gr -> GetXaxis() -> SetRangeUser(70,180);
gr -> GetXaxis() -> CenterTitle();
gr -> GetXaxis()->SetTitleOffset(1.3);
gr -> GetYaxis() -> CenterTitle();
gr -> GetYaxis()->SetTitleOffset(1.3);
gr -> SetTitle("Mean fragment kinetic energy as a function of fragment mass number");
// gr->SetMarkerSize(2);
gr->SetMarkerStyle(22);
gr->Draw("AeP");
//cEkin->SaveAs("../plot/mean_energy_vs_mass.pdf");


/////////////////////////////////////////////////////
// Neutron Multiplicities
/////////////////////////////////////////////////////



//Total neutron multiplicity
//Probabilty that a given nu occurs. eg. probability that 3 neutrons are emitted

TCanvas *c3 = new TCanvas("c3","Neutron Multiplicities",150,10,990,660);
mytree->Draw("n0>>hframe_n0");
mean = hframe_n0->GetMean();
cout << "\n" << endl;
cout << "nu_bar_PreFiss: " << mean << endl;
Double_t uncertainty_n0 = sqrt(mean*F)/F;
cout << "Uncertainty n0: " << uncertainty_n0 << "\n" << endl;

mytree->Draw("n1>>hframe_n1");
// hframe_n1->SaveAs("../plot/n_mult1.pdf");
//cout << "Total n1: " << total_n1 << endl; 
mean = hframe_n1->GetMean();
cout << "nu_bar_FF1: " << mean << endl;
Double_t uncertainty_n1 = sqrt(mean*F)/F;
cout << "Uncertainty n1: " << uncertainty_n1 << "\n" << endl;

mytree->Draw("n2>>hframe_n2");
// hframe_n2->SaveAs("../plot/n_mult2.pdf");
mean = hframe_n2->GetMean();
cout << "nu_bar_FF2: " << mean << endl;
Double_t uncertainty_n2 = sqrt(mean*F)/F;
cout << "Uncertainty n2: " << uncertainty_n2 << "\n" << endl;

mytree->Draw("n1:n2>>hframe_n_multi","","col");
mytree->Draw("n0:n2:n1>>hframe_n_multi3D","","lego");

//number of neutrons

//To compensate for that 0 can be emitted by the first, and 1 by the second, and then multiplicity is 1
Int_t lastbin = h_n_mult_total->GetNbinsX();
Double_t n_n0;
Double_t n_n1;
Double_t n_n2;
Double_t n_tot;
Double_t weight;
for (int k=0;k<lastbin;k++){
n_n0 = k-1;
hframe_n_multi3D->GetZaxis()->SetRange(k,k);
TH2D* pxy = (TH2D*) hframe_n_multi3D->Project3D("xy");
pxy->Draw("col");
  for (int i=0;i<lastbin;i++){
    n_n1 = i-1; // bin 1 contains evetns of multiplicity 0
    auto px = pxy->ProjectionX("px",i,i);
      for (int j=0;j<lastbin;j++){
      weight = px->GetBinContent(j);
      n_n2 = j-1; //bin 1 contains evetns of multiplicity 0
      if (n_n0<0){n_n0=0;} // bin 0 events (underflow) have to come in underflow
      if (n_n1<0){n_n1=0;} // bin 0 events (underflow) have to come in underflow
      if (n_n2<0){n_n2=0;} // bin 0 events (underflow) have to come in underflow
      n_tot = n_n0+n_n1+n_n2;
      h_n_mult_total->Fill(n_tot,weight);
      }
  }
}

// Int_t lastbin = h_n_mult_total->GetNbinsX();
// Double_t n_n1;
// Double_t n_n2;
// n_tot;
// Double_t weight;
// for (int i=0;i<lastbin;i++){
//   n_n1 = i-1; // bin 1 contains evetns of multiplicity 0
//   auto px = hframe_n_multi->ProjectionX("px",i,i);
//     for (int j=0;j<lastbin;j++){
//     weight = px->GetBinContent(j);
//     n_n2 = j-1; //bin 1 contains events of multiplicity 0
//     n_tot = n_n1+n_n2;
//     if (n_tot<0){n_tot=0;} // bin 0 events (underflow) have to come in underflow
//     h_n_mult_total->Fill(n_tot,weight);
//     }
// }




// normalize the multiplicites
norm = h_n_mult_total->Integral();
Double_t n_tot_2 = 0;
for(int i=2;i<nbins_h_n_mult_total+1;i++){
  n_tot_2 += h_n_mult_total->GetBinContent(i);
} 

h_n_mult_total->Scale(1/norm);

// factorial moments for nu_total
Double_t value;
Double_t moments[5];
std::fill_n(moments, 5, 0.);
// value = h_n_mult_total->GetBinContent(2);
// cout << value << endl;

// bin zero is underflow bin!
int nu;
//h_n_mult_total is P(nu), as it is already scaled!
Double_t unc_n_moments[5];
std::fill_n(unc_n_moments, 5, 0);

for(int i=1;i<nbins_h_n_mult_total+1;i++){
  nu = i-1; // bin 1 is 0 neutrons ...
  //value is P(given nu), for example P(nu=1), as it is the first bin in h_n_mult_total
  value = h_n_mult_total->GetBinContent(i);

  //moments[0] = 0; // dummy
  moments[1] += nu * value;
  moments[2] += nu * (nu-1) * value;
  moments[3] += nu * (nu-1) * (nu-2) * value;
  moments[4] += nu * (nu-1) * (nu-2) * (nu-3) * value;

  unc_n_moments[1] += value*F*nu;
  unc_n_moments[2] += value*F*nu*(nu-1);
  unc_n_moments[3] += value*F*nu*(nu-1)*(nu-2);
  unc_n_moments[4] += value*F*nu*(nu-1)*(nu-2)*(nu-3);

  // value = h_n_mult_total->GetBinContent(i);

}
//assume number of neutrons = moments[1]*number of fissions, because I don't find the correct value elsewhere
Double_t total_nr_neutrons = moments[1]*F;
Double_t uncertanity_n_avg = sqrt(total_nr_neutrons)/F;
cout << "Uncertainty in nu_avg: " << uncertanity_n_avg << endl;


Double_t n_multiplicity = moments[1];

cout << "\n Neutron Moments for nu" << endl;
// write to terminal
for(int i=1;i<5;i++){
  cout << "Moment" << i << " " << moments[i] << " Uncertainty: " << sqrt(unc_n_moments[i])/F << endl;
}


//factorial moments for nu_0
Double_t value_0;
Double_t moments_0[5];

std::fill_n(moments_0, 5, 0.);

Double_t unc_n0_moments[5];
std::fill_n(unc_n0_moments, 5, 0);

Double_t norm_0 = hframe_n0->Integral();
hframe_n0->Scale(1/norm_0);

int nu_0;
for(int i=1; i<nbins_h_n_mult_total+1;i++){
  nu_0 = i-1;

  value_0 = hframe_n0->GetBinContent(i);

  moments_0[1] += nu_0 * value_0;
  moments_0[2] += nu_0 * (nu_0-1) * value_0;
  moments_0[3] += nu_0 * (nu_0-1) * (nu_0-2) * value_0;
  moments_0[4] += nu_0 * (nu_0-1) * (nu_0-2) * (nu_0-3) * value_0;

  unc_n0_moments[1] += value_0*F*nu_0;
  unc_n0_moments[2] += value_0*F*nu_0*(nu_0-1);
  unc_n0_moments[3] += value_0*F*nu_0*(nu_0-1)*(nu_0-2);
  unc_n0_moments[4] += value_0*F*nu_0*(nu_0-1)*(nu_0-2)*(nu_0-3);
}


cout << "\n Neutron Moments for nu_0" << endl;
// write to terminal
for(int i=1;i<5;i++){
  cout << "Moment" << i << " " << moments_0[i] << " Uncertainty: " << sqrt(unc_n0_moments[i])/F << endl;
}


//factorial moments for nu_1
Double_t value_1;
Double_t moments_1[5];

Double_t norm_1 = hframe_n1->Integral();
hframe_n1->Scale(1/norm_1);

std::fill_n(moments_1, 5, 0.);

Double_t unc_n1_moments[5];
std::fill_n(unc_n1_moments, 5, 0);

int nu_1;
for(int i=1; i<nbins_h_n_mult_total+1;i++){
  nu_1 = i-1;

  value_1 = hframe_n1->GetBinContent(i);

  moments_1[1] += nu_1 * value_1;
  moments_1[2] += nu_1 * (nu_1-1) * value_1;
  moments_1[3] += nu_1 * (nu_1-1) * (nu_1-2) * value_1;
  moments_1[4] += nu_1 * (nu_1-1) * (nu_1-2) * (nu_1-3) * value_1;

  unc_n1_moments[1] += value_1*F*nu_1;
  unc_n1_moments[2] += value_1*F*nu_1*(nu_1-1);
  unc_n1_moments[3] += value_1*F*nu_1*(nu_1-1)*(nu_1-2);
  unc_n1_moments[4] += value_1*F*nu_1*(nu_1-1)*(nu_1-2)*(nu_1-3);
}


cout << "\n Neutron Moments for nu_1" << endl;
// write to terminal
for(int i=1;i<5;i++){
  cout << "Moment" << i << " " << moments_1[i] << " Uncertainty: " << sqrt(unc_n1_moments[i])/F << endl;
}

//factorial moments for nu_2
Double_t value_2;
Double_t moments_2[5];

std::fill_n(moments_2, 5, 0.);

Double_t unc_n2_moments[5];
std::fill_n(unc_n2_moments, 5, 0);

Double_t norm_2 = hframe_n2->Integral();
hframe_n2->Scale(1/norm_2);

int nu_2;
for(int i=1; i<nbins_h_n_mult_total+1;i++){
  nu_2 = i-1;

  value_2 = hframe_n2->GetBinContent(i);

  moments_2[1] += nu_2 * value_2;
  moments_2[2] += nu_2 * (nu_2-1) * value_2;
  moments_2[3] += nu_2 * (nu_2-1) * (nu_2-2) * value_2;
  moments_2[4] += nu_2 * (nu_2-1) * (nu_2-2) * (nu_2-3) * value_2;

  unc_n2_moments[1] += value_2*F*nu_2;
  unc_n2_moments[2] += value_2*F*nu_2*(nu_2-1);
  unc_n2_moments[3] += value_2*F*nu_2*(nu_2-1)*(nu_2-2);
  unc_n2_moments[4] += value_2*F*nu_2*(nu_2-1)*(nu_2-2)*(nu_2-3);
}


cout << "\n Neutron Moments for nu_2" << endl;
// write to terminal
for(int i=1;i<5;i++){
  cout << "Moment" << i << " " << moments_2[i] << " Uncertainty: " << sqrt(unc_n2_moments[i])/F << endl;
}

Double_t n_mult_variance;
Double_t n_mult_error;
// calculate the variance
n_mult_variance = moments[2] -  moments[1] * (moments[1] -1 );
n_mult_error = sqrt(n_mult_variance);
cout << "Std. Dev (?)" << n_mult_error << endl;

// x axis is such that [0,1] is multiplicity 0, ...
h_n_mult_total->GetXaxis()->SetTitle("Neutron multiplicity nu");
h_n_mult_total->GetYaxis()->SetTitle("P(nu)");
h_n_mult_total->SetTitle("Total neutron multiplicity");
h_n_mult_total->GetYaxis()->SetDecimals(3);
h_n_mult_total->GetXaxis()->SetRange(1,10);
h_n_mult_total->SetMarkerColor(4);
h_n_mult_total->SetMarkerStyle(8);
// h_n_mult_total->GetXaxis()->SetRangeUser(0.,nbins_h_n_mult_total);
h_n_mult_total->Draw("hist p");
c3->Update();
//c3->Print("../plot/n_mult.pdf");


///////////////
//Neutron multiplicity per fragment
TCanvas *cMults = new TCanvas("cMults","Neutron Multiplicities(2)",150,10,990,660);

norm = hframe_n2->Integral();
hframe_n2->SetMarkerStyle(8);
hframe_n2->SetMarkerColor(8);
hframe_n2->SetLineColor(8);
hframe_n2->Scale(1/norm);
hframe_n2->GetXaxis()->SetTitle("multiplicity nu");
hframe_n2->GetYaxis()->SetTitle("P(nu)");
hframe_n2->SetTitle("Neutron multiplicity, per fragment");
hframe_n2->GetXaxis()->SetRangeUser(0,10);
hframe_n2->GetYaxis()->SetRangeUser(0,0.6);
hframe_n2->Draw("hist pl");

norm = hframe_n1->Integral();
hframe_n1->Scale(1/norm);
hframe_n1->SetMarkerColor(2);
hframe_n1->SetLineColor(2);
hframe_n1->SetMarkerStyle(8);
//hframe_n1->SetTitle("Neutron multiplicity, per fragment");
hframe_n1->Draw("same, hist pl");

norm = hframe_n0->Integral();
hframe_n0->Scale(1/norm);
hframe_n0->SetMarkerColor(7);
hframe_n0->SetLineColor(7);
hframe_n0->SetMarkerStyle(8);

//hframe_n0->SetLineColor(3);
//hframe_n0->SetTitle("Neutron multiplicity, per fragment");
hframe_n0->Draw("same, hist pl");

//Plot total n multiplicity in same plot as above
h_n_mult_total->Draw("same, hist pl");

legend = new TLegend(0.7,0.75,0.9,0.9);
legend->SetTextSize(0.03);
// legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
legend->AddEntry(hframe_n0,     "Prefission");
legend->AddEntry(hframe_n1,     "Fragment 1");
legend->AddEntry(hframe_n2,     "Fragment 2");
legend->AddEntry(h_n_mult_total,"Total");
legend->Draw();

//cMults->SaveAs("../plot/n_mult_together.pdf");
;

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
    n_m1 = i-1; // bin 1 contains events of multiplicity 0
    auto px = pxy_m->ProjectionX("px",i,i);
      for (int j=0;j<lastbin1;j++){
      weight_m = px->GetBinContent(j);
      n_m2 = j-1; //bin 1 contains events of multiplicity 0
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

Double_t uncertainty_ph_avg = sqrt(p_multiplicity*F)/F;
cout << "Uncertainty photon mult: " << uncertainty_ph_avg << endl;

cout << "\n Photons Moments" << endl;
// write to terminal
for(int i=1;i<5;i++){
cout << "Moment" << i << " " << moments[i] << endl;
}

// x axis is such that [0,1] is multiplicity 0, ...
h_p_mult_total->GetXaxis()->SetTitle("Photon multiplicity PM");
h_p_mult_total->GetYaxis()->SetTitle("P(PM)");
h_p_mult_total->SetTitle("Total photon multiplicity");
h_p_mult_total->GetYaxis()->SetDecimals(3);
h_p_mult_total->GetXaxis()->SetRange(0,21);
h_p_mult_total->SetMarkerColor(4);
h_p_mult_total->SetMarkerStyle(8);
// h_n_mult_total->GetXaxis()->SetRangeUser(0.,nbins_h_n_mult_total);
h_p_mult_total->Draw("hist p");
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

legend = new TLegend(0.7,0.75,0.9,0.9);
legend->SetTextSize(0.03);
// legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
legend->AddEntry(hframe_p_mult_0,     "Prefission","l");
legend->AddEntry(hframe_p_mult_1,     "Fragment 1","l");
legend->AddEntry(hframe_p_mult_2,     "Fragment 2","l");
legend->AddEntry(h_p_mult_total,"Total","l");
legend->Draw();



/////////////////////////////////////////////////////
// Neutron Spectrum, neutron spectral shape
/////////////////////////////////////////////////////

TCanvas *c5 = new TCanvas("c5","Neutron Spectrum (spectral shape)",150,10,990,660);
mytree->Draw("P0>>hframe_n_E_0");
mytree->Draw("P1>>hframe_n_E_1");
mytree->Draw("P2>>hframe_n_E_2");
h_n_E_total->Add(hframe_n_E_0,1.0);
h_n_E_total->Add(hframe_n_E_1,1.0);
h_n_E_total->Add(hframe_n_E_2,1.0);


// at the moment if no neutron is emitted, this is tabulated in bin 0
cout << "Workaround as long as underflow bin is /WRONG/" << endl;
h_n_E_total->SetBinContent(1,0);
//

// normalize the multiplicites
h_n_E_total->Sumw2();
norm = h_n_E_total->Integral();
// norm = norm * n_multiplicity;
norm /= max_h_n_E;
h_n_E_total->Scale(1./norm);


// Boltzmann 
cout << "\n Neutron Spectrum" << endl;
TF1 *f_n_Boltz = new TF1("f_n_Boltz",n_Boltzman_fct,0,max_h_n_E,/*dimParams*/ 2);
f_n_Boltz->SetParNames("T","binWidth");
f_n_Boltz->SetParameter(0,2.); // Start value
f_n_Boltz->FixParameter(1,h_n_E_total->GetBinWidth(0)*max_h_n_E); // Default value
h_n_E_total->Fit("f_n_Boltz");
// Double_t chi2 = f_n_Boltz->GetChisquare();
// Double_t n_fitted_parameters = 2.;
// Double_t n_freedom = h_n_E_total->GetEntries() - n_fitted_parameters;
// Double_t chi2_red = chi2 / n_freedom;
// cout << "chi2: " << chi2 << endl;
// cout << "chi-red: " << chi2_red << endl;

c5->SetLogy();
h_n_E_total->GetXaxis()->SetTitle("Neutron Energy En [MeV]");
h_n_E_total->GetYaxis()->SetTitle("Neutrons per unit energy and nubar [MeV^(-1)]");
h_n_E_total->SetTitle("Neutron spectral shape");
h_n_E_total->GetYaxis()->SetDecimals(2);
h_n_E_total->GetYaxis()->SetRange(0,1);
//h_n_E_total->Draw("E");
//h_n_E_total->GetYaxis()->SetRange(0,1000);

   auto legend3 = new TLegend(0.7,0.8,0.3,0.9);
   legend3->SetTextSize(0.03);
   // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend3->AddEntry(f_n_Boltz,"Boltzmann distribution");
   legend3->AddEntry(h_n_E_total,"Neutron spectral shape");
   legend3->Draw();

//c5->Print("../plot/n_spectrum.pdf");


// //////////////////////////////////////////////////////
//Photon spectrum, photon spectral shape
// //////////////////////////////////////////////////////
TCanvas *c6 = new TCanvas("c6","Photon Spectrum (spectral shape)",150,10,990,660);
mytree->Draw("Q0>>hframe_ph_E_0");
mytree->Draw("Q1>>hframe_ph_E_1");
mytree->Draw("Q2>>hframe_ph_E_2");
h_ph_E_total->Add(hframe_ph_E_0,1.0);
h_ph_E_total->Add(hframe_ph_E_1,1.0);
h_ph_E_total->Add(hframe_ph_E_2,1.0);

// at the moment if no photon is emitted, this is tabulated in bin 0
cout << "Workaround as long as underflow bin is /WRONG/" << endl;
h_ph_E_total->SetBinContent(1,0);
//

// normalize the multiplicites
h_ph_E_total->Sumw2();
norm = h_ph_E_total->Integral();
// norm = norm * n_multiplicity;
norm /= max_h_ph_E;
h_ph_E_total->Scale(1./norm);

// Boltzmann 
cout << "\n Photon Spectrum" << endl;
TF1 *f_ph_Boltz = new TF1("f_ph_Boltz",ph_Boltzman_fct,0,max_h_ph_E,/*dimParams*/ 2);
f_ph_Boltz->SetParNames("T","binWidth");
f_ph_Boltz->SetParameter(0,2.); // Start value
f_ph_Boltz->FixParameter(1,h_ph_E_total->GetBinWidth(0)*max_h_ph_E); // Default value

h_ph_E_total->Fit("f_ph_Boltz");
// Double_t chi2 = f_n_Boltz->GetChisquare();
// Double_t n_fitted_parameters = 2.;
// Double_t n_freedom = h_n_E_total->GetEntries() - n_fitted_parameters;
// Double_t chi2_red = chi2 / n_freedom;
// cout << "chi2: " << chi2 << endl;
// cout << "chi-red: " << chi2_red << endl;

c6->SetLogy();
h_ph_E_total->GetXaxis()->SetTitle("Photon Energy En [MeV]");
h_ph_E_total->GetYaxis()->SetTitle("Photons per unit energy and mean_nr_of_photons [MeV^(-1)]");
h_ph_E_total->SetTitle("Photon spectral shape");
h_ph_E_total->GetYaxis()->SetDecimals(2);
h_ph_E_total->Draw("E");


   auto legend4 = new TLegend(0.7,0.8,0.3,0.9);
   legend4->SetTextSize(0.03);
   // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend4->AddEntry(f_n_Boltz,"Boltzmann distribution");
   legend4->AddEntry(h_ph_E_total,"Photon spectral shape");
   legend4->Draw();

/////////////////////////////////
// Neutron angular correlation
/////////////////////////////////



}

void create_frames() {
gStyle->SetOptStat(0);
int nbins;
int maxbin;

hframe_fragyield = new TH1D("hframe_fragyield","",250,0,249);
hframe_fragyield2 = new TH1D("hframe_fragyield2","",250,0,249);
htotal_fragyield = new TH1D("htotal_fragyield","",250,0,249);

//Numbers: #bins, startbin, endbin
hframe_n_yield = new TH1D("hframe_n_yield","",10,0,9);
hframe_n_yield2 = new TH1D("hframe_n_yield2","",10,0,9);
htotal_n_yield = new TH1D("htotal_n_yield","",10,0,9);

hframe_prodyield = new TH1D("hframe_prodyield","",250,0,249);
hframe_prodyield2 = new TH1D("hframe_prodyield2","",250,0,249);
htotal_prodyield = new TH1D("htotal_prodyield","",250,0,249);

hframe_Ekin1 = new TH2D("hframe_Ekin1","",200,0,199,250,0,249);
hframe_Ekin2 = new TH2D("hframe_Ekin2","",200,0,199,250,0,249);
h_Ekin = new TH2D("h_Ekin","",200,0,199,250,0,249);

nbins = nbins_h_n_mult_total;
maxbin = nbins;
hframe_n0 = new TH1D("hframe_n0","",nbins,-0.5,maxbin-0.5);
hframe_n1 = new TH1D("hframe_n1","",nbins,-0.5,maxbin-0.5);
hframe_n2 = new TH1D("hframe_n2","",nbins,-0.5,maxbin-0.5);
hframe_n_multi  = new TH2D("hframe_n_multi","",nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5);
hframe_n_multi3D  = new TH3F("hframe_n_multi3D","",nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5);
h_n_mult_total = new TH1D("h_n_mult_total","",nbins,-0.5,maxbin-.5);




nbins = nbins_h_p_mult_total;
maxbin = nbins;
hframe_p_mult_0 = new TH1D("hframe_p_mult_0","",nbins,-0.5,maxbin-.5);
hframe_p_mult_1 = new TH1D("hframe_p_mult_1","",nbins,-0.5,maxbin-.5);
hframe_p_mult_2 = new TH1D("hframe_p_mult_2","",nbins,-0.5,maxbin-0.5);
hframe_p_multi  = new TH2D("hframe_p_multi","",nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5);
hframe_p_multi3D  = new TH3F("hframe_p_multi3D","",nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5);
h_p_mult_total = new TH1D("h_p_mult_total","",nbins,-0.5,maxbin-.5);

nbins = nbins_h_n_E_total;
maxbin = max_h_n_E;
hframe_n_E_0 = new TH1D("hframe_n_E_0","",nbins,0,maxbin);
hframe_n_E_1 = new TH1D("hframe_n_E_1","",nbins,0,maxbin);
hframe_n_E_2 = new TH1D("hframe_n_E_2","",nbins,0,maxbin);
h_n_E_total = new TH1D("h_n_E_total","",nbins,0,maxbin);
// h_n_E_Boltzmann = new TH1D("h_n_E_Boltzmann","",nbins,0,maxbin);


nbins = nbins_h_ph_E_total;
maxbin = max_h_ph_E;
hframe_ph_E_0 = new TH1D("hframe_ph_E_0","",nbins,0,maxbin);
hframe_ph_E_1 = new TH1D("hframe_ph_E_1","",nbins,0,maxbin);
hframe_ph_E_2 = new TH1D("hframe_ph_E_2","",nbins,0,maxbin);
h_ph_E_total = new TH1D("h_ph_E_total","",nbins,0,maxbin);
// h_n_E_Boltzmann = new TH1D("h_n_E_Boltzmann","",nbins,0,maxbin);


// hframe_kin = new TH1D("hframe_kin","",140,0,139);
// hframe_kin2 = new TH1D("hframe_kin2","",140,0,139);
// htotal_kin = new TH1D("htotal_kin","",140,0,139);


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