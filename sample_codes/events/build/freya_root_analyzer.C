#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TROOT.h>
#include <iostream>
#include <cstdio>
#include "TFile.h"
#include "TChain.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

//////////////////////////////////////////////////////////////
// To be used in Dorthea's master, to analyze FREYA output files
//
// Includes:
// - mean number of photons emitted
// - mean photon energy
// - total photon energy, per fragment and for all three fragments
// - writes these quanties to a file data_as_func_of_excitation_energy.dat


// OBS: the error calculations for average and total gamma energies is not _entirely_ correct. 
// As the total number of gammas and the gamma multiplicities for certain gamma energies is correlated,
// the error term should include both error in total number of gammas AND the covariance
// See Master Journal pages 16-18 for more info


//OBS: Root places energies in bins. Therefore, the exact energies are lost, how much depends on the size of the bin
/////////////////////////////////////////////////////////////////


TFile *vetsex = new TFile("Pu240.dat.root", "READ");
TTree *mytree = (TTree *) gROOT->FindObject("FreyaTree");

//TTree* mytree_cut = mytree->CopyTree("Entry$<1000");

//
Double_t pi = 3.14159265359;

//need number of fissions
TH1D *hframe_fragyield;
TH1D *hframe_fragyield2;
TH1D *htotal_fragyield;

//Photon emission multiplicities
TH1D *hframe_p_mult_0;
TH1D *hframe_p_mult_1;
TH1D *hframe_p_mult_2;
TH2D *hframe_p_multi;
TH3F *hframe_p_multi3D;
TH1D *h_p_mult_total;
TH1D *hframe_p_mult_first;
TH1D *hframe_p_mult_second;
TH1D *hframe_p_mult_third;
int nbins_h_p_mult_total  = 20;

//Neutron emission multiplicities
TH1D *hframe_n_mult_first;
TH1D *hframe_n_mult_second;
TH1D *hframe_n_mult_third;

TH1D *hframe_ph_E_0;
TH1D *hframe_ph_E_1;
TH1D *hframe_ph_E_2;
TH1D *h_ph_E_total;
TH1D *hframe_n_E_1;
TH1D *hframe_n_E_2;
TH1D *h_n_E_total;
TH1D *hframe_ph_E_first;
TH1D *hframe_ph_E_second;
TH1D *hframe_ph_E_third;
TH1D *hframe_n_E_first;
TH1D *hframe_n_E_second;
TH1D *hframe_n_E_third;
TH1D *hframe_fragkin1;
TH1D *hframe_fragkin2;
// TH1D *h_n_E_Boltzmann;
int nbins_h_ph_E_total  = 500;
int max_h_ph_E  = 7; // in MeV

void create_frames();

/*
List of meanings:

iAf1: fission fragment 1, mass number
iAf2: fission fragment 2, mass number
iAp1: product fission fragment 1, mass number (after initial fragment has decayed)
iAp2: product fission fragment 2, mass number (after initial fragment has decayed)
E1kin: kinetic energy of 1 product
E2kin: kinetic energy of second product
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

m_first: # photons from first chance fission
m_second: # photons from second chance fission
m_second: # photons from third chance fission

Q_first: energy from PFG's from first chance
Q_second: energy from PFG's from second chance
Q_third: energy from PFG's from third chance
*/


void freya_root_analyzer() {

create_frames();
//Fission fragment distribution
TCanvas *c1 = new TCanvas("c1","Fragment Yield",150,10,990,660);
mytree->Draw("iAf1>>hframe_fragyield");
hframe_fragyield->SetLineColor(1);
mytree->Draw("iAf2>>hframe_fragyield2");
hframe_fragyield->SetLineColor(2);
hframe_fragyield->GetXaxis()->SetTitle("Mass [A]");
hframe_fragyield->GetYaxis()->SetTitle("Counts");
hframe_fragyield->SetTitle("Fission fragment mass distribution");
hframe_fragyield->Draw();
hframe_fragyield2->Draw("same");

auto legend_c1 = new TLegend(0.7,0.75,0.9,0.9);
legend_c1->SetTextSize(0.03);
// legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
legend_c1->AddEntry(hframe_fragyield,     "FF1, light,","l");
legend_c1->AddEntry(hframe_fragyield2,     "FF2, heavy","l");
legend_c1->Draw();


cout << "Initial mass FF1 [A]: " << hframe_fragyield->GetMean() << " Initial mass FF2 [A]: " << hframe_fragyield2->GetMean() << endl;
cout << "\n" << endl;

std::ofstream ofs20;
ofs20.open ("fragment_mass_distr.dat", std::ofstream::out | std::ofstream::app);
ofs20 << hframe_fragyield->GetMean() << "  " << hframe_fragyield2->GetMean() << endl;
ofs20.close();


//Find average kinetic energy of the two fission fragments
TCanvas *c18 = new TCanvas("c18", "Fragment kinetic energy",150,10,990,660);
mytree->Draw("E1kin>>hframe_fragkin1");
mytree->Draw("E2kin>>hframe_fragkin2");
hframe_fragkin1->SetLineColor(1);
hframe_fragkin2->SetLineColor(2);
//hframe_fragkin2->GetXaxis()->SetTitle("Mass [A]");
//hframe_fragkin2->GetYaxis()->SetTitle("Counts");
hframe_fragkin1->GetYaxis()->SetRangeUser(0,10000);
hframe_fragkin1->GetXaxis()->SetTitle("Energy [MeV]");
hframe_fragkin1->GetYaxis()->SetTitle("Counts");
hframe_fragkin1->SetTitle("Fission fragment kinE distribution");
hframe_fragkin1->Draw();
hframe_fragkin2->Draw("same");

auto legend_18 = new TLegend(0.7,0.75,0.9,0.9);
legend_18->SetTextSize(0.03);
// legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
legend_18->AddEntry(hframe_fragkin1,     "FF1, light","l");
legend_18->AddEntry(hframe_fragkin2,     "FF2, heavy","l");
legend_18->Draw();

std::ofstream ofs18;
ofs18.open ("fragment_kinE.dat", std::ofstream::out | std::ofstream::app);
ofs18 << hframe_fragkin1->GetMean() << "  " << hframe_fragkin2->GetMean() << endl;
ofs18.close();


///////////////////////////////////////
int sum_F = 0;

for(int i=0;i<300;i++){
  sum_F += hframe_fragyield->GetBinContent(i);
}

cout << "Sum: " << sum_F << endl;

int F = sum_F; //Number of fissions 

Double_t mean;
Double_t norm = 2;
Double_t nu;
Double_t value;


// /////////////////////////////////////////
// // Neutron multiplicities
// /////////////////////////////////////////
// //Don't use this to calculate anything: because of fact that I have to put variables to 0  in order for my mac to not put them to anything else, the value 0 is sent to the program, Must correct manually

// TCanvas *c17 = new TCanvas("c17", "Neutron multiplicities",150,10,990,660);
// mytree->Draw("n_first>>hframe_n_mult_first");
// mytree->Draw("n_second>>hframe_n_mult_second");
// hframe_n_mult_first->SetLineColor(1);
// //hframe_n_mult_first->GetYaxis()->SetRangeUser(0,350000);
// hframe_n_mult_first->Draw();
// hframe_n_mult_first->SetTitle("Neutron multiplicities");
// hframe_n_mult_second->SetLineColor(2);
// hframe_n_mult_second->Draw("same");
// //mytree->Draw("n_third>>hframe_n_mult_third");
// //hframe_n_mult_third->Draw("same");

// auto legend_17 = new TLegend(0.7,0.75,0.9,0.9);
// legend_17->SetTextSize(0.03);
// // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
// legend_17->AddEntry(hframe_n_mult_first,     "First","l");
// legend_17->AddEntry(hframe_n_mult_second,     "Second","l");
// //legend_17->AddEntry(hframe_p_mult_third,     "Third","l");
// legend_17->Draw();

// //cout << "N first: " << hframe_n_mult_first->GetMean() << " N second" << hframe_n_mult_second->GetMean() << " N third" << hframe_n_mult_third->GetMean() << endl;

// int sum_n_first = 0;
// for(int n=0;n<20;n++){
//   sum_n_first += hframe_n_mult_first->GetBinCenter(n)*hframe_n_mult_first->GetBinContent(n);
// }


// //cout << "n_first: " << sum_n_first << endl;

// //////////////////////////////////////////////////////
// //Photon multiplicities
// //////////////////////////////////////////////////////

// //First vs second chance fission
// //Debugged this, yields the correct multiplicity
// TCanvas *c20 = new TCanvas("c20", "Photon multiplicities",150,10,990,660);
// mytree->Draw("m_first>>hframe_p_mult_first");
// hframe_p_mult_first->SetTitle("Photon multiplicities");
// hframe_p_mult_first->Draw();
// mytree->Draw("m_second>>hframe_p_mult_second");
// hframe_p_mult_second->Draw();
// //mytree->Draw("m_third>>hframe_p_mult_third");
// //hframe_p_mult_third->Draw("same");

// int sum_m_first = 0;
// for(int n=0;n<20;n++){
//   sum_m_first += hframe_p_mult_first->GetBinCenter(n)*hframe_p_mult_first->GetBinContent(n);
// }

// cout << "Integral m_first: " << sum_m_first << endl;

// hframe_p_mult_first->SetLineColor(1);
// hframe_p_mult_first->Draw();
// hframe_p_mult_second->SetLineColor(2);
// hframe_p_mult_second->Draw("same");
// //hframe_p_mult_third->SetLineColor(3);
// //hframe_p_mult_third->SetBinContent(1,0);
// //hframe_p_mult_third->Draw("same");

// auto legend_20 = new TLegend(0.7,0.75,0.9,0.9);
// legend_20->SetTextSize(0.03);
// // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
// legend_20->AddEntry(hframe_p_mult_first,     "First","l");
// legend_20->AddEntry(hframe_p_mult_second,     "Second","l");
// legend_20->AddEntry(hframe_p_mult_third,     "Third","l");
// legend_20->Draw();

/////////////////////////////////////////////////////////////////////////

//Total photon multiplicity
TCanvas *c4 = new TCanvas("c4", "Photon multiplicities",150,10,990,660);

cout << "\n" << endl;
cout << "MEAN GAMMMA MULTIPLICITY" << endl;
mytree->Draw("m0>>hframe_p_mult_0");
mean = hframe_p_mult_0->GetMean();
//cout << "\n" << endl;
cout << "mean_number_of_photons_FF0: " << mean << " Uncertainty: " << sqrt(mean*F)/F << endl;

mytree->Draw("m1>>hframe_p_mult_1");
mean = hframe_p_mult_1->GetMean();
cout << "mean_number_of_photons_FF1: " << mean << " Uncertainty: " << sqrt(mean*F)/F << endl;

mytree->Draw("m2>>hframe_p_mult_2");
mean = hframe_p_mult_2->GetMean();
cout << "mean_number_of_photons_FF2: " << mean << " Uncertainty: " << sqrt(mean*F)/F << endl;

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

//Average photon multiplicity is the first moment
cout << "Average photon multiplicity: " << p_multiplicity << " Uncertainty: " << sqrt(p_multiplicity*F)/F << "\n" <<endl;

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
hframe_p_mult_2->SetTitle("Photon multiplicity, per fragment");
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
//Plot photon energy spectrum for first, second and third chance fission
TCanvas *c19 = new TCanvas("c19","Photon Energy Spectrum, Multichance fission ",150,10,990,660);
mytree->Draw("Q_first>>hframe_ph_E_first");
mytree->Draw("Q_second>>hframe_ph_E_second");
mytree->Draw("Q_third>>hframe_ph_E_third");

hframe_ph_E_first->GetXaxis()->SetTitle("Photon Energy Eg [MeV]");
hframe_ph_E_first->GetYaxis()->SetTitle("Number of Photons");
hframe_ph_E_first->SetTitle("Photon spectrum from multichance fission");
hframe_ph_E_first->SetLineColor(2);
hframe_ph_E_first->Draw();
hframe_ph_E_first->SetBinContent(1,0);
hframe_ph_E_second->SetLineColor(1);
hframe_ph_E_second->SetBinContent(1,0);
hframe_ph_E_second->Draw("same");
hframe_ph_E_third->SetLineColor(3);
c19->SetLogy();
hframe_ph_E_third->SetBinContent(1,0);
hframe_ph_E_third->Draw("same");

auto legend_19 = new TLegend(0.7,0.75,0.9,0.9);
legend_19->SetTextSize(0.03);
// legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
legend_19->AddEntry(hframe_ph_E_first,     "First","l");
legend_19->AddEntry(hframe_ph_E_second,     "Second","l");
legend_19->AddEntry(hframe_ph_E_third,     "Third","l");
legend_19->Draw();

double En;
int F_first, F_second, F_third, F_tot;

ifstream infile;
infile.open ("multichance_file_disp.dat");
infile >> En >> F_first >> F_second >> F_third >> F_tot;

Double_t p_multiplicity_pspec_first = 0;
Double_t p_total_energy_pspec_first = 0;
for(int i=0;i<nbins_h_ph_E_total+1;i++){
  p_multiplicity_pspec_first += hframe_ph_E_first->GetBinContent(i);
  p_total_energy_pspec_first += hframe_ph_E_first->GetBinContent(i)*hframe_ph_E_first->GetBinCenter(i);
}
p_multiplicity_pspec_first = p_multiplicity_pspec_first/F_first;
p_total_energy_pspec_first = p_total_energy_pspec_first/F_first;

cout << "Mg_first: " << p_multiplicity_pspec_first << " Etot_first: " << p_total_energy_pspec_first << endl;


Double_t p_multiplicity_pspec_second = 0;
Double_t p_total_energy_pspec_second = 0;
for(int i=0;i<nbins_h_ph_E_total+1;i++){
  p_multiplicity_pspec_second += hframe_ph_E_second->GetBinContent(i);
  p_total_energy_pspec_second += hframe_ph_E_second->GetBinContent(i)*hframe_ph_E_second->GetBinCenter(i);
}
p_multiplicity_pspec_second = p_multiplicity_pspec_second/F_second;
p_total_energy_pspec_second = p_total_energy_pspec_second/F_second;

cout << "Mg_second: " << p_multiplicity_pspec_second << " Etot_second: " << p_total_energy_pspec_second << endl;

Double_t p_multiplicity_pspec_third = 0;
Double_t p_total_energy_pspec_third = 0;
for(int i=0;i<nbins_h_ph_E_total+1;i++){
  p_multiplicity_pspec_third += hframe_ph_E_third->GetBinContent(i);
  p_total_energy_pspec_third += hframe_ph_E_third->GetBinContent(i)*hframe_ph_E_third->GetBinCenter(i);
}
p_multiplicity_pspec_third = p_multiplicity_pspec_third/F_third;
p_total_energy_pspec_third = p_total_energy_pspec_third/F_third;

cout << "Mg_third: " << p_multiplicity_pspec_third << " Etot_third: " << p_total_energy_pspec_third << endl;


///////////////////////////////////////////////////////////////////////////////////////////////////////////

//Plot of photon energy spectrum
TCanvas *c6 = new TCanvas("c6","Photon Energy Spectrum ",150,10,990,660);
//mytree->Draw("Q0>>hframe_ph_E_0");
mytree->Draw("Q1>>hframe_ph_E_1");
mytree->Draw("Q2>>hframe_ph_E_2");
//h_ph_E_total->Add(hframe_ph_E_0,1.0);
h_ph_E_total->Add(hframe_ph_E_1,1.0);
h_ph_E_total->Add(hframe_ph_E_2,1.0);


/*
// normalize the multiplicites
h_ph_E_total->Sumw2();
norm = h_ph_E_total->Integral();
// norm = norm * n_multiplicity;
norm /= max_h_ph_E;
h_ph_E_total->Scale(1./norm);
*/

// at the moment if no photon is emitted, this is tabulated in bin 0
//cout << "Workaround as long as underflow bin is /WRONG/" << endl;
h_ph_E_total->SetBinContent(1,0);

c6->SetLogy();
h_ph_E_total->GetXaxis()->SetTitle("Photon Energy En [MeV]");
h_ph_E_total->GetYaxis()->SetTitle("Number of Photons");
h_ph_E_total->SetTitle("Photon spectrum");
h_ph_E_total->Draw("E");
int nbins_h_ph_E_total= h_ph_E_total->GetNbinsX();


TH1F *h_ph_E_total_scaled = (TH1F*) h_ph_E_total->Clone();
//Must scale the low-energy gammas in order to reproduce experiment
// Double_t FREYA_scale = 3.7272612980187914*log(800)-18.139593352492675;
// for(int i=1;i<nbins_h_ph_E_total;i++){
//   if(h_ph_E_total_scaled->GetBinCenter(i) < 0.130){
//     h_ph_E_total_scaled->SetBinContent(i,0);
//   }

//   if(h_ph_E_total_scaled->GetBinCenter(i) < 0.450 && h_ph_E_total_scaled->GetBinCenter(i) > 0.130){
//     h_ph_E_total_scaled->SetBinContent(i,h_ph_E_total->GetBinContent(i)*(3.7272612980187914*log(h_ph_E_total_scaled->GetBinCenter(i)*1000)-18.139593352492675)/FREYA_scale);
//   }
// }

h_ph_E_total_scaled->SetLineColor(2);
h_ph_E_total_scaled->Draw("same");

Double_t p_multiplicity_pspec = 0;
Double_t p_total_energy_pspec = 0;

for(int i=0;i<nbins_h_ph_E_total+1;i++){
  p_multiplicity_pspec += h_ph_E_total_scaled->GetBinContent(i);
  p_total_energy_pspec += h_ph_E_total_scaled->GetBinContent(i)*h_ph_E_total_scaled->GetBinCenter(i);
}

//cout << "Mg: " << p_multiplicity_pspec << " Etot: " << p_total_energy_pspec << endl;

p_multiplicity_pspec = p_multiplicity_pspec/F;
p_total_energy_pspec = p_total_energy_pspec/F;

cout << "p_multiplicity_pspec: " << p_multiplicity_pspec << " " << "p_total_energy_pspec: " << p_total_energy_pspec << "  " << "p_avg_energy_pspec: " << p_total_energy_pspec/p_multiplicity_pspec << endl;; 
cout << "\n " << endl;

//std::ofstream ofs2;
//ofs2.open ("data_as_func_of_excitation_energy_scaled.dat", std::ofstream::out | std::ofstream::app);
//ofs2 << "         " << p_multiplicity_pspec << "       " << p_total_energy_pspec/p_multiplicity_pspec <<"             " << p_total_energy_pspec << endl;
//ofs2.close();


std::ofstream ofs1;
ofs1.open ("photon_spectrum.dat", std::ofstream::out | std::ofstream::app);
ofs1 << "Photon Energy [MeV]" << "    " << "Counts per fission" << endl;
for(int i=0; i<nbins_h_ph_E_total;i++)
  ofs1 << h_ph_E_total_scaled->GetBinCenter(i) << "        " << (h_ph_E_total_scaled->GetBinContent(i))/F << endl;
ofs1.close();


/////////////////////////////////////////////////////////////////
// Neutron energies
/////////////////////////////////////////////////////////////////
TCanvas *c16 = new TCanvas("c16","Neutron Energy Spectrum, Multichance fission ",150,10,990,660);
mytree->Draw("P_first>>hframe_n_E_first");
mytree->Draw("P_second>>hframe_n_E_second");
mytree->Draw("P_third>>hframe_n_E_third");
hframe_n_E_first->GetXaxis()->SetTitle("Neutron Energy Eg [MeV]");
hframe_n_E_first->GetYaxis()->SetTitle("Number of neutrons");
hframe_n_E_first->SetLineColor(2);
hframe_n_E_first->Draw();
hframe_n_E_first->SetBinContent(1,0);
hframe_n_E_second->SetLineColor(1);
hframe_n_E_second->SetBinContent(1,0);
hframe_n_E_second->Draw("same");
hframe_n_E_third->SetLineColor(3);
c19->SetLogy();
hframe_n_E_third->SetBinContent(1,0);
hframe_n_E_third->Draw("same");

Double_t n_multiplicity_pspec_first = 0;
Double_t n_total_energy_pspec_first = 0;
for(int i=0;i<2*nbins_h_ph_E_total+1;i++){
  n_multiplicity_pspec_first += hframe_n_E_first->GetBinContent(i);
  n_total_energy_pspec_first += hframe_n_E_first->GetBinContent(i)*hframe_n_E_first->GetBinCenter(i);
}
n_multiplicity_pspec_first = n_multiplicity_pspec_first/F_first;
n_total_energy_pspec_first = n_total_energy_pspec_first/F_first;

cout << "Mn_first: " << n_multiplicity_pspec_first << " Etot_n_first: " << n_total_energy_pspec_first << endl;

Double_t n_multiplicity_pspec_second = 0;
Double_t n_total_energy_pspec_second = 0;
for(int i=0;i<2*nbins_h_ph_E_total+1;i++){
  n_multiplicity_pspec_second += hframe_n_E_second->GetBinContent(i);
  n_total_energy_pspec_second += hframe_n_E_second->GetBinContent(i)*hframe_n_E_second->GetBinCenter(i);
}
n_multiplicity_pspec_second = n_multiplicity_pspec_second/F_second;
n_total_energy_pspec_second = n_total_energy_pspec_second/F_second;

cout << "Mn_second: " << n_multiplicity_pspec_second << " Etot_n_second: " << n_total_energy_pspec_second << endl;

Double_t n_multiplicity_pspec_third = 0;
Double_t n_total_energy_pspec_third = 0;
for(int i=0;i<2*nbins_h_ph_E_total+1;i++){
  n_multiplicity_pspec_third += hframe_n_E_third->GetBinContent(i);
  n_total_energy_pspec_third += hframe_n_E_third->GetBinContent(i)*hframe_n_E_third->GetBinCenter(i);
}
n_multiplicity_pspec_third = n_multiplicity_pspec_third/F_third;
n_total_energy_pspec_third = n_total_energy_pspec_third/F_third;

cout << "Mn_third: " << n_multiplicity_pspec_third << " Etot_n_third: " << n_total_energy_pspec_third << endl;


TCanvas *c15 = new TCanvas("15","Neutron Energy Spectrum ",150,10,990,660);
//mytree->Draw("Q0>>hframe_ph_E_0");
mytree->Draw("P1>>hframe_n_E_1");
mytree->Draw("P2>>hframe_n_E_2");
//h_ph_E_total->Add(hframe_ph_E_0,1.0);
h_n_E_total->Add(hframe_n_E_1,1.0);
h_n_E_total->Add(hframe_n_E_2,1.0);
h_n_E_total->SetBinContent(1,0);
c15->SetLogy();
h_n_E_total->GetXaxis()->SetTitle("Neutron Energy En [MeV]");
h_n_E_total->GetYaxis()->SetTitle("Number of neutrons");
h_n_E_total->SetTitle("Neutron spectrum");
h_n_E_total->Draw("E");
int nbins_h_n_E_total= h_n_E_total->GetNbinsX();

Double_t n_multiplicity_pspec = 0;
Double_t n_total_energy_pspec = 0;

for(int i=0;i<2*nbins_h_ph_E_total+1;i++){
  n_multiplicity_pspec += h_n_E_total->GetBinContent(i);
  n_total_energy_pspec += h_n_E_total->GetBinContent(i)*h_n_E_total->GetBinCenter(i);
}

//cout << "Mg: " << p_multiplicity_pspec << " Etot: " << p_total_energy_pspec << endl;

n_multiplicity_pspec = n_multiplicity_pspec/F;
n_total_energy_pspec = n_total_energy_pspec/F;

cout << "n_multiplicity_pspec: " << n_multiplicity_pspec << " " << "n_total_energy_pspec: " << n_total_energy_pspec << "  " << "n_avg_energy_pspec: " << n_total_energy_pspec/n_multiplicity_pspec << endl;; 


cout << "\n" << endl;

std::ofstream ofs16;
ofs16.open ("neutron_mult.dat", std::ofstream::out | std::ofstream::app);
ofs16 << "         " << n_multiplicity_pspec <<  "       " << n_total_energy_pspec/p_multiplicity_pspec << "       "<< n_total_energy_pspec;
ofs16 << "         " << n_multiplicity_pspec_first  << "       " << n_total_energy_pspec_first/p_multiplicity_pspec_first << "        " << n_total_energy_pspec_first;
ofs16 << "         " << n_multiplicity_pspec_second << "       " << n_total_energy_pspec_second/p_multiplicity_pspec_second << "        " << n_total_energy_pspec_second;
ofs16 << "         " << n_multiplicity_pspec_third  << "       " << n_total_energy_pspec_third/p_multiplicity_pspec_third  << "        " << n_total_energy_pspec_third << endl;

ofs16.close();

///////////////////////////////////////////////////////////////////////////////////////////////

//Total gamma ray energy

//Find total number of gamma rays emitted
//FABIO; here, if I change to int i=1, then the underflow bin is included...
Double_t sum;
Double_t total_ph_number_0 = 0;
Double_t total_ph_number_1 = 0;
Double_t total_ph_number_2 = 0;


Double_t un_mean_E_0 = 0;
Double_t un_mean_E_1 = 0;
Double_t un_mean_E_2 = 0;
Double_t un_mean = 0;

//Find energy emmited by all gammas per fragment, then divide on gammas per fragment to find mean gamma energy
Double_t mean_ph_E_0 = 0;
Double_t mean_ph_E_1 = 0;
Double_t mean_ph_E_2 = 0;




for (int i=2; i<nbins_h_ph_E_total+2;i++){
  total_ph_number_0 += hframe_ph_E_0->GetBinContent(i);	
  total_ph_number_1 += hframe_ph_E_1->GetBinContent(i);
  total_ph_number_2 += hframe_ph_E_2->GetBinContent(i);

  //The unaccurate way of calculating uncertainty for average gamma energy, uncertainty in number of gammas not included. Works for total gamma energy, where no division on number of gammas
  un_mean_E_0 += hframe_ph_E_0->GetBinContent(i)*hframe_ph_E_0->GetBinCenter(i)*hframe_ph_E_0->GetBinCenter(i);
  un_mean_E_1 += hframe_ph_E_1->GetBinContent(i)*hframe_ph_E_1->GetBinCenter(i)*hframe_ph_E_1->GetBinCenter(i);
  un_mean_E_2 += hframe_ph_E_2->GetBinContent(i)*hframe_ph_E_2->GetBinCenter(i)*hframe_ph_E_2->GetBinCenter(i);
  un_mean += h_ph_E_total->GetBinContent(i)*h_ph_E_total->GetBinCenter(i)*h_ph_E_total->GetBinCenter(i);
  
  mean_ph_E_0 += hframe_ph_E_0->GetBinContent(i)*hframe_ph_E_0->GetBinCenter(i);
  mean_ph_E_1 += hframe_ph_E_1->GetBinContent(i)*hframe_ph_E_1->GetBinCenter(i);
  mean_ph_E_2 += hframe_ph_E_2->GetBinContent(i)*hframe_ph_E_2->GetBinCenter(i);
}

cout << "Total ph numbers: " << total_ph_number_0 << " " << total_ph_number_1 << " " << total_ph_number_2 << endl;

//Average gamma ray energies, per fragment
if(total_ph_number_0==0){
	mean_ph_E_0 = 0;
}
else{
	mean_ph_E_0 = mean_ph_E_0/total_ph_number_0;
}
mean_ph_E_1 = mean_ph_E_1/total_ph_number_1;
mean_ph_E_2 = mean_ph_E_2/total_ph_number_2;


//Total number of photons
Double_t total_ph_number = total_ph_number_0+total_ph_number_1+total_ph_number_2;
//Total gamma ray energy
Double_t total_ph_energy = total_ph_number_0*mean_ph_E_0 + total_ph_number_1*mean_ph_E_1 + total_ph_number_2*mean_ph_E_2;
//Mean gamma ray energy, all fragments
Double_t mean_ph_E = total_ph_energy/total_ph_number;


//More accurate way of calculating uncertainty in average and total gamma energies
Double_t un_mean_ph_E_0_accurate = 0;
Double_t un_mean_ph_E_1_accurate = 0;
Double_t un_mean_ph_E_2_accurate = 0;
Double_t un_mean_ph_accurate = 0;

for(int i=2;i<nbins_h_ph_E_total+2;i++){
	un_mean_ph_E_0_accurate += hframe_ph_E_0->GetBinCenter(i)*hframe_ph_E_0->GetBinCenter(i)*hframe_ph_E_0->GetBinContent(i)/(total_ph_number_0*total_ph_number_0) + hframe_ph_E_0->GetBinCenter(i)*hframe_ph_E_0->GetBinContent(i)/total_ph_number_0;
	un_mean_ph_E_1_accurate += hframe_ph_E_1->GetBinCenter(i)*hframe_ph_E_1->GetBinCenter(i)*hframe_ph_E_1->GetBinContent(i)/(total_ph_number_1*total_ph_number_1) + hframe_ph_E_1->GetBinCenter(i)*hframe_ph_E_1->GetBinContent(i)/total_ph_number_1;
	un_mean_ph_E_2_accurate += hframe_ph_E_2->GetBinCenter(i)*hframe_ph_E_2->GetBinCenter(i)*hframe_ph_E_2->GetBinContent(i)/(total_ph_number_2*total_ph_number_2) + hframe_ph_E_2->GetBinCenter(i)*hframe_ph_E_2->GetBinContent(i)/total_ph_number_2;
	un_mean_ph_accurate += h_ph_E_total->GetBinCenter(i)*h_ph_E_total->GetBinCenter(i)*h_ph_E_total->GetBinContent(i)/(total_ph_number*total_ph_number) + h_ph_E_total->GetBinContent(i)*h_ph_E_total->GetBinCenter(i)/total_ph_number;
}


cout << "AVERAGE GAMMA RAY ENERGIES" << endl;
cout << "Mean photon energy FF0: " << mean_ph_E_0 << " MeV" << " Uncertainty: " << sqrt(un_mean_E_0)/total_ph_number_0<< endl;
cout << "Mean photon energy FF1: " << mean_ph_E_1 << " MeV" << " Uncertainty: "  << sqrt(un_mean_E_1)/total_ph_number_1 << endl;
cout << "Mean photon energy FF2: " << mean_ph_E_2 << " MeV" << " Uncertainty: " << sqrt(un_mean_E_2)/total_ph_number_2 << endl;
cout << "Mean photon energy: " << mean_ph_E  << " MeV" << " Uncertainty: " << sqrt(un_mean)/total_ph_number <<endl;
cout << "\n" << endl;

cout << "TOTAL GAMMA ENERGIES" << endl;
cout << "Total Photon Energy FF0: " << total_ph_number_0*mean_ph_E_0 << " MeV" << " Uncertainty: " << sqrt(un_mean_E_0) << endl;
cout << "Total Photon Energy FF1: " << total_ph_number_1*mean_ph_E_1 << " MeV" << " Uncertainty: " << sqrt(un_mean_E_1) << endl;
cout << "Total Photon Energy FF2: " << total_ph_number_2*mean_ph_E_2 << " MeV" << " Uncertainty: " << sqrt(un_mean_E_2)<< endl;
cout << "Total Photon energy all fragments: " << total_ph_energy<< " MeV" << " Uncertainty: " << sqrt(un_mean) << endl;
cout << "Total Photon energy per fragment: " << total_ph_energy/F << " MeV" << endl;


//Write to file
std::ofstream ofs;
ofs.open ("data_as_func_of_excitation_energy.dat", std::ofstream::out | std::ofstream::app);
//Without scaling of low-energy gammas
//ofs << "         " << p_multiplicity << "            " << sqrt(p_multiplicity*F)/F << "       " << mean_ph_E <<"             "<< sqrt(un_mean)/total_ph_number << "        " << total_ph_energy/F << "             " << sqrt(un_mean) <<  endl;

//With scaling of low-energy gammas
ofs << "         " << p_multiplicity_pspec << "            " << 0 << "       " << p_total_energy_pspec/p_multiplicity_pspec <<"             "<< 0 << "        " << p_total_energy_pspec << "             " << 0;
ofs << "         " << p_multiplicity_pspec_first << "            " << 0 << "       " << p_total_energy_pspec_first/p_multiplicity_pspec_first <<"             "<< 0 << "        " << p_total_energy_pspec_first << "             " << 0;
ofs << "         " << p_multiplicity_pspec_second << "            " << 0 << "       " << p_total_energy_pspec_second/p_multiplicity_pspec_second <<"             "<< 0 << "        " << p_total_energy_pspec_second << "             " << 0;
ofs << "         " << p_multiplicity_pspec_third << "            " << 0 << "       " << p_total_energy_pspec_third/p_multiplicity_pspec_third <<"             "<< 0 << "        " << p_total_energy_pspec_third << "             " << 0 << endl;

ofs.close();

}


void create_frames() {
gStyle->SetOptStat(0);
int nbins;
int maxbin;

hframe_fragyield = new TH1D("hframe_fragyield","",250,0,249);
hframe_fragyield2 = new TH1D("hframe_fragyield2","",250,0,249);
htotal_fragyield = new TH1D("htotal_fragyield","",250,0,249);


nbins = nbins_h_p_mult_total;
maxbin = nbins;
hframe_p_mult_0 = new TH1D("hframe_p_mult_0","",nbins,-0.5,maxbin-.5);
hframe_p_mult_1 = new TH1D("hframe_p_mult_1","",nbins,-0.5,maxbin-.5);
hframe_p_mult_2 = new TH1D("hframe_p_mult_2","",nbins,-0.5,maxbin-0.5);
hframe_p_multi  = new TH2D("hframe_p_multi","",nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5);
hframe_p_multi3D  = new TH3F("hframe_p_multi3D","",nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5,nbins,-0.5,maxbin-0.5);
h_p_mult_total = new TH1D("h_p_mult_total","",nbins,-0.5,maxbin-.5);
hframe_p_mult_first = new TH1D("hframe_p_mult_first","",nbins,-0.5,maxbin-.5);
hframe_p_mult_second = new TH1D("hframe_p_mult_second","",nbins,-0.5,maxbin-.5);
hframe_p_mult_third = new TH1D("hframe_p_mult_third","",nbins,-0.5,maxbin-.5);

hframe_n_mult_first = new TH1D("hframe_n_mult_first","",11,-1.5,10-0.5);
hframe_n_mult_second = new TH1D("hframe_n_mult_second","",11,-1.5,10-0.5);
hframe_n_mult_third = new TH1D("hframe_n_mult_third","",11,-1.5,10-0.5);

nbins = nbins_h_ph_E_total;
maxbin = max_h_ph_E;
hframe_ph_E_0 = new TH1D("hframe_ph_E_0","",nbins,0,maxbin);
hframe_ph_E_1 = new TH1D("hframe_ph_E_1","",nbins,0,maxbin);
hframe_ph_E_2 = new TH1D("hframe_ph_E_2","",nbins,0,maxbin);
h_ph_E_total = new TH1D("h_ph_E_total","",nbins,0,maxbin);

hframe_n_E_1 = new TH1D("hframe_n_E_1","",2*nbins,0,20);
hframe_n_E_2 = new TH1D("hframe_n_E_2","",2*nbins,0,20);
h_n_E_total = new TH1D("h_n_E_total","",2*nbins,0,20);

hframe_ph_E_first = new TH1D("hframe_ph_E_first","",nbins,0,maxbin);
hframe_ph_E_second = new TH1D("hframe_ph_E_second","",nbins,0,maxbin);
hframe_ph_E_third = new TH1D("hframe_ph_E_third","",nbins,0,maxbin);

hframe_n_E_first = new TH1D("hframe_n_E_first","",2*nbins,0,20);
hframe_n_E_second = new TH1D("hframe_n_E_second","",2*nbins,0,20);
hframe_n_E_third = new TH1D("hframe_n_E_third","",2*nbins,0,20);

hframe_fragkin1 = new TH1D("hframe_fragkin1","",1000,0,150);
hframe_fragkin2 = new TH1D("hframe_fragkin2","",1000,0,150);

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