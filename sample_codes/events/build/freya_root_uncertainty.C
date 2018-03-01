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
#include <iostream>

//////////////////////////////////////////////////////////////
// To be used in Dorthea's master, to calculate uncertainty in FREYA simulations
//
// Includes:
// - Cutting a long .dat.root file in smaller trees, alalyzing each tree induvidually, and writing average
// gamma energy, avg gamma multiplicity and total gamma energy to a file ...., where a script is used to find
// the uncertainties by formula
// -


// OBS: the error calculations for average and total gamma energies is not _entirely_ correct. 
// As the total number of gammas and the gamma multiplicities for certain gamma energies is correlated,
// the error term should include both error in total number of gammas AND the covariance
// See Master Journal pages 16-18 for more info


//OBS: Root places energies in bins. Therefore, the exact energies are lost, how much depends on the size of the bin
/////////////////////////////////////////////////////////////////

 //void create_frames();


void freya_root_uncertainty(){

  TFile *vetsex = new TFile("Cf252.dat.root", "READ");
  TTree *mytree_all = (TTree *) gROOT->FindObject("FreyaTree");

  for(int i=0;i<3;i++){
    
    int k1 = i * 1000;
    int k2 = (i+1) * 1000;
    string nEntry1 = to_string(int(k1));
    string nEntry2 = to_string(int(k2));
    auto my_cutStr = "Entry$>"  + nEntry1 + "&&" +  "Entry$<"  + nEntry2;
    cout << "Cut: " << my_cutStr << endl;
    const char *my_cut = my_cutStr.c_str();

    TTree* mytree = mytree_all->CopyTree(my_cut);


  //TTree* mytree = mytree_all->CopyTree("Entry$<1000");
  

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
    int nbins_h_p_mult_total  = 20;

    TH1D *hframe_ph_E_0;
    TH1D *hframe_ph_E_1;
    TH1D *hframe_ph_E_2;
    TH1D *h_ph_E_total;
    // TH1D *h_n_E_Boltzmann;
    int nbins_h_ph_E_total  = 500;
    int max_h_ph_E  = 7; // in MeV

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

    //create frames

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

    nbins = nbins_h_ph_E_total;
    maxbin = max_h_ph_E;
    hframe_ph_E_0 = new TH1D("hframe_ph_E_0","",nbins,0,maxbin);
    hframe_ph_E_1 = new TH1D("hframe_ph_E_1","",nbins,0,maxbin);
    hframe_ph_E_2 = new TH1D("hframe_ph_E_2","",nbins,0,maxbin);
    h_ph_E_total = new TH1D("h_ph_E_total","",nbins,0,maxbin);

    TCanvas *c1 = new TCanvas("c1","Fragment Yield",150,10,990,660);
    mytree->Draw("iAf1>>hframe_fragyield");
    mytree->Draw("iAf2>>hframe_fragyield2");


    int sum_F = 0;

    for(int i=0;i<300;i++){
      sum_F += hframe_fragyield->GetBinContent(i);
    }

    cout << "Sum: " << sum_F << endl;

    int F = sum_F; //Number of fissions -> Implement this automatically?

    Double_t mean;
    Double_t norm = 2;
    Double_t nu;
    Double_t value;

    //////////////////////////////////////////////////////
    //Photon multiplicities
    //////////////////////////////////////////////////////

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
    cout << "Mean photon energy FF0: " << mean_ph_E_0 << " MeV" << endl;
    cout << "Mean photon energy FF1: " << mean_ph_E_1 << " MeV" << endl;
    cout << "Mean photon energy FF2: " << mean_ph_E_2 << " MeV" << endl;;
    cout << "Mean photon energy: " << mean_ph_E  << " MeV" << endl;
    cout << "\n" << endl;

    cout << "TOTAL GAMMA ENERGIES" << endl;
    cout << "Total Photon Energy FF0: " << total_ph_number_0*mean_ph_E_0 << " MeV" << endl;
    cout << "Total Photon Energy FF1: " << total_ph_number_1*mean_ph_E_1 << " MeV" << endl;
    cout << "Total Photon Energy FF2: " << total_ph_number_2*mean_ph_E_2 << " MeV" << endl;
    cout << "Total Photon energy all fragments: " << total_ph_energy<< " MeV" << endl;

  //Delete previous canvas'es
  delete c1;
  delete c4;
  delete pMults;
  delete c6;
  delete hframe_fragyield;
  delete hframe_fragyield2;
  delete htotal_fragyield;
  delete hframe_p_mult_0;
  delete hframe_p_mult_1;
  delete hframe_p_mult_2;
  delete hframe_p_multi;
  delete hframe_p_multi3D;
  delete h_p_mult_total;
  delete hframe_ph_E_0;
  delete hframe_ph_E_1;
  delete hframe_ph_E_2;
  delete h_ph_E_total;


    //Write to file
  std::ofstream ofs;
  ofs.open ("uncertainties.dat", std::ofstream::out | std::ofstream::app);
  ofs << p_multiplicity << " " << mean_ph_E <<" " << total_ph_energy << endl;
  ofs.close();
  }
}
