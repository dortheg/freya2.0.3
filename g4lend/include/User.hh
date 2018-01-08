//
// global struct for storing user stuff
//
// must put the following in head of main program
//   #include "User.hh"
//   USER user;
//
// then can use global struct "user" like this
//   #include "User.hh"
//   ...
//   printf("user.blah= %d\n",user.blah);
//
// Doug Wright, LLNL

#include "TFile.h"
#include "TParameter.h"

struct USER
{
    //....variables
    TFile* file;
    double radius; // radius of target sphere (cm)
    
    //....initialize all variables in this struct
     USER(): file(0),radius(0)
            { open("sphere.root"); }
    
    //.....cleanup (and save things) when deleted
    ~USER(){
        file->cd();
        TParameter<int>("radius",radius).Write();
        file->Close();
    }

    //....functions
    TFile* open (TString name){
        file = new TFile(name,"recreate");
        return file;
    }
};

extern USER user; // declares instance "user" as global (must also be created in main)
