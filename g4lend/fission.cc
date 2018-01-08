//******************************************************************************
// fission.cc  GEANT4 user application for testing handling of
//             low-energy neutron induced fission
//
// fission                 # start interactive version
// fission test.mac        # read test.mac for geant commands
// fission -i test.mac     # read test.mac for geant commans and continue interactively
//
// 1.00 JMV, LLNL, APR-2017:  First version.
//******************************************************************************
//

// geant4 includes
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4ios.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

// misc includes
#include <fstream>
#include <math.h>

// package includes
#include "DetectorConstruction.hh"
#include "LENDPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "User.hh"
#include "options.hh"
USER user; // create global struct, see USER.hh
#include "End.hh"

//------------------------------------------------------------------------------
int main(int argc,char** argv) {
  string session_type = "";  //....use default

  //....parse command line options
  bool interactive          = options("i",argc,argv);
  bool interactive_terminal = options("t",argc,argv);
  bool vis                  = options("v",argc,argv);
  if (interactive_terminal){ interactive = true; session_type = "tcsh";}
  vector<string> files = options(argc,argv);

  //....construct default run manager
  G4RunManager* runManager = new G4RunManager;

  //....set mandatory initialization and user classes
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new LENDPhysicsList);
  runManager->SetUserAction(new PrimaryGeneratorAction());

  //....initialize G4 kernel
  runManager->Initialize();
  
  //....initialize visualization manager (if requested)
  G4VisManager* visManager=NULL;
  if (vis) {
    visManager = new G4VisExecutive;
    visManager->Initialize();
  }
   
  //....get user interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  //....get interactive user interface, will automatically try (in order): Qt, tcsh, Xm
  G4UIExecutive* UI = NULL;
  if (interactive) UI = new G4UIExecutive(argc, argv, session_type);

  //....startup macro files (optional -  will continue if file not found)
  //    must come after G4UIExecutive() but before SessionStart()
  if (vis) UImanager->ApplyCommand("/control/execute vis.mac");  // graphics definitions
  //....execute default macro file init.mac, unless other files given on the command line
  if( files.size() == 0 )
    UImanager->ApplyCommand("/control/execute init.mac");
  else
    for( unsigned int i=0; i<files.size(); i++)
      UImanager->ApplyCommand("/control/execute "+files[i]);

  //....start interactive user interface (if requested)
  if (interactive){ UI->SessionStart(); delete UI;}

  //....my personal end of job function
  End();

  //....cleanup required, otherwise you get warnings
  if (NULL!=visManager) delete visManager;
  delete runManager;
  return EXIT_SUCCESS;
}
