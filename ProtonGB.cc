//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This is the modified version of the example electromagnetic/TestEm7/TestEm7.cc
// This code simulates proton beam with realistic geometry incident on water cube.     

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4ScoringManager.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//#include "CrystalPhysics.hh"
//#include "G4ChannelingPhysics.hh"
 
int main(int argc,char** argv) {

  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
  
  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
    
  //Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  //set mandatory initialization classes
  DetectorConstruction*   det  = new DetectorConstruction();

  //set physics list
  PhysicsList* phys = new PhysicsList();
  //phys->RegisterPhysics(new CrystalPhysics);
  //phys->RegisterPhysics(new G4ChannelingPhysics);
 
  runManager->SetUserInitialization(det);
  runManager->SetUserInitialization(phys);
  
  //set UI-command base scorer
  G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
  scManager->SetVerboseLevel(1);

  //set user action classes
  PrimaryGeneratorAction* kin   = new PrimaryGeneratorAction(det);  
  RunAction*              run   = new RunAction(det,phys,kin);
  EventAction*            event = new EventAction();
  SteppingAction*         step  = new SteppingAction(det,run);              
  
  runManager->SetUserAction(kin); 
  runManager->SetUserAction(run); 
  runManager->SetUserAction(event); 
  runManager->SetUserAction(step);

  // Initialize visualization
  //

   G4VisManager* visManager = new G4VisExecutive;
   visManager->Initialize();

  //get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

 // Process macro or start UI session
 //
  
  if (! ui)   // batch mode  
    {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
    }
    
  else           //UI terminal for interactive mode
    { 
    UI->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
    }    

/*     
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);      
      ui->SessionStart();
      delete ui;
    }
*/
 
  delete visManager;   
    
  //job termination
  delete runManager;

  //return 0;
}
 
