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
// This is the modified version of the example electromagnetic/TestEm7/src/RunAction.cc
  
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "StepMax.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "Randomize.hh"
#include "Analysis.hh"

RunAction::RunAction(DetectorConstruction* det, PhysicsList* phys,
                     PrimaryGeneratorAction* kin)
 : G4UserRunAction(),
   fAnalysisManager(0), fDetector(det), fPhysics(phys), fKinematic(kin),
   fLayerEdep(new G4double[MaxLayer]), fEdeptot(0.), fNbPrimarySteps(0)
{ 
  // Book predefined histograms
  BookHisto();
}

RunAction::~RunAction()
{
  delete [] fLayerEdep;
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  if(!fAnalysisManager) { BookHisto(); }
  
  CLHEP::HepRandom::showEngineStatus();
     
  //initialize projected range, Ebeam, and book histograms
  fNbPrimarySteps = 0;
  fEdeptot = 0.;
  for (G4int j=0; j<MaxLayer; j++) fLayerEdep[j] = 0.;
  fKinematic->ResetEbeamCumul();
   
  if (fAnalysisManager->IsActive()) {
    fAnalysisManager->OpenFile(); 
    // histogram "1" is defined by the length of the target
    G4double length  = fDetector->GetAbsorSizeX();
    G4double stepMax = fPhysics->GetStepMaxProcess()->GetMaxStep();
    G4int nbmin = 100;
    G4int nbBins = (G4int)(0.5 + length/stepMax);
    if (nbBins < nbmin) nbBins = nbmin;
    fAnalysisManager->SetH1(1, nbBins, 0., length, "mm");
  }
  
   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   analysisManager->SetVerboseLevel(1);
   
   analysisManager->OpenFile("CrystalTrackInfo");	
   analysisManager->CreateNtuple("01","Track Information Al2O3 Crystal");
   analysisManager->CreateNtupleIColumn("Event ID");
   analysisManager->CreateNtupleDColumn("Position X (mm)");
   analysisManager->CreateNtupleDColumn("Position Y (mm)");
   analysisManager->CreateNtupleDColumn("Position Z (mm)");
   analysisManager->CreateNtupleDColumn("Energy Deposit (keV)");
   analysisManager->CreateNtupleDColumn("Initial Kinetic Energy (MeV)");   
   analysisManager->CreateNtupleDColumn("Step Length (um)");
   analysisManager->CreateNtupleIColumn("Particle Type (in PDG Format)");
   analysisManager->CreateNtupleSColumn("Origin Process");
   analysisManager->CreateNtupleSColumn("Last Volume");
   analysisManager->CreateNtupleDColumn("Parent ID");   
   analysisManager->FinishNtuple();	
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nbofEvents = aRun->GetNumberOfEvent();
  if (nbofEvents == 0) return;
    
  //run conditions 
  G4Material* material = fDetector->GetAbsorMaterial();
         
  //compute number of primary steps
  G4double nstep = G4double(fNbPrimarySteps)/G4double(nbofEvents);
  G4cout << " Mean number of primary steps = "<< nstep << G4endl;

  if (fAnalysisManager->IsActive() ) {        
    // normalize histograms  
    for (G4int j=1; j<3; j++) {  
      G4double binWidth = fAnalysisManager->GetH1Width(j);
      G4double fac = (mm/MeV)/(nbofEvents * binWidth);
      fAnalysisManager->ScaleH1(j, fac);
    }
    
    // save histograms
    fAnalysisManager->Write();
    fAnalysisManager->CloseFile();  
    delete fAnalysisManager;
    fAnalysisManager = 0;
  }
   
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

void RunAction::BookHisto()
{
  // Create or get analysis manager
  fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->SetFileName("ProtonGB");
  fAnalysisManager->SetVerboseLevel(1);
  fAnalysisManager->SetActivation(true);  // enable inactivation of histograms

  // Define histograms start values
  const G4int kMaxHisto = 3;
  const G4String id[] = { "0", "1", "2"};
  const G4String title[] = 
                { "dummy",                                      //0
                  "Edep (MeV/mm) along absorber ",              //1
                  "Edep (MeV/mm) along absorber zoomed"         //2
                 }; 

  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = fAnalysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    G4bool activ = false;
    if (k == 1) activ = true;
    fAnalysisManager->SetH1Activation(ih, activ);
  }
}


