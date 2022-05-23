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
// This is the modified version of the example electromagnetic/TestEm7/src/SteppingAction.cc

#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "Randomize.hh"

#include "Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
//#include "G4RunManager.hh"
#include "G4Event.hh"

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* RuAct)
:G4UserSteppingAction(),fDetector(det), fRunAction(RuAct)
{ }

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
   
   G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition();
   
   const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();		
   fEventID        = evt->GetEventID(); 
   fPType          = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
   fPosX           = pos.x();
   fPosY           = pos.y();
   fPosZ           = pos.z();  
   fEnergyDeposit  = aStep->GetTotalEnergyDeposit();  
   fKineticEnergy  = aStep->GetPreStepPoint()->GetKineticEnergy();
   fStepLength     = aStep->GetStepLength();
   fparentId       = aStep->GetTrack()->GetParentID();
   
   if (aStep->GetPostStepPoint()->GetStepStatus() == fWorldBoundary)
   	fLastVolumeName = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
   else
   	fLastVolumeName = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
   
   const G4VProcess* originProcess = aStep->GetTrack()->GetCreatorProcess();
   if (originProcess)
	fOriginProcessName = originProcess->GetProcessName();
   else
	fOriginProcessName = "Primary";

   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   analysisManager->OpenFile("CrystalTrackInfo");		
   analysisManager->FillNtupleIColumn(0,0, fEventID);
   analysisManager->FillNtupleDColumn(0,1, fPosX);
   analysisManager->FillNtupleDColumn(0,2, fPosY);
   analysisManager->FillNtupleDColumn(0,3, fPosZ);
   analysisManager->FillNtupleDColumn(0,4, fEnergyDeposit * 1000);
   analysisManager->FillNtupleDColumn(0,5, fKineticEnergy);   
   analysisManager->FillNtupleDColumn(0,6, fStepLength * 1000);
   analysisManager->FillNtupleIColumn(0,7, fPType);
   analysisManager->FillNtupleSColumn(0,8, fOriginProcessName);
   analysisManager->FillNtupleSColumn(0,9, fLastVolumeName);
   analysisManager->FillNtupleDColumn(0,10, fparentId);
   analysisManager->AddNtupleRow(0);  
   
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep <= 0.) return;
  
  fRunAction->FillEdep(edep);
  
  if (aStep->GetTrack()->GetTrackID() == 1) {
    fRunAction->AddPrimaryStep();
    /*
    G4cout << aStep->GetTrack()->GetMaterial()->GetName()
           << "  E1= " << aStep->GetPreStepPoint()->GetKineticEnergy()
           << "  E2= " << aStep->GetPostStepPoint()->GetKineticEnergy()
           << " Edep= " << edep << G4endl;
    */
  } 

  //Bragg curve
  G4StepPoint* prePoint  = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
   
  G4double x1 = prePoint->GetPosition().x();
  G4double x2 = postPoint->GetPosition().x();  
  G4double x  = x1 + G4UniformRand()*(x2-x1) + 0.5*(fDetector->GetAbsorSizeX());
  G4AnalysisManager* analman = G4AnalysisManager::Instance();
  analman->FillH1(1, x, edep); 
  analman->FillH1(2, x, edep);
  
  //fill layers
  G4int copyNb = prePoint->GetTouchableHandle()->GetCopyNumber();
  if (copyNb > 0) fRunAction->FillLayerEdep(copyNb, edep); 
}




