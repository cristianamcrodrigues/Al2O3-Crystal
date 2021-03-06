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
// This is the modified version of the example electromagnetic/TestEm7/include/SteppingAction.hh

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class RunAction;
class G4Event;

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(DetectorConstruction*, RunAction*);
    virtual ~SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
    
  private:
    DetectorConstruction* fDetector;
    RunAction*            fRunAction;
    //G4Event* aEvent;
    
	// Output variables
	G4int fEventID;
	G4float fPosX;
	G4float fPosY;
	G4float fPosZ;
	G4float fEnergyDeposit;
	G4float fKineticEnergy;
	G4float fStepLength;
	G4int fPType;
	G4String fOriginProcessName;
	G4String fLastVolumeName;
	G4double fparentId;
};

#endif
