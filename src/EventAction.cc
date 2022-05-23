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
// This is the modified version of the example electromagnetic/TestEm7/src/EventAction.cc
   
#include "EventAction.hh"
#include "EventActionMessenger.hh"
#include "G4Event.hh"

EventAction::EventAction()
:G4UserEventAction(),fDrawFlag("none"),fPrintModulo(10000),fEventMessenger(0)
{
  fEventMessenger = new EventActionMessenger(this);
}

EventAction::~EventAction()
{
  delete fEventMessenger;
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();

 //printing survey
 if (evtNb%fPrintModulo == 0)
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
}

void EventAction::EndOfEventAction(const G4Event*)
{

/*
 G4HCofThisEvent* HCE = event->GetHCofThisEvent();
 SDHitCollection* Hits = 0;
 Hits = GetHitCollection(HCE,"SDHitCollection"); // Particles detected at detector
 
 if (Hits) {
   const G4int nHits = Hits->entries();
   for (G4int iHit = 0; iHit<nHits; ++iHit) {
     outFile << (*Hits)[iHit] << G4endl;
   }
  }
 else
   G4Exception("EndOfEventAction","",JustWarning,"Hits collection SDHitCollection not
found."); 

*/
}

/*
SDHitCollection* EventAction::GetHitCollection(G4HCofThisEvent* HCE, const G4String & name)
{
 static G4SDManager* SDman = G4SDManager::GetSDMpointer();
 const G4int HitCollID = SDman->GetCollectionID(name);
 if (HCE && HitCollID > -1)
   return (SDHitCollection*)(HCE->GetHC(HitCollID));
 else
   return 0;
}
*/


