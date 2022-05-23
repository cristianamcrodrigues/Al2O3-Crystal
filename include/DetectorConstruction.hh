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
// This is the modified version of the example electromagnetic/TestEm7/include/DetectorConstruction.hh

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

const G4int MaxLayer = 60;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
       
     void SetSizeX         (G4double);
     void SetSizeYZ        (G4double);   
     void SetLayerNumber   (G4int);
     void SetLayerSizeYZ   (G4double);                          
     
     virtual     
     G4VPhysicalVolume* Construct();
     void UpdateGeometry();
     //void ConstructSDandField() override;
     
  public:  
  
     G4double     GetWorldSizeX()    {return world_sizeXY;};
     G4double     GetWorldSizeYZ()   {return world_sizeZ;};
     G4Material*  GetWorldMaterial() {return world_mat;};     
     G4double     GetAbsorSizeX()    {return cryst_dX;};
     G4double     GetAbsorSizeYZ()   {return cryst_dZ;};           
     G4Material*  GetAbsorMaterial() {return cryst_mat;};
      
     G4int            GetLayerNumber()         {return fLayerNumber;};
     G4double         GetLayerMass()           {return fLayerMass;}; 
     G4double         GetLayerSizeYZ()         {return fLayerSizeYZ;};

     void         PrintParameters();
                       
  private:
  
     G4double            fLayerSizeYZ;
     G4double            fLayerSizeX;
     
     G4int               fLayerNumber;
     G4double            fLayerMass; 
     DetectorMessenger*  fDetectorMessenger;
     
     G4double		  cryst_dX;
     G4double		  cryst_dY;
     G4double		  cryst_dZ;
     G4double		  world_sizeXY;     
     G4double		  world_sizeZ;
     
     G4Material*         world_mat; 
     G4Material*	  cryst_mat; 
     
     G4String fECfileName;
     
     //void ConstructSDandField();         


  private:
    
     void DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

#endif

