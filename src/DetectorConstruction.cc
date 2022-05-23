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
// This is the modified version of the example electromagnetic/TestEm7/src/DetectorConstruction.cc
   
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
//#include "SensitiveDetector.hh"
#include "G4SDManager.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"

#include "G4TransportationManager.hh"
#include "G4RunManager.hh" 

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4CrystalExtension.hh"
#include "G4ExtendedMaterial.hh"
#include "G4LogicalCrystalVolume.hh"
#include "MaterialExtensionData.hh"
//#include "G4ChannelingMaterialData.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  //fWorldMaterial(0),fAbsorMaterial(0),fLAbsor(0),
  fDetectorMessenger(0)
{

// Default Parameters

  // Aluminium Oxide Crystal
  cryst_dX = 1.0*mm;
  cryst_dY = 1.0*mm;
  cryst_dZ = 1.0*mm;
  
  // World
  world_sizeXY = 10.0*cryst_dX;
  world_sizeZ  = 10.0*cryst_dZ;

  fLayerSizeYZ = 0.1*mm;
  fLayerNumber = 1;
  fLayerSizeX = 0*mm;
  fLayerMass = 0;
  
 // fWorldMaterial = fAbsorMaterial = 0;
  //fLAbsor   = 0;
      
  DefineMaterials();
  
  // create commands for interactive definition of the detector  
  fDetectorMessenger = new DetectorMessenger(this);
}


DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}


void DetectorConstruction::DefineMaterials()
{ 

/*

  G4NistManager* man = G4NistManager::Instance();

  G4bool isotopes = false;

  G4Element*  O = man->FindOrBuildElement("O" , isotopes);
  G4Element* Al = man->FindOrBuildElement("Al", isotopes);

  G4Material* Al2O3 = new G4Material("Al2O3", 3.97*g/cm3, 2); // density from PubChem
  Al2O3->AddElement(Al, 2);
  Al2O3->AddElement(O , 3);
  
   G4bool isotopes = false;
  
   G4Material* Al2O3 = man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE Al_2O_3", isotopes); 
   
   */
   
}

  
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // Aluminium Oxide Crystal
  G4NistManager* nist = G4NistManager::Instance();
  cryst_mat = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
   
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // World                     
  //
  world_mat = nist->FindOrBuildMaterial("G4_AIR");
                             
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking   
                      
  // define crystal
  //
  
  G4Box* solidCryst = new G4Box("crystal", cryst_dX/2, cryst_dZ/2, cryst_dZ/2);
  
  /*

  G4LogicalVolume* logicCryst =
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicCryst,              //its logical volume
                    "CrystalLV",             //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  
    
    */
    
    G4ExtendedMaterial* CrystalMat = new G4ExtendedMaterial("crystal.material",cryst_mat);
    
    CrystalMat->RegisterExtension(std::unique_ptr<G4CrystalExtension>(new G4CrystalExtension(CrystalMat)));
    G4CrystalExtension* crystalExtension = (G4CrystalExtension*)CrystalMat->RetrieveExtension("crystal");
    
    crystalExtension->SetUnitCell(new G4CrystalUnitCell(4.75 * CLHEP::angstrom,
                                                        4.75 * CLHEP::angstrom,
                                                        12.99 * CLHEP::angstrom,
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi,
                                                        2/3 * CLHEP::pi,
                                                        167));
    
    CrystalMat->RegisterExtension(std::unique_ptr<MaterialExtensionData>(new MaterialExtensionData("ExtendedData")));
    MaterialExtensionData* materialExtension = (MaterialExtensionData*)CrystalMat->RetrieveExtension("ExtendedData");
    materialExtension->SetValue(57.);                                                     
                                                        
         /*

    CrystalMat->RegisterExtension(std::unique_ptr<G4ChannelingMaterialData>(new G4ChannelingMaterialData("channeling")));
    G4ChannelingMaterialData* crystalChannelingData = (G4ChannelingMaterialData*)CrystalMat->RetrieveExtension("channeling"); 
    crystalChannelingData->SetFilename(fECfileName);
    
    */  
    
    G4LogicalCrystalVolume* fBoxLogicCrystal = new G4LogicalCrystalVolume(solidCryst,
                                                                          CrystalMat,
                                                                          "crystal.logic");
    
    fBoxLogicCrystal->SetVerbose(1);
    
    new G4PVPlacement(0,
                      G4ThreeVector(),
                      fBoxLogicCrystal,
                      "crystal.physical",
                      logicWorld,
                      false,
                      0);  
                                                                                     

   // Layer (i.e., slices)
   if (fLayerNumber > 0) {
      fLayerSizeX = cryst_dX/fLayerNumber;
      G4Box* sLayer = new G4Box("Layer", fLayerSizeX/2,fLayerSizeYZ/2,fLayerSizeYZ/2); 
      G4LogicalVolume* fLLayer = new G4LogicalVolume(sLayer, cryst_mat, "Layer");	
      G4VPhysicalVolume* fPhysLayer = new G4PVReplica("Layer",		
      		                    fLLayer,		
      	                            fBoxLogicCrystal,		
                                    kXAxis,		
                                    fLayerNumber,		
                                    fLayerSizeX);	

      fLayerMass = fLayerSizeX*fLayerSizeYZ*fLayerSizeYZ*(cryst_mat->GetDensity());

    }                

  PrintParameters();
    
  //always return the World volume 
  return physWorld;
}

void DetectorConstruction::PrintParameters()
{
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Crystal is " << G4BestUnit(cryst_dX,"Length")
         << " of " << cryst_mat->GetName() << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
  
}

void DetectorConstruction::SetSizeX(G4double value)
{
  cryst_dX = value; world_sizeXY = 10.0*cryst_dX;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
  
void DetectorConstruction::SetSizeYZ(G4double value)
{
  cryst_dY = value; 
  world_sizeXY = 10.0*cryst_dY;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}  

void DetectorConstruction::SetLayerSizeYZ(G4double value)
{
  fLayerSizeYZ = value; 
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}  

void DetectorConstruction::SetLayerNumber(G4int value)
{
  fLayerNumber = value; 
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}  
  
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

/*
void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  
  // declare crystal as a sensitive scorer
  //  
  G4SDManager* SDmanager = G4SDManager::GetSDMpointer();
  SensitiveDetector* crystSD = new SensitiveDetector("/MycrystalSD");
  SDmanager->AddNewDetector(crystSD);
  crystal->SetSensitiveDetector(crystSD);
}

*/

