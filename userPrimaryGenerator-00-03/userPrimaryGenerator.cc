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
/// \file eventgenerator/userPrimaryGenerator/userPrimaryGenerator.cc
/// \brief Main program of the eventgenerator/userPrimaryGenerator example
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "G4PhysListFactory.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingVerbose.hh"
#include "PhysicsListMessenger.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"


#include "BesRunAction.hh"
#include "BesEventAction.hh"
#include "BesSteppingAction.hh"
#include "TROOT.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc,char** argv) {

gROOT->GetInterpreter();
  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);

  //construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  //set mandatory initialization classes
  //
  runManager->SetUserInitialization(new DetectorConstruction);

  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = nullptr;
  PhysicsListMessenger* mess = nullptr;
  G4String physName = "";
  G4String outfileName = "";

  //Physics List name defined via 3nd argument
  if (argc == 4) { outfileName = argv[2]; }
  if (argc == 3) { outfileName = argv[2]; }

  // Physics List name defined via environment variable
  if("" == physName) {
    char* path = getenv("PHYSLIST");
    if (path) { physName = G4String(path); }
  }

//  physName= "FTFP_BERT";
  //reference PhysicsList via its name
  if ("" != physName && factory.IsReferencePhysList(physName)) {
    phys = factory.GetReferencePhysList(physName);

    // instantiated messenger
    mess = new PhysicsListMessenger();
  } 

  //local Physics List
  if(!phys) { phys = new PhysicsList();
  
  }
	PhysicsList *phylist = new PhysicsList();
	phylist->AddPhysicsList("FTFP_BERT");

 


  runManager->SetUserInitialization(phylist);
	PrimaryGeneratorAction *PGA = new PrimaryGeneratorAction();
	PGA->SetRootFile(argv[3]);
  runManager->SetUserAction(PGA);
  BesRunAction *runAction = new BesRunAction;
  runManager->SetUserAction(runAction);
  BesEventAction *eventAction = new BesEventAction(runAction);
  runManager->SetUserAction(eventAction);
  BesSteppingAction *stepAction = new BesSteppingAction(eventAction);
  stepAction->SetOutFile(outfileName);
  runManager->SetUserAction(stepAction);

  //initialize visualization
  G4VisManager* visManager = nullptr;

  //get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (ui)  {
   //interactive mode
   visManager = new G4VisExecutive;
   visManager->Initialize();
   UImanager->ApplyCommand("/control/execute vis.mac");          
   ui->SessionStart();
   delete ui;
  }
  else  {
   //batch mode  
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   UImanager->ApplyCommand(command+fileName);
  }

  //job termination 
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 

