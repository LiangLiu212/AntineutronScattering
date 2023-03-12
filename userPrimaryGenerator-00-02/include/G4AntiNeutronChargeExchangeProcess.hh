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
//
//
//
// Geant4 Hadron Elastic Charge Exchange Process -- header file
//
// Created 21 April 2006 V.Ivanchenko
//
// Modified:
//

// Class Description
// Process for hadron nuclear elastic scattering using optimal
// combination of Geant4 models
// Class Description - End

#ifndef G4AntiNeutronChargeExchangeProcess_h
#define G4AntiNeutronChargeExchangeProcess_h 1

#include "globals.hh"
#include "G4HadronicProcess.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Nucleus.hh"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"


class G4ParticleDefinition;

class G4AntiNeutronChargeExchangeProcess : public  G4VDiscreteProcess
{
public:

  G4AntiNeutronChargeExchangeProcess(const G4String& procName = "antineutronchargeExchange");

  virtual ~G4AntiNeutronChargeExchangeProcess();

  virtual G4double GetMeanFreePath(const G4Track & aTrack, G4double previousStepSize, G4ForceCondition* condition);
  virtual G4VParticleChange* PostStepDoIt(const G4Track &, const G4Step&);
	virtual G4bool IsApplicable(const G4ParticleDefinition&);
	virtual void initialize();

private:


  const G4ParticleDefinition* theParticle;

  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theNeutron;
  const G4ParticleDefinition* theAProton;
  const G4ParticleDefinition* theANeutron;
  const G4ParticleDefinition* thePiPlus;
  const G4ParticleDefinition* thePiMinus;
  const G4ParticleDefinition* thePiZero;
  const G4ParticleDefinition* theKPlus;
  const G4ParticleDefinition* theKMinus;
  const G4ParticleDefinition* theK0S;
  const G4ParticleDefinition* theK0L;
  const G4ParticleDefinition* theL;
  const G4ParticleDefinition* theAntiL;
  const G4ParticleDefinition* theSPlus;
  const G4ParticleDefinition* theASPlus;
  const G4ParticleDefinition* theSMinus;
  const G4ParticleDefinition* theASMinus;
  const G4ParticleDefinition* theS0;
  const G4ParticleDefinition* theAS0;
  const G4ParticleDefinition* theXiMinus;
  const G4ParticleDefinition* theXi0;
  const G4ParticleDefinition* theAXiMinus;
  const G4ParticleDefinition* theAXi0;
  const G4ParticleDefinition* theOmega;
  const G4ParticleDefinition* theAOmega;
  const G4ParticleDefinition* theD;
  const G4ParticleDefinition* theT;
  const G4ParticleDefinition* theA;
  const G4ParticleDefinition* theHe3;

  G4ParticleDefinition *theSecondary1; // anti-proton
  G4ParticleDefinition *theSecondary2; // proton1
  G4ParticleDefinition *theSecondary3; // proton2
  G4ParticleDefinition *theSecondary4; // nucleus

	G4DynamicParticle *aDynaSec1;
	G4DynamicParticle *aDynaSec2;
	G4DynamicParticle *aDynaSec3;
	G4DynamicParticle *aDynaSec4;

	G4Track *aSecTrk1;
	G4Track *aSecTrk2;
	G4Track *aSecTrk3;
	G4Track *aSecTrk4;

  G4double        thEnergy;

  G4int    pPDG;
  G4bool   first;


  Float_t pp_rho, pm_rho, angle_pp, angle_pm;
  Float_t miss_mass, miss_rho, angle_miss, ppbar_mass;
  TMVA::Reader *reader0;
  TString method;
  TFile *fweight;
  TH1D *hS;
  TH1D *hB;


};

#endif
