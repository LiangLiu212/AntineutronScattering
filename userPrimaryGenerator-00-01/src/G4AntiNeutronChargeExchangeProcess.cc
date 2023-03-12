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
// Geant4 Hadron Charge Exchange Process -- source file
//
// Created 21 April 2006 V.Ivanchenko
//
// Modified:
// 24-Apr-06 V.Ivanchenko add neutron scattering on hydrogen from CHIPS
// 07-Jun-06 V.Ivanchenko fix problem of rotation of final state
// 25-Jul-06 V.Ivanchenko add 19 MeV low energy for CHIPS
// 23-Jan-07 V.Ivanchenko add cross section interfaces with Z and A
//                        and do not use CHIPS for cross sections
// 14-Sep-12 M.Kelsey -- Pass subType code to base ctor
// 06-Aug-15 A.Ribon migrating to G4Pow

#include "G4AntiNeutronChargeExchangeProcess.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4IsotopeVector.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4IonTable.hh"
#include "TMath.h"
#include "G4AtomicShells.hh"
#include "G4Element.hh"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "G4ParticleTable.hh"
#include "G4Pow.hh"
using namespace std;

G4AntiNeutronChargeExchangeProcess::G4AntiNeutronChargeExchangeProcess(const G4String& procName)
  : G4VDiscreteProcess(procName), first(true)
{
  thEnergy = 20.*MeV;
  pPDG = 0;
  verboseLevel= 1;
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theAProton  = G4AntiProton::AntiProton();
  theANeutron = G4AntiNeutron::AntiNeutron();
  thePiPlus   = G4PionPlus::PionPlus();
  thePiMinus  = G4PionMinus::PionMinus();
  thePiZero   = G4PionZero::PionZero();
  theKPlus    = G4KaonPlus::KaonPlus();
  theKMinus   = G4KaonMinus::KaonMinus();
  theK0S      = G4KaonZeroShort::KaonZeroShort();
  theK0L      = G4KaonZeroLong::KaonZeroLong();
  theL        = G4Lambda::Lambda();
  theAntiL    = G4AntiLambda::AntiLambda();
  theSPlus    = G4SigmaPlus::SigmaPlus();
  theASPlus   = G4AntiSigmaPlus::AntiSigmaPlus();
  theSMinus   = G4SigmaMinus::SigmaMinus();
  theASMinus  = G4AntiSigmaMinus::AntiSigmaMinus();
  theS0       = G4SigmaZero::SigmaZero();
  theAS0      = G4AntiSigmaZero::AntiSigmaZero();
  theXiMinus  = G4XiMinus::XiMinus();
  theXi0      = G4XiZero::XiZero();
  theAXiMinus = G4AntiXiMinus::AntiXiMinus();
  theAXi0     = G4AntiXiZero::AntiXiZero();
  theOmega    = G4OmegaMinus::OmegaMinus();
  theAOmega   = G4AntiOmegaMinus::AntiOmegaMinus();
  theD        = G4Deuteron::Deuteron();
  theT        = G4Triton::Triton();
  theA        = G4Alpha::Alpha();
  theHe3      = G4He3::He3();
}

G4AntiNeutronChargeExchangeProcess::~G4AntiNeutronChargeExchangeProcess()
{

	if(aSecTrk1) delete aSecTrk1;
	if(aSecTrk2) delete aSecTrk2;
	if(aSecTrk3) delete aSecTrk3;
	if(aSecTrk4) delete aSecTrk4;
	if(aDynaSec1) delete aDynaSec1;
	if(aDynaSec2) delete aDynaSec2;
	if(aDynaSec3) delete aDynaSec3;
	if(aDynaSec4) delete aDynaSec4;
	if(theSecondary1) delete theSecondary1;
	if(theSecondary2) delete theSecondary2;
	if(theSecondary3) delete theSecondary3;
	if(theSecondary4) delete theSecondary4;

}

G4double G4AntiNeutronChargeExchangeProcess::GetMeanFreePath(const G4Track &aTrack, G4double previousStepSize, G4ForceCondition* condition){

		int debug = 0;
		
		const G4Material* material = aTrack.GetMaterial();
		G4double density, A;
		G4double unitCrossSection = 200*microbarn;
	//	G4double myNa = 6.0221415;
		if(material->GetName() == "Beam") return 1e100;
		if(material->GetName() == "G4_AIR") return 1e100;
		if(material->GetName() == "Oil"){
				if(debug == 1){
						std::cout << "From G4AntiNeutronChargeExchangeProcess::GetMeanFreePath" << std::endl;
						std::cout << "G4Material : name : " << material->GetName() << std::endl;
						std::cout << "G4Material : Z    : 1 6" << std::endl;
						std::cout << "G4Material : A    : " << 14.13 << std::endl;
						std::cout << "G4Material : Density : " << material->GetDensity()/(g/cm3) << std::endl;
				}
				density = material->GetDensity()/(g/cm3);
				A = 14.13;
		}
		else if(material->GetName() == "Beryllium"){
				if(debug == 1){
						std::cout << "G4Material : name : " << material->GetName() << std::endl;
						std::cout << "G4Material : Z    : " << material->GetZ() << std::endl;
						std::cout << "G4Material : A    : " << material->GetA()/(g/mole) << std::endl;
						std::cout << "G4Material : Density : " << material->GetDensity()/(g/cm3) << std::endl;
				}
				density = material->GetDensity()/(g/cm3);
				A = material->GetA()/(g/mole);
		}
		else if(material->GetName() == "Gold"){
				if(debug == 1){
						std::cout << "G4Material : name : " << material->GetName() << std::endl;
						std::cout << "G4Material : Z    : " << material->GetZ() << std::endl;
						std::cout << "G4Material : A    : " << material->GetA()/(g/mole) << std::endl;
						std::cout << "G4Material : Density : " << material->GetDensity()/(g/cm3) << std::endl;
				}
				density = material->GetDensity()/(g/cm3);
				A = material->GetA()/(g/mole);

		}
		else 
				return 1e100;
		G4double trkmom = aTrack.GetMomentum().mag()/MeV;
		G4double lambdaI = A / ( density *  ( unitCrossSection / barn ) * TMath::Power(A, 2/3));
		if(debug == 1){
				std::cout << "Decay Length : " << lambdaI*cm << std::endl;
				std::cout << "Decay Length : " << lambdaI << std::endl;
				std::cout << " Momentum : " << trkmom << std::endl;
				std::cout << " Cross Section : " << unitCrossSection / barn << std::endl;
		}
		return lambdaI * 0.001;
}

G4VParticleChange* G4AntiNeutronChargeExchangeProcess::PostStepDoIt(const G4Track &aTrack,  const G4Step &aStep){


		int debug = 0;
		Double_t m_Proton = 0.938272081;
		Double_t m_Neutron = 0.939565413;
		Double_t m_Proton_amu = 1.00727646688;
		Double_t amu = m_Proton/m_Proton_amu;

		const G4Material* material = aTrack.GetMaterial();
		std::cout << material->GetName() << std::endl;
		G4double density, A, Z;
		TLorentzVector p4NucleusMother;
		TLorentzVector p4NucleusDaughter;
		G4double binding_energy, m_mass_daughter;
		if(material->GetName() == "Oil"){  // Carb
				G4double m_mass = 12 * amu;    // GeV
				p4NucleusMother.SetPxPyPzE(0,0,0, m_mass);
				binding_energy = 7680.144 * keV;
				m_mass_daughter = 11.011433563 * amu;
				A = 11;
				Z = 6;
		}
		else if(material->GetName() == "Beryllium"){
				G4double m_mass = 9.012183065 * amu;    // GeV
				p4NucleusMother.SetPxPyPzE(0,0,0, m_mass);
				binding_energy = 6462.668 * keV;
				m_mass_daughter = 8.005305102 * amu;
				A = 8;
				Z = 4;






		}
		else if(material->GetName() == "Gold"){
				G4double m_mass = 196.966568786 * amu;    // GeV
				p4NucleusMother.SetPxPyPzE(0,0,0, m_mass);
				binding_energy = 7680.144 * keV;
				m_mass_daughter = 195.966569908 * amu;
				A = 196;
				Z = 79;
		}
		else{
				aParticleChange.Initialize(aTrack);
				aParticleChange.ProposeEnergy(0.);
				aParticleChange.ProposeTrackStatus(fStopAndKill);
				return &aParticleChange;
		}
		if(aTrack.GetMomentum().mag() <  binding_energy){
				aParticleChange.Initialize(aTrack);
				aParticleChange.ProposeEnergy(0.);
				aParticleChange.ProposeTrackStatus(fStopAndKill);
				return &aParticleChange;
		}
		
		TLorentzVector p4an(aTrack.GetMomentum().x()/GeV,
						aTrack.GetMomentum().y()/GeV,
						aTrack.GetMomentum().z()/GeV,
						TMath::Sqrt((aTrack.GetMomentum().mag()/GeV)*(aTrack.GetMomentum().mag()/GeV) + m_Neutron*m_Neutron));

		TLorentzVector p4Nuclnbar = p4NucleusMother + p4an;


		TGenPhaseSpace event;
		Double_t masses2[3] = {m_mass_daughter, m_Proton, m_Proton};
		event.SetDecay(p4Nuclnbar, 3, masses2);
		int index  = 0;
		while(index < 1000){
				Double_t weight = event.Generate();  // weight < 0.4
				if(G4UniformRand() *0.7 < weight)
						break;
				index++;
		}
		TLorentzVector *p4dauNucl = event.GetDecay(0);
		TLorentzVector *p4pm12 = event.GetDecay(1);
		TLorentzVector *p4pm22 = event.GetDecay(2);

		theSecondary1 = G4AntiProton::AntiProton();
		theSecondary2 = G4Proton::Proton();
		theSecondary4 = nullptr;
		theSecondary4 = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z,A,0.0);

		G4ThreeVector p3pm12(p4pm12->Px()*GeV, 
						p4pm12->Py()*GeV,
						p4pm12->Pz()*GeV);
		G4ThreeVector p3pm22(p4pm22->Px()*GeV, 
						p4pm22->Py()*GeV,
						p4pm22->Pz()*GeV);
		G4ThreeVector p3dauNucl(p4dauNucl->Px()*GeV, 
						p4dauNucl->Py()*GeV,
						p4dauNucl->Pz()*GeV);


		aDynaSec1 = new G4DynamicParticle(theSecondary1, (p4pm12->E())*GeV, p3pm12);
		aDynaSec2 = new G4DynamicParticle(theSecondary2, (p4pm22->E())*GeV, p3pm22);
		aDynaSec4 = new G4DynamicParticle(theSecondary4, (p4dauNucl->E())*GeV, p3dauNucl);
		G4double  pathlenght = 0.00*mm;
		G4double  time0 = pathlenght/aTrack.GetVelocity();

		aSecTrk1  = new G4Track(aDynaSec1, aTrack.GetGlobalTime()+time0 , aTrack.GetPosition()  + p3pm12* time0);
		aSecTrk2  = new G4Track(aDynaSec2, aTrack.GetGlobalTime()+time0, aTrack.GetPosition() + p3pm22* time0 );
		//	aSecTrk3  = new G4Track(aDynaSec3, aTrack.GetGlobalTime()+time0, aTrack.GetPosition() + aTrack.GetPosition().unit() * +time0 );
		aSecTrk4  = new G4Track(aDynaSec4, aTrack.GetGlobalTime()+time0, aTrack.GetPosition() + p3dauNucl*time0);
	//	aSecTrk4->SetTrackStatus(fStopAndKill);
		if(debug == 1){
				std::cout << material->GetName() << std::endl;
				std::cout << "energy 1 : " << p4dauNucl->E() << std::endl;
				std::cout << "energy 2 : " << p4pm12->E() << std::endl;
				std::cout << "energy 3 : " << p4pm22->E() << std::endl;
				std::cout << std::fabs((aSecTrk1->GetMomentumDirection()).mag2()-1.0) << std::endl;
				std::cout << std::fabs((aSecTrk2->GetMomentumDirection()).mag2()-1.0) << std::endl;
				std::cout << std::fabs((aSecTrk4->GetMomentumDirection()).mag2()-1.0) << std::endl;
				std::cout << -1.0*(aSecTrk1->GetKineticEnergy())/MeV << std::endl;
				std::cout << -1.0*(aSecTrk2->GetKineticEnergy())/MeV << std::endl;
				std::cout << -1.0*(aSecTrk4->GetKineticEnergy())/MeV << std::endl;
				std::cout << (aTrack.GetGlobalTime() - aSecTrk1->GetGlobalTime())/ns << std::endl;
				std::cout << (aTrack.GetGlobalTime() - aSecTrk2->GetGlobalTime())/ns << std::endl;
				std::cout << (aTrack.GetGlobalTime() - aSecTrk4->GetGlobalTime())/ns << std::endl;
		}

		if(!(p4dauNucl->E() == p4dauNucl->E() && p4pm12->E() == p4pm12->E() && p4pm22->E() == p4pm22->E() &&
						std::fabs((aSecTrk1->GetMomentumDirection()).mag2()-1.0) == std::fabs((aSecTrk1->GetMomentumDirection()).mag2()-1.0) &&
						std::fabs((aSecTrk2->GetMomentumDirection()).mag2()-1.0) == std::fabs((aSecTrk2->GetMomentumDirection()).mag2()-1.0) &&
						std::fabs((aSecTrk4->GetMomentumDirection()).mag2()-1.0) == std::fabs((aSecTrk4->GetMomentumDirection()).mag2()-1.0) &&
						-1.0*(aSecTrk1->GetKineticEnergy())/MeV == -1.0*(aSecTrk1->GetKineticEnergy())/MeV &&
						-1.0*(aSecTrk2->GetKineticEnergy())/MeV == -1.0*(aSecTrk2->GetKineticEnergy())/MeV &&
						-1.0*(aSecTrk4->GetKineticEnergy())/MeV  == -1.0*(aSecTrk4->GetKineticEnergy())/MeV &&
						(aTrack.GetGlobalTime() - aSecTrk1->GetGlobalTime())/ns == (aTrack.GetGlobalTime() - aSecTrk1->GetGlobalTime())/ns &&
						(aTrack.GetGlobalTime() - aSecTrk2->GetGlobalTime())/ns == (aTrack.GetGlobalTime() - aSecTrk2->GetGlobalTime())/ns &&
						(aTrack.GetGlobalTime() - aSecTrk4->GetGlobalTime())/ns == (aTrack.GetGlobalTime() - aSecTrk4->GetGlobalTime())/ns)){

				aParticleChange.Initialize(aTrack);
				aParticleChange.ProposeEnergy(0.);
				aParticleChange.ProposeTrackStatus(fStopAndKill);
				return &aParticleChange;
		}


		aParticleChange.Initialize(aTrack);
		aParticleChange.ProposeEnergy(0.);
		aParticleChange.ProposeTrackStatus(fStopAndKill);
		aParticleChange.AddSecondary(aSecTrk1);
		aParticleChange.AddSecondary(aSecTrk2);
		aParticleChange.AddSecondary(aSecTrk4);
		return &aParticleChange;
}

G4bool G4AntiNeutronChargeExchangeProcess::IsApplicable(const G4ParticleDefinition& aParticleType){
		const G4ParticleDefinition* p = &aParticleType;
		return (p == theANeutron);
}









