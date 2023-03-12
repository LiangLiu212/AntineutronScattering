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
/// \file PrimaryGenerator.cc
/// \brief Implementation of the PrimaryGenerator1 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGenerator.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
using namespace CLHEP;
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::PrimaryGenerator()
: G4VPrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGenerator::SetRootFile(TString str){
		f1 = new TFile(str, "read");
		t1 = (TTree*)f1->Get("mctruth");

		t1->SetBranchAddress("mcpxpim", &mcpxpim);
		t1->SetBranchAddress("mcpypim", &mcpypim);
		t1->SetBranchAddress("mcpzpim", &mcpzpim);
		t1->SetBranchAddress("mcepim", &mcepim);
		t1->SetBranchAddress("inipxpim", &inipxpim);
		t1->SetBranchAddress("inipypim", &inipypim);
		t1->SetBranchAddress("inipzpim", &inipzpim);
		t1->SetBranchAddress("iniepim", &iniepim);

		t1->SetBranchAddress("mcpxnbar", &mcpxnbar);
		t1->SetBranchAddress("mcpynbar", &mcpynbar);
		t1->SetBranchAddress("mcpznbar", &mcpznbar);
		t1->SetBranchAddress("mcenbar", &mcenbar);
		t1->SetBranchAddress("inipxnbar", &inipxnbar);
		t1->SetBranchAddress("inipynbar", &inipynbar);
		t1->SetBranchAddress("inipznbar", &inipznbar);
		t1->SetBranchAddress("inienbar", &inienbar);

		t1->SetBranchAddress("mcpxp", &mcpxp);
		t1->SetBranchAddress("mcpyp", &mcpyp);
		t1->SetBranchAddress("mcpzp", &mcpzp);
		t1->SetBranchAddress("mcep", &mcep);
		t1->SetBranchAddress("inipxp", &inipxp);
		t1->SetBranchAddress("inipyp", &inipyp);
		t1->SetBranchAddress("inipzp", &inipzp);
		t1->SetBranchAddress("iniep", &iniep);

		fweight  = new TFile("/home/lliu/ustcfs/runarea/neutron/BesPip/myBesPip/AntiNeutronChargeExchange/run/EvtGen/reweight/weight.root", "read");
		hweight = (TH2D*)fweight->Get("hweight");


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{
		//vertex A uniform on a cylinder
		//

		Double_t m_Jpsi = 3.097;
		Double_t m_Proton = 0.938272081;
		Double_t m_Neutron = 0.939565413;
		Double_t m_Pion = 0.13957039;

		TLorentzVector ep, em, cms;
		em.SetPxPyPzE(0.011*m_Jpsi/2. ,0. , sqrt((pow(m_Jpsi,2)/4.)-(0.000511*0.000511)),( m_Jpsi/( sqrt(1-0.011*0.011) ) )/2. );
		ep.SetPxPyPzE(0.011*m_Jpsi/2. ,0. ,-sqrt((pow(m_Jpsi,2)/4.)-(0.000511*0.000511)),( m_Jpsi/( sqrt(1-0.011*0.011) ) )/2. );
		cms = ep+em;

		TGenPhaseSpace genphsp;

		Double_t masses[3] = {m_Proton, m_Neutron, m_Pion};
		genphsp.SetDecay(cms, 3, masses);

		int index = 0;
		while(index < 1000){
				Double_t weight = genphsp.Generate();  // weight < 0.4
				if(G4UniformRand() *0.7 < weight)
						break;
				index++;
		}
		TLorentzVector *p4pp = genphsp.GetDecay(0);
		TLorentzVector *p4an = genphsp.GetDecay(1);
		TLorentzVector *p4pim = genphsp.GetDecay(2);



		const G4double r = 2*mm;
		const G4double zmax = 5*mm;
		//
		G4double alpha = twopi*G4UniformRand();     //alpha uniform in (0, 2*pi)
		G4double ux = std::cos(alpha);
		G4double uy = std::sin(alpha);
		G4double z = zmax*(2*G4UniformRand() - 1);  //z uniform in (-zmax, +zmax)

		G4double beamPosX = 2.1022400000000001086;
		G4double beamPosY = -1.1219099999999999628;
		G4double beamPosZ  = 6.8130999999999994898;
		G4double beamSizeX = 0.38000000000000000444;
		G4double beamSizeY = 0.0057000000000000002054;
		G4double beamSizeZ = 7.4576099999999998502;

		double tmpPosX =  (beamPosX + beamSizeX* G4RandGauss::shoot() ) *mm;
		double tmpPosY =  (beamPosY + beamSizeY* G4RandGauss::shoot() ) *mm;
		double tmpPosZ =  (beamPosZ + beamSizeZ* G4RandGauss::shoot() ) *mm;

		G4ThreeVector positionA(tmpPosX, tmpPosY, tmpPosZ);
		G4double timeA = 0*s;
		// 
		G4PrimaryVertex* vertexA = new G4PrimaryVertex(positionA, timeA);

		//particle 1 at vertex A
		//
		//	G4ParticleDefinition* particleDefinition
		//			= G4ParticleTable::GetParticleTable()->FindParticle("geantino");
		G4ParticleDefinition* particleDefinition
				= G4ParticleTable::GetParticleTable()->FindParticle("proton");
		G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition);
		//		particle1->SetMomentumDirection(G4ThreeVector(ux,uy,0));    
		//		particle1->SetKineticEnergy(1*MeV);
		particle1->Set4Momentum(p4pp->Px()*GeV, p4pp->Py()*GeV, p4pp->Pz()*GeV, p4pp->E()*GeV);
		//
		vertexA->SetPrimary(particle1);

		//particles 2 and 3 at vertex B
		//
		G4ParticleDefinition* particleDefinition2
				= G4ParticleTable::GetParticleTable()->FindParticle("anti_neutron");
		G4PrimaryParticle* particle2 = new G4PrimaryParticle(particleDefinition2);
		//	particle2->SetMomentumDirection(G4ThreeVector(ux,uy,0));    
		//	particle2->SetKineticEnergy(1*keV);
		particle2->Set4Momentum(p4an->Px()*GeV, p4an->Py()*GeV, p4an->Pz()*GeV, p4an->E()*GeV);

		G4ParticleDefinition* particleDefinition3
				= G4ParticleTable::GetParticleTable()->FindParticle("pi-");
		particleDefinition3->SetPDGStable(true);
		G4PrimaryParticle* particle3 = new G4PrimaryParticle(particleDefinition3);
		//	particle3->SetMomentumDirection(G4ThreeVector(ux,uy,0));    
		//	particle3->SetKineticEnergy(1*GeV);
		particle3->Set4Momentum(p4pim->Px()*GeV, p4pim->Py()*GeV, p4pim->Pz()*GeV, p4pim->E()*GeV);
		//
		vertexA->SetPrimary(particle2);
		vertexA->SetPrimary(particle3);  
		event->AddPrimaryVertex(vertexA);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertexRootFromEvtGen(G4Event* event)
{
		int evt = event->GetEventID();
//		cout << "Event ID : " << evt << endl;
		t1->GetEntry(evt);

		TLorentzVector p4nbar;
		p4nbar.SetPxPyPzE(mcpxnbar, mcpynbar, mcpznbar, mcenbar);
		double nbar_rho  =p4nbar.Rho();
		double nbar_cos = p4nbar.CosTheta();
		int idx_rho = floor(nbar_rho / (1.2 / 100.));
		int idx_cos = floor( (nbar_cos +1) / (2. / 100.) );
		double weight =  hweight->GetBinContent(idx_cos + 1, idx_rho + 1);

		G4ParticleDefinition* particleDefinition1
				= G4ParticleTable::GetParticleTable()->FindParticle("proton");
		G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition1);
		particle1->Set4Momentum(mcpxp*GeV, mcpyp*GeV, mcpzp*GeV, mcep*GeV);

		G4ParticleDefinition* particleDefinition2
				= G4ParticleTable::GetParticleTable()->FindParticle("anti_neutron");
		G4PrimaryParticle* particle2 = new G4PrimaryParticle(particleDefinition2);
		particle2->Set4Momentum(mcpxnbar*GeV, mcpynbar*GeV, mcpznbar*GeV, mcenbar*GeV);

		G4ParticleDefinition* particleDefinition3
				= G4ParticleTable::GetParticleTable()->FindParticle("pi-");
		G4PrimaryParticle* particle3 = new G4PrimaryParticle(particleDefinition3);
		particle3->Set4Momentum(mcpxpim*GeV, mcpypim*GeV, mcpzpim*GeV, mcepim*GeV);


		G4ThreeVector positionA(inipxpim*cm, inipypim*cm, inipzpim*cm);

		G4PrimaryVertex* vertexA = new G4PrimaryVertex(positionA, 0);
		vertexA->SetPrimary(particle1);
		if(G4UniformRand() < weight){
				vertexA->SetPrimary(particle2);
		}
		vertexA->SetPrimary(particle3);
		event->AddPrimaryVertex(vertexA);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
