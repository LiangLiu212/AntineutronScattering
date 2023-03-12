/*****************************************************************************\
 *  BOSS (BES Offline Software System)                                       *
 *       Filename:  G4AntiNeutronChargeExchange.cc
 *                                                                           *
 *    Description:                                                           
 *                                                                           *
 *        Version:  1.0                                                      *
 *        Created:  10/12/2022 04:13:34 PM                                   *
 *       Revision:  none                                                     *
 *       Compiler:  gcc                                                      *
 *                                                                           *
 *         Author:  Liang Liu (USTC), <liangzy@mail.ustc.edu.cn>             *
 *   Organization:  BESIII Collaboration                                     *
 *        License:  GNU Lesser General Public License 3.0, see LICENSE.md.   *
\*****************************************************************************/


// G4AntiNeutronChargeExchange.hh 

// process anti-n + n -> anti-p + p
//
// process anti-n + n + p  -> anti-p + p + p


#include "G4AntiNeutronChargeExchange.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4NucleiProperties.hh"

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

#include "G4HadronicParameters.hh"


G4AntiNeutronChargeExchange::G4AntiNeutronChargeExchange() : G4HadronicInteraction("Anti-neutron Charge Exchange")
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );

  lowestEnergyLimit    = 1.*MeV;  

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

G4AntiNeutronChargeExchange::~G4AntiNeutronChargeExchange()
{}

G4HadFinalState* G4AntiNeutronChargeExchange::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();
  const G4HadProjectile* aParticle = &aTrack;
  G4double ekin = aParticle->GetKineticEnergy();

  G4int A = targetNucleus.GetA_asInt();
  G4int Z = targetNucleus.GetZ_asInt();

  if(ekin <= lowestEnergyLimit || A < 3) {
    theParticleChange.SetEnergyChange(ekin);
    theParticleChange.SetMomentumChange(0.0,0.0,1.0);
    return &theParticleChange;
  }

  G4double plab = aParticle->GetTotalMomentum();

  if (verboseLevel > 1)
    G4cout << "G4AntiNeutronChargeExchange::DoIt: Incident particle plab="
	   << plab/GeV << " GeV/c "
	   << " ekin(MeV) = " << ekin/MeV << "  "
	   << aParticle->GetDefinition()->GetParticleName() << G4endl;

  // Scattered particle referred to axis of incident particle
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  G4int N = A - Z;
  G4int projPDG = theParticle->GetPDGEncoding();
  if (verboseLevel > 1)
    G4cout << "G4AntiNeutronChargeExchange for " << theParticle->GetParticleName()
	   << " PDGcode= " << projPDG << " on nucleus Z= " << Z
	   << " A= " << A << " N= " << N
	   << G4endl;

  const G4ParticleDefinition* theDef = nullptr;

  G4double mass2 = G4NucleiProperties::GetNuclearMass(A, Z);
  G4LorentzVector lv1 = aParticle->Get4Momentum();  // anti-neutron beam momentum
  G4LorentzVector lv0(0.0,0.0,0.0,mass2); 			// particle at rest

	// try to use TGenPhaseSpace to generate the final states
	
  // Sample final particles
 			const G4ParticleDefinition* theRecoilNucleus = nullptr; // for both p pbar and p p pbar
 			const G4ParticleDefinition* theSecondary1 = nullptr;    // for both p pbar and p p pbar
 			const G4ParticleDefinition* theSecondary2 = nullptr;    // for both p pbar and p p pbar
 			const G4ParticleDefinition* theSecondary3 = nullptr;    // for p p pbar

  if(G4UniformRand() > 0.25){   // simulate process nnbar + N -> ppbar + N'

			

  }
  else {    // simulate process nnbar + N  -> pppbar + N'



  }
	
	
	
  	return &theParticleChange;
}

G4double G4AntiNeutronChargeExchange::SampleT(G4double tmax, G4int A)
{

}

