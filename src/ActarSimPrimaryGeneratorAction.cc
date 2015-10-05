////////////////////////////////////////////////////////////////
//*-- AUTHOR : Hector Alvarez-Pol  (hapol@fpddux.usc.es)
//*-- Date: 11/2004
//*-- Last Update: 05/05/08
//*-- Modified by
// --------------------------------------------------------------
// Description:
//   Actions to perform to generate a primary vertex
//
// --------------------------------------------------------------
// Comments:
//   - 27/01/05 Cleaning and improving calculations
//   - 25/11/04 Created based on example/novice/N01 structure
//
// --------------------------------------------------------------
/////////////////////////////////////////////////////////////////

#include "ActarSimPrimaryGeneratorAction.hh"

#include "ActarSimDetectorConstruction.hh"
#include "ActarSimGasDetectorConstruction.hh"
#include "ActarSimPrimaryGeneratorMessenger.hh"

#include "ActarSimROOTAnalysis.hh"
#include "ActarSimCinePrimGenerator.hh"
#include "ActarSimKinePrimGenerator.hh"
#include "ActarSimEulerTransformation.hh"

#include "ActarSimBeamInfo.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4Gamma.hh"
#include "G4Geantino.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "G4ParticleDefinition.hh"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

ActarSimPrimaryGeneratorAction::ActarSimPrimaryGeneratorAction(const char *inputfile)
  :gasDetector(0), incidentIon(0),targetIon(0),scatteredIon(0),recoilIon(0),
   beamInteractionFlag("off"),
   realisticBeamFlag("off"), reactionFromEvGenFlag("off"), reactionFromCrossSectionFlag("off"),
   reactionFromFileFlag("off"),reactionFromCineFlag("off"),
   randomThetaFlag("off"),reactionFile("He8onC12_100MeV_Elastic.dat"),reactionFromKineFlag("off"),
   reactionFromCRYFlag("off"),
   vertexPosition(0) { 
  //
  // Constructor: init values are filled
  //

  //Initial Values
  G4ThreeVector zero;

  //create a particleGun
  particleGun = new G4ParticleGun(1);

  // Read the cry input file
  std::ifstream inputFile;
  inputFile.open(inputfile,std::ios::in);
  char buffer[1000];

  if (inputFile.fail()) {
    if( *inputfile !=0)  //....only complain if a filename was given
      G4cout << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputfile << G4endl;
    InputState=-1;
  }else{
    std::string setupString("");
    while ( !inputFile.getline(buffer,1000).eof()) {
      setupString.append(buffer);
      setupString.append(" ");
    }

    CRYSetup *setup=new CRYSetup(setupString,"cosmic/data");

    gen = new CRYGenerator(setup);

    // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState=0;
  }
  // create a vector to store the CRY particle properties
  vect=new std::vector<CRYParticle*>;

  // Create the table containing all particle names
  particleTable = G4ParticleTable::GetParticleTable();

  //create a messenger for this class
  gunMessenger = new ActarSimPrimaryGeneratorMessenger(this);
  /*
  G4ParticleDefinition* pd = particleTable->FindParticle("proton");
  if(pd != 0)
    particleGun->SetParticleDefinition(pd);
  particleGun->SetParticlePosition(zero);
  particleGun->SetParticleTime(0.0);
  particleGun->SetParticlePolarization(zero);
  particleGun->SetParticleCharge(1.0);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.0));
  particleGun->SetParticleEnergy(1*MeV);

  //create a pointer that gives access to the tabulated xs
  pReadEvGen = new ActarSimEventGenerator();

  reactionQ = 0.0001;   //does 0 work? (QM)

  labEnergy = 100;      // 15MeV*numero de nucleones (EI)

  incidentEnergy=labEnergy;
  // (EN) and (ENI) are taken from the target and projectile ion definitions

  thetaLabAngle = 45 * deg;   // 45 degrees (TH)

  randomThetaMin = 0.*deg;
  randomThetaMax = 180.*deg;*/
}

ActarSimPrimaryGeneratorAction::~ActarSimPrimaryGeneratorAction() {
  //
  // Simple destructor
  //

  delete particleGun;
  delete gunMessenger;
}

void ActarSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  
  //this function is called at the begining of event
  // Generate most of the primary physics. See the comments on each possible case
  //
  //REMOVE IF NOT NEEDED!!!!!!
  const G4int verboseLevel = G4RunManager::GetRunManager()->GetVerboseLevel();
  //G4cout << G4endl << " ______ VerboseLevel  _______" <<verboseLevel<< G4endl;

  G4ThreeVector zero;
//   G4double energyLostInTargetGas = 0; //zero, to be calculated if realisticBeamFlag is on

  //Initial values for reactionFromEvGen
  G4double  LabParticleAngle = 85.0 * deg;
  G4double  LabParticleAngle_rec = 85.0 * deg;
  if(LabParticleAngle_rec){;}

  ActarSimDetectorConstruction* detector = (ActarSimDetectorConstruction*)
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  //reading this parameter from the geometry
  if(!gasDetector) {
    // ActarSimDetectorConstruction* detector = (ActarSimDetectorConstruction*)
    //   (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    gasDetector = (ActarSimGasDetectorConstruction*)(detector->GetGasDetector());
    if(gasDetector->GetDetectorGeometry()=="tube")
      lengthParameter = gasDetector->GetLengthGasTub();
    else
      lengthParameter = gasDetector->GetZGasBox();

    G4cout << "##################################################################" << G4endl
	   << "##### ActarSimPrimaryGeneratorAction::GeneratePrimaries()  #######" << G4endl
	   << "##### Information: Length of gas volume = " << lengthParameter << "    ######"
           << G4endl;
    G4cout << "##################################################################" << G4endl;
  }

  //CASE ?  reactionFromCRYFlag = on
  // Reaction products kinematics calculated using Kine
  if(reactionFromCRYFlag == "on")
    {
    if (InputState != 0) {
      G4String* str = new G4String("CRY library was not successfully initialized");
      //G4Exception(*str);
      G4Exception("ActarSimPrimaryGeneratorAction", "1",
		  RunMustBeAborted, *str);
    }
    G4String particleName;
    vect->clear();
    gen->genEvent(vect);
    /*
    //....debug output
    G4cout << "\nEvent=" << anEvent->GetEventID() << " "
	   << "CRY generated nparticles=" << vect->size()
	   << G4endl;
    */
    
    G4int Inside=0;

    for ( unsigned j=0; j<vect->size(); j++) {

      if((*vect)[j]->x()>0.134 || (*vect)[j]->x()<-0.134)Inside++;
    }

    while(Inside!=0){
      Inside=0;
      vect->clear();
      gen->genEvent(vect);
      for ( unsigned j=0; j<vect->size(); j++) {

	if((*vect)[j]->x()>0.134 || (*vect)[j]->x()<-0.134)Inside++;
      }

    }

    for ( unsigned j=0; j<vect->size(); j++) {
      particleName=CRYUtils::partName((*vect)[j]->id());
      /*
      //....debug output  
      cout << "  "          << particleName << " "
	   << "charge="      << (*vect)[j]->charge() << " "
	   << setprecision(4)
	   << "energy (MeV)=" << (*vect)[j]->ke()*MeV << " "
	   << "pos (m)"
	   << G4ThreeVector((*vect)[j]->x(), (*vect)[j]->y(), (*vect)[j]->z())
	   << " " << "direction cosines "
	   << G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w())
	   << " " << endl;
      */
      particleGun->SetParticleDefinition(particleTable->FindParticle((*vect)[j]->PDGid()));
      particleGun->SetParticleEnergy((*vect)[j]->ke()*MeV);
      // particleGun->SetParticlePosition(G4ThreeVector((*vect)[j]->x()*m, (*vect)[j]->y()*m, (*vect)[j]->z()*m));
      // particleGun->SetParticleMomentumDirection(G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w()));
      particleGun->SetParticlePosition(G4ThreeVector((*vect)[j]->x()*m, (*vect)[j]->z()*m + 200.*mm, (*vect)[j]->y()*m + 69.*mm));
      //particleGun->SetParticlePosition(G4ThreeVector(-(*vect)[j]->x()*m, (*vect)[j]->z()*m + 200.*mm, -(*vect)[j]->y()*m + 64.*mm));
      particleGun->SetParticleMomentumDirection(G4ThreeVector(-(*vect)[j]->u(), (*vect)[j]->w(), -(*vect)[j]->v()));
      particleGun->SetParticleTime((*vect)[j]->t());
      particleGun->GeneratePrimaryVertex(anEvent);
      delete (*vect)[j];

  }

 }


  // CASE 5
  // Particle selected manually (using the messenger commands)
 else{
    if(verboseLevel>0){
      G4cout << G4endl
	     << " *************************************************** " << G4endl
	     << " * ActarSimPrimaryGeneratorAction::GeneratePrimaries() " << G4endl
	     << " * No particular event generator or kinematics code... " << G4endl
	     << " * A single particle is thrown using messenger commands."  << G4endl;
      G4cout << " *************************************************** "<< G4endl;
    }
    
    if(realisticBeamFlag == "on")
      //particleGun->SetParticlePosition(vertexPosition);
      particleGun->SetParticlePosition(G4ThreeVector());
    //G4ThreeVector particlePosition;
    //particlePosition=particleGun->GetParticlePosition();
    
    
    //DPLoureiro adding random distribution for polar and azimuthal angles
    G4double cosTheta;
    G4double sinTheta;
    //G4double Theta;
    G4double y_coord;
    y_coord = -1 + 2.0*G4UniformRand();
    y_coord=10*y_coord/185.;
   
    if(randomThetaFlag == "on") {
      
      G4ParticleDefinition* pd = particleTable->FindParticle("alpha");
      //G4ParticleDefinition* pd = particleTable->FindParticle("proton");
      if(pd != 0)
	particleGun->SetParticleDefinition(pd);
     
      //G4cout<<cos(randomThetaMin*rad)<<" "<<cos(randomThetaMax*rad)<<G4endl;
      
      G4double CosRandomThetaMin=cos(randomThetaMin);     
      G4double CosRandomThetaMax=cos(randomThetaMax);     
      
      if(CosRandomThetaMin==1. && CosRandomThetaMax==0. ) {
	cosTheta = -1.0 + 2.0*G4UniformRand();
       sinTheta = sqrt(1 - cosTheta*cosTheta);
       //G4cout<<"UNIFORM"<<G4endl;     
      }
      else{
	//Theta = (randomThetaMin + (randomThetaMax-randomThetaMin) * G4UniformRand()) * rad;
	cosTheta = cos(randomThetaMin*rad)+(cos(randomThetaMax*rad)- cos(randomThetaMin*rad))*G4UniformRand();
	//cosTheta = cos(Theta);
	sinTheta = sqrt(1 - cosTheta*cosTheta);
       //G4cout<<"NON UNIFORM"<<G4endl;     
     }
    }
    else{ 
      sinTheta=0;
      cosTheta=1;
    }
    G4double phi;

    //G4cout<< "THETA"<<" "<<Theta<<G4endl;     
    
    if(randomPhiFlag == "on"){
      G4double CosRandomPhiMin=cos(randomPhiMin);     
      G4double CosRandomPhiMax=cos(randomPhiMax);     
      if(CosRandomPhiMin==1. && CosRandomPhiMax==1.) 
	phi=2*pi*G4UniformRand();
      else 
	phi=randomPhiMin + ((randomPhiMax-randomPhiMin) * G4UniformRand()) * rad;
    }
  
   //else phi=pi/2;
   else phi=0;
    
   //G4cout<< randomPhiMin<<" "<<randomPhiMax<<" "<<phi<<G4endl;     
   //G4cout<< randomThetaMin<<" "<<randomThetaMax<<" "<<Theta<<G4endl;     
    
    if(alphaSourceFlag == "on"){
      G4int i=rand() % 3;
      G4double alpha_energy[3]={5.16,5.49,5.81};
      //G4cout<<alpha_energy[i]<<G4endl;
      particleGun->SetParticleEnergy(alpha_energy[i]*MeV);
      
    }
    
   else{  ;
     //particleGun->SetParticleEnergy(5.5*MeV);
     particleGun->SetParticleEnergy(GetLabEnergy());
        //incidentIon =  (G4Ions*) particleTable->GetIon(50, 134, 0.);  // 134Sn
     //incidentIon =  (G4Ions*) particleTable->GetIon(30, 80, 0.);  // 80Zn
     //incidentIon =  (G4Ions*) particleTable->GetIon(38, 90, 0.);  // 90Sr
     //incidentIonCharge =  50;
     //incidentIonCharge =  38;
     //particleGun->SetParticleDefinition(incidentIon);
   }

    //if(incidentIon)particleGun->SetParticleDefinition(incidentIon);
    //else particleGun->SetParticleDefinition(particleTable->FindParticle("mu+")); 

    //these are the cosines for an isotropic direction
    //particleGun -> SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    //Selecting only positive Z
    //Circle of 2.5 mm radius
    Float_t radius=beamRadiusAtEntrance;
    //Float_t radius=32*mm;
    //Float_t radius=1*cm;
    //Float_t radius=2.5*mm;
    //Float_t radius=1*mm;
    Float_t X0=0;
    Float_t Y0=0;
    //Float_t Z0= -(detector->GetChamberZLength()-detector->GetZGasBoxPosition());
    //G4double Z0= -29*mm;

    do{
      X0=-1.0 + 2.0*G4UniformRand();
      Y0=-1.0 + 2.0*G4UniformRand();
    }while((pow(X0,2)+pow(Y0,2)>1));
    X0=radius*X0;
    Y0=radius*Y0;
    //G4cout<<X0<<" "<<Y0<<endl;
    //particleGun->SetParticlePosition(G4ThreeVector(0.,0.,Z0));
    //particleGun->SetParticlePosition(G4ThreeVector(X0,Y0,Z0));
    //particleGun->SetParticlePosition(G4ThreeVector(X0,Y0,0.));
    //particleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    //particleGun->SetParticlePosition(G4ThreeVector(0.,200.,50.));
    particleGun->SetParticlePosition(GetParticlePosition());
    if(cosTheta>0)
      //particleGun -> SetParticleMomentumDirection(G4ThreeVector(y_coord,sinTheta*sin(phi),cosTheta));
      //particleGun -> SetParticleMomentumDirection(G4ThreeVector(sinTheta*cos(phi),sinTheta*sin(phi),cosTheta));
      //particleGun -> SetParticleMomentumDirection(G4ThreeVector(0,-1,0));
      particleGun -> SetParticleMomentumDirection(GetParticleMomentumDirection());
    //particleGun -> SetParticleMomentumDirection(G4ThreeVector(sinTheta*cos(phi),y_coord,cosTheta));
    else
      //particleGun -> SetParticleMomentumDirection(G4ThreeVector(-y_coord,sinTheta*sin(phi),-cosTheta));
      //particleGun -> SetParticleMomentumDirection(G4ThreeVector(-sinTheta*cos(phi),sinTheta*sin(phi),-cosTheta));
    //particleGun -> SetParticleMomentumDirection(G4ThreeVector(-sinTheta*cos(phi),y_coord,-cosTheta));
    //particleGun -> SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    //particleGun -> SetParticleMomentumDirection(G4ThreeVector(0,-1,0));
    particleGun -> SetParticleMomentumDirection(GetParticleMomentumDirection());
    particleGun->GeneratePrimaryVertex(anEvent);
  
  }


  // G4ParticleDefinition *pdef=particleGun->GetParticleDefinition();   
    // G4cout << " *************************************************** "<< G4endl;
    // G4cout<<pdef->GetAtomicNumber()<<" "<<pdef->GetAtomicMass()<<" "<<particleGun->GetParticleEnergy()/MeV<<endl;
    // G4cout << " *************************************************** "<< G4endl;


  //TODO SOlve this assymetry No theta or energies should be in the argument!
  
  //G4cout<<"HERE I AM!!!!!!!"<<G4endl;  
  ActarSimBeamInfo *pBeamInfo = (ActarSimBeamInfo*) 0;
  if(gActarSimROOTAnalysis){
    //D. Perez
    pBeamInfo = gActarSimROOTAnalysis->GetBeamInfo();
    //G4cout << pBeamInfo << G4endl
    //	   << "Theta1 "<< pBeamInfo->GetThetaEntrance() / deg << " deg"<<G4endl
    //	   << "Theta2 "<< pBeamInfo->GetThetaVertex() / deg << " deg"<<G4endl;
    //Histogramming
    gActarSimROOTAnalysis->GeneratePrimaries(anEvent,pBeamInfo);
    //gActarSimROOTAnalysis->GeneratePrimaries(anEvent,
    //					     theta1,
    //					     theta2,
    //					     energy1,
    //					     energy2);
  }


}

//----------------------------------------------------------------------------//
void ActarSimPrimaryGeneratorAction::InputCRY()
{
  InputState=1;
}

//----------------------------------------------------------------------------//
void ActarSimPrimaryGeneratorAction::UpdateCRY(std::string* MessInput)
{
  CRYSetup *setup=new CRYSetup(*MessInput,"cosmic/data");

  gen = new CRYGenerator(setup);

  // set random number generator
  RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
  setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
  InputState=0;

}

//----------------------------------------------------------------------------//
void ActarSimPrimaryGeneratorAction::CRYFromFile(G4String newValue)
{
  // Read the cry input file
  std::ifstream inputFile;
  inputFile.open(newValue,std::ios::in);
  char buffer[1000];

  if (inputFile.fail()) {
    G4cout << "Failed to open input file " << newValue << G4endl;
    G4cout << "Make sure to define the cry library on the command line" << G4endl;
    InputState=-1;
  }else{
    std::string setupString("");
    while ( !inputFile.getline(buffer,1000).eof()) {
      setupString.append(buffer);
      setupString.append(" ");
    }

    CRYSetup *setup=new CRYSetup(setupString,"cosmic/data");

    gen = new CRYGenerator(setup);

  // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState=0;
  }
}

