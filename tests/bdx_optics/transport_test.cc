// Native Geant4 transport checks share the exact GEMC optics implementation.
#include "BdxCrystalOptics.h"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4GenericTrap.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4Electron.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include <iostream>
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4ParticleGun.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include <string>
#include <vector>

long firstDiffuse=0,firstSpecular=0,firstAbsorbed=0,birthInjections=0;
bool injectBirth=false,birthInjectedThisEvent=false;

class Detector final : public G4VUserDetectorConstruction {
    G4VPhysicalVolume* Construct() override {
        auto* nist=G4NistManager::Instance();
        auto* worldLV=new G4LogicalVolume(new G4Box("world",300*mm,300*mm,300*mm),
            nist->FindOrBuildMaterial("G4_AIR"),"world");
        auto* world=new G4PVPlacement(nullptr,{},worldLV,"world",nullptr,false,0);
        std::vector<G4TwoVector> vertices={{-10,-10},{-10,10},{10,10},{10,-10},
                                         {-10,-10},{-10,10},{10,10},{10,-10}};
        auto* bgoLV=new G4LogicalVolume(new G4GenericTrap("bgo",120*mm,vertices),
            nist->FindOrBuildMaterial("G4_BGO"),"bgo");
        auto* bgo=new G4PVPlacement(nullptr,{},bgoLV,"bdx_bgo",worldLV,false,0);
        auto* sipmLV=new G4LogicalVolume(new G4Box("sipm",1.5*mm,1.5*mm,0.05*mm),
            nist->FindOrBuildMaterial("G4_Galactic"),"sipm");
        auto* sipm=new G4PVPlacement(nullptr,{0,0,120.05*mm},sipmLV,"bdx_sipm",worldLV,false,0);
        bdxoptics::Construct(world,bgo,sipm);
        return world;
    }
};
class Physics final : public G4VModularPhysicsList {
public:
    Physics() { RegisterPhysics(new G4DecayPhysics); RegisterPhysics(new G4EmStandardPhysics);
        RegisterPhysics(new G4OpticalPhysics);bdxoptics::ConfigurePhysics(); }
    void ConstructProcess() override { G4VModularPhysicsList::ConstructProcess(); bdxoptics::InstallBoundaryProcess(); }
    void SetCuts() override {}
};
class Generator final : public G4VUserPrimaryGeneratorAction {
    G4ParticleGun gun{1};
public:
    explicit Generator(const std::string& mode) {
        gun.SetParticleDefinition(G4OpticalPhoton::Definition());
        gun.SetParticleEnergy(2.583004134*eV);
        G4ThreeVector direction(0,0,1),position(0,0,0);
        if(mode=="outside") position.setX(5*mm);
        if(mode=="one_reflection") direction={20,0,120};
        if(mode=="many_reflections") direction={2020,0,120};
        if(mode=="mixed") direction={10,0,30};
        if(mode=="far_side") {position={0,0,-100};direction={1,0,0.01};}
        if(mode=="secondary_birth") {
            injectBirth=true;
            gun.SetParticleDefinition(G4Electron::Definition());
            gun.SetParticleEnergy(1*MeV);
            direction={1,0,0};
        }
        gun.SetParticlePosition(position);
        gun.SetParticleMomentumDirection(direction.unit());
        gun.SetParticlePolarization({0,1,0});
    }
    void GeneratePrimaries(G4Event* event) override { gun.GeneratePrimaryVertex(event); }
};
class Event final : public G4UserEventAction {
    void BeginOfEventAction(const G4Event* event) override { birthInjectedThisEvent=false; bdxoptics::BeginEvent(event); }
    void EndOfEventAction(const G4Event* event) override { bdxoptics::EndEvent(event); }
};
class Step final : public G4UserSteppingAction {
    void UserSteppingAction(const G4Step* step) override {
        if(injectBirth && !birthInjectedThisEvent &&
           step->GetTrack()->GetDefinition()!=G4OpticalPhoton::Definition() &&
           step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="bdx_bgo") {
            for(const auto* child:*step->GetSecondaryInCurrentStep()) {
                if(child->GetDefinition()!=G4OpticalPhoton::Definition()) continue;
                // Emulate a newborn photon's rounded position and wrong BGO
                // touchable, before the production hook sees the secondary.
                auto* newborn=const_cast<G4Track*>(child);
                newborn->SetPosition({(10+1.28e-6)*mm,0,0});
                auto* worldHistory=new G4TouchableHistory(*static_cast<const G4TouchableHistory*>(step->GetPreStepPoint()->GetTouchable()));
                worldHistory->MoveUpHistory(1);
                G4TouchableHandle outside=worldHistory;
                newborn->SetTouchableHandle(outside);
                newborn->SetNextTouchableHandle(outside);
                birthInjectedThisEvent=true;++birthInjections;break;
            }
        }
        bdxoptics::Step(step);
        if(step->GetTrack()->GetCurrentStepNumber()==1 &&
           step->GetPostStepPoint()->GetStepStatus()==fGeomBoundary) {
            auto* processes=G4OpticalPhoton::Definition()->GetProcessManager()->GetProcessList();
            for(std::size_t i=0;i<processes->size();++i)
                if(auto* p=dynamic_cast<G4OpBoundaryProcess*>((*processes)[i])) {
                    if(p->GetStatus()==LambertianReflection) ++firstDiffuse;
                    if(p->GetStatus()==SpikeReflection) ++firstSpecular;
                    if(p->GetStatus()==Absorption) ++firstAbsorbed;
                }
        }
    }
};
int main(int argc,char** argv) {
    if(argc!=3) return 2;
    G4Random::setTheSeed(758629);
    auto* run=new G4RunManager;
    run->SetUserInitialization(new Detector);
    run->SetUserInitialization(new Physics);
    run->SetUserAction(new Generator(argv[1]));
    run->SetUserAction(new Event);
    run->SetUserAction(new Step);
    run->Initialize();
    run->BeamOn(std::stoi(argv[2]));
    std::cout << "FIRST_STATS " << firstDiffuse << " " << firstSpecular << " " << firstAbsorbed << "\n";
    std::cout << "BIRTH_INJECTIONS " << birthInjections << "\n";
    delete run;
}
