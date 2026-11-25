#ifndef BDX_CRYSTAL_OPTICS_H
#define BDX_CRYSTAL_OPTICS_H
#include "G4ThreeVector.hh"
class G4VSolid;
class G4GenericTrap;
class G4VPhysicalVolume;
class G4Step;
class G4Event;
namespace hipo { class dictionary; class event; }
namespace bdxoptics {
bool Enabled();
void RegisterHipoSchemas(hipo::dictionary&);
void AppendHipo(hipo::event&, hipo::dictionary&, int gemcEvent);
void Construct(G4VPhysicalVolume* world, G4VPhysicalVolume* crystal, G4VPhysicalVolume* sensor);
void ConfigurePhysics();
void InstallBoundaryProcess();
double SpecularFraction(double distance_mm);
// Bounded correction for numerical intersection errors; push uses Geant4 units.
G4ThreeVector RelocateReflectedPoint(const G4VSolid&, const G4ThreeVector&, double push);
G4VSolid* MakeConservativeBgoSolid(const G4GenericTrap&);
void BeginEvent(const G4Event*);
void EndEvent(const G4Event*);
// Returns true for optical tracks owned by this study: bypass legacy GEMC cuts.
bool Step(const G4Step*);
}
#endif
