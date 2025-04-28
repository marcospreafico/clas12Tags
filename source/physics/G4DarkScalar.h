/*
 * G4Scalar.h
 *
 *  Created on: April 28, 2025
 *      Author: msprea
 */

#ifndef PHYSICS_G4DARKSCALAR_H_
#define PHYSICS_G4DARKSCALAR_H_

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                       G4DarkScalar                                 ###
// ######################################################################

class G4DarkScalar: public G4ParticleDefinition {
private:
	static G4DarkScalar* theInstance;
	G4DarkScalar() {
	}
	~G4DarkScalar() {
	}

	/*Since this is a singleton, the mass is used only at the first call*/
public:
	static G4DarkScalar* Definition(G4double mass=0., G4double g_mu=0);
	static G4DarkScalar* DarkScalarDefinition(G4double mass=0., G4double g_mu=0);
	static G4DarkScalar* DarkScalar(G4double mass=0., G4double g_mu=0);
};

#endif /* PHYSICS_G4DARKSCALAR_H_ */
