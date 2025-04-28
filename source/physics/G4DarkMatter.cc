/*
 * G4DarkMatter.cc
 *
 *  Created on: March 15, 2020
 *      Author: celentan
 */

#include "G4DarkMatter.h"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"



// ######################################################################
// ###                         G4DarkMatter                               ###
// ######################################################################

G4DarkMatter* G4DarkMatter::theInstance = 0;
G4DarkMatter* G4DarkMatter::Definition(G4double mass,G4int twoJ)
{

  if (theInstance !=0) return theInstance;
  const G4String name = "DarkMatter";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);



  if (anInstance ==0)
  {
  // create particle
  //
  //    Arguments for constructor are as follows

      int PID = 911;
      std::string type = "";
      if(twoJ == 0) {
          PID = 811;
          type = "scalar";
      }
      if(twoJ == 1){
          PID = 911;
          type = "fermion";
      }
          
      
      
      
  anInstance = new G4ParticleDefinition(
                                        
//               name             mass          width         charge
                 name,	          mass,         	0,	           0,
//             2*spin           parity  C-conjugation
				 twoJ,			    1,	        	0,
//          2*Isospin       2*Isospin3       G-parity
                    0,			     0,             0,
//               type    lepton number  baryon number   PDG encoding
				 type,               0,             0,          PID,
//             stable         lifetime    decay table
				 true,              0,           NULL,
//             shortlived      subType    anti_encoding
                    false,         type
              );


  }
  theInstance = reinterpret_cast<G4DarkMatter*>(anInstance);
  return theInstance;
}

G4DarkMatter*  G4DarkMatter::DarkMatterDefinition(G4double mass,G4int twoJ)
{
  return Definition(mass,twoJ);
}

G4DarkMatter*  G4DarkMatter::DarkMatter(G4double mass,G4int twoJ)
{
  return Definition(mass,twoJ);
}
