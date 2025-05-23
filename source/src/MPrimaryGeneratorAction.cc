// G4 headers
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"

// gemc headers
#include "MPrimaryGeneratorAction.h"
#include "string_utilities.h"

// mlibrary
#include "gstring.h"

using namespace gstring;

// C++ headers
#include <sys/types.h>
#include <dirent.h>
#include <iostream>

//ROOT headers
#include "TH3D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TTree.h"

using namespace std;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"

using namespace CLHEP;

MPrimaryGeneratorAction::MPrimaryGeneratorAction(goptions *opts) {
    gemcOpt = opts;
    hd_msg = gemcOpt->optMap["LOG_MSG"].args + " Beam Settings >> ";
    input_gen = gemcOpt->optMap["INPUT_GEN_FILE"].args;
    background_gen = gemcOpt->optMap["MERGE_LUND_BG"].args;
    cosmics = gemcOpt->optMap["COSMICRAYS"].args;
    GEN_VERBOSITY = gemcOpt->optMap["GEN_VERBOSITY"].arg;
    ntoskip = gemcOpt->optMap["SKIPNGEN"].arg;
    PROPAGATE_DVERTEXTIME = gemcOpt->optMap["PROPAGATE_DVERTEXTIME"].arg;

    muonDecay = 0;
    muonDecay = gemcOpt->optMap["FORCE_MUON_RADIATIVE_DECAY"].arg; //M.S.
    
    particleTable = G4ParticleTable::GetParticleTable();

    beamPol = 0;

    setBeam();

    particleGun = new G4ParticleGun(1);

    if (input_gen == "gemc_internal") {
        vector <string> cvalues = get_info(gemcOpt->optMap["COSMICAREA"].args, string(",\""));

        if (cvalues.size() < 4) {
            cout << "  !!!  Warning:  COSMICAREA flag not set correctly. It should be 4 numbers: x,y,z and R." << endl;
        }

        cosmicTarget = G4ThreeVector(get_number(cvalues[0]), get_number(cvalues[1]), get_number(cvalues[2]));
        cosmicRadius = get_number(cvalues[3]);
        if (cvalues.size() == 5) {
            cosmicGeo = cvalues[4];
        } else {
            cosmicGeo = "sph";
        }

        if (cosmics == "no") {
            cout << endl << hd_msg << " Beam Type: " << Particle->GetParticleName() << endl;
            cout << hd_msg << " Beam Momentum: " << G4BestUnit(mom, "Energy");
            if (dmom > 0) cout << " +- " << G4BestUnit(dmom, "Energy");
            cout << endl;
            cout << hd_msg << " Beam Direction: (theta, phi) = (" << theta / deg << ", " << phi / deg << ") deg";
            if (dtheta > 0 || dphi > 0) cout << " +- (" << dtheta / deg << ", " << dphi / deg << ") deg ";
            cout << endl << hd_msg << " Beam Vertex: (" << vx / cm << ", " << vy / cm << ", " << vz / cm << ") cm";
            if (drdzOrdxdydz == 0) cout << " (radius, z-spread) = (" << dvr / cm << ", " << dvz / cm << ") cm, ";
            else cout << " (dvx, dvy, dvz) = (" << dvx / cm << ", " << dvy / cm << ", " << dvz / cm << ") cm, ";
            cout << (gaussOrFlatV ? "gaussian spread" : "flat spread");
            cout << endl << hd_msg << " Beam polarization: " << polDeg << "%" << endl;
            cout << hd_msg << " Polarization Direction: (theta, phi) = (" << polTheta / deg << ", " << polPhi / deg << ")";
            cout << endl;
        } else {
            cout << endl << hd_msg << " Beam Type: Cosmic rays." << endl;
            if(cosmicParticle != "muon") cout << hd_msg << "Beam parameters follow: A. Dar, Phys.Rev.Lett, 51,3,p.227 (1983) " << endl;
            else cout << hd_msg << "Beam parameters follow: 1509.06176 " << endl;
            cout << hd_msg << " Momentum Range: [" << cminp / GeV << " - " << cmaxp / GeV << "] GeV" << endl;
            cout << hd_msg << " Cosmic Area :" << cosmicTarget << endl;
            cout << hd_msg << " Cosmic Radius :" << cosmicRadius / cm << " cm " << endl;
            cout << hd_msg << " Cosmic Surface Type: " << cosmicGeo << endl;
            cout << hd_msg << " Cosmic Particle Type: " << cosmicParticle << endl;
            if (cosmicParticle == "muon" && muonDecay == 1) cout << hd_msg << " Radiative muon decay BR=100%" << endl; //M.S.
        }
    }


    if (NP > 0) {
        cout << endl << hd_msg << " Luminosity Particle Type: " << L_Particle->GetParticleName() << endl;
        cout << hd_msg << " Luminosity Particle Momentum: " << G4BestUnit(L_mom, "Energy");
        if (L_dmom > 0) cout << " +- " << G4BestUnit(L_dmom, "Energy");
        cout << endl;
        cout << hd_msg << " Luminosity Particle Direction: (theta, phi) = (" << L_theta / deg << ", " << L_phi / deg << ") deg";
        if (L_dtheta > 0 || L_dphi > 0) cout << " +- (" << L_dtheta / deg << ", " << L_dphi / deg << ") deg";
        cout << endl << hd_msg << " Luminosity Particle Vertex: (" << L_vx / cm << ", " << L_vy / cm << ", " << L_vz / cm << ") cm";
        if (L_dvr + L_dvz > 0) cout << " (radius, z-spread) = (" << L_dvr / cm << ", " << L_dvz / cm << ")";
        cout << endl << hd_msg << " Number of Luminosity Particles: " << NP << endl;
        cout << hd_msg << " Luminosity Time Window: " << TWINDOW / ns << " nanoseconds." << endl;
        cout << hd_msg << " Luminosity Time Between Bunches: " << TBUNCH / ns << " nanoseconds." << endl;
        cout << hd_msg << " Primary Beam Particle Bunch Time: " << TSIGNAL / ns << " nanoseconds." << endl;
    }

    if (NP2 > 0) {
        cout << endl << hd_msg << " Luminosity Particle 2 Type: " << L2_Particle->GetParticleName() << endl;
        cout << hd_msg << " Luminosity Particle 2 Momentum: " << G4BestUnit(L2_mom, "Energy");
        if (L2_dmom > 0) cout << " +- " << G4BestUnit(L2_dmom, "Energy");
        cout << endl;
        cout << hd_msg << " Luminosity Particle 2 Direction: (theta, phi) = (" << L2_theta / deg << ", " << L2_phi / deg << ") deg";
        if (L2_dtheta > 0 || L2_dphi > 0) cout << " +- (" << L2_dtheta / deg << ", " << L2_dphi / deg << ") deg";
        cout << endl << hd_msg << " Luminosity Particle Vertex: (" << L2_vx / cm << ", " << L2_vy / cm << ", " << L2_vz / cm << ") cm";
        if (L2_dvr + L2_dvz > 0) cout << " (radius, z-spread) = (" << L2_dvr / cm << ", " << L2_dvz / cm << ") cm";
        cout << endl << hd_msg << " Number of Luminosity Particles 2: " << NP2 << endl;
        cout << hd_msg << " Luminosity Time Between Bunches: " << TBUNCH2 / ns << " nanoseconds." << endl;
    }

    eventIndex = 1;

    // Set up to read saved RNGs

    string arg = gemcOpt->optMap["RERUN_SELECTED"].args;
    if (arg == "" || arg == "no") {
        rsp.enabled = false;
    } else {
        vector <string> values;
        string units;
        values = get_info(gemcOpt->optMap["RERUN_SELECTED"].args, string(",\""));
        if (values.size() <= 2) {
            rsp.enabled = true;
            rsp.run = get_number(values[0]);
            if (values.size() == 1) {
                rsp.dir = "./";
            } else {
                rsp.dir = values[1];
                if (rsp.dir[rsp.dir.size() - 1] != '/') rsp.dir += "/";
#ifdef WIN32
                std::replace(rsp.dir.begin(), rsp.dir.end(),'/','\\');
#endif
            }

            DIR *dirp = opendir(rsp.dir.c_str());
            struct dirent *dp;
            while ((dp = readdir(dirp)) != nullptr) {
                string dname(dp->d_name);
                size_t rpos = dname.find("run");
                size_t epos = dname.find("evt");
                size_t rnpos = dname.find(".rndm");
                if (rpos == string::npos || epos == string::npos
                    || rnpos != dname.size() - 5)
                    continue;

                unsigned rstring = get_number(dname.substr(rpos + 3, epos - rpos - 3));
                if (rstring == rsp.run) {
                    string estring = dname.substr(epos + 3, rnpos - epos - 3);
                    rsp.events.push_back(atoi(estring.c_str()));
                }
            }
            closedir(dirp);
            std::sort(rsp.events.begin(), rsp.events.end());
            rsp.currentevent = -1;
        }
    }
}


void MPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
    // Check first if event should be seeded
    if (rsp.enabled) {
        G4RunManager *runManager = G4RunManager::GetRunManager();;
        ++rsp.currentevent;
        if (rsp.currentevent < int(rsp.events.size())) {
            std::ostringstream os;
            os << "run" << rsp.run << "evt" << rsp.events[rsp.currentevent]
               << ".rndm" << '\0';
            G4String fileOut = rsp.dir + os.str();
            // Use restoreEngineStatus since RestoreRandomNumberStatus is
            // too verbose
            //	  runManager->RestoreRandomNumberStatus (fileOut);
            HepRandom::restoreEngineStatus(fileOut);
        } else {
            runManager->AbortRun();
            cout << " No more events to rerun." << endl;
            return;
        }
    }

    double raster_r = sqrt(G4UniformRand());
    double raster_phi = 2.0 * pi * G4UniformRand();
    double rasterx = rvdx * raster_r * cos(raster_phi);
    double rastery = rvdy * raster_r * sin(raster_phi);

    double beamSpot_phi = 2.0 * pi * G4UniformRand();
    double beamSpot_ox = G4RandGauss::shoot(0, bsdx) * cos(beamSpot_phi);
    double beamSpot_oy = G4RandGauss::shoot(0, bsdy) * cos(beamSpot_phi);

    double beamSpot_x = bssx + beamSpot_ox * cos(bsphi) - beamSpot_oy * sin(bsphi);
    double beamSpot_y = bssy + beamSpot_ox * sin(bsphi) + beamSpot_oy * cos(bsphi);

    double displaceZ = displaceZs + (2.0 * G4UniformRand() - 1) * displaceZd;

    // internal generator. Particle defined by command line
    if (input_gen == "gemc_internal") {
        if (cosmics == "no") { // internal, no cosmic
            // redefining particle if in graphic mode
            if (gemcOpt->optMap["USE_GUI"].arg > 0) {
                setBeam();
            }
            // Primary Particle
            // if one uses the /gun/particle /gun/ion command then we
            // make sure not to overwrite those values
            if (Particle->GetParticleName() != "GenericIon") {
                particleGun->SetParticleDefinition(Particle);
            }

            // this has to go in GeneratePrimaries
            // GetIonTable cannot be called in the constructor of PrimaryGeneratorAction
            if (gemcOpt->optMap["ION_P"].args != "no") {

                vector <string> values = get_info(gemcOpt->optMap["ION_P"].args, string(",\""));

                if (values.size() > 1) {

                    G4int Z = get_number(values[0]);
                    G4int A = get_number(values[1]);

                    G4int J = 0;

                    Particle = G4IonTable::GetIonTable()->GetIon(Z, A, J);
                    particleGun->SetParticleDefinition(Particle);

                    if (values.size() == 3) {
                        G4int Q = get_number(values[2]);
                        particleGun->SetParticleCharge(Q);
                    }
                }
            }


            // 4-momenta
            double Mom = mom / MeV + (2.0 * G4UniformRand() - 1.0) * dmom / MeV;
            double Theta = acos(G4UniformRand() * (cos(theta / rad - dtheta / rad) - cos(theta / rad + dtheta / rad)) + cos(theta / rad + dtheta / rad)) / rad;
            if (primaryFlat) {
                Theta = theta / rad + (2.0 * G4UniformRand() - 1.0) * dtheta / rad;
            }
            double Phi = phi / rad + (2.0 * G4UniformRand() - 1.0) * dphi / rad;
            double mass = Particle->GetPDGMass();


            // recalculating momentun, if isKINE its value is kinetic energy
            if (isKINE) {
                double kk = Mom;
                double p2 = kk * kk + 2 * mass * kk;
                if (p2 >= 0) Mom = sqrt(p2);
                else Mom = sqrt(-p2);
            }
            double akine = sqrt(Mom * Mom + mass * mass) - mass;

            // beam direction
            G4ThreeVector beam_dir;
            if (gemcOpt->optMap["ALIGN_ZAXIS"].args == "no") {
                beam_dir = G4ThreeVector(cos(Phi / rad) * sin(Theta / rad), sin(Phi / rad) * sin(Theta / rad), cos(Theta / rad));
            } else if (gemcOpt->optMap["ALIGN_ZAXIS"].args == "beamp") {
                beam_dir = G4ThreeVector(cos(phi / rad) * sin(theta / rad), sin(phi / rad) * sin(theta / rad), cos(theta / rad));
                const G4ThreeVector beam_axis(cos(phi / rad) * sin(theta / rad), sin(phi / rad) * sin(theta / rad), cos(theta / rad));
                const G4ThreeVector rotx1(cos(phi / rad), -sin(phi / rad), 0);
                beam_dir.rotate((2.0 * G4UniformRand() - 1.0) * dtheta / rad, rotx1);
                beam_dir.rotate((2.0 * G4UniformRand() - 1.0) * dphi / rad, beam_axis);
            } else {
                beam_dir = G4ThreeVector(cos(cphi / rad) * sin(ctheta / rad), sin(cphi / rad) * sin(ctheta / rad), cos(ctheta / rad));
                //		const G4ThreeVector beam_axis(cos(cphi/rad)*sin(ctheta/rad), sin(cphi/rad)*sin(ctheta/rad), cos(ctheta/rad));
                const G4ThreeVector beam_axis(beam_dir);
                const G4ThreeVector rotx1(cos(cphi / rad), -sin(cphi / rad), 0);
                beam_dir.rotate(Theta, rotx1);
                beam_dir.rotate(Phi, beam_axis);
            }

            particleGun->SetParticleEnergy(akine);
            particleGun->SetParticleMomentumDirection(beam_dir);
            G4ThreeVector beam_vrt;

            // vertex
            if (drdzOrdxdydz == 0) {

                double VR;

                if (gaussOrFlatV == 1) {
                    // r is gaussian
                    VR = G4RandGauss::shoot(sqrt(vx * vx + vy * vy), dvr);
                } else {
                    // r2 is flat
                    VR = sqrt(G4UniformRand()) * dvr;
                }

                double PHI = 2.0 * pi * G4UniformRand();
                double Vx = vx + VR * cos(PHI) + rasterx + beamSpot_x;
                double Vy = vy + VR * sin(PHI) + rastery + beamSpot_y;
                double Vz = vz + (2.0 * G4UniformRand() - 1.0) * dvz + displaceZ;

                if (resetRasterVertex || resetBeamSpot) {
                    Vx = rasterx + beamSpot_x;
                    Vy = rastery + beamSpot_y;
                }
                if (resetLundZvertex) {
                    Vz = displaceZ;
                }

                beam_vrt = G4ThreeVector(Vx, Vy, Vz);

            } else {

                double Vx, Vy, Vz;

                if (gaussOrFlatV == 1) {
                    Vx = G4RandGauss::shoot(vx, dvx) + rasterx + beamSpot_x;
                    Vy = G4RandGauss::shoot(vy, dvy) + rastery + beamSpot_y;
                    Vz = G4RandGauss::shoot(vz, dvz) + displaceZ;
                } else {
                    Vx = vx + (2.0 * G4UniformRand() - 1.0) * dvx + rasterx + beamSpot_x;
                    Vy = vy + (2.0 * G4UniformRand() - 1.0) * dvy + rastery + beamSpot_y;
                    Vz = vz + (2.0 * G4UniformRand() - 1.0) * dvz + displaceZ;
                }
                if (resetRasterVertex || resetBeamSpot) {
                    Vx = rasterx + beamSpot_x;
                    Vy = rastery + beamSpot_y;
                }
                if (resetLundZvertex) {
                    Vz = displaceZ;
                }

                beam_vrt = G4ThreeVector(Vx, Vy, Vz);

            }
            particleGun->SetParticlePosition(beam_vrt);


            // polarization
            double partPol = 0.0;
            double polCast = 100.0 * G4UniformRand();
            if (polCast <= polDeg) partPol = 1;
            double polX = partPol * sin(polTheta / rad) * cos(polPhi / rad);
            double polY = partPol * sin(polTheta / rad) * sin(polPhi / rad);
            double polZ = partPol * cos(polTheta / rad);
            particleGun->SetParticlePolarization(G4ThreeVector(polX, polY, polZ));

            // Primary particle generated at the chosen start time
            particleGun->SetParticleTime(TSIGNAL);
            particleGun->SetNumberOfParticles(1);
            particleGun->GeneratePrimaryVertex(anEvent);
            if (GEN_VERBOSITY > 3) {
                cout << hd_msg << " Particle id=" << Particle->GetParticleName()
                     << "  Vertex=" << beam_vrt / cm << "cm,  momentum=" << Mom / GeV << " GeV, Kinematic Energy=" << akine / GeV << " GeV, theta="
                     << Theta / deg << " degrees,   phi=" << Phi / deg << " degrees" << endl;
                if (partPol > 0)
                    cout << hd_msg << "   with polarization  angles: polar - " << polTheta / deg << " degrees, "
                         << "azimuthal - " << polPhi / deg << " degrees ";
                cout << endl;
            }
        } else { // internal, cosmic model

            // paper: A. Dar, Phys.Rev.Lett, 51,3,p.227 (1983)
            bool cosmicNeutrons = false;
            if (cosmicParticle != "muon") cosmicNeutrons = true;
            
            double thisMom;
            double thisKinE,KinEmin, KinEmax;
            double thisthe;
            double thisPhi;
            double akine;

            // first randomly pick a number inside the sphere
            double cosmicVX = 100000;
            double cosmicVY = 100000;
            double cosmicVZ = 100000;

            while ((cosmicVX - cosmicTarget.x()) * (cosmicVX - cosmicTarget.x()) +
                   (cosmicVY - cosmicTarget.y()) * (cosmicVY - cosmicTarget.y()) +
                   (cosmicVZ - cosmicTarget.z()) * (cosmicVZ - cosmicTarget.z()) >= cosmicRadius * cosmicRadius) {
                
                // point generated inside spherical or cylindrical volume
                if (cosmicGeo == "sph" || cosmicGeo == "sphere") {
                    // point inside spherical volume
                    cosmicVX = cosmicTarget.x() - cosmicRadius + 2 * cosmicRadius * G4UniformRand();
                    cosmicVY = cosmicTarget.y() - cosmicRadius + 2 * cosmicRadius * G4UniformRand();
                    cosmicVZ = cosmicTarget.z() - cosmicRadius + 2 * cosmicRadius * G4UniformRand();
                } else if(cosmicGeo == "rec"){ //M.S.
                    // Uniform over a rectangular surface Z = 6*Radius, X = 2*Radius, Y = 2*Radius
                    cosmicVX = (-cosmicRadius + 2 * cosmicRadius * G4UniformRand());
                    cosmicVZ = 3 * (-cosmicRadius + 2 * cosmicRadius * G4UniformRand());
                    cosmicVY = -cosmicRadius + 2 * cosmicRadius * G4UniformRand();
                } else {
                    // point inside a cylinder, height of the cylinder = cosmicRadius/2.
                    double h = cosmicRadius / 2.;
                    cosmicVX = -cosmicRadius + 2 * cosmicRadius * G4UniformRand();
                    double sig = 1.;
                    if ((2. * G4UniformRand() - 1) < 0) sig = -1.;
                    cosmicVY = h * (2. * G4UniformRand() - 1) * sig;
                    cosmicVZ = -cosmicRadius + 2 * cosmicRadius * G4UniformRand();
                }
            }

            if (cosmicNeutrons) {
                if (cminp < 0.1 || cmaxp > 10000) cout << "WARNING !!!! COSMIC NEUTRONS E (MeV) is OUT OF THE VALID RANGE !!!" << endl;
                // Model by Ashton (1973)
                // now generating random momentum, cos(theta)
                // the maximum of the distribution is the lowest momentum and 0 theta
                // normalizing by that number
                double cosmicProb = G4UniformRand() * cosmicNeutBeam(0., cminp / GeV);
                thisMom = cminp + (cmaxp - cminp) * G4UniformRand(); // momentum in MeV/c
                thisthe = pi * G4UniformRand() / 2.0; // [0,pi/2] zenith angle
                while (cosmicNeutBeam(thisthe, thisMom / GeV) < cosmicProb) {
                    thisMom = cminp + (cmaxp - cminp) * G4UniformRand();
                    thisthe = pi * G4UniformRand() / 2.0;
                }
            } else {
                // muons
                // now generating random momentum, cos(theta)
                // the maximum of the distribution is the lowest momentum and 0 theta
                // normalizing by that number
                
                KinEmin = sqrt(cminp / GeV * cminp / GeV + 0.104 * 0.104) - 0.104;
                KinEmax = sqrt(cmaxp / GeV * cmaxp / GeV + 0.104 * 0.104) - 0.104;
                
                double cosmicProbMax = cosmicMuBeam(1., KinEmin);              //max prob for theta vertical and low momentum
                double cosmicProbMin = cosmicMuBeam(0., KinEmax);              //min prob for theta horizontal and high momentum
                double cosmicProb = G4UniformRand() * (cosmicProbMax - cosmicProbMin) + cosmicProbMin;

                double thiscthe;
                
                thiscthe = G4UniformRand();              //ctheta uniform between 0 and 1
                thisthe = acos(thiscthe);

                thisKinE = (KinEmax - KinEmin) * G4UniformRand() + KinEmin;
                thisMom = sqrt((thisKinE + 0.104) * (thisKinE + 0.104) - 0.104 * 0.104);
                
                double cosmic = cosmicMuBeam(thiscthe, thisKinE);
                
                int Nextr = 0;
                while (cosmic < cosmicProb && Nextr < 1000000) {
                    thisMom = 0;
                    Nextr = Nextr + 1;

                    thisKinE = (KinEmax - KinEmin) * G4UniformRand() + KinEmin;
                    thisMom = sqrt((thisKinE + 0.104) * (thisKinE + 0.104) - 0.104 * 0.104);

                    thiscthe = G4UniformRand();
                    thisthe = acos(thiscthe);

                    cosmic = cosmicMuBeam(thiscthe, thisKinE);
                    cosmicProb = G4UniformRand() * (cosmicProbMax - cosmicProbMin) + cosmicProbMin;
                }
                if (Nextr > 999999) cout << " !!! LOOPING IN MU EXTRACTION !!! exceeded " << Nextr << " extractions !!!" << Nextr << endl;
                
                thisMom = thisMom * GeV;
            }

            // isotropic in phi
            thisPhi = -pi + 2 * pi * G4UniformRand();

            // now finding the vertex. Assuming twice the radius as starting point
            // attention:
            // axis transformation, z <> y,  x <> -x
            // cause the cosmics come from the sky
            double pvx = cosmicVX - 2 * cosmicRadius * sin(thisthe) * cos(thisPhi);
            double pvy = cosmicVY + 2 * cosmicRadius * cos(thisthe);
            double pvz = cosmicVZ + 2 * cosmicRadius * sin(thisthe) * sin(thisPhi);
            particleGun->SetParticlePosition(G4ThreeVector(pvx, pvy, pvz));

            if (cosmicNeutrons) {
                Particle = particleTable->FindParticle("neutron");
            } else {
                // choosing charge of the muons
                string muonType = "mu+";
                if (G4UniformRand() <= 0.5)
                    muonType = "mu-";
                Particle = particleTable->FindParticle(muonType);
            }
            double mass = Particle->GetPDGMass();
            akine = sqrt(thisMom * thisMom + mass * mass) - mass;

            particleGun->SetParticleDefinition(Particle);

            // when assigning momentum the direction is reversed
            G4ThreeVector beam_dir(cos(thisPhi) * sin(thisthe), -cos(thisthe), -sin(thisPhi) * sin(thisthe));

            if (GEN_VERBOSITY > 3) {
                cout << hd_msg << " Particle id=" << Particle->GetParticleName()
                     << "  Vertex=" << G4ThreeVector(pvx, pvy, pvz) / cm << "cm,  momentum=" << thisMom / GeV << " GeV, theta="
                     << thisthe / deg << " degrees,   phi=" << thisPhi / deg << " degrees" << endl;
                cout << endl;
            }

            particleGun->SetParticleEnergy(akine);
            particleGun->SetParticleMomentumDirection(beam_dir);
            particleGun->SetNumberOfParticles(1);
            particleGun->GeneratePrimaryVertex(anEvent);
        }
    } else { // external generator: input file


        // TODO: check beagle getLine is more elegant

        // LUND format:
        // Header (Event Info):
        // These are the original LUND variables, however after # particles, and except beam polarization, these can be user defined.
        // 1               2                     3                    4               5                 6  7  8   9   10
        // # of Particles, # of Target Nucleons, # of Target Protons, Pol. of Target, Pol. of Electron, x, y, W2, Q2, nu
        //
        // Body (Particle Info):
        // 1       2      3     4            5       6         7    8    9    10   11    12        13        14
        // index, charge, type, particle id, parent, daughter, p_x, p_y, p_z, p_t, mass, x vertex, y vertex, z vertex
        // type is 1 for particles in the detector
        if ((gformat == "LUND" || gformat == "lund") && !gif.eof()) {

            headerUserDefined.clear();

            string theWholeLine;

            if (rsp.enabled && eventIndex < int(rsp.events[rsp.currentevent])) {
                // Skip input file lines to find rerun event
                for (; eventIndex < int(rsp.events[rsp.currentevent]); ++eventIndex) {
                    // reading header
                    getline(gif, theWholeLine);
                    if (gif.eof()) {
                        return;
                    }
                    vector <string> headerStrings = getStringVectorFromString(theWholeLine);
                    headerUserDefined.clear();
                    for (auto &s: headerStrings) {
                        headerUserDefined.push_back(get_number(s));
                    }

                    int nparticles = headerUserDefined[0];
                    for (int p = 0; p < nparticles; p++) {
                        string theWholeLine;
                        getline(gif, theWholeLine);
                        if (gif.eof() && p != nparticles - 1) {
                            cout << " Input file " << gfilename << " appear to be truncated." << endl;
                            return;
                        }
                    }
                }
            }
            // reading header
            getline(gif, theWholeLine);
            if (gif.eof()) {
                return;
            }
            vector <string> headerStrings = getStringVectorFromString(theWholeLine);
            headerUserDefined.clear();
            for (auto &s: headerStrings) {
                headerUserDefined.push_back(get_number(s));
            }

            int nparticles = headerUserDefined[0];
            beamPol = headerUserDefined[4];
            if (beamPol > 1) {
                beamPol = 1;
            }

            userInfo.clear();
            for (int p = 0; p < nparticles; p++) {
                string theWholeLine;
                getline(gif, theWholeLine);
                if (gif.eof() && p != nparticles - 1) {
                    cout << " Input file " << gfilename << " appear to be truncated." << endl;
                    return;
                }

                vector <string> infoStrings = getStringVectorFromString(theWholeLine);

                userInforForParticle thisParticleInfo;
                for (auto &s: infoStrings) {
                    thisParticleInfo.infos.push_back(get_number(s));
                }
                userInfo.push_back(thisParticleInfo);

                // necessary geant4 info. Lund specifics:
                int pindex = thisParticleInfo.infos[0];
                int type = thisParticleInfo.infos[2];
                int pdef = thisParticleInfo.infos[3];
                double px = thisParticleInfo.infos[6];
                double py = thisParticleInfo.infos[7];
                double pz = thisParticleInfo.infos[8];

                // vertex is  received in cm from LUND
                double Vx = thisParticleInfo.infos[11] * cm + svx + rasterx + beamSpot_x;
                double Vy = thisParticleInfo.infos[12] * cm + svy + rastery + beamSpot_y;
                double Vz = thisParticleInfo.infos[13] * cm + svz + displaceZ;
                if (resetRasterVertex || resetBeamSpot) {
                    Vx = rasterx + beamSpot_x;
                    Vy = rastery + beamSpot_y;
                }
                if (resetLundZvertex) {
                    Vz = displaceZ;
                }


                bool isNan =
                        isnan(p) && isnan(pindex) && isnan(type) && isnan(pdef) && isnan(px) && isnan(py) && isnan(pz) && isnan(Vx) && isnan(Vy) && isnan(Vz);

                if (isNan) {
                    cout << " Error: NAN detected: " << endl;
                    cout << "p: " << isnan(p) << endl;
                    cout << "pindex: " << isnan(pindex) << endl;
                    cout << "type: " << isnan(type) << endl;
                    cout << "pdef: " << isnan(pdef) << endl;
                    cout << "px: " << isnan(px) << endl;
                    cout << "py: " << isnan(py) << endl;
                    cout << "pz: " << isnan(pz) << endl;
                    cout << "Vx: " << isnan(Vx) << endl;
                    cout << "Vy: " << isnan(Vy) << endl;
                    cout << "Vz: " << isnan(Vz) << endl;
                }

                if (PROPAGATE_DVERTEXTIME == 0) {
                    setParticleFromPars(p, pindex, type, pdef, px, py, pz, Vx, Vy, Vz, anEvent);
                } else {
                    // if this flag is set updated times are calculated for detached vertex events
                    setParticleFromParsPropagateTime(p, userInfo, anEvent);
                }
            }

            if (eventIndex <= ntoskip) {
                if (GEN_VERBOSITY > 3) {
                    cout << " This event will be skipped." << endl;
                }
            }
            eventIndex++;
        } else if ((gformat == "BEAGLE" || gformat == "beagle") && !gif.eof()) {
            // Format:
            // https://wiki.bnl.gov/eic/index.php/BeAGLE#Output_Data_Format
            //
            // BEAGLE EVENT FILE
            // ============================================
            // I, ievent, genevent, lepton, Atarg, Ztarg, pzlep, pztarg  pznucl, crang, crori, subprocess, nucleon,
            // targetparton,          xtargparton, beamparton, xbeamparton, thetabeamprtn, truey,       trueQ2, truex,
            // trueW2, trueNu, leptonphi, s_hat, t_hat, u_hat,    pt2_hat, Q2_hat, F2, F1, R, sigma_rad, SigRadCor, EBrems,
            // photonflux, b, Phib, Thickness, ThickScl, Ncollt, Ncolli,         Nwound, Nwdch, Nnevap, Npevap, Aremn,
            // NINC, NINCch, d1st, davg,   pxf, pyf, pzf, Eexc, RAevt, User1, User2, User3, nrTracks
            // ============================================
            // I  ISTHKK(I)  IDHKK(I)  JMOHKK(2,I)  JMOHKK(1,I)       JDAHKK(1,I)  JDAHKK(2,I)  PHKK(1,I)  PHKK(2,I)  PHKK(3,I)
            // PHKK(4,I)  PHKK(5,I)  VHKK(1,I) VHKK(2,I) VHKK(3,I) IDRES(I)      IDXRES(I) NOBAM(I)
            // ============================================

            string theWholeLine;

            // first 6 lines are file header
            if (beagleHeader == 0) {
                getline(gif, theWholeLine);
                getline(gif, theWholeLine);
                getline(gif, theWholeLine);
                getline(gif, theWholeLine);
                getline(gif, theWholeLine);
                getline(gif, theWholeLine);
                beagleHeader = 1;
            }

            // reading header
            getline(gif, theWholeLine);
            vector <string> headerStrings = getStringVectorFromString(theWholeLine);
            headerUserDefined.clear();
            for (auto &s: headerStrings) {
                headerUserDefined.push_back(get_number(s));
            }

            int nparticles = headerUserDefined.back();

            // header / particles separator
            getline(gif, theWholeLine);


            // getting info for each particle
            userInfo.clear();
            for (int p = 0; p < nparticles; p++) {

                getline(gif, theWholeLine);

                // reaching eof prematurely
                if (gif.eof()) {
                    cout << " Input file " << gfilename << " appear to be truncated." << endl;
                    return;
                }

                vector <string> infoStrings = getStringVectorFromString(theWholeLine);

                userInforForParticle thisParticleInfo;
                for (auto &s: infoStrings) {
                    thisParticleInfo.infos.push_back(get_number(s));
                }
                if (thisParticleInfo.infos.size() != 18) {
                    cout << " !!! Error: Beagle particle info size is " << thisParticleInfo.infos.size() << " instead of 18." << endl;
                }
                userInfo.push_back(thisParticleInfo);

                // necessary geant4 info. Lund specifics:
                int pindex = thisParticleInfo.infos[0];
                int type = thisParticleInfo.infos[1];
                int pdef = thisParticleInfo.infos[2];
                double px = thisParticleInfo.infos[7];
                double py = thisParticleInfo.infos[8];
                double pz = thisParticleInfo.infos[9];

                double A = thisParticleInfo.infos[15];
                double Z = thisParticleInfo.infos[16];

                // vertex is  received in cm from LUND
                double Vx = thisParticleInfo.infos[12] * cm + svx + rasterx + beamSpot_x;
                double Vy = thisParticleInfo.infos[13] * cm + svy + rastery + beamSpot_y;
                double Vz = thisParticleInfo.infos[14] * cm + svz;
                if (resetRasterVertex || resetBeamSpot) {
                    Vx = rasterx + beamSpot_x;
                    Vy = rastery + beamSpot_y;
                }
                if (resetLundZvertex) {
                    Vz = displaceZ;
                }

                setParticleFromPars(p, pindex, type, pdef, px, py, pz, Vx, Vy, Vz, anEvent, A, Z);
            }

            if (eventIndex <= ntoskip) {
                if (GEN_VERBOSITY > 3) {
                    cout << " This event will be skipped." << endl;
                }
                eventIndex++;
            }


            // reading end of event separator
            getline(gif, theWholeLine);
        } else if (gformat == "stdhep" || gformat == "STDHEP" || gformat == "StdHep" || gformat == "StdHEP") {
            //
            // StdHep is an (old like LUND) MC generator format in binary form.
            //

            long lerr = stdhep_reader->readEvent();  // Read the next event from the file.

            if (lerr == LSH_ENDOFFILE) {
                return;
            } else if (lerr != LSH_SUCCESS) {
                cout << hd_msg << " -- Error reading stdhep file: " << lerr << "\n";
                return;
            }

            int NPART = stdhep_reader->nTracks();

            for (int p = 0; p < NPART; p++) {
                if (stdhep_reader->daughter1(p) > 0) {
                    // The particle has daughters, so we do not want to generate this one.
                    continue;
                } else {

                    Particle = particleTable->FindParticle(stdhep_reader->pid(p));
                    if (!Particle) {
                        cout << hd_msg << " Particle id " << stdhep_reader->pid(p) << " not found in G4 table." << endl << endl;

                        cout << "Exiting !" << endl;
                        exit(1);
                    }

                    particleGun->SetParticleDefinition(Particle);

                    // 4-momenta
                    G4ThreeVector pmom(stdhep_reader->Px(p) * GeV, stdhep_reader->Py(p) * GeV, stdhep_reader->Pz(p) * GeV);
                    double Mom = pmom.mag();
                    double Phi = pmom.getPhi();
                    double Theta = pmom.getTheta();
                    double mass = Particle->GetPDGMass();
                    double akine = sqrt(Mom * Mom + mass * mass) - mass;

                    G4ThreeVector beam_dir(cos(Phi / rad) * sin(Theta / rad), sin(Phi / rad) * sin(Theta / rad), cos(Theta / rad));

                    if (gemcOpt->optMap["STEER_BEAM"].arg != 0) {
                        beam_dir.rotateY(theta);
                        beam_dir.rotateZ(phi);
                    }

                    particleGun->SetParticleEnergy(akine);
                    particleGun->SetParticleMomentumDirection(beam_dir);

                    G4ThreeVector beam_vrt;

                    // vertex is  received in cm from stdhep_reader
                    if (gemcOpt->optMap["STEER_BEAM"].arg == 0) {
                        beam_vrt = G4ThreeVector(stdhep_reader->X(p) * cm,
                                                 stdhep_reader->Y(p) * cm,
                                                 stdhep_reader->Z(p) * cm);
                    } else {
                        // vertex smear and offset
                        double VR = sqrt(G4UniformRand()) * dvr;
                        double PHI = 2.0 * pi * G4UniformRand();

                        // vertex is  received in cm from stdhep_reader
                        beam_vrt = G4ThreeVector(stdhep_reader->X(p) * cm + vx + VR * cos(PHI),
                                                 stdhep_reader->Y(p) * cm + vy + VR * sin(PHI),
                                                 stdhep_reader->Z(p) * cm + vz + (2.0 * G4UniformRand() - 1.0) * dvz);

                    }
                    particleGun->SetParticlePosition(beam_vrt);


                    // beam polarization only along the beam
                    // only for the first particle
                    if (p == 0) {
                        particleGun->SetParticlePolarization(G4ThreeVector(0, 0, beamPol));
                    }

                    // Primary particle generated at the chosen start time
                    particleGun->SetParticleTime(TSIGNAL);
                    particleGun->GeneratePrimaryVertex(anEvent);
                    if (GEN_VERBOSITY > 3)
                        cout << hd_msg << " Particle Number:  " << p << ", id=" << stdhep_reader->pid(p) << "(" << Particle->GetParticleName() << ")"
                             << "  Vertex=" << beam_vrt / cm << "cm,  momentum=" << pmom / GeV << " GeV" << endl;
                }
            }
        }
    }

    // merging (background) events from LUND format
    if (background_gen != "no") {
        int nparticles;
        double tmp;

        bgif >> nparticles;
        for (unsigned i = 0; i < 9; i++)
            bgif >> tmp;

        for (int p = 0; p < nparticles; p++) {
            double px, py, pz;
            int pdef, pindex;
            double time;
            double Vx, Vy, Vz;
            bgif >> pindex >> tmp >> tmp >> pdef >> tmp >> tmp >> px >> py >> pz >> time >> tmp >> Vx >> Vy >> Vz;
            if (pindex == p + 1) {
                // Primary Particle
                Particle = particleTable->FindParticle(pdef);
                if (!Particle) {
                    cout << hd_msg << " Particle id " << pdef << " not found in G4 table." << endl << endl;

                    return;
                }
                particleGun->SetParticleDefinition(Particle);

                // 4-momenta
                G4ThreeVector pmom(px * GeV, py * GeV, pz * GeV);
                double Mom = pmom.mag();
                double Phi = pmom.getPhi();
                double Theta = pmom.getTheta();
                double mass = Particle->GetPDGMass();
                double akine = sqrt(Mom * Mom + mass * mass) - mass;

                particleGun->SetParticleEnergy(akine);
                particleGun->SetParticleMomentumDirection(
                        G4ThreeVector(cos(Phi / rad) * sin(Theta / rad), sin(Phi / rad) * sin(Theta / rad), cos(Theta / rad)));

                // vertex is  received in cm from LUND
                G4ThreeVector beam_vrt(Vx * cm, Vy * cm, Vz * cm);
                particleGun->SetParticlePosition(beam_vrt);


                // Primary particle generated int the middle of Time window
                particleGun->SetParticleTime(time);
                particleGun->GeneratePrimaryVertex(anEvent);
                if (GEN_VERBOSITY > 3)
                    cout << hd_msg << " Merged Particle Number:  " << p + 1 << ", id=" << pdef << " (" << Particle->GetParticleName() << ")"
                         << "  Vertex=" << beam_vrt / cm << "cm,  momentum=" << pmom / GeV << " GeV" << endl;
            } else if (pindex != p + 1)
                if (GEN_VERBOSITY > 3)
                    cout << hd_msg << " Warning: file particle index " << tmp << " does not match read particle index " << p + 1 << endl;
        }
    }

    // Luminosity Particles
    int NBUNCHES = (int) floor(TWINDOW / TBUNCH);
    //int PBUNCH     = (int) floor((double)NP/NBUNCHES) - 1;
    vector<int> BUNCHES = setBunches(NP, NBUNCHES);

    if (NP > 0) {
        particleGun->SetParticleDefinition(L_Particle);

        // getting kinematics
        double L_mass = L_Particle->GetPDGMass();
        double L_Mom = L_mom;
        double L_Theta = L_theta;
        double L_Phi = L_phi;

        // all particles in a bunch are identical
        for (int b = 0; b < NBUNCHES; b++) {
            // spread momentum if requested
            if (L_dmom > 0) {
                L_Mom = L_mom + (2.0 * G4UniformRand() - 1.0) * L_dmom;
                L_Theta = acos(G4UniformRand() * (cos(L_theta / rad - L_dtheta / rad) - cos(L_theta / rad + L_dtheta / rad)) +
                               cos(L_theta / rad + L_dtheta / rad)) / rad;
                if (lumiFlat) {
                    L_Theta = L_theta + (2.0 * G4UniformRand() - 1.0) * L_dtheta;
                }

                L_Phi = L_phi + (2.0 * G4UniformRand() - 1.0) * L_dphi;
            }
            double L_akine = sqrt(L_Mom * L_Mom + L_mass * L_mass) - L_mass;
            particleGun->SetParticleEnergy(L_akine);
            particleGun->SetParticleMomentumDirection(
                    G4ThreeVector(cos(L_Phi / rad) * sin(L_Theta / rad), sin(L_Phi / rad) * sin(L_Theta / rad), cos(L_Theta / rad)));

            // luminosity vertex
            // vertex has uniform density across the cilinder
            double lvx = L_vx;
            double lvy = L_vy;

            if (L_dvr > 0.) {
                double tmp_sqrt_rho = sqrt(G4UniformRand());  // Square root gives uniform spread over circle.
                double tmp_phi = 2 * pi * G4UniformRand();

                // This *should* have an l_dvr_x and l_dvr_y for ellipsoidal beams!!!
                lvx = L_vx + L_dvr * tmp_sqrt_rho * cos(tmp_phi);
                lvy = L_vy + L_dvr * tmp_sqrt_rho * sin(tmp_phi);
            }

            double lvz = L_vz;

            // spread vertex if requested
            if (L_dvz > 0) {
                lvz = L_vz + (2.0 * G4UniformRand() - 1.0) * L_dvz;
            }

            particleGun->SetNumberOfParticles(BUNCHES[b]);

            // cout << " bunch " << b << " " << PBUNCH << endl;

            particleGun->SetParticleTime(TBUNCH * b);
            particleGun->SetParticlePosition(G4ThreeVector(lvx, lvy, lvz));
            particleGun->GeneratePrimaryVertex(anEvent);
        }
    }

    // Luminosity Particles2
    int NBUNCHES2 = (int) floor(TWINDOW / TBUNCH2);
    vector<int> BUNCHES2 = setBunches(NP2, NBUNCHES2);


    if (NP2 > 0) {
        particleGun->SetParticleDefinition(L2_Particle);
        //		particleGun->SetNumberOfParticles(PBUNCH);

        // getting kinematics
        double L2_mass = L2_Particle->GetPDGMass();
        double L2_Mom = L2_mom;
        double L2_Theta = L2_theta;
        double L2_Phi = L2_phi;


        // all particles in a bunch are identical
        for (int b = 0; b < NBUNCHES2; b++) {
            particleGun->SetNumberOfParticles(BUNCHES2[b]);
            particleGun->SetParticleTime(TBUNCH2 * b);
            // spread momentum if requested
            if (L2_dmom > 0) {
                L2_Mom += (2.0 * G4UniformRand() - 1.0) * L2_dmom;
                L2_Theta = acos(G4UniformRand() * (cos(L2_theta / rad - L2_dtheta / rad) - cos(L2_theta / rad + L2_dtheta / rad)) +
                                cos(L2_theta / rad + L2_dtheta / rad)) / rad;
                if (lumi2Flat)
                    L2_Theta += (2.0 * G4UniformRand() - 1.0) * L2_dtheta;
                L2_Phi += (2.0 * G4UniformRand() - 1.0) * L2_dphi;
            }
            double L2_akine = sqrt(L2_Mom * L2_Mom + L2_mass * L2_mass) - L2_mass;
            particleGun->SetParticleEnergy(L2_akine);
            particleGun->SetParticleMomentumDirection(
                    G4ThreeVector(cos(L2_Phi / rad) * sin(L2_Theta / rad), sin(L2_Phi / rad) * sin(L2_Theta / rad), cos(L2_Theta / rad)));

            // luminosity vertex 2
            double L2_VR = sqrt(G4UniformRand()) * L2_dvr;
            double L2_PHI = 2.0 * pi * G4UniformRand();
            L2_vx += L2_VR * cos(L2_PHI);
            L2_vy += L2_VR * sin(L2_PHI);

            // spread vertex if requested
            if (L2_dvz > 0) {
                L2_vz += (2.0 * G4UniformRand() - 1.0) * L2_dvz;
            }
            particleGun->SetParticlePosition(G4ThreeVector(L2_vx, L2_vy, L2_vz));

            particleGun->GeneratePrimaryVertex(anEvent);
        }

    }

    if (GEN_VERBOSITY > 5)
        cout << " Generation done " << endl;
}


void MPrimaryGeneratorAction::setBeam() {
    string hd_msg = gemcOpt->optMap["LOG_MSG"].args + " Beam Settings >> ";
    svx = 0;
    svy = 0;
    svz = 0;


    // vector of string - filled from the various option
    vector <string> values;
    string units;

    // getting vertex shifts
    values = get_info(gemcOpt->optMap["SHIFT_LUND_VERTEX"].args);
    string toshifts = trimSpacesFromString(values[0]);
    if (toshifts != "no") {
        units = trimSpacesFromString(values[3]);
        svx = get_number(values[0] + "*" + units);
        svy = get_number(values[1] + "*" + units);
        svz = get_number(values[2] + "*" + units);

        if (GEN_VERBOSITY > 3) {
            cout << hd_msg << " Vertex shift (" << svx << ", " << svy << ", " << svz << ")"
                 << units << endl;
        }

    }

    // Getting raster vertex from option value
    values = get_info(gemcOpt->optMap["RASTER_VERTEX"].args);

    rvdx = get_number(values[0]);
    rvdy = get_number(values[1]);
    resetRasterVertex = false;
    if (values.back().find("reset") != string::npos) {
        resetRasterVertex = true;
    }
    bool isNanVX = isnan(rvdx) && isnan(rvdy);

    if (isNanVX) {
        cout << " Error: NAN detected for RASTER_VERTEX parameters " << endl;
        cout << "rvdx: " << isnan(rvdx) << endl;
        cout << "rvdy: " << isnan(rvdy) << endl;
    }


    // Getting beamspot vertex from option value
    values = get_info(gemcOpt->optMap["BEAM_SPOT"].args);

    bssx = get_number(values[0]);
    bssy = get_number(values[1]);
    bsdx = get_number(values[2]);
    bsdy = get_number(values[3]);
    bsphi = get_number(values[4]);
    resetBeamSpot = false;
    if (values.back().find("reset") != string::npos) {
        resetBeamSpot = true;
    }

    bool isNanBS = isnan(bssx) && isnan(bssy) && isnan(bsdx) && isnan(bsdy) && isnan(bsphi);

    if (isNanBS) {
        cout << " Error: NAN detected for BEAM_SPOT parameters " << endl;
        cout << "bssx: " << isnan(bssx) << endl;
        cout << "bssy: " << isnan(bssy) << endl;
        cout << "bsdx: " << isnan(bsdx) << endl;
        cout << "bsdy: " << isnan(bsdy) << endl;
        cout << "bsphi: " << isnan(bsphi) << endl;
    }


    // Getting vertex z displacement from option value
    values = get_info(gemcOpt->optMap["RANDOMIZE_LUND_VZ"].args);

    displaceZs = get_number(values[0]);
    displaceZd = get_number(values[1]);
    resetLundZvertex = false;
    if (values.back().find("reset") != string::npos) {
        resetLundZvertex = true;
    }
    bool isNanVZD = isnan(displaceZs) && isnan(displaceZd);

    if (isNanVZD) {
        cout << " Error: NAN detected for RANDOMIZE_LUND_VZ parameters " << endl;
        cout << "displaceZs: " << isnan(displaceZs) << endl;
        cout << "displaceZd: " << isnan(displaceZd) << endl;
    }


    if (input_gen == "gemc_internal") {

        if (cosmics == "no") {
            // Getting particle name,  momentum from option value
            values = get_info(gemcOpt->optMap["BEAM_P"].args, string(",\""));
            string pname = trimSpacesFromString(values[0]);
            isKINE = false;
            if (values.size() == 5) {
                if (values[4].find("KE") != string::npos) {
                    isKINE = true;
                }
            }
            if (values.size() == 4 || (values.size() == 5 && isKINE)) {
                mom = get_number(values[1]);
                theta = get_number(values[2]);
                phi = get_number(values[3]);
            }

            // making sure the particle exists
            Particle = particleTable->FindParticle(pname);
            if (!Particle) {
                // it may be the "show_all" option. In this case print all available particle names
                if (pname == "show_all") {
                    for (int i = 0; i < particleTable->entries(); i++)
                        cout << hd_msg << " g4 particle: " << particleTable->GetParticleName(i)
                             << " pdg encoding: " << particleTable->GetParticle(i)->GetPDGEncoding() << endl;
                } else {
                    // otherwise it's not found. Need to exit here.
                    cout << hd_msg << " Particle " << pname << " not found in G4 table. Exiting" << endl << endl;
                }
                exit(1);
            }

            // Getting custom beam direction if it's set
            values = get_info(gemcOpt->optMap["ALIGN_ZAXIS"].args);
            string align = trimSpacesFromString(values[0]);
            if (align == "custom") {
                ctheta = get_number(values[1]);
                cphi = get_number(values[2]);
            }

            // Getting momentum spread from option value
            primaryFlat = 0;
            values = get_info(gemcOpt->optMap["SPREAD_P"].args);
            dmom = get_number(values[0]);
            dtheta = get_number(values[1]);
            dphi = get_number(values[2]);
            if (values.size() > 3) {
                if (trimSpacesFromString(values[3]) == "flat") {
                    primaryFlat = 1;
                }
            }

            // Getting vertex from option value
            values = get_info(gemcOpt->optMap["BEAM_V"].args);
            units = trimSpacesFromString(values[3]);

            vx = get_number(values[0] + "*" + units);
            vy = get_number(values[1] + "*" + units);
            vz = get_number(values[2] + "*" + units);

            // Getting vertex spread from option value
            values = get_info(gemcOpt->optMap["SPREAD_V"].args);
            // distribution type is not given. defaults to "flat"
            // number of argument can be 3 (drdz) or 4 (dxdydz)
            if (values.back().find("gauss") == string::npos && values.back().find("flat") == string::npos) {
                gaussOrFlatV = 0;

                // check if it's (dr, dz) or (dx, dy, dz)
                if (values.size() == 3) {
                    drdzOrdxdydz = 0;
                    units = trimSpacesFromString(values[2]);
                    dvr = get_number(values[0] + "*" + units);
                    dvz = get_number(values[1] + "*" + units);
                } else if (values.size() == 4) {
                    drdzOrdxdydz = 1;
                    units = trimSpacesFromString(values[3]);
                    dvx = get_number(values[0] + "*" + units);
                    dvy = get_number(values[1] + "*" + units);
                    dvz = get_number(values[2] + "*" + units);
                }
            } else {
                if (values.back().find("gauss") == string::npos) gaussOrFlatV = 0;
                else gaussOrFlatV = 1;

                // check if it's (dr, dz) or (dx, dy, dz)
                if (values.size() == 4) {
                    drdzOrdxdydz = 0;
                    units = trimSpacesFromString(values[2]);
                    dvr = get_number(values[0] + "*" + units);
                    dvz = get_number(values[1] + "*" + units);
                } else if (values.size() == 5) {
                    drdzOrdxdydz = 1;
                    units = trimSpacesFromString(values[3]);
                    dvx = get_number(values[0] + "*" + units);
                    dvy = get_number(values[1] + "*" + units);
                    dvz = get_number(values[2] + "*" + units);
                }
            }


            // Getting polarization from option value
            values = get_info(gemcOpt->optMap["POLAR"].args);
            polDeg = get_number(values[0]);
            polTheta = get_number(values[1]);
            polPhi = get_number(values[2]);
        } else {
            vector <string> csettings = get_info(cosmics, string(",\""));
            int len = csettings.size();
            // parsing information for COSMIC RAYS option
            if (csettings[0] == "default") {
                cosmicA = 55.6;
                cosmicB = 1.04;
                cosmicC = 64;

                cminp = get_number(csettings[1], 0) * GeV;
                cmaxp = get_number(csettings[2], 0) * GeV;

                // model is valid only starting at 1 GeV for now
                if (cminp < 1) cminp = 1;

                // select cosmic ray particle from data card
                if (len > 3) {
                    cosmicParticle = csettings[3];
                } else {
                    cosmicParticle = "muon";
                }
            } else {
                cosmicA = get_number(csettings[0], 0);
                cosmicB = get_number(csettings[1], 0);
                cosmicC = get_number(csettings[2], 0);

                cminp = get_number(csettings[3], 0) * GeV;
                cmaxp = get_number(csettings[4], 0) * GeV;

                // model is valid only starting at 1 GeV for now
                if (cminp < 1) cminp = 1;

                // select cosmic ray particle from data card
                if (len > 5) {
                    cosmicParticle = csettings[5];
                } else {
                    cosmicParticle = "muon";
                }
            }
        }
    } else if (input_gen.compare(0, 4, "LUND") == 0 || input_gen.compare(0, 4, "lund") == 0) {
        gformat.assign(input_gen, 0, input_gen.find(","));
        gfilename.assign(input_gen, input_gen.find(",") + 1, input_gen.size());
        cout << hd_msg << "LUND: Opening  " << gformat << " file: " << trimSpacesFromString(gfilename).c_str() << endl;
        gif.open(trimSpacesFromString(gfilename).c_str());
        if (!gif) {
            cerr << hd_msg << " Can't open LUND input file " << trimSpacesFromString(gfilename).c_str() << ". Exiting. " << endl;
            exit(201);
        }
    } else if (input_gen.compare(0, 6, "BEAGLE") == 0 || input_gen.compare(0, 6, "beagle") == 0) {
        gformat.assign(input_gen, 0, input_gen.find(","));
        gfilename.assign(input_gen, input_gen.find(",") + 1, input_gen.size());
        cout << hd_msg << "BEAGLE: Opening  " << gformat << " file: " << trimSpacesFromString(gfilename).c_str() << endl;
        gif.open(trimSpacesFromString(gfilename).c_str());
        if (!gif) {
            cerr << hd_msg << " Can't open BEAGLE input file " << trimSpacesFromString(gfilename).c_str() << ". Exiting. " << endl;
            exit(202);
        }
        beagleHeader = 0;
    } else if (input_gen.compare(0, 6, "stdhep") == 0 || input_gen.compare(0, 6, "STDHEP") == 0 ||
               input_gen.compare(0, 6, "StdHep") == 0 || input_gen.compare(0, 6, "StdHEP") == 0) {
        // StdHep is an (old like LUND) MC generator format in binary form.
        gformat.assign(input_gen, 0, input_gen.find(","));
        gfilename.assign(input_gen, input_gen.find(",") + 1, input_gen.size());
        cout << hd_msg << "StdHEP: Opening  " << gformat << " file: " << trimSpacesFromString(gfilename).c_str() << endl;
        stdhep_reader = new lStdHep(trimSpacesFromString(gfilename).c_str());

        if (!stdhep_reader) {
            cerr << hd_msg << " Can't open input file " << trimSpacesFromString(gfilename).c_str() << ". Exiting. " << endl;
            exit(203);
        }

        // For the STEER_BEAM option, we need to have the angles and vertex of the GCARD in BEAM_P and BEAM_V, SPREAD_V
        // Getting particle name,  momentum from option value
        values = get_info(gemcOpt->optMap["BEAM_P"].args);
        string pname = trimSpacesFromString(values[0]);
        mom = get_number(values[1]);
        theta = get_number(values[2]);
        phi = get_number(values[3]);

        // Getting vertex from option value
        values = get_info(gemcOpt->optMap["BEAM_V"].args);
        units = trimSpacesFromString(values[3]);
        vx = get_number(values[0] + "*" + units);
        vy = get_number(values[1] + "*" + units);
        vz = get_number(values[2] + "*" + units);

        // Getting vertex spread from option value
        values = get_info(gemcOpt->optMap["SPREAD_V"].args);
        units = trimSpacesFromString(values[2]);
        dvr = get_number(values[0] + "*" + units);
        dvz = get_number(values[1] + "*" + units);

    }

    // merging (background) events from LUND format
    if (background_gen != "no") {
        // file may be already opened cause setBeam is called again in graphic mode
        if (!bgif.is_open()) {
            bgif.open(trimSpacesFromString(background_gen).c_str());
            if (!bgif) {
                cerr << hd_msg << " Can't open background input file >" << trimSpacesFromString(background_gen).c_str() << "<. Exiting. " << endl;
                exit(204);
            }
        }
    }


    // %%%%%%%%%%%%%%%
    // Luminosity Beam
    // %%%%%%%%%%%%%%%

    // Getting particle name,  momentum from option value
    values = get_info(gemcOpt->optMap["LUMI_P"].args);
    string L_pname = trimSpacesFromString(values[0]);
    L_mom = get_number(values[1]);
    L_theta = get_number(values[2]);
    L_phi = get_number(values[3]);


    // Getting momentum spread from option value
    lumiFlat = 0;
    values = get_info(gemcOpt->optMap["LUMI_SPREAD_P"].args);
    L_dmom = get_number(values[0]);
    L_dtheta = get_number(values[1]);
    L_dphi = get_number(values[2]);
    if (values.size() == 4)
        if (trimSpacesFromString(values[3]) == "flat")
            lumiFlat = 1;


    // making sure the particle exists
    L_Particle = particleTable->FindParticle(L_pname);
    if (!L_Particle) {
        // it may be the "show_all" option. In this case print all available particle names
        if (L_pname == "show_all") {
            for (int i = 0; i < particleTable->entries(); i++)
                cout << hd_msg << " g4 particle: " << particleTable->GetParticleName(i) << endl;
        }
            // otherwise it's not found. Need to exit here.
        else
            cout << hd_msg << " Particle " << L_pname << " not found in G4 table. Exiting" << endl << endl;

        exit(1);
    }

    // Getting vertex from option value
    values = get_info(gemcOpt->optMap["LUMI_V"].args);
    units = trimSpacesFromString(values[3]);
    L_vx = get_number(values[0] + "*" + units);
    L_vy = get_number(values[1] + "*" + units);
    L_vz = get_number(values[2] + "*" + units);

    // Getting vertex spread from option value
    values = get_info(gemcOpt->optMap["LUMI_SPREAD_V"].args);
    units = trimSpacesFromString(values[2]);
    L_dvr = get_number(values[0] + "*" + units);
    L_dvz = get_number(values[1] + "*" + units);

    // Getting parameters from option value
    values = get_info(gemcOpt->optMap["LUMI_EVENT"].args);
    NP = (int) get_number(values[0]);
    TWINDOW = get_number(values[1]);
    TBUNCH = get_number(values[2]);
    TSIGNAL = TWINDOW / 2;
    if (values.size() == 4) TSIGNAL = get_number(values[3]);
    TSIGNAL = TBUNCH * floor(TSIGNAL / TBUNCH); // rounding to previous bunch


    // %%%%%%%%%%%%%%%%%
    // Luminosity Beam 2
    // %%%%%%%%%%%%%%%%%

    // Getting particle name,  momentum from option value
    values = get_info(gemcOpt->optMap["LUMI2_P"].args);
    string L2_pname = trimSpacesFromString(values[0]);
    L2_mom = get_number(values[1]);
    L2_theta = get_number(values[2]);
    L2_phi = get_number(values[3]);

    // Getting momentum spread from option value
    lumi2Flat = 0;
    values = get_info(gemcOpt->optMap["LUMI2_SPREAD_P"].args);
    L2_dmom = get_number(values[0]);
    L2_dtheta = get_number(values[1]);
    L2_dphi = get_number(values[2]);
    if (values.size() == 4)
        if (trimSpacesFromString(values[3]) == "flat")
            lumi2Flat = 1;


    // making sure the particle exists
    L2_Particle = particleTable->FindParticle(L2_pname);
    if (!L2_Particle) {
        // it may be the "show_all" option. In this case print all available particle names
        if (L_pname == "show_all") {
            for (int i = 0; i < particleTable->entries(); i++)
                cout << hd_msg << " g4 particle: " << particleTable->GetParticleName(i) << endl;
        } else
            // otherwise it's not found. Need to exit here.
            cout << hd_msg << " Particle " << L2_pname << " not found in G4 table. Exiting" << endl << endl;
        exit(1);
    }

    // Getting vertex from option value
    values = get_info(gemcOpt->optMap["LUMI2_V"].args);
    units = trimSpacesFromString(values[3]);
    L2_vx = get_number(values[0] + "*" + units);
    L2_vy = get_number(values[1] + "*" + units);
    L2_vz = get_number(values[2] + "*" + units);


    // Getting vertex spread from option value
    values = get_info(gemcOpt->optMap["LUMI2_SPREAD_V"].args);
    units = trimSpacesFromString(values[2]);
    L2_dvr = get_number(values[0] + "*" + units);
    L2_dvz = get_number(values[1] + "*" + units);

    // Getting parameters from option value
    values = get_info(gemcOpt->optMap["LUMI2_EVENT"].args);
    NP2 = (int) get_number(values[0]);
    TBUNCH2 = get_number(values[1]);

}


vector<int> MPrimaryGeneratorAction::setBunches(int nparticle, int nbunch) {
    vector<int> bunches;
    if (nparticle > 0) {
        for (int i = 0; i < nbunch; i++) bunches.push_back(0);
        for (int i = 0; i < nparticle; i++) {
            int ib = (int) floor(G4UniformRand() * nbunch);
            bunches[ib]++;
        }
    }
    return bunches;
}

MPrimaryGeneratorAction::~MPrimaryGeneratorAction() {
    delete particleGun;
    gif.close();
    bgif.close();
}


double MPrimaryGeneratorAction::cosmicMuBeam(double c, double e) { //M.S.
    // cosmic muons spectrum as a function of energy (GeV) and cosine of zenith angle
    // from 1509.06176: this is dN / dE dOmega
    // note that E is the kinetic energy

        double pa[5] = { 0.102573, -0.068287, 0.958633, 0.0407253, 0.817285 };
        double cst = sqrt((c * c + pa[0] * pa[0] + pa[1] * pow(c, pa[2]) + pa[3] * pow(c, pa[4])) / (1 + pa[0] * pa[0] + pa[1] + pa[3]));
        return 0.14 * pow((e * (1 + 3.64 / (e * pow(cst, 1.29)))), -2.7) * (1 / (1 + 1.1 * e * c / 115) + 0.054 / (1 + 1.1 * e * c / 850));
}

double MPrimaryGeneratorAction::cosmicNeutBeam(double t, double p) {
    // cosmic neutrons spectrum as a function of kinetic energy (GeV) and
    // zenith angle
    double massNeut = particleTable->FindParticle("neutron")->GetPDGMass() / GeV;
    double En = sqrt(p * p + massNeut * massNeut) - massNeut;
    double I0 = pow(En, -2.95);
    
    // From L.Faure Thesis
    // https://www.ge.infn.it/~batta/antinu/Tesi_Specialistica_Faure_Laura.pdf
    // Total normalization (neutrons/cm2/s):
    // 10-9 - 10-6:     2.0e-3
    // 10-6 - 2 MeV    30.0e-3
    // 2MeV - 1 GeV     4.6e-3
    // 1 GeV - 10 GeV   0.078e-3
    // 10 GeV -100 GeV  0.139 e-6
    double A = 1.006 * pow(10, -6);
    double B = 1.011 * pow(10, -3);
    double C = 1.53023e-7;
    
    if (En > 1e-9 && En <= 1e-6) I0 = 2e3;
    else if (En > 1e-6 && En <= 2) I0 = 2e-3 / En;
    else if (En > 2 && En <= 1000)
        I0 = A * exp(2.1451 * log(En) - 0.35 * log(En) * log(En)) +
        B * exp(-0.667 * log(En) - 0.4106 * log(En) * log(En));
    else if (En > 1000 && En <= 10000) I0 = C * pow(En / 1000., -2.95);
    
    return I0 * pow(cos(t), 3.5);
}


void MPrimaryGeneratorAction::setParticleFromPars(int p, int pindex, int type, int pdef, double px, double py, double pz, double Vx, double Vy, double Vz,
                                                  G4Event *anEvent, int A, int Z) {

    // PDG numbering scheme:
    // 2021 version: https://pdg.lbl.gov/2021/mcdata/mc_particle_id_contents.html
    // https://pdg.lbl.gov
    // Translating from geant3: https://www.star.bnl.gov/public/comp/simu/newsite/gstar/Manual/particle_id.html
    // 45 (deuteron), 46 (triton), 47 (alpha), 49 (helium3)

    if (type == 1 && pindex == p + 1) {
        // Primary Particle

        if (pdef != 80000) {
            int redefinedPdef = pdef;

            if (pdef == 45) {
                redefinedPdef = 1000010020;  // deuteron
            } else if (pdef == 46) {
                redefinedPdef = 1000010030;  // triton
            } else if (pdef == 47) {
                redefinedPdef = 1000020040;  // alpha
            } else if (pdef == 49) {
                redefinedPdef = 1000020030;  // He3
            }

            Particle = particleTable->FindParticle(redefinedPdef);
            if (!Particle) {
                cout << hd_msg << " Particle id " << redefinedPdef << " not found in G4 table." << endl << endl;

                cout << "Exiting !" << endl;
                exit(1);
            }
        } else {
            Particle = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
        }
        particleGun->SetParticleDefinition(Particle);

        // 4-momenta
        G4ThreeVector pmom(px * GeV, py * GeV, pz * GeV);
        double Mom = pmom.mag();
        double Phi = pmom.getPhi();
        double Theta = pmom.getTheta();
        double mass = Particle->GetPDGMass();
        double akine = sqrt(Mom * Mom + mass * mass) - mass;

        particleGun->SetParticleEnergy(akine);
        particleGun->SetParticleMomentumDirection(G4ThreeVector(cos(Phi / rad) * sin(Theta / rad), sin(Phi / rad) * sin(Theta / rad), cos(Theta / rad)));

        // vertex
        G4ThreeVector beam_vrt(Vx, Vy, Vz);
        particleGun->SetParticlePosition(beam_vrt);


        // beam polarization only along the beam
        // only for the first particle
        if (p == 0) {
            particleGun->SetParticlePolarization(G4ThreeVector(0, 0, beamPol));
        }

        if (GEN_VERBOSITY > 3)
            cout << hd_msg << " Particle n. " << p + 1 << ", id=" << pdef << " (" << Particle->GetParticleName() << ")"
                 << "  Vertex=" << beam_vrt / cm << "cm,  momentum=" << pmom / GeV << " GeV" << endl;

        // Primary particle generated at the chosen start time
        if (eventIndex > ntoskip) {
            particleGun->SetParticleTime(TSIGNAL);
            particleGun->SetNumberOfParticles(1);
            particleGun->GeneratePrimaryVertex(anEvent);
        }
    } else if (pindex != p + 1) {
        if (GEN_VERBOSITY > 3)
            cout << hd_msg << " Warning: file particle index " << pindex << " does not match read particle index " << p + 1 << endl;
    }
}

void MPrimaryGeneratorAction::setParticleFromParsPropagateTime(int p, vector <userInforForParticle> Userinfo, G4Event *anEvent, int A, int Z) {

    int pindex = Userinfo[p].infos[0];
    int type = Userinfo[p].infos[2];
    int pdef = Userinfo[p].infos[3];
    double px = Userinfo[p].infos[6];
    double py = Userinfo[p].infos[7];
    double pz = Userinfo[p].infos[8];
    double Vx = Userinfo[p].infos[11];
    double Vy = Userinfo[p].infos[12];
    double Vz = Userinfo[p].infos[13];
    //Make a list of particles parents
    vector<int> parentindex;
    parentindex.push_back(Userinfo[p].infos[4]);


    if (type == 1 && pindex == p + 1) {
        if (pdef != 80000) {

            int redefinedPdef = pdef;

            if (pdef == 45) {
                redefinedPdef = 1000010020;  // deuteron
            } else if (pdef == 46) {
                redefinedPdef = 1000010030;  // triton
            } else if (pdef == 47) {
                redefinedPdef = 1000020040;  // alpha
            } else if (pdef == 49) {
                redefinedPdef = 1000020030;  // He3
            }

            Particle = particleTable->FindParticle(redefinedPdef);
            if (!Particle) {
                cout << hd_msg << " Particle id " << redefinedPdef << " not found in G4 table." << endl << endl;

                cout << "Exiting !" << endl;
                exit(1);
            }
        } else {
            Particle = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
        }
        particleGun->SetParticleDefinition(Particle);

        // 4-momenta
        G4ThreeVector pmom(px * GeV, py * GeV, pz * GeV);
        double Mom = pmom.mag();
        double Phi = pmom.getPhi();
        double Theta = pmom.getTheta();
        double mass = Particle->GetPDGMass();
        double akine = sqrt(Mom * Mom + mass * mass) - mass;

        particleGun->SetParticleEnergy(akine);
        particleGun->SetParticleMomentumDirection(G4ThreeVector(cos(Phi / rad) * sin(Theta / rad), sin(Phi / rad) * sin(Theta / rad), cos(Theta / rad)));

        // vertex
        G4ThreeVector beam_vrt(Vx, Vy, Vz);
        particleGun->SetParticlePosition(beam_vrt);


        // beam polarization only along the beam
        // only for the first particle
        if (p == 0) {
            particleGun->SetParticlePolarization(G4ThreeVector(0, 0, beamPol));
        }
        if (GEN_VERBOSITY > 3) {
            cout << hd_msg << " Particle Number:  " << p + 1 << ", id=" << pdef << " (" << Particle->GetParticleName() << ")"
                 << "  Vertex=" << beam_vrt / cm << "cm,  momentum=" << pmom / GeV << " GeV" << endl;
        }

        // Primary particle generated in the middle of Time windowat the chosen time, while non primary particles have a time offset
        if (eventIndex > ntoskip) {
            double timeoffset = 0;
            //determine if the particle has an inactive parent
            if (parentindex[0] != 0 && Userinfo[parentindex[0] - 1].infos[2] != 1) {

                double px_parent = Userinfo[parentindex[0] - 1].infos[6];
                double py_parent = Userinfo[parentindex[0] - 1].infos[7];
                double pz_parent = Userinfo[parentindex[0] - 1].infos[8];
                double E_parent = Userinfo[parentindex[0] - 1].infos[9];
                double Vx_parent = Userinfo[parentindex[0] - 1].infos[11] * cm;
                double Vy_parent = Userinfo[parentindex[0] - 1].infos[12] * cm;
                double Vz_parent = Userinfo[parentindex[0] - 1].infos[13] * cm;
                //Push back the parents parent to the list
                parentindex.push_back(Userinfo[parentindex[0] - 1].infos[4]);

                //vertex difference between particle and its parent
                double vertex_diff = sqrt(pow(Vx - Vx_parent, 2) + pow(Vy - Vy_parent, 2) + pow(Vz - Vz_parent, 2));
                //set parent momentum
                G4ThreeVector pmom_parent(px_parent * GeV, py_parent * GeV, pz_parent * GeV);
                double Mom_parent = pmom_parent.mag();
                //calculate beta of parent
                double beta_parent = Mom_parent / (E_parent * GeV);
                double speedoflight = 29.979246 * (cm / ns);
                //calculate time between particle and its parent
                timeoffset += vertex_diff / (beta_parent * speedoflight);

                //checking for any further parents
                for (int i = 1;; i++) {
                    if (parentindex[i] != 0 && userInfo[parentindex[i - 1] - 1].infos[2] != 1) {

                        double Vx_parent1 = Userinfo[parentindex[i - 1] - 1].infos[11];
                        double Vy_parent1 = Userinfo[parentindex[i - 1] - 1].infos[12];
                        double Vz_parent1 = Userinfo[parentindex[i - 1] - 1].infos[13];

                        double px_parent2 = Userinfo[parentindex[i] - 1].infos[6];
                        double py_parent2 = Userinfo[parentindex[i] - 1].infos[7];
                        double pz_parent2 = Userinfo[parentindex[i] - 1].infos[8];
                        double E_parent2 = Userinfo[parentindex[i] - 1].infos[9];
                        double Vx_parent2 = Userinfo[parentindex[i] - 1].infos[11];
                        double Vy_parent2 = Userinfo[parentindex[i] - 1].infos[12];
                        double Vz_parent2 = Userinfo[parentindex[i] - 1].infos[13];
                        //push back further Parents to the list
                        parentindex.push_back(Userinfo[parentindex[i] - 1].infos[4]);
                        //vertex difference between the two parents
                        double Parents_vertex_diff =
                                sqrt(pow(Vx_parent1 - Vx_parent2, 2) +
                                     pow(Vy_parent1 - Vy_parent2, 2) +
                                     pow(Vz_parent1 - Vz_parent2, 2));

                        //set the second parents momentum
                        G4ThreeVector pmom_parent2(px_parent2 * GeV, py_parent2 * GeV, pz_parent2 * GeV);
                        double Mom_parent2 = pmom_parent2.mag();
                        //calculate beta of the second parent
                        double beta_parent2 = Mom_parent2 / (E_parent2 * GeV);
                        //calculate additional time offset
                        double timeoffset2 = Parents_vertex_diff / (beta_parent2 * speedoflight);
                        timeoffset += timeoffset2;

                    } else {
                        break;
                    }
                }


            }

            particleGun->SetParticleTime(TSIGNAL + timeoffset);
            particleGun->SetNumberOfParticles(1);
            particleGun->GeneratePrimaryVertex(anEvent);

        }

    } else if (pindex != p + 1) {
        if (GEN_VERBOSITY > 3)
            cout << hd_msg << " Warning: file particle index " << pindex << " does not match read particle index " << p + 1 << endl;
    }
}
