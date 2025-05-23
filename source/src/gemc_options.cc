// gemc headers
#include "gemcOptions.h"

void goptions::setGoptions()
{
	// Initialize all the options in the map<string, aopt>
	//
	// The "string" of the pair in the map is the argument name: -"name"=
	// args = the string type argument
	// arg  = the numeric type argument
	// help = Long explanation.
	// name = Short description.
	// type = 1 for argumenst that are strings, 0 for numbers.
	// argsJSONTypes: S = string, F = float, VS = vector of strings
	//
	//
	//
	//
	//
	//
	//
	//
	//

	// Generator
	// ---------

	optMap["BEAM_P"].args  = "e-, 11*GeV, 0*deg, 0*deg";
	optMap["BEAM_P"].help  = "Beam particle, momentum, angles (in respect of z-axis). \n";
	optMap["BEAM_P"].help += "      Example: -BEAM_P=\"e-, 6*GeV, 15*deg, 20*deg\" sets 6 GeV electrons 15 degrees in theta, 20 degrees in phi. \n";
	optMap["BEAM_P"].help += "      Use -BEAM_P=\"show_all\" to print the list of G4 supported particles.\n";
	optMap["BEAM_P"].help += "      An optional string \"KE\" can be appended to this option to specify Kinetic energy instead of momentum.\n";
	optMap["BEAM_P"].name  = "Primary particle, Energy, Theta, Phi";
	optMap["BEAM_P"].type  = 1;
	optMap["BEAM_P"].ctgr  = "generator";
	optMap["BEAM_P"].argsJSONDescription  = "particleID, momentum, theta, phi";
	optMap["BEAM_P"].argsJSONTypes        = "S, F, F, F";

	optMap["SPREAD_P"].args  = "0*GeV, 0*deg, 0*deg";
	optMap["SPREAD_P"].help  = "Spread Primary Particle energy and angles (in respect of z-axis). \n";
	optMap["SPREAD_P"].help += "      Example: -SPREAD_P=\"0*GeV, 10*deg, 20*deg\" spreads 10 degrees in theta, 20 degrees in phi. \n";
	optMap["SPREAD_P"].help += "      If the word flat is added as last argument to the option, theta is flat instead of cos(theta). \n";
	optMap["SPREAD_P"].help += "      An optional argument \"KE\" can be appended to this option to specify Kinetic energy instead of momentum.\n";
	optMap["SPREAD_P"].name  = "delta_Energy, delta_Theta, delta_phi";
	optMap["SPREAD_P"].type  = 1;
	optMap["SPREAD_P"].ctgr  = "generator";
	optMap["SPREAD_P"].argsJSONDescription  = "denergy, dtheta, dphi";
	optMap["SPREAD_P"].argsJSONTypes        = "F, F, F";

	optMap["ALIGN_ZAXIS"].args  = "no";
	optMap["ALIGN_ZAXIS"].help  = "Align z axis to a custom direction. Options:\n";
	optMap["ALIGN_ZAXIS"].help += "      - \"beamp\"  aligns z axis to the beam directions specified by BEAM_P.\n";
	optMap["ALIGN_ZAXIS"].help += "      - \"custom, theta*unit, phi*unit\" aligns z axis to a custom direction, changes BEAM_P reference frame.";
	optMap["ALIGN_ZAXIS"].name  = "Align z axis to a custom direction.";
	optMap["ALIGN_ZAXIS"].type  = 1;
	optMap["ALIGN_ZAXIS"].ctgr  = "generator";
	optMap["ALIGN_ZAXIS"].argsJSONDescription  = "to:";
	optMap["ALIGN_ZAXIS"].argsJSONTypes        = "S";

	optMap["BEAM_V"].args = "(0, 0, 0)cm";
	optMap["BEAM_V"].help = "Primary Particle Vertex. Example: -BEAM_V=\"(0, 0, -20)cm\". ";
	optMap["BEAM_V"].name = "Primary Particle Vertex";
	optMap["BEAM_V"].type = 1;
	optMap["BEAM_V"].ctgr = "generator";
	optMap["BEAM_V"].argsJSONDescription  = "vx vy vz unit";
	optMap["BEAM_V"].argsJSONTypes        = "F F F S";

	optMap["SPREAD_V"].args = "(0, 0)cm";
	optMap["SPREAD_V"].help = "Spread Primary Particle Radius, Z position. Example: -SPREAD_V=\"(0.1, 10)cm\". ";
	optMap["SPREAD_V"].name = "Primary Particle Vertex Spread";
	optMap["SPREAD_V"].type = 1;
	optMap["SPREAD_V"].ctgr = "generator";
	optMap["SPREAD_V"].argsJSONDescription  = "dr dz unit";
	optMap["SPREAD_V"].argsJSONTypes  = "F F S";

	optMap["POLAR"].args  = "100, 0*deg, 0*deg";
	optMap["POLAR"].help  = "Primary Particle polarization percentage and angles  (in respect of z-axis). \n";
	optMap["POLAR"].help += "      Example: -POLAR=\"90, 90*deg, 270*deg\" sets 90% polarization 90 degrees in theta, 270 degrees in phi. \n";
	optMap["POLAR"].help += "      Use -POLAR=\"show_all\" to print the list of G4 supported particles.\n";
	optMap["POLAR"].name  = "Primary Particle polarization in %, Theta, Phi";
	optMap["POLAR"].type  = 1;
	optMap["POLAR"].ctgr  = "generator";
	optMap["POLAR"].argsJSONDescription  = "percentage theta phi";
	optMap["POLAR"].argsJSONTypes  = "F F F";

	optMap["N"].arg  = 0;
	optMap["N"].help = "Number of events to be simulated.";
	optMap["N"].name = "Number of events to be simulated";
	optMap["N"].type = 0;
	optMap["N"].ctgr = "generator";
	optMap["N"].argsJSONDescription  = "nevents";
	optMap["N"].argsJSONTypes  = "F";

	optMap["EVTN"].arg  = 1;
	optMap["EVTN"].help = "Starting Event Number";
	optMap["EVTN"].name = "Starting Event Number";
	optMap["EVTN"].type = 0;
	optMap["EVTN"].ctgr = "generator";
	optMap["EVTN"].argsJSONDescription  = "startEventNumber";
	optMap["EVTN"].argsJSONTypes  = "F";

	optMap["ION_P"].args  = "no";
	optMap["ION_P"].help  = "Primary particle is an ion. Sets Z, A and Charge of the primary particle. Setting the charge is optional. \n";
	optMap["ION_P"].help += "      Example 1: -ION_P=\"6 12\" sets the Z to 6, A to 12. Charge will be Z. \n";
	optMap["ION_P"].help += "      Example 2: -ION_P=\"6 12 -4\" sets the Z to 6, A to 12, and charge to -4. \n";
	optMap["ION_P"].help += "      Notice: BEAM_P still sets the particle momentum, and its type must be GenericIon\n";
	optMap["ION_P"].name  = "Primary particle is an ion. Sets A, Z and Charge of the primary particle. ";
	optMap["ION_P"].type  = 1;
	optMap["ION_P"].ctgr  = "generator";
	optMap["ION_P"].argsJSONDescription  = "isPrimaryAPion";
	optMap["ION_P"].argsJSONTypes  = "S";

	optMap["INPUT_GEN_FILE"].args = "gemc_internal";
	optMap["INPUT_GEN_FILE"].help = "Generator Input. Current availables file formats:\n";
	optMap["INPUT_GEN_FILE"].help += "      LUND. \n";
	optMap["INPUT_GEN_FILE"].help += "      example: -INPUT_GEN_FILE=\"LUND, input.dat\" or -INPUT_GEN_FILE=\"StdHEP, darkphoton.stdhep\" \n";
	optMap["INPUT_GEN_FILE"].name = "Generator Input File";
	optMap["INPUT_GEN_FILE"].type = 1;
	optMap["INPUT_GEN_FILE"].ctgr = "generator";
	optMap["INPUT_GEN_FILE"].argsJSONDescription  = "type, filename";
	optMap["INPUT_GEN_FILE"].argsJSONTypes  = "S S";

	optMap["SHIFT_LUND_VERTEX"].args = "(0, 0, 0)cm";
	optMap["SHIFT_LUND_VERTEX"].help = "Shift Generator File tracks vertices.\n";
	optMap["SHIFT_LUND_VERTEX"].help += "      example: -SHIFT_LUND_VERTEX=\"(0, 0, -3)cm\" \n";
	optMap["SHIFT_LUND_VERTEX"].name = "Shift Generator File tracks vertices.";
	optMap["SHIFT_LUND_VERTEX"].type = 1;
	optMap["SHIFT_LUND_VERTEX"].ctgr = "generator";
	optMap["SHIFT_LUND_VERTEX"].argsJSONDescription  = "svx svy svz unit";
	optMap["SHIFT_LUND_VERTEX"].argsJSONTypes  = "F F F S";

	optMap["RANDOMIZE_LUND_VZ"].args = "0*cm, 0*cm, default";
	optMap["RANDOMIZE_LUND_VZ"].help = " Randomizes the z vertexes using, in order: Z shift, DZ sigma\n";
	optMap["RANDOMIZE_LUND_VZ"].help += "     By default the randomization is relative to the LUND vertex values. \n";
	optMap["RANDOMIZE_LUND_VZ"].help += "     If the third argument \"reset\" is given, the vertexes are relative to VZ=0 \n";
	optMap["RANDOMIZE_LUND_VZ"].help += "      example 1: -RANDOMIZE_LUND_VZ=\"-3*cm, 5*cm\" \n";
	optMap["RANDOMIZE_LUND_VZ"].help += "       This randomizes the z vertex by plus-minus 5cm around the original LUND values and shift it by -3cm \n";
	optMap["RANDOMIZE_LUND_VZ"].help += "      example 1: -RANDOMIZE_LUND_VZ=\"-3*cm, 5*cm, reset\" \n";
	optMap["RANDOMIZE_LUND_VZ"].help += "       This randomizes the z vertex by 5cm around vz = 0 and shift it by -3cm \n";
	optMap["RANDOMIZE_LUND_VZ"].name = "Randomizes the z vertex according to a shift and a flat value";
	optMap["RANDOMIZE_LUND_VZ"].type = 1;
	optMap["RANDOMIZE_LUND_VZ"].ctgr = "generator";
	optMap["RANDOMIZE_LUND_VZ"].argsJSONDescription  = "lundvzs lundvzd default/reset";
	optMap["RANDOMIZE_LUND_VZ"].argsJSONTypes  = "F F S";

	optMap["RASTER_VERTEX"].args = "0*cm, 0*cm, default";
	optMap["RASTER_VERTEX"].help = "Randomizes the x, y generated partice vertexes within an ellipse defined by the x, y radii, around their values.\n";
	optMap["RASTER_VERTEX"].help += "     If the third argument \"reset\" is given, the vertexes are relative to (VX, VY) = (0, 0) \n";
	optMap["RASTER_VERTEX"].help += "      example 1: -RASTER_VERTEX=\"2*cm, 3*cm\" \n";
	optMap["RASTER_VERTEX"].help += "       This randomizes the vertexes around the original LUND values. \n";
	optMap["RASTER_VERTEX"].help += "      example 2: -RASTER_VERTEX=\"2*cm, 3*cm, reset\" \n";
	optMap["RASTER_VERTEX"].help += "       This randomizes the vertexes around (VX, VY) = (0, 0). \n";
	optMap["RASTER_VERTEX"].name = "Randomizes the x, y generated partice vertexes";
	optMap["RASTER_VERTEX"].type = 1;
	optMap["RASTER_VERTEX"].ctgr = "generator";
	optMap["RASTER_VERTEX"].argsJSONDescription  = "rvdx rvdy default/reset";
	optMap["RASTER_VERTEX"].argsJSONTypes  = "F F S";

	optMap["BEAM_SPOT"].args = "0*cm, 0*cm, 0*cm, 0*cm, 0*deg, default";
	optMap["BEAM_SPOT"].help = " Randomizes the vertexes using, in order: X, Y shifts, DX DY sigmas, elippse phi orientation\n";
	optMap["BEAM_SPOT"].help += "     By default the randomization is relative to the LUND vertex values. \n";
	optMap["BEAM_SPOT"].help += "     If the sixth argument \"reset\" is given, the vertexes are relative to (VX, VY) = (0, 0) \n";
	optMap["BEAM_SPOT"].help += "      example 1: -BEAM_SPOT=\"2*cm, 3*cm, 0.2*cm, 0.1*cm, 22*deg\" \n";
	optMap["BEAM_SPOT"].help += "       This randomizes the vertexes around the original LUND values, but shifted by (VX, VY) = (2, 3)cm \n";
	optMap["BEAM_SPOT"].help += "       A gaussian with sigmas (SX, SY) = (0.2, 0.1)cm are used, rotated 22 degrees around z.\n";
	optMap["BEAM_SPOT"].help += "      example 2: -BEAM_SPOT=\"2*cm, 3*cm, 0.2*cm, 0.1*cm, 22*deg, reset\" \n";
	optMap["BEAM_SPOT"].help += "       This randomizes the vertexes around (VX, VY) = (2, 3)cm \n";
	optMap["BEAM_SPOT"].help += "       A gaussian with sigmas (SX, SY) = (0.2, 0.1)cm are used, rotated 22 degrees around z.\n";
	optMap["BEAM_SPOT"].name = "Randomizes the x, y generated partice vertexes according to  gaussian";
	optMap["BEAM_SPOT"].type = 1;
	optMap["BEAM_SPOT"].ctgr = "generator";
	optMap["BEAM_SPOT"].argsJSONDescription  = "bssx bssy bsdx bsdy bsphi default/reset";
	optMap["BEAM_SPOT"].argsJSONTypes  = "F F F F F S";
	

	optMap["MERGE_LUND_BG"].args = "no";
	optMap["MERGE_LUND_BG"].help = "LUND Input file to merge background events\n";
	optMap["MERGE_LUND_BG"].help += "      example: -MERGE_LUND_BG=\"background.dat\" \n";
	optMap["MERGE_LUND_BG"].name = "LUND Input file to merge background events";
	optMap["MERGE_LUND_BG"].type = 1;
	optMap["MERGE_LUND_BG"].ctgr = "generator";
	optMap["MERGE_LUND_BG"].argsJSONDescription  = "bgLundfilename";
	optMap["MERGE_LUND_BG"].argsJSONTypes  = "S";


	optMap["MERGE_BGHITS"].args = "no";
	optMap["MERGE_BGHITS"].help = "ASCII file to merge background hits\n";
	optMap["MERGE_BGHITS"].help += "      example: -MERGE_BGHITS=\"background.dat\" \n";
	optMap["MERGE_BGHITS"].name = "ASCII file to merge background hits";
	optMap["MERGE_BGHITS"].type = 1;
	optMap["MERGE_BGHITS"].ctgr = "generator";
	optMap["MERGE_BGHITS"].argsJSONDescription  = "bgfilename";
	optMap["MERGE_BGHITS"].argsJSONTypes  = "S";

	optMap["NGENP"].arg  = 10;
	optMap["NGENP"].help = "Max Number of Generated Particles to save in the Output.";
	optMap["NGENP"].name = "Max Number of Generated Particles to save in the Output";
	optMap["NGENP"].type = 0;
	optMap["NGENP"].ctgr = "generator";
	optMap["NGENP"].argsJSONDescription  = "maxGenToSave";
	optMap["NGENP"].argsJSONTypes  = "F";

	optMap["STEER_BEAM"].arg = 0;
	optMap["STEER_BEAM"].type = 0;
	optMap["STEER_BEAM"].ctgr = "generator";
	optMap["STEER_BEAM"].name = "STEER_BEAM";
	optMap["STEER_BEAM"].help = "Steer the beam, and translate the vertex, of an StdHep file by the amount specified in Beam_P, Beam_V, Spread_V \n";
	optMap["STEER_BEAM"].argsJSONDescription  = "steerBeam";
	optMap["STEER_BEAM"].argsJSONTypes  = "F";

	optMap["SKIPNGEN"].arg  = 0;
	optMap["SKIPNGEN"].help = "Skip N events";
	optMap["SKIPNGEN"].name = "Skip N events";
	optMap["SKIPNGEN"].type = 0;
	optMap["SKIPNGEN"].ctgr = "generator";
	optMap["SKIPNGEN"].argsJSONDescription  = "nEventsToSkip";
	optMap["SKIPNGEN"].argsJSONTypes  = "F";

	optMap["PROPAGATE_DVERTEXTIME"].arg  = 0;
	optMap["PROPAGATE_DVERTEXTIME"].help = "Calculate propogation time of detached vertex events and fire them at this later time. \n";
	optMap["PROPAGATE_DVERTEXTIME"].help += "         0: Off (default)\n";
	optMap["PROPAGATE_DVERTEXTIME"].help += "         1: On\n";
	optMap["PROPAGATE_DVERTEXTIME"].name = "Calculate Propogation Time";
	optMap["PROPAGATE_DVERTEXTIME"].type = 0;
	optMap["PROPAGATE_DVERTEXTIME"].ctgr = "generator";
	optMap["PROPAGATE_DVERTEXTIME"].argsJSONDescription  = "propagateTimeFromDetachedVertex";
	optMap["PROPAGATE_DVERTEXTIME"].argsJSONTypes  = "F";
	//
	//
	//
	//
	//
	//
	//
	//
	//
	// Cosmic Rays
	// ----------_

	optMap["COSMICRAYS"].args = "no";
	optMap["COSMICRAYS"].help = "Cosmic Generator. The model has a (cos(theta), p) probability function:\n\n";
	optMap["COSMICRAYS"].help += "              a^(b*cos(theta))/(c*p^2). \n\n";
	optMap["COSMICRAYS"].help += "      The COSMICRAYS option sets the parameters and the momentum range in the last two numbers. \n";
	optMap["COSMICRAYS"].help += "      By default the parameters are: \n";
	optMap["COSMICRAYS"].help += "       a = 55.6: \n";
	optMap["COSMICRAYS"].help += "       b = 1.04: \n";
	optMap["COSMICRAYS"].help += "       c = 64: \n";
	optMap["COSMICRAYS"].help += "      One can use the defaults or set the pars with the options: \n";
	optMap["COSMICRAYS"].help += "      example 1: -COSMICRAYS=\"default, 1, 10\" will use the default parameterization, and momentum range [1-10] GeV \n";
	optMap["COSMICRAYS"].help += "      example 2: -COSMICRAYS=\"55, 2, 66, 3, 4\" will set the parameterization, and momentum range [3-4] GeV \n";
	optMap["COSMICRAYS"].name = "Cosmic Generator";
	optMap["COSMICRAYS"].type = 1;
	optMap["COSMICRAYS"].ctgr = "generator";
	optMap["COSMICRAYS"].argsJSONDescription  = "parameterization, minP, maxP";
	optMap["COSMICRAYS"].argsJSONTypes  = "S F F";

	optMap["COSMICAREA"].args = "0*cm, 0*cm, 0*cm, 50*cm";
	optMap["COSMICAREA"].help = "Target (x,y,z) location and radius of area of interest";
	optMap["COSMICAREA"].name = "Target (x,y,z) location and radius of area of interest";
	optMap["COSMICAREA"].type = 1;
	optMap["COSMICAREA"].ctgr = "generator";
	optMap["COSMICAREA"].argsJSONDescription  = "x, y, z, radius";
	optMap["COSMICAREA"].argsJSONTypes  = "F F F F";
	//
	//
	//
	//
	//
	//
	//
	//
	//
	// Luminosity Beam
	// ---------------
 	optMap["LUMI_EVENT"].args = "0, 0*ns, 2*ns, 0*ns";
 	optMap["LUMI_EVENT"].help = "Luminosity Particle Parameters: number of Particles/Event, Time Window, Time Between Bunches, Primary Beam Particle Time (optional)\n";
 	optMap["LUMI_EVENT"].help += "            Example: -LUMI_EVENT=\"10000, 120*ns, 2*ns, 45*ns\" simulate 10K particles per event distributed over 120 ns, at 2ns intervals, with beam particle time set to last bunch before 45 ns. \n";
	optMap["LUMI_EVENT"].name = "Luminosity Particle Parameters";
	optMap["LUMI_EVENT"].type = 1;
	optMap["LUMI_EVENT"].ctgr = "luminosity";
	optMap["LUMI_EVENT"].argsJSONDescription  = "nparticles, timeWindow, bunchTime, signalTime";
	optMap["LUMI_EVENT"].argsJSONTypes  = "S F F F";

	optMap["LUMI_P"].args  = "e-, 11*GeV, 0*deg, 0*deg";
	optMap["LUMI_P"].help  = "Luminosity Particle, momentum, angles (in respect of z-axis). \n";
	optMap["LUMI_P"].help += "            Example: -LUMI_P=\"proton, 1*GeV, 25*deg, 2*deg\" sets 1 GeV protons, 25 degrees in theta, 2 degrees in phi. \n";
	optMap["LUMI_P"].help += "            Use -LUMI_P=\"show_all\" to print the list of G4 supported particles.\n";
	optMap["LUMI_P"].name  = "Luminosity Particle, Energy, Theta, Phi";
	optMap["LUMI_P"].type  = 1;
	optMap["LUMI_P"].ctgr = "luminosity";
	optMap["LUMI_P"].argsJSONDescription  = "particleID, momentum, theta, phi";
	optMap["LUMI_P"].argsJSONTypes        = "S, F, F, F";

	optMap["LUMI_SPREAD_P"].args  = "0*GeV, 0*deg, 0*deg";
	optMap["LUMI_SPREAD_P"].help  = "Spread Luminosity Particle energy and angles (in respect of z-axis). \n";
	optMap["LUMI_SPREAD_P"].help += "      Example: -LUMI_SPREAD_P=\"0*GeV, 10*deg, 20*deg\" spreads 10 degrees in theta, 20 degrees in phi. \n";
	optMap["LUMI_SPREAD_P"].help += "      If the word flat is added to the option, theta is flat instead of cos(theta). \n";
	optMap["LUMI_SPREAD_P"].name  = "delta_Energy, delta_Theta, delta_phi";
	optMap["LUMI_SPREAD_P"].type  = 1;
	optMap["LUMI_SPREAD_P"].ctgr  = "generator";
	optMap["LUMI_SPREAD_P"].argsJSONDescription  = "denergy, dtheta, dphi";
	optMap["LUMI_SPREAD_P"].argsJSONTypes        = "F, F, F";

	optMap["LUMI_V"].args = "(0, 0, -20)cm";
	optMap["LUMI_V"].help = "Luminosity Particle Vertex. Example: -LUMI_V=\"(0, 0, -20)cm\". ";
	optMap["LUMI_V"].name = "Luminosity Particle Vertex";
	optMap["LUMI_V"].type = 1;
	optMap["LUMI_V"].ctgr = "luminosity";
	optMap["LUMI_V"].argsJSONDescription  = "vx vy vz unit";
	optMap["LUMI_V"].argsJSONTypes        = "F F F S";

	optMap["LUMI_SPREAD_V"].args = "(0, 0)cm";
	optMap["LUMI_SPREAD_V"].help = "Spread Luminosity Particle Radius, Z position. Example: -SPREAD_V=\"(0.1, 10)cm\". ";
	optMap["LUMI_SPREAD_V"].name = "Luminosity Particle Vertex Spread";
	optMap["LUMI_SPREAD_V"].type = 1;
	optMap["LUMI_SPREAD_V"].ctgr = "luminosity";
	optMap["LUMI_SPREAD_V"].argsJSONDescription  = "dr dz unit";
	optMap["LUMI_SPREAD_V"].argsJSONTypes  = "F F S";

	optMap["LUMI2_EVENT"].args = "0, 2*ns";
	optMap["LUMI2_EVENT"].help = "Luminosity Particle 2 Parameters: number of Particles/Event, Time Between Bunches. The Time Window is specified with the LUMI_EVENT flag\n";
	optMap["LUMI2_EVENT"].help += "            Example: -LUMI2_EVENT=\"10000, 2*ns\" simulate 10K particles per event at 2ns intervals. \n";
	optMap["LUMI2_EVENT"].name = "Luminosity Particle 2 Parameters";
	optMap["LUMI2_EVENT"].type = 1;
	optMap["LUMI2_EVENT"].ctgr = "luminosity";
	optMap["LUMI2_EVENT"].argsJSONDescription  = "nparticles, timeWindow, bunchTime";
	optMap["LUMI2_EVENT"].argsJSONTypes  = "S F F";

	optMap["LUMI2_P"].args  = "proton, 50*GeV, 175*deg, 180*deg";
	optMap["LUMI2_P"].help  = "Luminosity Particle 2, momentum, angles (in respect of z-axis). \n";
	optMap["LUMI2_P"].help += "            Example: -LUMI2_P=\"proton, 1*GeV, 25*deg, 2*deg\" sets 1 GeV protons, 25 degrees in theta, 2 degrees in phi. \n";
	optMap["LUMI2_P"].help += "            Use -LUMI2_P=\"show_all\" to print the list of G4 supported particles.\n";
	optMap["LUMI2_P"].name  = "Luminosity Particle 2, Energy, Theta, Phi";
	optMap["LUMI2_P"].type  = 1;
	optMap["LUMI2_P"].ctgr = "luminosity";
	optMap["LUMI2_P"].argsJSONDescription  = "particleID, momentum, theta, phi";
	optMap["LUMI2_P"].argsJSONTypes        = "S, F, F, F";

	optMap["LUMI2_SPREAD_P"].args  = "0*GeV, 0*deg, 0*deg";
	optMap["LUMI2_SPREAD_P"].help  = "Spread Luminosity Particle 2 energy and angles (in respect of z-axis). \n";
	optMap["LUMI2_SPREAD_P"].help += "      Example: -LUMI2_SPREAD_P=\"0*GeV, 10*deg, 20*deg\" spreads 10 degrees in theta, 20 degrees in phi. \n";
	optMap["LUMI2_SPREAD_P"].help += "      If the word flat is added to the option, theta is flat instead of cos(theta). \n";
	optMap["LUMI2_SPREAD_P"].name  = "delta_Energy, delta_Theta, delta_phi";
	optMap["LUMI2_SPREAD_P"].type  = 1;
	optMap["LUMI2_SPREAD_P"].ctgr  = "generator";
	optMap["LUMI2_SPREAD_P"].argsJSONDescription  = "denergy, dtheta, dphi";
	optMap["LUMI2_SPREAD_P"].argsJSONTypes        = "F, F, F";

	optMap["LUMI2_V"].args = "(4, 0, 50)cm";
	optMap["LUMI2_V"].help = "Luminosity Particle 2 Vertex. Example: -LUMI2_V=\"(0, 0, -20)cm\". ";
	optMap["LUMI2_V"].name = "Luminosity Particle 2 Vertex";
	optMap["LUMI2_V"].type = 1;
	optMap["LUMI2_V"].ctgr = "luminosity";
	optMap["LUMI2_V"].argsJSONDescription  = "vx vy vz unit";
	optMap["LUMI2_V"].argsJSONTypes        = "F F F S";

	optMap["LUMI2_SPREAD_V"].args = "(0, 0)cm";
	optMap["LUMI2_SPREAD_V"].help = "Spread Luminosity Particle 2 Radius, Z position. Example: -SPREAD_V=\"(0.1, 10)cm\". ";
	optMap["LUMI2_SPREAD_V"].name = "Luminosity Particle Vertex 2 Spread";
	optMap["LUMI2_SPREAD_V"].type = 1;
	optMap["LUMI2_SPREAD_V"].ctgr = "luminosity";
	optMap["LUMI2_SPREAD_V"].argsJSONDescription  = "denergy, dtheta, dphi";
	optMap["LUMI2_SPREAD_V"].argsJSONTypes        = "F, F, F";
	//
	//
	//
	//
	//
	//
	//
	//
	//
	// MySQL Database
	// --------------

	optMap["DBHOST"].args = "no";
	optMap["DBHOST"].help = "Selects mysql server host name.";
	optMap["DBHOST"].name = "mysql server host name";
	optMap["DBHOST"].type = 1;
	optMap["DBHOST"].ctgr = "mysql";
	optMap["DBHOST"].argsJSONDescription  = "mysqlServer";
	optMap["DBHOST"].argsJSONTypes        = "S";

	optMap["DATABASE"].args = "no";
	optMap["DATABASE"].help = "Selects mysql Database.";
	optMap["DATABASE"].name = "mysql Database";
	optMap["DATABASE"].type = 1;
	optMap["DATABASE"].ctgr = "mysql";
	optMap["DATABASE"].argsJSONDescription  = "mysqlDatabase";
	optMap["DATABASE"].argsJSONTypes        = "S";

	optMap["DBUSER"].args = "gemc";
	optMap["DBUSER"].help = "Select mysql user name";
	optMap["DBUSER"].name = "Select mysql user name";
	optMap["DBUSER"].type = 1;
	optMap["DBUSER"].ctgr = "mysql";
	optMap["DBUSER"].argsJSONDescription  = "mysqlUserName";
	optMap["DBUSER"].argsJSONTypes        = "S";

	optMap["DBPSWD"].args = "no";
	optMap["DBPSWD"].help = "mysql password";
	optMap["DBPSWD"].name = "Select mysql password";
	optMap["DBPSWD"].type = 1;
	optMap["DBPSWD"].ctgr = "mysql";
	optMap["DBPSWD"].argsJSONDescription  = "mysqlUserPwd";
	optMap["DBPSWD"].argsJSONTypes        = "S";

	optMap["DBPORT"].arg = 0;
	optMap["DBPORT"].help = "Select mysql server port.";
	optMap["DBPORT"].name = "Select mysql server port";
	optMap["DBPORT"].type = 0;
	optMap["DBPORT"].ctgr = "mysql";
	optMap["DBPORT"].argsJSONDescription  = "mysqlPort";
	optMap["DBPORT"].argsJSONTypes        = "F";
	//
	//
	//
	//
	//
	//
	//
	//
	//
	// Verbosity
	// ---------

	optMap["G4P_VERBOSITY"].arg  = 1;
	optMap["G4P_VERBOSITY"].help = "Controls Physical Volumes Construction Log Output.";
	optMap["G4P_VERBOSITY"].name = "Logical Volume Verbosity";
	optMap["G4P_VERBOSITY"].type = 0;
	optMap["G4P_VERBOSITY"].ctgr = "verbosity";
	optMap["G4P_VERBOSITY"].type = 0;
	optMap["G4P_VERBOSITY"].ctgr = "verbosity";
	optMap["G4P_VERBOSITY"].argsJSONDescription  = "G4P_VERBOSITY";
	optMap["G4P_VERBOSITY"].argsJSONTypes        = "F";

	optMap["GEO_VERBOSITY"].arg  = 1;
	optMap["GEO_VERBOSITY"].help = "Controls Geometry Construction Log Output.";
	optMap["GEO_VERBOSITY"].name = "Geometry Verbosity";
	optMap["GEO_VERBOSITY"].type = 0;
	optMap["GEO_VERBOSITY"].ctgr = "verbosity";
	optMap["GEO_VERBOSITY"].argsJSONDescription  = "GEO_VERBOSITY";
	optMap["GEO_VERBOSITY"].argsJSONTypes        = "F";

	optMap["GUI_VERBOSITY"].arg  = 1;
	optMap["GUI_VERBOSITY"].help = "Controls GUI Construction Log Output.";
	optMap["GUI_VERBOSITY"].name = "GUI Verbosity";
	optMap["GUI_VERBOSITY"].type = 0;
	optMap["GUI_VERBOSITY"].ctgr = "verbosity";
	optMap["GUI_VERBOSITY"].argsJSONDescription  = "GUI_VERBOSITY";
	optMap["GUI_VERBOSITY"].argsJSONTypes        = "F";

	optMap["HIT_VERBOSITY"].arg  = 1;
	optMap["HIT_VERBOSITY"].help = "Controls Hits Log Output. ";
	optMap["HIT_VERBOSITY"].name = "Hit Verbosity";
	optMap["HIT_VERBOSITY"].type = 0;
	optMap["HIT_VERBOSITY"].ctgr = "verbosity";
	optMap["HIT_VERBOSITY"].argsJSONDescription  = "HIT_VERBOSITY";
	optMap["HIT_VERBOSITY"].argsJSONTypes        = "F";

	optMap["CATCH"].args = "Maurizio";
	optMap["CATCH"].help = "Catch volumes matching the given string.";
	optMap["CATCH"].name = "Volume catcher";
	optMap["CATCH"].type = 1;
	optMap["CATCH"].ctgr = "verbosity";
	optMap["CATCH"].argsJSONDescription  = "CATCH";
	optMap["CATCH"].argsJSONTypes        = "S";

	optMap["FIELD_VERBOSITY"].arg   = 0;
	optMap["FIELD_VERBOSITY"].help  = "Controls Electro-Magnetic Fields Log Output:\n";
	optMap["FIELD_VERBOSITY"].help += "  0: no log";
	optMap["FIELD_VERBOSITY"].help += "  1: field definitions log";
	optMap["FIELD_VERBOSITY"].help += "  2: max field details";
	optMap["FIELD_VERBOSITY"].name  = "Electro-Magnetic Fields Verbosity";
	optMap["FIELD_VERBOSITY"].type  = 0;
	optMap["FIELD_VERBOSITY"].ctgr  = "verbosity";
	optMap["FIELD_VERBOSITY"].argsJSONDescription  = "FIELD_VERBOSITY";
	optMap["FIELD_VERBOSITY"].argsJSONTypes        = "F";

	optMap["PRINT_EVENT"].arg  = 1000;
	optMap["PRINT_EVENT"].help = "-PRINT_EVENT=N: Print Event Number every N events.";
	optMap["PRINT_EVENT"].name = "Print Event Modulus";
	optMap["PRINT_EVENT"].type = 0;
	optMap["PRINT_EVENT"].ctgr = "verbosity";
	optMap["PRINT_EVENT"].argsJSONDescription  = "PRINT_EVENT";
	optMap["PRINT_EVENT"].argsJSONTypes        = "F";

	optMap["BANK_VERBOSITY"].arg  = 1;
	optMap["BANK_VERBOSITY"].help = "Controls Bank Log Output.";
	optMap["BANK_VERBOSITY"].name = "Bank Output Verbosity";
	optMap["BANK_VERBOSITY"].type = 0;
	optMap["BANK_VERBOSITY"].ctgr = "verbosity";
	optMap["BANK_VERBOSITY"].argsJSONDescription  = "BANK_VERBOSITY";
	optMap["BANK_VERBOSITY"].argsJSONTypes        = "F";

	optMap["EVENT_VERBOSITY"].arg  = 0;
	optMap["EVENT_VERBOSITY"].help = "Controls Event Log Output.";
	optMap["EVENT_VERBOSITY"].name = "Event Output Verbosity";
	optMap["EVENT_VERBOSITY"].type = 0;
	optMap["EVENT_VERBOSITY"].ctgr = "verbosity";
	optMap["EVENT_VERBOSITY"].argsJSONDescription  = "EVENT_VERBOSITY";
	optMap["EVENT_VERBOSITY"].argsJSONTypes        = "F";

	optMap["PHY_VERBOSITY"].arg  = 1;
	optMap["PHY_VERBOSITY"].help = "Controls Physics List Log Output.";
	optMap["PHY_VERBOSITY"].name = "Physics List Verbosity";
	optMap["PHY_VERBOSITY"].type = 0;
	optMap["PHY_VERBOSITY"].ctgr = "verbosity";
	optMap["PHY_VERBOSITY"].argsJSONDescription  = "PHY_VERBOSITY";
	optMap["PHY_VERBOSITY"].argsJSONTypes        = "F";

	optMap["GEN_VERBOSITY"].arg  = 0;
	optMap["GEN_VERBOSITY"].help = "Controls Geant4 Generator Verbosity.";
	optMap["GEN_VERBOSITY"].name = "Geant4 Generator Verbosity";
	optMap["GEN_VERBOSITY"].type = 0;
	optMap["GEN_VERBOSITY"].ctgr = "verbosity";
	optMap["GEN_VERBOSITY"].argsJSONDescription  = "GEN_VERBOSITY";
	optMap["GEN_VERBOSITY"].argsJSONTypes        = "F";

	optMap["G4TRACK_VERBOSITY"].arg  = 0;
	optMap["G4TRACK_VERBOSITY"].help = "Controls Geant4 Track Verbosity.";
	optMap["G4TRACK_VERBOSITY"].name = "Geant4 Track Verbosity";
	optMap["G4TRACK_VERBOSITY"].type = 0;
	optMap["G4TRACK_VERBOSITY"].ctgr = "verbosity";
	optMap["G4TRACK_VERBOSITY"].argsJSONDescription  = "G4TRACK_VERBOSITY";
	optMap["G4TRACK_VERBOSITY"].argsJSONTypes        = "F";

	optMap["MATERIAL_VERBOSITY"].arg  = 0;
	optMap["MATERIAL_VERBOSITY"].help = "Controls Geant4 Material Verbosity.";
	optMap["MATERIAL_VERBOSITY"].name = "Geant4 Material Verbosity";
	optMap["MATERIAL_VERBOSITY"].type = 0;
	optMap["MATERIAL_VERBOSITY"].ctgr = "verbosity";
	optMap["MATERIAL_VERBOSITY"].argsJSONDescription  = "MATERIAL_VERBOSITY";
	optMap["MATERIAL_VERBOSITY"].argsJSONTypes        = "F";

	optMap["PARAMETER_VERBOSITY"].arg  = 0;
	optMap["PARAMETER_VERBOSITY"].help = "Controls Parameters Verbosity.";
	optMap["PARAMETER_VERBOSITY"].name = "Parameters Verbosity";
	optMap["PARAMETER_VERBOSITY"].type = 0;
	optMap["PARAMETER_VERBOSITY"].ctgr = "verbosity";
	optMap["PARAMETER_VERBOSITY"].argsJSONDescription  = "PARAMETER_VERBOSITY";
	optMap["PARAMETER_VERBOSITY"].argsJSONTypes        = "F";

	optMap["MIRROR_VERBOSITY"].arg  = 0;
	optMap["MIRROR_VERBOSITY"].help = "Controls Mirrors Verbosity.";
	optMap["MIRROR_VERBOSITY"].name = "Mirrors Verbosity";
	optMap["MIRROR_VERBOSITY"].type = 0;
	optMap["MIRROR_VERBOSITY"].ctgr = "verbosity";
	optMap["MIRROR_VERBOSITY"].argsJSONDescription  = "MIRROR_VERBOSITY";
	optMap["MIRROR_VERBOSITY"].argsJSONTypes        = "F";

	//
	//
	//
	//
	//
	//
	//
	//
	//
	// Run Control
	// -----------

	optMap["EXEC_MACRO"].args = "no";
	optMap["EXEC_MACRO"].help = "Executes commands in macro file.";
	optMap["EXEC_MACRO"].name = "Executes commands in macro file";
	optMap["EXEC_MACRO"].type = 1;
	optMap["EXEC_MACRO"].ctgr = "control";
	optMap["EXEC_MACRO"].repe  = 1;

	optMap["CHECK_OVERLAPS"].arg  = 0;
	optMap["CHECK_OVERLAPS"].help  = "Checks Overlapping Volumes:\n";
	optMap["CHECK_OVERLAPS"].help += "      1.  Check Overlaps at Construction Time\n";
	optMap["CHECK_OVERLAPS"].help += "      2.  Check Overlaps based on standard lines grid setup\n";
	optMap["CHECK_OVERLAPS"].help += "      3.  Check Overlaps by shooting lines according to a cylindrical pattern\n";
	optMap["CHECK_OVERLAPS"].name = "Checks Overlapping Volumes";
	optMap["CHECK_OVERLAPS"].type = 0;
	optMap["CHECK_OVERLAPS"].ctgr = "control";
	
	optMap["USE_GUI"].arg   = 1;
	optMap["USE_GUI"].help  = " GUI switch\n";
	optMap["USE_GUI"].help += "      0.  Don't use the graphical interface\n";
	optMap["USE_GUI"].help += "      1.  OGLSQt (OGL, OGLS, OGLSQt)\n";
	optMap["USE_GUI"].help += "      2.  OGLIQt (OGLI, OGLIQt) \n";
	optMap["USE_GUI"].name  = "QT Gui";
	optMap["USE_GUI"].type  = 0;
	optMap["USE_GUI"].ctgr  = "control";
	
	optMap["geometry"].args="620x620";
	optMap["geometry"].help = "Specify the size of the QT display window. Default '600x600' ";
	optMap["geometry"].name="geometry";
	optMap["geometry"].type=1;
	optMap["geometry"].ctgr = "control";
	
	optMap["GUIPOS"].args="(150, 150)";
	optMap["GUIPOS"].help = "Specify the position of the QT display window. Default is at 50, 50 ";
	optMap["GUIPOS"].name="geometry";
	optMap["GUIPOS"].type=1;
	optMap["GUIPOS"].ctgr = "control";
	
	optMap["QTSTYLE"].args  = "no";
	optMap["QTSTYLE"].name  = "Sets the GUI Style";
	optMap["QTSTYLE"].help  = "Sets the GUI Style. Available options: \n";
	optMap["QTSTYLE"].help += "      - QCleanlooksStyle \n";
	optMap["QTSTYLE"].help += "      - QMacStyle \n";
	optMap["QTSTYLE"].help += "      - QPlastiqueStyle \n";
	optMap["QTSTYLE"].help += "      - QWindowsStyle \n";
	optMap["QTSTYLE"].help += "      - QMotifStyle";
	optMap["QTSTYLE"].type  = 1;
	optMap["QTSTYLE"].ctgr  = "control";
	
	optMap["RANDOM"].args = "TIME";
	optMap["RANDOM"].help = "Random Engine Initialization. The argument (seed) can be an integer or the string TIME.";
	optMap["RANDOM"].name = "Random Engine Initialization";
	optMap["RANDOM"].type = 1;
	optMap["RANDOM"].ctgr = "control";
	
	optMap["gcard"].args = "no";
	optMap["gcard"].help = "gemc card file.";
	optMap["gcard"].name = "gemc card file";
	optMap["gcard"].type = 1;
	optMap["gcard"].ctgr = "control";
	
	optMap["EVN"].arg  = 1;
	optMap["EVN"].help = "Initial Event Number.";
	optMap["EVN"].name = "Initial Event Number";
	optMap["EVN"].type = 0;
	optMap["EVN"].ctgr = "control";

	optMap["ENERGY_CUT"].arg  = -1.;
	optMap["ENERGY_CUT"].help = "Set an energy cut in MeV below which no particle will be tracked further. -1. turns this off.";
	optMap["ENERGY_CUT"].name = "Tracking Energy cut (in MeV)";
	optMap["ENERGY_CUT"].type = 0;
	optMap["ENERGY_CUT"].ctgr = "control";
	
	optMap["MAX_X_POS"].arg  = 20000.;
	optMap["MAX_X_POS"].help = "Max X Position in millimeters. Beyond this the track will be killed";
	optMap["MAX_X_POS"].name = "Max X Position in millimeters. Beyond this the track will be killed.";
	optMap["MAX_X_POS"].type = 0;
	optMap["MAX_X_POS"].ctgr = "control";
	
	optMap["MAX_Y_POS"].arg  = 20000.;
	optMap["MAX_Y_POS"].help = "Max Y Position in millimeters. Beyond this the track will be killed";
	optMap["MAX_Y_POS"].name = "Max Y Position in millimeters. Beyond this the track will be killed.";
	optMap["MAX_Y_POS"].type = 0;
	optMap["MAX_Y_POS"].ctgr = "control";
	
	optMap["MAX_Z_POS"].arg  = 20000.;
	optMap["MAX_Z_POS"].help = "Max Z Position in millimeters. Beyond this the track will be killed";
	optMap["MAX_Z_POS"].name = "Max Z Position in millimeters. Beyond this the track will be killed.";
	optMap["MAX_Z_POS"].type = 0;
	optMap["MAX_Z_POS"].ctgr = "control";
	
	optMap["DAWN_N"].arg = 0;
	optMap["DAWN_N"].name = "Number of events to be displayed with the DAWN driver (also activate the DAWN driver)";
	optMap["DAWN_N"].help = "Number of events to be displayed with the DAWN driver (also activate the DAWN driver).";
	optMap["DAWN_N"].type = 0;
	optMap["DAWN_N"].ctgr = "control";
	
	optMap["HIT_PROCESS_LIST"].args = "clas12";
	optMap["HIT_PROCESS_LIST"].name = "Registers Hit Process Routines.";
	optMap["HIT_PROCESS_LIST"].help = "Registers Hit Process Routines. Can register multiple experiments, separated by space, e.v. \"clas12 aprime\"\n";
	optMap["HIT_PROCESS_LIST"].help += "      clas12.  CLAS12 hit process routines (default)\n";
	optMap["HIT_PROCESS_LIST"].help += "      aprime.  aprime hit process routines\n";
	optMap["HIT_PROCESS_LIST"].help += "      gluex.   GlueX  hit process routines\n";
	optMap["HIT_PROCESS_LIST"].type = 1;
	optMap["HIT_PROCESS_LIST"].ctgr = "control";
	
	optMap["SAVE_ALL_MOTHERS"].arg = 0;
	optMap["SAVE_ALL_MOTHERS"].name = "Set to 1 to save mother vertex and pid infos in output. High Memory Usage";
	optMap["SAVE_ALL_MOTHERS"].help  = "Set to 1 to save mother vertex and pid infos in output. High Memory Usage. Default is 0.\n";
	optMap["SAVE_ALL_MOTHERS"].help += "                       2: saves in LUND format any particle that generate hits\n";
	optMap["SAVE_ALL_MOTHERS"].help += "                       3: same as 2 but if a track’s mother also produced hits, only the mother is saved so the hits are not double counted\n";
	optMap["SAVE_ALL_MOTHERS"].type = 0;
	optMap["SAVE_ALL_MOTHERS"].ctgr = "control";
	
	optMap["SAVE_ALL_ANCESTORS"].arg = 0;
	optMap["SAVE_ALL_ANCESTORS"].name = "Set to 1 to save all ancestors of hits in output. High Memory Usage";
	optMap["SAVE_ALL_ANCESTORS"].help  = "Set to 1 to save all ancestors of hits. High Memory Usage. Default is 0.\n";
	optMap["SAVE_ALL_ANCESTORS"].type = 0;
	optMap["SAVE_ALL_ANCESTORS"].ctgr = "control";
	
	optMap["HIGH_RES"].arg = 1;
	optMap["HIGH_RES"].name = "Use High Resolution Graphics";
	optMap["HIGH_RES"].help = "Use High Resolution Graphics\n";
	optMap["HIGH_RES"].type = 0;
	optMap["HIGH_RES"].ctgr = "control";
	
	optMap["RECORD_PASSBY"].arg = 0;
	optMap["RECORD_PASSBY"].name = "Set to one if you want to save zero energy hits in the output";
	optMap["RECORD_PASSBY"].help = "Set to one if you want to save zero energy hits in the output. Default is 0.\n";
	optMap["RECORD_PASSBY"].type = 0;
	optMap["RECORD_PASSBY"].ctgr = "control";
	
	optMap["RECORD_OPTICALPHOTONS"].arg = 0;
	optMap["RECORD_OPTICALPHOTONS"].name = "Set to one if you want the optical photons to produce hits";
	optMap["RECORD_OPTICALPHOTONS"].help = "Set to one if you want the optical photons to produce hits. Default is 0.\n";
	optMap["RECORD_OPTICALPHOTONS"].type = 0;
	optMap["RECORD_OPTICALPHOTONS"].ctgr = "control";

	optMap["RECORD_MIRRORS"].arg = 0;
	optMap["RECORD_MIRRORS"].name = "Set to one if you want to save mirror hits in the output";
	optMap["RECORD_MIRRORS"].help = "Set to one if you want to save mirror hits in the output. Default is 0.\n";
	optMap["RECORD_MIRRORS"].type = 0;
	optMap["RECORD_MIRRORS"].ctgr = "control";

	optMap["RUNNO"].arg  = -1;
	optMap["RUNNO"].name = "Run Number. Controls the geometry and calibration parameters";
	optMap["RUNNO"].help = "Run Number. Controls the geometry and calibration parameters. Default is 1\n";
	optMap["RUNNO"].type = 0;
	optMap["RUNNO"].ctgr = "control";
	
	optMap["RUN_WEIGHTS"].args  = "no";
	optMap["RUN_WEIGHTS"].name  = "Simulate events based on run based conditions table (text file)";
	optMap["RUN_WEIGHTS"].help  = "Simulate events based on run based conditions table (text file)\n";
	optMap["RUN_WEIGHTS"].help += "      The text file must have two columns, run# and weight.\n";
	optMap["RUN_WEIGHTS"].help += "      For example:\n\n";
	optMap["RUN_WEIGHTS"].help += "        11 0.1\n";
	optMap["RUN_WEIGHTS"].help += "        12 0.7\n";
	optMap["RUN_WEIGHTS"].help += "        13 0.2\n\n";
	optMap["RUN_WEIGHTS"].help += "      Will simulate 10% of events with run number 11 conditions, 70% run 12 and 20% run 13.\n";
	optMap["RUN_WEIGHTS"].type  = 1;
	optMap["RUN_WEIGHTS"].ctgr  = "control";

	optMap["RFSETUP"].args = "no";
	optMap["RFSETUP"].name = "Radio-frequency signal";
	optMap["RFSETUP"].help = "Radio-frequency signal. This are a minium of 2 parameters for one given RF signal:\n";
	optMap["RFSETUP"].help += "      1. radioFrequency period T (ns). The frequency in GHz is 1/T \n";
	optMap["RFSETUP"].help += "      2. radioInterval (in number of bunches): distance between RF signals, or pre-scale factor.\n";
	optMap["RFSETUP"].help += "     Any additional parameter (in ns) will add an RF signal with that time distance from the original.\n";
	optMap["RFSETUP"].help += "     Example1:  4.008 (0.2495 GHz) RF signal, 1 output, output frequency is 50 bunches: \n";
	optMap["RFSETUP"].help += "      -RFSETUP=\"4.008, 50\" \n";
	optMap["RFSETUP"].help += "     Example2: two (0.5 GHz) and 2 RF signals, separated by 30 ns and both output frequency is 80 bunches: \n";
	optMap["RFSETUP"].help += "      -RFSETUP=\"2.004, 80, 30\" \n";
	optMap["RFSETUP"].help += "     By default the RFSETUP is set to 'clas12': the above constants are to be read from the CCDB database, using as run: RUNNO and as variation the one specified in the option DIGITIZATION_VARIATION.\n";
	optMap["RFSETUP"].type = 1;
	optMap["RFSETUP"].ctgr = "control";

	optMap["RFSTART"].args = "eventVertex, 0, 0, 0";
	optMap["RFSTART"].name = "Radio-frequency time model";
	optMap["RFSTART"].help = "Radio-frequency time model. Available options are:\n";
	optMap["RFSTART"].help += "      \"eventVertex, 0, 0, 0\" (default): the RF time is the event start time + the light time-distance of the first particle from the point (0,0,0)\n";
	optMap["RFSTART"].help += "      \"eventTime\".....................: the RF time is identical to the event start time\n";
	optMap["RFSTART"].ctgr = "control";

	optMap["DIGITIZATION_VARIATION"].args = "default";
	optMap["DIGITIZATION_VARIATION"].name = "Common variation to be used in the digitization routines";
	optMap["DIGITIZATION_VARIATION"].help = "Common variation to be used in the digitization routines\n";
	optMap["DIGITIZATION_VARIATION"].type = 1;
	optMap["DIGITIZATION_VARIATION"].ctgr = "control";

	optMap["DIGITIZATION_TIMESTAMP"].args = "no";
	optMap["DIGITIZATION_TIMESTAMP"].name = "Timestamp to be used in digitization routines";
	optMap["DIGITIZATION_TIMESTAMP"].help = "Timestamp to be used in digitization routines\n";
	optMap["DIGITIZATION_TIMESTAMP"].type = 1;
	optMap["DIGITIZATION_TIMESTAMP"].ctgr = "control";

	optMap["HARDWARESTATUS"].arg  = 0;
	optMap["HARDWARESTATUS"].name = "Accounts for hardware status";
	optMap["HARDWARESTATUS"].help = "Accounts for hardware status\n";
	optMap["HARDWARESTATUS"].help += "This flag can be used by the digitization routines to account for hardware status\n";
	optMap["HARDWARESTATUS"].type = 0;
	optMap["HARDWARESTATUS"].ctgr = "control";

	optMap["DETECTOR_INEFFICIENCY"].arg  = 0;
	optMap["DETECTOR_INEFFICIENCY"].name = "Accounts for detector inefficiency";
	optMap["DETECTOR_INEFFICIENCY"].help = "Accounts for detector inefficiency\n";
	optMap["DETECTOR_INEFFICIENCY"].help += "This flag can be used by the digitization routines to account for hardware status\n";
	optMap["DETECTOR_INEFFICIENCY"].type = 0;
	optMap["DETECTOR_INEFFICIENCY"].ctgr = "control";

	optMap["APPLY_THRESHOLDS"].arg  = 0;
	optMap["APPLY_THRESHOLDS"].name = "Accounts for detector thresholds";
	optMap["APPLY_THRESHOLDS"].help = "Accounts for detector thresholds\n";
	optMap["APPLY_THRESHOLDS"].help += "This flag can be used by the digitization routines to account for hardware status\n";
	optMap["APPLY_THRESHOLDS"].type = 0;
	optMap["APPLY_THRESHOLDS"].ctgr = "control";
	


	// Output
	// ------
	optMap["OUTPUT"].args = "no, output";
	optMap["OUTPUT"].help = "Type of output, output filename. Supported output: evio, txt. Example: -OUTPUT=\"evio, out.ev\"";
	optMap["OUTPUT"].name = "Type of output, output filename. ";
	optMap["OUTPUT"].type = 1;
	optMap["OUTPUT"].ctgr = "output";
	
	optMap["INTEGRATEDRAW"].args = "no";
	optMap["INTEGRATEDRAW"].help = "Activates integrated geant4 true info output for system(s). Example: -INTEGRATEDRAW=\"dc, ftof\"";
	optMap["INTEGRATEDRAW"].help = "To acticate the true info for all sensitive detectors:  -INTEGRATEDRAW=\"*\"";
	optMap["INTEGRATEDRAW"].name = "Activates integrated geant4 true info output for system(s)";
	optMap["INTEGRATEDRAW"].type = 1;
	optMap["INTEGRATEDRAW"].ctgr = "output";
	
	optMap["INTEGRATEDDGT"].args = "yes";
	optMap["INTEGRATEDDGT"].help = "De-activates integrated digitized output for system(s). Example: -INTEGRATEDDGT=\"DC, TOF\"";
	optMap["INTEGRATEDDGT"].name = "De-activates integrated digitized output for system(s). ";
	optMap["INTEGRATEDDGT"].type = 1;
	optMap["INTEGRATEDDGT"].ctgr = "output";

	optMap["SKIPREJECTEDHITS"].arg  = 1;
	optMap["SKIPREJECTEDHITS"].help = "Skips hits that are rejected by digitization. Default: yes (1)";
	optMap["SKIPREJECTEDHITS"].name = "Skips hits that are rejected by digitization. Default: yes";
	optMap["SKIPREJECTEDHITS"].type = 0;
	optMap["SKIPREJECTEDHITS"].ctgr = "output";

	optMap["ALLRAWS"].args = "no";
	optMap["ALLRAWS"].help = "Activates step-by-step output for system(s). Example: -ALLRAWS=\"DC, TOF\"";
	optMap["ALLRAWS"].name = "Activates step-by-step output for system(s). ";
	optMap["ALLRAWS"].type = 1;
	optMap["ALLRAWS"].ctgr = "output";
	
	optMap["SIGNALVT"].args = "no";
	optMap["SIGNALVT"].help = "Activates voltage (t) output for system(s). Example: -SIGNALVT=\"DC, TOF\"";
	optMap["SIGNALVT"].name = "Activates voltage (t) output for system(s). ";
	optMap["SIGNALVT"].type = 1;
	optMap["SIGNALVT"].ctgr = "output";

	optMap["ELECTRONICNOISE"].args = "no";
	optMap["ELECTRONICNOISE"].help = "Activates electronic noise routines for system(s). Example: -ELECTRONICNOISE=\"dc, ftof\"";
	optMap["ELECTRONICNOISE"].name = "Activates electronic noise routines for system(s). ";
	optMap["ELECTRONICNOISE"].type = 1;
	optMap["ELECTRONICNOISE"].ctgr = "output";

	optMap["VTRESOLUTION"].arg = 0.1;
	optMap["VTRESOLUTION"].help = "Voltage versus time resolution, in ns";
	optMap["VTRESOLUTION"].name = "Voltage versus time resolution, in ns.";
	optMap["VTRESOLUTION"].type = 0;
	optMap["VTRESOLUTION"].ctgr = "output";
	
	optMap["FILTER_NULL_VARIABLES"].arg = 0;
	optMap["FILTER_NULL_VARIABLES"].help = "If set to 1, filters out variables with no valuable info in the output.";
	optMap["FILTER_NULL_VARIABLES"].name = "If set to 1, filters out variables with no valuable info in the output.";
	optMap["FILTER_NULL_VARIABLES"].type = 0;
	optMap["FILTER_NULL_VARIABLES"].ctgr = "output";

	optMap["FILTER_HITS"].arg = 0;
	optMap["FILTER_HITS"].help = "If set to 1, do not write output if there are no hits in the detectors";
	optMap["FILTER_HITS"].name = "If set to 1, do not write output if there are no hits in the detectors";
	optMap["FILTER_HITS"].type = 0;
	optMap["FILTER_HITS"].ctgr = "output";

	optMap["FILTER_HADRONS"].arg = 0;
	optMap["FILTER_HADRONS"].help = "If set to 1, do not write events if there are no hadrons. Otherwise if \n";
	optMap["FILTER_HADRONS"].help += "nonzero write only events having a hadron with matching ID. For example\n";
	optMap["FILTER_HADRONS"].help += " -FILTER_HADRONS=2212 for protons";
	optMap["FILTER_HADRONS"].name = "If set to 1 (or >1), do not write output if there are no (matching) hadrons in the detectors";
	optMap["FILTER_HADRONS"].type = 0;
	optMap["FILTER_HADRONS"].ctgr = "output";

	// No output if no high mom hit
	optMap["FILTER_HIGHMOM"].arg = 0;
	optMap["FILTER_HIGHMOM"].help = "If set to non-0, do not write events if there are no high mom hit. Otherwise if \n";
	optMap["FILTER_HIGHMOM"].help += "nonzero write only events having a hit with mom > FILTER_HIGHMOM. For example\n";
	optMap["FILTER_HIGHMOM"].help += " -FILTER_HIGHMOM=1 for mom > 1MeV";
	optMap["FILTER_HIGHMOM"].name = "If set to non-0 (or >1), do not write output if there are no  high mom hit";
	optMap["FILTER_HIGHMOM"].type = 0;
	optMap["FILTER_HIGHMOM"].ctgr = "output";
	
	// sampling time of electronics (typically FADC), and number of sampling / event
	// the VT output is sampled every TSAMPLING nanoseconds to produce a ADC
	// the default number of samples is 500 ADC points, at 4ns intervals (total electronic event time = 2 microseconds)
	optMap["TSAMPLING"].args = "4, 250";
	optMap["TSAMPLING"].help = "Sampling time of electronics (typically FADC)";
	optMap["TSAMPLING"].name = "Sampling time of electronics (typically FADC)";
	optMap["TSAMPLING"].type = 1;
	optMap["TSAMPLING"].ctgr = "output";

	// Activates RNG saving for selected events
	optMap["SAVE_SELECTED"].args  = "";
	optMap["SAVE_SELECTED"].help  = "Save events with selected hit types\n";
	optMap["SAVE_SELECTED"].help  += "  arg is list of id, pid, low limit, high limit, variable[, directory]\n";
	optMap["SAVE_SELECTED"].help  += "  e.g. 7xx10000, 11, 0.0*MeV, 2000*MeV, trackE, /.\n";
	optMap["SAVE_SELECTED"].name  = "Save events with selected hit types";
	optMap["SAVE_SELECTED"].type  = 1;
	optMap["SAVE_SELECTED"].ctgr  = "output";
	
	// Reruns saved events
	optMap["RERUN_SELECTED"].args  = "";
	optMap["RERUN_SELECTED"].help  = "Rerun saved events";
	optMap["RERUN_SELECTED"].help  += "  arg is list of run #[, directory]\n";
	optMap["RERUN_SELECTED"].name  = "Rerun saved events";
	optMap["RERUN_SELECTED"].type  = 1;
	optMap["RERUN_SELECTED"].ctgr  = "control";
	






	// Physics
	// -------
	optMap["PHYSICS"].args = "STD + FTFP_BERT";
	optMap["PHYSICS"].help =  "  Physics List. The list is modular. \n\n";
	optMap["PHYSICS"].help +=  "     For example, 'STD' would only activate the standard electromagnetic processes, while\n";
	optMap["PHYSICS"].help +=  "     'QGSC_BERT + STD + HP + Optical' would also activate high precision hadronic and optical physics.\n\n";
	optMap["PHYSICS"].help +=  "     For a complete list, see http://geant4.cern.ch/support/physicsLists/referencePL/referencePL.shtml'\n\n";
	optMap["PHYSICS"].help +=  "     Available modules:\n\n";
	optMap["PHYSICS"].help +=  "      Hadronic: \n";
	optMap["PHYSICS"].help +=  "      The FTF model is based on the FRITIOF description of string excitation and fragmentation.\n";
	optMap["PHYSICS"].help +=  "      - FTFP_BERT\n";
	optMap["PHYSICS"].help +=  "      - FTFP_BERT_HP\n";
	optMap["PHYSICS"].help +=  "      - FTFP_BERT_TRV\n";
	optMap["PHYSICS"].help +=  "      - FTF_BIC\n\n";
	optMap["PHYSICS"].help +=  "      QGSP is the basic physics list applying the quark gluon string model for high ";
	optMap["PHYSICS"].help +=  " energy interactions of protons, neutrons, pions, and Kaons and nuclei.\n";
	optMap["PHYSICS"].help +=  "      - QGSC_BERT\n";
	optMap["PHYSICS"].help +=  "      - QGSP_BERT_HP\n";
	optMap["PHYSICS"].help +=  "      - QGSP_BIC\n";
	optMap["PHYSICS"].help +=  "      - QGSP_BIC_HP\n";
	optMap["PHYSICS"].help +=  "      - QGSP_FTFP_BERT\n\n";
	optMap["PHYSICS"].help +=  "      Electromagnetic: \n";
	optMap["PHYSICS"].help +=  "      - STD\n";
	optMap["PHYSICS"].help +=  "      - EMV: Designed for HEP productions. Production thresholds are enabled on secondary particles for all EM processes.\n";
	optMap["PHYSICS"].help +=  "      - EMX: Designed for HEP productions. Alternative G4Generator2BS angular generator is used for the bremsstrahlung process.\n";
	optMap["PHYSICS"].help +=  "      - EMY: Higher accuracy of electrons, hadrons and ion tracking without magnetic field\n";
	optMap["PHYSICS"].help +=  "      - EMZ: Higher accuracy of electrons, hadrons and ion tracking. Use the most accurate standard and low-energy models.\n";
	optMap["PHYSICS"].help +=  "      - LIV: Higher accuracy of electrons, hadrons and ion tracking without magnetic field. Livermore model on top of emstandard_opt4.\n";
	optMap["PHYSICS"].help +=  "      - PEN: Higher accuracy of electrons, hadrons and ion tracking without magnetic field. Penelope model on top of emstandard_opt4.\n";
	optMap["PHYSICS"].name = "Choice of Physics List";
	optMap["PHYSICS"].type = 1;
	optMap["PHYSICS"].ctgr = "physics";

	optMap["SYNRAD"].arg  = 0;
	optMap["SYNRAD"].help = " Adds Synchrotron Radiation";
	optMap["SYNRAD"].help +=  "     Possible values:\n";
	optMap["SYNRAD"].help +=  "     0: no Synchrotron Radiation (default)\n";
	optMap["SYNRAD"].help +=  "     1: vacuum Synchrotron Radiation \n";
	optMap["SYNRAD"].help +=  "     2: material Synchrotron Radiation \n";
	optMap["SYNRAD"].name = " Adds Synchrotron radiation";
	optMap["SYNRAD"].type = 0;
	optMap["SYNRAD"].ctgr = "physics";

	optMap["FORCE_MUON_RADIATIVE_DECAY"].arg = 0;
	optMap["FORCE_MUON_RADIATIVE_DECAY"].help = "Force muon radiative decay";
	optMap["FORCE_MUON_RADIATIVE_DECAY"].name = "Muon rad decay BR 100%";
	optMap["FORCE_MUON_RADIATIVE_DECAY"].type = 0;
	optMap["FORCE_MUON_RADIATIVE_DECAY"].ctgr = "physics";



	// materials
	optMap["HALL_MATERIAL"].args = "G4_AIR";
	optMap["HALL_MATERIAL"].help = "Composition of the Experimental Hall. \n";
	optMap["HALL_MATERIAL"].help += "            Air normal simulation\n";
	optMap["HALL_MATERIAL"].help += "            Air_Opt Simulation with Optical Physics\n";
	optMap["HALL_MATERIAL"].help += "            Vacuum (default)\n";
	optMap["HALL_MATERIAL"].name = "Composition of the Experimental Hall";
	optMap["HALL_MATERIAL"].type = 1;
	optMap["HALL_MATERIAL"].ctgr = "materials";
	
	optMap["DEFAULT_MATERIAL"].args = "none";
	optMap["DEFAULT_MATERIAL"].help = "Default material for missing material field.\n";
	optMap["DEFAULT_MATERIAL"].name = "Default material for missing material field";
	optMap["DEFAULT_MATERIAL"].type = 1;
	optMap["DEFAULT_MATERIAL"].ctgr = "materials";

	optMap["SWITCH_MATERIALTO"].args = "none";
	optMap["SWITCH_MATERIALTO"].help = "Switches a material for another. Usage: -SWITCH_MATERIALTO=\"old, new\" \n";
	optMap["SWITCH_MATERIALTO"].name = "Switches a material for another. ";
	optMap["SWITCH_MATERIALTO"].type = 1;
	optMap["SWITCH_MATERIALTO"].ctgr = "materials";
	optMap["SWITCH_MATERIALTO"].repe  = 1;

	optMap["CHANGEVOLUMEMATERIALTO"].args = "none";
	optMap["CHANGEVOLUMEMATERIALTO"].help = "Changes a volume material. Usage: -CHANGEVOLUMEMATERIALTO=\"volume, new material\" \n";
	optMap["CHANGEVOLUMEMATERIALTO"].name = "Changes a volume material. ";
	optMap["CHANGEVOLUMEMATERIALTO"].type = 1;
	optMap["CHANGEVOLUMEMATERIALTO"].ctgr = "materials";
	optMap["CHANGEVOLUMEMATERIALTO"].repe  = 1;


	// EM Fields
	optMap["HALL_FIELD"].args = "no";
	optMap["HALL_FIELD"].help = "Magnetic Field of the Hall. \n";
	optMap["HALL_FIELD"].name = "Magnetic Field of the Hall";
	optMap["HALL_FIELD"].type = 1;
	optMap["HALL_FIELD"].ctgr = "fields";

	optMap["FIELD_DIR"].args = "env";
	optMap["FIELD_DIR"].help = "Magnetic Field Maps Location. \n";
	optMap["FIELD_DIR"].name = "Magnetic Field Maps Location";
	optMap["FIELD_DIR"].type = 1;
	optMap["FIELD_DIR"].ctgr = "fields";
	
	optMap["NO_FIELD"].args = "none";
	optMap["NO_FIELD"].help = "Sets Magnetic Field of a volume to zero. \"all\" means no magnetic field at all. \n";
	optMap["NO_FIELD"].name = "Sets Magnetic Field of a volume to zero. \"all\" means no magnetic field at all ";
	optMap["NO_FIELD"].type = 1;
	optMap["NO_FIELD"].ctgr = "fields";
	optMap["NO_FIELD"].repe  = 0;

	optMap["G4FIELDCACHESIZE"].arg  = 3;
	optMap["G4FIELDCACHESIZE"].help = "Sets Geant4 Magnetic Field Cache Size (in mm) \n";
	optMap["G4FIELDCACHESIZE"].name = "Sets Geant4 Magnetic Field Cache Size (in mm). ";
	optMap["G4FIELDCACHESIZE"].type = 1;
	optMap["G4FIELDCACHESIZE"].ctgr = "fields";
	optMap["G4FIELDCACHESIZE"].repe  = 0;

	optMap["PHYS_VERBOSITY"].arg = 0;
	optMap["PHYS_VERBOSITY"].help = "Physics List Verbosity";
	optMap["PHYS_VERBOSITY"].name = "Physics List Verbosity";
	optMap["PHYS_VERBOSITY"].type = 0;
	optMap["PHYS_VERBOSITY"].ctgr = "fields";
	
	// by default set max field step in geant4 is 100m
	optMap["MAX_FIELD_STEP"].arg =  0;
	optMap["MAX_FIELD_STEP"].help = "Sets Maximum Acceptable Step in Magnetic Field (in mm).\n";
	optMap["MAX_FIELD_STEP"].name = "Sets Maximum Acceptable Step in Magnetic Field (in mm) ";
	optMap["MAX_FIELD_STEP"].type = 0;
	optMap["MAX_FIELD_STEP"].ctgr = "fields";
	
	optMap["SCALE_FIELD"].args  = "no, 1";
	optMap["SCALE_FIELD"].help  = "Scales Magnetic Field by a factor.\n";
	optMap["SCALE_FIELD"].help += "      Usage:\n";
	optMap["SCALE_FIELD"].help += "      -SCALE_FIELD=\"fieldname, scalefactor\"\n";
	optMap["SCALE_FIELD"].help += "      Example: -SCALE_FIELD=\"srr-solenoid, 0.5\"\n";
	optMap["SCALE_FIELD"].name  = "Electro-Magnetic Field scale";
	optMap["SCALE_FIELD"].type  = 1;
	optMap["SCALE_FIELD"].ctgr  = "fields";
	optMap["SCALE_FIELD"].repe  = 1;
	optMap["SCALE_FIELD"].argsJSONDescription = "filename, scale";
	optMap["SCALE_FIELD"].argsJSONTypes  = "S, F";

	optMap["DISPLACE_FIELDMAP"].args  = "no";
	optMap["DISPLACE_FIELDMAP"].help  = "Displace Magnetic Field map from the lab origin.\n";
	optMap["DISPLACE_FIELDMAP"].help += "      Usage:\n";
	optMap["DISPLACE_FIELDMAP"].help += "      -DISPLACE_FIELDMAP=\"fieldname, dx*unit, dy*unit, dz*unit\"\n";
	optMap["DISPLACE_FIELDMAP"].help += "      Example: -DISPLACE_FIELDMAP=\"srr-solenoid, 3.5*mm, 0*mm, 0*mm\"\n";
	optMap["DISPLACE_FIELDMAP"].name  = "Displace Magnetic Field map from the lab origin.";
	optMap["DISPLACE_FIELDMAP"].type  = 1;
	optMap["DISPLACE_FIELDMAP"].ctgr  = "fields";
	optMap["DISPLACE_FIELDMAP"].repe  = 1;
	
	optMap["ROTATE_FIELDMAP"].args  = "no";
	optMap["ROTATE_FIELDMAP"].help  = "Rotate Magnetic Field map.\n";
	optMap["ROTATE_FIELDMAP"].help += "      Usage:\n";
	optMap["ROTATE_FIELDMAP"].help += "      -ROTATE_FIELDMAP=\"fieldname, alpha*unit, beta*unit, gamma*unit\"\n";
	optMap["ROTATE_FIELDMAP"].help += "      Example: -ROTATE_FIELDMAP=\"srr-solenoid, 0*deg, 0*deg, 2*deg\"\n";
	optMap["ROTATE_FIELDMAP"].name  = "Rotate Magnetic Field map.";
	optMap["ROTATE_FIELDMAP"].type  = 1;
	optMap["ROTATE_FIELDMAP"].ctgr  = "fields";
	optMap["ROTATE_FIELDMAP"].repe  = 1;
	
	optMap["ACTIVEFIELDS"].args  = "none";
	optMap["ACTIVEFIELDS"].help  = "List of activated fields\n";
	optMap["ACTIVEFIELDS"].name  = "List of activated fields";
	optMap["ACTIVEFIELDS"].type  = 1;
	optMap["ACTIVEFIELDS"].ctgr  = "fields";
	optMap["ACTIVEFIELDS"].repe  = 0;
	optMap["ACTIVEFIELDS"].argsJSONDescription  = "ACTIVEFIELDS";
	optMap["ACTIVEFIELDS"].argsJSONTypes        = "VS";

	optMap["FIELD_PROPERTIES"].args  = "no";
	optMap["FIELD_PROPERTIES"].help  = "Mapped field minimum step, integration method, interpolation\n\n";
	optMap["FIELD_PROPERTIES"].help += "      Usage:\n";
	optMap["FIELD_PROPERTIES"].help += "      -FIELD_PROPERTIES=\"fieldname, minStep, integralAlgorithm, (interpolationMethod)\"\n\n";
	optMap["FIELD_PROPERTIES"].help += "      Example: -FIELD_PROPERTIES=\"srr-solenoid, 1*mm, G4ClassicalRK4, linear\"\n\n";
	optMap["FIELD_PROPERTIES"].help += "      Available Integration Methods:\n";
	optMap["FIELD_PROPERTIES"].help += "       - G4CashKarpRKF45: Fift Order Range Kutta, for very smooth fields\n";
	optMap["FIELD_PROPERTIES"].help += "       - G4ClassicalRK4: Fourth Order Range Kutta. Robust for every field.\n";
	optMap["FIELD_PROPERTIES"].help += "       - G4SimpleHeum: Third order stepper.\n";
	optMap["FIELD_PROPERTIES"].help += "       - G4SimpleRunge: Simplified (second order) Range Kutta (faster).\n";
	optMap["FIELD_PROPERTIES"].help += "       - G4ImplicitEuler: Second order stepper, for faster varying fields.\n";
	optMap["FIELD_PROPERTIES"].help += "       - G4ExplicitEuler: First order stepper, for rough fields.\n";
	optMap["FIELD_PROPERTIES"].help += "       - G4HelixImplicitEuler: Second order, specialized for helix-like trajectories.\n";
	optMap["FIELD_PROPERTIES"].help += "       - G4HelixExplicitEuler: First order, specialized for helix-like trajectories.\n";
	optMap["FIELD_PROPERTIES"].help += "       - G4HelixSimpleRunge: Second order Range Kutta, specialized for helix-like trajectories.\n";
	optMap["FIELD_PROPERTIES"].help += "       - G4NystromRK4: provides accuracy near that of G4ClassicalRK4 with a significantly reduced cost in field evaluation.\n\n";
	optMap["FIELD_PROPERTIES"].help += "       Available Interpolation Methods:\n";
	optMap["FIELD_PROPERTIES"].help += "       - none: closest grid point.\n";
	optMap["FIELD_PROPERTIES"].help += "       - linear: linear interpolation.\n\n";
	optMap["FIELD_PROPERTIES"].help += "       Note: specifying interpolation method is optional. \"linear\" is the default.\n";
	optMap["FIELD_PROPERTIES"].name  = "Mapped field minimum step, integration method, interpolation";
	optMap["FIELD_PROPERTIES"].type  = 1;
	optMap["FIELD_PROPERTIES"].ctgr  = "fields";
	optMap["FIELD_PROPERTIES"].repe  = 1;

	optMap["FIELDS_FILENAMES"].args  = "none";
	optMap["FIELDS_FILENAMES"].help  = "List of activated fields\n";
	optMap["FIELDS_FILENAMES"].name  = "List of activated fields";
	optMap["FIELDS_FILENAMES"].type  = 1;
	optMap["FIELDS_FILENAMES"].ctgr  = "fields";
	optMap["FIELDS_FILENAMES"].repe  = 0;
	optMap["FIELDS_FILENAMES"].argsJSONDescription  = "label filename";
	optMap["FIELDS_FILENAMES"].argsJSONTypes        = "S S";

	// transportation
	optMap["REMOVESENSITIVITY"].args = "none";
	optMap["REMOVESENSITIVITY"].help = "Removes a volumes sensitivity. Usage: -REMOVESENSITIVITY=volumename1, volumename2, ... \n";
	optMap["REMOVESENSITIVITY"].name = "Changes a volumes sensitivity. ";
	optMap["REMOVESENSITIVITY"].type = 1;
	optMap["REMOVESENSITIVITY"].ctgr = "transportation";

	optMap["PRODUCTIONCUT"].arg   = 10;
	optMap["PRODUCTIONCUT"].help  = "Production cut for root, in mm\n";
	optMap["PRODUCTIONCUT"].name  = "Production cut for root, in mm";
	optMap["PRODUCTIONCUT"].type  = 0;
	optMap["PRODUCTIONCUT"].ctgr  = "transportation";

	optMap["PRODUCTIONCUTFORVOLUMES"].args   = "no";
	optMap["PRODUCTIONCUTFORVOLUMES"].help  = "Production cut for volumes (separated by commas), in mm\n";
	optMap["PRODUCTIONCUTFORVOLUMES"].help += "      Example: \"outerMount, taggerInnerShield, 50\" will set the production cut for the volumes outerMount and taggerInnerShield to 5cm\n";
	optMap["PRODUCTIONCUTFORVOLUMES"].name  = "Production cut for volumes, in mm";
	optMap["PRODUCTIONCUTFORVOLUMES"].type  = 1;
	optMap["PRODUCTIONCUTFORVOLUMES"].ctgr  = "transportation";
	optMap["PRODUCTIONCUTFORVOLUMES"].repe  = 1;

	// voltage versus time resolution, in ns
	optMap["FASTMCMODE"].arg   = 0;   // 1: production cut set to 5m. 2: turn off all physics but transportation
	optMap["FASTMCMODE"].help  = "Activate FASTMC Mode. Possible values: \n\n";
	optMap["FASTMCMODE"].help += "       0: No FASTMC (default).\n";
	optMap["FASTMCMODE"].help += "       1: Disable secondaries, disable hit processes.\n";
	optMap["FASTMCMODE"].help += "       2: Disable all physics but transportation, disable hit process.\n";
	optMap["FASTMCMODE"].help += "       10: Disable secondaries, hit processes are enabled.\n";
	optMap["FASTMCMODE"].help += "       20: Disable all physics but transportation, hit processes are enabled.\n";
	optMap["FASTMCMODE"].name  = "Activate FASTMC Mode if non zero.";
	optMap["FASTMCMODE"].type  = 0;
	optMap["FASTMCMODE"].ctgr  = "transportation";

	optMap["HALL_DIMENSIONS"].args = "20*m, 20*m, 20*m";
	optMap["HALL_DIMENSIONS"].help = "(x,y,z) semi-dimensions of the experimental Hall.\n";
	optMap["HALL_DIMENSIONS"].name = "(x,y,z) semi-dimensions of the experimental Hall.";
	optMap["HALL_DIMENSIONS"].type = 1;
	optMap["HALL_DIMENSIONS"].ctgr = "transportation";



	// General
	optMap["DF"].args  = "no";
	optMap["DF"].help  = "Selects Detector System and Factory type. \n";
	optMap["DF"].help += "      Example:\n";
	optMap["DF"].help += "        -DF=\"CTOF, MYSQL\"  selects the MYSQL factory for the detector CTOF\n";
	optMap["DF"].name  = "Detector System and Factory type.";
	optMap["DF"].type  = 1;
	optMap["DF"].ctgr  = "general";
	optMap["DF"].repe  = 1;
    
    //
    //
    //
    //
    //
    //
    //
    //
    //
    // Dark Matter
    // --------------
    
    /* Vanilla Dark Photon model */
    optMap["DARK_PHOTON"].args = "no";
    optMap["DARK_PHOTON"].help = "Dark Photon parameters: mass, spin*2, parity\n";
    optMap["DARK_PHOTON"].help += "example: -DARK_PHOTON=\"200*MeV,2,-1\"";
    optMap["DARK_PHOTON"].type = 1;
    optMap["DARK_PHOTON"].ctgr = "dark_photon";

    optMap["DARK_MATTER"].args = "no";
    optMap["DARK_MATTER"].help = "Dark Matter parameters: mass, spin*2\n";
    optMap["DARK_MATTER"].help += "example: -DARK_MATTER=\"200*MeV,1\"";
    optMap["DARK_MATTER"].type = 1;
    optMap["DARK_MATTER"].ctgr = "dark_photon";

    optMap["DARK_COUPLINGS"].args = "no";
    optMap["DARK_COUPLINGS"].help = "Dark photon couplings: eps, alphaD, min_kin_energy\n";
    optMap["DARK_COUPLINGS"].help += "example: -DARK_COUPLINGS=\"0.001,0.5,0.1*GeV\"";
    optMap["DARK_COUPLINGS"].type = 1;
    optMap["DARK_COUPLINGS"].ctgr = "dark_photon";
    /* end of Dark Photon model */
    
    /* Dark Scalar model */
    optMap["DARK_SCALAR"].args = "no";
    optMap["DARK_SCALAR"].help = "Dark Scalar parameters: mass, coupling\n";
    optMap["DARK_SCALAR"].help += "example: -DARK_PHOTON=\"100*MeV,3.87e-4\"";
    optMap["DARK_SCALAR"].type = 1;
    optMap["DARK_SCALAR"].ctgr = "dark_scalar";
    
    optMap["DARK_SCALAR_BIAS"].arg = 0;
    optMap["DARK_SCALAR_BIAS"].help = "Flat increase Dark Scalar production rate\n";
    optMap["DARK_SCALAR_BIAS"].name = "Set bias for dark scalar production (e.g. 1e7)";
    optMap["DARK_SCALAR_BIAS"].type = 0;
    optMap["DARK_SCALAR_BIAS"].ctgr = "dark_scalar_bias";
    /* end of Dark Scalar*/
}
