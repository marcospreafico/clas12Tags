 // G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "scint_hitprocess.h"

// Root utils
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <fstream>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> scint_HitProcess::integrateDgt(MHit* aHit, int hitn) {
    
    map<string, double> dgtz;
    vector<identifier> identity = aHit->GetId();
    
    int sector = identity[0].id;    // module
    int channel    = identity[1].id;
    int module    = identity[2].id;

    // Parameters for plastic scintillator
    double light_yield       = 9200;            // photons / MeV
    double attenuation_length = 100 * cm;       // attenuation length in plastic
    double veff              = 13 * cm / ns;    // effective velocity of light
    double sipm_pde          = 0.5;             // SiPM photon detection efficiency

    // Dimensions of plastic scintillator
    double dx = aHit->GetDetector().dimensions[0] * 2; // full width  (x)
    double dy = aHit->GetDetector().dimensions[1] * 2; // full height (y)
    double dz = aHit->GetDetector().dimensions[2] * 2; // full length (z) - readout axis

    double z_readout = dz / 2.; // readout at one end of the scintillator

    // Extract Birks constant from material properties
    double birks_constant = aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
    
    // Get hit information
    vector<G4ThreeVector> Lpos  = aHit->GetLPos();   // local positions of steps
    vector<G4double>      Edep  = aHit->GetEdep();   // energy deposited per step
    vector<G4double>      Dx    = aHit->GetDx();     // step lengths
    vector<int>           charge = aHit->GetCharges(); // particle charge per step
    vector<G4double>      times  = aHit->GetTime();  // time per step

    unsigned int nsteps = Edep.size();

    // Sum total energy deposit
    double Etot = 0;
    for (unsigned int s = 0; s < nsteps; s++) Etot += Edep[s];

    double Etot_B  = 0;   // Birks-attenuated total energy
    double etot    = 0;   // light collected at readout after attenuation
    double time_hit = 0;  // average hit time at readout
    double ADC_scint = 0;
    double TDC     = 4096; 

    double eff_ly = 22.04; // effective LY ???? is different to light_yield by such a big factor? 

    if (Etot > 0) {
        for (unsigned int s = 0; s < nsteps; s++) {
            
            // Apply Birks attenuation
            double Edep_B = BirksAttenuation(Edep[s], Dx[s], charge[s], birks_constant);
            Etot_B += Edep_B;

            // Distance from hit position to readout along z
            double zhit = Lpos[s].z();
            double dist = fabs(zhit - z_readout);

            // Light collected: attenuated by propagation to readout
            etot     += Edep_B * exp(-dist / attenuation_length);

            // Time: hit time + propagation time to readout
            time_hit += (times[s] + dist / veff) / nsteps;
        }

        // Convert energy to photoelectrons
        double npe = etot * eff_ly;  // number of detected photoelectrons 

        // Apply Poisson fluctuations to the number of photoelectrons (Ask M. Spreafico if this is needed when using WF)
        npe = G4Poisson(npe);

        // Gaussian smearing for SiPM excess noise factor
        // target: sqrt(3.70^2 - 0.45^2) ≈ 3.67 phe (fixed value to match expected resolution at 22 phe)
        double sigma_sipm = 3.67;
        npe = G4RandGauss::shoot(npe, sigma_sipm);
        if (npe < 0) npe = 0;

        // ADC: energy from effective LY
        ADC_scint = npe;
        if (ADC_scint < 0) ADC_scint = 0;
        TDC = time_hit;
    }

    dgtz["hitn"]      = hitn;
    dgtz["sector"]    = sector;
    dgtz["layer"]     = channel;
    dgtz["component"] = module;
    dgtz["ADC_order"] = 0;
    dgtz["ADC_ADC"]   = 1000 * ADC_scint;     // number of photoelectrons
    dgtz["ADC_time"]  = 1000 * Etot_B;       
    dgtz["ADC_ped"]   = 0;

    return dgtz;
}


vector<identifier> scint_HitProcess::processID(vector<identifier> id, G4Step *step, detector Detector) {
	id[id.size() - 1].id_sharing = 1;
	return id;
}

double scint_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks) {
	//Example of Birk attenuation law in organic scintillators.
	//adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
	//
	// Taken from GEANT4 examples advanced/amsEcal and extended/electromagnetic/TestEm3
	//
	double response = destep;
	if (birks * destep * stepl * charge != 0.) {
		response = destep / (1. + birks * destep / stepl);
	}
	return response;
}

double scint_HitProcess::BirksAttenuation2(double destep, double stepl, int charge, double birks) {
	//Extension of Birk attenuation law proposed by Chou
	// see G.V. O'Rielly et al. Nucl. Instr and Meth A368(1996)745
	//
	double C = 9.59 * 1E-4 * mm * mm / MeV / MeV;
	double response = destep;
	if (birks * destep * stepl * charge != 0.) {
		response = destep / (1. + birks * destep / stepl + C * pow(destep / stepl, 2.));
	}
	return response;
}

map<string, vector<int> > scint_HitProcess::multiDgt(MHit* aHit, int hitn) {
	map<string, vector<int> > MH;

	return MH;
}

double* scint_HitProcess::WaveForm(double npe, double* time, string crs_material){
    double c = exp(-2);
    
    int Nsamp_WF = 2500; // number of samples to generate WF; 2500 = 10 us
    
    static double* WFsample = new double[Nsamp_WF]; // object to save the WF
    for(unsigned int s = 0; s < Nsamp_WF; s++){ WFsample[s] = 0; }
    
    double smp_t = 4./1000.; // assuming fADC sampling at 250 MHz 1 sample = 4 ns
    
    array<double, 6> p = {0., 0., 0., 0., 0., 0.}; // crs scintillation time parameters
    if(crs_material == "G4_PbWO4"){
        p = { 0., 0.00680, 0.64, 0.0334, 0.36, 0. }; // PbWO4
    }
    else if(crs_material == "CsI_Tl"){
        p = { 0., 0.680, 0.64, 3.34, 0.36, 0. }; // BaBar CsI
    } else if(crs_material == "G4_BGO"){
        p = { 0., 0.060, 0.10, 0.230, 0.90, 0. }; // PbWO4
    }
    else{
        cout << "ERROR: no waveform parametrization for this material ( "+crs_material+" )" << endl;
        return 0;
    }
    TF1* tdistrib = new TF1("tdistrib", "([2]/[1]*exp(-x/[1])+[4]/[3]*exp(-x/[3]))/([2]/[1]+[4]/[3])", 0, Nsamp_WF*smp_t);
    for(int ii = 0; ii < 6; ii ++){ tdistrib->SetParameter(ii, p[ii]); }
    
    
    //Definition of phe shape
    double tau = 15.; // ampli response time constant (in ns)
    double t0 = 0.01; // t0 starting time (in ns)
    double area = (tau / c / 2.);
    double A = 1. / area; // amplitude at mnax (55.41 to have it normalized to integral=1, otherwise the max is at 1)
    
    double t_spread = 1. * 0.000; // pream time spread in us
    double A_spread = 1. * 0.4 * A; // pream amp spread (in fraction of 1pe amplitude = A)
    if(crs_material == "CsI_Tl"){
         t_spread = 1. * 0.020; // pream time spread in us
         A_spread = 1. * 0.05 * A; // pream amp spread (in fraction of 1pe amplitude = A)
    }
    
    static double AmpWF[80]; // phe WF
    for (unsigned int s = 0; s < 80; s++) {
        double t = 1000. * s * smp_t;
        double func = (t - t0) * (t - t0) * exp(-(t - t0) / tau) * A / (4 * tau * tau * c) * 0.5 * (abs(t - t0) / (t - t0) + 1);
        AmpWF[s] = smp_t * 1000. * func;
    }
    
    static double frac = 1 - ((p[2] * exp(-smp_t * Nsamp_WF / p[1]) + p[4] * exp(-smp_t * Nsamp_WF / p[3])));    // fraction of pe in Nsamp_WF
    
    
    // generate waveform sample
    double t; int it; // time variable
    int mNpe = G4Poisson(frac * npe); // number of phe in Nsamp_WF
    for(unsigned int s = 1; s <= mNpe; s++){
        t = tdistrib->GetRandom();
        t = G4RandGauss::shoot(t, t_spread);
        if(t < 0.) t = 0.;
        if(t > smp_t * Nsamp_WF) t = smp_t * Nsamp_WF;
        it = t / smp_t;
        for(unsigned int s = 0; s < 80; s++){ // sum the phe WF at the phe time
            double func = G4RandGauss::shoot(AmpWF[s], A_spread);
            if((it + s) > Nsamp_WF) break;
            WFsample[it + s] += func;
        }
    }
    
    // mimicking a CF discriminator at 1/3 of the max signal
    *time = 0.;
    double time_max = -100;
    int s = 0, s_time_max = 0;
    while(time_max < WFsample[s]){// search for max WF time
        time_max = 1/2. * (WFsample[s + 1] + WFsample[s]);
        s_time_max = s;
        *time = 1000. * smp_t * s_time_max / 3.;
        s++;
    }
    
    return WFsample;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> scint_HitProcess::electronicNoise() {
	vector<MHit*> noiseHits;

	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)

	return noiseHits;
}

// - charge: returns charge/time digitized information / step
map<int, vector<double> > scint_HitProcess::chargeTime(MHit* aHit, int hitn) {
	map<int, vector<double> > CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double scint_HitProcess::voltage(double charge, double time, double forTime) {
	return 0.0;
}

