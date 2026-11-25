#include "BdxCrystalOptics.h"
#include "BdxOpticalSchema.h"
#include "hipo4/event.h"
#include "hipo4/bank.h"
#include "G4VUserTrackInformation.hh"
#include <vector>
#include "G4Event.hh"
#include "G4Box.hh"
#include "G4GenericTrap.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicalConstants.hh"
#include "G4ParticleChange.hh"
#include "G4PrimaryVertex.hh"
#include "G4ProcessManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4VProcess.hh"
#include "G4IonisParamMat.hh"
#include <algorithm>
#include "G4VPhysicalVolume.hh"
#include "G4Version.hh"
#include "G4Threading.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4TouchableHistory.hh"
#include <array>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>
#include <sstream>

namespace bdxoptics {
namespace {
struct Config {
    double light_yield_per_MeV=8500, refractive_index=2.15;
    double absorption_length_mm=10000, decay_time_ns=300, wavelength_nm=480;
    double wrap_reflectivity=0.98;
    double diffuse_length_mm=40, transition_length_mm=100, specular_length_mm=100;
    double boundary_push_mm=1.e-6;
    std::string wrap_model="paint";
    bool cherenkov=false;
};
Config cfg;
G4VPhysicalVolume *crystalPV=nullptr, *sensorPV=nullptr, *worldPV=nullptr;
double sensorZ=0, sensorY=0;
bool fullBase=false;
std::array<double,7> dimensions;
std::ofstream events, arrivals;
G4MaterialPropertyVector* specularProbability=nullptr;
struct PhotonInfo final : G4VUserTrackInformation {
    int source=-1, reflections=0;
    G4ThreeVector birth;
    double birthTime=0;
    PhotonInfo(int sourceStep,const G4Track& t): source(sourceStep),birth(t.GetPosition()),birthTime(t.GetGlobalTime()) {}
    void Print() const override {}
};
struct Deposit {
    int track,parent,pid,step,scint,cerenkov;
    G4ThreeVector start,end;
    double t0,t1,energy;
};
struct Arrival {
    int track,parent,origin,source,steps,reflections;
    G4ThreeVector position,direction,birth;
    double time,energy,path,birthTime,weight;
};
struct Event {
    long generated=0, generatedC=0, injected=0;
    long arrived=0, arrivedC=0, arrivedI=0;
    long bulk=0, surface=0, escaped=0, reflections=0, diffuse=0, specular=0, maxSteps=0, boundaryRepairs=0, birthRepairs=0;
    double maxBoundaryShift=0, maxBirthShift=0;
    double edep=0, weightedZ=0, weightedD=0, primaryZ=0;
    int index=-1;
    std::vector<Deposit> deposits;
    std::vector<Arrival> flux;
};
thread_local Event ev;
std::string Env(const char* key) {
    const auto* value=std::getenv(key);
    if(!value || !*value) throw std::runtime_error(std::string("BDX optics: missing ")+key);
    return value;
}
void Fail(const std::string& message) { throw std::runtime_error("BDX optics: "+message); }
std::string Trim(const std::string& s) {
    const auto a=s.find_first_not_of(" \t\r\n");
    return a==std::string::npos ? "" : s.substr(a,s.find_last_not_of(" \t\r\n")-a+1);
}
Config ReadConfig() {
    Config c;
    std::ifstream input(Env("BDX_OPTICS_CONFIG"));
    if(!input) Fail("cannot open optical configuration");
    std::map<std::string,double*> numbers={
        {"light_yield_per_MeV",&c.light_yield_per_MeV}, {"refractive_index",&c.refractive_index},
        {"absorption_length_mm",&c.absorption_length_mm}, {"decay_time_ns",&c.decay_time_ns},
        {"wavelength_nm",&c.wavelength_nm}, {"wrap_reflectivity",&c.wrap_reflectivity},
        {"diffuse_length_mm",&c.diffuse_length_mm}, {"transition_length_mm",&c.transition_length_mm},
        {"specular_length_mm",&c.specular_length_mm}, {"boundary_push_mm",&c.boundary_push_mm}};
    std::map<std::string,bool> seen;
    std::string line;
    while(std::getline(input,line)) {
        line=Trim(line.substr(0,line.find('#')));
        if(line.empty()) continue;
        const auto p=line.find('=');
        if(p==std::string::npos) Fail("expected key=value: "+line);
        const auto k=Trim(line.substr(0,p)), v=Trim(line.substr(p+1));
        if(seen[k]) Fail("duplicate key "+k);
        seen[k]=true;
        if(k=="wrap_model") c.wrap_model=v;
        else if(k=="cherenkov") {
            if(v!="0" && v!="1") Fail("cherenkov must be 0 or 1");
            c.cherenkov=v=="1";
        } else {
            const auto it=numbers.find(k);
            if(it==numbers.end()) Fail("unknown key "+k);
            std::size_t n=0; const auto x=std::stod(v,&n);
            if(n!=v.size() || !std::isfinite(x)) Fail("invalid number for "+k);
            *it->second=x;
        }
    }
    if(c.light_yield_per_MeV<=0 || c.refractive_index<1 || c.absorption_length_mm<0 ||
       c.decay_time_ns<=0 || c.wavelength_nm<=250 || c.wavelength_nm>=800 ||
       c.wrap_reflectivity<0 || c.wrap_reflectivity>1 || c.diffuse_length_mm<0 ||
       c.transition_length_mm<0 || c.specular_length_mm<0 || c.boundary_push_mm<0 || c.boundary_push_mm>1.e-3) Fail("configuration value out of range");
    if(c.wrap_model!="paint" && c.wrap_model!="contact" && c.wrap_model!="bare")
        Fail("wrap_model must be paint, contact or bare");
    if(c.wrap_reflectivity==1 && c.absorption_length_mm==0 && c.wrap_model!="bare")
        Fail("perfect lossless reflector can trap rays indefinitely; use finite absorption");
    return c;
}
void Constant(G4MaterialPropertiesTable* m, const char* key, double x) {
    double e[2]={h_Planck*c_light/(800*nm),h_Planck*c_light/(250*nm)}, v[2]={x,x};
    m->AddProperty(key,e,v,2);
}
G4Material* Material(G4LogicalVolume* lv,const std::string& name,double index,bool scintillating) {
    auto* base=lv->GetMaterial();
    auto* m=G4Material::GetMaterial(name,false);
    if(!m) m=new G4Material(name,base->GetDensity(),base);
    auto* table=new G4MaterialPropertiesTable;
    Constant(table,"RINDEX",index);
    if(scintillating) {
        if(cfg.absorption_length_mm>0) Constant(table,"ABSLENGTH",cfg.absorption_length_mm*mm);
        // Narrow peak proxy. All transport coefficients are wavelength independent
        // in the geometry study, so integrated collection does not depend on it.
        double e[3]={h_Planck*c_light/((cfg.wavelength_nm+0.01)*nm),
                     h_Planck*c_light/(cfg.wavelength_nm*nm),
                     h_Planck*c_light/((cfg.wavelength_nm-0.01)*nm)};
        double intensity[3]={0,1,0};
#if G4VERSION_NUMBER < 1100
        table->AddProperty("FASTCOMPONENT",e,intensity,3);
        table->AddConstProperty("FASTTIMECONSTANT",cfg.decay_time_ns*ns);
        table->AddConstProperty("YIELDRATIO",1.0);
#else
        table->AddProperty("SCINTILLATIONCOMPONENT1",e,intensity,3);
        table->AddConstProperty("SCINTILLATIONTIMECONSTANT1",cfg.decay_time_ns*ns);
        table->AddConstProperty("SCINTILLATIONYIELD1",1.0);
#endif
        table->AddConstProperty("SCINTILLATIONYIELD",cfg.light_yield_per_MeV/MeV);
        table->AddConstProperty("RESOLUTIONSCALE",1.0);
        // No plastic-scintillator Birks coefficient; BGO proportional reference.
        m->GetIonisation()->SetBirksConstant(0.0);
    }
    m->SetMaterialPropertiesTable(table);
    lv->SetMaterial(m);
    return m;
}
// Preserve GenericTrap's exact ruled geometry, using stable polynomial roots
// and conservative scalar safety bounds on its twisted faces.
// A safety must be a LOWER bound on distance to a boundary; using a tangent
// plane anchored at another z can otherwise let multiple scattering leave BGO.
class ConservativeBgoTrap final : public G4GenericTrap {
    double length,h0,h1,a0,a1,c0,c1,xmax,ymax,sideScale,yScale;
    std::array<double,6> Margins(const G4ThreeVector& p) const {
        const double t=(p.z()+length)/(2*length);
        const double h=h0+t*(h1-h0),a=a0+t*(a1-a0),c=c0+t*(c1-c0);
        return {{length-p.z(),length+p.z(),(h-p.y())/yScale,(h+p.y())/yScale,
                 (h*(a-p.x())+c*p.y())/sideScale,(h*(a+p.x())+c*p.y())/sideScale}};
    }
    using Polynomial=std::array<double,3>;
    std::array<Polynomial,6> RayPolynomials(const G4ThreeVector& p,const G4ThreeVector& v) const {
        const double t=(static_cast<double>(p.z())+length)/(2*length);
        const double h=h0+t*(h1-h0),a=a0+t*(a1-a0),c=c0+t*(c1-c0);
        const double dh=(h1-h0)*static_cast<double>(v.z())/(2*length),
                          da=(a1-a0)*static_cast<double>(v.z())/(2*length),
                          dc=(c1-c0)*static_cast<double>(v.z())/(2*length);
        const auto side=[&](int sign) -> Polynomial {
            const double x=a+sign*static_cast<double>(p.x()),dx=da+sign*static_cast<double>(v.x());
            return {{(h*x+c*p.y())/sideScale,(h*dx+dh*x+c*v.y()+dc*p.y())/sideScale,
                     (dh*dx+dc*v.y())/sideScale}};
        };
        return {{{{length-p.z(),-v.z(),0}},{{length+p.z(),v.z(),0}},
                 {{(h-p.y())/yScale,(dh-v.y())/yScale,0}},
                 {{(h+p.y())/yScale,(dh+v.y())/yScale,0}},side(-1),side(1)}};
    }
    G4double Crossing(const G4ThreeVector& p,const G4ThreeVector& v,bool entering) const {
        const auto polynomials=RayPolynomials(p,v);
        std::array<double,12> roots{};
        std::size_t count=0;
        const double tolerance=1.e-9*mm;
        const auto add=[&](double root) {
            if(std::isfinite(root) && root>=-tolerance && root<kInfinity)
                roots[count++]=std::max(0.0,root);
        };
        for(const auto& f:polynomials) {
            const double c=f[0],b=f[1],a=f[2];
            if(a==0) {if(b!=0) add(-c/b);continue;}
            double discriminant=b*b-4*a*c;
            const double error=64*std::numeric_limits<double>::epsilon()*(b*b+std::abs(4*a*c));
            if(discriminant<0 && discriminant>=-error) discriminant=0;
            if(discriminant<0) continue;
            const double q=-0.5*(b+std::copysign(std::sqrt(discriminant),b));
            if(q==0) add(-b/(2*a));
            else {add(q/a);add(c/q);}
        }
        std::sort(roots.begin(),roots.begin()+count);
        for(std::size_t rootIndex=0;rootIndex<count;++rootIndex) {
            const double distance=roots[rootIndex];
            bool before=true,after=true;
            for(const auto& f:polynomials) {
                const double value=(f[2]*distance+f[1])*distance+f[0];
                if(value < -tolerance) {before=after=false;break;}
                if(value>tolerance) continue;
                const double derivative=2*f[2]*distance+f[1];
                if(derivative>0) before=false;
                else if(derivative<0) after=false;
                else if(f[2]<0) before=after=false;
            }
            if((entering && !before && after) || (!entering && before && !after))
                return static_cast<G4double>(distance);
        }
        return kInfinity;
    }
public:
    explicit ConservativeBgoTrap(const G4GenericTrap& source)
        :G4GenericTrap(source.GetName()+"_conservative_safety",source.GetZHalfLength(),source.GetVertices()) {
        const auto& v=source.GetVertices();length=source.GetZHalfLength();
        for(int k:{0,4}) {
            const double tolerance=1.e-9*mm;
            if(std::abs(v[k].x()+v[k+3].x())>tolerance || std::abs(v[k+1].x()+v[k+2].x())>tolerance ||
               std::abs(v[k].y()-v[k+3].y())>tolerance || std::abs(v[k+1].y()-v[k+2].y())>tolerance ||
               std::abs(v[k].y()+v[k+1].y())>tolerance)
                Fail("conservative BGO safety requires the documented centred trapezoid vertices");
        }
        h0=v[1].y();h1=v[5].y();
        a0=(v[3].x()+v[2].x())/2;a1=(v[7].x()+v[6].x())/2;
        c0=(v[2].x()-v[3].x())/2;c1=(v[6].x()-v[7].x())/2;
        xmax=0;ymax=0;
        for(const auto& point:v) {xmax=std::max(xmax,std::abs(point.x()));ymax=std::max(ymax,std::abs(point.y()));}
        const double hmax=std::max(h0,h1),amax=std::max(std::abs(a0),std::abs(a1));
        const double cmax=std::max(std::abs(c0),std::abs(c1));
        const double dh=(h1-h0)/(2*length),da=(a1-a0)/(2*length),dc=(c1-c0)/(2*length);
        // Global bound on the norm of grad[h(z)*(a(z) +/- x)+c(z)*y]
        // throughout the bounding box. F/|grad F|max is a conservative safety.
        const double dz=std::abs(dh)*(amax+xmax)+hmax*std::abs(da)+std::abs(dc)*ymax;
        sideScale=std::sqrt(hmax*hmax+cmax*cmax+dz*dz);
        yScale=std::sqrt(1+dh*dh);
    }
    using G4GenericTrap::DistanceToIn;
    using G4GenericTrap::DistanceToOut;
    G4double DistanceToOut(const G4ThreeVector& p) const override {
        const auto distances=Margins(p);
        return std::max(0.0,*std::min_element(distances.begin(),distances.end()));
    }
    G4double DistanceToIn(const G4ThreeVector& p) const override {
        // Outside the bounding box its distance is already a valid lower bound.
        const double dx=std::max(0.0,std::abs(p.x())-xmax),dy=std::max(0.0,std::abs(p.y())-ymax),
                     dz=std::max(0.0,std::abs(p.z())-length);
        if(dx>0 || dy>0 || dz>0) return std::sqrt(dx*dx+dy*dy+dz*dz);
        const auto distances=Margins(p);
        return std::max(0.0,-*std::min_element(distances.begin(),distances.end()));
    }
    G4double DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& v) const override {
        return Crossing(p,v,true);
    }
    G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
                         const G4bool calcNormal=false,G4bool* validNormal=nullptr,G4ThreeVector* normal=nullptr) const override {
        const double distance=Crossing(p,v,false);
        if(calcNormal) {
            if(normal && distance<kInfinity) *normal=G4GenericTrap::SurfaceNormal(p+distance*v);
            // False is conservative for curved/twisted faces and edges.
            if(validNormal) *validNormal=false;
        }
        return distance;
    }
    G4VSolid* Clone() const override { return new ConservativeBgoTrap(*this); }
};

// GEMC and these checks use the sequential G4RunManager. Update the UNIFIED
// mixture before the native boundary process samples it at the collision point.
// No coating segmentation, seams, overlapping solids or post-step direction edits.
class SpatialBoundary final : public G4OpBoundaryProcess {
public:
    G4VParticleChange* PostStepDoIt(const G4Track& track,const G4Step& step) override {
        if(cfg.wrap_model=="paint" && specularProbability &&
           step.GetPostStepPoint()->GetStepStatus()==fGeomBoundary &&
           step.GetPreStepPoint()->GetPhysicalVolume()==crystalPV &&
           step.GetPostStepPoint()->GetPhysicalVolume()==worldPV) {
            const double fraction=SpecularFraction((sensorZ-step.GetPostStepPoint()->GetPosition().z())/mm);
            specularProbability->PutValue(0,fraction);
            specularProbability->PutValue(1,fraction);
        }
        auto* result=G4OpBoundaryProcess::PostStepDoIt(track,step);
        const auto status=GetStatus();
        const bool reflected=status==LambertianReflection || status==SpikeReflection ||
            status==LobeReflection || status==BackScattering || status==FresnelReflection ||
            status==TotalInternalReflection;
        if(reflected && cfg.boundary_push_mm>0 && step.GetPreStepPoint()->GetPhysicalVolume()==crystalPV) {
            // G4GenericTrap in Geant4 10.7.4 can locate a reflected photon in the
            // world at the next finite step. Repair intersection roundoff and
            // place it strictly inside the same solid with the configured clearance;
            // no direction, energy or survival reweighting; convergence is tested.
            auto* solid=crystalPV->GetLogicalVolume()->GetSolid();
            const auto point=step.GetPostStepPoint()->GetPosition();
            const auto inside=RelocateReflectedPoint(*solid,point,cfg.boundary_push_mm*mm);
            if(solid->Inside(point)==kOutside) ++ev.boundaryRepairs;
            ev.maxBoundaryShift=std::max(ev.maxBoundaryShift,(inside-point).mag()/mm);
            auto* change=dynamic_cast<G4ParticleChange*>(result);
            if(!change) Fail("unexpected boundary particle change type");
            change->ProposePosition(inside);
        }
        return result;
    }
};

bool Transmitted(G4OpBoundaryProcessStatus status) {
    return status==FresnelRefraction || status==Transmission || status==SameMaterial;
}
}

bool Enabled() {
    const auto* c=std::getenv("BDX_OPTICS_CONFIG");
    return c && *c;
}
G4VSolid* MakeConservativeBgoSolid(const G4GenericTrap& source) { return new ConservativeBgoTrap(source); }
G4ThreeVector RelocateReflectedPoint(const G4VSolid& solid,const G4ThreeVector& point,double push) {
    // GenericTrap's nearly degenerate ray/surface intersection can land a few
    // nanometres outside the mathematical surface. A fixed 1 nm push is then
    // insufficient. Keep the normal push as the fast path, then find the
    // closest classified-inside point along a short inward segment. Do not use
    // GenericTrap's scalar DistanceToIn as a bound: for twisted faces it can
    // report micrometres for a point only picometres outside the exact section.
    const double maxRoundoffRepair=1.e-3*mm; // hard bound: never conceal a macroscopic geometry error
    if(!std::isfinite(push) || push<=0 || push>maxRoundoffRepair ||
       !std::isfinite(point.x()) || !std::isfinite(point.y()) || !std::isfinite(point.z()))
        Fail("invalid reflected-point relocation input");
    const double limit=maxRoundoffRepair+push;
    const auto normal=solid.SurfaceNormal(point);
    G4ThreeVector candidate=point;
    if(std::isfinite(normal.mag2()) && normal.mag2()>0) {
        candidate=point-normal.unit()*push;
        if((candidate-point).mag()<=limit && solid.Inside(candidate)==kInside) return candidate;
        // A normal at an edge, or a warped face classified within its tolerance
        // band, need not cross every neighbouring face. The BDX solids are
        // star-shaped about their bounding-box centre. Search towards that known
        // interior point, always within the same small numerical bound.
        G4ThreeVector minimum,maximum;
        solid.BoundingLimits(minimum,maximum);
        const auto centre=(minimum+maximum)*0.5;
        if(solid.Inside(centre)==kInside && (centre-point).mag2()>0) {
            const auto inward=(centre-point).unit();
            candidate=point+inward*limit;
            if(solid.Inside(candidate)==kInside) {
                double low=0,high=limit;
                for(int iteration=0;iteration<48 && high-low>std::max(1.e-12*mm,push*1.e-3);++iteration) {
                    const double middle=(low+high)*0.5;
                    if(solid.Inside(point+inward*middle)==kInside) high=middle;
                    else low=middle;
                }
                candidate=point+inward*std::min(high+push,limit);
                if(solid.Inside(candidate)==kInside) return candidate;
            }
        }
    }
    std::ostringstream detail;
    detail << std::setprecision(17) << "cannot relocate reflected point within numerical bound"
           << "; solid=" << solid.GetName() << "; point_mm=" << point/mm << "; candidate_mm=" << candidate/mm
           << "; surface_normal=" << solid.SurfaceNormal(point)
           << "; outside_distance_mm=" << solid.DistanceToIn(point)/mm
           << "; displacement_mm=" << (candidate-point).mag()/mm
           << "; allowed_mm=" << (maxRoundoffRepair+push)/mm;
    Fail(detail.str());
    return point; // Fail always throws.
}

double SpecularFraction(double distance_mm) {
    if(distance_mm<=cfg.diffuse_length_mm) return 0;
    if(cfg.transition_length_mm==0) return 1;
    return std::min(1.0,(distance_mm-cfg.diffuse_length_mm)/cfg.transition_length_mm);
}
void InstallBoundaryProcess() {
    if(!Enabled()) return;
    if(G4Threading::IsMultithreadedApplication()) Fail("spatial coating requires sequential G4RunManager");
    auto* manager=G4OpticalPhoton::Definition()->GetProcessManager();
    auto* list=manager->GetProcessList();
    G4VProcess* previous=nullptr;
    for(std::size_t i=0;i<list->size();++i) {
        if(dynamic_cast<SpatialBoundary*>((*list)[i])) return;
        if(dynamic_cast<G4OpBoundaryProcess*>((*list)[i])) previous=(*list)[i];
    }
    if(!previous) Fail("native optical boundary process missing: enable Optical physics");
    manager->RemoveProcess(previous);
    delete previous;
    manager->AddDiscreteProcess(new SpatialBoundary);
}
void ConfigurePhysics() {
    if(!Enabled()) return;
    cfg=ReadConfig();
    auto* p=G4OpticalParameters::Instance();
    p->SetProcessActivation("Scintillation",true);
    p->SetProcessActivation("Cerenkov",cfg.cherenkov);
    p->SetProcessActivation("OpAbsorption",true);
    p->SetProcessActivation("OpBoundary",true);
    p->SetProcessActivation("OpWLS",false);
    p->SetProcessActivation("OpWLS2",false);
    p->SetProcessActivation("OpRayleigh",false);
    p->SetProcessActivation("OpMieHG",false);
    p->SetScintByParticleType(false);
    p->SetScintFiniteRiseTime(false);
    p->SetScintTrackSecondariesFirst(false);
#if G4VERSION_NUMBER < 1100
    p->SetScintEnhancedTimeConstants(false);
#endif
}
void Construct(G4VPhysicalVolume* world, G4VPhysicalVolume* crystal, G4VPhysicalVolume* sensor) {
    if(!Enabled()) return;
    cfg=ReadConfig();
    if(!world || !crystal || !sensor) Fail("load exactly one bdx_bgo and bdx_flux geometry");
    crystalPV=crystal; sensorPV=sensor; worldPV=world;
    auto* shape=dynamic_cast<G4GenericTrap*>(crystal->GetLogicalVolume()->GetSolid());
    auto* pad=dynamic_cast<G4Box*>(sensor->GetLogicalVolume()->GetSolid());
    auto* fluxShape=dynamic_cast<G4GenericTrap*>(sensor->GetLogicalVolume()->GetSolid());
    if(!shape || (!pad && !fluxShape)) Fail("unexpected crystal or collector solid");
    fullBase=fluxShape!=nullptr;
    if(!dynamic_cast<ConservativeBgoTrap*>(shape)) {
        shape=static_cast<G4GenericTrap*>(MakeConservativeBgoSolid(*shape));
        crystal->GetLogicalVolume()->SetSolid(shape);
    }
    if(crystal->GetMotherLogical()!=world->GetLogicalVolume() ||
       sensor->GetMotherLogical()!=world->GetLogicalVolume()) Fail("crystal and collector must be world daughters");
    if(crystal->GetTranslation().mag()>1.e-9*mm || crystal->GetRotation() || sensor->GetRotation()) {
        // GEMC assigns identity matrices even for zero rotations.
        if(crystal->GetTranslation().mag()>1.e-9*mm ||
           (crystal->GetRotation() && !crystal->GetRotation()->isIdentity()) ||
           (sensor->GetRotation() && !sensor->GetRotation()->isIdentity())) Fail("study requires local unrotated crystal at origin");
    }
    const auto& v=shape->GetVertices();
    dimensions={{shape->GetZHalfLength()/mm, v[6].x()/mm, v[7].x()/mm, v[6].y()/mm,
                 v[2].x()/mm, v[3].x()/mm, v[2].y()/mm}};
    sensorZ=sensor->GetTranslation().z()-(pad?pad->GetZHalfLength():fluxShape->GetZHalfLength());
    sensorY=dimensions[3]*(dimensions[1]-dimensions[2])/(3*(dimensions[1]+dimensions[2]))*mm;
    if(fullBase) {
        if(sensor->GetTranslation().x()!=0 || sensor->GetTranslation().y()!=0)
            Fail("full-base collector must be centred on the crystal axis");
        const auto& f=fluxShape->GetVertices();
        for(int i=0;i<8;++i) if((f[i]-v[4+i%4]).mag()>1.e-7*mm)
            Fail("collector must reproduce the exact large-base quadrilateral");
    }
    if(std::abs(sensorZ-shape->GetZHalfLength())>1.e-7*mm) Fail("collector must touch the +z large base without overlap or gap");
    const double length=2*shape->GetZHalfLength()/mm;
    if(std::abs(cfg.diffuse_length_mm+cfg.transition_length_mm+cfg.specular_length_mm-length)>1.e-6)
        Fail("diffuse + transition + specular lengths must equal crystal length");
    specularProbability=nullptr;
    G4LogicalBorderSurface::CleanSurfaceTable();
    G4LogicalSkinSurface::CleanSurfaceTable();
    Material(world->GetLogicalVolume(),"BDX_optical_air",1.0,false);
    Material(crystal->GetLogicalVolume(),"BDX_optical_BGO",cfg.refractive_index,true);
    // Ideal collection: same refractive index at the entrance, no PDE, package,
    // grease, noise, saturation or electronics. Each transmitted photon ends once.
    Material(sensor->GetLogicalVolume(),"BDX_ideal_readout",cfg.refractive_index,false);
    if(cfg.wrap_model!="bare") {
        auto* wrap=new G4OpticalSurface("BDX_white_paint");
        wrap->SetModel(unified);
        wrap->SetType(dielectric_metal); // opaque effective paint: reflection or absorption
        wrap->SetFinish(cfg.wrap_model=="paint" ? ground : polished);
        wrap->SetSigmaAlpha(0);
        auto* props=new G4MaterialPropertiesTable;
        Constant(props,"REFLECTIVITY",cfg.wrap_reflectivity);
        Constant(props,"EFFICIENCY",0);
        Constant(props,"SPECULARSPIKECONSTANT",1);
        Constant(props,"SPECULARLOBECONSTANT",0);
        Constant(props,"BACKSCATTERCONSTANT",0);
        if(cfg.wrap_model=="paint") specularProbability=props->GetProperty("SPECULARSPIKECONSTANT");
        wrap->SetMaterialPropertiesTable(props);
        // Ordered border: it covers all crystal faces except the collector face,
        // whose neighbour is a separate volume rather than the world.
        new G4LogicalBorderSurface("BDX_crystal_wrap",crystal,world,wrap);
    }
    if(crystal->CheckOverlaps(20000,0,true)) Fail("crystal geometry overlaps another volume");
    if(sensor->CheckOverlaps(20000,0,true)) Fail("collector geometry overlaps another volume");
    if(!events.is_open()) {
        events.open(Env("BDX_OPTICS_CSV"));
        if(!events) Fail("cannot open event CSV");
        events << "event,primary_z_mm,nominal_distance_mm,edep_MeV,edep_axial_distance_mm,edep_euclidean_distance_mm,generated_scint,arrived_scint,generated_cherenkov,arrived_cherenkov,injected_optical,arrived_injected,bulk_absorbed,surface_absorbed,escaped,reflections,diffuse_reflections,specular_reflections,max_optical_steps,boundary_repairs,max_boundary_shift_mm,secondary_birth_repairs,max_birth_shift_mm\n";
        if(const auto* path=std::getenv("BDX_OPTICS_ARRIVALS")) {
            arrivals.open(path);
            if(!arrivals) Fail("cannot open arrivals CSV");
            arrivals << "event,track,parent,origin,time_ns,wavelength_nm,path_mm,steps\n";
        }
    }
    G4cout << "BDX_OPTICS_V2: " << G4Version << "; conservative BGO safety; wrap=" << cfg.wrap_model
           << "; R=" << cfg.wrap_reflectivity << "; absorption=" << cfg.absorption_length_mm
           << " mm; collector=" << (fullBase?"entire large base":"test aperture")
           << "; collection plane z=" << sensorZ/mm << " mm" << G4endl;
}
void BeginEvent(const G4Event* event) {
    if(!Enabled()) return;
    ev=Event{};
    ev.index=event->GetEventID();
    const auto* v=event->GetPrimaryVertex();
    if(!v) Fail("event has no primary vertex");
    ev.primaryZ=v->GetZ0()/mm;
}
void EndEvent(const G4Event* event) {
    if(!Enabled()) return;
    if(event->IsAborted()) Fail("aborted event is not a valid optical measurement");
    const long generated=ev.generated+ev.generatedC+ev.injected;
    const long terminated=ev.arrived+ev.arrivedC+ev.arrivedI+ev.bulk+ev.surface+ev.escaped;
    if(generated!=terminated) Fail("photon budget does not close: generated="+std::to_string(generated)+
                                  " terminated="+std::to_string(terminated));
    events << std::setprecision(12) << event->GetEventID() << ',' << ev.primaryZ << ','
           << sensorZ/mm-ev.primaryZ << ',' << ev.edep << ',';
    if(ev.edep>0) events << sensorZ/mm-ev.weightedZ/ev.edep << ',' << ev.weightedD/ev.edep;
    else events << "nan,nan";
    events << ',' << ev.generated << ',' << ev.arrived << ',' << ev.generatedC << ',' << ev.arrivedC
           << ',' << ev.injected << ',' << ev.arrivedI << ',' << ev.bulk << ',' << ev.surface << ','
           << ev.escaped << ',' << ev.reflections << ',' << ev.diffuse << ',' << ev.specular << ',' << ev.maxSteps
           << ',' << ev.boundaryRepairs << ',' << ev.maxBoundaryShift
           << ',' << ev.birthRepairs << ',' << ev.maxBirthShift << '\n';
    events.flush();
    if(arrivals.is_open()) arrivals.flush();
    if(!events || (arrivals.is_open() && !arrivals)) Fail("CSV write failure");
}
bool Step(const G4Step* step) {
    if(!Enabled()) return false;
    auto* track=step->GetTrack();
    const auto* pre=step->GetPreStepPoint();
    auto* post=step->GetPostStepPoint();
    const auto* from=pre->GetPhysicalVolume();
    const auto* to=post->GetPhysicalVolume();
    const bool optical=track->GetDefinition()==G4OpticalPhoton::Definition();
    if(!optical) {
        if(from==crystalPV) {
            const double edep=step->GetTotalEnergyDeposit()/MeV;
            const auto midpoint=(pre->GetPosition()+post->GetPosition())*0.5;
            ev.edep+=edep;
            ev.weightedZ+=edep*midpoint.z()/mm;
            ev.weightedD+=edep*(midpoint-G4ThreeVector(0,sensorY,sensorZ)).mag()/mm;
            const int source=static_cast<int>(ev.deposits.size());
            int scint=0,cerenkov=0;
            for(const auto* child:*step->GetSecondaryInCurrentStep()) {
                if(child->GetDefinition()!=G4OpticalPhoton::Definition()) continue;
                // Prepare newborn tracks before the stack starts transporting
                // them. Scintillation close to a charged-particle boundary can
                // inherit its candidate outside touchable or a rounded position
                // just outside BGO. Repair within the same numerical bound;
                // never transport such a photon through the world and call it loss.
                auto* newborn=const_cast<G4Track*>(child);
                if(newborn->GetUserInformation()) Fail("optical child already has user track information");
                newborn->SetUserInformation(new PhotonInfo(source,*newborn));
                const auto point=newborn->GetPosition();
                bool repaired=false;
                auto* solid=crystalPV->GetLogicalVolume()->GetSolid();
                if(cfg.boundary_push_mm>0 && solid->Inside(point)!=kInside) {
                    const auto inside=RelocateReflectedPoint(*solid,point,cfg.boundary_push_mm*mm);
                    newborn->SetPosition(inside);
                    ev.maxBirthShift=std::max(ev.maxBirthShift,(inside-point).mag()/mm);
                    repaired=true;
                }
                if(!newborn->GetTouchable() || newborn->GetVolume()!=crystalPV) {
                    newborn->SetTouchableHandle(pre->GetTouchableHandle());
                    newborn->SetNextTouchableHandle(pre->GetTouchableHandle());
                    repaired=true;
                }
                if(repaired) ++ev.birthRepairs;
                const auto* creator=child->GetCreatorProcess();
                if(creator && creator->GetProcessName()=="Scintillation") { ++ev.generated; ++scint; }
                else if(creator && creator->GetProcessName()=="Cerenkov") { ++ev.generatedC; ++cerenkov; }
                else Fail("unexpected optical production process");
            }
            if(edep>0 || scint || cerenkov) ev.deposits.push_back({track->GetTrackID(),track->GetParentID(),
                track->GetDefinition()->GetPDGEncoding(),track->GetCurrentStepNumber(),scint,cerenkov,
                pre->GetPosition(),post->GetPosition(),pre->GetGlobalTime(),post->GetGlobalTime(),edep});
        }
        return false;
    }
    const auto* creator=track->GetCreatorProcess();
    const auto origin=creator ? creator->GetProcessName() : "primary";
    if(track->GetCurrentStepNumber()==1 && !creator) {
        if(from!=crystalPV) Fail("optical test primary must start inside BGO");
        ++ev.injected;
    }

    auto* info=dynamic_cast<PhotonInfo*>(track->GetUserInformation());
    if(!info) {
        if(creator) Fail("optical secondary lacks its BGO source record");
        info=new PhotonInfo(-1,*track);
        info->birth=track->GetVertexPosition();
        info->birthTime=pre->GetGlobalTime();
        track->SetUserInformation(info);
    }
    ev.maxSteps=std::max(ev.maxSteps,static_cast<long>(track->GetCurrentStepNumber()));
    G4OpBoundaryProcessStatus status=Undefined;
    if(post->GetStepStatus()==fGeomBoundary) {
        const auto* processes=track->GetDefinition()->GetProcessManager()->GetProcessList();
        for(std::size_t i=0;i<processes->size();++i)
            if(auto* boundary=dynamic_cast<G4OpBoundaryProcess*>((*processes)[i])) status=boundary->GetStatus();
        if(from!=crystalPV && status==StepTooSmall && step->GetStepLength()<1.e-6*mm) {
            // Geant4 10.7 boundary relocation after a reflection temporarily has
            // the world as the pre-volume. This zero-length step is not escape.
            track->UseGivenVelocity(true);
            track->SetVelocity(c_light/cfg.refractive_index);
            return true;
        }
        if(status==NoRINDEX) Fail("NoRINDEX is a configuration error, not physical attenuation");
        if(to==sensorPV && Transmitted(status)) {
            if(origin=="Scintillation") ++ev.arrived;
            else if(origin=="Cerenkov") ++ev.arrivedC;
            else if(origin=="primary") ++ev.arrivedI;
            else Fail("arrival from unexpected optical process");
            ev.flux.push_back({track->GetTrackID(),track->GetParentID(),
                origin=="Scintillation"?1:origin=="Cerenkov"?2:0,info->source,
                track->GetCurrentStepNumber(),info->reflections,post->GetPosition(),pre->GetMomentumDirection(),
                info->birth,post->GetGlobalTime(),pre->GetKineticEnergy(),track->GetTrackLength(),
                info->birthTime,track->GetWeight()});
            if(arrivals.is_open()) {
                const auto* event=G4RunManager::GetRunManager()->GetCurrentEvent();
                arrivals << std::setprecision(12) << event->GetEventID() << ',' << track->GetTrackID() << ','
                         << track->GetParentID() << ',' << origin << ',' << post->GetGlobalTime()/ns << ','
                         << h_Planck*c_light/pre->GetKineticEnergy()/nm << ',' << track->GetTrackLength()/mm
                         << ',' << track->GetCurrentStepNumber() << '\n';
            }
            track->SetTrackStatus(fStopAndKill);
            return true;
        }
        if(to==worldPV && Transmitted(status)) {
            ++ev.escaped; track->SetTrackStatus(fStopAndKill); return true;
        }
        if(status==Absorption) { ++ev.surface; return true; }
        if(status==Detection) Fail("unexpected boundary detection: SiPM is an entrance counter");
        if(status!=StepTooSmall && status!=Undefined && status!=NotAtBoundary) ++ev.reflections;
        const bool reflected=status==LambertianReflection || status==SpikeReflection ||
            status==LobeReflection || status==BackScattering || status==FresnelReflection ||
            status==TotalInternalReflection;
        if(reflected) ++info->reflections;
        if(reflected && cfg.boundary_push_mm>0) {
            // Transportation's post-step touchable still names the candidate
            // outside volume. Synchronise it with the strictly internal proposed
            // position before the next step; otherwise an entire BGO flight can
            // incorrectly use the world's material and Fresnel index.
            auto* navigator=G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
            G4TouchableHandle corrected=new G4TouchableHistory;
            navigator->LocateGlobalPointAndUpdateTouchableHandle(post->GetPosition(),
                track->GetMomentumDirection(),corrected,false);
            if(corrected->GetVolume()!=crystalPV) Fail("reflected track relocation selected the wrong volume");
            post->SetTouchableHandle(corrected);
            post->SetMaterial(crystalPV->GetLogicalVolume()->GetMaterial());
            post->SetMaterialCutsCouple(pre->GetMaterialCutsCouple());
            post->SetSensitiveDetector(pre->GetSensitiveDetector());
            track->SetTouchableHandle(corrected);
            track->SetNextTouchableHandle(corrected);
        }
        if(status==LambertianReflection) ++ev.diffuse;
        if(status==SpikeReflection || status==LobeReflection) ++ev.specular;
    }
    if(from!=crystalPV) {
        std::ostringstream detail;
        detail << "unexpected optical propagation outside BGO; from=" << (from?from->GetName():"null")
               << " to=" << (to?to->GetName():"null") << " origin=" << origin
               << " step=" << track->GetCurrentStepNumber() << " length_mm=" << step->GetStepLength()/mm
               << " pre=" << pre->GetPosition()/mm << " post=" << post->GetPosition()/mm
               << " boundary_status=" << status << " step_status=" << post->GetStepStatus();
        Fail(detail.str());
    }
    track->UseGivenVelocity(true);
    track->SetVelocity(c_light/cfg.refractive_index);
    const auto* process=post->GetProcessDefinedStep();
    if(track->GetTrackStatus()==fStopAndKill) {
        if(process && process->GetProcessName()=="OpAbsorption") ++ev.bulk;
        else Fail("optical track terminated without physical absorption, escape or SiPM arrival");
    }
    // No energy, time or step-count truncation: survival is a physical result.
    return true;
}
void RegisterHipoSchemas(hipo::dictionary& dictionary) {
    if(Enabled()) AddOpticalSchemas(dictionary);
}
void AppendHipo(hipo::event& output,hipo::dictionary& dictionary,int gemcEvent) {
    if(!Enabled()) return;
    if(ev.index<0) Fail("HIPO output requested before event scoring");
    auto append=[&](hipo::bank& b) {
        const auto required=static_cast<std::size_t>(output.getSize())+b.getSize()+8+1024;
        if(output.getEventBuffer().size()<required) output.getEventBuffer().resize(required);
        output.addStructure(b);
    };
    auto ints=[](hipo::bank& b,int row,std::initializer_list<std::pair<const char*,int>> values) {
        for(const auto& v:values) b.putInt(v.first,row,v.second);
    };
    auto doubles=[](hipo::bank& b,int row,std::initializer_list<std::pair<const char*,double>> values) {
        for(const auto& v:values) b.putDouble(v.first,row,v.second);
    };
    hipo::bank c(dictionary.getSchema("BDX::Config"),1);
    const auto envInt=[](const char* name) { const auto* p=std::getenv(name); return p?std::stoi(p):0; };
    ints(c,0,{{"schema_version",2},{"crystal_type",envInt("BDX_OPTICS_TYPE")},
        {"seed",envInt("BDX_OPTICS_SEED")},{"collector",fullBase?1:0},
        {"wrap_model",cfg.wrap_model=="paint"?1:cfg.wrap_model=="contact"?2:0},{"cherenkov",cfg.cherenkov?1:0}});
    doubles(c,0,{{"light_yield_per_MeV",cfg.light_yield_per_MeV},{"refractive_index",cfg.refractive_index},
        {"absorption_length_mm",cfg.absorption_length_mm},{"decay_time_ns",cfg.decay_time_ns},
        {"wavelength_nm",cfg.wavelength_nm},{"wrap_reflectivity",cfg.wrap_reflectivity},
        {"diffuse_length_mm",cfg.diffuse_length_mm},{"transition_length_mm",cfg.transition_length_mm},
        {"specular_length_mm",cfg.specular_length_mm},{"boundary_push_mm",cfg.boundary_push_mm},
        {"half_length_mm",dimensions[0]},{"Bg_mm",dimensions[1]},{"Bm_mm",dimensions[2]},
        {"H_mm",dimensions[3]},{"bg_mm",dimensions[4]},{"bm_mm",dimensions[5]},{"h_mm",dimensions[6]},
        {"center_y_mm",sensorY/mm},{"plane_z_mm",sensorZ/mm}});
    append(c);
    hipo::bank e(dictionary.getSchema("BDX::Event"),1);
    ints(e,0,{{"event",gemcEvent},{"event_index",ev.index},{"generated_scint",int(ev.generated)},
        {"arrived_scint",int(ev.arrived)},{"generated_cherenkov",int(ev.generatedC)},{"arrived_cherenkov",int(ev.arrivedC)},
        {"injected_optical",int(ev.injected)},{"arrived_injected",int(ev.arrivedI)},{"bulk_absorbed",int(ev.bulk)},
        {"surface_absorbed",int(ev.surface)},{"escaped",int(ev.escaped)},{"reflections",int(ev.reflections)},
        {"diffuse_reflections",int(ev.diffuse)},{"specular_reflections",int(ev.specular)},
        {"max_optical_steps",int(ev.maxSteps)},{"boundary_repairs",int(ev.boundaryRepairs)},
        {"secondary_birth_repairs",int(ev.birthRepairs)}});
    const auto nan=std::numeric_limits<double>::quiet_NaN();
    doubles(e,0,{{"primary_z_mm",ev.primaryZ},{"nominal_distance_mm",sensorZ/mm-ev.primaryZ},
        {"edep_MeV",ev.edep},{"edep_axial_distance_mm",ev.edep>0?sensorZ/mm-ev.weightedZ/ev.edep:nan},
        {"edep_euclidean_distance_mm",ev.edep>0?ev.weightedD/ev.edep:nan},
        {"max_boundary_shift_mm",ev.maxBoundaryShift},{"max_birth_shift_mm",ev.maxBirthShift}});
    append(e);
    hipo::bank d(dictionary.getSchema("BDX::Edep"),ev.deposits.size());
    for(std::size_t i=0;i<ev.deposits.size();++i) {
        const auto& v=ev.deposits[i];
        ints(d,i,{{"source_step",int(i)},{"track",v.track},{"parent",v.parent},{"pid",v.pid},{"step",v.step},
            {"generated_scint",v.scint},{"generated_cherenkov",v.cerenkov}});
        doubles(d,i,{{"x0_mm",v.start.x()/mm},{"y0_mm",v.start.y()/mm},{"z0_mm",v.start.z()/mm},
            {"x1_mm",v.end.x()/mm},{"y1_mm",v.end.y()/mm},{"z1_mm",v.end.z()/mm},
            {"t0_ns",v.t0/ns},{"t1_ns",v.t1/ns},{"edep_MeV",v.energy}});
    }
    append(d);
    hipo::bank f(dictionary.getSchema("BDX::Flux"),ev.flux.size());
    for(std::size_t i=0;i<ev.flux.size();++i) {
        const auto& v=ev.flux[i];
        ints(f,i,{{"track",v.track},{"parent",v.parent},{"origin",v.origin},{"source_step",v.source},
            {"steps",v.steps},{"reflections",v.reflections}});
        doubles(f,i,{{"x_mm",v.position.x()/mm},{"y_mm",v.position.y()/mm},{"z_mm",v.position.z()/mm},
            {"lx_mm",v.position.x()/mm},{"ly_mm",(v.position.y()-sensorY)/mm},
            {"ux",v.direction.x()},{"uy",v.direction.y()},{"uz",v.direction.z()},
            {"time_ns",v.time/ns},{"energy_eV",v.energy/eV},{"wavelength_nm",h_Planck*c_light/v.energy/nm},
            {"path_mm",v.path/mm},{"birth_x_mm",v.birth.x()/mm},{"birth_y_mm",v.birth.y()/mm},
            {"birth_z_mm",v.birth.z()/mm},{"birth_time_ns",v.birthTime/ns},{"weight",v.weight}});
    }
    append(f);
}

}
