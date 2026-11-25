// Regression for the actual BGO-1 grazing intersection that stopped the user's run.
#include "BdxCrystalOptics.h"
#include "G4GenericTrap.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <fstream>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
void Require(bool condition,const char* message) {
    if(!condition) throw std::runtime_error(message);
}
std::vector<G4TwoVector> Vertices(const std::string& root,int type,double& halfLength) {
    std::ifstream input(root+"/experiments/bdx_optics/bdx_optics__geometry_bgo"+std::to_string(type)+".txt");
    Require(bool(input),"missing crystal geometry fixture");
    std::string line,field;std::getline(input,line);std::stringstream row(line);
    for(int i=0;i<=7;++i) std::getline(row,field,'|');
    std::string::size_type p;
    while((p=field.find("*mm"))!=std::string::npos) field.erase(p,3);
    std::stringstream dimensions(field);dimensions>>halfLength;
    std::vector<G4TwoVector> vertices;
    for(int i=0;i<8;++i) {
        double x,y;dimensions>>x>>y;
        Require(bool(dimensions),"invalid crystal dimensions");
        vertices.emplace_back(x*mm,y*mm);
    }
    halfLength*=mm;
    return vertices;
}
}

int main(int argc,char** argv) {
    if(argc!=2) return 2;
    G4Random::setTheSeed(2010314138);
    long checked=0;double maximum=0;
    for(int type=1;type<=8;++type) {
        double halfLength;auto vertices=Vertices(argv[1],type,halfLength);
        std::cout << "RELOCATION_TYPE " << type << std::endl;
        G4GenericTrap solid("regression_bgo"+std::to_string(type),halfLength,vertices);
        if(type==1) {
            // Captured with the original run's exact seed and input. DistanceToOut
            // put this hit 1.276 nm outside: the previous fixed 1 nm push fails.
            const G4ThreeVector point(-14.554468366128804*mm,-22.489342449694302*mm,75.987450886933132*mm);
            Require(solid.Inside(point-solid.SurfaceNormal(point)*1.e-6*mm)!=kInside,
                    "fixture no longer reproduces the fixed-push failure");
            for(const double push:{1.e-7*mm,1.e-6*mm,1.e-5*mm}) {
                const auto repaired=bdxoptics::RelocateReflectedPoint(solid,point,push);
                Require(solid.Inside(repaired)==kInside,"regression point was not repaired");
                Require((repaired-point).mag()<2.e-5*mm,"regression correction is too large");
                maximum=std::max(maximum,(repaired-point).mag()/mm);++checked;
            }
        }
        auto* conservative=bdxoptics::MakeConservativeBgoSolid(solid);
        for(double z:{100.0,40.0,-20.0,-80.0}) for(double vz:{0.0,6.123233995736766e-17,-1.e-8,1.e-8,-1.e-5,1.e-5,-1.e-3,1.e-3}) {
            const G4ThreeVector start(-100*mm,0,z*mm),direction=G4ThreeVector(1,0,vz).unit();
            const double entry=conservative->DistanceToIn(start,direction);
            Require(entry>0 && entry<100*mm,"transverse beam missed BGO");
            const auto point=start+direction*(entry+1.e-6*mm);
            Require(conservative->Inside(point)==kInside,"transverse beam entry is not internal");
            const double chord=conservative->DistanceToOut(point,direction);
            Require(chord>10*mm && chord<80*mm,"transverse beam has a zero or invalid chord");
            Require(conservative->Inside(point+direction*(chord-1.e-5*mm))!=kOutside,"exit occurs before material ends");
            Require(conservative->Inside(point+direction*(chord+1.e-5*mm))!=kInside,"exit remains inside material");
        }
        for(int sample=0;sample<2000;++sample) {
            const auto boundary=solid.GetPointOnSurface();
            const auto inward=boundary-solid.SurfaceNormal(boundary)*1.e-4*mm;
            if(solid.Inside(inward)==kInside) {
                Require(conservative->Inside(inward)==kInside,"replacement changed the solid");
                Require(conservative->DistanceToOut(inward)<=(boundary-inward).mag()+1.e-8*mm,
                        "inside safety is not conservative");
            }
            const auto outward=boundary+solid.SurfaceNormal(boundary)*1.e-4*mm;
            if(solid.Inside(outward)==kOutside)
                Require(conservative->DistanceToIn(outward)<=(outward-boundary).mag()+1.e-8*mm,
                        "outside safety is not conservative");

            const auto surface=solid.GetPointOnSurface();
            const auto origin=surface*0.5;
            if(conservative->Inside(origin)==kInside) {
                const auto direction=G4ThreeVector(2*G4UniformRand()-1,2*G4UniformRand()-1,2*G4UniformRand()-1).unit();
                const double distance=conservative->DistanceToOut(origin,direction);
                Require(distance>0 && distance<1000*mm,"random ray has invalid exit");
                Require(conservative->Inside(origin+direction*(distance-1.e-5*mm))!=kOutside,"random exit precedes the boundary");
                const auto outside=origin+direction*(distance+1.e-5*mm);
                Require(conservative->Inside(outside)!=kInside,"random exit remains in the crystal");
                const double reverse=conservative->DistanceToIn(outside,-direction);
                Require(reverse>=0 && reverse<2.e-5*mm,"reversed ray missed the same boundary");
            }
            const auto normal=solid.SurfaceNormal(surface);
            // Cover exact surfaces and points displaced outside by rounding
            // error, including both trapezoidal ends and the lateral faces.
            for(const double rounding:{0.0,1.3e-6*mm,1.e-5*mm}) {
                const auto point=surface+normal*rounding;
                const auto repaired=bdxoptics::RelocateReflectedPoint(solid,point,1.e-6*mm);
                Require(solid.Inside(repaired)==kInside,"surface sample was not repaired");
                Require((repaired-point).mag()<=1.001e-3*mm,"surface correction exceeded the numerical bound");
                maximum=std::max(maximum,(repaired-point).mag()/mm);++checked;
            }
        }
        for(int i=0;i<8;++i) {
            const auto v=vertices[i];
            const G4ThreeVector point(v.x(),v.y(),i<4 ? -halfLength : halfLength);
            const auto repaired=bdxoptics::RelocateReflectedPoint(solid,point,1.e-6*mm);
            Require(solid.Inside(repaired)==kInside,"corner was not repaired");
            Require((repaired-point).mag()<=1.001e-3*mm,"corner correction exceeded the numerical bound");
            maximum=std::max(maximum,(repaired-point).mag()/mm);++checked;
        }
        if(type==1) {
            const G4ThreeVector p(-20.994472375460159,-0.027083894775988219,100.03453154064427);
            const double t=(p.z()+halfLength)/(2*halfLength);
            const double h=vertices[1].y()+t*(vertices[5].y()-vertices[1].y());
            const double bm=vertices[3].x()+t*(vertices[7].x()-vertices[3].x());
            const double bg=vertices[2].x()+t*(vertices[6].x()-vertices[2].x());
            const double left=-(bm+(bg-bm)*(p.y()+h)/(2*h));
            const double distance=p.x()-left;
            Require(solid.DistanceToOut(p)>distance,"fixture no longer exposes non-conservative native safety");
            Require(conservative->DistanceToOut(p)<=distance,"captured step safety is not conservative");
            std::cout << "SAFETY_REGRESSION native_mm=" << solid.DistanceToOut(p)
                      << " conservative_mm=" << conservative->DistanceToOut(p)
                      << " boundary_upper_bound_mm=" << distance << std::endl;
        }
        bool rejected=false;
        try { bdxoptics::RelocateReflectedPoint(solid,{1000*mm,0,0},1.e-6*mm); }
        catch(const std::runtime_error&) { rejected=true; }
        Require(rejected,"macroscopic geometry error was silently repaired");
    }
    std::cout << std::setprecision(12) << "RELOCATION_CHECKS " << checked
              << " MAX_SHIFT_MM " << maximum << "\n";
}
