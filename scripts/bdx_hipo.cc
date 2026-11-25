// Offline export and integrity checks. No Geant4 dependency.
#include "hipo4/reader.h"
#include "hipo4/bank.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>

void require(bool ok,const std::string& why) { if(!ok) throw std::runtime_error(why); }
bool close(double a,double b) { return std::abs(a-b)<=1.e-8*std::max(1.,std::max(std::abs(a),std::abs(b))); }
void jsonBank(std::ostream& out,const hipo::bank& bank,int row) {
    // The schema names are fixed identifier strings; every field is numerical.
    auto& schema=const_cast<hipo::bank&>(bank).getSchema();
    out << '{';
    for(int col=0;col<schema.getEntries();++col) {
        if(col) out << ',';
        out << '"' << schema.getEntryName(col) << "\":";
        const int type=schema.getEntryType(col);
        if(type==hipo::kLong) out << bank.getLong(col,row);
        else if(type==hipo::kDouble || type==hipo::kFloat) {
            const double value=type==hipo::kDouble?bank.getDouble(col,row):bank.getFloat(col,row);
            if(std::isfinite(value)) out << value; else out << "null";
        } else out << bank.getInt(col,row);
    }
    out << '}';
}
int main(int argc,char** argv) { try {
    if(argc<4 || argc>5) throw std::runtime_error("Usage: bdx_hipo FILE --check REPORT.json | --export DIRECTORY | --reduce REPORT.json [square_mm]");
    const std::string mode=argv[2];
    const bool exporting=mode=="--export", reducing=mode=="--reduce";
    require(exporting || reducing || mode=="--check","unknown operation");
    require(argc==4 || reducing,"square selection is only valid with --reduce");
    const double square=argc==5?std::stod(argv[4]):0;
    require(argc==4 || (std::isfinite(square) && square>0),"invalid square size");
    const std::string destination=argv[3];
    std::ofstream reduced;
    if(reducing) {
        reduced.open(destination+".tmp");require(bool(reduced),"cannot create reduction");
        reduced << std::setprecision(17) << "{\"format\":\"BDX_OPTICS_REDUCTION_V1\",\"events\":[";
    }
    hipo::reader reader; reader.open(argv[1]);
    require(reader.is_open(),"cannot open HIPO input");
    hipo::dictionary dictionary; reader.readDictionary(dictionary);
    std::map<std::string,std::unique_ptr<hipo::bank>> banks;
    std::map<std::string,long long> totals;
    std::map<std::string,std::unique_ptr<std::ofstream>> streams;
    const std::set<std::string> exports={"BDX::Config","BDX::Event","BDX::Edep","BDX::Flux"};
    for(const auto& name:dictionary.getSchemaList()) {
        banks[name].reset(new hipo::bank(dictionary.getSchema(name.c_str())));
        totals[name]=0;
        if(exporting && exports.count(name)) {
            const auto base=name.substr(5);
            auto& stream=streams[name]; stream.reset(new std::ofstream(destination+"/"+base+".csv"));
            require(bool(*stream),"cannot create CSV export");
            *stream << "record_index";
            auto& schema=banks[name]->getSchema();
            for(int c=0;c<schema.getEntries();++c) *stream << ',' << schema.getEntryName(c);
            *stream << '\n' << std::setprecision(17);
        }
    }
    for(const auto* name:{"BDX::Config","BDX::Event","BDX::Edep","BDX::Flux","RUN::config","MC::Particle"})
        require(banks.count(name),std::string("missing schema ")+name);
    int index=0; hipo::event event;
    while(reader.next()) {
        reader.read(event);
        for(auto& item:banks) {
            auto& bank=*item.second; bank.setRows(0); event.getStructure(bank);
            totals[item.first]+=bank.getRows();
            if(streams.count(item.first)) {
                auto& out=*streams[item.first]; const auto& schema=bank.getSchema();
                for(int row=0;row<bank.getRows();++row) {
                    out << index;
                    for(int col=0;col<schema.getEntries();++col) {
                        out << ',';
                        switch(schema.getEntryType(col)) {
                            case hipo::kFloat: out << bank.getFloat(col,row); break;
                            case hipo::kDouble: out << bank.getDouble(col,row); break;
                            case hipo::kLong: out << bank.getLong(col,row); break;
                            default: out << bank.getInt(col,row);
                        }
                    }
                    out << '\n';
                }
                require(bool(out),"CSV write failure");
            }
        }
        const auto& e=*banks.at("BDX::Event"), &c=*banks.at("BDX::Config"), &f=*banks.at("BDX::Flux"), &d=*banks.at("BDX::Edep");
        require(e.getRows()==1 && c.getRows()==1,"missing event/config row");
        require(c.getInt("schema_version",0)==2,"unsupported BDX schema version");
        require(e.getInt("event_index",0)==index,"missing or reordered optical event");
        require(banks.at("RUN::config")->getRows()==1,"missing RUN::config");
        require(e.getInt("event",0)==banks.at("RUN::config")->getInt("event",0),"event id mismatch");
        const int generated=e.getInt("generated_scint",0)+e.getInt("generated_cherenkov",0)+e.getInt("injected_optical",0);
        const int arrived=e.getInt("arrived_scint",0)+e.getInt("arrived_cherenkov",0)+e.getInt("arrived_injected",0);
        require(generated==arrived+e.getInt("bulk_absorbed",0)+e.getInt("surface_absorbed",0)+e.getInt("escaped",0),"photon budget mismatch");
        require(f.getRows()==arrived,"lost or duplicated flux rows");
        double energy=0; int scint=0,cerenkov=0;
        for(int i=0;i<d.getRows();++i) {
            require(d.getInt("source_step",i)==i,"invalid source step");
            energy+=d.getDouble("edep_MeV",i); scint+=d.getInt("generated_scint",i); cerenkov+=d.getInt("generated_cherenkov",i);
        }
        require(close(energy,e.getDouble("edep_MeV",0)) && scint==e.getInt("generated_scint",0) && cerenkov==e.getInt("generated_cherenkov",0),"source/event mismatch");
        std::set<int> tracks; int origins[3]={0,0,0};
        std::vector<int> sourceScint(d.getRows(),0),sourceCerenkov(d.getRows(),0);
        for(int i=0;i<f.getRows();++i) {
            require(tracks.insert(f.getInt("track",i)).second,"duplicate photon first arrival");
            const int origin=f.getInt("origin",i),source=f.getInt("source_step",i);
            require(origin>=0 && origin<=2,"invalid photon origin"); ++origins[origin];
            if(origin) {
                require(source>=0 && source<d.getRows(),"missing photon source step");
                (origin==1?sourceScint:sourceCerenkov)[source]++;
            } else require(source==-1,"injected photon has source step");
            const double x=f.getDouble("x_mm",i),y=f.getDouble("y_mm",i),z=f.getDouble("z_mm",i);
            require(std::isfinite(x) && std::isfinite(y) && close(z,c.getDouble("plane_z_mm",0)),"arrival outside collection plane");
            const double H=c.getDouble("H_mm",0),Bg=c.getDouble("Bg_mm",0),Bm=c.getDouble("Bm_mm",0);
            const double width=Bm+(Bg-Bm)*(y+H)/(2*H);
            require(std::abs(y)<=H+1e-5 && std::abs(x)<=width+1e-5,"arrival outside large-base trapezoid");
            require(close(f.getDouble("lx_mm",i),x) && close(f.getDouble("ly_mm",i),y-c.getDouble("center_y_mm",0)),"invalid local coordinates");
            const double ux=f.getDouble("ux",i),uy=f.getDouble("uy",i),uz=f.getDouble("uz",i);
            require(uz>0 && close(ux*ux+uy*uy+uz*uz,1),"invalid incident direction");
            require(f.getDouble("path_mm",i)>0 && f.getDouble("time_ns",i)>=f.getDouble("birth_time_ns",i),"invalid photon history");
            require(close(f.getDouble("energy_eV",i)*f.getDouble("wavelength_nm",i),1239.841984332),"invalid photon wavelength");
            require(f.getDouble("weight",i)==1,"unsupported weighted production");
        }
        require(origins[1]==e.getInt("arrived_scint",0) && origins[2]==e.getInt("arrived_cherenkov",0) && origins[0]==e.getInt("arrived_injected",0),"origin count mismatch");
        for(int i=0;i<d.getRows();++i) require(sourceScint[i]<=d.getInt("generated_scint",i) && sourceCerenkov[i]<=d.getInt("generated_cherenkov",i),"source arrival excess");
        if(reducing) {
            int selected=0;
            for(int row=0;row<f.getRows();++row) if(f.getInt("origin",row)==1 &&
                (argc==4 || (std::abs(f.getDouble("lx_mm",row))<=square/2 && std::abs(f.getDouble("ly_mm",row))<=square/2))) ++selected;
            if(index) reduced << ',';
            reduced << "{\"config\":"; jsonBank(reduced,c,0);
            reduced << ",\"event\":"; jsonBank(reduced,e,0);
            reduced << ",\"selected_scint\":" << selected << ",\"primaries\":[";
            const auto& primaries=*banks.at("MC::Particle");
            for(int row=0;row<primaries.getRows();++row) { if(row) reduced << ',';jsonBank(reduced,primaries,row); }
            reduced << "]}";
        }
        ++index;
    }
    require(index==reader.getEntries(),"HIPO index/event count mismatch");
    if(reducing) {
        reduced << "]}\n";reduced.flush();require(bool(reduced),"reduction write failure");reduced.close();
        require(std::rename((destination+".tmp").c_str(),destination.c_str())==0,"cannot publish reduction");
        return 0;
    }
    const auto report=exporting?destination+"/banks.json":destination;
    std::ofstream json(report); require(bool(json),"cannot create bank inventory");
    json << "{\"format\":\"BDX_OPTICS_HIPO_V2\",\"events\":" << index << ",\"banks\":{";
    bool first=true;
    for(const auto& item:totals) { if(!first) json << ','; first=false; json << '"' << item.first << "\":" << item.second; }
    json << "}}\n"; require(bool(json),"inventory write failure");
    std::cout << "Validated " << index << " events; " << totals["BDX::Flux"] << " optical first arrivals\n";
    return 0;
} catch(const std::exception& e) { std::cerr << "BDX HIPO: " << e.what() << '\n'; return 1; } }
