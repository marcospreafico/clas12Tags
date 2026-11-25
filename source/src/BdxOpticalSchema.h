#ifndef BDX_OPTICAL_SCHEMA_H
#define BDX_OPTICAL_SCHEMA_H
#include "hipo4/dictionary.h"
#include "hipo4/writer.h"
#include <stdexcept>
// Independent of Geant4 so readers and schema tests share the file contract.
namespace bdxoptics {
inline void AddOpticalSchemas(hipo::dictionary& d) {
    auto add=[&](const char* name,int item,const char* fields) {
        hipo::schema s(name,4601,item); s.parse(fields); d.addSchema(s);
    };
    add("BDX::Config",1,"schema_version/I,crystal_type/I,seed/I,collector/I,wrap_model/I,cherenkov/I,light_yield_per_MeV/D,refractive_index/D,absorption_length_mm/D,decay_time_ns/D,wavelength_nm/D,wrap_reflectivity/D,diffuse_length_mm/D,transition_length_mm/D,specular_length_mm/D,boundary_push_mm/D,half_length_mm/D,Bg_mm/D,Bm_mm/D,H_mm/D,bg_mm/D,bm_mm/D,h_mm/D,center_y_mm/D,plane_z_mm/D");
    add("BDX::Event",2,"event/I,event_index/I,primary_z_mm/D,nominal_distance_mm/D,edep_MeV/D,edep_axial_distance_mm/D,edep_euclidean_distance_mm/D,generated_scint/I,arrived_scint/I,generated_cherenkov/I,arrived_cherenkov/I,injected_optical/I,arrived_injected/I,bulk_absorbed/I,surface_absorbed/I,escaped/I,reflections/I,diffuse_reflections/I,specular_reflections/I,max_optical_steps/I,boundary_repairs/I,max_boundary_shift_mm/D,secondary_birth_repairs/I,max_birth_shift_mm/D");
    add("BDX::Edep",3,"source_step/I,track/I,parent/I,pid/I,step/I,x0_mm/D,y0_mm/D,z0_mm/D,x1_mm/D,y1_mm/D,z1_mm/D,t0_ns/D,t1_ns/D,edep_MeV/D,generated_scint/I,generated_cherenkov/I");
    add("BDX::Flux",4,"track/I,parent/I,origin/I,source_step/I,x_mm/D,y_mm/D,z_mm/D,lx_mm/D,ly_mm/D,ux/D,uy/D,uz/D,time_ns/D,energy_eV/D,wavelength_nm/D,path_mm/D,steps/I,reflections/I,birth_x_mm/D,birth_y_mm/D,birth_z_mm/D,birth_time_ns/D,weight/D");
}
inline void WriteOpticalRecord(hipo::writer& writer,hipo::event& event) {
    // HIPO uses a strict capacity comparison, leaving one index slot unused.
    hipo::recordbuilder record(2,event.getSize()+1024);
    if(!record.addEvent(event)) throw std::runtime_error("BDX optical HIPO record insertion failed");
    writer.writeRecord(record);
}

}
#endif
