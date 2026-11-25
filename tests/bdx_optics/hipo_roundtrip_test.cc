#include "BdxOpticalSchema.h"
#include "hipo4/writer.h"
#include "hipo4/reader.h"
#include <stdexcept>
#include <iostream>
int main(int argc,char** argv) {
    if(argc!=2) return 2;
    hipo::writer writer; bdxoptics::AddOpticalSchemas(writer.getDictionary());
    hipo::schema runSchema("RUN::config",10000,11);runSchema.parse("event/I");writer.getDictionary().addSchema(runSchema);
    hipo::schema mcSchema("MC::Particle",40,2);mcSchema.parse("pid/I");writer.getDictionary().addSchema(mcSchema);
    writer.open(argv[1]);
    const int photons=100000;
    for(int index=0;index<2;++index) {
        const int count=index==0?photons:0;
        hipo::event event(count*180+4096);
        hipo::bank run(runSchema,1);run.putInt("event",0,index+1);event.addStructure(run);
        hipo::bank mc(mcSchema,1);mc.putInt("pid",0,13);event.addStructure(mc);
        hipo::bank c(writer.getDictionary().getSchema("BDX::Config"),1);
        c.putInt("schema_version",0,2);c.putInt("collector",0,1);
        for(auto key:{"H_mm","Bg_mm","Bm_mm"}) c.putDouble(key,0,10);
        c.putDouble("plane_z_mm",0,120);c.putDouble("absorption_length_mm",0,10000);event.addStructure(c);
        hipo::bank e(writer.getDictionary().getSchema("BDX::Event"),1);
        e.putInt("event",0,index+1);e.putInt("event_index",0,index);
        e.putDouble("edep_MeV",0,1);e.putDouble("edep_axial_distance_mm",0,120);
        e.putInt("generated_scint",0,index==0?photons:1);e.putInt("arrived_scint",0,count);
        e.putInt("bulk_absorbed",0,index==0?0:1);event.addStructure(e);
        hipo::bank d(writer.getDictionary().getSchema("BDX::Edep"),1);
        d.putInt("source_step",0,0);d.putInt("generated_scint",0,index==0?photons:1);d.putDouble("edep_MeV",0,1);event.addStructure(d);
        hipo::bank f(writer.getDictionary().getSchema("BDX::Flux"),count);
        for(int i=0;i<count;++i) {
            f.putInt("track",i,i+2);f.putInt("parent",i,1);f.putInt("origin",i,1);f.putInt("source_step",i,0);
            f.putDouble("x_mm",i,i%2?5:0);f.putDouble("lx_mm",i,i%2?5:0);
            f.putDouble("z_mm",i,120);f.putDouble("uz",i,1);f.putDouble("path_mm",i,120);f.putDouble("time_ns",i,1);
            f.putDouble("energy_eV",i,1239.841984332/480);f.putDouble("wavelength_nm",i,480);f.putDouble("weight",i,1);
        }
        event.addStructure(f);
        if(index==0 && event.getSize()<8*1024*1024) throw std::runtime_error("test event too small");
        bdxoptics::WriteOpticalRecord(writer,event);
        std::cout << "Wrote test event of " << event.getSize() << " bytes\n";
    }
    writer.close();
    return 0;
}
