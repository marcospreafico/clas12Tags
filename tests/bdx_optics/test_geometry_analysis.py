import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from bdx_geometry import load_types, geometry_text, generate, SOURCE
from analyze_bdx_optics import summarize
from run_bdx_optics import config_values, validate_events


class GeometryTests(unittest.TestCase):
    def test_exact_match_to_executed_upstream_perl(self):
        # Execute the upstream geometry routines with a print_det capture; this
        # independently checks extraction, units, vertex ordering and the active
        # branches (including the asymmetric trapezoid faces).
        program=r'''
use strict; use warnings; use JSON::PP;
our %configuration=(large_crs=>$ARGV[1],small_crs=>$ARGV[1]);
our @saved;
sub init_det { return (); }
sub print_det { my ($config,$detector)=@_; push @saved,{%$detector} if $detector->{name} eq "crs_0_0"; }
do $ARGV[0] or die $@ || $!;
make_large_col_crs(0);
print encode_json($saved[0]);
'''
        for kind,shape in load_types().items():
            detector=json.loads(subprocess.check_output(['perl','-e',program,str(SOURCE),str(kind)],text=True))
            self.assertEqual(detector['type'],'G4GenericTrap')
            upstream=[float(x.removesuffix('*cm'))*10 for x in detector['dimensions'].split()]
            row=geometry_text(shape).splitlines()[0].split('|')
            ours=[float(x.removesuffix('*mm')) for x in row[7].split()]
            self.assertEqual(len(ours),17)
            for x,y in zip(upstream,ours): self.assertAlmostEqual(x,y,places=10)

    def test_centroid_and_sensor_contact(self):
        types=load_types()
        self.assertEqual(len(types),8)
        self.assertEqual(types[8]['sensor_center_mm'],[0,0,120])
        for shape in types.values():
            self.assertEqual(shape['half_length_mm'],120)
            rows=[row.split('|') for row in geometry_text(shape).splitlines()]
            self.assertTrue(all(len(row)==18 for row in rows))
            z=float(rows[1][3].split()[2].removesuffix('*mm'))
            self.assertAlmostEqual(z-0.05,120)

    def test_regeneration_and_full_base(self):
        with tempfile.TemporaryDirectory() as path:
            generate(Path(path))
            self.assertEqual(len(list(Path(path).glob('*.gcard'))),8)
        for shape in load_types().values():
            row=geometry_text(shape).splitlines()[1].split('|')
            self.assertEqual(row[0].strip(),'bdx_flux')
            self.assertEqual(row[6].strip(),'G4GenericTrap')
            dims=[float(x.removesuffix('*mm')) for x in row[7].split()]
            self.assertEqual(dims[0],0.05)
            face=sum(shape['vertices_mm'][4:],[])
            for x,y in zip(dims[1:],face+face): self.assertAlmostEqual(x,y,places=10)


class AnalysisTests(unittest.TestCase):
    def test_ratio_of_totals_and_single_event(self):
        rows=[dict(generated_scint=100,arrived_scint=10,edep_MeV=1,
                   edep_axial_distance_mm=50,max_optical_steps=150),
              dict(generated_scint=300,arrived_scint=60,edep_MeV=3,
                   edep_axial_distance_mm=60,max_optical_steps=200)]
        result=summarize(rows,100)
        self.assertAlmostEqual(result['collection_efficiency'],70/400)
        self.assertAlmostEqual(result['arrived_per_MeV'],17.5)
        self.assertAlmostEqual(result['edep_axial_distance_mm'],57.5)
        self.assertGreater(result['efficiency_event_sem'],0)
        self.assertIsNone(summarize(rows[:1],100)['efficiency_event_sem'])

    def test_zero_arrivals_has_nonzero_upper_interval(self):
        rows=[dict(generated_scint=1000,arrived_scint=0,edep_MeV=1,
                   edep_axial_distance_mm=50,max_optical_steps=1)]
        result=summarize(rows,1000)
        self.assertGreater(result['efficiency_photon_wilson68_high'],0)
        self.assertEqual(result['collection_efficiency'],0)

    def test_input_validation(self):
        config=ROOT/'experiments/bdx_optics/geometry_only.conf'
        for override in ('light_yield_per_MeV=0','wrap_reflectivity=1.1',
                         'absorption_length_mm=-5','unknown=5','refractive_index=nan'):
            with self.assertRaises(ValueError): config_values(config,[override])


if __name__=='__main__': unittest.main()
