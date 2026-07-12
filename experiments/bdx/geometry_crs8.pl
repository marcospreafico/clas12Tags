use strict;
use warnings;

use Math::Trig;

#TO DO: change veto identifiers to use sector to discriminate between IV and OV
our %configuration;

my $degrad = 0.01745329252;
my $cic    = 2.54;

#BEGIN general geometry parameters

my $d = 1./2; #thickness of aluminum crs cap
# BEGIN BGO parameters
my %bgo_types = (
                 "1" => {
                 b => 2.38/2.,
                 B => 4.82/2.,
                 h => 2.44/2.,
                 H => 4.96/2.,
                 L => 24./2.,
                 theta => 6./2,
                 },
"2" => {
    b => 2.34/2.,
    B => 5.26/2.,
    h => 2.38/2.,
    H => 5.31/2.,
    L => 24./2.,
    theta => 7./2,
},
"3" => {
    b => 2.27/2.,
    B => 5.64/2.,
    h => 2.13/2.,
    H => 5.27/2.,
    L => 24./2.,
    theta => 7.5/2,
},
"4" => {
    b => 2.23/2.,
    B => 6.17/2.,
    h => 2.45/2.,
    H => 6.65/2.,
    L => 24./2.,
    theta => 10./2,
},
"5" => {
    b => 2.16/2.,
    B => 6.47/2.,
    h => 2.08/2.,
    H => 6.28/2.,
    L => 24./2.,
    theta => 10./2,
},
"6" => {
    b => 2.08/2.,
    B => 6.67/2.,
    h => 1.88/2.,
    H => 6.08/2.,
    L => 24./2.,
    theta => 10./2,
},
"7" => {
    b => 2.01/2.,
    B => 6.74/2.,
    h => 1.74/2.,
    H => 5.80/2.,
    L => 24./2.,
    theta => 9.667/2,
},
"8" => {
    b => 1.96/2.,
    B => 6.68/2.,
    h => 1.69/2.,
    H => 5.75/2.,
    L => 24./2.,
    theta => 9.667/2,
},
);


my $alv_x = 7.4/2;
my $alv_y = 9.4/2;
my $alv_z = 26;

my $alv_t = 0.2;

my $alv_in_x = $alv_x-$alv_t;
my $alv_in_y = $alv_y-$alv_t;
my $alv_in_z = $alv_z;

sub make_main
{
    my %detector = init_det();
    $detector{"name"}        = "main_volume";
    $detector{"mother"}      = "root";
    $detector{"description"} = "World";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    
    my $X = 0.;
    my $Y = 0.;
    my $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    my $par1 = 800;
    my $par2 = 400.;
    my $par3 = 400.;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
}


sub make_8_crs{
    my %detector = init_det();

    my $ix = $_[1];
    my $iy = $_[2];
    my $iz = $_[3];

    my $X = 0;
    my $Y = 0;
    my $Z = 0;

    # Compute crystal X positions: 0.1 cm gap between adjacent crystal surfaces
    my $gap = 0.01;
    my @Bv  = map { $bgo_types{$_}{B} } (1..8);

    my @Xpos = (0);
    for my $i (1..7) {
        push @Xpos, $Xpos[-1] + $Bv[$i-1] + $gap + $Bv[$i];
    }
    my $Xshift = ($Xpos[0] + $Xpos[-1]) / 2;
    @Xpos = map { $_ - $Xshift } @Xpos;

    # Alveolus half-width: rightmost edge + wall thickness + 0.5 cm margin
    my $large_alv_x    = $Xpos[-1] + $Bv[-1] + $alv_t + 0.5;
    my $large_alv_in_x = $large_alv_x - $alv_t;

    %detector = init_det();
    $detector{"name"}        = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$large_alv_x*cm $alv_y*cm $alv_z*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "alveolus air";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$large_alv_in_x*cm $alv_in_y*cm $alv_in_z*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    # ============================ Crs type1 ============================
    my $Xcrs  = $Xpos[0];
    my $b     = $bgo_types{"1"}{b};
    my $B     = $bgo_types{"1"}{B};
    my $h     = $bgo_types{"1"}{h};
    my $H     = $bgo_types{"1"}{H};
    my $L     = $bgo_types{"1"}{L};
    my $theta = $bgo_types{"1"}{theta};
    $Y = -$H/2.-$h/2.;
    $Z = 0;

    %detector = init_det();
    $detector{"name"}        = "crs1_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs1 type1";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
    $detector{"material"}    = "G4_BGO";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 1 xch manual $ix ych manual $iy zch manual $iz SiPM manual 6025";
    print_det(\%configuration, \%detector);

    $Y = $Y+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z+($L+$d)*cos($theta * pi /180.0);

    %detector = init_det();
    $detector{"name"}        = "cap_crs1_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs1 type1 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "cap_hole_crs1_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cap_crs1_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs1 type1 cap hole";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    # ============================ Crs type2 ============================
    $Xcrs  = $Xpos[1];
    $b     = $bgo_types{"2"}{b};
    $B     = $bgo_types{"2"}{B};
    $h     = $bgo_types{"2"}{h};
    $H     = $bgo_types{"2"}{H};
    $L     = $bgo_types{"2"}{L};
    $theta = $bgo_types{"2"}{theta};
    $Y = -$H/2.-$h/2.;
    $Z = 0;

    %detector = init_det();
    $detector{"name"}        = "crs2_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs2 type2";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "-$theta*deg 180*deg 0*deg";
    $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
    $detector{"material"}    = "G4_BGO";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 2 xch manual $ix ych manual $iy zch manual $iz SiPM manual 6025";
    print_det(\%configuration, \%detector);

    $Y = $Y+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z-($L+$d)*cos($theta * pi /180.0);

    %detector = init_det();
    $detector{"name"}        = "cap_crs2_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs2 type2 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "-$theta*deg 180*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "cap_hole_crs2_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cap_crs2_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs2 type2 cap hole";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    # ============================ Crs type3 ============================
    $Xcrs  = $Xpos[2];
    $b     = $bgo_types{"3"}{b};
    $B     = $bgo_types{"3"}{B};
    $h     = $bgo_types{"3"}{h};
    $H     = $bgo_types{"3"}{H};
    $L     = $bgo_types{"3"}{L};
    $theta = $bgo_types{"3"}{theta};
    $Y = -$H/2.-$h/2.;
    $Z = 0;

    %detector = init_det();
    $detector{"name"}        = "crs3_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs3 type3";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
    $detector{"material"}    = "G4_BGO";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 3 xch manual $ix ych manual $iy zch manual $iz SiPM manual 6025";
    print_det(\%configuration, \%detector);

    $Y = $Y+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z+($L+$d)*cos($theta * pi /180.0);

    %detector = init_det();
    $detector{"name"}        = "cap_crs3_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs3 type3 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "cap_hole_crs3_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cap_crs3_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs3 type3 cap hole";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    # ============================ Crs type4 ============================
    $Xcrs  = $Xpos[3];
    $b     = $bgo_types{"4"}{b};
    $B     = $bgo_types{"4"}{B};
    $h     = $bgo_types{"4"}{h};
    $H     = $bgo_types{"4"}{H};
    $L     = $bgo_types{"4"}{L};
    $theta = $bgo_types{"4"}{theta};
    $Y = -$H/2.-$h/2.;
    $Z = 0;

    %detector = init_det();
    $detector{"name"}        = "crs4_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs4 type4";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "-$theta*deg 180*deg 0*deg";
    $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
    $detector{"material"}    = "G4_BGO";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 4 xch manual $ix ych manual $iy zch manual $iz SiPM manual 6025";
    print_det(\%configuration, \%detector);

    $Y = $Y+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z-($L+$d)*cos($theta * pi /180.0);

    %detector = init_det();
    $detector{"name"}        = "cap_crs4_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs4 type4 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "-$theta*deg 180*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "cap_hole_crs4_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cap_crs4_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs4 type4 cap hole";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    # ============================ Crs type5 ============================
    $Xcrs  = $Xpos[4];
    $b     = $bgo_types{"5"}{b};
    $B     = $bgo_types{"5"}{B};
    $h     = $bgo_types{"5"}{h};
    $H     = $bgo_types{"5"}{H};
    $L     = $bgo_types{"5"}{L};
    $theta = $bgo_types{"5"}{theta};
    $Y = -$H/2.-$h/2.;
    $Z = 0;

    %detector = init_det();
    $detector{"name"}        = "crs5_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs5 type5";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
    $detector{"material"}    = "G4_BGO";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 5 xch manual $ix ych manual $iy zch manual $iz SiPM manual 6025";
    print_det(\%configuration, \%detector);

    $Y = $Y+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z+($L+$d)*cos($theta * pi /180.0);

    %detector = init_det();
    $detector{"name"}        = "cap_crs5_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs5 type5 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "cap_hole_crs5_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cap_crs5_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs5 type5 cap hole";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    # ============================ Crs type6 ============================
    $Xcrs  = $Xpos[5];
    $b     = $bgo_types{"6"}{b};
    $B     = $bgo_types{"6"}{B};
    $h     = $bgo_types{"6"}{h};
    $H     = $bgo_types{"6"}{H};
    $L     = $bgo_types{"6"}{L};
    $theta = $bgo_types{"6"}{theta};
    $Y = -$H/2.-$h/2.;
    $Z = 0;

    %detector = init_det();
    $detector{"name"}        = "crs6_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs6 type6";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "-$theta*deg 180*deg 0*deg";
    $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
    $detector{"material"}    = "G4_BGO";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 6 xch manual $ix ych manual $iy zch manual $iz SiPM manual 6025";
    print_det(\%configuration, \%detector);

    $Y = $Y+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z-($L+$d)*cos($theta * pi /180.0);

    %detector = init_det();
    $detector{"name"}        = "cap_crs6_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs6 type6 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "-$theta*deg 180*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "cap_hole_crs6_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cap_crs6_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs6 type6 cap hole";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    # ============================ Crs type7 ============================
    $Xcrs  = $Xpos[6];
    $b     = $bgo_types{"7"}{b};
    $B     = $bgo_types{"7"}{B};
    $h     = $bgo_types{"7"}{h};
    $H     = $bgo_types{"7"}{H};
    $L     = $bgo_types{"7"}{L};
    $theta = $bgo_types{"7"}{theta};
    $Y = -$H/2.-$h/2.;
    $Z = 0;

    %detector = init_det();
    $detector{"name"}        = "crs7_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs7 type7";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
    $detector{"material"}    = "G4_BGO";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 7 xch manual $ix ych manual $iy zch manual $iz SiPM manual 6025";
    print_det(\%configuration, \%detector);

    $Y = $Y+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z+($L+$d)*cos($theta * pi /180.0);

    %detector = init_det();
    $detector{"name"}        = "cap_crs7_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs7 type7 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "cap_hole_crs7_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cap_crs7_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs7 type7 cap hole";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    # ============================ Crs type8 ============================
    $Xcrs  = $Xpos[7];
    $b     = $bgo_types{"8"}{b};
    $B     = $bgo_types{"8"}{B};
    $h     = $bgo_types{"8"}{h};
    $H     = $bgo_types{"8"}{H};
    $L     = $bgo_types{"8"}{L};
    $theta = $bgo_types{"8"}{theta};
    $Y = -$H/2.-$h/2.;
    $Z = 0;

    %detector = init_det();
    $detector{"name"}        = "crs8_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs8 type8";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "-$theta*deg 180*deg 0*deg";
    $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
    $detector{"material"}    = "G4_BGO";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 8 xch manual $ix ych manual $iy zch manual $iz SiPM manual 6025";
    print_det(\%configuration, \%detector);

    $Y = $Y+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z-($L+$d)*cos($theta * pi /180.0);

    %detector = init_det();
    $detector{"name"}        = "cap_crs8_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs8 type8 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$Xcrs*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "-$theta*deg 180*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "cap_hole_crs8_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cap_crs8_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs8 type8 cap hole";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

}

sub make_detector{
    my %detector = init_det();
    $detector{"name"}        = "detector_volume";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "detector_volume";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    my $X = 0;
    my $Y = 0;
    my $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "30*cm 30*cm 30*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
    
    # scint_down: top surface flush with alveolus bottom
    my $scint_ht    = 0.5;
    my $scint_bot_Y = -($alv_y + $scint_ht);   
    my $scint_top_Y = $scint_bot_Y + 13;        

    # 4 top scints equally spaced 3 cm apart in Z, centered on the array
    my @Ztop = (-9, -3, 3, 9);

    %detector = init_det();
    $detector{"name"}        = "scint_1";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic scint 1";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $scint_top_Y*cm $Ztop[0]*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "1*cm 1*cm 28*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 1 xch manual 0 ych manual 0 zch manual 0";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "scint_2";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic scint 2";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $scint_top_Y*cm $Ztop[1]*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "1*cm 1*cm 28*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 2 xch manual 0 ych manual 0 zch manual 0";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "scint_3";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic scint 3";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $scint_top_Y*cm $Ztop[2]*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "1*cm 1*cm 28*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 3 xch manual 0 ych manual 0 zch manual 0";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "scint_4";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic scint 4";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $scint_top_Y*cm $Ztop[3]*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "1*cm 1*cm 28*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 4 xch manual 0 ych manual 0 zch manual 0";
    print_det(\%configuration, \%detector);

    make_8_crs("1", 0, 0, 0);

    # Scintillator under the crystal array
    %detector = init_det();
    $detector{"name"}        = "scint_down";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic scint down";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $scint_bot_Y*cm 0*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "15*cm 0.5*cm 28*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 5 xch manual 0 ych manual 0 zch manual 0";
    print_det(\%configuration, \%detector);

}

sub make_all
{
    make_main;
    make_detector;
}


1;
