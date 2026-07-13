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
my $alv_z = 25./2;
my $alv_t = 0.2;

# Alveolus Y sized to fit the largest crystal face (type 4 has the largest H).
# Crystal is centred at y=0 in the alveolus; 0.1 cm gap between adjacent layer faces.
my $H_max     = $bgo_types{4}{H};   # 3.325 cm half-height of widest face
my $layer_gap = 0.1;                # gap between crystal faces of neighbouring layers
my $alv_in_y  = $H_max + $layer_gap / 2;   # 3.375 cm
my $alv_y     = $alv_in_y + $alv_t;        # 3.575 cm  (layer step = 2*alv_y = 7.15 cm)

my $alv_in_x = $alv_x - $alv_t;
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


sub make_65431_crs{
    my %detector = init_det();

    my $ix = $_[1];
    my $iy = $_[2];   # layer index 0..7
    my $iz = $_[3];

    # Only crystal types used in this geometry
    my @types = (6, 5, 4, 3, 1);
    my @Bv    = map { $bgo_types{$_}{B} } @types;

    my $gap  = 0.01;
    my @Xpos = (0);
    for my $i (1..$#types) {
        push @Xpos, $Xpos[-1] + $Bv[$i-1] + $gap + $Bv[$i];
    }
    my $Xshift = ($Xpos[0] + $Xpos[-1]) / 2;
    @Xpos = map { $_ - $Xshift } @Xpos;

    my $large_alv_x    = $Xpos[-1] + $Bv[-1] + $alv_t;
    my $large_alv_in_x = $large_alv_x - $alv_t;

    # Y centre of this layer: 8 layers centred around 0
    my $Y_layer = ($iy - 3.5) * 2 * $alv_y;

    %detector = init_det();
    $detector{"name"}        = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $Y_layer*cm 0*cm";
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

    # Loop over crystal types to create each crystal and its cap
    for my $ti (0..$#types) {
        my $type  = $types[$ti];
        my $Xcrs  = $Xpos[$ti];
        my $b     = $bgo_types{$type}{b};
        my $B     = $bgo_types{$type}{B};
        my $h     = $bgo_types{$type}{h};
        my $H     = $bgo_types{$type}{H};
        my $L     = $bgo_types{$type}{L};
        my $theta = $bgo_types{$type}{theta};

        my $is_even  = ($iy % 2 == 0) ? 1 : 0;
        my $rot      = $is_even ? "0*deg 0*deg 0*deg" : "0*deg 180*deg 0*deg";
        my $cap_sign = $is_even ? 1 : -1;

        my $cY   = 0;    # crystal centred at alveolus y=0
        my $cZ   = 0;
        my $capY = $cY;
        my $capZ = $cZ + $cap_sign * ($L+$d)*cos($theta * pi /180.0);

        %detector = init_det();
        $detector{"name"}        = "crs${type}_${ix}_${iy}_${iz}";
        $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
        $detector{"description"} = "crs type$type";
        $detector{"color"}       = "00ffff";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Trd";
        $detector{"pos"}         = "$Xcrs*cm $cY*cm $cZ*cm";
        $detector{"rotation"}    = $rot;
        $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
        $detector{"material"}    = "G4_BGO";
        $detector{"sensitivity"} = "crs";
        $detector{"hit_type"}    = "crs";
        $detector{"identifiers"} = "sector manual $type xch manual $ix ych manual $iy zch manual $iz SiPM manual 6025";
        print_det(\%configuration, \%detector);

        %detector = init_det();
        $detector{"name"}        = "cap_crs${type}_${ix}_${iy}_${iz}";
        $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
        $detector{"description"} = "crs type$type cap";
        $detector{"color"}       = "555555";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "$Xcrs*cm $capY*cm $capZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);

        %detector = init_det();
        $detector{"name"}        = "cap_hole_crs${type}_${ix}_${iy}_${iz}";
        $detector{"mother"}      = "cap_crs${type}_${ix}_${iy}_${iz}";
        $detector{"description"} = "crs type$type cap hole";
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
}

sub make_87722_crs{
    my %detector = init_det();

    my $ix = $_[1];
    my $iy = $_[2];   # layer index (8 or 9, placed above the 0..7 stack)
    my $iz = $_[3];

    # Crystal types for this row — includes duplicates, so $ti is used in names
    my @types = (8, 7, 7, 2, 2);
    my @Bv    = map { $bgo_types{$_}{B} } @types;

    my $gap  = 0.01;
    my @Xpos = (0);
    for my $i (1..$#types) {
        push @Xpos, $Xpos[-1] + $Bv[$i-1] + $gap + $Bv[$i];
    }
    my $Xshift = ($Xpos[0] + $Xpos[-1]) / 2;
    @Xpos = map { $_ - $Xshift } @Xpos;

    my $large_alv_x    = $Xpos[-1] + $Bv[-1] + $alv_t;
    my $large_alv_in_x = $large_alv_x - $alv_t;

    # Same Y_layer formula as make_8_crs — iy=8,9 naturally extend above the stack
    my $Y_layer = ($iy - 3.5) * 2 * $alv_y;

    %detector = init_det();
    $detector{"name"}        = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $Y_layer*cm 0*cm";
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

    for my $ti (0..$#types) {
        my $type  = $types[$ti];
        my $Xcrs  = $Xpos[$ti];
        my $b     = $bgo_types{$type}{b};
        my $B     = $bgo_types{$type}{B};
        my $h     = $bgo_types{$type}{h};
        my $H     = $bgo_types{$type}{H};
        my $L     = $bgo_types{$type}{L};
        my $theta = $bgo_types{$type}{theta};

        my $is_even  = ($iy % 2 == 0) ? 1 : 0;
        my $rot      = $is_even ? "${theta}*deg 0*deg 0*deg" : "-${theta}*deg 180*deg 0*deg";
        my $cap_sign = $is_even ? 1 : -1;

        my $cY   = 0;
        my $cZ   = 0;
        my $capY = $cY + ($L+$d)*sin($theta * pi /180.0);
        my $capZ = $cZ + $cap_sign * ($L+$d)*cos($theta * pi /180.0);

        # Include $ti in name to handle duplicate type numbers within the same row
        %detector = init_det();
        $detector{"name"}        = "crs${type}_t${ti}_${ix}_${iy}_${iz}";
        $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
        $detector{"description"} = "crs type$type pos$ti";
        $detector{"color"}       = "00ffff";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Trd";
        $detector{"pos"}         = "$Xcrs*cm $cY*cm $cZ*cm";
        $detector{"rotation"}    = $rot;
        $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
        $detector{"material"}    = "G4_BGO";
        $detector{"sensitivity"} = "crs";
        $detector{"hit_type"}    = "crs";
        $detector{"identifiers"} = "sector manual $type xch manual $ix ych manual $iy zch manual $iz SiPM manual 6025";
        print_det(\%configuration, \%detector);

        %detector = init_det();
        $detector{"name"}        = "cap_crs${type}_t${ti}_${ix}_${iy}_${iz}";
        $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
        $detector{"description"} = "crs type$type pos$ti cap";
        $detector{"color"}       = "555555";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "$Xcrs*cm $capY*cm $capZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);

        %detector = init_det();
        $detector{"name"}        = "cap_hole_crs${type}_t${ti}_${ix}_${iy}_${iz}";
        $detector{"mother"}      = "cap_crs${type}_t${ti}_${ix}_${iy}_${iz}";
        $detector{"description"} = "crs type$type pos$ti cap hole";
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
    $detector{"dimensions"}  = "30*cm 65*cm 30*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
    
    # Bottom edge of lowest crystal layer (iy=0): (0-3.5)*2*alv_y - alv_y = -8*alv_y
    my $crystal_stack_bot = -8 * $alv_y;
    my $scint_ht    = 0.5;                              # half-thickness of scint
    my $scint_bot_Y = $crystal_stack_bot - 2 - $scint_ht;   # 2 cm gap below bottom layer
    my $scint_top_Y = $scint_bot_Y + 80;                # 80 cm above bottom scint

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

    for my $layer (0..7) { make_65431_crs("1", 0, $layer, 0); }

    # Two additional rows of (8,7,7,2,2) above the main stack
    for my $layer (8..9) { make_87722_crs("1", 0, $layer, 0); }

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
    $detector{"dimensions"}  = "15*cm 0.5*cm 20*cm";
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
