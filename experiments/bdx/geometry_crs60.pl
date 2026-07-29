use strict;
use warnings;

use Math::Trig;
# New 60crs geometry using column loops rather than rows.
# Proper placement using 7x7 and 6x6 masks to get accurate dimensions.
# No Alveolus volume in this version. Use masks to define the crystal placement.


#TO DO: change veto identifiers to use sector to discriminate between IV and OV
our %configuration;

my $degrad = 0.01745329252;
my $cic    = 2.54;

my $d = 1./2;  # half-thickness of Al cap

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


my $layer_gap        = 0.1;    # clearance between adjacent layer faces (cm)
my $col_gap          = 0.25;   # face-to-face gap between columns — support rods (cm)
my $Z_shift          = -0.5;   # global Z offset for all crystal volumes (cm)
my $mask_large_half  = 3.5;    # Y half-size of 7 cm large mask (types 4-8)
my $mask_small_half  = 5.8/2.; # Y half-size of 5.8 cm small mask (types 1-3)
my $cap_large_half   = 3.5;    # Y half-height of large plastic mask
my $cap_large_half_x = 3.4;    # X half-width  of large plastic mask
my $cap_small_half   = 5.8/2.; # half-dimension of small (6x6 cm) plastic mask
my $layer_half_y     = $mask_large_half + $layer_gap / 2;  # 3.55 cm

# Row centre-to-centre pitch: pitch = mask_half + mask_gap/2
my $mask_gap          = 2.98;
my $large_layer_pitch = $mask_large_half  + $mask_gap / 2;  # 5.0 cm
my $small_layer_pitch = $mask_small_half  + $mask_gap / 2;  # ~4.4 cm

# Y shift so the bottom edge of small-col row 0 aligns with large-col row 0.
my $Y_small_shift = 5.5 * ($small_layer_pitch - $large_layer_pitch)
                  + $layer_half_y - $mask_small_half
                  - ($cap_large_half - $cap_small_half);

# Fixed column X centres: c-to-c = hw[i] + col_gap + hw[i+1], then centred.
my @col_mask_hw = ($cap_large_half_x, $cap_large_half_x, $cap_large_half_x,
                   $cap_small_half,   $cap_small_half);
my @col_centres = (0);
for my $i (1..4) {
    push @col_centres, $col_centres[-1] + $col_mask_hw[$i-1] + $col_gap + $col_mask_hw[$i];
}
my $Xcshift = ($col_centres[0] + $col_centres[-1]) / 2;
@col_centres = map { $_ - $Xcshift } @col_centres;

# Full crystal serial number table: $crs_id{$iy}[$ti], ti=0=leftmost, ti=N=rightmost
# iy=0 bottom row, iy=10 top row
my %crs_id = (
    #  iy  |  col=0    col=1    col=2    col=3    col=4   | row type
     0 => [  604,   1116,    407,    311,   1532],  # 65431
     1 => [  627,    507,    418,    307,    124],  # 65431
     2 => [  625,    504,    430,    330,   1513],  # 65431
     3 => [  607,    503,    419,    305,   1521],  # 65431
     4 => [ 1002,   1119,    415,   1303,    111],  # 65431
     5 => [ 1009,    533,   1228,    327,   1524],  # 65431
     6 => [ 1026,   1102,    433,    322,   1531],  # 65431
     7 => [ 1001,   1109,   1212,   1421,    132],  # 65421
     8 => [  825,    916,    911,    220,   1420],  # 87722
     9 => [  819,    912,    733,   1425,   1414],  # 87722
    10 => [  812,    924,    926,    217,    212],  # 87722
);

# Per-row crystal type layouts — drives X-position computation
my @row_layouts = (
    [6, 5, 4, 3, 1],  # iy=0
    [6, 5, 4, 3, 1],  # iy=1
    [6, 5, 4, 3, 1],  # iy=2
    [6, 5, 4, 3, 1],  # iy=3
    [6, 5, 4, 3, 1],  # iy=4
    [6, 5, 4, 3, 1],  # iy=5
    [6, 5, 4, 3, 1],  # iy=6
    [6, 5, 4, 2, 1],  # iy=7  (65421)
    [8, 7, 7, 2, 2],  # iy=8  (87722)
    [8, 7, 7, 2, 2],  # iy=9
    [8, 7, 7, 2, 2],  # iy=10
);

sub make_main{    my %detector = init_det();
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

# Per crystal: BGO Trd → Al cap → air hole (child of cap) → plastic cap mask → plastic back mask.
# Caps and crystals follow $Y_crs; masks follow $Y_layer so row pitch is uniform.

# --- Large crystal columns (types 4-8: 7x7 cm masks) ---
sub make_large_col_crs {
    my $ti     = shift;
    my $mask_t = 0.5;  # half-thickness of plastic masks (full = 1 cm)

    for my $iy (0..10) {
        my @row = @{ $row_layouts[$iy] };

        my $type  = $row[$ti];
        my $Xcrs  = $col_centres[$ti];
        my $b     = $bgo_types{$type}{b};
        my $B     = $bgo_types{$type}{B};
        my $h     = $bgo_types{$type}{h};
        my $H     = $bgo_types{$type}{H};
        my $L     = $bgo_types{$type}{L};

        my $Y_layer  = ($iy - 5.5) * $large_layer_pitch;
        my $Y_crs    = $Y_layer + (-$layer_half_y + $mask_large_half);

        my $is_even  = ($iy % 2 == 0) ? 1 : 0;
        my $rot      = $is_even ? "0*deg 0*deg 0*deg" : "0*deg 180*deg 0*deg";
        my $cap_sign = $is_even ? 1 : -1;
        my $capZ     = $cap_sign * ($L + $d)                    + $cap_sign * $Z_shift;
        my $cmaskZ   = $cap_sign * ($L + 2*$d + $mask_t + 0.01) + $cap_sign * $Z_shift;
        my $backZ    = -$cap_sign * ($L + $mask_t + 0.035)      + $cap_sign * $Z_shift;

        my $serial = $crs_id{$iy}[$ti] // 0;

        # 1. BGO crystal
        my %detector = init_det();
        $detector{"name"}        = "crs_${ti}_${iy}";
        $detector{"mother"}      = "detector_volume";
        $detector{"description"} = "crs col$ti row$iy type$type";
        $detector{"color"}       = "00ffff";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Trd";
        $detector{"pos"}         = "$Xcrs*cm $Y_crs*cm " . ($cap_sign * $Z_shift) . "*cm";
        $detector{"rotation"}    = $rot;
        $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
        $detector{"material"}    = "G4_BGO";
        $detector{"sensitivity"} = "crs";
        $detector{"hit_type"}    = "crs";
        $detector{"identifiers"} = "sector manual $type xch manual $serial ych manual $iy zch manual $ti SiPM manual 6025";
        print_det(\%configuration, \%detector);

        # 2. Cap — follows crystal, dimensions B x H
        %detector = init_det();
        $detector{"name"}        = "cap_${ti}_${iy}";
        $detector{"mother"}      = "detector_volume";
        $detector{"description"} = "Al cap col$ti row$iy";
        $detector{"color"}       = "555555";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "$Xcrs*cm $Y_crs*cm $capZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);

        # 3. Hole through cap
        %detector = init_det();
        $detector{"name"}        = "cap_hole_${ti}_${iy}";
        $detector{"mother"}      = "cap_${ti}_${iy}";
        $detector{"description"} = "cap hole col$ti row$iy";
        $detector{"color"}       = "ffffff";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "0*cm 0*cm 0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
        $detector{"material"}    = "G4_AIR";
        print_det(\%configuration, \%detector);

        # 4. Plastic cap mask 7x7 cm — centred on ROW ($Y_layer)
        %detector = init_det();
        $detector{"name"}        = "cmask_${ti}_${iy}";
        $detector{"mother"}      = "detector_volume";
        $detector{"description"} = "plastic cap mask col$ti row$iy";
        $detector{"color"}       = "fffdd0";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "$Xcrs*cm $Y_layer*cm $cmaskZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$cap_large_half_x*cm $cap_large_half*cm $mask_t*cm";
        $detector{"material"}    = "G4_POLYSTYRENE";
        print_det(\%configuration, \%detector);

        # 5. Plastic back mask 6.8x3 cm — centred on ROW ($Y_layer)
        %detector = init_det();
        $detector{"name"}        = "bmask_${ti}_${iy}";
        $detector{"mother"}      = "detector_volume";
        $detector{"description"} = "plastic back mask col$ti row$iy";
        $detector{"color"}       = "f2d5aa";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "$Xcrs*cm $Y_layer*cm $backZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$cap_large_half_x*cm 1.475*cm $mask_t*cm";
        $detector{"material"}    = "G4_POLYSTYRENE";
        print_det(\%configuration, \%detector);
    }
}

# --- Small crystal columns (types 1-3: 6x6 cm masks) ---
sub make_small_col_crs {
    my $ti     = shift;
    my $mask_t = 0.5;  # half-thickness of plastic masks (full = 1 cm)

    for my $iy (0..10) {
        my @row = @{ $row_layouts[$iy] };

        my $type  = $row[$ti];
        my $Xcrs  = $col_centres[$ti];
        my $b     = $bgo_types{$type}{b};
        my $B     = $bgo_types{$type}{B};
        my $h     = $bgo_types{$type}{h};
        my $H     = $bgo_types{$type}{H};
        my $L     = $bgo_types{$type}{L};

        my $Y_layer  = ($iy - 5.5) * $small_layer_pitch;
        my $Y_crs    = $Y_layer + (-$layer_half_y + $mask_small_half) + $Y_small_shift;

        my $is_even  = ($iy % 2 == 0) ? 1 : 0;
        my $rot      = $is_even ? "0*deg 0*deg 0*deg" : "0*deg 180*deg 0*deg";
        my $cap_sign = $is_even ? 1 : -1;
        my $capZ     = $cap_sign * ($L + $d)                    + $cap_sign * $Z_shift;
        my $cmaskZ   = $cap_sign * ($L + 2*$d + $mask_t + 0.01) + $cap_sign * $Z_shift;
        my $backZ    = -$cap_sign * ($L + $mask_t + 0.035)      + $cap_sign * $Z_shift;

        my $serial = $crs_id{$iy}[$ti] // 0;

        # 1. BGO crystal
        my %detector = init_det();
        $detector{"name"}        = "crs_${ti}_${iy}";
        $detector{"mother"}      = "detector_volume";
        $detector{"description"} = "crs col$ti row$iy type$type";
        $detector{"color"}       = "00ffff";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Trd";
        $detector{"pos"}         = "$Xcrs*cm $Y_crs*cm " . ($cap_sign * $Z_shift) . "*cm";
        $detector{"rotation"}    = $rot;
        $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
        $detector{"material"}    = "G4_BGO";
        $detector{"sensitivity"} = "crs";
        $detector{"hit_type"}    = "crs";
        $detector{"identifiers"} = "sector manual $type xch manual $serial ych manual $iy zch manual $ti SiPM manual 6025";
        print_det(\%configuration, \%detector);

        # 2. Al cap — follows crystal, dimensions B x H
        %detector = init_det();
        $detector{"name"}        = "cap_${ti}_${iy}";
        $detector{"mother"}      = "detector_volume";
        $detector{"description"} = "Al cap col$ti row$iy";
        $detector{"color"}       = "555555";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "$Xcrs*cm $Y_crs*cm $capZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);

        # 3. Hole through Al cap
        %detector = init_det();
        $detector{"name"}        = "cap_hole_${ti}_${iy}";
        $detector{"mother"}      = "cap_${ti}_${iy}";
        $detector{"description"} = "cap hole col$ti row$iy";
        $detector{"color"}       = "ffffff";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "0*cm 0*cm 0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
        $detector{"material"}    = "G4_AIR";
        print_det(\%configuration, \%detector);

        # 4. Plastic cap mask 6x6 cm — centred on ROW ($Y_layer)
        %detector = init_det();
        $detector{"name"}        = "cmask_${ti}_${iy}";
        $detector{"mother"}      = "detector_volume";
        $detector{"description"} = "plastic cap mask col$ti row$iy";
        $detector{"color"}       = "fffdd0";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "$Xcrs*cm $Y_crs*cm $cmaskZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$cap_small_half*cm $cap_small_half*cm $mask_t*cm";
        $detector{"material"}    = "G4_POLYSTYRENE";
        print_det(\%configuration, \%detector);

        # 5. Plastic back mask 6x3 cm — centred on ROW ($Y_layer)
        %detector = init_det();
        $detector{"name"}        = "bmask_${ti}_${iy}";
        $detector{"mother"}      = "detector_volume";
        $detector{"description"} = "plastic back mask col$ti row$iy";
        $detector{"color"}       = "f2d5aa";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "$Xcrs*cm $Y_crs*cm $backZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$cap_small_half*cm 1.475*cm $mask_t*cm";
        $detector{"material"}    = "G4_POLYSTYRENE";
        print_det(\%configuration, \%detector);
    }
}

#Column col=0 (604,627,...). Large crystals, 7x7 cm masks.
sub make_col0_crs { make_large_col_crs(0);}

#Column col=1 (1116,507,...). Large crystals, 7x7 cm masks.
sub make_col1_crs { make_large_col_crs(1);}

#Column col=2 (407,418,...). Large crystals, 7x7 cm masks.
sub make_col2_crs { make_large_col_crs(2);}

#Column col=3 (311,307,...). Small crystals, 6x6 cm masks.
sub make_col3_crs { make_small_col_crs(3);}

#Column col=4 (1532,124,...). Small crystals, 6x6 cm masks.
sub make_col4_crs { make_small_col_crs(4);}

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
    $detector{"dimensions"}  = "30*cm 50*cm 30*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);

    my $crystal_stack_bot = (0 - 5.5) * $large_layer_pitch - $layer_half_y;
    my $scint_ht    = 1.0;   # scint half-thickness
    my $scint_bot_Y = $crystal_stack_bot - 0.5 - $scint_ht;  # 0.5 cm gap below crystal stack
    my $scint_top_Y = $scint_bot_Y + 60;                      # 60 cm above bottom scint

    # 4 top scints at fixed Z positions matching the 4 bottom scints
    my @Ztop = (-7.55, -3.05, 2.45, 7.45);
    my $top_length = 35./2.;

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
    $detector{"dimensions"}  = "1*cm 1*cm $top_length*cm";
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
    $detector{"dimensions"}  = "1*cm 1*cm $top_length*cm";
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
    $detector{"dimensions"}  = "1*cm 1*cm $top_length*cm";
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
    $detector{"dimensions"}  = "1*cm 1*cm $top_length*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 4 xch manual 0 ych manual 0 zch manual 0";
    print_det(\%configuration, \%detector);

    make_col0_crs();
    make_col1_crs();
    make_col2_crs();
    make_col3_crs();
    make_col4_crs();

    # 4 bottom scints — same Z positions as top
    my @Zbot = (-7.55, -3.05, 2.45, 7.45);
    my $bot_length = 45./2.;

    %detector = init_det();
    $detector{"name"}        = "scint_5";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic scint 5";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $scint_bot_Y*cm $Zbot[0]*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "1*cm 1*cm $bot_length*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 5 xch manual 0 ych manual 0 zch manual 0";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "scint_6";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic scint 6";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $scint_bot_Y*cm $Zbot[1]*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "1*cm 1*cm $bot_length*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 6 xch manual 0 ych manual 0 zch manual 0";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "scint_7";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic scint 7";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $scint_bot_Y*cm $Zbot[2]*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "1*cm 1*cm $bot_length*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 7 xch manual 0 ych manual 0 zch manual 0";
    print_det(\%configuration, \%detector);

    %detector = init_det();
    $detector{"name"}        = "scint_8";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic scint 8";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $scint_bot_Y*cm $Zbot[3]*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "1*cm 1*cm $bot_length*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 8 xch manual 0 ych manual 0 zch manual 0";
    print_det(\%configuration, \%detector);

}

sub make_all
{
    make_main;
    make_detector;
}


1;
