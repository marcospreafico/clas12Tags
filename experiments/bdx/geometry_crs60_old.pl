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

# Alveolus Y sized by the large-crystal mask (7x7 cm cap-end mask → half = 3.5 cm).
# Small crystals (types 1-3) use a 6x6 cm mask → placed 2.5 cm from alveolus bottom,
# giving a ~1 cm physical misalignment between large and small crystal rows.
my $layer_gap       = 0.1;   # small clearance between adjacent layer faces
my $mask_large_half = 3.5;   # half of 7 cm cap mask — drives layer pitch for types 4-8
my $mask_small_half = 2.5;   # placement offset from alveolus bottom for types 1-3
my $alv_in_y  = $mask_large_half + $layer_gap / 2;  # 3.55 cm
my $alv_y     = $alv_in_y + $alv_t;                 # 3.75 cm  (layer step = 2*alv_y = 7.5 cm)

my $alv_in_x = $alv_x - $alv_t;
my $alv_in_z = $alv_z;

# Centre-to-centre pitch between adjacent alveolus layers.
# Even rows (iy=0,2,4...) share the same cap orientation — the face-to-face distance
# in Y between their masks (half-size = mask_large_half) separated by 2 layers is:
#   mask_gap = 2*layer_pitch - 2*mask_large_half
# Solving: layer_pitch = mask_large_half + mask_gap/2
my $mask_gap    = 3.0;   # desired face-to-face gap in Y between same-orientation masks (cm)
my $layer_pitch = $mask_large_half + $mask_gap / 2;  # = 3.5 + 1.5 = 5.0 cm

# Full crystal serial number table: $crs_id{$iy}[$ti], ti=0=leftmost, ti=N=rightmost
# iy=0 bottom row, iy=10 top row
my %crs_id = (
    #  iy  |  ti=0    ti=1    ti=2    ti=3    ti=4   | row type
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

sub make_main{
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

    my $gap  = 0.05;
    my @Xpos = (0);
    for my $i (1..$#types) {
        push @Xpos, $Xpos[-1] + $Bv[$i-1] + $gap + $Bv[$i];
    }
    my $Xshift = ($Xpos[0] + $Xpos[-1]) / 2;
    @Xpos = map { $_ - $Xshift } @Xpos;

    my $large_alv_x    = $Xpos[-1] + $Bv[-1] + $alv_t;
    my $large_alv_in_x = $large_alv_x - $alv_t;

    # Y centre of this layer: layers centred around 0
    my $Y_layer = ($iy - 5.5) * $layer_pitch;

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

        # Large crystals (types 4-8): cap mask 7x7 cm → center 3.5 cm from alveolus bottom
        # Small crystals (types 1-3): cap mask 6x6 cm → center 2.5 cm from alveolus bottom (~1 cm lower)
        my $mask_offset = ($type >= 4) ? $mask_large_half : $mask_small_half;
        my $cY   = -$alv_in_y + $mask_offset;
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
        my $serial   = $crs_id{$iy}[$ti] // 0;
        $detector{"identifiers"} = "sector manual $type xch manual $serial ych manual $iy zch manual $ti SiPM manual 6025";
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

    # Same Y_layer formula — iy=8,9,10 extend above the main stack
    my $Y_layer = ($iy - 5.5) * $layer_pitch;

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
        my $rot      = $is_even ? "0*deg 0*deg 0*deg" : "0*deg 180*deg 0*deg";
        my $cap_sign = $is_even ? 1 : -1;

        # Large crystals (types 4-8): cap mask 7x7 cm → center 3.5 cm from alveolus bottom
        # Small crystals (types 1-3): cap mask 6x6 cm → center 2.5 cm from alveolus bottom (~1 cm lower)
        my $mask_offset = ($type >= 4) ? $mask_large_half : $mask_small_half;
        my $cY   = -$alv_in_y + $mask_offset;
        my $cZ   = 0;
        my $capY = $cY;
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
        my $serial   = $crs_id{$iy}[$ti] // 0;
        $detector{"identifiers"} = "sector manual $type xch manual $serial ych manual $iy zch manual $ti SiPM manual 6025";
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

sub make_65421_crs{
    my %detector = init_det();

    my $ix = $_[1];
    my $iy = $_[2];
    my $iz = $_[3];

    my @types = (6, 5, 4, 2, 1);
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

    my $Y_layer = ($iy - 5.5) * $layer_pitch;

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
        my $rot      = $is_even ? "0*deg 0*deg 0*deg" : "0*deg 180*deg 0*deg";
        my $cap_sign = $is_even ? 1 : -1;

        # Large crystals (types 4-8): cap mask 7x7 cm → center 3.5 cm from alveolus bottom
        # Small crystals (types 1-3): cap mask 6x6 cm → center 2.5 cm from alveolus bottom (~1 cm lower)
        my $mask_offset = ($type >= 4) ? $mask_large_half : $mask_small_half;
        my $cY   = -$alv_in_y + $mask_offset;
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
        my $serial   = $crs_id{$iy}[$ti] // 0;
        $detector{"identifiers"} = "sector manual $type xch manual $serial ych manual $iy zch manual $ti SiPM manual 6025";
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
    
    # Bottom edge of lowest layer (iy=0): centre = (0-5.5)*layer_pitch, bottom face = centre - alv_y
    my $crystal_stack_bot = (0 - 5.5) * $layer_pitch - $alv_y;
    my $scint_ht    = 1.0;                              # half-thickness of scint (2x2 cm)
    my $scint_bot_Y = $crystal_stack_bot - 0.5 - $scint_ht;   # 0.5cm gap below bottom layer
    my $scint_top_Y = $scint_bot_Y + 60;                # 1 m above bottom scint

    # 4 top scints equally spaced 3 cm apart in Z, centered on the array
    # Gaps between scint faces (left→right): 2.5, 3, 3 cm  (each scint is 2 cm wide in Z)
    # center-to-center: 4.5, 5, 5 cm  →  centered: -7.25, -2.75, +2.25, +7.25
    my @Ztop = (-7.25, -2.75, 2.25, 7.25);
    my $top_length = 35./2.;   # length of top scint in cm

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

    # iy 0-6:  7 rows of (6,5,4,3,1)
    for my $layer (0..6)  { make_65431_crs("1", 0, $layer, 0); }
    # iy 7:    1 row  of (6,5,4,2,1)
    make_65421_crs("1", 0, 7, 0);
    # iy 8-10: 3 rows of (8,7,7,2,2)
    for my $layer (8..10) { make_87722_crs("1", 0, $layer, 0); }

    # 4 bottom scints — same Z positions as top, at scint_bot_Y
    my @Zbot = (-7.25, -2.75, 2.25, 7.25);
    my $bot_length = 45./2.;   # length of bottom scint in cm

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
