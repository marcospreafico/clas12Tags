use strict;
use warnings;

our %configuration;

my $ncol = 3; # Number of columns (X width)      #             ||
my $nrow = 3; # Number of rows (Y width)         #             \/
my $ndep = 2; # Number of (ncol x nrow) matrixes

# INNER VETO dimensions
# This is the BDX-PROTO veto
#   ______
#  |______|
#  | |  | |
#  |------|
#   ------
# Front and back panel are inserted inside the side

my $iv_t_x = 38.0/2.0;
my $iv_t_y = 1.0/2.0;
my $iv_t_z = 74.0/2.0;

my $iv_s_x = 1.0/2.0;
my $iv_s_y = 35.0/2.0;
my $iv_s_z = 74.0/2.0;

my $iv_f_x = 34.8/2.0;
my $iv_f_y = 34.5/2.0;
my $iv_f_z = 1.0/2.0;

# OUTER VETO dimensions
# Mounting of sides:
#   ______
#  | |____|
#  | |  | |
#  |----| |
#   ------
# Front and back panel are placed inside the side panels

my $ov_t_x = 40.0/2.0;
my $ov_t_y = 2.0/2.0;
my $ov_t_z = 80.0/2.0;

my $ov_s_x = 2.0/2.0;
my $ov_s_y = 40.0/2.0;
my $ov_s_z = 80.0/2.0;

my $ov_f_x = 38.0/2.0;
my $ov_f_y = 38.0/2.0;
my $ov_f_z = 2.0/2.0;


# Cal center in X, Y, Z
my $X0=0; my $Y0=0; my $Z0=0;

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

my $d = 1./2; #thickness of aluminum crs cap

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
    
    my $par1 = 200;
    my $par2 = 200.;
    my $par3 = 400.;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}

sub make_bgo_crs{
    
    my %detector = init_det();
    
    my $ix = $_[1];
    my $iy = $_[2];
    my $iz = $_[3];
    
    my $X = -2.0*$alv_x*5+$alv_x+$ix*$alv_x*2.0;
    my $Y = -2.0*$alv_y*4+$alv_y+$iy*$alv_y*2.0;
    my $Z = ($iz-1)*$alv_z*2.0;
    
    $X=0; $Y=0; $Z=0;
    
    $detector{"name"}        = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "crypt";
    $detector{"description"} = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    
    $detector{"dimensions"}  = "$alv_x*cm $alv_y*cm $alv_z*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
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
    
    $detector{"dimensions"}  = "$alv_in_x*cm $alv_in_y*cm $alv_in_z*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
    
    
    
    
    my $type = $_[0];
    if ($type eq "3") {
        $type = "4";
    }
    
    my $b = $bgo_types{$type}{b};
    my $B = $bgo_types{$type}{B};
    my $h = $bgo_types{$type}{h};
    my $H = $bgo_types{$type}{H};
    my $L = $bgo_types{$type}{L};
    my $theta = $bgo_types{$type}{theta};
    
    
    $Y = -$H/2.-$h/2.;
    $Z = 0;
    if($type eq "3" or $type eq "4"){
        $Y = -0.5*$bgo_types{"4"}{H}-0.5*$bgo_types{"3"}{h};
        $Z = 0.5;
    }
    
    %detector = init_det();
    $detector{"name"}        = "crs0_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs0";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "0*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
    $detector{"material"}    = "G4_BGO";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 0 xch manual 0 ych manual 0 zch manual 0 SiPM manual 6025";
    print_det(\%configuration, \%detector);
    
    
    $Y = $Y+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z+($L+$d)*cos($theta* pi / 180.0);
    
    %detector = init_det();
    $detector{"name"}        = "cap_crs0_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs0 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);
    
    %detector = init_det();
    $detector{"name"}        = "cap_hole_crs0_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cap_crs0_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs0 cap";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm 2*cm $d*cm 0*deg 360*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
    $Z = 0;
    $Y = $H/2.+$h/2.+0.2;
    
    if($type eq "3" or $type eq "4"){
        $Y = 0.5*$bgo_types{"4"}{H}+0.5*$bgo_types{"3"}{h}+0.1;
        $b = $bgo_types{"3"}{b};
        $B = $bgo_types{"3"}{B};
        $h = $bgo_types{"3"}{h};
        $H = $bgo_types{"3"}{H};
        $L = $bgo_types{"3"}{L};
        #$theta = $bgo_types{"3"}{theta};
        $Z = -0.5;
    }
    
    %detector = init_det();
    $detector{"name"}        = "crs1_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs1";
    $detector{"color"}       = "00bbbb";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "0*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $b*cm $H*cm $h*cm $L*cm";
    $detector{"material"}    = "G4_BGO";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 0 xch manual 0 ych manual 0 zch manual 0 SiPM manual 6025";
    print_det(\%configuration, \%detector);
    
    $Y = $Y-+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z-($L+$d)*cos($theta* pi / 180.0);
    
    %detector = init_det();
    $detector{"name"}        = "cap_crs1_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs0 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);
    
    %detector = init_det();
    $detector{"name"}        = "cap_hole_crs1_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cap_crs1_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs0 cap";
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

sub make_lead{
    my %detector = init_det();
    $detector{"name"}        = "lead";
    $detector{"mother"}      = "proto_mother";
    $detector{"description"} = "lead";
    $detector{"color"}       = "222222";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    my $pb_x = $iv_f_x;
    my $pb_y = $iv_f_y;
    my $pb_z = $iv_t_z-2*$iv_f_z;
    $detector{"dimensions"}  = "$pb_x*cm $pb_y*cm $pb_z*cm";
    $detector{"material"}    = "G4_Pb";
    #$detector{"material"}   = "KryptoniteLight";
    print_det(\%configuration, \%detector);
    
    %detector = init_det();
    $detector{"name"}        = "crypt";
    $detector{"mother"}      = "lead";
    $detector{"description"} = "crypt";
    $detector{"color"}       = "222222";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    my $crypt_x = $pb_x - 5.; #one lead block per side
    my $crypt_y = $pb_y - 2.5-0.5; #block below and foil above
    my $crypt_z = $pb_z - 1.; #just foil
    
    my $Y = 2;
    
    $detector{"pos"}         = "0*cm $Y*cm 0*cm";
    $detector{"dimensions"}  = "$crypt_x*cm $crypt_y*cm $crypt_z*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}   = "KryptoniteLight";
    print_det(\%configuration, \%detector);
}

#BEGIN inner veto
sub make_iveto{
        # Assuming fixed parameters for lead
        # Dimensions defined as follows
        #     _ _________
        #  y |_|_________|
        #     x     z
        my %detector = init_det();
        
        $detector{"mother"}      = "proto_mother";
        $detector{"material"}    = "ScintillatorB";
        
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        
        # UPSTREAM
        $detector{"name"}        = "iveto_upstream";
        $detector{"description"} = "inner veto upstream";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = 0;
        my $Y = 0;
        my $Z = -$iv_t_z+$iv_f_z;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_f_x*cm $iv_f_y*cm $iv_f_z*cm ";
        my $ch_id = 9;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 3";
        # WARNING!
        # Need to change the parameters of the veto to accomodate a new hitprocess that uses the real parameters of this prototype
        print_det(\%configuration, \%detector);
        
        # DOWNSTREAM
        $detector{"name"}        = "iveto_downstream";
        $detector{"description"} = "inner veto downstream";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = 0;
        $Y = 0;
        $Z = -$Z;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_f_x*cm $iv_f_y*cm $iv_f_z*cm ";
        $ch_id = 10;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 4";
        print_det(\%configuration, \%detector);
        
        # TOP
        $detector{"name"}        = "iveto_top";
        $detector{"description"} = "inner veto top";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = 0;
        $Y = $iv_s_y+$iv_t_y;
        $Z = 0;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_t_x*cm $iv_t_y*cm $iv_t_z*cm ";
        $ch_id = 1;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 1";
        print_det(\%configuration, \%detector);
        
        # BOTTOM
        $detector{"name"}        = "iveto_bottom";
        $detector{"description"} = "inner veto bottom";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = 0;
        $Y = -$Y;
        $Z = 0;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_t_x*cm $iv_t_y*cm $iv_t_z*cm ";
        $ch_id = 2;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 2";
        print_det(\%configuration, \%detector);
        
        # LEFT
        $detector{"name"}        = "iveto_left";
        $detector{"description"} = "inner veto left";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = -$iv_t_x+$iv_s_x;
        $Y = 0;
        $Z = 0;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_s_x*cm $iv_s_y*cm $iv_s_z*cm ";
        $ch_id = 5;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 6";
        print_det(\%configuration, \%detector);
        
        # RIGHT
        $detector{"name"}        = "iveto_right";
        $detector{"description"} = "inner veto right";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = -$X;
        $Y = 0;
        $Z = 0;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_s_x*cm $iv_s_y*cm $iv_s_z*cm ";
        $ch_id = 5;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 5";
        print_det(\%configuration, \%detector);
}
#END inner veto

#BEGIN outer veto
sub make_oveto{
    # Assuming fixed parameters for lead
    # Dimensions defined as follows
    #     _ _________
    #  y |_|_________|
    #     x     z
    
    my %detector = init_det();
    
    $detector{"mother"}      = "proto_mother";
    $detector{"material"}    = "ScintillatorB";
    
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    
    # UPSTREAM
    $detector{"name"}        = "oveto_upstream";
    $detector{"description"} = "inner veto upstream";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $X = 0;
    my $Y = 0;
    my $Z = -$ov_t_z+$ov_f_z;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ov_f_x*cm $ov_f_y*cm $ov_f_z*cm ";
    my $ch_id = 9;
    $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 3";
    # WARNING!
    # Need to change the parameters of the veto to accomodate a new hitprocess that uses the real parameters of this prototype
    print_det(\%configuration, \%detector);
    
    # DOWNSTREAM
    $detector{"name"}        = "oveto_downstream";
    $detector{"description"} = "inner veto downstream";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = 0;
    $Y = 0;
    $Z = -$Z;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ov_f_x*cm $ov_f_y*cm $ov_f_z*cm ";
    $ch_id = 10;
    $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 4";
    print_det(\%configuration, \%detector);
    
    # TOP
    $detector{"name"}        = "oveto_top";
    $detector{"description"} = "inner veto top";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = -$ov_s_x;
    $Y = $ov_s_y;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ov_t_x*cm $ov_t_y*cm $ov_t_z*cm ";
    $ch_id = 1;
    $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 1";
    print_det(\%configuration, \%detector);
    
    # BOTTOM
    $detector{"name"}        = "oveto_bottom";
    $detector{"description"} = "inner veto bottom";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = -$X;
    $Y = -$Y;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ov_t_x*cm $ov_t_y*cm $ov_t_z*cm ";
    $ch_id = 2;
    $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 2";
    print_det(\%configuration, \%detector);
    
    # LEFT
    $detector{"name"}        = "oveto_left";
    $detector{"description"} = "inner veto left";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $Y = $ov_t_y;
    $X = $ov_t_x;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ov_s_x*cm $ov_s_y*cm $ov_s_z*cm ";
    $ch_id = 5;
    $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 6";
    print_det(\%configuration, \%detector);
    
    # RIGHT
    $detector{"name"}        = "oveto_right";
    $detector{"description"} = "inner veto right";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = -$X;
    $Y = -$Y;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ov_s_x*cm $ov_s_y*cm $ov_s_z*cm ";
    $ch_id = 5;
    $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 5";
    print_det(\%configuration, \%detector);
}
#END outer veto

sub make_crs{
    make_crystal(0, 0, 0, 0, 0, 0);
}

sub make_proto{
    my $X = 0;
    my $Y = 0;
    my $Z = 0;
    
    my $proto_X = 42.0/2.0;
    my $proto_Y = 42.0/2.0;
    my $proto_Z = 80.0/2.0;
    
    my %detector = init_det();
    $detector{"name"}        = "proto_mother";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "mother of prototype";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$proto_X*cm $proto_Y*cm $proto_Z*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
    make_iveto;
    make_oveto;
    make_lead;
    make_bgo_crs("3", 0, 0, 0);;
}

sub make_all
{
    make_main;
    make_proto;
}


1;
