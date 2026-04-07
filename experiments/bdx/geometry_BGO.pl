#CAMBIARE
#lati a svastica!
# per i lati 9 fibre di lettura
# per il sopra e sotto 11 (OV) e 10 (IV)
# spazio tra esterno IV e interno OV 45 mm
# lastra di acciaio sopra per tenere pb sopra
# mettere cristalli di prad nella parte di pbwo

use strict;
use warnings;

use Math::Trig;

#TO DO: change veto identifiers to use sector to discriminate between IV and OV
our %configuration;

my $degrad = 0.01745329252;
my $cic    = 2.54;

#BEGIN general geometry parameters


# crypt volume (crytal space)
my $crypt_x = 80./2;
my $crypt_y = 80./2;
my $crypt_z = 130./2;

my $bgo_z = 85./2; # space allocated for BGO part of calorimeter
my $pbwo_z = $crypt_z-$bgo_z; # remaining space for PbWO


my $air_gap = 1.; #air grap between lead and veto

# lead shielding
my $pb_t = 5./2; #lead thikness
my $pb_x = $crypt_x+2.*$pb_t+$air_gap;
my $pb_y = $crypt_y+2.*$pb_t+$air_gap;
my $pb_z = $crypt_z+2.*$pb_t+$air_gap;


my $veto_t = 1.5/2; #veto thickness
# veto

my $overhang = 2.; #extra length of veto for overhanging portions

# General shape of the veto from the top.
# Front face is between sides and top
# Side panel closes caps and is between two layers of top
# Overhanging parts to increase hermeticity
#    _______________________________________________
#   |   _________________________________________   |
#   |  |_________________________________________|  |
#   |    |  |                               |  |    |
#   |    |  |                               |  |    |
#   |    |  |                               |  |    |
#   |    |  |                               |  |    |
#   |    |  |                               |  |    |
#   |   _|__|_______________________________|__|_   |
#   |  |_________________________________________|  |
#   |_______________________________________________|
#

# BEGIN IV

# size (mm x mm)
# 930x965 IV front x
# 930x1465 IV side x
# 1500x1000 IV top x

#front panel - between top and side paddle
my $iv_f_x = 96.5/2;
my $iv_f_y = 93.0/2;
my $iv_f_z = $veto_t;

#side panel - closes the top and is between up and down
my $iv_s_x = $veto_t;
my $iv_s_y =  93.0/2;
my $iv_s_z = 146.5/2;

#top ppanel - closes front and side panels
my $iv_t_x = 100.0/2;
my $iv_t_y = $veto_t;
my $iv_t_z = 150.0/2;

#general dimensions of iv
my $iv_x = ($iv_f_x > $iv_s_x) ? $iv_f_x : $iv_s_x;
$iv_x = ($iv_x > $iv_t_x) ? $iv_x : $iv_t_x;
my $iv_y = ($iv_f_y > $iv_s_y) ? $iv_f_y : $iv_s_y;
$iv_y = ($iv_y > $iv_t_y) ? $iv_y : $iv_t_y;
my $iv_z = ($iv_f_z > $iv_s_z) ? $iv_f_z : $iv_s_z;
$iv_z = ($iv_z > $iv_t_z) ? $iv_z : $iv_t_z;

$iv_x = $iv_x+2*$veto_t+$air_gap;
$iv_y = $iv_y+2*$veto_t+$air_gap;
$iv_z = $iv_z+2*$veto_t+$air_gap;
#END IV

# BEGIN OV

# size (mm x mm)
# 970x1065 OV front x
# 1575x970 OV side x
# 1610x1100 OV top x

#front panel - between top and side paddle
my $ov_f_x = 106.5/2;
my $ov_f_y =  97.0/2;
my $ov_f_z = $veto_t;

#side panel - closes the top and is between up and down
my $ov_s_x = $veto_t;
my $ov_s_y =  97.0/2;
my $ov_s_z = 157.5/2;

#top ppanel - closes front and side panels
my $ov_t_x = 110.0/2;
my $ov_t_y = $veto_t;
my $ov_t_z = 161.0/2;

#general dimensions of ov
my $ov_x = ($ov_f_x > $ov_s_x) ? $ov_f_x : $ov_s_x;
$ov_x = ($ov_x > $ov_t_x) ? $ov_x : $ov_t_x;
my $ov_y = ($ov_f_y > $ov_s_y) ? $ov_f_y : $ov_s_y;
$ov_y = ($ov_y > $ov_t_y) ? $ov_y : $ov_t_y;
my $ov_z = ($ov_f_z > $ov_s_z) ? $ov_f_z : $ov_s_z;
$ov_z = ($ov_z > $ov_t_z) ? $ov_z : $ov_t_z;

$ov_x = $ov_x+2*$veto_t+$air_gap;
$ov_y = $ov_y+2*$veto_t+$air_gap;
$ov_z = $ov_z+2*$veto_t+$air_gap;
#END OV

#vessel
my $vessel_t = 1.; #thickness of vessel

my $vessel_x = $ov_x+$air_gap+$vessel_t;
my $vessel_y = $ov_y+$air_gap+$vessel_t;
my $vessel_z = $ov_z+$air_gap+$vessel_t;

#overall detector volume

my $detector_x = $vessel_x+$air_gap;
my $detector_y = $vessel_y+$air_gap;
my $detector_z = $vessel_z+$air_gap;



#END general geometry parameters


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
my $alv_z = $bgo_z/3;

my $alv_t = 0.2;

my $alv_in_x = $alv_x-$alv_t;
my $alv_in_y = $alv_y-$alv_t;
my $alv_in_z = $alv_z;
#END BGO parameters

#BEGIN PbWO parameters
my %pbwo_types = (
                 "EC" => {
                 b => 4.8/2.,
                 B => 5.2/2.,
                 h => 4.8/2.,
                 H => 5.2/2.,
                 L => 20./2.,
                 theta => atan((5.2-4.8)/(2.*20))/2*180./pi
                 },
"PRad" => {
    b => 4.10/2.,
    B => 4.10/2.,
    h => 4.10/2.,
    H => 4.10/2.,
    L => 18./2.,
    theta => 0.
    }
);


my $alv_pbwo_x = 6./2;
my $alv_pbwo_y = 11./2;
my $alv_pbwo_z = $pbwo_z/2.;

my $x_pbwo = int($crypt_x/$alv_pbwo_x);
my $y_pbwo = int($crypt_y/$alv_pbwo_y);

my $alv_pbwo_in_x = $alv_pbwo_x-$alv_t;
my $alv_pbwo_in_y = $alv_pbwo_y-$alv_t;
my $alv_pbwo_in_z = $alv_pbwo_z;
#END PbWO parametrers

my $d = 1./2; #thickness of aluminum crs cap

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
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
}



sub make_lead{
    my %detector = init_det();
    $detector{"name"}        = "lead_volume";
    $detector{"mother"}      = "iv_volume";
    $detector{"description"} = "lead";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    $detector{"dimensions"}  = "$pb_x*cm $pb_y*cm $pb_z*cm";
    $detector{"material"}    = "G4_Pb";
    print_det(\%configuration, \%detector);
    
    %detector = init_det();
    $detector{"name"}        = "crypt";
    $detector{"mother"}      = "lead_volume";
    $detector{"description"} = "crypt";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    $detector{"dimensions"}  = "$crypt_x*cm $crypt_y*cm $crypt_z*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
}



#BEGIN inner veto
sub make_iv{
    my %detector = init_det();
    $detector{"mother"}      = "ov_volume";
    $detector{"material"}    = "G4_AIR";
    $detector{"name"}        = "iv_volume";
    $detector{"description"} = "inner veto volume";
    $detector{"color"}       = "0000FF";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$iv_x*cm $iv_y*cm $iv_z*cm";
    print_det(\%configuration, \%detector);
    
    %detector = init_det();
    
    
    %detector = init_det();
    
    # general parameters
    $detector{"mother"}      = "iv_volume";
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
    my $X = 2*$veto_t;
    my $Y = 0;
    my $Z = -$pb_z-$iv_f_z;
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
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = -$X;
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
    $Y = $pb_y+$iv_t_y;
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
    $detector{"style"}       = 1;
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
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $Y = 0;
    $X = $pb_x+$iv_s_x;
    $Z = -2*$veto_t;
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
    $Z = $Z;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$iv_s_x*cm $iv_s_y*cm $iv_s_z*cm ";
    $ch_id = 5;
    $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 5";
    print_det(\%configuration, \%detector);
}
#END inner veto

#BEGIN outer veto
sub make_ov{
    my %detector = init_det();
    $detector{"mother"}      = "vessel_air";
    $detector{"material"}    = "G4_AIR";
    $detector{"name"}        = "ov_volume";
    $detector{"description"} = "outer veto volume";
    $detector{"color"}       = "0000FF";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ov_x*cm $ov_y*cm $ov_z*cm";
    print_det(\%configuration, \%detector);
    
    %detector = init_det();
    
    
    %detector = init_det();
    
    # general parameters
    $detector{"mother"}      = "ov_volume";
    $detector{"material"}    = "ScintillatorB";
    
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    
    # UPSTREAM
    $detector{"name"}        = "oveto_upstream";
    $detector{"description"} = "outer veto upstream";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $X = 0;
    my $Y = 0;
    my $Z = -$iv_z-$ov_f_z;
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
    $detector{"description"} = "outer veto downstream";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 1;
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
    $detector{"description"} = "outer veto top";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = 0;
    $Y = $iv_y+$iv_t_y;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ov_t_x*cm $ov_t_y*cm $ov_t_z*cm ";
    $ch_id = 1;
    $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 1";
    print_det(\%configuration, \%detector);
    
    # BOTTOM
    $detector{"name"}        = "oveto_bottom";
    $detector{"description"} = "outer veto bottom";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = 0;
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
    $detector{"description"} = "outer veto left";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $Y = 0;
    $X = $iv_x+$iv_s_x;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ov_s_x*cm $ov_s_y*cm $ov_s_z*cm ";
    $ch_id = 5;
    $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 6";
    print_det(\%configuration, \%detector);
    
    # RIGHT
    $detector{"name"}        = "oveto_right";
    $detector{"description"} = "outer veto right";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = -$X;
    $Y = 0;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ov_s_x*cm $ov_s_y*cm $ov_s_z*cm ";
    $ch_id = 5;
    $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 5";
    print_det(\%configuration, \%detector);
}
#END outer veto

sub make_vessel{
    my %detector = init_det();
    $detector{"mother"}      = "detector_volume";
    $detector{"material"}    = "G4_Al";
    $detector{"name"}        = "vessel";
    $detector{"description"} = "vessel";
    $detector{"color"}       = "333333";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$vessel_x*cm $vessel_y*cm $vessel_z*cm ";
    print_det(\%configuration, \%detector);
    
    %detector = init_det();
    $detector{"mother"}      = "vessel";
    $detector{"material"}    = "G4_AIR";
    $detector{"name"}        = "vessel_air";
    $detector{"description"} = "vessel air";
    $detector{"color"}       = "333333";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    my $par1 = $vessel_x-$vessel_t;
    my $par2 = $vessel_y-$vessel_t;
    my $par3 = $vessel_z-$vessel_t;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm ";
    print_det(\%configuration, \%detector);
    
    %detector = init_det();
}

sub make_bgo_crs{
    
    my %detector = init_det();
    
    my $ix = $_[1];
    my $iy = $_[2];
    my $iz = $_[3];
    
    my $X = -2.0*$alv_x*5+$alv_x+$ix*$alv_x*2.0;
    my $Y = -2.0*$alv_y*4+$alv_y+$iy*$alv_y*2.0;
    my $Z = ($iz-1)*$alv_z*2.0;
    
    $detector{"name"}        = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "bgo_volume";
    $detector{"description"} = "alveolus_$ix"."_"."$iy"."_"."$iz";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    
    $detector{"dimensions"}  = "$alv_x*cm $alv_y*cm $alv_z*cm";
    $detector{"material"}    = "G4_Al";
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

sub make_bgo_cal{
    
    my %detector = init_det();
    $detector{"name"}        = "bgo_volume";
    $detector{"mother"}      = "crypt";
    $detector{"description"} = "bgo_volume";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    my $X = 0;
    my $Y = 0;
    my $Z = 0; #-$crypt_z+$bgo_z;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    $detector{"dimensions"}  = "$crypt_x*cm $crypt_y*cm $bgo_z*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
    
    my @counter = (0, 0, 0, 0, 0, 0, 0, 0);
    my $id = 0;
    for(my $iz = 0; $iz < 3; $iz++){
        for(my $iy = 0; $iy < 8; $iy++){
            for(my $ix = 0; $ix < 10; $ix++){
                if($ix == 0 or $iy == 0 or $ix == 9 or $iy == 7){
                    $id = 1;
                    if($counter[0] == 32) {
                        $id = 2;
                    }
                    if($counter[0] == 32 and $counter[1] ==32){
                        for(my $ii = 3; $ii<9 ; $ii++){
                            if($counter[$ii-1] != 32){
                                $id = $ii;
                                last;
                            }
                        }
                    }
                    make_bgo_crs("$id", $ix, $iy, $iz);
                    $counter[$id-1] = $counter[$id-1]+1;
                } else {
                    for(my $ii = 3; $ii<9 ; $ii++){
                        if($counter[$ii-1] != 32){
                            $id = $ii;
                            last;
                        }
                    }
                    make_bgo_crs("$id", $ix, $iy, $iz);
                    $counter[$id-1] = $counter[$id-1]+1;
                }
            }
        }
    }
    
    for(my $ii = 0; $ii < 8; $ii++){
        print "$counter[$ii] , "
    }
}

sub make_pbwo_cal{
    
    my %detector = init_det();
    $detector{"name"}        = "prad_cal";
    $detector{"mother"}      = "crypt";
    $detector{"description"} = "pbwo_volume";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    my $X = 0;
    my $Y = 0;
    my $Z = -$crypt_z+$pbwo_z/2;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    $detector{"dimensions"}  = "$crypt_x*cm $crypt_y*cm $pbwo_z*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
    
    my $ncrs = 0;
    
    
        for(my $iy = 0; $iy < 17; $iy++){
            for(my $ix = 0; $ix < 17; $ix++){
                make_pbwo_crs("PRad", $ix, $iy, 0);
            }
        }
    
    %detector = init_det();
    $detector{"name"}        = "pbwo_volume";
    $detector{"mother"}      = "crypt";
    $detector{"description"} = "pbwo_volume";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    my $X = 0;
    my $Y = 0;
    my $Z = $crypt_z-$pbwo_z/2;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    $detector{"dimensions"}  = "$crypt_x*cm $crypt_y*cm $pbwo_z*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
        for(my $iy = 0; $iy < $y_pbwo; $iy++){
            for(my $ix = 0; $ix < $x_pbwo; $ix++){
                make_pbwo_crs("EC", $ix, $iy, 1);
                $ncrs = $ncrs+4;
            }
        }
    
    
    
    print "Used $ncrs pbwo crystals \n";
    
}

sub make_pbwo_crs{
    
    my %detector = init_det();
    
    my $type = $_[0];
    my $ix = $_[1];
    my $iy = $_[2];
    my $iz = $_[3];
    
    my $b = $pbwo_types{$type}{b};
    my $B = $pbwo_types{$type}{B};
    my $h = $pbwo_types{$type}{h};
    my $H = $pbwo_types{$type}{H};
    my $L = $pbwo_types{$type}{L};
    my $theta = $pbwo_types{$type}{theta};
    
    my $X = -$alv_pbwo_x*int($x_pbwo)+$alv_pbwo_x+$ix*2.0*$alv_pbwo_x;
    my $Y = -$alv_pbwo_y*int($y_pbwo)+$alv_pbwo_y+$iy*2.0*$alv_pbwo_y;
    my $Z = 0; #1.*$pbwo_z/2;
    
    if($type eq "PRad"){
     $X = -$B*17+$B+$ix*2.0*$B;
    $Y = -$H*17+$H+$iy*2.0*$H;
        $Z = 0; #-1.*$pbwo_z/2;
    }
    
    print("$pbwo_z $Z \n");
    
    $detector{"name"}        = "alveolus_pbwo_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "pbwo_volume";
    if($type eq "PRad"){
        $detector{"mother"}     = "prad_cal";
    }
    $detector{"description"} = "alveolus_pbwo_$ix"."_"."$iy"."_"."$iz";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    
    $detector{"dimensions"}  = "$alv_pbwo_x*cm $alv_pbwo_y*cm $alv_pbwo_z*cm";
    if($type eq "PRad"){
        my $L1 = $L+1;
        $detector{"dimensions"}  = "$B*cm $H*cm $L1*cm";
    }
    $detector{"material"}    = "G4_Al";
    #$detector{"material"}    = "KryptoniteLight";
        print_det(\%configuration, \%detector);
    
    %detector = init_det();
    $detector{"name"}        = "alveolus_pbwo_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_pbwo_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "alveolus air";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    
    $detector{"dimensions"}  = "$alv_pbwo_in_x*cm $alv_pbwo_in_y*cm $alv_pbwo_in_z*cm";
    if($type eq "PRad"){
        my $L1 = $L+1;
        $detector{"dimensions"}  = "$B*cm $H*cm $L1*cm";
    }
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
        print_det(\%configuration, \%detector);
    
    
    $Y = -$H/2.-$h/2.;
    $Z = 0;
    if($type eq "PRad"){
     $Y=0;
    }
    
    %detector = init_det();
    $detector{"name"}        = "crs_pbwo0_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_pbwo_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs0";
    if($type eq "PRad"){
        $detector{"color"} = "ffd700";
    }
    else{
        $detector{"color"}       = "d7263d";
    }
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "0*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$b*cm $B*cm $h*cm $H*cm $L*cm";
    $detector{"material"}    = "G4_PbWO4";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 0 xch manual 0 ych manual 0 zch manual 0 SiPM manual 6025";
    print_det(\%configuration, \%detector);
    
    
    $Y = $Y+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z+($L+$d)*cos($theta* pi / 180.0);
    
    %detector = init_det();
    $detector{"name"}        = "cap_pbwo_crs0_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_pbwo_air_$ix"."_"."$iy"."_"."$iz";
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
    
    
    $Z = 0;
    $Y = $H/2.+$h/2.+0.2;
    
    %detector = init_det();
    $detector{"name"}        = "crs1_pbwo_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_pbwo_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs1";
    if($type eq "PRad"){
        $detector{"color"} = "ffd700";
    }
    else{
        $detector{"color"}       = "d7263d";
    }
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "0*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $b*cm $H*cm $h*cm $L*cm";
    $detector{"material"}    = "G4_PbWO4";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 0 xch manual 0 ych manual 0 zch manual 0 SiPM manual 6025";
    if($type ne "PRad"){
        print_det(\%configuration, \%detector);
    }
    
    $Y = $Y-+($L+$d)*sin($theta * pi /180.0);
    $Z = $Z-($L+$d)*cos($theta* pi / 180.0);
    
    %detector = init_det();
    $detector{"name"}        = "cap_pbwo_crs1_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "alveolus_pbwo_air_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "crs0 cap";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$theta*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$B*cm $H*cm $d*cm";
    $detector{"material"}    = "G4_Al";
    if($type ne "PRad"){
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
    
    $detector{"dimensions"}  = "$detector_x*cm $detector_y*cm $detector_z*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
    make_vessel;
    make_ov;
    make_iv;
    make_lead;
    make_bgo_cal;
    make_pbwo_cal;
}

sub make_all
{
    make_main;
    make_detector;
}


1;
