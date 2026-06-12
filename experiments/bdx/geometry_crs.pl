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
    my $par1 = 200;
    my $par2 = 200.;
    my $par3 = 400.;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "KryptoniteLight";
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
    $detector{"mother"}      = "detector_volume";
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
    $detector{"identifiers"} = "sector manual 0 xch manual 0 ych manual 1 zch manual 0 SiPM manual 6025";
    #print_det(\%configuration, \%detector);
    
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
    #print_det(\%configuration, \%detector);
    
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
    #print_det(\%configuration, \%detector);
    
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
    
    my %detector = init_det();
    $detector{"name"}        = "plastic_up";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic_up";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 6*cm 10*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "1*cm 0.25*cm 10*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 0 scint manual 0 channel manual 0 module manual 0";
    print_det(\%configuration, \%detector);
    
    make_bgo_crs("3", 0, 0, 0);
    
    %detector = init_det();
    $detector{"name"}        = "plastic_down";
    $detector{"mother"}      = "detector_volume";
    $detector{"description"} = "plastic_down";
    $detector{"color"}       = "555555";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm -7*cm 10*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "1*cm 0.25*cm 10*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "scint";
    $detector{"hit_type"}    = "scint";
    $detector{"identifiers"} = "sector manual 0 scint manual 0 channel manual 1 module manual 0";
    print_det(\%configuration, \%detector);
}

sub make_all
{
    make_main;
    make_detector;
}


1;
