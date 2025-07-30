use strict;
use warnings;

our %configuration;

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
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    my $par1 = 200;
    my $par2 = 200.;
    my $par3 = 400.;
    $detector{"dimensions"}  = "500*cm 1000*cm 1500*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}

sub make_dirt{
    my %detector = init_det();
    $detector{"mother"}      = "main_volume";
    $detector{"name"}        = "dirt_box";
    $detector{"description"} = "Dirt";
    $detector{"color"}       = "614022";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "400*cm 750*cm 1400*cm";
    $detector{"material"}    = "BDX_Dirt";
    print_det(\%configuration, \%detector);
}


sub make_vault{
    my %detector = init_det();
    $detector{"mother"}      = "dirt_box";
    $detector{"name"}        = "vault_wall";
    $detector{"description"} = "Vault Concrete";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 316*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "151.92*cm 434*cm 558.5*cm";
    $detector{"material"}    = "BDX_Dirt";
    print_det(\%configuration, \%detector);
    
    $detector{"mother"}      = "vault_wall";
    $detector{"name"}        = "vault_air";
    $detector{"description"} = "Vault air";
    $detector{"color"}       = "614022";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 15*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "120*cm 419*cm 528.5*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
    
    $detector{"mother"}      = "vault_air";
    $detector{"name"}        = "vault_shielding_volume";
    $detector{"description"} = "Shielding volume";
    $detector{"color"}       = "614022";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm -343*cm -114.1*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "120*cm 76*cm 414.4*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}

sub make_shielding{
    my %detector = init_det();
    $detector{"mother"}      = "vault_shielding_volume";
    $detector{"name"}        = "steel_front";
    $detector{"description"} = "Steel front";
    $detector{"color"}       = "777777";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm -156.9*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "76*cm 76*cm 257.5*cm";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    
    $detector{"mother"}      = "vault_shielding_volume";
    $detector{"name"}        = "steel_det";
    $detector{"description"} = "Steel around detector";
    $detector{"color"}       = "777777";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 257.5*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "76*cm 76*cm 156.9*cm";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    
    $detector{"mother"}      = "steel_det";
    $detector{"name"}        = "det_air";
    $detector{"description"} = "Detector space";
    $detector{"color"}       = "dddddd";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "52*cm 52*cm 156.9*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
    $detector{"mother"}      = "vault_shielding_volume";
    $detector{"name"}        = "steel_l";
    $detector{"description"} = "Steel left";
    $detector{"color"}       = "777777";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "ITrd";
    $detector{"pos"}         = "87.1*cm 0*cm -233.4*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "21.5*cm 0.5*cm 76*cm 76*cm  181*cm -3.32*deg 0*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    
    $detector{"mother"}      = "vault_shielding_volume";
    $detector{"name"}        = "steel_r";
    $detector{"description"} = "Steel right";
    $detector{"color"}       = "777777";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "ITrd";
    $detector{"pos"}         = "-87.1*cm 0*cm -233.4*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "21.5*cm 0.5*cm 76*cm 76*cm  181*cm 3.32*deg 0*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    
    $detector{"mother"}      = "vault_air";
    $detector{"name"}        = "steel_t";
    $detector{"description"} = "Steel top";
    $detector{"color"}       = "777777";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "ITrd";
    $detector{"pos"}         = "0*cm -255.9*cm -347.5*cm";
    $detector{"rotation"}    = "0*deg 0*deg 90*deg";
    $detector{"dimensions"}  = "21.5*cm 0.5*cm 76*cm 76*cm  181*cm 3.32*deg 0*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    
    $detector{"mother"}      = "steel_front";
    $detector{"name"}        = "lead_1";
    $detector{"description"} = "Lead block 1";
    $detector{"color"}       = "333333";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm -227*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "37.5*cm 37.5*cm 28.5*cm";
    $detector{"material"}    = "G4_Pb";
    print_det(\%configuration, \%detector);
    
    $detector{"mother"}      = "steel_front";
    $detector{"name"}        = "lead_2";
    $detector{"description"} = "Lead block 2";
    $detector{"color"}       = "333333";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm -170*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "25*cm 25*cm 28.5*cm";
    $detector{"material"}    = "G4_Pb";
    print_det(\%configuration, \%detector);
    
    $detector{"mother"}      = "steel_front";
    $detector{"name"}        = "lead_3";
    $detector{"description"} = "Lead block 3";
    $detector{"color"}       = "333333";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm -113*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "12.5*cm 12.5*cm 28.5*cm";
    $detector{"material"}    = "G4_Pb";
    print_det(\%configuration, \%detector);
    
}




sub make_all
{
    make_main;
    make_dirt;
    make_vault;
    make_shielding;
}


1;
