use strict;
use warnings;

our %configuration;

# Full detector and or Hodo in the center of the hall
my $flag_detcentered = 0;
my $flag_inner_lead = 1;
my $flag_outer_lead = 0;
# to control Prototype into the cave
# $flag_JlabCT=0; -> LNS/CT config
# $flag_JlabCT=1; -> proto at JLab config
my $flag_JlabCT = 0 ;
if ($configuration{"variation"} eq "CT") {$flag_JlabCT = 0 ;}

# Nov16: to control test in front of Hall-A beam dump
# to measure muon flux
# $flag_mutest = 0 -> standard Proposal
# $flag_mutest = 1 -> Hall-A mu tests
### $flag_minibdx=0 standard BDX-Hodo
### $flag_minibdx!=0 BDX-Mini
my $flag_mutest =0;
my $flag_minibdx =0;

### JLab_paddle_setup=0 use bottom inner veto lid
### JLab_paddle_setup=0 use external paddle in place
my $JLab_paddle_setup = 1;
$JLab_paddle_setup =$JLab_paddle_setup *$flag_detcentered; # checking that the BDX-MINI is NOT in the pipe
###########################################################################################
###########################################################################################
# Define the relevant parameters of BDX geometry
#
# the BDX geometry will be defined starting from these parameters 
#
# all dimensions are in mm
#


my $degrad = 0.01745329252;
my $cic    = 2.54;


my $BDXmini_externalBox_X = 100;
my $BDXmini_externalBox_Y = 100;
my $BDXmini_externalBox_Z = 100;


###########################################################################################
# Build Crystal Volume and Assemble calorimeter
###########################################################################################



#################################################################################################
#
# Begin: Hall A Beam Dump
#
#################################################################################################
#
# Note: These numbers assume that the origin is at the upstream side of the BEAM DUMP.
#
# Generate particles at (0,0,-200)


sub make_bdx_main_volume
    {
     my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
        {

        $detector{"name"}        = "bdx_real_main_volume";
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
            my $wallthk=1; # now it's 15cm or 470cm
            my $par1 = 600.+$wallthk;
            my $par2 = 400.+$wallthk;
            my $par3 = 600.+$wallthk;

        if (($flag_detcentered eq "1"))
        {
            $par1 = 300.;
            $par2 = 300.;
            $par3 = 300.;
        }
            
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
            $detector{"material"}    = "G4_CONCRETE";
            # $detector{"material"}    = "G4_AIR";
        print_det(\%configuration, \%detector);
 

        my $p1= $par1-$wallthk  ;
        my $p2 =$par2-$wallthk ;
        my $p3 =$par3-$wallthk ;
        $detector{"name"}        = "main_volume";
        $detector{"mother"}      = "bdx_real_main_volume";
        $detector{"description"} = "concrete walls";
        $detector{"color"}       = "f00000";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$p1*cm $p2*cm $p3*cm";
        $detector{"material"}    = "G4_Galactic";
        print_det(\%configuration, \%detector);
    }
    else
    {
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
    
    my $par1 = 1000.;
    my $par2 = 2000.;
    my $par3 = 4000.;
    if (($configuration{"variation"} eq "BDXmini")){
	$par1=$BDXmini_externalBox_X;
	$par2=$BDXmini_externalBox_Y;
	$par3=$BDXmini_externalBox_Z;
    }

    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_Galactic";
    print_det(\%configuration, \%detector);
    }
}



my $Dirt_xmin = -800 ;
my $Dirt_xmax = +800. ;
my $Dirt_ymin = -762. ;
my $Dirt_ymax = +762. ;  # This is "xgrade" the depth of the beamline underground.
my $Dirt_zmin =-1000. ;
my $Dirt_zmax = 3200. ;

sub make_mutest_bdxMiniDirt{
    my %detector = init_det();

    my $par1 = 100;
    my $par2 = 100;
    my $par3 = 100;
    
    $detector{"name"}        = "dirt";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "Dirt";
    $detector{"color"}       = "f00000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"material"}    = "BDX_Dirt"; # if not defined use G4_SiO2
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    print_det(\%configuration, \%detector);
}


# To keep the beamdump at (0,0,0) we have symmetric x,y. For z we use two joined volumes.
sub make_dirt_u
{
    my %detector = init_det();
    
    my $X = ($Dirt_xmax+$Dirt_xmin) ;
    my $Y = ($Dirt_ymax+$Dirt_ymin) ;
    my $Z = 0;
    my $par1 = ($Dirt_xmax-$Dirt_xmin)/2.;
    my $par2 = ($Dirt_ymax-$Dirt_ymin)/2.;
    my $par3 = -($Dirt_zmin) ;

    $detector{"name"}        = "dirt_u";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "Upstream side of Mother volume of Earth/dirt, below grade level";
    $detector{"color"}       = "f00000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "Component";
    print_det(\%configuration, \%detector);
}
sub make_dirt_d
{
    my %detector = init_det();
    
    my $X = ($Dirt_xmax+$Dirt_xmin) ;
    my $Y = ($Dirt_ymax+$Dirt_ymin) ;
    my $Z = -$Dirt_zmin + ($Dirt_zmax+$Dirt_zmin)/2.;
    my $par1 = ($Dirt_xmax-$Dirt_xmin)/2.;
    my $par2 = ($Dirt_ymax-$Dirt_ymin)/2.;
    my $par3 = ($Dirt_zmax+$Dirt_zmin)/2. ;
    
    $detector{"name"}        = "dirt_d";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "Downpstream side of Mother volume of Earth/dirt, below grade level";
    $detector{"color"}       = "d00000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "Component";
    print_det(\%configuration, \%detector);
}
sub make_dirt
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = 0.;
    
    $detector{"name"}        = "dirt";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "Mother volume of Earth/dirt, below grade level";
    $detector{"color"}       = "D0A080";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Operation: dirt_u+dirt_d";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"material"}    = "BDX_Dirt"; # if not defined use G4_SiO2
    # $detector{"material"}    = "G4_AIR"; # if not defined use G4_SiO2
    print_det(\%configuration, \%detector);
}

sub make_dirt_top
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y =$Dirt_ymax+300./2.;
    my $Z = 400. ;
    my $par1 = 800./2;
    my $par2 = 300./2.;
    my $par3 = -$Dirt_zmin+800/2 ;
    
    $detector{"name"}        = "dirt_top";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "Hill protecting dump area";
    $detector{"color"}       = "D0A080";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Dirt";
    # $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}



my $Bunker_L_inner = 914. ;
my $Bunker_DZ_end = 548. ;
my $Bunker_Z_upstream = -600. ;
my $Bunker_zmin = $Bunker_Z_upstream;
my $Bunker_zmax = $Bunker_Z_upstream + $Bunker_L_inner + $Bunker_DZ_end ;
my $Bunker_dx = 564. ;
my $Bunker_dy = 564. ;

my $Bunker_cutout_l = 914 ;
my $Bunker_cutout_r = 213./2. ;
my $Bunker_cutout_shim = 1. ;



my $Bunker_end_dc = 91. ;
my $Bunker_end_dx = $Bunker_cutout_r + $Bunker_end_dc;
my $Bunker_end_dz = ($Bunker_dx - ($Bunker_cutout_r + $Bunker_end_dc))/2. ;

my $Bunker_dz = ($Bunker_L_inner + $Bunker_DZ_end)/2. - $Bunker_end_dz*2.+100 ;

my $Bunker_main_rel_z = $Bunker_Z_upstream + ($Bunker_zmax-$Bunker_zmin)/2.-100;

sub make_bunker_main
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = $Bunker_main_rel_z;
    my $par1 = $Bunker_dx;
    my $par2 = $Bunker_dy;
    my $par3 = $Bunker_dz;
    
    $detector{"name"}        = "Bunker_main";
    $detector{"mother"}      = "dirt";
    $detector{"description"} = "Main block volume of concrete bunker";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
}
sub make_bunker_tunnel
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = -90.;
    my $par1 = $Bunker_cutout_r;
    my $par2 = $Bunker_cutout_r;
    my $par3 = $Bunker_dz-90.;
    
    $detector{"name"}        = "Bunker_tunnel";
    $detector{"mother"}      = "Bunker_main";
    $detector{"description"} = "Cutout of bunker for tunnel";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}
sub make_bunker
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = $Bunker_main_rel_z;
    
    $detector{"name"}        = "Bunker";
    $detector{"mother"}      = "dirt";
    $detector{"description"} = "Cutout of bunker for tunnell";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Operation: Bunker_main-Bunker_tunnel";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"material"}    = "BDX_Concrete";
    #    print_det(\%configuration, \%detector);
}
sub make_bunker_end
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = $Bunker_main_rel_z + $Bunker_dz + $Bunker_end_dz;
    
    $detector{"name"}        = "Bunker_end";
    $detector{"mother"}      = "dirt";
    $detector{"description"} = "End cone of bunker";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$Bunker_dx*cm $Bunker_end_dx*cm $Bunker_dy*cm $Bunker_end_dx*cm $Bunker_end_dz*cm";
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
    
    # add flux detectors
    my $bunker_nflux=10;
    my $bunker_flux_dz=40.;
    my $bunker_flux_lx=130.;
    my $bunker_flux_ly=130.;
    my $bunker_flux_lz=0.01;
    for(my $iz=0; $iz<$bunker_nflux; $iz++)
    {
        my %detector = init_det();
        $detector{"name"}        = "bunker_flux_$iz";
        $detector{"mother"}      = "Bunker_end";
        $detector{"description"} = "bunker flux detector $iz ";
        $detector{"color"}       = "cc00ff";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = 0.;
        $Y = 0.;
        $Z=($iz-4.5)*$bunker_flux_dz+$bunker_flux_lz+1.;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$bunker_flux_lx*cm $bunker_flux_ly*cm $bunker_flux_lz*cm";
        $detector{"material"}    = "BDX_Concrete";
        $detector{"sensitivity"} = "flux";
        $detector{"hit_type"}    = "flux";
        my $nflux=101+$iz;
        $detector{"identifiers"} = "id manual $nflux";
        # print_det(\%configuration, \%detector);
    }

}

    my $beamdump_zdelta = 374;
    my $beamdump_zmin   = 0.;
    my $beamdump_zmax   = $beamdump_zmin + $beamdump_zdelta;
    my $beamdump_radius = 54.5;

my $beamdump_z = -($Bunker_dz+$beamdump_zdelta/2.)+2*$Bunker_dz-$beamdump_zdelta+145-91;

sub make_hallaBD_org
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = $beamdump_z;
    my $par1 = 0.;
    my $par2 = $beamdump_radius;
    my $par3 = $beamdump_zdelta/2.;
    my $par4 = 0.;
    my $par5 = 360.;
   
    $detector{"name"}        = "hallaBD";
    $detector{"mother"}      = "Bunker_tunnel";
    $detector{"description"} = "BD vessel";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_WATER";
    print_det(\%configuration, \%detector);
     $X = 0. ;
     $Y = 0. ;
     $Z = 83.5;
    
     $par1 = 0.;
     $par2 = 25.1;
     $par3 = 106/2.;
     $par4 = 0.;
     $par5 = 360.;
    $detector{"name"}        = "hallaBD-Al-solid";
    $detector{"mother"}      = "hallaBD";
    $detector{"description"} = "Al solid dump";
    $detector{"color"}       = "0000FF";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);
    
    my $nplates= 80;
    my $z1=-178.6;
    my $z2=$z1+198.;
    my $zfil=($z2-$z1)*80/100;
    my $zinc=0.01;
    
    my $zsum=0;
    my $zi=1;
    
    for(my $iz=0; $iz<$nplates; $iz++)
    {
        $zsum=$zsum+$zi;
        $zi= $zi*($zinc+1);
    }
        $zi=$zfil/$zsum;
        my $zw=($z2-$z1-$zfil)/$zsum;
        $Z=$z1;
        for(my $iz=0; $iz<$nplates; $iz++)
    {
        $Z=$Z+$zi/2+$zw/2;
        $par3=$zi/2;
        my %detector = init_det();
        $detector{"name"}        = "hallaBD-Al-plates_$iz";
        $detector{"mother"}      = "hallaBD";
        $detector{"description"} = "Al plates $iz ";
        $detector{"color"}       = "0000ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
        $Z=$Z +$zi/2+$zw/2;
        $zi= $zi*($zinc+1);
        $zw= $zw*($zinc+1);

    }
  
    
    

}

$beamdump_zdelta = 320;
$beamdump_zmin   = 0.;
$beamdump_zmax   = $beamdump_zmin + $beamdump_zdelta;
$beamdump_radius = 54.5;

my $water_at_the_end=29;
my $Al_at_the_end=76.2;

my $beamdump_start;
#my $beamdump_z = -($Bunker_dz+$beamdump_zdelta/2.)+2*$Bunker_dz-$beamdump_zdelta+145-91;

$beamdump_z = -229.5;


sub make_hallaBD
{
    
    
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = $beamdump_z;
    my $par1 = 0.;
    my $par2 = $beamdump_radius;
    my $par3 = $beamdump_zdelta/2.;
    my $par4 = 0.;
    my $par5 = 360.;
    
    $detector{"name"}  = "hallaBDextvessel";
    $detector{"mother"} = "Bunker_tunnel";
    $detector{"description"} = "BD ext vessel";
    $detector{"color"} = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Polycone";
    
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    my $dimen = "0.0*deg 360*deg 6*counts";
    $dimen = $dimen . " 0*cm 0*cm 0*cm 0*cm 0*cm 0*cm";
    $dimen = $dimen . " 29*cm 29*cm 53*cm 53*cm 2.28*inch 2.28*inch";
    $dimen = $dimen . " -8*cm 27*cm 47*cm 351*cm 380*cm 402*cm";
    
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = $dimen;
    $detector{"material"}    = "G4_WATER";
    print_det(\%configuration, \%detector);
    
    $X=0;
    $Y=0;
    $Z=402+$beamdump_z+0.2;
    
    $par2=2.28;
    $par3=0.1;
   

    %detector = init_det();
    $detector{"name"}        = "mutest_flux_bdvessel";
    $detector{"mother"}      = "Bunker_tunnel";
    $detector{"description"} = "Mu Test flux detector to sample the muon flux just out of the bd vessel";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*inch $par3*cm $par4*deg $par5*deg";
    #$detector{"sensitivity"} = "flux";
    #$detector{"hit_type"}    = "flux";
    $detector{"material"}    = "G4_WATER";
    #$detector{"identifiers"} = "id manual 4000";
    print_det(\%configuration, \%detector);
    
    
    
    
    
    
    
    
    
    %detector = init_det();
    
    $par3=$beamdump_zdelta/2.;
    $Z=$par3;
    $par2=25.2;
    $detector{"name"}        = "hallaBD";
    $detector{"mother"}      = "hallaBDextvessel";
    $detector{"description"} = "BD vessel";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_WATER";
    print_det(\%configuration, \%detector);
    $X = 0. ;
    $Y = 0. ;
    $Z = ($beamdump_zdelta/2.-$water_at_the_end-$Al_at_the_end)+$Al_at_the_end/2;
    
    %detector = init_det();
    $par1 = 0.;
    $par2 = 25.1;
    $par3 = $Al_at_the_end/2.;
    $par4 = 0.;
    $par5 = 360.;
    $detector{"name"}        = "hallaBD-Al-solid";
    $detector{"mother"}      = "hallaBD";
    $detector{"description"} = "Al solid dump";
    $detector{"color"}       = "0000FF";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);
    
    my $totAlPlates=0;
    my $first_plate_Z_center=-154.18;
    my $thisZ=$first_plate_Z_center;
    my $iz=0;
    my $plateID=0;
    #first 3 plates are .5inch
    my $thick=0.5;
    for ($iz=0;$iz<3;$iz++){
        %detector = init_det();
        $par3=$thick/2;
        $detector{"name"}        = "hallaBD-Al-plates_$plateID";$plateID++;
        $detector{"mother"}      = "hallaBD";
        $detector{"description"} = "Al plates $iz ";
        $detector{"color"}       = "0000ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $thisZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*inch $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
        $thisZ=$thisZ+2*$thick*2.54;
        $totAlPlates=$totAlPlates+$thick;
    }
    #next 3 plates are .375inch. The spacer before the first is .375inch
    $thisZ=$thisZ+(-0.5*2+0.5/2+0.375+0.375/2)*2.54;
    $thick=0.375;
    for ($iz=0;$iz<3;$iz++){
        %detector = init_det();
        $par3=$thick/2;
        $detector{"name"}        = "hallaBD-Al-plates_$plateID";$plateID++;
        $detector{"mother"}      = "hallaBD";
        $detector{"description"} = "Al plates $iz ";
        $detector{"color"}       = "0000ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $thisZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*inch $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
        $thisZ=$thisZ+$thick*2*2.54;
        $totAlPlates=$totAlPlates+$thick;
    }
    #next 8 plates are .25 inch. The spacer before the first is .25 inch
    $thisZ=$thisZ+(-.375*2+.375/2+.25+.25/2)*2.54;
    $thick=0.25;
    for ($iz=0;$iz<8;$iz++){
        %detector = init_det();
        $par3=$thick/2;
        $detector{"name"}        = "hallaBD-Al-plates_$plateID";$plateID++;
        $detector{"mother"}      = "hallaBD";
        $detector{"description"} = "Al plates $iz ";
        $detector{"color"}       = "0000ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $thisZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*inch $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
        $thisZ=$thisZ+$thick*2*2.54;
        $totAlPlates=$totAlPlates+$thick;
    }
    #next 50 plates are .188 inch. The spacer before the first is .188 inch
    $thisZ=$thisZ+(-.25*2+.25/2+.188+.188/2)*2.54;
    $thick=0.188;
    for ($iz=0;$iz<50;$iz++){
        %detector = init_det();
        $par3=$thick/2;
        $detector{"name"}        = "hallaBD-Al-plates_$plateID";$plateID++;
        $detector{"mother"}      = "hallaBD";
        $detector{"description"} = "Al plates $iz ";
        $detector{"color"}       = "0000ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $thisZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*inch $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
        $thisZ=$thisZ+$thick*2*2.54;
        $totAlPlates=$totAlPlates+$thick;
    }
    #next 16 plates are .25 inch. The spacer before the first is .25 inch
    $thisZ=$thisZ+(-.188*2+.188/2+.25+.25/2)*2.54;
    $thick=0.25;
    for ($iz=0;$iz<16;$iz++){
        %detector = init_det();
        $par3=$thick/2;
        $detector{"name"}        = "hallaBD-Al-plates_$plateID";$plateID++;
        $detector{"mother"}      = "hallaBD";
        $detector{"description"} = "Al plates $iz ";
        $detector{"color"}       = "0000ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $thisZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*inch $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
        $thisZ=$thisZ+$thick*2*2.54;
        $totAlPlates=$totAlPlates+$thick;
    }
    
    #next 9 plates are .375 inch. The spacer before the first is .375 inch
    $thisZ=$thisZ+(-.25*2+.25/2+.375+.375/2)*2.54;
    $thick=0.375;
    for ($iz=0;$iz<9;$iz++){
        %detector = init_det();
        $par3=$thick/2;
        $detector{"name"}        = "hallaBD-Al-plates_$plateID";$plateID++;
        $detector{"mother"}      = "hallaBD";
        $detector{"description"} = "Al plates $iz ";
        $detector{"color"}       = "0000ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $thisZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*inch $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
        $thisZ=$thisZ+$thick*2*2.54;
        $totAlPlates=$totAlPlates+$thick;
    }
    #next 12 plates are .5 inch. The spacer before the first is .5 inch
    $thisZ=$thisZ+(-.375*2+.375/2+.5+.5/2)*2.54;
    $thick=0.5;
    for ($iz=0;$iz<12;$iz++){
        %detector = init_det();
        $par3=$thick/2;
        $detector{"name"}        = "hallaBD-Al-plates_$plateID";$plateID++;
        $detector{"mother"}      = "hallaBD";
        $detector{"description"} = "Al plates $iz ";
        $detector{"color"}       = "0000ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $thisZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*inch $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
        $thisZ=$thisZ+$thick*2*2.54;
        $totAlPlates=$totAlPlates+$thick;
    }
    #next 12 plates are .75 inch. The spacer before the first is .5 inch, as well as the spacers in the middle
    $thisZ=$thisZ+(-0.5*2+0.5/2+0.5+0.75/2)*2.54;
    $thick=0.75;
    for ($iz=0;$iz<12;$iz++){
        %detector = init_det();
        $par3=$thick/2;
        $detector{"name"}        = "hallaBD-Al-plates_$plateID";$plateID++;
        $detector{"mother"}      = "hallaBD";
        $detector{"description"} = "Al plates $iz ";
        $detector{"color"}       = "0000ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $thisZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*inch $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
        $thisZ=$thisZ+$thick*2.54+0.5*2.54;
        $totAlPlates=$totAlPlates+$thick;
    }
    #next 6 plates are .1 inch. The spacer before the first is .5 inch, as well as the spacers in the middle
    $thisZ=$thisZ+(-.75/2+1/2)*2.54;
    $thick=1.;
    for ($iz=0;$iz<6;$iz++){
        %detector = init_det();
        $par3=$thick/2;
        $detector{"name"}        = "hallaBD-Al-plates_$plateID";$plateID++;
        $detector{"mother"}      = "hallaBD";
        $detector{"description"} = "Al plates $iz ";
        $detector{"color"}       = "0000ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $thisZ*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*inch $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Al"; 
        print_det(\%configuration, \%detector); 
        $thisZ=$thisZ+$thick*2.54+0.5*2.54;
        $totAlPlates=$totAlPlates+$thick;
    }
    
    
    # print "Total Al plates: $plateID. Tot Al (inches): $totAlPlates","\n";
    
    
} 




sub make_hallaBD_flux_barrel
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = 0.;
    my $par1 = $beamdump_radius*1.1;
    my $par2 = $beamdump_radius*1.1+0.1;
    my $par3 = $beamdump_zdelta/2.+2.5;
    my $par4 = 0.;
    my $par5 = 360.;
    
    $detector{"name"}        = "hallaBD_flux_barrel";
    $detector{"mother"}      = "Bunker_tunnel";
    $detector{"description"} = "Beamdump flux detector";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "Component";
    print_det(\%configuration, \%detector);
}
sub make_hallaBD_flux_endcup
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z =    $beamdump_zdelta/2.+2.5;
    # my $Z = $beamdump_z+$beamdump_zdelta/2.+2.5+0.05;

    my $par1 = 0.;
    my $par2 = $beamdump_radius*1.1+0.1;
    my $par3 = 0.1;
    my $par4 = 0.;
    my $par5 = 360.;
    
    $detector{"name"}        = "hallaBD_flux_endcup";
    $detector{"mother"}      = "Bunker_tunnel";
    $detector{"description"} = "Beamdump flux detector";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "Component";
    print_det(\%configuration, \%detector);
}
sub make_hallaBD_flux
{
    my %detector = init_det();
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = $beamdump_z;
    
    $detector{"name"}        = "hallaBD_flux";
    $detector{"mother"}      = "Bunker_tunnel";
    $detector{"description"} = "Beamdump flux detector";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Operation: hallaBD_flux_barrel + hallaBD_flux_endcup";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"material"}    = "G4_AIR";
    $detector{"identifiers"} = "id manual 0";
    #print_det(\%configuration, \%detector);
}

sub make_hallaBD_flux_sample
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z =    $beamdump_zdelta/2.+2.5;
    $Z = -$beamdump_zdelta/2.+4.;
    
    my $par1 = 0.;
    my $par2 = $beamdump_radius*0.95;
    my $par3 = 0.1;
    my $par4 = 0.;
    my $par5 = 360.;
    
    $detector{"name"}        = "hallaBD_flux_sample";
    $detector{"mother"}      = "hallaBD";
    $detector{"description"} = "Beamdump flux detector to sample em shower in BD";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"material"}    = "G4_AIR";
    $detector{"identifiers"} = "id manual 5000";
    #   print_det(\%configuration, \%detector);
}


my $mutest_pipe_x = 0.;
my $mutest_pipe_y = 0.;
#my $mutest_pipe_z = 960.;
my $mutest_pipe_InRad = 24./2;
my $mutest_pipe_thick = 1./2;
my $mutest_pipe_length = 1500./2;



## WELL1 - NOMINAL##
my $mutest_pipe_z = 2257.1;
if (($configuration{"variation"} eq "BDXmini")){
    $mutest_pipe_z=0;
    $mutest_pipe_length=$BDXmini_externalBox_Y;
}


## WELL2 - NOMINAL##
my $mutest_pipe_z2 = 2561.9;



sub make_mutest_pipe
{
    my %detector = init_det();
    
     my $X = 0. ;
     my $Y = 0. ;
     my $Z = $mutest_pipe_z;
    #my $Z = 0.;# to center the pipe
    print  " Pipe beam-height position (cm) X1=",$X,"\n";
    print  " Pipe beam-height position position (cm) Y=",$Y,"\n";
    print  " Pipe beam-height position position (cm) Z=",$Z,"\n";

    my $par1 = $mutest_pipe_InRad;
    my $par2 = $mutest_pipe_InRad+$mutest_pipe_thick;
    my $par3  =$mutest_pipe_length;
    my $par4 = 0.;
    my $par5 = 360.;
    
    $detector{"name"}        = "mutest_pipe";
    $detector{"mother"}      = "dirt";
    $detector{"description"} = "pipe to host muon flux detector";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "mutest_pipe_air";
    $detector{"mother"}      = "mutest_pipe";
    $detector{"description"} = "air in pipe ";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par1*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}

sub make_mutest_pipe2
{
    my %detector = init_det();
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = $mutest_pipe_z2;
    #my $Z = 0.;# to center the pipe
    
    my $par1 = $mutest_pipe_InRad;
    my $par2 = $mutest_pipe_InRad+$mutest_pipe_thick;
    my $par3  =$mutest_pipe_length;
    my $par4 = 0.;
    my $par5 = 360.;
    
    $detector{"name"}        = "mutest_pipe2";
    $detector{"mother"}      = "dirt";
    $detector{"description"} = "pipe 2 to host muon flux detector";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Al";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "mutest_pipe_air2";
    $detector{"mother"}      = "mutest_pipe2";
    $detector{"description"} = "air in pipe 2";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par1*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
}


sub make_mutest_flux
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = 0. ;

    
    my $par1 = 0.;
    my $par2 = 23./2;
    my $par3 = 0.5/2;
    my $par4 = 0.;
    my $par5 = 360.;
    
    $detector{"name"}        = "mutest_flux";
    $detector{"mother"}      = "mutest_pipe_air";
    $detector{"description"} = "Mu Test flux detector to sample the muon flux in the pipe";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"material"}    = "G4_AIR";
    $detector{"identifiers"} = "id manual 6000";
    print_det(\%configuration, \%detector);
}


my  $Muon_absorber_dz = 810./2.; # Half length on muon iron
my  $Muon_absorber_dx = 132.;    # Half width of muon absorber iron
my  $Muon_absorber_zmax = $Bunker_zmax + $Muon_absorber_dz*2.;
my  $Muon_absorber_c1_dz = 40./2.; # half length of upstream concrete block
my  $Muon_absorber_iron_dz = 660./2.; # half length of iron
my  $Muon_absorber_c2_dz = $Muon_absorber_dz-$Muon_absorber_c1_dz-$Muon_absorber_iron_dz; # half length of upstream concrete block

sub make_muon_absorber
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = $Bunker_zmax+$Muon_absorber_dz;
    my $par1 = $Muon_absorber_dx;
    my $par2 = $Muon_absorber_dx;
    my $par3 = $Muon_absorber_dz;
    
    $detector{"name"}        = "muon_absorber";
    $detector{"mother"}      = "dirt";
    $detector{"description"} = "Muon absorber";
    $detector{"color"}       = "a0a0a0";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Dirt";
    print_det(\%configuration, \%detector);

    $par3 = $Muon_absorber_c1_dz;
    $Z    =-$Muon_absorber_dz + $Muon_absorber_c1_dz;
    $detector{"name"}        = "muon_absorber_concrete1";
    $detector{"mother"}      = "muon_absorber";
    $detector{"description"} = "Muon absorber upstream concrete";
    $detector{"color"}       = "cccccc";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
    #    print "#####################################","\n";
    
    # print $par3,"\n";
    #print "#####################################","\n";

    $par3 = $Muon_absorber_iron_dz;
    $Z    =-$Muon_absorber_dz + $Muon_absorber_c1_dz*2. + $Muon_absorber_iron_dz;
    $detector{"name"}        = "muon_absorber_iron";
    $detector{"mother"}      = "muon_absorber";
    $detector{"description"} = "Muon absorber iron";
    $detector{"color"}       = "A05030";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Iron";
    print_det(\%configuration, \%detector);
    #  print "#####################################","\n";

    #  print $par3,"\n";
    #  print "#####################################","\n";
 
    
    $par3 = $Muon_absorber_c2_dz;
    $Z    = $Muon_absorber_dz - $Muon_absorber_c2_dz;
    $detector{"name"}        = "muon_absorber_concrete2";
    $detector{"mother"}      = "muon_absorber";
    $detector{"description"} = "Muon absorber upstream concrete";
    $detector{"color"}       = "cccccc";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
    # print "#####################################","\n";
    
    # print $par3,"\n";
    #print "#####################################","\n";
   
    # add flux detectors
    my $absorber_nflux=20;
    my $absorber_flux_dz=40.;
    my $absorber_flux_lx=130.;
    my $absorber_flux_ly=130.;
    my $absorber_flux_lz=0.01;
    for(my $iz=0; $iz<$absorber_nflux; $iz++)
    {
	$X = 0.;
        $Y = 0.;
        $Z=$iz*$absorber_flux_dz-$Muon_absorber_dz+$absorber_flux_lz+1.;
        my %detector = init_det();
        $detector{"name"}        = "absorber_flux_$iz";
	if($Z>=-$Muon_absorber_dz && $Z<=-$Muon_absorber_dz+$Muon_absorber_c1_dz*2.) { # flux detector is in upstream concrete block
	  $Z = $Z-(-$Muon_absorber_dz+$Muon_absorber_c1_dz);
	  $detector{"mother"}      = "muon_absorber_concrete1";
	  $detector{"material"}    = "BDX_Concrete";
        }
	elsif($Z>-$Muon_absorber_dz+$Muon_absorber_c1_dz*2. && $Z<$Muon_absorber_dz-$Muon_absorber_c2_dz*2.) { # flux detector is in iron block
	  $Z = $Z-(-$Muon_absorber_dz+$Muon_absorber_c1_dz*2+$Muon_absorber_iron_dz);
	  $detector{"mother"}      = "muon_absorber_iron";
	  $detector{"material"}    = "BDX_Iron";
        }
	else {
	  $Z = $Z-($Muon_absorber_dz-$Muon_absorber_c2_dz);
	  $detector{"mother"}      = "muon_absorber_concrete2";
	  $detector{"material"}    = "BDX_Concrete";
	}
        $detector{"description"} = "absorber flux detector $iz ";
        $detector{"color"}       = "cc00ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$absorber_flux_lx*cm $absorber_flux_ly*cm $absorber_flux_lz*cm";
        $detector{"sensitivity"} = "flux";
        $detector{"hit_type"}    = "flux";
        my $nflux=201+$iz;
        $detector{"identifiers"} = "id manual $nflux";
        # print_det(\%configuration, \%detector);
    }

}
    my $Building_dx = 900/2. ;
    my $Building_dy = 300/2. ;     # Headroom in detector building is 3m?
    my $Building_dz = 900/2. ;

    my $Building_cc_thick = 30 ; # Concrete walls are 30cm thick?

    my $Building_x_offset = $Building_dx - 200 ;

    my $Building_shaft_dx = 300/2. ; # 3x4.5 m shaft
    my $Building_shaft_dz = 450/2. ;
    my $Building_shaft_dy = ($Dirt_ymax/2.  - $Building_dy/2.) ;

    my $Building_shaft_offset_x = $Building_dx - $Building_shaft_dx ;
    my $Building_shaft_offset_y = $Building_dy + $Building_shaft_dy ;
    my $Building_shaft_offset_z = -$Building_dz + $Building_shaft_dz ;

sub make_det_house_outer
{
    my %detector = init_det();
    
    my $X = $Building_x_offset ;
    my $Y = 0. ;
    my $Z = $Muon_absorber_zmax + ($Building_dz + $Building_cc_thick);
    my $par1 = $Building_dx+$Building_cc_thick;
    my $par2 = $Building_dy+$Building_cc_thick;
    my $par3 = $Building_dz+$Building_cc_thick;
    
    $detector{"name"}        = "Det_house_outer";
    $detector{"mother"}      = "dirt";
    $detector{"description"} = "Outer envelope of detector house";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
}
my $roomx = $Building_dx;
my $roomy = $Building_dy+$Building_cc_thick/2.;
my $roomz = $Building_dz;

sub make_det_house_inner
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0.+$Building_cc_thick/2. ;
    my $Z = 0.;
    
    $detector{"name"}        = "Det_house_inner";
    $detector{"mother"}      = "Det_house_outer";
    $detector{"description"} = "Inner envelope of detector house";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$roomx*cm $roomy*cm $roomz*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
}




sub make_det_shaft_outer
{
    my %detector = init_det();

    
    my $X = $Building_x_offset+ $Building_shaft_offset_x;
    my $Y = 0. +$Building_shaft_offset_y+$Building_cc_thick/2.;
    my $Z = $Muon_absorber_zmax + ($Building_dz + $Building_cc_thick)+$Building_shaft_offset_z;


    my $par1 = $Building_shaft_dx+$Building_cc_thick;
    my $par2 = $Building_shaft_dy-$Building_cc_thick/2.;
    my $par3 = $Building_shaft_dz+$Building_cc_thick;
    
    $detector{"name"}        = "Det_shaft_outer";
    $detector{"mother"}      = "dirt";
    $detector{"description"} = "Outer envelope of shaft";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Concrete";
     print_det(\%configuration, \%detector);
}
sub make_det_shaft_inner
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = 0.;
    
    my $par1 = $Building_shaft_dx;
    my $par2 = $Building_shaft_dy-$Building_cc_thick/2.;
    my $par3 = $Building_shaft_dz;
    
    $detector{"name"}        = "Det_shaft_inner";
    $detector{"mother"}      = "Det_shaft_outer";
    $detector{"description"} = "Inner envelope of shaft";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}





my  $Building_stair_dz=$Building_dz-$Building_shaft_dz ;
my $Building_stair_dx =$Building_dz;
my $Building_stair_dy = ($Dirt_ymax/2.  - $Building_dy/2.) ;

my $Building_stair_offset_x = $Building_dx - $Building_stair_dx ;
my $Building_stair_offset_y = $Building_dy + $Building_stair_dy ;
my $Building_stair_offset_z = -$Building_dz + $Building_stair_dz+2*$Building_shaft_dz ;



sub make_stair_outer
{
    my %detector = init_det();
    
    
    my $X = $Building_x_offset+ $Building_stair_offset_x;
    my $Y = 0. +$Building_stair_offset_y +$Building_cc_thick/2.;
    my $Z = $Muon_absorber_zmax + ($Building_dz + 2*$Building_cc_thick)+$Building_stair_offset_z;
   
    
    my $par1 = $Building_stair_dx+$Building_cc_thick;
    my $par2 = $Building_stair_dy-$Building_cc_thick/2.;
    my $par3 = $Building_stair_dz;
    
    $detector{"name"}        = "Stair_outer";
    $detector{"mother"}      = "dirt";
    $detector{"description"} = "Outer stair";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
}
sub make_stair_inner
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = 0.;
    
    my $par1 = $Building_stair_dx;
    my $par2 = $Building_stair_dy-$Building_cc_thick/2.;
    my $par3 = $Building_stair_dz-$Building_cc_thick;
   
    $detector{"name"}        = "Stair_inner";
    $detector{"mother"}      = "Stair_outer";
    $detector{"description"} = "Inner stair";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}

my $strwallposz=(-$Building_dz/2+$Building_stair_dz+3/2.*$Building_cc_thick);

#print $strwallposz;

sub make_stair_wall
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = $strwallposz;
    
    my $par1 = $Building_stair_dx;
    my $par2 = $Building_dy+$Building_cc_thick/2.;
    my $par3 = $Building_cc_thick/2;
    
    $detector{"name"}        = "Stair_wall";
    $detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Stair wall ";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
}
sub make_stair_wall_door
{
    my %detector = init_det();
    
    my $door_l = 80. ;
    my $door_h = 250. ;
    my $door_p = -50. ;
    my $X =  $door_p;
    my $Y = -(2*($Building_dy+$Building_cc_thick/2.)-$door_h)/2;
    my $Z = 0.;
    
    my $par1 = $door_l/2.;
    my $par2 = $door_h/2;
    my $par3 = $Building_cc_thick/2;
    
    $detector{"name"}        = "Stair_wall_door";
    $detector{"mother"}      = "Stair_wall";
    $detector{"description"} = "Stair gate";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}
sub make_shaft_wall
{
    my %detector = init_det();
    
    my $X = +$Building_shaft_offset_x-$Building_shaft_dx-$Building_cc_thick/2.;
    my $Y = 0. ;
    my $Z = 0.-$Building_dz/2.+$Building_cc_thick/2.;
    
    my $par1 = $Building_cc_thick/2;
    my $par2 = $Building_dy+$Building_cc_thick/2.;
    my $par3 = $Building_shaft_dz+$Building_cc_thick/2.;
    
    $detector{"name"}        = "Shaft_wall";
    $detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Shaft wall";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
}
sub make_shaft_wall_door
{
    my %detector = init_det();
    
    my $door_l = 300. ;
    my $door_h = 150. ;
    my $door_p = -50. ;
    my $X = 0.;
    my $Y =  -(2*($Building_dy+$Building_cc_thick/2.)-$door_h)/2;
    my $Z = $door_p;
    
    my $par1 = $Building_cc_thick/2;
    my $par2 = $door_h/2;
    my $par3 = $door_l/2.;
    
    $detector{"name"}        = "Shaft_wall_door";
    $detector{"mother"}      = "Shaft_wall";
    $detector{"description"} = "Shaft gate";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}


my $ext_house_shift_dz=$Muon_absorber_zmax + ($Building_dz + $Building_cc_thick)-$Dirt_zmin - ($Dirt_zmax+$Dirt_zmin)/2.+100;
my $ext_house_dy=300./2.;
sub make_ext_house_outer
{
    my %detector = init_det();
    
    my $X = $Building_x_offset ;
    my $Y = ($Dirt_ymax-$Dirt_ymin)/2.+ $ext_house_dy+$Building_cc_thick;
    my $Z = $ext_house_shift_dz;
    my $par1 = $Building_dx+$Building_cc_thick;
    my $par2 = $ext_house_dy+$Building_cc_thick;
    my $par3 = $Building_dz+$Building_cc_thick;
    
    $detector{"name"}        = "ext_house_outer";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "Outer envelope of externalr building";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
}
sub make_ext_house_inner
{
    my %detector = init_det();
    
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = 0.;
    my $par1 = $Building_dx;
    my $par2 = $ext_house_dy;
    my $par3 = $Building_dz;
    
    $detector{"name"}        = "ext_house_inner";
    $detector{"mother"}      = "ext_house_outer";
    $detector{"description"} = "Inner envelope of external building";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}
sub make_ext_house_shaft_hole
{
    my %detector = init_det();
    
    my $X = $Building_dx-$Building_shaft_dx ;
    my $Y = -$ext_house_dy-$Building_cc_thick/2.;
    my $Z = -($Building_stair_dz-$Building_cc_thick);
    my $par1 = $Building_shaft_dx;
    my $par2 = $Building_cc_thick/2.;
    my $par3 = $Building_dz-($Building_stair_dz-$Building_cc_thick);
    
    $detector{"name"}        = "ext_house_shaft_hole";
    $detector{"mother"}      = "ext_house_outer";
    $detector{"description"} = "Shaft hall in the ext building floor";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}
sub make_ext_house_stair_hole
{
    my %detector = init_det();
    
    my $X =  $Building_stair_offset_x;
    my $Y = -$ext_house_dy-$Building_cc_thick/2.;
    my $Z =  $Building_dz-$Building_stair_dz+$Building_cc_thick;
    my $par1 = $Building_stair_dx;
    my $par2 = $Building_cc_thick/2.;
    my $par3 = $Building_stair_dz-$Building_cc_thick;
    
    $detector{"name"}        = "ext_house_stair_hole";
    $detector{"mother"}      = "ext_house_outer";
    $detector{"description"} = "Shaft hall in the ext building floor";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}
sub make_stair_steps_1
{
    my %detector = init_det();
    
    my $X =  0.;
    my $Y = 50.;
    my $Z =  -($Building_stair_dz-$Building_cc_thick)/2.;
    my $par1 = $Building_stair_dx;
    my $par2 = $Building_cc_thick/2.;
    my $par3 = ($Building_stair_dz-$Building_cc_thick)/2.;
    my $par4 = 0. ;
    my $par5 = 0. ;
    my $par6 = 0. ;
    
    $detector{"name"}        = "stair_steps_1";
    $detector{"mother"}      = "Stair_inner";
    $detector{"description"} = "Shaft hall in the ext building floor";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Parallelepiped";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0.*deg 0.*deg 30*deg";
   
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg $par6*deg";

    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
}
sub make_stair_steps_2
{
    my %detector = init_det();
    
    my $X =  300.;
    my $Y = -220.;
    my $Z =  ($Building_stair_dz-$Building_cc_thick)/2.;
    my $par1 = 0.25*$Building_stair_dx;
    my $par2 = $Building_cc_thick/2.;
    my $par3 = ($Building_stair_dz-$Building_cc_thick)/2.;
    my $par4 = 0. ;
    my $par5 = 0. ;
    my $par6 = 0. ;
    
    $detector{"name"}        = "stair_steps_2";
    $detector{"mother"}      = "Stair_inner ";
    $detector{"description"} = "Shaft hall in the ext building floor";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Parallelepiped";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0.*deg 0.*deg -30*deg";
    
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg $par6*deg";
    
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
}
sub make_stair_steps_3
{
    my %detector = init_det();
    
    my $X =  -80.;
    my $Y =  0.;
    my $Z =  2.*($Building_stair_dz-$Building_cc_thick)-38.;
    my $par1 = 0.67*$Building_stair_dx;
    my $par2 = $Building_cc_thick/2.;
    my $par3 = ($Building_stair_dz-$Building_cc_thick)/2.;
    my $par4 = 0. ;
    my $par5 = 0. ;
    my $par6 = 0. ;
    
    $detector{"name"}        = "stair_steps_3";
    $detector{"mother"}      = "Det_house_inner ";
    $detector{"description"} = "Shaft hall in the ext building floor";
    $detector{"color"}       = "A0A0A0";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Parallelepiped";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0.*deg 0.*deg -30*deg";
    
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg $par6*deg";
    
    $detector{"material"}    = "BDX_Concrete";
    print_det(\%configuration, \%detector);
}

#################################################################################################
#
# End: Hall-A beam dump
#
#################################################################################################

#################################################################################################
#
# Start: BDX-p veto
#
#################################################################################################


# inner veto
my $cormo_iveto_gap=1.;
my $cormo_iveto_tn=1./2.;
my $cormo_iveto_lx=20.;
my $cormo_iveto_ly=20.;
my $cormo_iveto_lz=52.9;

# lead shield
my $cormo_leadshield_tn=5./2.;
my $cormo_leadshield_gap=1.;
my $cormo_leadshield_lx=$cormo_iveto_lx+$cormo_leadshield_gap+2*$cormo_leadshield_tn;
my $cormo_leadshield_ly=$cormo_iveto_ly+$cormo_leadshield_gap+2*$cormo_iveto_tn;
my $cormo_leadshield_lz=$cormo_iveto_lz+$cormo_leadshield_gap+2*$cormo_leadshield_tn;

# outer veto
my $cormo_oveto_gap=1.;
my $cormo_oveto_tn=1.;
my $cormo_oveto_lx=$cormo_leadshield_lx+$cormo_oveto_gap+2*$cormo_oveto_tn;
my $cormo_oveto_ly=$cormo_leadshield_ly+$cormo_oveto_gap+2*$cormo_leadshield_tn;
my $cormo_oveto_lz=$cormo_leadshield_lz+$cormo_oveto_gap+2*$cormo_oveto_tn;
my $cormo_oveto_dz=$cormo_oveto_lz/2.;

my $cormo_z=0.;
my $cormo_box_lx=$cormo_iveto_lx-$cormo_iveto_tn;
my $cormo_box_ly=$cormo_iveto_ly-$cormo_iveto_tn;
my $cormo_box_lz=$cormo_iveto_lz-$cormo_iveto_tn;




# Start inner veto
# UP=1
# Bottom=3
# Downstream (Z bigger)=0
# Upstream (Z smaller or negative)=2
# Right (looking at the beam from the front - from Z positive)=5
# Left (looking at the beam from the front - from Z positive)=4

# shift for jlab location (wrt the hall)
my $jshiftx=-$Building_x_offset;
my $jshifty=0;
my $jshiftz=-$Building_cc_thick/2.-($Muon_absorber_zmax + ($Building_dz + $Building_cc_thick));



# starting from 0,0,0 define the different shifts
 my $shX=0.;
 my $shY=0.;
 my $shZ=0.;
if (($flag_JlabCT) eq ("1"))
{   $shX=-$Building_x_offset;
    $shY=0.;
    $shZ=-($Building_dz-$Building_cc_thick/2)+240;}


sub make_cormo_iveto
{

    my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
    
    $detector{"name"}        = "cormo_iveto_top";
    $detector{"description"} = "inner veto top";
    $detector{"color"}       = "0000FF";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $cormo_iveto_lx=42.8/2;;#RUN1=40.1/2; RUN2=42.8/2
    my $cormo_iveto_ly=1.0/2.;
    my $cormo_iveto_lz=98.5/2;#RUN1=105.8/2; RUN2=98.5/2
    my $X = $shX + 0.;
    my $Y = $shY + (35.1+1.0)/2.+0.5;
    my $Z = $shZ + (105.8-98.5)/2 + 7;#RUN1=0.; RUN2=+(105.8-98.5)/2
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$cormo_iveto_lx*cm $cormo_iveto_ly*cm $cormo_iveto_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 1 channel manual 1";
    print_det(\%configuration, \%detector);

    $detector{"name"}        = "cormo_iveto_bottom";
    $detector{"description"} = "inner veto bottom";
    $cormo_iveto_lx=40.1/2;#RUN1=42.8/2; RUN2=40.1/2
    $cormo_iveto_ly=1.0/2.;
    $cormo_iveto_lz=105.8/2.;#RUN1=98.5/2; RUN2=105.8/2
    $X = $shX + 0.;
    $Y = $shY -((35.1+1.0)/2+0.2);
    $Z = $shZ +0. + 7;#RUN1=-(105.8-98.5)/2; RUN2=0.
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_iveto_lx*cm $cormo_iveto_ly*cm $cormo_iveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 1 channel manual 2";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "cormo_iveto_upstream";
    $detector{"description"} = "inner veto upstream";
    $cormo_iveto_lx=40.8/2;
    $cormo_iveto_ly=34.6/2.;
    $cormo_iveto_lz=1.0/2.;
    $X = $shX +0.;
    $Y = $shY -(35.1-34.6)/2;
    $Z = $shZ -(105.8-1.0)/2.+10.+ 7;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_iveto_lx*cm $cormo_iveto_ly*cm $cormo_iveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 1 channel manual 3";
    print_det(\%configuration, \%detector);

    #$detector{"name"}        = "cormo_iveto_downstream-full";
    #$detector{"description"} = "inner veto downstream";
    #my $cormo_iveto_lx=40.8/2;
    #my $cormo_iveto_ly=34.6/2.;
    #my $cormo_iveto_lz=1.0/2.;
    #my $X = 0.;
    #my $Y = -(35.1-34.6)/2;
    #my $Z = +(105.8-1.0-(105.8-103.2))/2.-15.;
    #$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    #$detector{"dimensions"}  = "$cormo_iveto_lx*cm $cormo_iveto_ly*cm $cormo_iveto_lz*cm";
    #$detector{"identifiers"} = "sector manual 0 veto manual 1 channel manual 2";
    #print_det(\%configuration, \%detector);
    $detector{"name"}        = "cormo_iveto_downstream";
    $detector{"description"} = "inner veto downstream";
    $cormo_iveto_lx=40.8/2;
    $cormo_iveto_ly=(34.6-2.00)/2;
    $cormo_iveto_lz=1.0/2.;
    $X = $shX + 0.;
    $Y = $shY -(35.1-34.6)/2+2.00/2.;
    $Z =  $shZ +(105.8-1)/2.-1+ 7;#RUN1=-(105.8-103.2)/2.+103.2/2.+0.5-4.5; RUN2=+(105.8-1)/2.-1
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_iveto_lx*cm $cormo_iveto_ly*cm $cormo_iveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 1 channel manual 4";
    print_det(\%configuration, \%detector);
    #$detector{"name"}        = "cormo_iveto_downstream2";
    #$detector{"description"} = "inner veto downstream";
    #my $cormo_iveto_lx=(40.8-12.)/2;
    #my $cormo_iveto_ly=34.6/2.-(34.6-2.00)/2.;
    #my $cormo_iveto_lz=1.0/2.;
    #my $X = 0.+12./2.;
    #my $Y = -(35.1-34.6)/2-(34.6-2.00)/2.;
    #my $Z = +(105.8-1.0-(105.8-103.2))/2.-15.;
    #$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    #$detector{"dimensions"}  = "$cormo_iveto_lx*cm $cormo_iveto_ly*cm $cormo_iveto_lz*cm";
    #$detector{"identifiers"} = "sector manual 0 veto manual 1 channel manual 2";
    #print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "cormo_iveto_right";
    $detector{"description"} = "inner veto right";
    $cormo_iveto_lx=1.0/2;
    $cormo_iveto_ly=35.1/2.;
    $cormo_iveto_lz=103.2/2.;
    $X = $shX -(42.8-1.0)/2.;
    $Y = $shY + 0.;
    $Z = $shZ +(105.8-103.2)/2.+ 7;# RUN1 - (105.8-103.2)/2.
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_iveto_lx*cm $cormo_iveto_ly*cm $cormo_iveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 1 channel manual 5";
    print_det(\%configuration, \%detector);

    $detector{"name"}        = "cormo_iveto_left";
    $detector{"description"} = "inner veto left";
    $X = $shX +(42.8-1.0)/2.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_iveto_tn*cm $cormo_iveto_ly*cm $cormo_iveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 1 channel manual 6";
    print_det(\%configuration, \%detector);
}
# END inner veto
# Lead shield
sub make_cormo_lead
{
    my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}

    $detector{"name"}        = "cormo_lead_upstream";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $cormo_leadshield_lx=50./2;
    my $cormo_leadshield_ly=50./2.;
    my $cormo_leadshield_lz=5.0/2.;
    my $X = $shX + 0.;
    my $Y = $shY - ((35.1+1.0+5.0+1.0)/2+0.2)-5/2+50/2.;
    my $Z = $shZ -(105.8+5.0)/2.-2.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$cormo_leadshield_lx*cm $cormo_leadshield_ly*cm $cormo_leadshield_lz*cm";
    $detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "G4_WATER";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "cormo_lead_downstream";
    $X = $shX + 0.;
    $Z = $shZ +120+5.-(105.8+5.0)/2.+2.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_leadshield_lx*cm $cormo_leadshield_ly*cm $cormo_leadshield_lz*cm";
    print_det(\%configuration, \%detector);

    
    $detector{"name"}        = "cormo_lead_bottom";
    $detector{"style"}       = 0;
    $cormo_leadshield_lx=40./2;
    $cormo_leadshield_ly=5.0/2.;
    $cormo_leadshield_lz=120/2.;
    $X = $shX + 0.;
    $Y = $shY -((35.1+1.0+5.0+1.0)/2+0.2);
    $Z = $shZ -(105.8-120)/2.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_leadshield_lx*cm $cormo_leadshield_ly*cm $cormo_leadshield_lz*cm";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "cormo_lead_top";
    $cormo_leadshield_lx=50./2;
    $cormo_leadshield_ly=5.0/2.;
    $cormo_leadshield_lz=130/2.;
    $Y = $shY -((35.1+1.0+5.0+1.0)/2+0.2)+50.0+1.3;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_leadshield_lx*cm $cormo_leadshield_ly*cm $cormo_leadshield_lz*cm";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "cormo_lead_right";
    $detector{"style"}       = 0;
    $cormo_leadshield_lx=5.0/2.;
    $cormo_leadshield_ly=50./2.;
    $cormo_leadshield_lz=120/2.;
    $X = $shX  -(40+5.0)/2.-2.2;
    $Y = $shY -((35.1+1.0+5.0+1.0)/2+0.2)-5/2+50/2.;
    $Z = $shZ -(105.8-120)/2.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_leadshield_lx*cm $cormo_leadshield_ly*cm $cormo_leadshield_lz*cm";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "cormo_lead_left";
    $X = $shX  +(40+5.0)/2.+2.2;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_leadshield_lx*cm $cormo_leadshield_ly*cm $cormo_leadshield_lz*cm";
     print_det(\%configuration, \%detector);
}
# END lead shield

sub make_sarc_lead
{
    my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
 
    
    $detector{"name"}        = "sarc_lead_top";
    $detector{"description"} = "sarc lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $cormo_leadshield_lx=40./2;
    my $cormo_leadshield_ly=50./2.;
    my $cormo_leadshield_lz=120.0/2.;
    my $X = $shX + 0.;
    my $Y = $shY -((35.1+1.0+5.0+1.0)/2+0.2)+75.+30;
    my $Z = $shZ -(105.8-120)/2.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$cormo_leadshield_lx*cm $cormo_leadshield_ly*cm $cormo_leadshield_lz*cm";
    #$detector{"material"}    = "Borotron";
    $detector{"material"}    = "G4_AIR";
    #$detector{"material"}    = "G4_GRAPHITE";
    print_det(\%configuration, \%detector);

    my $sarc_x=0.;
    my $sarc_y=0.;
    my $sarc_z=0.;
    my $sarc_dx=5.5;
    my $sarc_dy=5;
    my $sarc_dz=20;
    my $sarc_th=5;

    $detector{"name"}        = "sarc_lead_out";
    $detector{"description"} = "sarc lead shield out";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    my $sarc_leadshield_lx=($sarc_dx+$sarc_th);
    my $sarc_leadshield_ly=($sarc_dy+$sarc_th);
    my $sarc_leadshield_lz=($sarc_dz+$sarc_th);
    $detector{"pos"}         = "$sarc_x*cm $sarc_y*cm $sarc_z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$sarc_leadshield_lx*cm $sarc_leadshield_ly*cm $sarc_leadshield_lz*cm";
    $detector{"material"}    = "Component";
    print_det(\%configuration, \%detector);

    
    $detector{"name"}        = "sarc_lead_in";
    $detector{"description"} = "sarc lead shield in";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    
    $sarc_leadshield_lx=$sarc_dx;
    $sarc_leadshield_ly=$sarc_dy;
    $sarc_leadshield_lz=$sarc_dz;
    
    $detector{"pos"}         = "$sarc_x*cm $sarc_y*cm $sarc_z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$sarc_leadshield_lx*cm $sarc_leadshield_ly*cm $sarc_leadshield_lz*cm";
    $detector{"material"}    = "Component";
    print_det(\%configuration, \%detector);
    
    $sarc_x=0.;
    $sarc_y=-3.8;
    $sarc_z=3.;

    $detector{"name"}        = "sarc_lead";
    $detector{"description"} = "sarc lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Operation:sarc_lead_out-sarc_lead_in";
    $detector{"pos"}         = "$sarc_x*cm $sarc_y*cm $sarc_z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "Borotron";
    print_det(\%configuration, \%detector);


}
# END lead shield


sub make_cormo_oveto
{
    my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}

    # Run 1 - fall 2015 - fall 2016
    ## Begin
    $detector{"name"}        = "cormo_oveto_top_upstream";
    $detector{"description"} = "outer veto top upstream";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $cormo_oveto_lx = 40./2;
    my $cormo_oveto_ly =2.0/2 ;
    my $cormo_oveto_lz =80./2 ;# MB May16 should be 80
    my $X = $shX + 0;
    my $Y = $shY -((35.1+1.0+5.0+1.0)/2+0.2)+50.0+1.3+5./2.+2/2;
    my $Z = $shZ -(105.8-120)/2.-80./2.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 1";
    #print_det(\%configuration, \%detector);
    $detector{"name"}        = "cormo_oveto_top_downstream";
    $detector{"description"} = "outer veto top downstream";
    $detector{"color"}       = "ff8000";
    $Z = $shZ -(105.8-120)/2.+80./2.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 2";
    #print_det(\%configuration, \%detector);
    ## End
    
    # Run 2 - fall 2016 - current
    ## Begin
    ## Only one scintillator replaces cormo_oveto_top_upstream
    ## sec 0 veto 2 ch 2 not used in this config
    $detector{"name"}        = "cormo_oveto_top_upstream";
    $detector{"description"} = "outer veto top upstream";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    
    $cormo_oveto_lx = 58./2;
    $cormo_oveto_ly =2.5/2 ;
    $cormo_oveto_lz =181./2 ;# MB May16 should be 80
    $X = $shX + 0;
    $Y = $shY -((35.1+1.0+5.0+1.0)/2+0.2)+50.0+1.3+5./2.+$cormo_oveto_ly;
    $Z = $shZ -(105.8-120)/2.-80./2.+40.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 1";
    print_det(\%configuration, \%detector);
    ## End
    
    $detector{"name"}        = "cormo_oveto_bottom_upstream";
    $detector{"color"}       = "ff8000";
    $detector{"description"} = "outer veto bottom upstream";
   
    $cormo_oveto_lx = 40./2;
    $cormo_oveto_ly =2.0/2 ;
    $cormo_oveto_lz =80./2 ;# MB May16 should be 80
    
    $X = $shX + 0;
    $Y = $shY -((35.1+1.0+5.0+1.0)/2+0.2)-5/2-7-$cormo_oveto_ly-3.5;#RUN1=-0; RUN2=-3.5
    $Z = $shZ -(105.8-120)/2.-80./2.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 3";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "cormo_oveto_bottom_downstream";
    $detector{"color"}       = "ff8000";
    $detector{"description"} = "outer veto bottom downstream";
    $Z = $shZ -(105.8-120)/2.+80./2.;#RUN1=-0; RUN2=+66./2
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 4";
    print_det(\%configuration, \%detector);

    

    $detector{"style"}       = 0;
    $detector{"name"}        = "cormo_oveto_upstream";
    $detector{"color"}       = "ff8000";
    $detector{"description"} = "outer veto upstream";
    $cormo_oveto_lx = 58./2; #RUN1=50./2; RUN2=58./2
    $cormo_oveto_ly =66./2 ; #RUN1=56./2; RUN2=66./2
    $cormo_oveto_lz =2.5/2 ;#RUN1=2; RUN2=-2.5
    $X = $shX + 0.;
    $Y = $shY -((35.1+1.0+5.0+1.0)/2+0.2)-5/2+66/2-10; #RUN1=+56/2; RUN2=+66./2
    $Z = $shZ -(105.8+5.0)/2.-2.5-2.5-$cormo_oveto_lz-12.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 5";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "cormo_oveto_downstream";
    $detector{"color"}       = "ff8000";
    $detector{"description"} = "outer veto downstream";
    $Z = $shZ + 120+5.-(105.8+5.0)/2.+2.5+2.5+$cormo_oveto_lz+12.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 6";
    print_det(\%configuration, \%detector);
    $detector{"style"}       = 0;
    
    $detector{"name"}        = "cormo_oveto_right1";
    $detector{"color"}       = "ff8000";
    $detector{"description"} = "outer veto Right";
    $cormo_oveto_lx = 2./2;
    $cormo_oveto_ly =80./2 ;
    $cormo_oveto_lz =40./2 ;
    $X = $shX  -(40+5.0)/2.-2.2-2.5-2/2-3.6;
    $Y = $shY -((35.1+1.0+5.0+1.0)/2+0.2)-5/2-11+80/2.-6.;#RUN1=-0; RUN2=-6
    $Z = $shZ -20 -(105.8+5.0)/2.-2.5-2.5+40/2+7.5;# 7.5 (and not 8.5) to make it suymmetric)
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 7";
    print_det(\%configuration, \%detector);
        $detector{"style"}       = 0;
    $detector{"name"}        = "cormo_oveto_right2";
    $Z = $shZ -20+ 40.-(105.8+5.0)/2.-2.5-2.5+40/2+7.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 8";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "cormo_oveto_right3";
    $Z = $shZ -20+ 80.-(105.8+5.0)/2.-2.5-2.5+40/2+7.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 9";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "cormo_oveto_right4";
    $Z = $shZ -20+ 120.-(105.8+5.0)/2.-2.5-2.5+40/2+7.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 10";
    print_det(\%configuration, \%detector);

    $detector{"name"}        = "cormo_oveto_left1";
    $detector{"color"}       = "ff8000";
    $detector{"description"} = "outer veto left";
    $cormo_oveto_lx = 2./2;
    $cormo_oveto_ly =80./2 ;
    $cormo_oveto_lz =40./2 ;
    $X = $shX + -(-(40+5.0)/2.-2.2-2.5-2/2-3.6);
    $Y = $shY -((35.1+1.0+5.0+1.0)/2+0.2)-5/2-11+80/2.-6.;#RUN1=-0; RUN2=-6
    $Z = $shZ -20 -(105.8+5.0)/2.-2.5-2.5+40/2+7.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 11";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "cormo_oveto_left2";
    $Z = $shZ -20 +40.-(105.8+5.0)/2.-2.5-2.5+40/2+7.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 12";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "cormo_oveto_left3";
    $Z = $shZ -20 +80.-(105.8+5.0)/2.-2.5-2.5+40/2+7.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 13";
    print_det(\%configuration, \%detector);
        $detector{"name"}        = "cormo_oveto_left4";
    $Z = $shZ -20 +120.-(105.8+5.0)/2.-2.5-2.5+40/2+7.5;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$cormo_oveto_lx*cm $cormo_oveto_ly*cm $cormo_oveto_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 2 channel manual 14";
    print_det(\%configuration, \%detector);

}



# define cormorino shield geometry
# all sizes are in cm
my $cormo_hole_r=50.;
my $cormo_hole_h=50.;

my $cormo_shield_tn_r=100.;
my $cormo_shield_tn_h=100.;

my $cormo_shield_r=$cormo_hole_r+$cormo_shield_tn_r;
my $cormo_shield_h=$cormo_hole_h+$cormo_shield_tn_h;

my $cormo_shield_nplanes=6;
my @cormo_shield_ir = (              0.,              0.,   $cormo_hole_r,   $cormo_hole_r,              0.,              0.);
my @cormo_shield_or = ( $cormo_shield_r, $cormo_shield_r, $cormo_shield_r, $cormo_shield_r, $cormo_shield_r, $cormo_shield_r);
my @cormo_shield_z  = (-$cormo_shield_h,  -$cormo_hole_h,  -$cormo_hole_h,   $cormo_hole_h,   $cormo_hole_h, $cormo_shield_h);

sub make_cormo_shield
{
    my %detector = init_det();
    $detector{"name"}        = "cormo_shield";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "cormorad shield";
    $detector{"color"}       = "BDBDBD";
    $detector{"style"}       = 0;
    $detector{"type"}        = "Polycone";
    my $X = 0.;
    my $Y = 0.;
    my $Z = $cormo_z;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    my $dimen = "0.0*deg 360*deg $cormo_shield_nplanes*counts";
    for(my $i = 0; $i <$cormo_shield_nplanes ; $i++)
    {
        $dimen = $dimen ." $cormo_shield_ir[$i]*cm";
    }
    for(my $i = 0; $i <$cormo_shield_nplanes ; $i++)
    {
        $dimen = $dimen ." $cormo_shield_or[$i]*cm";
    }
    for(my $i = 0; $i <$cormo_shield_nplanes ; $i++)
    {
        $dimen = $dimen ." $cormo_shield_z[$i]*cm";
    }
    $detector{"dimensions"} = $dimen;
    $detector{"material"}    = "BDX_Iron";
    print_det(\%configuration, \%detector);
}
#############
# COSMIC Sphere
# CT fulldet     cosmicradius=110.
# CT cryst only  cosmicradius=50.
# Proposal       cosmicradius=180. Z = +1950
# Proposal       single crystal = +1830
# Proposal       mutest = +1830
sub make_flux_cosmic_sph
{
    my %detector = init_det();
    my $cosmicradius=50.;
    if (($configuration{"variation"} eq "CT") or ($flag_mutest eq 1) or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
    my $X = 0.;
    my $Y = 0.;
    my $Z = 0.;
    
     if (($configuration{"variation"} eq "Proposal") and ($flag_JlabCT eq "0"))
    {$X=-$Building_x_offset;
     $Y=0.;
     $Z=-($Building_dz-$Building_cc_thick/2)+240;}
    
    if (($flag_minibdx eq "1"))
    {#print $flag_detcentered,"\n";
        $X = $shX + $mutest_pipe_x;
        $Y = $shY + $mutest_pipe_y;
        $Z = $shZ + $mutest_pipe_z;
        
        if (($flag_detcentered eq "1"))
        {#print $flag_detcentered,"\n";
            $X = $shX + 0. ;
            $Y = $shY + 0.;
            $Z = $shZ + 0.;
        }
        
        print "Place the cosmic shere at X=",$X," Y=",$Y," Z=",$Z," with Radius R=",$cosmicradius,"(cm)","\n";
    }
    
 
    my $par1 = $cosmicradius;
    my $par2 = $cosmicradius+0.01;
    my $parz3 = 0.;
    my $par4 = 360.;
    my $par5 = 0.;
    my $par6 = 90.;
    
    
    $detector{"name"}        = "flux_cosmic_sph";
    $detector{"description"} = "Beamdump flux detector";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Sphere";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $parz3*cm $par4*deg $par5*deg $par6*deg";
    $detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 1001";
    print_det(\%configuration, \%detector);
    
}




#############
# COSMIC cylinder (CT)
# CT fulldet     cosmicradius=110.
# CT cryst only  cosmicradius=50.
# Proposal       cosmicradius=180. Z = +1950
# Proposal       single crystal = +1830
# Proposal       mutest = +1830
sub make_flux_cosmic_cyl
{
    my %detector = init_det();
    my $cosmicradius=110.1;
    $cosmicradius=20.1;#mutest
    if (($configuration{"variation"}) eq ("CT")  or ($flag_mutest eq 1) or ($flag_detcentered eq "1") )
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
    
    if ((($configuration{"variation"}) eq ("CT")) || (($configuration{"variation"} eq "Proposal") and ($flag_JlabCT eq "1")))
    {$cosmicradius=110.1;} #fixed cosmic radius for proto conf both jlab/ct

    my $X = $shX + 0. ;
    my $Y = $shY + 0. ;
    my $Z = $shZ +  1822;#130.
    my $par1 = $cosmicradius;
    my $par2 = $cosmicradius+0.01;
    my $parz3 = ($cosmicradius/2);
    my $par4 = 0.;
    my $par5 = 360.;


    $detector{"name"}        = "flux_cosmic_cyl";
    $detector{"description"} = "Beamdump flux detector";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $parz3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 1001";
    print_det(\%configuration, \%detector);

    $X = $shX + 0. ;
    $Y = $shY + $parz3 ;

    $par1 =0;
    $par2 = $cosmicradius;
    my $par3 = 0.01;
    $par4 = 0.;
    $par5 = 360.;

    $detector{"name"}        = "flux_cosmic_cyl_top";
    $detector{"description"} = "Beamdump flux detector";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 1002";
    print_det(\%configuration, \%detector);
    
    $X = $shX +  0. ;
    $Y = $shY -$parz3;

    $detector{"name"}        = "flux_cosmic_cyl_bot";
    $detector{"description"} = "Beamdump flux detector";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 1003";
    print_det(\%configuration, \%detector);

}


#############
# COSMIC rec
# CT fulldet     cosmicradius=50.
# CT cryst only  cosmicradius=20.
# Proposal       cosmicradius=65. Z = +240; Z = +1950 from the hall
# Proposal       single crystal = +1830
sub make_flux_cosmic_rec
{

    my %detector = init_det();
    my $cosmicradius=65.;
    
   
    
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
    if ((($configuration{"variation"}) eq ("CT")) || (($configuration{"variation"} eq "Proposal") and ($flag_JlabCT eq "1")))
    {$cosmicradius=44.9;} #fixed cosmic radius for proto conf both jlab/ct



    my $X = $shX +  0. ;
    my $Y = $shY + $cosmicradius;
    my $Z = $shZ + 0;#130.
    
    if (($configuration{"variation"} eq "Proposal") and ($flag_JlabCT eq "0"))
    {$X=-$Building_x_offset;
    $Z=-($Building_dz-$Building_cc_thick/2)+240;}

    my $par1 = $cosmicradius;
    my $par2 = 0.01;
    my $par3 = 3*$cosmicradius;
    
    
    $detector{"name"}        = "flux_cosmic_rec_top";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "cc00aa";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 2001";
    print_det(\%configuration, \%detector);
    
    $Y = $shY +  -$cosmicradius;
    $detector{"name"}        = "flux_cosmic_rec_bottom";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "cc00ff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "flux";
    $detector{"hit_type"}    = "flux";
    $detector{"identifiers"} = "id manual 2002";
    print_det(\%configuration, \%detector);


}

########
##
##  BaBar Crystals
# ALL in cm
# Averaging the crystal size
# Endcap: short side X (4.3+3.9)/2=4.1cm
my $cr_ssx=4.1/2 ;
# Endcap: short side Y 4.7
my $cr_ssy=4.7/2 ;
# Endcap: long side X (5+4.6)/2=4.8cm
my $cr_lsx=4.8/2 ;
# Endcap: long side Y 5.4
my $cr_lsy=5.4/2 ;
# Endcap: lenght side Y 32.5
my $cr_lgt=32.5/2.;
# Mylar wrapping thikness
my $cr_mylar=0.005/2;
# Wrapping thikness
my $cr_airgap=0.1/2;
# Alveolus thikness (0.03cm if Cfiber, 0.2 cm if Al)
my $cr_alv=0.2/2;
######
# G.Ottonello measurement
#  short side X
$cr_ssx=4.7/2 ;
# short side Y
$cr_ssy=4.8/2 ;
# long side X
$cr_lsx=5.8/2 ;
#  long side Y
$cr_lsy=6.0/2 ;
#  lenght side Y
$cr_lgt=31.6/2.;
# Mylar wrapping thikness
$cr_mylar=0.005/2;
# Wrapping thikness
$cr_airgap=0.1/2;
# Alveolus thikness (0.03cm if Cfiber, 0.2 cm if Al)
$cr_alv=0.2/2;

# Wrapped crystals
my $wr_cr_ssx=$cr_ssx+$cr_mylar;
my $wr_cr_ssy=$cr_ssy+$cr_mylar;
my $wr_cr_lsx=$cr_lsx+$cr_mylar;
my $wr_cr_lsy=$cr_lsy+$cr_mylar;
my $wr_cr_lgt=$cr_lgt+$cr_mylar;
# Air around crystals
my $ar_wr_cr_ssx=$wr_cr_ssx+$cr_airgap;
my $ar_wr_cr_ssy=$wr_cr_ssy+$cr_airgap;
my $ar_wr_cr_lsx=$wr_cr_lsx+$cr_airgap;
my $ar_wr_cr_lsy=$wr_cr_lsy+$cr_airgap;
my $ar_wr_cr_lgt=$wr_cr_lgt+$cr_airgap;
# Crystal alveolus
my $al_ar_wr_cr_ssx=$ar_wr_cr_ssx+$cr_alv;
my $al_ar_wr_cr_ssy=$ar_wr_cr_ssy+$cr_alv;
my $al_ar_wr_cr_lsx=$ar_wr_cr_lsx+$cr_alv;
my $al_ar_wr_cr_lsy=$ar_wr_cr_lsy+$cr_alv;
my $al_ar_wr_cr_lgt=$ar_wr_cr_lgt+$cr_alv;

if (($configuration{"variation"}) eq ("CT"))
{
# fixed Al alveoles
$al_ar_wr_cr_ssx=8.2;
$al_ar_wr_cr_ssy=8./2;
$al_ar_wr_cr_lsx=8./2.;
$al_ar_wr_cr_lsy=8./2.;
$al_ar_wr_cr_lgt=(29.5+5+5)/2.;
}

#distance betweens modules (arbitraty fixed at XXcm)
my $blocks_distance=$cr_lgt*2.+2.;


# $fg=1 Flipped crystal positioning
# $fg=0 Unflipped crystal positioning
my $fg=0;
# Number of modules (blocks or sectors)
my $nblock=8;
# Nuumber of columns (vert or X)
my $ncol=10;
# Number of rows (horiz or Y)
my $nrow=10;
#  <----- X Y ||
#             ||
#             ||
#             \/

# makes the alveoles parallelepipedal in shape (assuming that short sides are both < long sides)
my $irectalv=1 ;

if($irectalv==1) {
  $al_ar_wr_cr_ssy=$al_ar_wr_cr_lsy;
    #$al_ar_wr_cr_lsx=$al_ar_wr_cr_lsy;
  $al_ar_wr_cr_ssx=$al_ar_wr_cr_lsx;
  $ar_wr_cr_ssy=$ar_wr_cr_lsy;
    #$ar_wr_cr_lsx=$ar_wr_cr_lsy;
  $ar_wr_cr_ssx=$ar_wr_cr_lsx;
}


# to place it in in Hall-A detector house
# Size of the whole calorimeter
my $cal_sz_x=2*$al_ar_wr_cr_lsx*$ncol;
my $cal_sz_y=2*$al_ar_wr_cr_lsy*$nrow;
my $cal_sz_z=($nblock-1)*$blocks_distance+2*$al_ar_wr_cr_lgt;



my $det_dist_from_wall=110.;
my $shiftx=-$Building_x_offset;
my $shifty=-$Building_cc_thick/2.;
my $shiftz=-($Building_dz-$Building_cc_thick/2)+$det_dist_from_wall;

    if (($flag_detcentered) eq "1")
    {
        $shiftx=$shX;
        $shifty=$shY;
        $shiftz=-($cal_sz_z/2.-$al_ar_wr_cr_lgt);
    }

# CT configuration or # Prototype at JLab
if ((($configuration{"variation"}) eq ("CT")) || (($configuration{"variation"} eq "Proposal") and ($flag_JlabCT eq "1")))
{   $fg=0;
    $nblock=1;
    $ncol=4;#RUN1=1; RUN2=4
    $nrow=4;#RUN1=1; RUN2=4
    $shiftx=$shX - 1.;#RUN1=+0.; RUN2=+1.
    $shifty=$shY -((35.1+1.0)/2+0.2)+0.5+17.-$cr_lsx+2.;#RUN1=+0.; RUN2=+2.
    $shiftz=$shZ -(105.8-1.0)/2.+1./2+54+5;#RUN1=+0.; RUN2=-5.
}


my $tocntx=($ncol-1)/2.*($al_ar_wr_cr_lsx+$al_ar_wr_cr_lsx);
my $tocnty=($nrow-1)/2.*($al_ar_wr_cr_lsy+$al_ar_wr_cr_ssy);
#print $tocnty;

my $tol=0.;
if (($flag_inner_lead) eq "1") {$tol=15/2};

## IV
my $iv_tol=$tol; # 0 in case of outer lead
my $iv_thk=1.0/2.;
my $iv_lgt=100./2.+$iv_tol/4.;
my $n_iv=int($cal_sz_z/($iv_lgt*2)+1);
my $iv_wdtx=int($cal_sz_x/10.+1.)*10/2.+$iv_tol;
my $iv_wdty=int($cal_sz_y/10.+1.)*10/2.+$iv_tol;
# IV top and bottom width is increased by the thikness to have an hermetic arrangment
my $iv_wdtxH=$iv_wdtx+2*$iv_thk;
my $iv_nsipmx=2*int($iv_wdtxH/10.+1);
my $iv_nsipmy=2*int($iv_wdty/10.+1);
my $iv_nsipmu=2.;

# Cal center in X, Y, Z
my $cal_centx=($tocntx+$shiftx+$tocntx-($ncol-1)*2.*($al_ar_wr_cr_lsx+$al_ar_wr_cr_lsx)/2.+$shiftx)/2;
my $cal_centy=($tocnty+$shifty+$tocnty-($nrow-1)*2.*($al_ar_wr_cr_ssy+$al_ar_wr_cr_lsy)/2+$shifty)/2.;
my $cal_centz=$cal_sz_z/2.-$al_ar_wr_cr_lgt+$shiftz;

## Lead
my $ld_tol=5.0;
my $ld_thk=5./2.;
my $ld_lgt=$n_iv*$iv_lgt+$ld_tol/2.;
my $ld_wdx=$iv_wdtxH+$ld_thk+$ld_tol/2.;
my $ld_wdy=$iv_wdty+2*$iv_thk+$ld_tol/2.;
my $ld_lgtH=$ld_lgt+2*$ld_thk;

## Lead vault inside INNER veto
my $ld_in_tol=-15.0;
my $ld_in_thk=5./2.;
my $ld_in_lgt=$n_iv*$iv_lgt+$ld_in_tol/2;
my $ld_in_wdx=$iv_wdtxH+($ld_in_thk+$ld_in_tol/2.);
my $ld_in_wdy=$iv_wdty+2*$iv_thk+$ld_in_tol/2.;
my $ld_in_lgtH=$ld_in_lgt+2*$ld_in_thk;


## OV
my $ov_thk=2.0/2.;
my $ov_lgt=170./2.;
my $ov_wdtL=32./2.;
my $ov_szL=0;
my $n_ov=0;
my $ov_wdtx=0;
my $ov_wdty=0;
my $ov_wdtxH=0;
my $n_ovL=0;
my $ov_lgtL=0;
my $ov_wdxU=0;
my $ov_wdyU=0;
my $ov_zU=0;

if (($flag_outer_lead) eq "1")
    {
        $n_ov=int($ld_lgtH*2./($ov_lgt*2)+1);
        $ov_wdtx=int($ld_wdx/5.+1.)*5.;
        $ov_wdty=int($ld_wdy/5.+1.)*5.;
        $ov_szL=$ov_wdty+10.;

        # OV top and bottom width is increased by the thikness to have an hermetic arrangment
        $ov_wdtxH=$ov_wdtx;
        # OV laterals
        $n_ovL=int($n_ov*$ov_lgt/$ov_wdtL+1);
        $ov_lgtL=$n_ovL*$ov_wdtL;
        # OV Upstream/Downstream
        $ov_wdxU=$ov_wdtx;
        $ov_wdyU=$ld_wdy+2*$ld_thk;
        $ov_zU=($ld_lgtH+$ov_thk);
}
if (($flag_outer_lead) eq "0")
{
    $n_ov=int($n_iv*$iv_lgt*2./($ov_lgt*2)+1);
    $ov_wdtx=int($iv_wdtx/5.+1.)*5.;
    $ov_wdty=int($iv_wdty/5.+1.)*5.;
    $ov_szL=$ov_wdty+10.;
    # OV top and bottom width is increased by the thikness to have an hermetic arrangment
    $ov_wdtxH=$ov_wdtx;
    # OV laterals
    $n_ovL=int($n_ov*$ov_lgt/$ov_wdtL+1);
    $ov_lgtL=$n_ovL*$ov_wdtL;
    # OV Upstream/Downstream
    $ov_wdxU=$ov_wdtx;
    $ov_wdyU=$iv_wdty+2*$iv_thk;
    $ov_zU=($n_iv*$iv_lgt+$ov_thk);
}





#print " ",$iv_posz," ",$blocks_distance+$shiftz,"\n";
#print " ",$iv_posx," ",$tocntx+$shiftx," ",($tocntx+$shiftx)-$ncol*($al_ar_wr_cr_lsx+$al_ar_wr_cr_lsx)/2.,"\n";
#print " ",$iv_posy," ",$tocnty+$shifty," ",$tocnty-($nrow-1)*2.*($al_ar_wr_cr_ssy+$al_ar_wr_cr_lsy)/2+$shifty,"\n";
print "#####################################","\n";
print "###  Room size ############","\n";
print "###  ","\n";
print "###  Room size (m) (lg X Y): ",($roomz+$strwallposz)/100.," x ",2*$roomx/100.," x ",2*$roomy/100.,"\n";
print "#####################################","\n";
print "###  DETECTOR parameters ############","\n";
print "###  ","\n";
print "###       CALORIMETER","\n";
print "###  Module size and number (X,Y,Z): ",$ncol," ",$nrow," ",$nblock,"\n";
print "###  Total crystals: ",$ncol*$nrow*$nblock,"\n";
print "###  Distance between modules: ",$blocks_distance-$cr_lgt*2.,"\n";
print "###  Cal size (cm) (lgt X Y)= ",$cal_sz_z," x ",$cal_sz_x," x ",$cal_sz_y,"\n";
print "###  Cal center (cm) (X Y Z)= ",$cal_centx," ",$cal_centy," ",$cal_centz,"\n";
print "###       IV","\n";
print "###  IV overall lgt=",$n_iv*$iv_lgt*2.," in ",$n_iv," pieces","\n";
print "###  IV TOP/BOTTOM (cm) : ",2*$iv_lgt," x ",2*$iv_wdtxH," x ",2*$iv_thk,," , ",$n_iv," pieces","\n";
print "###  IV sides (cm) : ",2*$iv_lgt," x ",2*$iv_wdty," x ",2*$iv_thk,," , ",$n_iv," pieces","\n";
print "###  IV U/D (cm) : ",2*$iv_wdtx," x ",2*$iv_wdty," x ",2*$iv_thk,," , 2 pieces","\n";
if (($flag_outer_lead) eq "1")
{
print "###       LEAD","\n";
print "###  Lead TOP/BOTTOM (cm) : ",2*$ld_lgtH," x ",2*$ld_wdx," x ",2*$ld_thk,"\n";
print "###  Lead sides (cm) : ",2*$ld_lgt," x ",2*$ld_wdy," x ",2*$ld_thk,"\n";
print "###  Lead U/D (cm) : ",2*$ld_wdx," x ",2*$ld_wdy," x ",2*$ld_thk,"\n";
print "###  Lead tolerance (cm): ",$ld_tol,"\n";
print "###  Lead weight (ton): ",2*$ld_thk*(8*($ld_wdx*$ld_lgtH)+8*$ld_wdy*($ld_wdx+$ld_lgt))*11.3/1e6,"\n";
}
if (($flag_outer_lead) eq "1")
{
    print "###       Inner LEAD","\n";
    print "###  Lead TOP/BOTTOM (cm) : ",2*$ld_in_lgtH," x ",2*$ld_in_wdx," x ",2*$ld_in_thk,"\n";
    print "###  Lead sides (cm) : ",2*$ld_in_lgt," x ",2*$ld_in_wdy," x ",2*$ld_in_thk,"\n";
    print "###  Lead U/D (cm) : ",2*$ld_in_wdx," x ",2*$ld_in_wdy," x ",2*$ld_in_thk,"\n";
    print "###  Lead tolerance (cm): ",$ld_in_tol,"\n";
    print "###  Lead weight (ton): ",2*$ld_in_thk*(8*($ld_in_wdx*$ld_in_lgtH)+8*$ld_in_wdy*($ld_in_wdx+$ld_in_lgt))*11.3/1e6,"\n";
}
print "###       OV","\n";
print "###  OV overall lgt=",$n_ov*$ov_lgt*2.," in ",$n_ov," pieces","\n";
print "###  OV TOP/BOTTOM (cm) : ",2*$ov_lgt,"x",2*$ov_wdtx,"x",2*$ov_thk," , ",$n_ov,"x2 pieces","\n";
print "###  OV sides (cm) : ",2*$ov_szL,"x",2*$ov_wdtL,"x",2*$ov_thk," , ",$n_ovL,"x2 pieces","\n";
print "###  OV U/D (cm) : ",2*$ov_wdxU,"x",2*$ov_wdyU,"x",2*$ov_thk," , 2 pieces","\n";
print "###      Sizes","\n";
print "###  Cal x= ",$cal_sz_x," IV x=",2*$iv_wdtx*2.," Lead x=",2*$ld_wdx*2.," OV x=",2*$ov_wdtx*2.,"\n";
print "###  Cal y= ",$cal_sz_y," IV y=",2*$iv_wdty*2, " Lead y=",2*$ld_wdy*2.," OV y=",2*$ov_wdty*2.,"\n";
print "###  Cal z= ",$cal_sz_z," IV z=",$n_iv*2*$iv_lgt,  " Lead z=",2*$ld_lgtH, " OV z=",$n_ov*2*$ov_lgt,"\n";
print "###      Electronics","\n";
print "###  SiPM Cal (single readout 6x6 mm): ",$ncol*$nrow*$nblock,"\n";
print "###  SiPM IV (2 sides readout 3x3 mm): ",2*($n_iv*($iv_nsipmx+$iv_nsipmy)+$iv_nsipmu),"\n";
print "###  PMTs OV (1 side readout):",2*($n_ov+$n_ovL+$ov_thk),"\n";
print "###  ","\n";
print "###  END DETECTOR parameters ############","\n";
print "#########################################","\n";

#print "Lead x= ",$cal_sz_x," OVeto x=",$iv_wdtx*2.,"\n";
#print "Lead y= ",$cal_sz_y," OVeto y=",$iv_wdty*2,"\n";


#print "Shiftz ",$cal_sz_z,"\n";

 sub make_cry_module
{

    if ($configuration{"variation"} eq "Proposal")
    {
    # add flux detectors
    my $det_nflux=1;
    my $det_flux_dz=40.;
    my $det_flux_lx=130.;
    my $det_flux_ly=130.;
    my $det_flux_lz=0.01;
    for(my $iz=0; $iz<$det_nflux; $iz++)
    {
        my %detector = init_det();
        if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
        {$detector{"mother"}      = "main_volume";}
        else
        {$detector{"mother"}      = "Det_house_inner";}

        $detector{"name"}        = "det_flux_$iz";
        $detector{"description"} = "det flux detector $iz ";
        $detector{"color"}       = "cc00ff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = $cal_centx;
        my $Y = $cal_centy;
        my $Z=($iz+1)*$det_flux_dz+$det_flux_lz+1.-$Building_dz;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$det_flux_lx*cm $det_flux_ly*cm $det_flux_lz*cm";
        $detector{"material"}    = "G4_AIR";
        #$detector{"sensitivity"} = "flux";
        #$detector{"hit_type"}    = "flux";
        my $nflux=301+$iz;
        #$detector{"identifiers"} = "id manual $nflux";
        #print_det(\%configuration, \%detector);
        }
    }
    

    #    print $tocntx;
    # FLAG to rotate a plane of crystals
    # $rotated=1 rotate
    # $rotated=0 standard
    # change also $blocks_distance=$cr_lgt*2.+1.;
    my $rotated=1;

    if ($rotated == 1)
    {
        ######### Init rotated crystal
        for(my $im=0; $im<($nblock); $im++)
        {
            for(my $ib=0; $ib<($ncol); $ib++)
            {
                for(my $ir=0; $ir<($nrow); $ir++)
                {
                    my $rot90=90*($ir/2-int($ir/2))/0.5;
                    #print "AAA",$im," ",$ib," ",$ir,"\n";
                    my $rot=$fg*180.*((int(($ib+1.)/2.)-int(($ib)/2.))-(int(($ir+1.)/2.)-int(($ir)/2.)))+$rot90;#RUN1=+0; #RUN2=+180.
                    # Carbon/Aluminum alveols
                    my %detector = init_det();
                    $detector{"name"}        = "cry_alveol_$ib"."_"."$ir"."_"."$im";
                    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
                    {$detector{"mother"}      = "main_volume";}
                    else
                    {$detector{"mother"}      = "Det_house_inner";}
                    
                    
                    $detector{"description"} = "Carbon/Al container_$ib"."_"."$ir"."_"."$im";
                    $detector{"color"}       = "00ffff";
                    $detector{"style"}       = 0;
                    $detector{"visible"}     = 1;
                    $detector{"type"}        = "Trd";
                    my $X = $tocntx-$ib*2.*($al_ar_wr_cr_lsx+$al_ar_wr_cr_lsx)/2.+$shiftx;
                    my $Y = $tocnty-$ir*2.*($al_ar_wr_cr_ssy+$al_ar_wr_cr_lsy)/2+$shifty ;
                    my $Z = $im*$blocks_distance+$shiftz;
                    
                    if ($rot90 != 0)
                    {
                      
                    my $ik = ($im/2-int($im/2))/0.5;
                        #  $X=$ik*$blocks_distance+$shiftz-($al_ar_wr_cr_lgt)-7.+80;
                        #$Z=$tocntx-$ib*2.*($al_ar_wr_cr_lsx+$al_ar_wr_cr_lsx)/2.+$shiftx+($al_ar_wr_cr_lgt)+8+int($im/2)*2*($blocks_distance)-80;
                        $X=$ik*$blocks_distance+$shiftz-($al_ar_wr_cr_lgt)+117;
                        $Z=$tocntx-$ib*2.*($al_ar_wr_cr_lsx+$al_ar_wr_cr_lsx)/2.+$shiftx+($al_ar_wr_cr_lgt)+8+int($im/2)*2*($blocks_distance)-125;
                    }
                    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
                    $detector{"rotation"}    = "0*deg $rot*deg 0*deg";
                    my $par1 =$al_ar_wr_cr_ssx ;
                    my $par2 =$al_ar_wr_cr_lsx;
                    my $par3 =$al_ar_wr_cr_ssy  ;
                    my $par4 =$al_ar_wr_cr_lsy  ;
                    my $par5 =$al_ar_wr_cr_lgt ;
                    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
                    $detector{"material"}    = "G4_Al";
                    #print_det(\%configuration, \%detector);
                    #print  " cryst",$Z," ","\n";
                    
                    #print $detector{"name"},"\n";
                    #print $detector{"pos"},"\n";
                    
                    # Air layer
                    %detector = init_det();
                    $detector{"name"}        = "cry_air_$ib"."_"."$ir"."_"."$im";
                    $detector{"mother"}      = "cry_alveol_$ib"."_"."$ir"."_"."$im";
                    #$detector{"mother"}      = "Det_house_inner";
                    $detector{"description"} = "Air $ib"."_"."$ir"."_"."$im";
                    $detector{"color"}       = "00fff1";
                    $detector{"style"}       = 0;
                    $detector{"visible"}     = 1;
                    $detector{"type"}        = "Trd";
                    $X = 0.;
                    $Y = 0.;
                    $Z = 0.;
                    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
                    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
                    $par1 =$ar_wr_cr_ssx ;
                    $par2 =$ar_wr_cr_lsx;
                    $par3 =$ar_wr_cr_ssy  ;
                    $par4 =$ar_wr_cr_lsy  ;
                    $par5 =$ar_wr_cr_lgt ;
                    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
                    $detector{"material"}    = "G4_AIR";
                    #print_det(\%configuration, \%detector);
                    
                    # Mylar wrapping
                    %detector = init_det();
                    $detector{"name"}        = "cry_mylar_$ib"."_"."$ir"."_"."$im";
                    $detector{"mother"}      = "cry_air_$ib"."_"."$ir"."_"."$im";
                    #$detector{"mother"}      = "Det_house_inner";
                    $detector{"description"} = "Mylar wrapping_$ib"."_"."$ir"."_"."$im";
                    $detector{"color"}       = "00fff2";
                    $detector{"style"}       = 0;
                    $detector{"visible"}     = 1;
                    $detector{"type"}        = "Trd";
                    $X = 0.;
                    $Y = 0.;
                    $Z = 0.;
                    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
                    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
                    $par1 =$wr_cr_ssx ;
                    $par2 =$wr_cr_lsx;
                    $par3 =$wr_cr_ssy  ;
                    $par4 =$wr_cr_lsy  ;
                    $par5 =$wr_cr_lgt ;
                    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
                    $detector{"material"}    = "bdx_mylar";
                    #print_det(\%configuration, \%detector);
                    
                    
                    # Crystals
                    %detector = init_det();
                    $detector{"name"}        = "crystal_$ib"."_"."$ir"."_"."$im";
                    $detector{"mother"}      = "cry_mylar_$ib"."_"."$ir"."_"."$im";
                    #$detector{"mother"}      = "Det_house_inner";
                    $detector{"description"} = "Crystal_$ib"."_"."$ir"."_"."$im";
                    $detector{"color"}       = "00ffff";
                    $detector{"style"}       = 1;
                    $detector{"visible"}     = 1;
                    $detector{"type"}        = "Trd";
                    $X = 0.;
                    $Y = 0.;
                    $Z = 0.;
                    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
                    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
                    $par1 =$cr_ssx ;
                    $par2 =$cr_lsx;
                    $par3 =$cr_ssy  ;
                    $par4 =$cr_lsy  ;
                    $par5 =$cr_lgt ;
                    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
                    $detector{"material"}    = "CsI_Tl";
                    #$detector{"material"}    = "G4_AIR";
                    $detector{"sensitivity"} = "crs";
                    $detector{"hit_type"}    = "crs";
                    my $i_im=$im+1;
                    my $i_ir=$ir+1;
                    my $i_ib=$ib+1;
                    $detector{"identifiers"} = "sector manual $i_im xch manual $i_ir ych manual $i_ib";
                    #print_det(\%configuration, \%detector);
                    
                    
                    
                }
            }
        }
        
        ######### End rortated crystal
        

    }
    else{
    for(my $im=0; $im<($nblock); $im++)
    {
    for(my $ib=0; $ib<($ncol); $ib++)
    {
         for(my $ir=0; $ir<($nrow); $ir++)
        {
            #print "AAA",$im," ",$ib," ",$ir,"\n";
            my $rot=$fg*180.*((int(($ib+1.)/2.)-int(($ib)/2.))-(int(($ir+1.)/2.)-int(($ir)/2.)))+180.;#RUN1=+0; #RUN2=+180.
            # Carbon/Aluminum alveols
            my %detector = init_det();
            $detector{"name"}        = "cry_alveol_$ib"."_"."$ir"."_"."$im";
            if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
            {$detector{"mother"}      = "main_volume";}
            else
            {$detector{"mother"}      = "Det_house_inner";}
            

            $detector{"description"} = "Carbon/Al container_$ib"."_"."$ir"."_"."$im";
            $detector{"color"}       = "00ffff";
            $detector{"style"}       = 0;
            $detector{"visible"}     = 1;
            $detector{"type"}        = "Trd";
            my $X = $tocntx-$ib*2.*($al_ar_wr_cr_lsx+$al_ar_wr_cr_lsx)/2.+$shiftx;
            my $Y = $tocnty-$ir*2.*($al_ar_wr_cr_ssy+$al_ar_wr_cr_lsy)/2+$shifty ;
            my $Z = $im*$blocks_distance+$shiftz;
            $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
            $detector{"rotation"}    = "0*deg $rot*deg 0*deg";
            my $par1 =$al_ar_wr_cr_ssx ;
            my $par2 =$al_ar_wr_cr_lsx;
            my $par3 =$al_ar_wr_cr_ssy  ;
            my $par4 =$al_ar_wr_cr_lsy  ;
            my $par5 =$al_ar_wr_cr_lgt ;
            $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
            $detector{"material"}    = "G4_Al";
            print_det(\%configuration, \%detector);
            #print  " cryst",$Z," ","\n";
            
            #print $detector{"name"},"\n";
            #print $detector{"pos"},"\n";
            
            # Air layer
            %detector = init_det();
            $detector{"name"}        = "cry_air_$ib"."_"."$ir"."_"."$im";
            $detector{"mother"}      = "cry_alveol_$ib"."_"."$ir"."_"."$im";
            #$detector{"mother"}      = "Det_house_inner";
            $detector{"description"} = "Air $ib"."_"."$ir"."_"."$im";
            $detector{"color"}       = "00fff1";
            $detector{"style"}       = 0;
            $detector{"visible"}     = 1;
            $detector{"type"}        = "Trd";
            $X = 0.;
            $Y = 0.;
            $Z = 0.;
            $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
            $detector{"rotation"}    = "0*deg 0*deg 0*deg";
            $par1 =$ar_wr_cr_ssx ;
            $par2 =$ar_wr_cr_lsx;
            $par3 =$ar_wr_cr_ssy  ;
            $par4 =$ar_wr_cr_lsy  ;
            $par5 =$ar_wr_cr_lgt ;
            $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
            $detector{"material"}    = "G4_AIR";
            print_det(\%configuration, \%detector);
            
            # Mylar wrapping
            %detector = init_det();
            $detector{"name"}        = "cry_mylar_$ib"."_"."$ir"."_"."$im";
            $detector{"mother"}      = "cry_air_$ib"."_"."$ir"."_"."$im";
            #$detector{"mother"}      = "Det_house_inner";
            $detector{"description"} = "Mylar wrapping_$ib"."_"."$ir"."_"."$im";
            $detector{"color"}       = "00fff2";
            $detector{"style"}       = 0;
            $detector{"visible"}     = 1;
            $detector{"type"}        = "Trd";
            $X = 0.;
            $Y = 0.;
            $Z = 0.;
            $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
            $detector{"rotation"}    = "0*deg 0*deg 0*deg";
            $par1 =$wr_cr_ssx ;
            $par2 =$wr_cr_lsx;
            $par3 =$wr_cr_ssy  ;
            $par4 =$wr_cr_lsy  ;
            $par5 =$wr_cr_lgt ;
            $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
            $detector{"material"}    = "bdx_mylar";
            print_det(\%configuration, \%detector);
            
            
            # Crystals
            %detector = init_det();
            $detector{"name"}        = "crystal_$ib"."_"."$ir"."_"."$im";
            $detector{"mother"}      = "cry_mylar_$ib"."_"."$ir"."_"."$im";
            #$detector{"mother"}      = "Det_house_inner";
            $detector{"description"} = "Crystal_$ib"."_"."$ir"."_"."$im";
            $detector{"color"}       = "00ffff";
            $detector{"style"}       = 1;
            $detector{"visible"}     = 1;
            $detector{"type"}        = "Trd";
            $X = 0.;
            $Y = 0.;
            $Z = 0.;
            $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
            $detector{"rotation"}    = "0*deg 0*deg 0*deg";
            $par1 =$cr_ssx ;
            $par2 =$cr_lsx;
            $par3 =$cr_ssy  ;
            $par4 =$cr_lsy  ;
            $par5 =$cr_lgt ;
            $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
            $detector{"material"}    = "CsI_Tl";
            #$detector{"material"}    = "G4_AIR";
            $detector{"sensitivity"} = "crs";
            $detector{"hit_type"}    = "crs";
            my $i_im=$im+1;
            my $i_ir=$ir+1;
            my $i_ib=$ib+1;
            $detector{"identifiers"} = "sector manual $i_im xch manual $i_ir ych manual $i_ib";
            print_det(\%configuration, \%detector);
            

        
        }
         }
    }
    }
    
    
}

sub make_cry_module_single
{
    $fg=0;
    $nblock=1;
    $ncol=1;
    $nrow=1;
    $al_ar_wr_cr_lgt=(29.5+5)/2.;

                my $rot=$fg*180+90 ;                # Carbon/Aluminum alveols
                my %detector = init_det();
                $detector{"name"}        = "cry_alveol";
                if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
                {$detector{"mother"}      = "main_volume";}
                else
                {$detector{"mother"}      = "Det_house_inner";}
                $detector{"description"} = "Carbon/Al container";
                $detector{"color"}       = "00ffff";
                $detector{"style"}       = 0;
                $detector{"visible"}     = 1;
                $detector{"type"}        = "Trd";
                my $X = +1;
                my $Y = -13.5 ;
                my $Z = 38;
                $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
                $detector{"rotation"}    = "0*deg $rot*deg 0*deg";
                my $par1 =$al_ar_wr_cr_ssx ;
                my $par2 =$al_ar_wr_cr_lsx;
                my $par3 =$al_ar_wr_cr_ssy  ;
                my $par4 =$al_ar_wr_cr_lsy  ;
                my $par5 =$al_ar_wr_cr_lgt ;
                $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
                $detector{"material"}    = "G4_AIR";
                print_det(\%configuration, \%detector);
                #print  " cryst",$Z," ","\n";

    # Air layer
                %detector = init_det();
                $detector{"name"}        = "cry_air";
                $detector{"mother"}      = "cry_alveol";
                #$detector{"mother"}      = "Det_house_inner";
                $detector{"description"} = "Air ";
                $detector{"color"}       = "00fff1";
                $detector{"style"}       = 0;
                $detector{"visible"}     = 1;
                $detector{"type"}        = "Trd";
                $X = 0.;
                $Y = 0.;
                $Z = 0.;
                $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
                $detector{"rotation"}    = "0*deg 0*deg 0*deg";
                $par1 =$ar_wr_cr_ssx ;
                $par2 =$ar_wr_cr_lsx;
                $par3 =$ar_wr_cr_ssy  ;
                $par4 =$ar_wr_cr_lsy  ;
                $par5 =$ar_wr_cr_lgt ;
                $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
                $detector{"material"}    = "G4_AIR";
                print_det(\%configuration, \%detector);
                
                # Mylar wrapping
                %detector = init_det();
                $detector{"name"}        = "cry_mylar";
                $detector{"mother"}      = "cry_air";
                #$detector{"mother"}      = "Det_house_inner";
                $detector{"description"} = "Mylar wrapping";
                $detector{"color"}       = "00fff2";
                $detector{"style"}       = 0;
                $detector{"visible"}     = 1;
                $detector{"type"}        = "Trd";
                $X = 0.;
                $Y = 0.;
                $Z = 0.;
                $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
                $detector{"rotation"}    = "0*deg 0*deg 0*deg";
                $par1 =$wr_cr_ssx ;
                $par2 =$wr_cr_lsx;
                $par3 =$wr_cr_ssy  ;
                $par4 =$wr_cr_lsy  ;
                $par5 =$wr_cr_lgt ;
                $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
                $detector{"material"}    = "bdx_mylar";
                print_det(\%configuration, \%detector);
                
                
                # Crystals
                %detector = init_det();
                $detector{"name"}        = "crystal";
                $detector{"mother"}      = "cry_mylar";
                #$detector{"mother"}      = "Det_house_inner";
                $detector{"description"} = "Crystal";
                $detector{"color"}       = "00ffff";
                $detector{"style"}       = 1;
                $detector{"visible"}     = 1;
                $detector{"type"}        = "Trd";
                $X = 0.;
                $Y = 0.;
                $Z = 0.;
                $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
                $detector{"rotation"}    = "0*deg 0*deg 0*deg";
                $par1 =$cr_ssx ;
                $par2 =$cr_lsx;
                $par3 =$cr_ssy  ;
                $par4 =$cr_lsy  ;
                $par5 =$cr_lgt ;
                $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
                $detector{"material"}    = "CsI_Tl";
                #$detector{"material"}    = "G4_AIR";
                $detector{"sensitivity"} = "crs";
                $detector{"hit_type"}    = "crs";
                $detector{"identifiers"} = "sector manual 100 xch manual 0 ych manual 0";
                print_det(\%configuration, \%detector);

}




sub make_iveto
{
    my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
    
    for(my $ir=0; $ir<($n_iv); $ir++)
    {

    $detector{"name"}        = "iveto_top_$ir";
    $detector{"description"} = "inner veto top";
    $detector{"color"}       = "0000FF";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $X = $cal_centx;;
    my $Y = $cal_centy+$iv_wdty/2.+$iv_wdty/2.+$iv_thk;
    my $Z = $cal_centz-$n_iv*$iv_lgt+(2*$ir+1)*$iv_lgt;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$iv_wdtxH*cm $iv_thk*cm $iv_lgt*cm ";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual $ir veto manual 4 channel manual 1";
    print_det(\%configuration, \%detector);

        $detector{"name"}        = "iveto_bottom_$ir";
        #$detector{"mother"}      = "Det_house_inner";
        $detector{"description"} = "inner veto bottom";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = $cal_centx;;
        $Y = $cal_centy-$iv_wdty/2.-$iv_wdty/2.-$iv_thk;
        $Z = $cal_centz-$n_iv*$iv_lgt+(2*$ir+1)*$iv_lgt;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_wdtxH*cm $iv_thk*cm $iv_lgt*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $detector{"identifiers"} = "sector manual $ir veto manual 4 channel manual 2";
        print_det(\%configuration, \%detector);

        $detector{"name"}        = "iveto_right_$ir";
        #$detector{"mother"}      = "Det_house_inner";
        $detector{"description"} = "inner veto right";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = $cal_centx-$iv_wdtx/2.-$iv_wdtx/2.-$iv_thk;
        $Y = $cal_centy;
        $Z = $cal_centz-$n_iv*$iv_lgt+(2*$ir+1)*$iv_lgt;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_thk*cm $iv_wdty*cm $iv_lgt*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $detector{"identifiers"} = "sector manual $ir veto manual 4 channel manual 5";
        print_det(\%configuration, \%detector);

        $detector{"name"}        = "iveto_left_$ir";
        #$detector{"mother"}      = "Det_house_inner";
        $detector{"description"} = "inner veto left";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X =$cal_centx+$iv_wdtx/2.+$iv_wdtx/2.+$iv_thk;
        $Y =$cal_centy;
        $Z = $cal_centz-$n_iv*$iv_lgt+(2*$ir+1)*$iv_lgt;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_thk*cm $iv_wdty*cm $iv_lgt*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $detector{"identifiers"} = "sector manual $ir veto manual 4 channel manual 6";
        print_det(\%configuration, \%detector);
    
    }
    
        $detector{"name"}        = "iveto_downstream";
        $detector{"description"} = "inner veto downstream";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = $cal_centx;;
        my $Y = $cal_centy;
        my $Z = $cal_centz+$n_iv*$iv_lgt-$iv_thk;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_wdtx*cm $iv_wdty*cm $iv_thk*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $detector{"identifiers"} = "sector manual 0 veto manual 4 channel manual 4";
        print_det(\%configuration, \%detector);

    $detector{"name"}        = "iveto_upstream";
    $detector{"description"} = "inner veto upstream";
    $detector{"color"}       = "0000FF";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = $cal_centx;;
    $Y = $cal_centy;
    $Z = $cal_centz-$n_iv*$iv_lgt+$iv_thk;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$iv_wdtx*cm $iv_wdty*cm $iv_thk*cm ";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 4 channel manual 3";
    print_det(\%configuration, \%detector);

}
# END Proposal inner veto
# BEGIN Porposal LEAD
sub make_lead
{
    my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
    
    my $X = $cal_centx;
    my $Y = ($cal_centy+$iv_wdty/2.+$iv_wdty/2.+$iv_thk)+$iv_thk + $ld_thk+$ld_tol/2.;
    my $Z = $cal_centz;
    $detector{"name"}        = "lead_top";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_wdx*cm $ld_thk*cm $ld_lgtH*cm";
    $detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    $Y = ($cal_centy-$iv_wdty/2.-$iv_wdty/2.-$iv_thk-$iv_thk) - $ld_thk-$ld_tol/2.;
    $detector{"name"}        = "lead_bottom";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_wdx*cm $ld_thk*cm $ld_lgtH*cm";
    #$detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    $X = $cal_centx;
    $Y = $cal_centy;
    $Z = $cal_centz-$ld_lgtH+$ld_thk;
    $detector{"name"}        = "lead_upstream";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_wdx*cm  $ld_wdy*cm $ld_thk*cm";
    #$detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    $Z = $cal_centz+$ld_lgtH-$ld_thk;
    $detector{"name"}        = "lead_downstream";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_wdx*cm  $ld_wdy*cm $ld_thk*cm";
    #$detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    $X = $cal_centx-$ld_wdx+$ld_thk;
    $Y = $cal_centy;
    $Z = $cal_centz;
    $detector{"name"}        = "lead_right";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_thk*cm  $ld_wdy*cm $ld_lgt*cm";
    #$detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    $X = $cal_centx+$ld_wdx-$ld_thk;
    $Y = $cal_centy;
    $Z = $cal_centz;
    $detector{"name"}        = "lead_left";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_thk*cm  $ld_wdy*cm $ld_lgt*cm";
    #$detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}

# END Proposal LEAD

# BEGIN Porposal INNER LEAD
sub make_in_lead
{
    my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
    
    my $X = $cal_centx;
    my $Y = ($cal_centy+$iv_wdty/2.+$iv_wdty/2.+$iv_thk)+$iv_thk + $ld_in_thk + $ld_in_tol/2.;
    my $Z = $cal_centz;
    $detector{"name"}        = "lead_in_top";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_in_wdx*cm $ld_in_thk*cm $ld_in_lgtH*cm";
    $detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    $Y = ($cal_centy-$iv_wdty/2.-$iv_wdty/2.-$iv_thk)-$iv_thk - $ld_thk - $ld_in_tol/2.;
    $detector{"name"}        = "lead_in_bottom";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_in_wdx*cm $ld_in_thk*cm $ld_in_lgtH*cm";
    #$detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    $X = $cal_centx;
    $Y = $cal_centy;
    $Z = $cal_centz-$ld_in_lgtH+$ld_in_thk;
    $detector{"name"}        = "lead_in_upstream";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_in_wdx*cm  $ld_in_wdy*cm $ld_in_thk*cm";
    #$detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    $Z = $cal_centz+$ld_in_lgtH-$ld_in_thk;
    $detector{"name"}        = "lead_in_downstream";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_in_wdx*cm  $ld_in_wdy*cm $ld_in_thk*cm";
    #$detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    $X = $cal_centx-$ld_in_wdx+$ld_in_thk;
    $Y = $cal_centy;
    $Z = $cal_centz;
    $detector{"name"}        = "lead_in_right";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_in_thk*cm  $ld_in_wdy*cm $ld_in_lgt*cm";
    #$detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    $X = $cal_centx+$ld_in_wdx-$ld_in_thk;
    $Y = $cal_centy;
    $Z = $cal_centz;
    $detector{"name"}        = "lead_in_left";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$ld_in_thk*cm  $ld_in_wdy*cm $ld_in_lgt*cm";
    #$detector{"material"}    = "G4_Pb";
    #$detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
}

# END Proposal INNER LEAD



# BEGIN Porposal OV
sub make_oveto
{
    my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
    
    for(my $ir=0; $ir<($n_ov); $ir++)
    {
        
        $detector{"name"}        = "oveto_top_$ir";
        $detector{"description"} = "Outer veto top";
        $detector{"color"}       = "088A4B";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = $cal_centx;;
        my $Y = $cal_centy+$ov_wdyU+$ov_thk;
        my $Z = $cal_centz-$n_ov*$ov_lgt+(2*$ir+1)*$ov_lgt;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$ov_wdtx*cm $ov_thk*cm $ov_lgt*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $detector{"identifiers"} = "sector manual $ir veto manual 5 channel manual 1";
        print_det(\%configuration, \%detector);
        
        $detector{"name"}        = "oveto_bottom_$ir";
        $detector{"description"} = "Outer veto bottom";
        $detector{"color"}       = "088A4B";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = $cal_centx;;
        $Y = $cal_centy-$ov_wdyU-$ov_thk;
        $Z = $cal_centz-$n_ov*$ov_lgt+(2*$ir+1)*$ov_lgt;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$ov_wdtx*cm $ov_thk*cm $ov_lgt*cm ";  
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $detector{"identifiers"} = "sector manual $ir veto manual 5 channel manual 2";
        print_det(\%configuration, \%detector);
        
        
    }
       for(my $ir=0; $ir<($n_ovL); $ir++)
    {
        my $ir_rev=$n_ovL-$ir-1;
        $detector{"name"}        = "oveto_left_$ir_rev";
        $detector{"description"} = "Outer veto Right";
        $detector{"color"}       = "088A4B";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = $cal_centx+$ov_wdtx+$ov_thk;
        my $Y = $cal_centy;
        my $Z = $cal_centz+($n_ovL-1)*$ov_wdtL-2*$ov_wdtL*$ir;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$ov_thk*cm $ov_szL*cm $ov_wdtL*cm";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $detector{"identifiers"} = "sector manual $ir_rev veto manual 5 channel manual 4";
        print_det(\%configuration, \%detector);
        
        #print  "OV  SECTOR L=",$ir_rev," Position in Z = ",$Z,"\n";

        
        $detector{"name"}        = "oveto_right_$ir_rev";
        $detector{"description"} = "Outer veto Left";
        $detector{"color"}       = "088A4B";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = $cal_centx-$ov_wdtx-$ov_thk;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$ov_thk*cm $ov_szL*cm $ov_wdtL*cm";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $detector{"identifiers"} = "sector manual $ir_rev veto manual 5 channel manual 3";
        print_det(\%configuration, \%detector);
        #print  "OV  SECTOR R=",$ir_rev," Position in Z = ",$Z,"\n";
        
        
    }
    $detector{"name"}        = "oveto_downstream";
    $detector{"description"} = "Outer veto Downstream";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $X = $cal_centx;
    my $Y = $cal_centy;
    my $Z = $cal_centz+$ov_zU;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = " $ov_wdxU*cm $ov_wdyU*cm $ov_thk*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 5 channel manual 5";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "oveto_upstream";
    $detector{"description"} = "Outer veto Upstream";
    $detector{"color"}       = "088A4B";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $Z = $cal_centz-$ov_zU;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = " $ov_wdxU*cm $ov_wdyU*cm $ov_thk*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 5 channel manual 6";
    print_det(\%configuration, \%detector);
    

}
# END Proposal OV
sub make_csi_pad
{
    my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
    
    $detector{"name"}        = "csi_pad_up";
    $detector{"description"} = "paddle over the crystal";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $csi_pad_lx =12./2.;
    my $csi_pad_ly =1.0/2 ;
    my $csi_pad_lz =12./2 ;
    my $X = $shiftx+ 0;
    my $Y = $shifty-$cr_lsy-1.5-0.5+0.5+20.;
    my $Z = $shiftz+$cr_lgt-$csi_pad_lz-8.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 4 channel manual 1";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "csi_pad_down";
    $detector{"description"} = "baddle below the crystal";
    $detector{"type"}        = "Box";
    $csi_pad_lx =12./2.;
    $csi_pad_ly =1.0/2 ;
    $csi_pad_lz =12./2 ;
    $X = $shiftx + 0;
    $Y = $shifty-$cr_lsy-1.5-0.5;
    $Z = $shiftz+$cr_lgt-$csi_pad_lz-8.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 4 channel manual 2";
    print_det(\%configuration, \%detector);
}
sub make_cal_pad
{
    my %detector = init_det();
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {$detector{"mother"}      = "main_volume";}
    else
    {$detector{"mother"}      = "Det_house_inner";}
    
    $detector{"name"}        = "cal_pad_up";
    $detector{"description"} = "paddle over the calorimeter";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $csi_pad_lx =5/2.;
    my $csi_pad_ly =1.0/2 ;
    my $csi_pad_lz =32./2 ; # TO BE CHECKED
    my $X = $shiftx;
    my $Y = $shifty+37.5;
    my $Z = $shiftz;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 3 channel manual 1";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "cal_pad_down";
    $detector{"description"} = "paddle below the calorimeter";
    $detector{"type"}        = "Box";
    $X = $shiftx ;
    $Y = $shifty-35.5;
    $Z = $shiftz;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 3 channel manual 2";
    print_det(\%configuration, \%detector);
}

sub make_mutest_detector
{
    my %detector = init_det();
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = 100.; # depth in the pipe
    my $rotX=0.;

    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {
        $detector{"mother"}      = "main_volume";
        $rotX=0.; #0 for H 90deg for V
        $Z = 0.
    }
    else
    {$detector{"mother"}      = "mutest_pipe_air";}

    
         my $vs_top_tk=0.5;
         my $vs_side_tk=1.0;
         my $vs_lg=26.0*2.;
         my $vs_ir=20.;
        my $par1 = $vs_ir/2; #InRad
        my $par2 = ($vs_ir+$vs_side_tk)/2; #InRad+thick
    
        my $par3  =$vs_lg/2.;#length
        my $par4 = 0.;
        my $par5 = 360.;
        
        $detector{"name"}        = "mutest_vessel";
        # $detector{"mother"}      = "mutest_pipe_air1";
        $detector{"description"} = "Vessel of BDX-Hodo";
        $detector{"color"}       = "A05070";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "$rotX*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
        $detector{"material"}    = "G4_Fe";
        print_det(\%configuration, \%detector);
 
    
    
        $detector{"name"}        = "mutest_vessel_air";
        $detector{"mother"}      = "mutest_vessel";
        $detector{"description"} = "air in BDX-Hodo vessel ";
        $detector{"color"}       = "A05070";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Tube";
        $detector{"pos"}         = "0*cm 0*cm 0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "0*cm $par1*cm $par3*cm $par4*deg $par5*deg";
        $detector{"material"}    = "G4_AIR";
        print_det(\%configuration, \%detector);
   

    $Z = $vs_lg/2-$vs_top_tk/2 ;
    $par2 = $vs_ir/2; #InRad
    $par3  =$vs_top_tk/2;#length
    $par4 = 0.;
    $par5 = 360.;

    
    $detector{"name"}        = "mutest_vessel_top";
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    $Z = -($vs_lg/2-$vs_top_tk/2) ;

    $detector{"name"}        = "mutest_vessel_bottom";
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    
    
    %detector = init_det();
    # cystal
    $fg=0;
    $nblock=1;
    $ncol=1;
    $nrow=1;
    $al_ar_wr_cr_lgt=(29.5+5)/2.;
    
    
    my $rot=0;                # Carbon/Aluminum alveols
    %detector = init_det();
    $detector{"name"}        = "cry_alveol";
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "Carbon/Al container";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0;
    $Y = 0;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg $rot*deg 90*deg";
    $par1 =$al_ar_wr_cr_ssx ;
    $par2 =$al_ar_wr_cr_lsx;
    $par3 =$al_ar_wr_cr_ssy  ;
    $par4 =$al_ar_wr_cr_lsy  ;
    $par5 =$al_ar_wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "G4_AIR";
    #print_det(\%configuration, \%detector);# decomment if you want external box
    #print  " cryst",$Z," ","\n";
    
    # Air layer
    %detector = init_det();
    $detector{"name"}        = "cry_air";
    $detector{"mother"}      = "mutest_vessel_air";
    # $detector{"mother"}      = "cry_alveol"; # decomment if you want external box
    $detector{"description"} = "Air ";
    $detector{"color"}       = "00fff1";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$ar_wr_cr_ssx ;
    $par2 =$ar_wr_cr_lsx;
    $par3 =$ar_wr_cr_ssy  ;
    $par4 =$ar_wr_cr_lsy  ;
    $par5 =$ar_wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
    # Mylar wrapping
    %detector = init_det();
    $detector{"name"}        = "cry_mylar";
    $detector{"mother"}      = "cry_air";
    #$detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Mylar wrapping";
    $detector{"color"}       = "00fff2";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$wr_cr_ssx ;
    $par2 =$wr_cr_lsx;
    $par3 =$wr_cr_ssy  ;
    $par4 =$wr_cr_lsy  ;
    $par5 =$wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "bdx_mylar";
    print_det(\%configuration, \%detector);
    
    # Crystals
    %detector = init_det();
    $detector{"name"}        = "crystal in BDX-Hodo";
    $detector{"mother"}      = "cry_mylar";
    #$detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Crystal";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$cr_ssx ;
    $par2 =$cr_lsx;
    $par3 =$cr_ssy  ;
    $par4 =$cr_lsy  ;
    $par5 =$cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "CsI_Tl";
    #$detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 200 xch manual 0 ych manual 0";
    print_det(\%configuration, \%detector);
    
    
    
    # Scintillators
    # Ref system has x == Z, y == X z ==Y
    my $hodo_sc_thk = 1.9;
    $detector{"name"}        = "bdx-hodo-top";# top
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "bdx-hodo paddle on top";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $csi_pad_lx =11.0/2.;#short side (readout)
    my $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    my $csi_pad_lz =11.0/2 ;#long side (readout)
    $X = 0.;
    $Y = 0;
    $Z = 31.6/2+2.9+$hodo_sc_thk/2.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 13";
    print_det(\%configuration, \%detector);
    $Z = -(31.6/2+1.5+$hodo_sc_thk/2.);
    $detector{"name"}        = "bdx-hodo-bot";# bottom
    $detector{"description"} = "bdx-hodo paddle on bottom";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 12";
    print_det(\%configuration, \%detector);
    
    $hodo_sc_thk = .95;
    $detector{"name"}        = "bdx-hodo-left";# sideL
    $detector{"description"} = "bdx-hodo paddle on left side";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =7.6/2.;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =31.4/2 ;#long side (readout)
    $X = 4.7+0.5-$hodo_sc_thk/2.;
    $Y = 0;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 90*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 1";
    print_det(\%configuration, \%detector);
    $X = -(4.7+0.5-$hodo_sc_thk/2.);
    $detector{"name"}        = "bdx-hodo-right";# Side R
    $detector{"description"} = "bdx-hodo paddle on right side";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 2";
    print_det(\%configuration, \%detector);
    
    
    
    $detector{"name"}        = "bdx-hodo-front-ffl";# front-front-large (ffl)
    $detector{"description"} = "bdx-hodo paddle on front front large";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =8.0/2.;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =19.2/2 ;#long side (readout)
    $X = $csi_pad_lx-2.7;
    $Y = (4.1+0.9)+$hodo_sc_thk/2.;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 6";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "bdx-hodo-front-fbl";# front-back-large (fbl)
    $detector{"description"} = "bdx-hodo paddle on front back large";
    $X = -($csi_pad_lx-2.7);
    $Y = (4.1+0.8)-$hodo_sc_thk/2.;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 5";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "bdx-hodo-front-ffs";#  front-front-small (ffs)
    $detector{"description"} = "bdx-hodo paddle on front front small";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =2.5/2.;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =19.2/2 ;#long side (readout)
    $X = -($csi_pad_lx+2.8);
    $Y = (4.1+0.9)+$hodo_sc_thk/2.;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 7";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "bdx-hodo-front-fbs";# front-back-small (fbs)
    $detector{"description"} = "bdx-hodo paddle on front back small";
    $X =($csi_pad_lx+2.8);
    $Y = (4.1+0.8)-$hodo_sc_thk/2.;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 4";
    print_det(\%configuration, \%detector);

    
    
    
    $detector{"name"}        = "bdx-hodo-back-large";# back large (bl)
    $detector{"description"} = "bdx-hodo paddle on back large";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =14.4/2;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =19.2/2 ;#long side (readout)
    $X = 0.;
    $Y = -(4.0+$hodo_sc_thk/2.);
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 3";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "bdx-hodo-back-bfl";# back-front-large (bfl)
    $detector{"description"} = "bdx-hodo paddle on back front large";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =5.0/2.;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =10.6/2 ;#long side (readout)
    $X = 0.;
    $Y = -(4.0+1.6*$hodo_sc_thk);
    $Z = -($csi_pad_lx-1.2);
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "180*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 11";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "bdx-hodo-back-bbl";#  back-back-large bbl)
    $detector{"description"} = "bdx-hodo paddle on back back large";
    $X = 0;
    $Y = -(4.0+2.7*$hodo_sc_thk);
    $Z = ($csi_pad_lx-1.2);
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 8";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "bdx-hodo-front-bfs";# back-front-small (bfs)
    $detector{"description"} = "bdx-hodo paddle on back front small";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =2.5/2.;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =10.6/2 ;#long side (readout)
    $X = 0.;
    $Y = -(4.0+1.6*$hodo_sc_thk);
    $Z = ($csi_pad_lx+1.3);
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 10";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "bdx-hodo-front-bbs";# back-back-small (bbs)
    $detector{"description"} = "bdx-hodo paddle on back back small";
    $X = 0.;
    $Y = -(4.0+2.7*$hodo_sc_thk);
    $Z = -($csi_pad_lx+1.3);
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 9";
    print_det(\%configuration, \%detector);


}

sub make_mutest_minibdx
{
    print  " &&&&&&&&&  BDX-MINI PARAMTERS &&&&&&&&","\n";

    my %detector = init_det();
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = 0.; # depth in the pipe
    my $rotX=0.;
    my $rotY=0.;
    #nominal:
    #my $rotZ=180.;

    my $rotZ=160.;
    
    if (($configuration{"variation"} eq "CT") or ($flag_detcentered eq "1"))
    {
        $detector{"mother"}      = "main_volume";
        $rotX=90.; #0 for H 90deg for V
        $Z = 0.
    }
    else
    {
	$detector{"mother"}      = "mutest_pipe_air";
    }
    
    
    my $vs_top_tk=0.5;
    my $vs_side_tk=1.0;
    my $vs_lg=89.;
    my $vs_ir=20.;
    my $par1 = $vs_ir/2; #InRad
    my $par2 = ($vs_ir+$vs_side_tk)/2; #InRad+thick
    
    my $par3  =$vs_lg/2.;#length
    my $par4 = 0.;    my $par5 = 360.;
    
    print  " BDX-MINI position (cm) X=",$X,"\n";
    print  " BDX-MINI position (cm) Y=",$Y,"\n";
    print  " BDX-MINI position (cm) Z=",$Z,"\n";
    $detector{"name"}        = "mutest_vessel_mother";
    $detector{"description"} = "Mother volume for the vessel";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$rotX*deg $rotY*deg $rotZ*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    $X=0.;
    $Y=0.;
    $Z=0.;
    
    $detector{"name"}        = "mutest_vessel";
    $detector{"mother"}      = "mutest_vessel_mother";
    $detector{"description"} = "Vessel of BDX-Hodo";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    

    print  " Vessel INNER radius = ",$par1,"\n";
    print  " Vessel OUTER radius = ",$par2,"\n";

    $detector{"name"}        = "mutest_vessel_air";
    $detector{"mother"}      = "mutest_vessel_mother";
    $detector{"description"} = "air in BDX-Hodo vessel ";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par1*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    

    $Z = $vs_lg/2-$vs_top_tk/2 ;
    $par2 = $vs_ir/2; #InRad
    $par3  =$vs_top_tk/2;#length
    $par4 = 0.;
    $par5 = 360.;
    
    
    $detector{"name"}        = "mutest_vessel_top";
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    $Z= -($vs_lg/2-$vs_top_tk/2) ;

    
    $detector{"name"}        = "mutest_vessel_bottom";
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    


    my $ov_top_tk=0.8;
    my $ov_side_tk=0.8*2;
    my $ov_lg=53.;
    my $ov_ir=17.7;
    my $ov_or= ($ov_ir+$ov_side_tk); #InRad+thick
    $par1 = $ov_ir/2; #InRad
    $par2 = $ov_or/2; #InRad
    $par3  =$ov_lg/2.;#length
    $par4 = 0.;
    $par5 = 360.;
    $Z=0;
    ###
    
    my $NchOV = 1; # MAX 8 otherwise overlaps with OV_TOP/BOTTOM (8,9)
    my $DeltaAng = 360./$NchOV;
    my $DeltaToll = 0.;#Tolerance (in cm) between adiacent sectors
    my $DeltaAngShift = 0.; #fixed shift
    $DeltaAngShift = 90. -$DeltaAng/2; #to have #1 perfectly back
    ###
    for(my $ib=1; $ib<($NchOV+1); $ib++)
    #for(my $ib=1; $ib<(2); $ib++)
    {

        $par4=($ib-1)*$DeltaAng+$DeltaAngShift+$DeltaToll/2;
        $par5=$DeltaAng-$DeltaToll;# -0.1deg tolerance
    $detector{"name"}        = "mutest_OV_$ib";
    $detector{"description"} = "OV in the pipe";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 7 channel manual $ib";
    print_det(\%configuration, \%detector);
        #print  " OV segments = ",$par4,"-",$par5,"\n";
    }
    print  " OV INNER radius = ",$par1,"\n";
    print  " OV OUTER radius = ",$par2,"\n";
    print  " OV segmentations Nch = ",$NchOV,"\n";
    
    $detector{"name"}        = "mutest_OV_air";
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "air in BDX-Mini OV ";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par1*cm $par3*cm 0.*deg 360.*deg";
    $detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "no";

    print_det(\%configuration, \%detector);
    
    $Z = $ov_lg/2+$ov_top_tk/2 ;
    $par2 = $ov_or/2; #InRad
    $par3  =$ov_top_tk/2;#thickness
    $par4 = 0.;
    $par5 = 360.;


    $detector{"name"}        = "mutest_OV_top";
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "cc6804";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 7 channel manual 9";
    print_det(\%configuration, \%detector);
    
    $Z = -($ov_lg/2-$ov_top_tk/2) ;
    $par2 = $ov_ir/2; #InRad
    
    $detector{"name"}        = "mutest_OV_bottom";
    $detector{"mother"}      = "mutest_OV_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "cc6804";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 7 channel manual 10";
   print_det(\%configuration, \%detector);

    my $iv_top_tk=0.8;
    my $iv_side_tk=0.8*2;
    my $iv_lg=49.4;
    my $iv_ir=(16.6+0.8-$iv_side_tk);
    my $iv_or=($iv_ir+$iv_side_tk);

    
    $par1 = $iv_ir/2; #InRad
    $par2 = $iv_or/2; #InRad+thick
    my $par1a = $par1 * 0.92388;
    my $par2a = $par2 * 0.92388;

    $par3  =$iv_lg/2.;#length
    my $par3a = -$iv_lg/2.;
    $par4 = 0.;
    $par5 = 360.;
    $Z=0;
    ###
    
    my $NchIV = 1; #MAX 8 otherwise overlasps with IV TOP/Bottom (9/10)
    $DeltaAng = 360./$NchIV;
    $DeltaToll = 0.; #Tolerance (in cm) between adiacent sectors
    #$DeltaAngShift = 0.; #fixed shift
    $DeltaAngShift = 90. -$DeltaAng/2; #to have #1 perfectly back
    ###
    for(my $ib=1; $ib<($NchIV+1); $ib++)
    #for(my $ib=1; $ib<(2); $ib++)
    {
        
        $par4=($ib-1)*$DeltaAng+$DeltaAngShift+$DeltaToll/2;
        $par5=$DeltaAng-$DeltaToll;# -0.1deg tolerance

    
	$detector{"name"}        = "mutest_IV_$ib";
	$detector{"mother"}      = "mutest_OV_air";
	$detector{"description"} = "IV in the pipe";
	$detector{"color"}       = "ffff02";
	$detector{"style"}       = 0;
	$detector{"visible"}     = 1;

#	$detector{"type"}        = "Tube";
#	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
#	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
#	$detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";

	$detector{"type"}        = "Pgon";
	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	$detector{"rotation"}    = "0*deg 0*deg 22.5*deg";
	$detector{"dimensions"}  = "0*deg 360*deg 8 2 $par1a*cm $par2a*cm $par3a*cm $par1a*cm $par2a*cm $par3*cm";

	$detector{"material"}    = "ScintillatorB";
	$detector{"sensitivity"} = "veto";
	$detector{"hit_type"}    = "veto";
	$detector{"identifiers"} = "sector manual 0 veto manual 8 channel manual $ib";
	print_det(\%configuration, \%detector);
    }
    print  " IV INNER radius = ",$par1,"\n";
    print  " IV OUTER radius = ",$par2,"\n";
    print  " IV segmentations Nch = ",$NchIV,"\n";
    

#    $detector{"name"}        = "mutest_IV_air";
#    $detector{"mother"}      = "mutest_OV_air";
#    $detector{"description"} = "air in BDX-Mini IV ";
#    $detector{"color"}       = "ffff02";
#    $detector{"style"}       = 0;
#    $detector{"visible"}     = 1;
##    $detector{"type"}        = "Tube";
##    $detector{"pos"}         = "0*cm 0*cm 0*cm";
##    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
##    $detector{"dimensions"}  = "0*cm $par1*cm $par3*cm 0.*deg 360.*deg";

#    $detector{"type"}        = "Pgon";
#    $detector{"pos"}         = "0*cm 0*cm 0*cm";
#    $detector{"rotation"}    = "0*deg 0*deg 22.5*deg";
#    $detector{"dimensions"}  = "0*deg 360*deg 8 2 0*cm $par1a*cm $par3a*cm 0*cm $par1a*cm $par3*cm";

#    $detector{"material"}    = "G4_AIR";
#    $detector{"sensitivity"} = "no";
#    print_det(\%configuration, \%detector);

    $Z = $iv_lg/2+$iv_top_tk/2 ;
    $par2 = $iv_or/2; #InRad
    $par3  =$iv_top_tk/2;#thickness
    $par2a = $par2 * 0.92388;
    $par3a = -$par3;

    $par4 = 0.;
    $par5 = 360.;
    
    $detector{"name"}        = "mutest_IV_top";
    $detector{"mother"}      = "mutest_OV_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "c4c401";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    #$detector{"type"}        = "Tube";
    #$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    #$detector{"rotation"}    = "0*deg 0*deg 0*deg";
    #$detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";

    $detector{"type"}        = "Pgon";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 22.5*deg";
    $detector{"dimensions"}  = "0*deg 360*deg 8 2 0*cm $par2a*cm $par3a*cm 0*cm $par2a*cm $par3*cm";
    

    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 8 channel manual 9";
    print_det(\%configuration, \%detector);
    $Z = -($iv_lg/2-$iv_top_tk/2) ;
    $par2 = $iv_ir/2; #InRad
    $par2a = $par2 * 0.92388;
 
    $detector{"name"}        = "mutest_IV_bottom";
    $detector{"mother"}      = "mutest_OV_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "c4c401";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;

#    $detector{"type"}        = "Tube";
#    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
#    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
#    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";

    $detector{"type"}        = "Pgon";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg -22.5*deg";
    $detector{"dimensions"}  = "0*deg 360*deg 8 2 0*cm $par2a*cm $par3a*cm 0*cm $par2a*cm $par3*cm";

    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 8 channel manual 10";
   if (($JLab_paddle_setup) eq "0")
   {print_det(\%configuration, \%detector);}



    
    my $w_top_tk=0.8;
    my $w_side_tk=0.8*2;
    my $w_lg=45.0;
    my $w_ir=(14.8+0.8-$w_side_tk);
     my $w_or=($w_ir+$w_side_tk);
    $par1 = $w_ir/2; #InRad
    $par2 = $w_or/2; #InRad+thick
    
    $par1a = $par1*0.92388;
    $par2a = $par2*0.92388;

    $par3  =$w_lg/2.;#thickness
    $par3a =-$par3;
    $par4 = 0.;
    $par5 = 360.;
        $Z=0;
    ###
    $detector{"name"}        = "mutest_W_vault";
    $detector{"mother"}      = "mutest_OV_air";
    $detector{"description"} = "W vault in the pipe";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    
#    $detector{"type"}        = "Tube";
#    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
#    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
#    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
 
    $detector{"type"}        = "Pgon";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg -22.5*deg";
    $detector{"dimensions"}  = "0*deg 360*deg 8 2 $par1a*cm $par2a*cm $par3a*cm $par1a*cm $par2a*cm $par3*cm";

    $detector{"material"}    = "G4_W";
    $detector{"sensitivity"} = "no";
    print_det(\%configuration, \%detector);
    
    print  " W shield radius = ",$par1,"\n";
    print  " W shield radius = ",$par2,"\n";
    
#    $detector{"name"}        = "mutest_W_air";
#    $detector{"mother"}      = "mutest_IV_air";
#    $detector{"description"} = "air in BDX-Mini W ";
#    $detector{"color"}       = "A1D8F0";
#    $detector{"style"}       = 0;
#    $detector{"visible"}     = 1;

##    $detector{"type"}        = "Tube";
##    $detector{"pos"}         = "0*cm 0*cm 0*cm";
##    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
##    $detector{"dimensions"}  = "0*cm $par1*cm $par3*cm $par4*deg $par5*deg";

#   $detector{"type"}        = "Pgon";
#    $detector{"pos"}         = "0*cm 0*cm 0*cm";
#    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
#    $detector{"dimensions"}  = "0*deg 360*deg 8 2 0*cm $par1a*cm $par3a*cm 0*cm $par1a*cm $par3*cm";
#    $detector{"material"}    = "G4_AIR";
#    print_det(\%configuration, \%detector);

    $Z = $w_lg/2+$w_top_tk/2 ;
    $par2 = $w_or/2; #InRad
    $par2a = $par2 * 0.92388;
    $par3  =$w_top_tk/2;#length
    $par3a = -$par3;
    $par4 = 0.;
    $par5 = 360.;
    
    $detector{"name"}        = "mutest_W_top";
    $detector{"mother"}      = "mutest_OV_air";
    $detector{"description"} = "top W ";
    $detector{"color"}       = "85a4c1";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
#    $detector{"type"}        = "Tube";
#    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
#    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
#    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";

    $detector{"type"}        = "Pgon";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg -22.5*deg";
    $detector{"dimensions"}  = "0*deg 360*deg 8 2 0*cm $par2a*cm $par3a*cm 0*cm $par2a*cm $par3*cm";

    $detector{"material"}    = "G4_W";


    print_det(\%configuration, \%detector);
    $Z = -($w_lg/2-$w_top_tk/2) ;
    $par2 = $w_ir/2; #InRad
    $par2a = $par2 * 0.92388;
    $detector{"name"}        = "mutest_W_bottom";
    $detector{"mother"}      = "mutest_OV_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "85a4c1";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;

#    $detector{"type"}        = "Tube";
#    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
#    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
#    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";

    $detector{"type"}        = "Pgon";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg -22.5*deg";
    $detector{"dimensions"}  = "0*deg 360*deg 8 2 0*cm $par2a*cm $par3a*cm 0*cm $par2a*cm $par3*cm";

    $detector{"material"}    = "G4_W";
    print_det(\%configuration, \%detector);
    ####
    if (($JLab_paddle_setup) eq "1")
    {
     ### SApecial run at Jlab in March using small paddle
    $detector{"name"}        = "bdx-hodo-front-ffs";#  front-front-small (ffs)
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "bdx-hodo paddle on front front small";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $csi_pad_lx =2.5/2.;#short side (readout)
    my $csi_pad_ly = .95/2.; #Thikness
    my $csi_pad_lz =19.2/2 ;#long side (readout)
    $X =0.;
    $Y= ($vs_lg/2.)-54.5-$csi_pad_lz/2-1.2;
    $Z = -($vs_ir/2+$vs_side_tk);
    print $Z,"\n";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 8 channel manual 10";
    print_det(\%configuration, \%detector);
    }
    ####

    
    # START PbWO4e crystals 400 1 15 TOP
    my $cr_Panda_Lp=2.0;;
    my $cr_Panda_lgt=20;
    my $cr_FT_Lf=1.5;;
    my $cr_FT_lgt=20;
    my $cr_airgap=0.05;
    my $cr_airgapz=0.2;
    my $Delta=0.;
    my $nx=3;
    my $Delta_pos=($cr_airgap+$cr_Panda_Lp);
    my $Delta_pos1=($cr_airgap+0.5*$cr_Panda_Lp);
    
    
    # START PANDA crystals
    # icc  : LEFT=icc=-1 Right=icc=+1 (along the beam)
    # il  :  Upstream= negative ; Downstream = positive (following the beam direction)
    # irr : Top=irr=1 BOTTOM=irr=-1 (below/above the beamn height)
    
    for(my $ir=0; $ir<2; $ir++)
    {
     for(my $is=0; $is<2; $is++) #
       {
        my $irr=2*$ir-1;
        my $iss=2*$is-1;
        my $irot=180*$ir;
        for(my $il=0; $il<4; $il++) # PANDA Lines 0 to +3
        {
            if (($il) == 3) {$nx=2}# only 1 crystal in the first line
            else {$nx=3}
            for(my $ic=1; $ic<$nx; $ic++) # PANDA column 1-2
            {
                my $icc=$ic*$iss;
                
                if (($il) == 0 || ($il) == 1)
                {$Delta=$Delta_pos1;
                    $X=$iss*($cr_FT_Lf+$Delta+($ic-1)*$Delta_pos);
                }# using FT shift for rows 0 and 1
                else
                {$Delta=$Delta_pos;
                    $X=$iss*(($ic-0.5)*$Delta_pos);
                }
                


         $Y=(-0.5+$il)*$Delta_pos;
         $Z = ($cr_Panda_lgt+$cr_airgapz)/2.*$irr;
         $detector{"name"}        = "CrsP"."_"."$icc"."_"."$il"."_"."$irr";
         $detector{"mother"}      = "mutest_OV_air";
            $detector{"description"} = "Panda PbWO4 crystal"."$icc"."$il"."$irr";
         $detector{"color"}       = "1c86ea";
         $detector{"style"}       = 1;
         $detector{"visible"}     = 1;
         $detector{"type"}        = "Trd";
         #print 'Ind all: ',$IndxAss,' ',$irr,"\n";
                #print 'Ind all: ',$Delta,"\n";
         $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
         $detector{"rotation"}    = "0*deg $irot*deg 0*deg";
         $par1 =$cr_Panda_Lp/2;
         $par2 =$cr_Panda_Lp/2;
         $par3 =$cr_Panda_Lp/2 ;
         $par4 =$cr_Panda_Lp/2 ;
         $par5 =$cr_Panda_lgt/2 ;
         $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
         $detector{"material"}    = "G4_PbWO4";
         $detector{"sensitivity"} = "crs";
         $detector{"hit_type"}    = "crs";
                my $man=400+$irr+1;
         $detector{"identifiers"} = "sector manual $man xch manual $icc ych manual $il";
                print_det(\%configuration, \%detector);
            }
        }
       }
    }
    # END PANDA crystals
    # START FT crystals
    for(my $ir=0; $ir<2; $ir++) # Top=irr=1 BOTTOM=irr=-1)
    {
            my $irr=2*$ir-1;
            my $irot=180*$ir;
            my $orn=0;
        my $suba=0;
        my $icc=0.;
            for(my $il=-2; $il<2; $il++) # FT  Lines -2 to +1
            {
                if ($il eq -2)
                {$nx=5;
                    $orn=0.;
                 $suba=2;
                }
                else
                {$nx=1;
                    $orn=90.;
                $suba=0;
                }
              
                for(my $is=0; $is<$nx; $is++)
                {
                    $icc=$is-$suba;
                    #print 'Ind : ',$icc,"\n";
                    $X=$icc*($cr_FT_Lf+$cr_airgap);
                    $Y=($il-1.5)*($cr_FT_Lf+$cr_airgap)+$cr_Panda_Lp;
                    if ($il == -2)
                    {
                        $Y=-3.5*$cr_airgap-4*$cr_FT_Lf+$cr_Panda_Lp;
                        if ($is == 0 || $is == 4) {$Y=-1.5*$cr_airgap-1*$cr_FT_Lf-$cr_Panda_Lp};
                    }
                    
                    #print 'Y: ',$Y,"\n";
                    $Z = ($cr_Panda_lgt+$cr_airgapz)/2.*$irr;
                    $detector{"name"}        = "CrsFT"."_"."$icc"."_"."$il"."_"."$irr";
                    $detector{"mother"}      = "mutest_OV_air";
                    $detector{"description"} = "FT PbWO4 crystal"."$icc"."$il"."$irr";
                    $detector{"color"}       = "2d50ff";
                    $detector{"style"}       = 1;
                    $detector{"visible"}     = 1;
                    $detector{"type"}        = "Trd";
                    #print 'Ind all: ',$IndxAss,' ',$irr,"\n";
                    #print 'Ind all: ',$Delta,"\n";
                    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
		    $detector{"rotation"}    = "0*deg $irot*deg $orn*deg";
                    $par1 =$cr_FT_Lf/2;
                    $par2 =$cr_FT_Lf/2;
                    $par3 =$cr_FT_Lf ;
                    $par4 =$cr_FT_Lf ;
                    $par5 =$cr_FT_lgt/2 ;
                    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
                    $detector{"material"}    = "G4_PbWO4";
                    $detector{"sensitivity"} = "crs";
                    $detector{"hit_type"}    = "crs";
                     my $man=500+$irr+1;
                    $detector{"identifiers"} = "sector manual $man xch manual $icc ych manual $il";
                    print_det(\%configuration, \%detector);
                }
            }
        }
    

    
    # END PbWO4e crystals
=pod

    ### CsI Crystals
    ###  TR (300 0 0) TL (300 0 1)
    ###  BR (300 0 2) BL (300 0 3)
    ##############
    %detector = init_det();
    # cystal
    $fg=0;
    $nblock=1;
    $ncol=1;
    $nrow=1;
    $al_ar_wr_cr_lgt=(29.5+5)/2.;
    my $z_shift=$al_ar_wr_cr_lgt-1;
    #$z_shift=0;
    
    #CRYSTAL 1
    # Air layer
    my $rot=180;
    %detector = init_det();
    $detector{"name"}        = "cry_air_TR";
    $detector{"mother"}      = "mutest_W_air";
    # $detector{"mother"}      = "cry_alveol"; # decomment if you want external box
    $detector{"description"} = "Air ";
    $detector{"color"}       = "00fff1";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = -3.075;
    $Y = 0.;
    $Z = 0.+$z_shift;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg $rot*deg 0*deg";
    $par1 =$ar_wr_cr_ssx ;
    $par2 =$ar_wr_cr_lsx;
    $par3 =$ar_wr_cr_ssy  ;
    $par4 =$ar_wr_cr_lsy  ;
    $par5 =$ar_wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
    # Mylar wrapping
    %detector = init_det();
    $detector{"name"}        = "cry_mylar_TR";
    $detector{"mother"}      = "cry_air_TR";
    #$detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Mylar wrapping";
    $detector{"color"}       = "00fff2";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$wr_cr_ssx ;
    $par2 =$wr_cr_lsx;
    $par3 =$wr_cr_ssy  ;
    $par4 =$wr_cr_lsy  ;
    $par5 =$wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "bdx_mylar";
    print_det(\%configuration, \%detector);
    

    %detector = init_det();
    $detector{"name"}        = "crystal_TR";
    $detector{"mother"}      = "cry_mylar_TR";
    #$detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Crystal TR";
    $detector{"color"}       = "dd0404";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$cr_ssx ;
    $par2 =$cr_lsx;
    $par3 =$cr_ssy  ;
    $par4 =$cr_lsy  ;
    $par5 =$cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "CsI_Tl";
    #$detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 300 xch manual 2 ych manual 1";
    print_det(\%configuration, \%detector);
    
    #CRYSTAL 2
    
    # Air layer
    $rot=0;
    %detector = init_det();
    $detector{"name"}        = "cry_air_TL";
    $detector{"mother"}      = "mutest_W_air";
    # $detector{"mother"}      = "cry_alveol"; # decomment if you want external box
    $detector{"description"} = "Air ";
    $detector{"color"}       = "00fff1";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = +3.075;
    $Y = 0;
    $Z = 0.+$z_shift;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg $rot*deg 0*deg";
    $par1 =$ar_wr_cr_ssx ;
    $par2 =$ar_wr_cr_lsx;
    $par3 =$ar_wr_cr_ssy  ;
    $par4 =$ar_wr_cr_lsy  ;
    $par5 =$ar_wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
    # Mylar wrapping
    %detector = init_det();
    $detector{"name"}        = "cry_mylar_TL";
    $detector{"mother"}      = "cry_air_TL";
    #$detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Mylar wrapping";
    $detector{"color"}       = "00fff2";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$wr_cr_ssx ;
    $par2 =$wr_cr_lsx;
    $par3 =$wr_cr_ssy  ;
    $par4 =$wr_cr_lsy  ;
    $par5 =$wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "bdx_mylar";
    print_det(\%configuration, \%detector);
    
    # Crystals
    %detector = init_det();
    $detector{"name"}        = "crystal_TL";
    $detector{"mother"}      = "cry_mylar_TL";
    #$detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Crystal TL";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$cr_ssx ;
    $par2 =$cr_lsx;
    $par3 =$cr_ssy  ;
    $par4 =$cr_lsy  ;
    $par5 =$cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "CsI_Tl";
    #$detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 300 xch manual 1 ych manual 1";
    print_det(\%configuration, \%detector);

    #CRYSTAL 3
    # Air layer
    $rot=180;
    %detector = init_det();
    $detector{"name"}        = "cry_air_BR";
    $detector{"mother"}      = "mutest_W_air";
    # $detector{"mother"}      = "cry_alveol"; # decomment if you want external box
    $detector{"description"} = "Air ";
    $detector{"color"}       = "00fff1";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = -3.075;
    $Y = 0.;
    $Z = 0.-$z_shift;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg $rot*deg 0*deg";
    $par1 =$ar_wr_cr_ssx ;
    $par2 =$ar_wr_cr_lsx;
    $par3 =$ar_wr_cr_ssy  ;
    $par4 =$ar_wr_cr_lsy  ;
    $par5 =$ar_wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
    # Mylar wrapping
    %detector = init_det();
    $detector{"name"}        = "cry_mylar_BR";
    $detector{"mother"}      = "cry_air_BR";
    #$detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Mylar wrapping";
    $detector{"color"}       = "00fff2";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$wr_cr_ssx ;
    $par2 =$wr_cr_lsx;
    $par3 =$wr_cr_ssy  ;
    $par4 =$wr_cr_lsy  ;
    $par5 =$wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "bdx_mylar";
    print_det(\%configuration, \%detector);
    
    
    %detector = init_det();
    $detector{"name"}        = "crystal_BR";
    $detector{"mother"}      = "cry_mylar_BR";
    #$detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Crystal BR";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$cr_ssx ;
    $par2 =$cr_lsx;
    $par3 =$cr_ssy  ;
    $par4 =$cr_lsy  ;
    $par5 =$cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "CsI_Tl";
    #$detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 300 xch manual 2 ych manual -1";
    print_det(\%configuration, \%detector);
    
    #CRYSTAL 4
    
    # Air layer
     $rot=0;
    %detector = init_det();
    $detector{"name"}        = "cry_air_BL";
    $detector{"mother"}      = "mutest_W_air";
    # $detector{"mother"}      = "cry_alveol"; # decomment if you want external box
    $detector{"description"} = "Air ";
    $detector{"color"}       = "00fff1";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = +3.075;
    $Y = 0;
    $Z = 0.-$z_shift;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg $rot*deg 0*deg";
    $par1 =$ar_wr_cr_ssx ;
    $par2 =$ar_wr_cr_lsx;
    $par3 =$ar_wr_cr_ssy  ;
    $par4 =$ar_wr_cr_lsy  ;
    $par5 =$ar_wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
    # Mylar wrapping
    %detector = init_det();
    $detector{"name"}        = "cry_mylar_BL";
    $detector{"mother"}      = "cry_air_BL";
    #$detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Mylar wrapping";
    $detector{"color"}       = "00fff2";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$wr_cr_ssx ;
    $par2 =$wr_cr_lsx;
    $par3 =$wr_cr_ssy  ;
    $par4 =$wr_cr_lsy  ;
    $par5 =$wr_cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "bdx_mylar";
    print_det(\%configuration, \%detector);
    
    # Crystals
    %detector = init_det();
    $detector{"name"}        = "crystal_BL";
    $detector{"mother"}      = "cry_mylar_BL";
    #$detector{"mother"}      = "Det_house_inner";
    $detector{"description"} = "Crystal BL";
    $detector{"color"}       = "00ffff";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Trd";
    $X = 0.;
    $Y = 0.;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $par1 =$cr_ssx ;
    $par2 =$cr_lsx;
    $par3 =$cr_ssy  ;
    $par4 =$cr_lsy  ;
    $par5 =$cr_lgt ;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*cm $par5*cm";
    $detector{"material"}    = "CsI_Tl";
    #$detector{"material"}    = "G4_AIR";
    $detector{"sensitivity"} = "crs";
    $detector{"hit_type"}    = "crs";
    $detector{"identifiers"} = "sector manual 300 xch manual 1 ych manual -1";
    print_det(\%configuration, \%detector);

   

    
    # Scintillators
    # Ref system has x == Z, y == X z ==Y
    my $hodo_sc_thk = 1.9;
    $detector{"name"}        = "bdx-hodo-top";# top
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "bdx-hodo paddle on top";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $csi_pad_lx =11.0/2.;#short side (readout)
    my $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    my $csi_pad_lz =11.0/2 ;#long side (readout)
    $X = 0.;
    $Y = 0;
    $Z = 31.6/2+2.9+$hodo_sc_thk/2.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "90*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 13";
    print_det(\%configuration, \%detector);
    $Z = -(31.6/2+1.5+$hodo_sc_thk/2.);
    $detector{"name"}        = "bdx-hodo-bot";# bottom
    $detector{"description"} = "bdx-hodo paddle on bottom";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 12";
    print_det(\%configuration, \%detector);
    
    $hodo_sc_thk = .95;
    $detector{"name"}        = "bdx-hodo-left";# sideL
    $detector{"description"} = "bdx-hodo paddle on left side";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =7.6/2.;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =31.4/2 ;#long side (readout)
    $X = 4.7+0.5-$hodo_sc_thk/2.;
    $Y = 0;
    $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 90*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 1";
    print_det(\%configuration, \%detector);
    $X = -(4.7+0.5-$hodo_sc_thk/2.);
    $detector{"name"}        = "bdx-hodo-right";# Side R
    $detector{"description"} = "bdx-hodo paddle on right side";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 2";
    print_det(\%configuration, \%detector);
    
    
    
    $detector{"name"}        = "bdx-hodo-front-ffl";# front-front-large (ffl)
    $detector{"description"} = "bdx-hodo paddle on front front large";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =8.0/2.;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =19.2/2 ;#long side (readout)
    $X = $csi_pad_lx-2.7;
    $Y = (4.1+0.9)+$hodo_sc_thk/2.;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 6";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "bdx-hodo-front-fbl";# front-back-large (fbl)
    $detector{"description"} = "bdx-hodo paddle on front back large";
    $X = -($csi_pad_lx-2.7);
    $Y = (4.1+0.8)-$hodo_sc_thk/2.;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 5";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "bdx-hodo-front-ffs";#  front-front-small (ffs)
    $detector{"description"} = "bdx-hodo paddle on front front small";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =2.5/2.;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =19.2/2 ;#long side (readout)
    $X = -($csi_pad_lx+2.8);
    $Y = (4.1+0.9)+$hodo_sc_thk/2.;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 7";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "bdx-hodo-front-fbs";# front-back-small (fbs)
    $detector{"description"} = "bdx-hodo paddle on front back small";
    $X =($csi_pad_lx+2.8);
    $Y = (4.1+0.8)-$hodo_sc_thk/2.;
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 4";
    print_det(\%configuration, \%detector);
    
    
    
    
    $detector{"name"}        = "bdx-hodo-back-large";# back large (bl)
    $detector{"description"} = "bdx-hodo paddle on back large";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =14.4/2;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =19.2/2 ;#long side (readout)
    $X = 0.;
    $Y = -(4.0+$hodo_sc_thk/2.);
    $Z = 0;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 3";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "bdx-hodo-back-bfl";# back-front-large (bfl)
    $detector{"description"} = "bdx-hodo paddle on back front large";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =5.0/2.;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =10.6/2 ;#long side (readout)
    $X = 0.;
    $Y = -(4.0+1.6*$hodo_sc_thk);
    $Z = -($csi_pad_lx-1.2);
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "180*deg 90*deg 0*deg";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 11";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "bdx-hodo-back-bbl";#  back-back-large bbl)
    $detector{"description"} = "bdx-hodo paddle on back back large";
    $X = 0;
    $Y = -(4.0+2.7*$hodo_sc_thk);
    $Z = ($csi_pad_lx-1.2);
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 8";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "bdx-hodo-front-bfs";# back-front-small (bfs)
    $detector{"description"} = "bdx-hodo paddle on back front small";
    $detector{"color"}       = "ff8000";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $csi_pad_lx =2.5/2.;#short side (readout)
    $csi_pad_ly = $hodo_sc_thk/2.; #Thikness
    $csi_pad_lz =10.6/2 ;#long side (readout)
    $X = 0.;
    $Y = -(4.0+1.6*$hodo_sc_thk);
    $Z = ($csi_pad_lx+1.3);
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"dimensions"}  = "$csi_pad_lx*cm $csi_pad_ly*cm $csi_pad_lz*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 10";
    print_det(\%configuration, \%detector);
    $detector{"name"}        = "bdx-hodo-front-bbs";# back-back-small (bbs)
    $detector{"description"} = "bdx-hodo paddle on back back small";
    $X = 0.;
    $Y = -(4.0+2.7*$hodo_sc_thk);
    $Z = -($csi_pad_lx+1.3);
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"identifiers"} = "sector manual 0 veto manual 6 channel manual 9";
    print_det(\%configuration, \%detector);
=cut
    print  " &&&&&&&&&  BDX-MINI PARAMTERS END &&&&&&&&","\n";

}



#################################################################################################
#
# End: BDX-p veto
#
#################################################################################################
sub make_hallA_bdx
{
    if (($flag_detcentered) eq "0")
        {
            if(($flag_mutest) eq "1")
            {
                make_bdx_main_volume();
                make_dirt_u();
                make_dirt_d();
                make_dirt();
                make_bunker_main();
                make_bunker_tunnel();
                make_bunker();
                make_bunker_end();
                make_hallaBD();
                make_hallaBD_flux_barrel();
                make_hallaBD_flux_endcup();
                make_hallaBD_flux();
                make_hallaBD_flux_sample();
                make_dirt_top();
            }
            else
            {
                make_bdx_main_volume();
                make_dirt_u();
                make_dirt_d();
                make_dirt();
                make_bunker_main();
                make_bunker_tunnel();
                make_bunker();
                make_bunker_end();
                make_hallaBD();
                #make_hallaBD_flux_barrel();
                #make_hallaBD_flux_endcup();
                #make_hallaBD_flux();
                #make_hallaBD_flux_sample();
                make_muon_absorber();
                make_det_house_outer();
                make_det_house_inner();
                make_det_shaft_outer();
                make_det_shaft_inner();
                make_stair_outer();
                make_stair_inner();
                make_stair_wall();
                make_stair_wall_door();
                make_shaft_wall();
                make_shaft_wall_door();
                make_ext_house_outer();
                make_ext_house_inner();
                make_ext_house_shaft_hole();
                make_ext_house_stair_hole();
                make_stair_steps_1();
                make_stair_steps_2();
                make_stair_steps_3();
                #make_cry_module();
            }
        }
    else
        {
            make_bdx_main_volume();
        }
}


sub make_hallA_bdxMini
{
    make_mutest_pipe();
    make_mutest_minibdx();
    make_mutest_bdxMiniDirt();
    make_bdx_main_volume();
}



sub make_detector_bdx
 {
     if (($configuration{"variation"} eq "Proposal") and ($flag_JlabCT eq "1"))
     {
         make_bdx_main_volume();
         #make_cormo_flux();
         #make_cormo_det();
         make_cormo_iveto;
         make_cormo_oveto;
         make_cormo_lead();
         # make_csi_pad();
         # make_flux_cosmic_cyl;
         # make_flux_cosmic_rec;
         #make_cormo_shield();
         #  make_babar_crystal();
         #make_cry_module_up();
         #make_cry_module_down();
         make_cry_module();
         make_cal_pad();
         #make_sarc_lead();
         #  make_cry_module_II();1')
     }
     else
     {
         if(($flag_mutest) eq '1')
         {
             if(($flag_minibdx) eq '0')
              # Standard BDX-Hodo
             {make_mutest_detector();
             }
             else
             {make_mutest_minibdx();
                 #make_flux_cosmic_sph();
             }
             if (($flag_detcentered) eq '0')
             {
                 make_mutest_pipe();
                 make_mutest_pipe2();
                 #make_mutest_flux();
                 #make_flux_cosmic_cyl;
             }
         }
         else
         {
         #make_cormo_flux();
         #make_cormo_det();
         #make_cormo_iveto;
         #make_cormo_oveto;
         #make_cormo_lead();
         #make_csi_pad();
         #make_cormo_shield();
         #  make_babar_crystal();
         #make_cry_module_up();
         #make_cry_module_down();
         make_cry_module();
             make_iveto;
          if (($flag_outer_lead) eq "1") {make_lead};
          if (($flag_inner_lead) eq "1") {make_in_lead};
             make_oveto;
             # make_flux_cosmic_sph();

         # make_flux_cosmic_cyl;
         # make_flux_cosmic_rec;
         
         #  make_cry_module_II();
        }
     
        }
 }









1;

