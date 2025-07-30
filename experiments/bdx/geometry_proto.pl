use strict;
use warnings;

our %configuration;

my $degrad = 0.01745329252;
my $cic    = 2.54;

my $crystal = "CsI";

#
my $ncol = 3; # Number of columns (X width)      #             ||
my $nrow = 3; # Number of rows (Y width)         #             \/
my $ndep = 2; # Number of (ncol x nrow) matrixes

# parameters for vertical crystal configuration (only considered if flag on)
my $nplane = 6;         # number of vertical planes
my $plane_side = 3;     # number of crystals side by side in the plane
my $plane_depth = 1;    # each plane is made by just a row of crystals
my $matrix_side = 33.0/2.0;     # the side of the "plane" is the longest dimension i.e. the crystal lenght

#CRYSTAL PARAMETERS
# spacings
my $cr_mylar = 0.005/2.0; # Mylar wrapping thikness

# alveolus shape
my $alv_s = 8.0/2.0; #side
my $alv_l = 33.0/2.0; #length
my $alv_t = 0.2/2.0; #thickness

# crystal parameters -- all in cm -- CsI crystals
my $crs_x=4.7/2 ; # Endcap: short side X (4.3+3.9)/2=4.1cm
my $crs_y=4.8/2 ; # Endcap: short side Y 4.7
my $crs_X=5.8/2 ; # Endcap: long side X (5+4.6)/2=4.8cm
my $crs_Y=6.0/2 ; # Endcap: long side Y 5.4
my $crs_l=31.6/2.; # Endcap: lenght side Y 32.5

#plastic cap dimension:
my $cap_tk = ($alv_l - $crs_l - $cr_mylar)/2.0; #cap fills the missing space

# LEAD shielding dimensions
# Mounting of sides:
#   ______
#  | |____|
#  | |  | |
#  |----| |
#   ------
# Front and back panel are placed in front

my $pb_t_x = 34.8/2.0;
my $pb_t_y = 1.0/2.0;
my $pb_t_z = 68.6/2.0;

my $pb_s_x = 1.0/2.0;
my $pb_s_y = 33.8/2.0;
my $pb_s_z = 68.6/2.0;

my $pb_f_x = 35.8/2.0;
my $pb_f_y = 34.8/2.0;
my $pb_f_z = 1.0/2.0;

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
    #$detector{"material"}   = "KryptoniteLight"; 
    print_det(\%configuration, \%detector);
}

# Open a file to write translation table
my $tt_title = "TT.txt";
if($configuration{"vertical_crystals"} == 1){
    $tt_title = "TT_vertical.txt";
}
open(my $traslation_table, '>', $tt_title) or die "Can't open output file, sorry";

sub make_crystal{
    my $ix = $_[0];
    my $iy = $_[1];
    my $iz = $_[2];
    
    my $rotx = 0;
    my $roty = 0;
    my $rotz = 0;
    
    my $rot=0;
    
    # All parallelepiped are defined as
    #     _ _________
    #  y |_|_________|
    #     x     l
    
    # For crystals the capital letter refers to the larger value of the face side
    
    # Wrapped crystals
    my $wr_x = $crs_x+$cr_mylar;
    my $wr_X = $crs_X+$cr_mylar;
    my $wr_y = $crs_y+$cr_mylar;
    my $wr_Y = $crs_Y+$cr_mylar;
    my $wr_l = $crs_l+$cr_mylar;
    
    # Air around crystals
    my $air_X = $alv_s-$alv_t; #$wr_X+$cr_air;
    my $air_Y = $alv_s-$alv_t; #$wr_Y+$cr_air;
    my $air_l = $alv_l-2.0 * $cap_tk;
    
    # Crystal alveolus
    my $al_X = $alv_s; #$air_X+$cr_alv;
    my $al_Y = $alv_s; #$air_Y+$cr_alv;
    my $al_l = $alv_l;
    
    # Shift between alveoli
    my $alv_dx = 2.0*$alv_s+0.5; # !!! Here I am assuming a gap of 1 mm between the alveoli
    my $alv_dy = 2.0*$alv_s+0.5;
    my $alv_dz = 2.0*$alv_l+0.5;
    
    
    my %detector = init_det();
        
        $detector{"name"}        = "cry_alveol_$ix"."_"."$iy"."_"."$iz";
        $detector{"mother"}      = "proto_mother";
        $detector{"description"} = "Al container_$ix"."_"."$iy"."_"."$iz";
        $detector{"color"}       = "00ffff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        # pos =  center of module - half size of block + alv. side
        my $X = -($ncol-1)*$alv_dx/2.0 + $alv_dx*$ix;#-($ncol-1)*$al_X + 2*$al_X*$ix;
        my $Y = -($ncol-1)*$alv_dy/2.0 + $alv_dy*$iy;#-($nrow-1)*$al_Y + 2*$al_Y*$iy;
        my $Z = $iz * $alv_dz - 0.5*$alv_dz*($ndep-1) ;
        
        if($configuration{"vertical_crystals"} eq 1){   #in this case ix is the number of the colunm/row, iy the position along the row and iz the plane
            if($iz % 2 == 0){
                $rotx = 90;
                $X = -($plane_side-1)*$al_X+ 2*$al_X*$ix;
                $Y = -($plane_depth-1)*0.5*$alv_dz + $alv_dz*$iy;
                $Z = -0.5*$alv_dy*($nplane-1) + $iz*$alv_dy;
            }elsif($iz % 2 != 0){
                $roty = 90;
                $X = -($plane_depth-1)*0.5*$alv_dz + $alv_dz*$iy;
                $Y = -($plane_side-1)*$al_X+ 2*$al_X*$ix;
                $Z = -0.5*$alv_dy*($nplane-1) + $iz*$alv_dy;
                print("$X $Y $Z $plane_depth \n");
                
            }
        }
        
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "$rotx*deg $roty*deg $rotz*deg";
        $detector{"dimensions"}  = "$alv_s*cm $alv_s*cm $alv_l*cm";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
    
    # Plastic caps
    # !!! NOTE that the material has to be changed
    $detector{"name"}        = "cry_front_cap_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cry_alveol_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "Front cap crystal_$ix"."_"."$iy"."_"."$iz";
    $detector{"color"}       = "000000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = 0;
    $Y = 0;
    $Z = $alv_l-$cap_tk;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$alv_s*cm $alv_s*cm $cap_tk*cm";
    $detector{"material"}    = "G4_POLYSTYRENE";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "cry_back_cap_$ix"."_"."$iy"."_"."$iz";
    $detector{"mother"}      = "cry_alveol_$ix"."_"."$iy"."_"."$iz";
    $detector{"description"} = "Back cap crystal_$ix"."_"."$iy"."_"."$iz";
    $detector{"color"}       = "000000";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = 0;
    $Y = 0;
    $Z = -$Z;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$alv_s*cm $alv_s*cm $cap_tk*cm";
    $detector{"material"}    = "G4_POLYSTYRENE";
    print_det(\%configuration, \%detector);
        
        # Air layer
        %detector = init_det();
        $detector{"name"}        = "cry_air_$ix"."_"."$iy"."_"."$iz";
        $detector{"mother"}      = "cry_alveol_$ix"."_"."$iy"."_"."$iz";
        $detector{"description"} = "Air $ix"."_"."$iy"."_"."$iz";
        $detector{"color"}       = "00fff1";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "0*cm 0*cm 0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$air_X*cm $air_Y*cm $air_l*cm";
        $detector{"material"}    = "G4_AIR";
        #$detector{"material"}   = "KryptoniteLight";
        print_det(\%configuration, \%detector);
        
        # Mylar wrapping
        %detector = init_det();
        $detector{"name"}        = "cry_mylar_$ix"."_"."$iy"."_"."$iz";
        $detector{"mother"}      = "cry_air_$ix"."_"."$iy"."_"."$iz";
        $detector{"description"} = "Mylar wrapping_$ix"."_"."$iy"."_"."$iz";
        $detector{"color"}       = "848484";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Trd";
        $detector{"pos"}         = "0*cm 0*cm 0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$wr_x*cm $wr_X*cm $wr_y*cm $wr_Y*cm $wr_l*cm";
        $detector{"material"}    = "bdx_mylar";
        print_det(\%configuration, \%detector);
        
        # Crystals
        %detector = init_det();
        $detector{"name"}        = "crystal_$ix"."_"."$iy"."_"."$iz";
        $detector{"mother"}      = "cry_mylar_$ix"."_"."$iy"."_"."$iz";
        $detector{"description"} = "Crystal_$ix"."_"."$iy"."_"."$iz";
        $detector{"color"}       = "00fff2";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Trd";
        $detector{"pos"}         = "0*cm 0*cm 0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$crs_x*cm $crs_X*cm $crs_y*cm $crs_Y*cm $crs_l*cm";
        $detector{"material"}    = "CsI_Tl";
        $detector{"sensitivity"} = "crs";
        $detector{"hit_type"}    = "crs";
        
        my $i_ix = $ix+1;
        my $i_iy = $iy+1;
        my $i_iz = $iz+1;
        
        # for geometry with rotated crystals I need to correct the labels
        #in this case ix is the number of the colunm/row, iy the position along the row and iz the plane
        if($configuration{"vertical_crystals"} eq 1){
            if($iz % 2 != 0){
                $i_ix = $iy+1;
                $i_iy = $ix+1;
            }
        }
        $detector{"identifiers"} = "sector manual 0 xch manual $i_ix ych manual $i_iy zch manual $i_iz";
        #$detector{"identifiers"} = "sector manual $i_im xch manual $X ych manual $Y zch manual $Z";
        print_det(\%configuration, \%detector);
        
        # write translation table
        # entry = module, idx, idy, idz, xc, yc, zc, lx, ly, lz
        if($configuration{"vertical_crystals"} eq 1){   #in this case ix is the number of the colunm/row, iy the position along the row and iz the plane
            if($iz % 2 == 0){
                print $traslation_table "$i_ix $i_iy $i_iz $X $Y $Z $crs_X $crs_l $crs_Y\n";
            }elsif($iz % 2 != 0){
                print $traslation_table "$i_ix $i_iy $i_iz $X $Y $Z $crs_l $crs_Y $crs_X \n";
            }
        }else{
            print $traslation_table " $i_ix $i_iy $i_iz $X $Y $Z $crs_X $crs_Y $crs_l \n";
            print "printing data \n"
        }
    
    
}

# BEGIN inner lead
sub make_lead
{
    # Assuming fixed parameters for lead
    # Dimensions defined as follows
    #     _ _________
    #  y |_|_________|
    #     x     z
    
    my $im = $_[0];
    
    my %detector = init_det();
    
    $detector{"mother"}      = "proto_mother";
    $detector{"material"}    = "G4_Pb";
    
    # Top/bottom
    $detector{"name"}        = "lead_top";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $X = -$pb_s_x;
    my $Y = $pb_s_y;
    $detector{"pos"}         = "$X*cm $Y*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_t_x*cm $pb_t_y*cm $pb_t_z*cm";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "lead_bottom";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = -$X;
    $Y = -$Y;
    $detector{"pos"}         = "$X*cm $Y*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_t_x*cm $pb_t_y*cm $pb_t_z*cm";
    print_det(\%configuration, \%detector);
    
    # Front/back
    
    #holes in the front piece of lead
    my $h_w = 2.2/2.0; # width of hole
    my $h_h = 1.0/2.0; #height of hole
    
    #upstream lead has two holes so I split it
    $detector{"name"}        = "lead_upstream_top";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $Z = $pb_s_z+$pb_f_z;
    $detector{"pos"}         = "0*cm $h_h*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    my $up_l = $pb_f_y-$h_h;
    $detector{"dimensions"}  = "$pb_f_x*cm  $up_l*cm $pb_f_z*cm";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "lead_upstream_bottom";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $low_y = (-$pb_f_y+$h_h);
    $detector{"pos"}         = "0*cm $low_y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    my $low_l = $pb_f_x-2.0*$h_w;
    $detector{"dimensions"}  = "$low_l*cm $h_h*cm $pb_f_z*cm";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "lead_downstream";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $Z = -$Z;
    $detector{"pos"}         = "0*cm 0*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_f_x*cm $pb_f_y*cm $pb_f_z*cm";
    print_det(\%configuration, \%detector);
    
    #right and left
    $detector{"name"}        = "lead_right";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $Y = $pb_t_y;
    $X = $pb_t_x;
    $detector{"pos"}         = "$X*cm $Y*cm 0*cm";
    $detector{"pos"}         = "$X*cm $Y*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_s_x*cm $pb_s_y*cm $pb_s_z*cm";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "lead_left";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $Y = -$Y;
    $X = -$X;
    $detector{"pos"}         = "$X*cm $Y*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_s_x*cm $pb_s_y*cm $pb_s_z*cm";
    print_det(\%configuration, \%detector);
}
#END inner lead

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

my $geom_file_name = "geom.txt";
open(my $output, '>', $geom_file_name) or die "Can't open output file, sorry";

sub make_crs{
    if($configuration{"vertical_crystals"} eq 1){
        for(my $iz=0; $iz< $nplane; $iz++){
            for(my $iy=0; $iy<$plane_depth; $iy++){
            for(my $ix = 0; $ix< $plane_side; $ix++){
                    make_crystal($ix, $iy, $iz, 0, 0, 0);
                }
            }
        }
    }else{
        
        for(my $iz=0; $iz< $ndep; $iz++){
            for(my $ix = 0; $ix< $ncol; $ix++){
                for(my $iy=0; $iy<$nrow; $iy++){
                    make_crystal($ix, $iy, $iz, 0, 0, 0);
                }
            }
        }
    }
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
    #$detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
    
    
    
    
    make_lead;
    make_iveto;
    make_oveto;
    make_crs;
}

sub make_all
{
    make_main;
    #make_proto;
}


1;
