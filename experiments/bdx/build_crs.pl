#!/usr/bin/perl -w

use strict;
use lib ("$ENV{GEMC}/api/perl");
use utils;
use parameters;
use geometry;
use hit;
use bank;
use math;
use materials;

use Math::Trig;

# Help Message
sub help()
{
	print "\n Usage: \n";
	print "   bdx.pl <configuration filename> [crs_zrot]\n";
 	print "   Will create the BDX using the variation specified in the configuration file\n";
 	print "   Optional 2nd argument overrides crs_zrot (deg), must be 0 or 90, e.g. ./build_crs.pl config_1.dat 90\n";
 	print "   Note: The passport and .visa files must be present to connect to MYSQL. \n\n";
	exit;
}

# Make sure the argument list is correct
if( scalar @ARGV < 1 || scalar @ARGV > 2)
{
	help();
	exit;
}


# Loading configuration file and paramters
our %configuration = load_configuration($ARGV[0],[0]);

print "CONFIGURATION IS ::: ",$configuration{"variation"},"\n";

# Optional 2nd argument: crystal Z-rotation (deg). Only 0 or 90 allowed.
if( scalar @ARGV == 2 )
{
	if( $ARGV[1] != 0 && $ARGV[1] != 90 )
	{
		print "\n ERROR: crs_zrot must be either 0 or 90 (got '$ARGV[1]').\n\n";
		exit;
	}
	$configuration{"crs_zrot"} = $ARGV[1];
	print "CRS_ZROT OVERRIDE ::: $ARGV[1] deg\n";
}

# One can change the "variation" here if one is desired different from the config.dat
# $configuration{"variation"} = "myvar";

# To get the parameters proper authentication is needed.
# our %parameters    = get_parameters(%configuration);

# Loading FT specific subroutines
require "./hit.pl";
require "./bank.pl";
#require "./geometry_with_dump.pl";
#require "./geometry_crs.pl";
#require "./geometry_crs8.pl";
require "./geometry_crs60.pl";
#require "./geometry_floors.pl";
require "./materials.pl";

define_banks();
define_bdx_hits();
define_bdx_materials();
make_all(); 
    
    
