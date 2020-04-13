package MyModule::GlobalVar;

use warnings;
use strict;
use File::Spec;
use File::Basename;
use Cwd qw(abs_path);
use Config::IniFiles;
use Exporter;
use vars qw(@ISA @EXPORT);


@ISA = qw(Exporter);
@EXPORT = qw($PATH_BIN $PATH_CONF $PATH_DB $PATH_TOOLS $CONFIG_FILE $DB_LIST_KRAKEN
             $REF_H_CONFIG $PYTHON3 $MAX_AVAILEBLE_MEM $REF_H_DBPATH_KRAKEN $REF_H_DBMEM_KRAKEN);


# Path
my $modulePath = dirname(File::Spec->rel2abs(__FILE__));
my $path_pipeline  = abs_path("$modulePath/../../..");
our $PATH_BIN   = "$path_pipeline/bin";
our $PATH_CONF  = "$path_pipeline/conf";
our $PATH_DB    = "$path_pipeline/db";
our $PATH_TOOLS = "$path_pipeline/tools";

# File
our $CONFIG_FILE    = "$PATH_CONF/config.ini";
our $DB_LIST_KRAKEN = "$PATH_DB/db_Kraken/db.list";


# Config file
&getCfg($CONFIG_FILE);


## Kraken db info
#&getKrakenDbInfo($DB_LIST_KRAKEN);


#===============================================================================
#   Subroutine
#===============================================================================
sub getCfg {
    my ($cfg_file, ) = @_;
    tie my %ini, "Config::IniFiles", (-file, => $cfg_file);

    our $REF_H_CONFIG      = \%ini;
    our $PYTHON3           = $ini{"software"}{"python3"};
    our $MAX_AVAILEBLE_MEM = $ini{"resource"}{"MaxAvailableMem"};
}

sub getKrakenDbInfo {
    my ($dbList, ) = @_;
    $dbList = abs_path($dbList);
    my $dbDir = dirname($dbList);
    
    my %h_dbPath = ();
    my %h_dbMem = ();
    open IN, "$dbList" or die $!;
    while (<IN>) {
        chomp;
        
        my($dbName, $dbPath, $dbMem) = split /\s+/;
        $dbPath = "$dbDir/$dbPath"; 
        $h_dbPath{$dbName} = $dbPath;
        $h_dbMem{$dbName}  = $dbMem;
    }
    close IN;
    
    our $REF_H_DBPATH_KRAKEN = \%h_dbPath;
    our $REF_H_DBMEM_KRAKEN  = \%h_dbMem;
}


1;
