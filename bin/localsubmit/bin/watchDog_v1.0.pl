#!/usr/bin/perl

#===============================================================================
#   FileName: watchDog_v1.0.pl
#   Author  : mayubin@genomics.cn
#   Version : 1.0
#   Date    : 2017-10-17
#   Description: The monitor program to run task parallel on local machine
#===============================================================================

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use Cwd qw(abs_path);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../lib/perl5";
use Parallel::ForkManager;
use MyModule::GlobalVar qw($MAX_AVAILEBLE_MEM);


my($taskList, $mem, $num_paral, $num_line, $help);
GetOptions(
    "mem=s"    => \$mem,
    "maxjob=i" => \$num_paral,
    "lines=i"  => \$num_line,
    "help|h"   => \$help,
);
pod2usage(-exitval => 1, -verbose => 1) if ($help || @ARGV < 1);
($taskList, ) = @ARGV;
$taskList = abs_path($taskList);
$mem       ||= 2;
$mem = &convMemToM($mem);
$num_paral ||= 3;
$num_line  ||= 1;

# Get the max parallel number based on the max available memory
my $max_paral = int(&convMemToM($MAX_AVAILEBLE_MEM) / $mem);
$num_paral = $num_paral < $max_paral ? $num_paral : $max_paral;


#===============================================================================
#   Check the total memory
#===============================================================================
my($total_mem, $valid_mem) = &getMemInfo();
if ($mem > $total_mem) {
    die "Error: No enough total memory, Total_Mem: $total_mem M, Request_Mem: $mem M\n";
}


#===============================================================================
#   Global variables
#===============================================================================
# Log file for watchDog and log directory for all tasks
my $logFile  = "$taskList.$$.log";
my $scptDir  = "$taskList.$$.watchDog";
mkdir $scptDir;


#===============================================================================
#   The main code
#===============================================================================
printLog2File("Task list: $taskList", $logFile);
printLog2File("Request memory for each task : $mem M", $logFile);
printLog2File("Number of parallel: $num_paral", $logFile);

# Get the shell scipts of the tasks
my @arr_taskList = &readTaskList($taskList, $num_line, $scptDir);
printLog2File("Total number of task: " . scalar(@arr_taskList), $logFile);

# Execute the tasks parrallel
my $pm = new Parallel::ForkManager($num_paral);
for (my $i = 0; $i < @arr_taskList; $i++)
{
    my $taskSh = $arr_taskList[$i];
    
    while(1) {
        ($total_mem, $valid_mem) = &getMemInfo();
        last if ($valid_mem > $mem);
        printLog2File("Waiting for memory, avaliable_Mem: $valid_mem M, Request_Mem: $mem M", $logFile);
        sleep(300);
    }
    
    $pm->start and next;
    &runTask($taskSh);
    $pm->finish;
}
$pm->wait_all_children;

printLog2File("Finish all tasks", $logFile);


#===============================================================================
#   Subroutine
#===============================================================================
# Convert the value of --mem to Mbp
sub convMemToM {
    my($mem, ) = @_;
    
    if ($mem =~ /^\d+$/) {
        return $mem * 1024;
    } elsif ($mem =~ /^(\d+)(g|G)$/) {
        return $1 * 1024;
    } elsif ($mem =~ /^(\d+)(m|M)$/) {
        return $1;
    } else {
        die "Error: Incorrect value for --mem\n";
    }
}

# Get the absolute path in the task
sub getAbsPath {
    my($oldTask, $path) = @_;
    
    my @task = split(/\s+/, $oldTask);
    for (my $i = 0; $i < @task; $i++) {
        if ($task[$i] =~ /^(\.\/|\.\.\/)/) {
            $task[$i] = "$path/$task[$i]";
        }
    }
    return join(" ", @task);
}

# Read the tasks in task.list into a array
sub readTaskList
{
    my($list, $num_line, $scptDir) = @_;
    
    my @arr_taskSh = ();
    my $listDir = dirname($list);
    
    open IN, "$list" or die "Can't open $list: $!\n";
    my $numInSh = 0;
    my $taskSeq = 1;
    while (<IN>)
    {
        chomp;
        
        my $newTask = &getAbsPath($_, $listDir);
        
        my $taskSh = sprintf("$scptDir/tast_%04d.sh", $taskSeq);
        unlink $taskSh if (-e $taskSh && $numInSh == 0);
        
        open SH, ">>$taskSh" or die "Can't open $taskSh: $!\n";
        print SH "$newTask\n";
        close SH;
        
        $numInSh++;
        if ($numInSh == $num_line) {
            push @arr_taskSh, $taskSh;
            $taskSeq++;
            $numInSh = 0;
        }
    }
    close IN;
    
    return @arr_taskSh;
}

# Print log message to file
sub printLog2File {
    my($msg, $file, ) = @_;

    my $time_log = strftime("[%Y-%m-%d %H:%M:%S]", localtime());
    system("echo '$time_log: $msg' >> $file");
}


# Execute the task
sub runTask
{
    my($taskSh, ) = @_;
    
    printLog2File("Start: $taskSh", $logFile);
    
    my $outLog  = "$taskSh.$$.o";
    my $errLog  = "$taskSh.$$.e";
    
    my $exit_val = system("sh $taskSh >$outLog 2>$errLog");
    
    unless ($exit_val == 0) {
        printLog2File("Error: There is error when execute task: $taskSh", $logFile);
    }
    
    printLog2File("Finish: $taskSh", $logFile);
}

# Get the memory information of the machine
sub getMemInfo {
    my $file_memInfo = "/proc/meminfo";
    
    my($total_mem, $valid_mem) = (0, 0);
    open INFO, "$file_memInfo" or die "Can't open $file_memInfo: $!\n";
    while (<INFO>) {
        chomp;
        
        $total_mem =  $1 if (/MemTotal:\s+(\d+)\s+kB/);
        $valid_mem += $1 if (/MemFree:\s+(\d+)\s+kB/);
        $valid_mem += $1 if (/Buffers:\s+(\d+)\s+kB/);
        $valid_mem += $1 if (/Cached:\s+(\d+)\s+kB/);
    }
    close INFO;
    
    $total_mem /= 1024 unless ($total_mem == 0);
    $valid_mem /= 1024 unless ($valid_mem == 0);
    
    return (int($total_mem), int($valid_mem));
}


__END__


=pod

=head1 NAME

watchDog_v1.0.pl - The monitor program to run task parallel on local machine

=head1 VERSION

v1.0

=head1 SYNOPSIS

perl watchDog.pl <task.list> [options]

=head1 ARGUMENTS

=over 8

=item B<task.list> <file>

The list of tasks which need to be executed parallel. The path in tasks must be 
the absolute path or the relative path start with './' or '../'.

=back

=head1 OPTIONS

=over 8

=item B<--mem> <Str>

The request memory for each task, the default unit is G, you can also us M. E.g: 
2, 2g, 2G, 2000m, 2000M. [2G]

=item B<--maxjob> <INT>

The max parallel number. [3]

=item B<--lines> <INT>

The number of lines for one tast. [1]

=item B<--help> | B<-h>

Print this information.

=back

=cut
