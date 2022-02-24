#!/usr/bin/perl -w 
use strict;
use warnings;
use List::Util qw( max min );
use Config;
use threads;
use threads::shared;
#use Data::Dumper;
#
#  Author: Adam Grofe
#  Email: grofe001@umn.edu
#  Created: 3-12-2015
#
#  Purpose: To streamline the process of evaluating spectra from the xdip files
#           
#  Input files needed: xdip files from the simulations
#
#  Execution: ./spectral-analysis.pl
#
#
#==============================================================================
# DECLARATIONS
# -----------------------------------------------------------------------------
# THINGS TO EDIT
my $output_file :shared = "spectral-analysis.log";
my $xdip_chunks :shared = 200000; # LENGTH OF TIME SEGMENT IN ACF CALC
my $acf_points :shared = 65536;   # NUMBER OF POINTS IN ACF 
my $nprocessors = 7;      # NUMBER OF PROCESSES GENERATED DURING ACF CALC
my $band_width = 6;       # BANDWIDTH FOR SMOOTHING FUNCTION IN R
my $norm_min = 2000;      # RANGE THROUGH WHICH TO NORMALIZE PEAK
my $norm_max = 2400;      # NORMALIZED PEAK IS IN THE WHOLE SPEC DATA FILE
my $remove_charge_transfer = "yes"; # EITHER "yes" or "no"
my $windows = 1;          # Time resolve the spectra using several windows

#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
# THINGS NOT TO EDIT
#------------------------------------------------------------------------------
my $i;
my $j;
my @xdip_files;
my $NTFR = 0;           # FREQ TO SKIP FRAMES
my $basename :shared;
my @dipole_x1;
my @dipole_x2;
my @dipole_y1;
my @dipole_y2;
my @dipole_z1;
my @dipole_z2;
my $timesteps;
my $hostname;
my @temp_files;
my $max = 0;
my $success = 0;
my $basename_spectrum;
my $start;
my $end;
my $dir_name;

#==============================================================================
# EVALUATE ARGUMENTS
# -----------------------------------------------------------------------------
my $early;
my $jump_to_r;
my $jump_to_splice;
my $jump_to_acf;

if( @ARGV ){
  my $nargs = @ARGV;
    if( $ARGV[0] =~ m/\s*-t\s*/ ){
      $early = 1;
    }
    elsif( $ARGV[0] =~ m/\s*-s\s*/){
       $jump_to_splice = 1;
    }
    elsif( $ARGV[0] =~ m/\s*-R\s*/){
         $jump_to_r = 1;
    }
    elsif( $ARGV[0] =~ m/\s*-acf\s*/){
        $jump_to_acf = 1;
        $basename = $ARGV[1];
        $start = $ARGV[2];
        $end = $ARGV[3];
        $dir_name = $ARGV[4];
    }
}
unless( $jump_to_acf ){
    unlink "$output_file";
}
#=============================================================================
# In order to conserve memory the parallel acf calculation will invoked 
# separately from the rest of the program.
if( $jump_to_acf ){
   close LOG;
   calc_acf();
   exit(0);
}
#==============================================================================
# INITIALIZATION AND GENERIC SANITY CHECKS
# -----------------------------------------------------------------------------

open (LOG, '>>', "$output_file") or #<----No semi-colon by design
     die "ERROR: Unable to open $output_file for output\n";

my $date = `date`;
print LOG "$date\n";
print LOG "PARAMETERS: XDIP-Segments: $xdip_chunks \n ACF-Points: $acf_points \n PROCESSORS: $nprocessors \n";
print LOG " BAND WIDTH: $band_width \n NORMALIZATION Min: $norm_min \n NORMALIZATION Max: $norm_max \n";
print LOG " REMOVE CT: $remove_charge_transfer \n NUMBER of Windows: $windows \n \n";

#$hostname = `hostname`;
#unless( $hostname =~ m/\s*node.*/ ){
#  print LOG "Warning: This program needs to be executed on Itasca or a computer\n";
#  print LOG "         that has R installed.\n";
#  print "Warning: This program needs to be executed on Itasca or a computer\n";
#  print "         that has R installed.\n";
#}
my $r_test = `R --vanilla < /dev/null `;
unless( $r_test ){
    print LOG "ERROR: R is not loaded \n";
    print LOG "       Run \"module load R\" or install R \n";
    print "ERROR: R is not loaded \n";
    die   "       Run \"module load R\" or install R \n";
}
unless( -e "acf.exe" ) {
   print LOG "ERROR: acf.exe is not located in this directory\n";
   die "ERROR: acf.exe is not located in this directory\n";
}

$Config{useithreads} or die "ERROR: ithreads are not enabled on this version of Perl\n       Recompile perl with ithreads enabled\n";

if( $norm_min > $norm_max){
    die "ERROR: Norm_min is greater than Norm_max\n       Norm_min: $norm_min\n       Norm_max: $norm_max\n";
}
elsif( $norm_min == $norm_max ){
    print "WARNING: Norm_min is equal to Norm_max:\n";
    print "         Norm_min = $norm_min\n";
    print "         Norm_max = $norm_max\n";
}

#==============================================================================
# OPEN THE XDIP FILES AND DETERMINE THE BASENAME OF THE FILES
# -----------------------------------------------------------------------------

@xdip_files = glob '*.xdip';

unless( @xdip_files ) {
   die "ERROR: There are no .xdip files in this directory \n";
}
print LOG "XDIP FILES:\n";

my @basename_array;
my @sim_number;
foreach $i ( @xdip_files ){
   if( $i =~ m/\s*(\S+?)([0-9]+)\.xdip/ ){
     push @basename_array, $1;
     push @sim_number, $2; 
   }
   else{
      print "ERROR: XDIP glob failure \n";
      print LOG  "ERROR: XDIP glob failure \n";
      print LOG  "       please rename xdip files to the following form:\n";
      print LOG  "       basename_1.xdip basename_2.xdip\n";
      print "       please rename xdip files to the following form:\n";
      die "       basename_1.xdip basename_2.xdip\n";
   }
}

foreach $i ( @basename_array ){
  foreach $j ( @basename_array ){
     if( $i eq $j ){
        next;
     }
     else{
        die "ERROR: There are different kinds of xdip files in this directory: $i  and $j \n";
     }
  }
}
$basename = $basename_array[0];
my $sim_min = min @sim_number;
my $sim_max = max @sim_number;

#==============================================================================================
# PERFORM JUMP TO ANOTHER SECTION OF THE PROGRAM
#----------------------------------------------------------------------------------------------
if( $jump_to_r ){
    goto R;
}
if( $jump_to_splice ){
    $basename_spectrum = $basename . "_spectrum";
    goto SPLICE;
}
#==============================================================================================
# DETERMINE THE FREQUENCY TO REMOVE DOUBLES
#----------------------------------------------------------------------------------------------
open (XDIP, '<', $basename . $sim_min . '.xdip') or 
     die "ERROR: Unable to open $basename $sim_min .xdip file\n";
my $oldx = 0;
my $oldy = 0;
my $oldz = 0;
my $frame = 1;
my $last = 1;
my %frame_hist;
my $diff;

while(<XDIP>){
    if(m/\s*(-*[0-9]+\.[0-9]{4})[0-9]{1}\s*(-*[0-9]+\.[0-9]{4})[0-9]{1}\s*(-*[0-9]+\.[0-9]{4})[0-9]{1}\s*(-*[0-9]+\.[0-9]{4})[0-9]{1}\s*(-*[0-9]+\.[0-9]{4})[0-9]{1}\s*(-*[0-9]+\.[0-9]{4})[0-9]{1}\s*/){
      if($1 == $oldx and $2 == $oldy and $3 == $oldz){
          $diff = $frame - $last;
          $frame_hist{$diff}++;
          $last = $frame;
      }
   }
   $frame++;
   $oldx = $1;
   $oldy = $2;
   $oldz = $3;
 
}
close XDIP;

my @sorted_hist = sort{ $frame_hist{$b} <=> $frame_hist{$a} } keys %frame_hist;

$NTFR = shift @sorted_hist;

if( %frame_hist ){
  if( $frame_hist{$NTFR} < 5 ){
    print "WARNING: The count for the frequency histogram is less then 5\n";
    print "         Freq: $NTFR    Count: $frame_hist{$NTFR} \n";
    print "         Frequency is being set to 0\n";
    $NTFR = 0;
 }
}
else{
   $NTFR = 0;
}
#==============================================================================================
# FILTER OUT THE DUPLICATES
#----------------------------------------------------------------------------------------------

print LOG "Analyzing simulations $sim_min through $sim_max\n";
print LOG "Simulation basenames are $basename \n\n";

my $basename_filtered = $basename . "filtered";

print LOG "---------------------------------------------------------------------------------\n";
print LOG "    FILTER XDIP FILES FOR DUPLICATES AND SPLICE TOGETHER\n";
print "    FILTER XDIP FILES FOR DUPLICATES AND SPLICE TOGETHER\n";

unless( -d "xdip" ){
    mkdir "xdip", 0770 or die "ERROR: Cannot create xdip directory\n";
 }

if ( -e "xdip/$basename_filtered.xdip" ){
    unlink "xdip/$basename_filtered.xdip";
}
open (XOUT, '>', "xdip/$basename_filtered.xdip") or 
     die "ERROR: Unable to open xdip/$basename_filtered.xdip for output.\n";

for($i=$sim_min; $i<=$sim_max; $i++ ){
    my $filename = $basename . $i . ".xdip";
    open ( XDIP, '<', "$filename" ) or die "ERROR: Unable to open $filename \n";
    $j = 0;
    my $skip = $NTFR;
    while(<XDIP>){
      $j++;
      if( $j == $skip){
          $skip = $skip + $NTFR;
          next;
      } 
      print XOUT $_; 
    }
    print LOG "$filename added to xdip/$basename_filtered.xdip file\n";
}
close XOUT;

print LOG "\n";
print LOG "XDIP files were filtered to remove duplicate frames at every $NTFR frames\n";
print LOG "XDIP Frequency data:\n";
print LOG "  Freq   Count \n";
print LOG "   $_    $frame_hist{$_} \n" for( keys %frame_hist );


open (XIN, '<', "xdip/$basename_filtered.xdip" ) or 
     die "ERROR: Unable to open xdip/$basename_filtered.xdip for input.\n";

while(<XIN>){
   $timesteps++;
}

#==============================================================================
# PRINT THE NUMBER OF TIMESTEPS AND DETERMING THE ANALYSIS WINDOWS
# -----------------------------------------------------------------------------
print LOG "The number of timesteps in the xdip/$basename_filtered.xdip file is: $timesteps\n";
if( $timesteps < $xdip_chunks) { 
    $xdip_chunks = $timesteps;
    print LOG "Number of Timesteps is less than Xdip Chunks \n";
    print LOG "    Xdip chunks being assigned to the number of timesteps \n";
}

my $dt;
my @start;
my @end;
if( $windows > 1 ){
   print LOG "Time Resolution requested\n Windows:\n";
   $j = 0;
   $dt = int( $timesteps/($windows) ) ;
   for( $i=0; $i<$windows; $i++){
      $start[$i] = $j;
      $end[$i] = $j + $dt;
      print LOG "  $start[$i]   $end[$i] \n";
      $j += $dt;
   }
   $start[$windows] = 0;
   $end[$windows] = $timesteps;
   print LOG "  $start[$windows]  $end[$windows]\n";
}
else{
   $start[0] = 0;
   $end[0] = $timesteps;
}

# IF THE -t FLAG IS USED ON INVOCATION THEN THE SCRIPT TERMINATES HERE
if( $early ){
    print LOG "\nNORMAL TERMINATION\n";
    exit(0);
}
#==============================================================================
# OPEN THE FILTERED XDIP FILE AND LOAD THE DIPOLE DATA
# -----------------------------------------------------------------------------

$i = 0;
seek XIN, 0, 0; # REWINDS THE XIN UNIT TO THE BEGINNING
while(<XIN>){
   if( m/\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*/ ) {
       $dipole_x1[$i] = $1;
       $dipole_y1[$i] = $2;
       $dipole_z1[$i] = $3;
       $dipole_x2[$i] = $4;
       $dipole_y2[$i] = $5;
       $dipole_z2[$i] = $6;
       $i++;
       $success = 1;
   }
}

close XIN;
if( $success != 1 ){
   print "ERROR: Filtered xdip file read unsuccessful \n";
   print "       Make sure the input xdip files are populated\n";
   print "       OR Contact the maintainer of this program\n";
   print '       Maintainer: grofe001@umn.edu';
   print "\n";
   die "\n";
}
print LOG "Filtered XDIP file loaded into memory\n";

#==============================================================================
#==============================================================================
# SPECTRA CALCULATIONS BEGIN HERE
#==============================================================================
#==============================================================================
my $t;
my $tmax = @start;
for($t=0; $t<$tmax; $t++){


print LOG "============================================================================\n";
print LOG "             WINDOW:  $start[$t]    $end[$t] \n";
print LOG "----------------------------------------------------------------------------\n";
print LOG "   GENERATING TEMP FILES FOR ACF PROGRAM\n";
print "    GENERATING TEMP FILES FOR ACF PROGRAM\n";

#==============================================================================
# PRINT OUT THE TEMP FILE FOR THE ACF CALCULATION AND EXECUTE 
# -----------------------------------------------------------------------------

$i = $start[$t];
my $stop = $start[$t] + $xdip_chunks;

$dir_name = "acf-".$start[$t]."-".$end[$t];

 if( -e $dir_name ) {
     unlink glob "$dir_name/*.acf";
 }
 else{
   mkdir $dir_name, 0770 or die "ERROR: Unable to make directory $dir_name \n  $!";
 }

while( $stop <= $end[$t] ){
  open (TEMP, '>', "temp.$basename.$i.xdip") or 
         die "ERROR: Unable to open temp.$basename.$i.xdip \n";
  if( $remove_charge_transfer eq "yes" ){
    print LOG "**Using Dipole vectors with charge transfer removed\n";
    for( $j=$i; $j<$stop; $j++){
      printf TEMP "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n", $dipole_x2[$j],$dipole_y2[$j],
                  $dipole_z2[$j],$dipole_x1[$j],$dipole_y1[$j],$dipole_z1[$j];
    }
  }
  elsif( $remove_charge_transfer eq "no"){
    for( $j=$i; $j<$stop; $j++){
      printf TEMP "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n", $dipole_x1[$j],$dipole_y1[$j],
                  $dipole_z1[$j],$dipole_x2[$j],$dipole_y2[$j],$dipole_z2[$j];
    }
  }
  else{
    print LOG "**Evaluation of \$remove_charge_transfer did not succeed \n";
    print LOG "**Using the default dipole (unaltered)\n";
    for( $j=$i; $j<$stop; $j++){
      printf TEMP "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n", $dipole_x1[$j],$dipole_y1[$j],
                  $dipole_z1[$j],$dipole_x2[$j],$dipole_y2[$j],$dipole_z2[$j];
    }
  }
  close TEMP; 
  print LOG "temp.$basename.$i.xdip successfully written\n";
  push @temp_files, "temp.$basename.$i.xdip";
  $i += $xdip_chunks / 2;
  $stop += $xdip_chunks / 2;
}

print LOG "----------------------------------------------------------------------------\n";
print LOG "   CALLING ACF PROGRAM\n";
print "    CALLING ACF PROGRAM\n";

close LOG; # Do this so that the acf.exe can write output to the log file

system "./spectral-analysis.pl -acf $basename  $start[$t]  $end[$t] $dir_name";
wait;

sub calc_acf{ 
 my @temp_files = glob 'temp.*.*.xdip';

 foreach my $i ( @temp_files ){
   my @thrds;
   my $cnt;
   
   unless( $cnt ){
     $cnt = 0;
   }
   $thrds[$cnt] = threads -> create(\&call_acf, $i );
   $cnt++;

   my $n_thrds = threads->list();
   if( $n_thrds >= $nprocessors ){
       foreach my $k ( threads->list() ){
          $k -> join();
       }
   }   
 }
 my $n_thrds = threads->list();
 if( $n_thrds > 0 ){
    foreach my $k ( threads->list() ){
       $k -> join();
    }
 }
}

sub call_acf {
      my $file = $_[0];
      $file =~ m/temp\.\S+\.([0-9]+)\.xdip/;
      my $j = $1;
      open (ACFIN, '>', "acf-input.$file.temp" ) or die "ERROR: Unable to open acf-input.$file.temp\n";
      print ACFIN "temp.$basename.$j.xdip, $basename.$j.acf, $xdip_chunks, $acf_points\n";
      close ACFIN;
      system "cat acf-input.$file.temp >> $output_file";
      system "cat acf-input.$file.temp | ./acf.exe >> $output_file ";
      wait;
      unlink "acf-input.$file.temp";
      unlink "temp.$basename.$j.xdip";
      system "mv $basename.$j.acf $dir_name/$basename.$j.acf";
}


# Find the Max value for the temp files
my $analysis_min = 1000000000000000000000000000000;
foreach $i (@temp_files){
  if( $i =~ m/temp\.\S+\.([0-9]*)\.xdip/){
    $j = $1;
    if( $j > $max ){
       $max = $j;
    }
    if( $j < $analysis_min ){
       $analysis_min = $j;
    }
  }
}
@temp_files = qw();


open LOG, '>>', "$output_file";


#==============================================================================
# GENERATE R INPUT FILE AND ANALYZE ACF FILES
# -----------------------------------------------------------------------------
R:
my $interval = $xdip_chunks / 2;
$basename_spectrum = $basename . "spectrum-".$start[$t]."-".$end[$t];

print LOG "----------------------------------------------------------------------------\n";
print LOG "    ANALYZE ACF FILES USING R \n";
print "    ANALYZE ACF FILES USING R \n";

open (RIN, '>', "r-input.temp" ) or die "ERROR: Unable to open r-input.temp \n $! \n";

print RIN "bw <- $band_width;\n";
print RIN "xlim <- c(600, 3600);\n";
print RIN "ylim <- c(0, 3);\n";
print RIN "options('help_type'='text');\n";
print RIN "del  <- 1.0/1000; # 1 fs interval\n";
print RIN "\n";
print RIN "# MeOH PM3 Z\n";
print RIN "avgft <- complex($acf_points);\n";
print RIN "count <- 0;\n";
print RIN "for (init in seq($analysis_min, $max, $interval)) {\n";
print RIN "  file <- paste(\'$dir_name/$basename.\', format(init, scientific=F), \'.acf\', sep=\"\");\n";
print RIN "  v    <- read.table(file);\n";
print RIN "  v1   <- v\$V2[1:(length(v\$V1)-1)];\n";
print RIN "  len  <- length(v1);\n";
print RIN "  ft1   <- fft(v1);\n";
print RIN "  lft1 <- length(ft1);\n";
print RIN "  k <- (1:((lft1+1)/2));\n";
print RIN "  lk <- length(k);\n";
print RIN "  nuc <- 1/(2*del);\n";
print RIN "  del.nu <- nuc/lk;\n";
print RIN "  avgft <- avgft + ft1;\n";
print RIN "  count <- count + 1;\n";
print RIN "}\n";
print RIN "avgft  <- avgft/count;\n";
print RIN "m.meoh.pm3.z <- cbind((0:(lk-1))*del.nu/0.03, Re(avgft[k])^2);\n";
print RIN "s.meoh.pm3.z <- ksmooth(m.meoh.pm3.z[,1], m.meoh.pm3.z[,2], \"normal\", bandwidth = bw);\n";
print RIN "\n";
print RIN "min.meoh.pm3.z <- min(s.meoh.pm3.z\$y);\n";
print RIN "\n";
print RIN "write(s.meoh.pm3.z\$x, file=\"temp.x-data.dat\", ncolumns = 1, append = \"FALSE\", sep=\"\\n\")\n";
print RIN "write(s.meoh.pm3.z\$y-min.meoh.pm3.z, file=\"temp.y-data.dat\", ncolumns = 1, append = \"FALSE\", sep=\"\\n\")\n";
print RIN "write(m.meoh.pm3.z[,1], file=\"temp.rough-xdata.dat\", ncolumns = 2, append = \"FALSE\", sep=\"\\n\")\n";
print RIN "write( m.meoh.pm3.z[,2], file=\"temp.rough-ydata.dat\", ncolumns = 1, append = \"FALSE\", sep=\"\\n\")\n";
print RIN "\n";
print RIN "pdf(\'$basename_spectrum.pdf\');\n";
print RIN "plot(0, 0, xlim=c(1000, 3500), ylim=c(0, 1.0), xlab='Wavenumber (cm-1)', ylab='Intensity (arbitrary)');\n";
print RIN "lines(s.meoh.pm3.z\$x, s.meoh.pm3.z\$y - min.meoh.pm3.z, lwd=2, col='blue' );\n";
print RIN "dev.off();\n";

close RIN;

print LOG "\n";
print LOG "   R INPUT:\n\n";
close LOG;
system "cat r-input.temp >> $output_file ";

system "echo -e \"\n\n    R OUTPUT:\n\n\" >> $output_file";

system "R --vanilla < r-input.temp >> $output_file 2>&1";

unlink "r-input.temp";

#==============================================================================
# READ IN THE SPECTRAL DATA AND PLACE THEM IN A SINGLE FILE
# -----------------------------------------------------------------------------

SPLICE:

my @rx;
my @ry;
my @x;
my @y;
my @temp_y;

open RXDATA, '<', "temp.rough-xdata.dat";
$i = 0;
while(<RXDATA>){
   $rx[$i] = $_;
   $i++;
}
close RXDATA;
open RYDATA, '<', "temp.rough-ydata.dat";
$i = 0;
while(<RYDATA>){
   $ry[$i] = $_;
   $i++;
}
close RYDATA;
open XDATA, '<', "temp.x-data.dat";
$i = 0;
while(<XDATA>){
   $x[$i] = $_;
   $i++;
}
close XDATA;
open YDATA, '<', "temp.y-data.dat";
$i = 0;
while(<YDATA>){
   $y[$i] = $_;
   $i++;
}
close YDATA;


open DATAOUT, '>', "$basename_spectrum.dat";

my $length_smooth = @x;
my $length_rough = @rx;
my $length_print;

for( $i=0; $i<$length_smooth; $i++){
  if( $x[$i] >= $norm_min ){
    if ($x[$i] <= $norm_max ){
        push  @temp_y, $y[$i];
    }
  }
}
my $norm_const = max( @temp_y );

if( $length_smooth > $length_rough ){
   $length_print = $length_smooth;
}
else{
   $length_print = $length_rough;
}

my $min_rough = min( @ry );

foreach $i ( @ry ){
   $i = $i - $min_rough;
}

print DATAOUT "Frequency, Smoothed Intensity, Frequency, Rough Intensity\n";

foreach $j ( @y ){
   $j = $j / $norm_const;
}

foreach $j ( @ry ){
   $j = $j / $norm_const;
}
   
for( $i=50; $i<$length_print; $i++ ){
   printf DATAOUT "%16.5f%16.5f%16.5f%16.5f\n", $x[$i], $y[$i], $rx[$i], $ry[$i];
}

close DATAOUT;

unlink "temp.rough-xdata.dat";
unlink "temp.rough-ydata.dat";
unlink "temp.x-data.dat";
unlink "temp.y-data.dat";

open LOG, '>>', "$output_file";
print LOG "\n";
print LOG "----------------------------------------------------------------------------\n";
print LOG "============================================================================\n";
}
#==============================================================================
#==============================================================================
# SPECTRA CALCULATIONS END HERE
# =============================================================================
# ============================================================================
print LOG "    NORMAL TERMINATION \n";


exit 0;

