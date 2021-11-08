use strict;
use warnings;
use threads;
use FindBin;
use lib "$FindBin::Bin/lib";
use Getopt::Long;
use File::Basename qw<basename dirname>;


#########################################################################################################################################
# Before using, Concat_Seq.pl needs to pre-install several Perl Modules and NCBI BLAST+ programs
# Required Perl Modules: threads File::Basename Getopt::Long Bio::SeqIO; (part of BioPerl, http://bioperl.org)
#########################################################################################################################################

my $now_time = localtime;
print "\n$now_time : Concat_Seq start...\n\n";

my $usage = <<USAGE; 

=NAME

Concat_Seq.pl

=DESCRIPTION

    Run this command to enble users to concatenate alignment files into a pseudo-DNA fasta file. 

=USAGE

Concat_Seq.pl -dir alignment_file_directory

FOR EXAMPLE: 
# perl /media/xiangyang/X128/test/127/Orthogroups/Concat_Seq.pl -dir /media/xiangyang/X128/test/127/Orthogroups/alignment_dir

=ARGUMENTS
=======================
    REQUIRED ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -dir, --alignment_file_directory
           A directory containing alignment files as FASTA format file. To avoid a mistake, sequence ids cannot use special character,
           such as space, equal.
       -h, --help
           Show this message.

=AUTHOR

Dr. Xiangyang Li (E-mail: lixiangyang\@fudan.edu.cn), Fudan university; Kaili University; Bacterial Genome Data mining & Bioinformatic Analysis (www.microbialgenomic.com/).

=COPYRIGHT

Copyright 2021, Xiangyang Li. All Rights Reserved.

USAGE


my %options = (
    'alignment_file_directory=s'                 => undef,   
    'help'                                     => undef

);

GetOptions(
    'dir|alignment_file_directory=s'             => \$options{alignment_file_directory},    
    'h|help'                                   => \$options{help}

);

if (defined $options{help}) {
    print "$usage\n";
    exit(0);
}

my $home_directory = $FindBin::Bin;           # obtaining the home directory where concatenation sequence file located


    print "\nConcatenate alignment files into a pseudo-DNA fasta file\n";
    my $SNP_dir =  
    my $snp_fasta_out = concatenation_sequence ($options{alignment_file_directory}, $home_directory);

$now_time = localtime;
print "$now_time : Finished!\n";


sub concatenation_sequence {

my ($path, $out_dir) = @_;
my $snp_fasta_out = "$out_dir/concatenation.fasta";
open (SNP_OUT, ">$snp_fasta_out");

my $Equnce_number=0;
my $number_files     = 0;    # count number of files
my @hash_ref_array   = ();   # array with hash references
my %nseq_hash        = ();   # key:infile, val:nseq

#  Read all infiles to count sequences
opendir DIR, $path or die $!;
my @dir = readdir DIR; 
closedir DIR;
foreach my $each_file (@dir) {
    my $infile       = "$path/$each_file";
    my %seq_hash     = parse_fasta($infile, "T"); # key: seqid, value:sequence
    ## Save sequences in array with hash references.
    push(@hash_ref_array, \%seq_hash);

} # Done with file

######################################concat sequences of identy strain from multiple file############   
    ## First, concatenate all sequences from hashes
    my %print_hash = (); # key:label, value:sequence
    foreach my $h_ref (@hash_ref_array) {
        $number_files++;
        foreach my $seqid (sort keys %$h_ref) { 
            $Equnce_number++;                   # %$h_ref 为哈希引用
            $print_hash{$seqid} .= $h_ref->{$seqid};
        }
    }
    foreach my $label (sort keys  %print_hash) {
        $print_hash{$label} =~ s/\S{60}/$&\n/gs; ## replace word of size $lwidth with itself and "\n"
            print SNP_OUT ">$label\n";
            print SNP_OUT $print_hash{$label}, "\n";
    }

    return $snp_fasta_out;
}


sub parse_fasta {

    my ($infile, $unify_id) = @_;
    my $term     = $/; # input record separator;
    my %seq_hash = (); # key:seqid, val:seq
    open my $INFILE, $infile or die "could not open infile '$infile' : $! \n"; 
    $/ = ">"; # input record separator;
    while(<$INFILE>) {

        chomp;

        next if($_ eq '');
        my ($id, @sequencelines) = split /\n/;
        if ($unify_id eq "T"){
            my @ww = split ";", $id;                                 #使每个文件的ID变为一致
            $id = $ww[2];
        }

        foreach my $line (@sequencelines) {
            $seq_hash{$id} .= $line;
        }
    }
    $/ = $term;
    return(%seq_hash);

} # end of parse_fasta

