#!/usr/bin/perl -w

use strict;
use warnings;
use threads;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename qw<basename dirname>;
use FindBin;
use File::Spec;

#########################################################################################################################################
# Before using, extract_protein_dir.pl needs to pre-install several Perl Modules.
# Required Perl Modules: threads File::Basename Getopt::Long FindBin File::Spec Bio::SeqIO; (part of BioPerl, http://bioperl.org)
#########################################################################################################################################

my $usage = <<USAGE; 

=NAME

extract_protein_dir.pl

=DESCRIPTION

    Run this command to enble users to extract protein sequences for genomes using multiple threads. 

=USAGE

extract_protein_dir.pl -dir genbank_file_directory [options]

FOR EXAMPLE: 
# perl /media/xiangyang/X128/test/127/Orthogroups/extract_protein_dir.pl -dir /media/xiangyang/X128/stutzeri_123 -m 3

=ARGUMENTS
=======================
    REQUIRED ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -dir, --genbank_file_directory
           A directory containing annotated genomes as Genbank format file. To avoid a mistake, genome names cannot use special character,
           such as space, equal. For large number of genomes, users are recommended to download using Aspera, a high-speed file transfer
           tool (https://downloads.asperasoft.com/).                           
 
    OPTIONAL ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -o, --output_directory
           An output directory holding all the generated files by extract_protein_dir.pl. if this option is not set, extract_protein_dir.pl will create a directory named "interested_gene_workplace" in the bin directory from where extract_protein_dir.pl was invoked.
       -m, --multiple_threads
           set thread number (Default: 1)
       -h, --help
           Show this message.

=AUTHOR

Dr. Xiangyang Li (E-mail: lixiangyang\@fudan.edu.cn), Fudan university; Kaili University; Bacterial Genome Data mining & Bioinformatic Analysis (www.microbialgenomic.com/).

=COPYRIGHT

Copyright 2021, Xiangyang Li. All Rights Reserved.

USAGE


my %options = (
    'genbank_file_directory'                   => undef,   
    'output_directory'                         => undef,  
    'multiple_threads'                         => "1",
    'help'                                     => undef

);

GetOptions(
    'dir|genbank_file_directory=s'             => \$options{genbank_file_directory},    
    'o|output_directory=s'                     => \$options{output_directory},     
    'm|multiple_threads=s'                     => \$options{multiple_threads},     
    'h|help'                                   => \$options{help}

);

if (defined $options{help}) {
    print "$usage\n";
    exit(0);
}

#check for required options
if ( !( defined( $options{genbank_file_directory} ) ) ) {
    print $usage;
    exit(1);
}


my $now_time = localtime;
print "\n$now_time: extract_protein_dir.pl start...\n\n";

my $path_genbank = File::Spec->rel2abs($options{genbank_file_directory});
my $thread_number = $options{multiple_threads};

my $home_directory = $FindBin::Bin;           # obtaining the home directory where  protein sequences of each genome located

#check for extract_protein_dir.pl workplace options
my $workplace;
if ( defined( $options{output_directory} ) ) {
    $workplace = File::Spec->rel2abs($options{output_directory});
    mkdir $workplace;
}else {

    $workplace = "$home_directory/Whole_proteome";
    $workplace =~ s/\/\//\//g;
    mkdir $workplace;
}

     
     batch_genbank_sequence_TFT_extract($path_genbank, $workplace, $thread_number);

sub batch_genbank_sequence_TFT_extract {

    my ($path_genbank, $path_protein, $thread_number)=@_;

    opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";
    my @path_genbank_temp = readdir PATH_GENBANK;
    foreach (@path_genbank_temp){
        if (/[ =]/) {
            my $new_name = $_;
            $new_name =~ s/ /_/g;   #repalce the space with "_"  for all GenBank files
            $new_name =~ s/=/_/g;   #repalce "=" with "_"  for all GenBank files
            $new_name =~ s/_+/_/g;  
            system ("mv '$path_genbank/$_' '$path_genbank/$new_name'");  
        }

    }
    closedir PATH_GENBANK;

    opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";
    my @path_genbank = readdir PATH_GENBANK;
    @path_genbank =grep ($_!~/^\./ ,@path_genbank);  #delete hidden file . .. 
    #@path_genbank =grep ($_!~/^\~/ ,@path_genbank);  #delete temp file 
    closedir PATH_GENBANK;

    my $sub_file_number = scalar @path_genbank;

    my @input;
    my @outfile_1;
 
    foreach my $file_genbank(@path_genbank){ 

        my $input="$path_genbank/$file_genbank"; 
        push (@input,$input);  
        @input = sort @input;  

        my $output_1="$path_protein/$file_genbank.fasta";
        push (@outfile_1,$output_1);
        @outfile_1 = sort @outfile_1;

    }

    my $thread;
    my @threads;
    my $job_number=0;

    while(){ 

        last if ($job_number eq $sub_file_number);  

        while(scalar(threads->list())<$thread_number) {     #set thread number
            $job_number++;    
            my $input_file = $input[$job_number-1];
            my $output_file_1 = $outfile_1[$job_number-1];
            print "Genbank_extraction: $job_number\n";
            $threads[$job_number-1]=threads->new(\&genbank_protein_sequence_extract, $input_file, $output_file_1); 

            last if ($job_number eq $sub_file_number);  
        }

        foreach $thread(threads->list(threads::all)){
            if($thread->is_joinable()) {              
                $thread->join();                     
            }
        }
    }

    foreach my $th (threads->list(threads::all)) {
        $th->join();
    }

}


sub genbank_protein_sequence_extract {
####### Part of script in this subroutine comes from CGView Comparison Tool. Reference: Grant JR, Arantes AS, Stothard P (2012) Comparing thousands of circular genomes using the CGView Comparison Tool. BMC Genomics 13:202.

    my ($input, $output_1)=@_;
    my $file_name = basename($input);
    open(OUTPUT_1, ">$output_1");

        my $in = new Bio::SeqIO( -file => $input, -format => 'genbank' );

        while ( my $seqObject = $in->next_seq() ) {
            my $acc = $seqObject->accession;
            my @features = $seqObject->get_SeqFeatures();
            my $count = 0;       
            
            foreach (@features) {
                my $feat = $_;	        
                unless ($feat->primary_tag =~ /^CDS/) {
                    next;
                }
                    if ( $feat->has_tag('pseudo') ) {

                    }else {$count++;}
                    #print "genome_$count\n";

                    my $locus_tag;
                    if ( $feat->has_tag('locus_tag') ) {
                        ($locus_tag) = $feat->get_tag_values('locus_tag');

                    }else {$locus_tag = "CDS_$acc"."_".$count;}
               
                    my $product;
                    if ( $feat->has_tag('product') ) {
                        ($product) = $feat->get_tag_values('product');
                        $product =~ s/ /_/g;
                    }
 
                    my $protein_seqeunce;
                    if ( $feat->has_tag('translation')) {
                        ($protein_seqeunce) = $feat->get_tag_values('translation');
                      
                            print OUTPUT_1 ">$locus_tag;$product;$file_name\n$protein_seqeunce\n";
                        
                    }
            }
        }

    close OUTPUT_1;

}


