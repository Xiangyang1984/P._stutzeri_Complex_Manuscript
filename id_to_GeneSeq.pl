#!usr/bin/perl -w

use strict;
use warnings;


# perl /media/xiangyang/X128/AMRG_TT.v-1.02/extract_gene_seq/id_to_GeneSeq.pl /media/xiangyang/X128/AMRG_TT.v-1.02/extract_gene_seq/proteins_fold /media/xiangyang/X128/AMRG_TT.v-1.02/extract_gene_seq/id.out

my $protein_sequence_fold = $ARGV[0];   # protein sequences for each genome under a fold
my $id_list = $ARGV[1]; # a gene id list file;

&extract_sequence($protein_sequence_fold, $id_list);

sub extract_sequence {

    my ($protein_folder, $parsed_list) = @_;

    opendir (PROTEIN_FOLDER, $protein_folder);
    my @total_protein = readdir PROTEIN_FOLDER;
    closedir PROTEIN_FOLDER;

    chdir $protein_folder;

    my %seq_hash;

    foreach my $eachfile (@total_protein) {
        my $term     = $/; # input record separator;
        $/ = ">";
        open (EACHFILE, $eachfile) or die "could not open infile $eachfile\n"; 

        while(<EACHFILE>) {
            chomp;
            next if($_ eq '');
            my ($id, @sequencelines) = split /\n/;
            foreach my $line (@sequencelines) {
                $seq_hash{$id} .= $line;
            }
        }
        
        $/ = $term;
    }
 
    my $folder_path = go_up_folder($protein_folder);

    my $out_protein_file = "$folder_path/gene_sequence.fasta";
    open (OUT, ">$out_protein_file");
 
    open (LIST, "$parsed_list") or die "can not open $parsed_list";

   
        while (<LIST>){
            chomp;
            print OUT ">$_\n$seq_hash{$_}\n";

        }

        close OUT;


#return $out_protein_file;
}



sub go_up_folder {

    my $current_path = shift;
    $current_path =~ s/(.*)\///; 
    my $up_folder = $1;
    return $up_folder;
}
