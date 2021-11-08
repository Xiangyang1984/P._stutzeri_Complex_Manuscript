use strict;
use warnings;


# Usage
############################################################################## 
# perl /media/xiangyang/X128/test/127/Orthogroups/format_modification.pl /media/xiangyang/X128/test/127/Orthogroups/genomic_name.txt /media/xiangyang/X128/test/127/Orthogroups/merged_Orthogroups.txt /media/xiangyang/X128/test/127/Orthogroups/OUT.txt
############################################################################## 

my $strain_list = $ARGV[0];
my $orthogroup = $ARGV[1];
my $ortho_to_PGPA = $ARGV[2];

open ("LIST", $strain_list);
my @strain;
while(<LIST>){
    chomp;
    @strain = split "\t", $_;
}
close LIST;

open ("OUT", ">$ortho_to_PGPA");
print OUT "ClutserID\t", join ("\t", @strain), "\n"; # print all the strain nams 

open ("ORTH", $orthogroup);


my $cluster_number =0;
while(<ORTH>){
    chomp;
    $cluster_number++;
    print OUT "$cluster_number\t";
    my @homologs = split " ", $_;
    foreach my $strain (@strain){
        my @cluster = grep($_ =~ /$strain/, @homologs); #obtaining all homologous genes from each genome
        print OUT join (", ", @cluster) , "\t" if (scalar @cluster) >0;
        print OUT "-\t" if(scalar @cluster) eq 0;
        print OUT "\t";
    }
    print OUT "\n";
}

close ORTH;
close OUT;
