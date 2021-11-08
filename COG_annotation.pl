#! usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/lib";
use threads;
use File::Basename qw<basename dirname>;
use Bio::SeqIO;

# perl /media/xiangyang/X128/COG_database/COG_annotation.pl cog-20.fa cog-20.cog.csv cog-20.def.tab fun-20.tab merged_cogs.fa Pan_repersnet_seq.fasta COG_out 16


my $prot_file     =  $ARGV[0];			# /media/xiangyang/X128/COG_database/cog-20.fa  test_lst_1
my $cog_file      =  $ARGV[1];			# /media/xiangyang/X128/COG_database/cog-20.cog.csv
my $cog_def_tab   =  $ARGV[2];		        # /media/xiangyang/X128/COG_database/cog-20.def.tab
my $fun_tab       =  $ARGV[3];		        # /media/xiangyang/X128/COG_database/fun-20.tab
my $merged_cogs   =  $ARGV[4];                  # /media/xiangyang/X128/COG_database/merged_cogs.fa
my $inputfile     =  $ARGV[5];                  # /media/xiangyang/X128/COG_database/Pan_repersnet_seq.fasta
my $cog_out       =  $ARGV[6];                  # /media/xiangyang/X128/COG_database/outfile
my $thread_number =  $ARGV[7];
my $bacth_blast   = "T";


my $e_value = 0.00001;
my $blastp        = `which blastp`;
$blastp =~ s/\n//g;

my $makeblastdb   = `which makeblastdb`;
$makeblastdb =~ s/\n//g;
####################################################################################
my %Func_category_ID_des;
open (FUN_TAB, $fun_tab);
while(<FUN_TAB>){
    chomp;
    my ($Func_category_ID, $color, $Func_category_des) = split "\t", $_;
    $Func_category_ID_des{$Func_category_ID} = $Func_category_des;
    #print "$Func_category_ID\t$Func_category_des\n";
}


####################################################################################
my %COG_ID_Func_category_ID;
open (COG_DEF_TAB, $cog_def_tab);
while(<COG_DEF_TAB>){
    chomp;
    my @cog_def_tab = split "\t", $_;
    $COG_ID_Func_category_ID{$cog_def_tab[0]} = $cog_def_tab[1];
    #print "$cog_def_tab[0]\t$cog_def_tab[1]\n";
}


####################################################################################
my %pID_Func_category_ID;
open (COG_FILE, $cog_file);
while(<COG_FILE>){
    chomp;
#print "$_\n";
    #print "$_\n" if $_ =~ /WP_015254381.1/;
    my @cog_file = split ",", $_;
    my $pID = $cog_file[2];
    $pID =~ s/\./_/;
    #print "$cog_file[2]\t$pID\n";
    $pID_Func_category_ID{$pID} = $COG_ID_Func_category_ID{$cog_file[6]};
}


open (MERGED_COGS, ">$merged_cogs");
my %seq_hash      = parse_fasta($prot_file); # key: seqid, value:sequence
foreach (keys %seq_hash){

    chomp;
    my $new_id = $_;
    $new_id =~ s/ .*//g;
    $new_id =~ s/\./_/;

    print MERGED_COGS ">$new_id", "|", $pID_Func_category_ID{$new_id},"\n", $seq_hash{$_},"\n";
    #print MERGED_COGS ">$new_id", "|", $pID_Func_category_ID{$new_id}, "|", $Func_category_ID_des{$pID_Func_category_ID{$new_id}}, "\n", $seq_hash{$_},"\n" if defined  $Func_category_ID_des{$pID_Func_category_ID{$new_id}};
#print "$new_id\n" if !defined ($pID_Func_category_ID{$new_id});
}

close MERGED_COGS;

format_sequence($makeblastdb, $merged_cogs);
#diamond_format_sequence($diamond, $new_db_file);
 
#chdir $path_protein;
my @input;
my @outfile;

my $home_directory = $FindBin::Bin;   # obtaining the home directory where AMRG_Anno.pl located
my $path_protein = "$home_directory/split_diretory";
my $blast_fold = "$home_directory/blast_diretory";
mkdir $blast_fold;

if ($bacth_blast eq "T"){

    mkdir $path_protein;
    &split_to_subfiles ($inputfile, $thread_number, $path_protein);

}

opendir PATH_PROTEIN, $path_protein or die "could not open $path_protein";
my @PATH_PROTEIN = readdir PATH_PROTEIN; 
@PATH_PROTEIN =grep ($_!~/^\./ ,@PATH_PROTEIN);  #delete hidden file . ..
closedir PATH_PROTEIN;

foreach my $file_protein(@PATH_PROTEIN){
 
    my $input="$path_protein/$file_protein"; 
    push (@input,$input);  
    my $output="$blast_fold/$file_protein.blastout";
    push (@outfile,$output);

}

my $sub_file_number = scalar @input;


my @sub_function_parameters = (\@input, $merged_cogs, \@outfile, $e_value, $blastp);
&bacth_blast($sub_file_number, $thread_number, "do_blastP", \@sub_function_parameters);

#my @sub_function_parameters = (\@input, $merged_cogs, \@outfile, $e_value, $diamond);###edit
#&bacth_blast($sub_file_number, $thread_number, "diamond_blastP", \@sub_function_parameters);

#chdir $blast_fold;
system ("cat $blast_fold/* > $blast_fold/all_vs_all.blast");

open(ALL_ALL, "$blast_fold/all_vs_all.blast");
open(COG_OUT, ">$cog_out");
while(<ALL_ALL>){

    chomp;
    my @blast_out_all = split ("\t", $_);
    my $COG_cat_ture = $blast_out_all[1];  #WP_014357680_1|T
    $COG_cat_ture =~ s/.*\|//g;
    my @COG_cat_ture = split (//, $COG_cat_ture);
    foreach (@COG_cat_ture){

        print COG_OUT "$blast_out_all[0]\t$COG_cat_ture\t$_\n";

    }
}


sub parse_fasta {

    my $infile = shift;
    my $term     = $/; # input record separator;
    my %seq_hash = (); # key:seqid, val:seq
    open my $INFILE, $infile or die "could not open infile '$infile' : $! \n"; 
    $/ = ">"; # input record separator;
    while(<$INFILE>) {

        chomp;

        next if ($_ eq '');
        my ($id, @sequencelines) = split /\n/;

        foreach my $line (@sequencelines) {
            $seq_hash{$id} .= $line;
        }
    }
    $/ = $term;
    return(%seq_hash);

} # end of parse_fasta


sub splitfile{
    my ($ref_input_line, $everyfile_num) = @_;
    
    my $total_line_num = @$ref_input_line; 

    unlink "temp_lst";
    open(TEMPA, ">> temp_lst") || die ("Could not open file temp_lst! \n");

    my $sq_num = 0;
    my $line_index;
    foreach my $line (@$ref_input_line){
        $line_index ++;

        printf TEMPA ("$line");

        if($line_index%$everyfile_num == 0){   #取余运算或取模？？？ #Set partition lines
            $sq_num = $line_index/$everyfile_num;
            my $lst_name = "test_"."lst_".$sq_num;
            rename "temp_lst",  "$lst_name";
            close TEMPA;

            unlink "temp_lst";
            open(TEMPA, ">> temp_lst") || die ("Could not open file temp_lst! \n");
        }
    }

    rename "temp_lst",  "test_"."lst_".($sq_num+1);
    close TEMPA;

    return $sq_num+1;
} #end of splitfile


sub split_to_subfiles{

    my ($file,$split,$output_dir,$log_file);
    open(my $out_fh, '>-') or die "Could not open stdout for writing";

    if ( @_ && $_[0] eq __PACKAGE__)
    {
        GetOptions('i|input-file=s' => \$file,
                   'o|output-dir=s' => \$output_dir,
                   'n|split-number=i' => \$split,
                   'l|log-file=s' => \$log_file) or die "Invalid options\n";
    }
    else
    {
        ($file,$split,$output_dir,$log_file) = @_;
    }

       
    
    #this will find records with a quickness
    (my $total) = `grep -c ">" $file`;
    chomp $total;

    die "no fasta records identified, exiting.\n" unless $total;
    
    my $records_per_file = int ($total / $split);
    my $leftover_records = $total % $split;
    
    if ($split > $total) {
       warn "Total records in file ($total) less than splitnum ($split); adjusting toatal\n";
       $records_per_file = 1;
    }
    
    my $input_file_name = basename($file);
    my $output_base_path = (defined $output_dir) ? "$output_dir/$input_file_name" : $input_file_name;
    my $in = new Bio::SeqIO (-file=>$file, -format=>"fasta");
    my $x;
    my @outs;
    my $records;
    for my $x (1..$split) { 
      #adjust # of records per file for leftover records
      my $adjusted_records_per_file; 
      $adjusted_records_per_file = $records_per_file+1 if $x <= $leftover_records;
      push @outs, new Bio::SeqIO (-file=>">$output_base_path.$x", -format=>"fasta")
    }
    
    my $out = shift @outs;
    my $filecounter = 1;
    my $recordcounter =1;
    while (my $seq = $in->next_seq) {

      my $adjusted_records_per_file = $filecounter<=$leftover_records?$records_per_file+1:$records_per_file;
      $out->write_seq($seq);
      if (++$records>=$adjusted_records_per_file) {
        $out = shift @outs; $records =0; $filecounter++;   
      }
    }
}


#########################################################
#########################################################
sub format_sequence {
    my $makeblastdb = shift;
    my $db_file = shift;
     
    system ("$makeblastdb -in $db_file -dbtype prot -logfile $db_file.makeblast.log");
}



#########################################################
#########################################################
sub do_blastP {
    my ($input, $db_file, $outfile, $e_value, $blastp) =@_;
    system("$blastp -outfmt '6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore' -evalue $e_value -max_hsps 1 -max_target_seqs 1 -num_threads 1 -query $input -db $db_file -out $outfile");  #-max_hsps -num_descriptions -num_alignments

}


sub diamond_format_sequence {
    my $diamond = shift;
    my $db_file = shift;

    system ("$diamond makedb --in $db_file -d $db_file --quiet"); #--sensitive或--more-senstive提高敏感度
}


sub diamond_blastP {
    my ($input, $db_file, $outfile, $e_value, $diamond) =@_;
    system("$diamond blastp --outfmt 6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore --evalue $e_value --max-target-seqs 1 -q $input -d $db_file -o $outfile --quiet --more-sensitive"); #--more-sensitive

}



#########################################################
#########################################################
sub bacth_blast {

    my ($work_number, $thread_number, $subfunction, $sub_function_parameters) = @_;
    my ($input_file, $db_file, $output_file, $e_value, $blastp) = @$sub_function_parameters;
 
    my @input = @$input_file;
    my @outfile = @$output_file;
    my $thread;
    my @threads;
    my $job_number=0;

    while(){ 
        last if ($job_number>=$work_number);                         
        while(scalar(threads->list())<$thread_number) {     #set threadnumber；
            $job_number++;                                 
            my $input_file = $input[$job_number-1];
            my $output_file = $outfile[$job_number-1];
            print "Blast_perform: $job_number\n";
            $threads[$job_number-1]=threads->new(\&$subfunction, $input_file, $db_file, $output_file, $e_value, $blastp);   
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

