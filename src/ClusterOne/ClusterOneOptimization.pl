#!/usr/bin/perl -w
use strict;
use Data::Dumper;

#a function to return the matched protein complexes with CORUM and GO complexes database.
#they take two hashes, and mapped the second hash to the first hash.
#overlap(hash1, hash2)---indicates how many complexes in hash2 are mapped to hash1
#cutoff: mached_proteins/min(no_proteins_hash1, no_proteins_hash2) >= 0.5

sub overlap{
    
    my %hash1 = %{shift()};
    my %hash2 = %{shift()};

    my $no_matched_complexes = 0;
    
    foreach my $key (sort keys %hash1) {
        
        my @temp_hash1 = split("\t", $hash1{$key});
        
        my $no_hash1 = scalar @temp_hash1;
        
        foreach my $key (sort keys %hash2) {
            
            my @temp_hash2 = split("\t", $hash2{$key});
            
            if (scalar @temp_hash2 > 2) {
                
            my %array1 = map { $_ => 1 } @temp_hash2;
            my @intersect = grep { $array1{$_} } @temp_hash1;

            my $no_proteins;
            
            if (scalar @temp_hash1 < scalar @temp_hash2) {
                $no_proteins = scalar @temp_hash1;
            } else {
                $no_proteins = scalar @temp_hash2;
            }
            
            if (((scalar @intersect) / $no_proteins) >= 0.5) {

                $no_matched_complexes += 1;
                last;
            } 
            }
        }
}
    return $no_matched_complexes;
}

#Read the CORUM complex database, so the ClusterONE predicted ones can be used to compared to the CORUM complex database...
open(IN, "temporal_output.txt") or die("Can't open IN CORUM files for reading\n");
my %Corum_complexes;

my $no_corum_complex = 0;

while (my $lines = <IN>) {
    chomp $lines;
    $no_corum_complex += 1;
    my @temp = split(/\t/, $lines);
    $Corum_complexes{$no_corum_complex} = join("\t", @temp);
}
close(IN);
#print Dumper(\%Corum_complexes);

#Read the GO complex database, so the ClusterONE predicted ones can be used to compared to the GO complex database...
open(IN, "Go_elegans_complex_experiments_mapped.txt") or die("Can't open IN GO files for reading\n");
my %Go_complexes;

my $no_Go_complex = 0;

while (my $lines = <IN>) {
    chomp $lines;
    $no_Go_complex += 1;
    my @temp = split(/\t/, $lines);
    $Go_complexes{$no_Go_complex} = join("\t", @temp);
}
close(IN);
#print Dumper(\%Go_complexes);

#optimize the parameters by two for loops for trying different sizes and densities...
for my $size (2 .. 3) {
    for my $density_interger (1 .. 9) {
        my $density = $density_interger/10;
        system("java -jar cluster_one-1.0.jar TESTFILE_final_PPI_ratio5.txt -s ".$size." -d ".$density." -F plain > temporal_output.txt");
        
        #read the predicted complexes file and store it into a hash...
        open(IN, "temporal_output.txt") or die("Can't open IN files for reading\n");
        my $no_complex = 0;
        my %temp_complex;
        while (my $lines = <IN>) {
            chomp $lines;
            $no_complex += 1;
            my @temp = split(/\t/, $lines);
            $temp_complex{$no_complex} = join("\t", @temp);
        }
        close(IN);
        
         
        my $overlap_predicted_with_corum = overlap(\%temp_complex, \%Corum_complexes);
        my $overlap_corum_with_predicted = overlap(\%Corum_complexes, \%temp_complex);
        
        my $overlap_predicted_with_GO = overlap(\%temp_complex, \%Go_complexes);
        my $overlap_GO_with_predicted = overlap(\%Go_complexes, \%temp_complex);
        
        print ("the number of overlap_predicted_with_corum is ", $overlap_predicted_with_corum, "\n");
        print ("the number of overlap_corum_with_predicted is ", $overlap_corum_with_predicted, "\n");
        print ("the number of overlap_predicted_with_GO is ", $overlap_predicted_with_GO, "\n");
        print ("the number of overlap_GO_with_predicted is ", $overlap_GO_with_predicted, "\n");


        
    }
}























