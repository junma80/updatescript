use strict;
use List::MoreUtils qw(firstidx);

my $outfolder="/space1/databaseUpdates/";
my $date =`date +"%Y-%m-%d %k:%M:%S"`;
chomp($date);
my @datecols=split(/-/,$date);
my %datehash=("01" => "12","02" => "01","03" => "02", "04" => "03", "05" => "04", "06" => "05", "07" => "06", "08" => "07", "09" => "08", "10" => "09", "11" => "10", "12" =>"11");
my %dateyearhash=("2015" => "2014", "2016" => "2015", "2017" => "2016");
my $datestring="$datecols[1]_$datecols[0]";
my $year=$datecols[0];
if ($datecols[1] eq "01"){
	$year=$dateyearhash{$datecols[0]};
}


my $snpfile="HGMD_2013_SNP_V4.txt";
my $indelfile="HGMD_2013_Indel_V4.txt";
my $outfile="$datestring.hgmd.insert_new.sql";

unless (-e "$outfolder$snpfile"){
	`cp $snpfile $outfolder`;
}
unless (-e "$outfolder$indelfile"){
	`cp $indelfile $outfolder`;
}

open (SNP,"<$outfolder$snpfile");
open (INDEL,"<$outfolder$indelfile");
open (OUT,">$outfolder$outfile");

print OUT "DROP TABLE IF EXISTS HGMD_update;\n";
print OUT "CREATE TABLE HGMD_update LIKE HGMD;\n";

my $header = <SNP>;
chomp($header);
my @headercols=split(/\t/,$header);

my $typepos=firstidx {$_ eq 'Variant_class'} @headercols;
my $positionpos=firstidx {$_ eq 'genomic_coordinates_hg19'} @headercols;
my $genepos=firstidx {$_ eq 'gene'} @headercols;
my $diseasepos=firstidx {$_ eq 'disease'} @headercols;
my $sequencepos=firstidx {$_ eq 'sequence_context_hg19'} @headercols;
my $pmidpos=firstidx {$_ eq 'pmid'} @headercols;
my $hgmdIDpos=firstidx {$_ eq 'ACC_NUM'} @headercols;
my $dnapos=firstidx {$_ eq 'HGVS_cdna'} @headercols;


my $count=0;
my @linearray=();
my $newline;
while (my $line =<SNP>){
	chomp($line);
	my @cols=split(/\t/,$line);
	my $hgmdID=$cols[$hgmdIDpos];
	my $gene=$cols[$genepos];
	my $disease=$cols[$diseasepos];
	$disease =~ s/^"//;
	$disease =~ s/"$//;
	$disease =~ s/\'/\\'/g;
	my $pmid=$cols[$pmidpos];
	my $type=$cols[$typepos];
	my $position=$cols[$positionpos];
	my $refvar;
	my $ref;
	my $var;
	my $chr;
	my $pos;
	if ($position ne "null"){
		my @positioncols=split(/:/,$position);
		$chr=$positioncols[0];
		$pos=$positioncols[1];
		my $strand=$positioncols[2];
		$chr =~ s/chr//;
		my $dna=$cols[$dnapos];
		if ($dna ne "null"){
			my @dnacols=split(/:/,$dna);
			$refvar=$dnacols[1];
			$refvar =~ /c.(\d+)(\w+)-(\w+)/;
			$ref=$2;
			$var=$3;
		
		}else{
			my $sequence=$cols[$sequencepos];
			$sequence =~ /(\w+)\[(\w+)\/(\w+)\](\w+)/;
			$ref=$2;
			$var=$3;
		}
		if ($strand =~ /-/){
			$ref =~ tr/ATGC/TACG/;
			$var =~ tr/ATGC/TACG/;
		}
		$count=$count+1;
		push (@linearray,"(\'$chr\',$pos,$pos,\'$gene\',\'$disease\',\'$type\',\'$ref\',\'$var\',\'$pmid\',\'$hgmdID\')");
		#print "$chr\t$pos\t$gene\t$disease\t$type\t$ref\t$var\t$pmid\t$hgmdID\n";
		if ($count >10000){
			$newline="INSERT INTO HGMD_update VALUES ";
			foreach my $ele (@linearray){
				$newline.="$ele,";
			}
			$newline =~ s/,$/;\n/;
			print OUT "$newline";
			@linearray=();
			$count=0;
		}
	}
}
$newline="INSERT INTO HGMD_update VALUES ";
foreach my $ele (@linearray){
	$newline.="$ele,";
}
$newline =~ s/,$/;\n/;
print OUT "$newline";
@linearray=();
$count=0;


$header=<INDEL>;
chomp($header);
@headercols=split("\t",$header);

$typepos=firstidx {$_ eq 'tag'} @headercols;
$positionpos=firstidx {$_ eq 'genomic_coords_hg19'} @headercols;
$genepos=firstidx {$_ eq 'gene'} @headercols;
$diseasepos=firstidx {$_ eq 'disease'} @headercols;
$sequencepos=firstidx {$_ eq 'sequence_context_hg19'} @headercols;
$pmidpos=firstidx {$_ eq 'pmid'} @headercols;
$hgmdIDpos=firstidx {$_ eq 'acc_num'} @headercols;
$dnapos=firstidx {$_ eq 'hgvs'} @headercols;

my $outvcf ="$outfolder$datestring.hgmd_indel.vcf";
open(OUTVCF,">$outvcf");
print OUTVCF  "##fileformat=VCFv4.1
##fileDate=20120810
##source=GenerateReportDataAndVCFv2.0.0.0
##reference=HumanNCBI37_UCSC_XX
##phasing=none
##FILTER=<ID=q20,Description=\"Quality below 20\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Reference Read Depth\">
##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Major Variant Read Depth\">
##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO\n";
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO	GT:VR:RR:DP:GQ	clia_MAXGT\n";


while (my $line =<INDEL>){
	chomp($line);
	my @cols=split(/\t/,$line);
	my $hgmdID=$cols[$hgmdIDpos];
	my $gene=$cols[$genepos];
	my $disease=$cols[$diseasepos];
	$disease =~ s/\'/\\'/g;
	$disease =~ s/^"//;
	$disease =~ s/"$//;	
	my $pmid=$cols[$pmidpos];
	my $type=$cols[$typepos];
	my $position=$cols[$positionpos];
	my $refvar;
	my $ref;
	my $var;
	my $chr;
	my $pos;
	my $origpos;
	if ($position ne "null"){
		my @positioncols=split(/:/,$position);
		$chr=$positioncols[0];
		my @pos1cols=split(/\s/,$positioncols[1]);
		my @pos2cols=split(/-/,$pos1cols[0]);
		$origpos=$pos2cols[0];
		my $strand=$pos1cols[1];
		$chr =~ s/chr//;
		$pos=$origpos;
		my $dna=$cols[$dnapos];
		if ($dna ne "null"){
			my @dnacols=split(/:/,$dna);
			$refvar=$dnacols[1];
			my $sequence=$cols[$sequencepos];
			$sequence =~ /(\w+)\[((\w+|-))\/((-|\w+))\](\w+)/;
			my $var1=$2;
			my $var2=$4;
			#print "$refvar\t$1\t$2\t$4\t$6\n";
			my $first=$1;
			my $last=$6;
			if ($strand =~ /-/) {
				$ref=substr($last,0,1);
				$var1= reverse $var1;
				$var2= reverse $var2;
			}else{
				$ref=chop($first);
			}
			#print "$refvar\t$var1\t$var2\t$ref\n";
			if ($refvar =~ /del/) {
				$pos=$origpos-1;
				if ($refvar =~ /ins/) {
					$var=$ref.$var2;
					$ref=$ref.$var1;
				}else{				
					$var=$ref;
					$ref=$ref.$var1;
				}
			}elsif ($refvar =~ /dup/){
				$var=$ref.$var2;
			}elsif ($refvar =~ /ins/){
				$var=$ref.$var2;				
			}
		}	
		if ($strand =~ /-/){
			$ref =~ tr/ATGC/TACG/;
			$var =~ tr/ATGC/TACG/;
		}
		$count=$count+1;
		print OUTVCF "$chr\t$pos\t.\t$ref\t$var\t60\tPASS\tAA=$chr|$origpos|$gene|$disease|$type|$pmid|$hgmdID;\n"
	}
}

my $normvcf="$outfolder$datestring.hgmd_indel.norm.vcf";

`/home/ma/htslib-master/htscmd vcfnorm -f /home/matthew/projects/databases/humanReference/hg19.fa $outvcf > $normvcf`;


open(INFILE, "<$normvcf");

$count=0;
@linearray=();
$newline="";
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO";
while (my $line=<INFILE>){
    chomp($line);
    if ($line !~ /^#/) {
        my @cols=split("\t",$line);
        my $chr=$cols[0];
        my $pos=$cols[1];
        my $ref=$cols[3];
        my $var=$cols[4];
        my $info=$cols[7];
        my @infocols=split(/\|/,$info);
        my $origpos=$infocols[1];
        my $gene=$infocols[2];
        my $disease=$infocols[3];
        my $type=$infocols[4];
        my $pmid=$infocols[5];
        my $hgmdID=$infocols[6];
        $count=$count+1;
	push (@linearray,"(\'$chr\',$pos,$origpos,\'$gene\',\'$disease\',\'$type\',\'$ref\',\'$var\',\'$pmid\',\'$hgmdID\')");
        if ($count >10000){
            $newline="INSERT INTO HGMD_update VALUES ";
            foreach my $ele (@linearray){
                    $newline.="$ele,";
            }
            $newline =~ s/,$/;\n/;
            print OUT "$newline";
            @linearray=();
            $count=0;
        }
    } 
}
$newline="INSERT INTO HGMD_update VALUES ";
foreach my $ele (@linearray){
	$newline.="$ele,";
}
$newline =~ s/,$/;\n/;
print OUT "$newline";
print OUT "RENAME TABLE HGMD to HGMD_$datestring;\n";
print OUT "RENAME TABLE HGMD_update to HGMD;\n";
print OUT "UPDATE METADATA SET tableName='HGMD_$datestring' where tableName='HGMD';\n";
 my $md5sumsnp = `md5sum $outfolder$snpfile`;
 my $md5sumindel= `md5sum $outfolder$indelfile`;
 #$pwd = `pwd`;chomp($pwd);
 my $snpfileinfo = `ls -lah $outfolder$snpfile`;
 my $indelfileinfo= `ls -lah $outfolder$indelfile`;

chomp($md5sumsnp);
chomp($md5sumindel);
chomp($snpfileinfo);
chomp($indelfileinfo);
chomp($date);

print OUT "INSERT INTO METADATA VALUES ('HGMD','$snpfileinfo','$md5sumsnp','$date');\n";
print OUT "INSERT INTO METADATA VALUES ('HGMD','$indelfileinfo','$md5sumindel','$date');\n";


