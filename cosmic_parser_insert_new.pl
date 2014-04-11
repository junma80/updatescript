use strict;
use List::MoreUtils qw(firstidx);
use List::MoreUtils qw(uniq);
#my $mutantfile="mutantsample.txt";

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

my $mutantfilezip="$datestring.CosmicMutantExport.tsv.gz";
my $indelfilezip="$datestring.CosmicInsMutExport.tsv.gz";
my $vcffilezip="$datestring.CosmicCodingMuts.vcf.gz";

my $mutantfile="$datestring.CosmicMutantExport.tsv";
my $indelfile="$datestring.CosmicInsMutExport.tsv";
my $vcffile="$datestring.CosmicCodingMuts.vcf";

unless (-e "$outfolder$mutantfile"){
	`wget ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicMutantExport_v68.tsv.gz -O $mutantfilezip`;
	`mv $mutantfilezip $outfolder`;
        `gzip -d $outfolder$mutantfilezip`;
}
unless (-e "$outfolder$indelfile"){
	`wget ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicInsMutExport_v68.tsv.gz -O $indelfilezip`;
	`mv $indelfilezip $outfolder`;
	`gzip -d $outfolder$indelfilezip`;
}
unless (-e "$outfolder$vcffile"){
	`wget ftp://ngs.sanger.ac.uk/production/cosmic/CosmicCodingMuts_v68.vcf.gz -O $vcffilezip`;
	`mv $vcffilezip $outfolder`;
	`gzip -d $outfolder$vcffilezip`;
}

#my $mutantfile="CosmicMutantExport_v67_241013.tsv";
#my $indelfile="CosmicInsMutExport_v67_241013.tsv";
#my $vcffile="CosmicCodingMuts_v67_20131024.vcf";

my $outfile="$datestring.COSMIC_insert_new.sql";

open(VCF,"<$outfolder$vcffile");
open(COSMUT,"<$outfolder$mutantfile");
open(INDEL,"<$outfolder$indelfile");
open(OUTFILE,">$outfolder$outfile");
open(OUTRENAME,">$outfolder$datestring.COSMIC_rename.sql");


print OUTFILE "DROP TABLE IF EXISTS COSMIC_update;\n"; 
print OUTFILE "CREATE table COSMIC_update LIKE COSMIC;\n";
print OUTFILE "DROP INDEX idx_name on COSMIC_update;\n";
my %vcfhash;
while (my $line=<VCF>){
	chomp($line);
	if ($line !~ /^#/){
		my @cols=split(/\t/,$line);
		my $chr=$cols[0];
		my $pos=$cols[1];
		my $ref=$cols[3];
		my $var=$cols[4];
		my $info=$cols[7];
		my @infocols=split(/;/,$info);
		my $genename =$infocols[0];
		$genename =~ s/GENE=//g;
		#my @genenamecols= split(/_/,$genename);
		#$genename=$genenamecols[0];
		my $AA =$infocols[3];
		$AA =~ s/AA=//g;
		my $count =$infocols[4];
		$count =~ s/CNT=//g;
		my $strand =$infocols[1];
		$strand =~ s/STRAND=//g;
		my $cds = $infocols[2];
		$cds =~ s/CDS=//g;
		if ($strand =~ /-/){
			$ref =~ tr/ATGC/TACG/;
			$var =~ tr/ATGC/TACG/;
		}
		my $mutantstring="\'$chr\',$pos,\'$cds\'";
		if (length($ref)<100 && length($var)<100) {
			$vcfhash{$mutantstring}="\'$ref\',\'$var\',\'$genename\',\'$AA\',$count";
		}
	}
}	

#my $size = keys %vcfhash;
#print "$size\n";


my $header = <COSMUT>;
chomp($header);
my @headercols=split(/\t/,$header);
my $mutantIDpos= firstidx {$_ eq 'Mutation ID'} @headercols;
my $sampleIDpos= firstidx {$_ eq 'ID_tumour'} @headercols;
my $genenamepos= firstidx {$_ eq 'Gene name'} @headercols;
my $acnumberpos= firstidx {$_ eq 'Accession Number'} @headercols;
my $primarysitepos= firstidx {$_ eq 'Primary site'} @headercols;
my $sitesubtypepos= firstidx {$_ eq 'Site subtype'} @headercols;
my $primaryhistpos= firstidx {$_ eq 'Primary histology'} @headercols;
my $histsubtypepos= firstidx {$_ eq 'Histology subtype'} @headercols;
my $cdspos=firstidx {$_ eq 'Mutation CDS'} @headercols;
my $aapos=firstidx {$_ eq 'Mutation AA'} @headercols;
my $mutdesppos= firstidx {$_ eq 'Mutation Description'} @headercols;
my $zygositypos= firstidx {$_ eq 'Mutation zygosity'} @headercols;
my $positionpos= firstidx {$_ eq 'Mutation GRCh37 genome position'} @headercols;
my $pubmedpos= firstidx {$_ eq 'Pubmed_PMID'} @headercols;
my $somaticstatpos= firstidx {$_ eq 'Mutation somatic status'} @headercols;

my %mutanthash;
my %mutantdetailhash;
while (my $line = <COSMUT>){
	chomp($line);
	my @cols=split(/\t/,$line);
	my $mutantID=$cols[$mutantIDpos];
	my $sampleID=$cols[$sampleIDpos];
	my $position=$cols[$positionpos];
	my @positioncols=split(/:/,$position);
	my $chr=$positioncols[0];
	my @positioncols1=split(/-/,$positioncols[1]);
	my $pos=$positioncols1[0];
	my $cds=$cols[$cdspos];
	$cds =~ /c.(\d+)(\w+)>(\w+)/;
	my $ref=$2;
	my $var=$3;
	my $genename=$cols[$genenamepos];
	my $acnumber=$cols[$acnumberpos];
	my $aa=$cols[$aapos];
	my $mutdesp=$cols[$mutdesppos];
	my $zygosity=$cols[$zygositypos];
	my $pubmed=$cols[$pubmedpos];
	my $somaticstat=$cols[$somaticstatpos];

	my $primarysite=$cols[$primarysitepos];
	my $sitesubtype=$cols[$sitesubtypepos];
	my $primaryhist=$cols[$primaryhistpos];
	my $histsubtype=$cols[$histsubtypepos];

	#my $samplestring="$pubmed:$primarysite:$sitesubtype:$primaryhist:$histsubtype";
	#my $samplestring ="$primarysite;$sitesubtype&$primaryhist;$histsubtype";
	my $samplestring="$primarysite;$sitesubtype";
	my $mutantdetail= "\'$acnumber\',\'$mutdesp\',\'$zygosity\',";
	if ($chr ne "" && $ref ne ""){
		my $mutantstring="\'$chr\',$pos,\'$cds\'";
		#$mutanthash{$mutantstring}{$sampleID}=$samplestring;
		$mutantdetailhash{$mutantstring}=$mutantdetail;
		push (@{$mutanthash{$mutantstring}},$samplestring);
	}
}

#$size = keys %mutanthash;
#print "$size\n";

$header =<INDEL>;
chomp($header);
@headercols=split(/\t/,$header);
$mutantIDpos= firstidx {$_ eq 'Mutation ID'} @headercols;
$genenamepos= firstidx {$_ eq 'Gene name'} @headercols;
$acnumberpos= firstidx {$_ eq 'Accession Number'} @headercols;
$primarysitepos= firstidx {$_ eq 'Primary site'} @headercols;
$sitesubtypepos= firstidx {$_ eq 'Site subtype'} @headercols;
$primaryhistpos= firstidx {$_ eq 'Primary histology'} @headercols;
$histsubtypepos= firstidx {$_ eq 'Histology subtype'} @headercols;
$cdspos=firstidx {$_ eq 'Mutation CDS'} @headercols;
$aapos=firstidx {$_ eq 'Mutation AA'} @headercols;
$mutdesppos= firstidx {$_ eq 'Mutation Description'} @headercols;
$zygositypos= firstidx {$_ eq 'Mutation zygosity'} @headercols;
$positionpos= firstidx {$_ eq 'GRCh37 genome position'} @headercols;
$pubmedpos= firstidx {$_ eq 'Pubmed_PMID'} @headercols;



my %indelhash;
my %indeldetailhash;
while (my $line = <INDEL>){
	chomp($line);
	my @cols=split(/\t/,$line);
	my $mutantID=$cols[$mutantIDpos];
	my $position=$cols[$positionpos];
	my @positioncols=split(/:/,$position);
	my $chr=$positioncols[0];
	my @positioncols1=split(/-/,$positioncols[1]);
	my $pos=$positioncols1[0];
	my $cds=$cols[$cdspos];
	$cds =~ /c.(\d+)(\w+)>(\w+)/;
	my $ref=$2;
	my $var=$3;
	my $genename=$cols[$genenamepos];
	my $acnumber=$cols[$acnumberpos];
	my $aa=$cols[$aapos];
	my $mutdesp=$cols[$mutdesppos];
	my $zygosity=$cols[$zygositypos];
	my $pubmed=$cols[$pubmedpos];
	my $primarysite=$cols[$primarysitepos];
	my $sitesubtype=$cols[$sitesubtypepos];
	my $primaryhist=$cols[$primaryhistpos];
	my $histsubtype=$cols[$histsubtypepos];

	#my $samplestring="$pubmed:$primarysite:$sitesubtype:$primaryhist:$histsubtype";
	#my $samplestring ="$primarysite:$sitesubtype&$primaryhist;$histsubtype";
	my $samplestring ="$primarysite;$sitesubtype";
	my $mutantdetail= "\'$acnumber\',\'$mutdesp\',\'$zygosity\',";
	if ($chr ne ""){
		my $mutantstring="\'$chr\',$pos,\'$cds\'";
		#$mutanthash{$mutantstring}{$sampleID}=$samplestring;
		$indeldetailhash{$mutantstring}=$mutantdetail;
		push (@{$indelhash{$mutantstring}},$samplestring);
	}
}

#$size = keys %indelhash;
#print "$size\n";


#my @uniqsamples=uniq @samplearray;
#my $samplesize=@uniqsamples;
my $count=0;
my $newline;
my @linearray=();
foreach my $mutantID (keys %vcfhash){
	if (exists $mutanthash{$mutantID}){
		my $line=$vcfhash{$mutantID};
		$line.=",".$mutantdetailhash{$mutantID}."\'";
		#my $count=0;
		my @records=@{$mutanthash{$mutantID}};
		my @uniquerecords= uniq @records;
		foreach my $record (@uniquerecords){
			#$line.="$record|";
			$line.="$record<br>";	
			#$count=$count+1;
		}
		#chop($line);
		$line =~ s/<br>$//;
		$line.="\'";
		#if ($count > 2){
		#	print "$mutantID,$line;\n";
		#}
		$count=$count+1;
		push (@linearray,"($mutantID,$line)");
		#print OUTFILE "INSERT INTO COSMIC (chr,pos,cds,ref,var,genename,aachanges,numcnt,acnumber,mutdescription,muttype,tumortype) VALUES ($mutantID,$line);\n";
	}elsif (exists $indelhash{$mutantID}){
		my $line=$vcfhash{$mutantID};
		$line.=",".$indeldetailhash{$mutantID}."\'";
		#my $count=0;
		my @records=@{$indelhash{$mutantID}};
		my @uniquerecords= uniq @records;
		foreach my $record (@uniquerecords){
			#$line.="$record|";
			$line.="$record<br>";	
			#$count=$count+1;
		}
		$line =~ s/<br>$//;
		#chop($line);
		$line.="\'";
		#if ($count > 2){
		#	print "$mutantID,$line;\n";
		#}
		$count=$count+1;
		push (@linearray,"($mutantID,$line)");
		#print OUTFILE "INSERT INTO COSMIC (chr,pos,cds,ref,var,genename,aachanges,numcnt,acnumber,mutdescription,muttype,tumortype) VALUES ($mutantID,$line);\n";
	}else{
		$count=$count+1;
		push (@linearray,"($mutantID,$vcfhash{$mutantID},'','','','')");
		#print OUTFILE "INSERT INTO COSMIC (chr,pos,cds,ref,var,genename,aachanges,numcnt,acnumber,mutdescription,muttype,tumortype) VALUES ($mutantID,$vcfhash{$mutantID});\n";
	}
	if ($count > 10000) {
		$newline="INSERT INTO COSMIC_update (chr,pos,cds,ref,var,genename,aachanges,numcnt,acnumber,mutdescription,muttype,tumortype) VALUES ";
		foreach my $ele (@linearray){
			$newline.="$ele,";
		}
		$newline =~ s/,$/;\n/;
		print OUTFILE "$newline";
		@linearray=();
		$count=0;
	}
	
}
$newline="INSERT INTO COSMIC_update (chr,pos,cds,ref,var,genename,aachanges,numcnt,acnumber,mutdescription,muttype,tumortype) VALUES ";
foreach my $ele (@linearray){
	$newline.="$ele,";
}
$newline =~ s/,$/;\n/;
print OUTFILE "$newline";
print OUTFILE "ALTER IGNORE TABLE COSMIC_update add unique index idx_name(chr,pos,ref,var,genename);";

my $md5sum = `md5sum $outfolder$vcffile`;
my $file = `ls -lah $outfolder$vcffile`;
chomp($file);
chomp($md5sum);

print OUTRENAME "RENAME TABLE COSMIC to COSMIC_$datehash{$datecols[1]}_$year,COSMIC_update to COSMIC;\n";
print OUTRENAME "UPDATE METADATA SET tableName='COSMIC_$datehash{$datecols[1]}_$year' where tableName='COSMIC';\n";
print OUTRENAME "INSERT INTO METADATA VALUES ('COSMIC','$file','$md5sum','$date');\n";
close OUTFILE;
close OUTRENAME;

