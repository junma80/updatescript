# This script takes the omim.txt and genemap2.txt as input to map omimID with genes and add clinical features for diseases. The omimID/genes is mainly  from omim title and disorders column in genemap2.txt. 
#wget ftp://anonymous:junm%40bcm.edu@ftp.omim.org/OMIM/genemap2.txt -p /space1/databasepdates 
#wget ftp://anonymous:junm%40bcm.edu@ftp.omim.org/OMIM/omim.txt.Z -p /space1/databaseupdates
#`gzip -d omim.txt.Z`

my $outfolder="/space1/databaseUpdates/";

my $date =`date +"%Y-%m-%d %k:%M:%S"`;
chomp($date);
@datecols=split(/-/,$date);
%datehash=("01" => "12","02" => "01","03" => "02", "04" => "03", "05" => "04", "06" => "05", "07" => "06", "08" => "07", "09" => "08", "10" => "09", "11" => "10", "12" =>"11");
%dateyearhash=("2015" => "2014", "2016" => "2015", "2017" => "2016");
$datestring="$datecols[1]_$datecols[0]";
$year=$datecols[0];
if ($datecols[1] eq "01"){
	$year=$dateyearhash{$datecols[0]};
}


my $omimtext="omim.$datestring.txt";
my $omimzip="omim.$datestring.txt.Z";
my $genemap="omim.genemap2.$datestring.txt";

unless (-e "$outfolder$genemap"){
	`wget ftp://anonymous:junm%40bcm.edu\@ftp.omim.org/OMIM/genemap2.txt -O $genemap`;
	`mv $genemap $outfolder`;
}
unless (-e "$outfolder$omimtext"){
	`wget ftp://anonymous:junm%40bcm.edu\@ftp.omim.org/OMIM/omim.txt.Z -O $omimzip`;
	`mv $omimzip $outfolder`;
	`gzip -d $outfolder$omimzip`;
}

open(OMIM, "<$outfolder$omimtext");
open(GENEMAP, "<$outfolder$genemap");
open(OUT, ">$outfolder$datestring.omim.insert_new.sql");
open(OUTRENAME, ">$outfolder$datestring.omim.rename.sql");

my @tempout=();
my @updateout=();

#write SQL to update METADAT table

 $tempTable = "OMIM_$datehash{$datecols[1]}_$year";
 $destTable="OMIM";
 $updateTable="OMIM_update";
#print SQLOUT "TRUNCATE $destTable;\n";
 $sourcefile = "$outfolder$omimtext";
#print SQLOUT "DELETE FROM METADATA where tableName='$destTable';\n";
 $md5sum = `md5sum $sourcefile`;
 #$pwd = `pwd`;chomp($pwd);
 $file = `ls -lah $sourcefile`;
 
 $md5summap = `md5sum $outfolder$genemap`;
 $mapfile=`ls -lah $outfolder$genemap`;

chomp($md5sum);
chomp($file);
chomp($date);
chomp($md5summap);
chomp($mapfile);
#print SQLOUT "INSERT INTO METADATA VALUES ('$destTable','$file','$md5sum','$date');\n";
#
print OUT "DROP TABLE IF EXISTS $updateTable;\n"; 
print OUT "CREATE table $updateTable LIKE $destTable;\n";
#print OUT "INSERT $tempTable SELECT * FROM $destTable;\n";
#print OUT "TRUNCATE $destTable;\n";

#creast hash for genemap(omimID/genes)
my %idhash;
while (my $line = <GENEMAP>){
	chomp($line);
        @A = split(/\|/,$line);
	@B = split(/[ ,]/,$A[5]); #gene symbols
	$gene=uc($B[0]);
	$gene =~s/;//g;
        @omimdisorders = ();
	#add in relevent omimid
	push (@omimdisorders,$A[8]);
	#add in omimid in disorder
        while($A[11] =~ m/ [0-9]{6} /g){
       		push(@omimdisorders, $&);
    	}
        foreach my $d (@omimdisorders){
		$d =~ s/^\s+|\s+$//g;
		push(@{$idhash{$d}},$gene);
       	}	
}
foreach my $key  (keys %idhash){
       @{$idhash{$key}}=uniq(@{$idhash{$key}});
}


#write insert sql for OMIM table
my $cnt=0;
#Process omim file
#push(@updateout, "LOCK TABLES OMIM WRITE;\n");
while($line = <OMIM>){
    chomp($line);
    if ($line eq "*RECORD*"){
	$cnt++;
	if ($cnt >1){
		procRec(@record);
		@record=();
	}
    }
    push(@record,$line);
}
procRec(@record);
#print UPDATEOUT "DROP TABLE omim_remove;";
$newline = "INSERT INTO $updateTable VALUES ";
$count=0;
foreach my $line(@tempout){
	$count+=1;
	if ($count <=5000){
		$newline.=$line;
	}else{
		$newline =~ s/.$//;
		print OUT "$newline;\n";
		$count=0;
		$newline ="INSERT INTO $updateTable VALUES ";
	}
}
$newline =~ s/,$//;
print OUT "$newline;\n";


foreach my $line(@updateout){
	print OUT $line;
}
print OUTRENAME "RENAME TABLE $destTable to $tempTable,$updateTable to $destTable;\n";
print OUTRENAME "UPDATE METADATA SET tableName='$tempTable' where tableName='$destTable';\n";
print OUTRENAME "INSERT INTO METADATA VALUES ('$destTable','$file','$md5sum','$date');\n";
print OUTRENAME "INSERT INTO METADATA VALUES ('$destTable','$mapfile','$md5summap','$date');\n";
close OUT;
close OUTRENAME;


sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}

sub procRec(){
	my $omimID;
	my $title;
	my @genes;
	my $desc;
	my $clinicfeats;
	my $cs;
	my @record=@_;
        my $omimIDpos=0,$titlepos=0,$descpos=0,$cfpos=0,$cspos=0;
	# find positions of tag
	for (my $i=0; $i <=$#record; $i++){
     		if ($record[$i] eq "*FIELD* NO"){
			$omimIDpos=$i+1;
		}		
        	if ($record[$i] eq "*FIELD* TI"){
            		$titlepos=$i+1;
       		}
		if ($record[$i] eq "DESCRIPTION"){
			$descpos=$i+2;
		}
		if ($record[$i] eq "CLINICAL FEATURES"){
			$cfpos=$i+2;
		}
		if ($record[$i] eq "*FIELD* CS"){
			$cspos=$i+1;
		}
   	 }
	 #find OMIMID
	 my $omimID = $record[$omimIDpos];
	 #find Title
	 for(my $i = $titlepos; $i <= $#record; $i++){
         	#last if($record[$i]=~/^;;/);
         	last if($record[$i]=~/^\*FIELD\*/);
         	$title.=$record[$i]." ";
         }
	 $title =~ s/\'/\\'/g;
         chop($title);
         return if($title=~/MOVED TO/);
	 #find more Genes in title
	 my @cols=split(/;;/,$title);
	 foreach my $col (@cols){
		$col  =~ s/^\s+|\s+$//g;
		my @subcols=split(/;/,$col);
		if ($#subcols>0 && $col !~ /INCLUDED|FORMERLY/){
			$gene=uc($subcols[1]);
			$gene =~ s/^\s+|\s+$//g;
			push(@{$idhash{$omimID}},$gene);
			@{$idhash{$omimID}}=uniq(@{$idhash{$omimID}});
		}
		if (length($col) < 6){
			$gene = uc($col);
			push(@{$idhash{$omimID}},$gene);
                        @{$idhash{$omimID}}=uniq(@{$idhash{$omimID}});
		}
	 }
	 #find Description
	 if ($descpos != 0){
	 	for ($i = $descpos; $i <= $#record; $i++){
			$uc=uc($record[$i]);
			last if($record[$i] eq $uc && $record[$i] ne "");
			last if($record[$i]=~/^\*FIELD\*/);
			$desc.=$record[$i]." ";
	 	} 
		$desc =~ s/\'/\\'/g;
	}
	else {
		$desc="";
	}
	# find Clinical Features
	if ($cfpos != 0){
		for ($i = $cfpos; $i <= $#record; $i++){
                	$uc=uc($record[$i]);
                	last if($record[$i] eq $uc && $record[$i] ne "");
                	last if($record[$i]=~/^\*FIELD\*/);
			if ($record[$i] eq ""){
				$record[$i]="<br>";
			}
                	$clinicfeats.=$record[$i]." ";
         	}
	 	$clinicfeats =~ s/\'/\\'/g;
	}
	else{
		$clinicfeats="";
	}
	# find CS
	if ($cspos != 0){
		for ($i = $cspos; $i <= $#record; $i++){
			last if($record[$i]=~/^\*FIELD\*/);
			$cs.=$record[$i]." ";
		}
		$cs =~ s/\'/\\'/g;
	}
	else{
		$cs="";
	}
	# find GENE based on hash from genemap
        if (exists $idhash{$omimID}){
		@genes=@{$idhash{$omimID}};
	}
	if ($#genes != -1){
		foreach $gene (@genes){
#			print OUT "$omimID\t$gene\t$title\t$desc\t$clinicfeats\t$cs\n";
		        push(@tempout, "('$omimID','$gene','$title','$desc','$clinicfeats','$cs'),");
			#push(@updateout, "INSERT INTO $destTable (omimID,gene,title,description,clinFeat,clinFeatSum) VALUES ('$omimID','$gene','$title','$desc','$clinicfeats','$cs') ON DUPLICATE KEY UPDATE omimID='$omimID',gene='$gene',title='$title',description='$desc',clinFeat='$clinicfeats',clinFeatSum='$cs';\n");
		}
	}else{
		push(@tempout,"('$omimID','','$title','$desc','$clinicfeats','$cs'),");
	#	print OUT "$omimID\tnull\t$title\t$desc\t$clincfeats\t$cs\n";
	#push(@updateout, "INSERT INTO $destTable (omimID,gene,title,description,clinFeat,clinFeatSum) VALUES ('$omimID','NULL','$title','$desc','$clinicfeats','$cs') ON DUPLICATE KEY UPDATE omimID='$omimID',gene='NULL',title='$title',description='$desc',clinFeat='$clinicfeats',clinFeatSum='$cs';\n");
	}
}

