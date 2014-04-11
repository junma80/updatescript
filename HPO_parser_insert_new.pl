use strict;
my $date =`date +"%Y-%m-%d %k:%M:%S"`;
chomp($date);
my @datecols=split(/-/,$date);
my %datehash=("01" => "12","02" => "01","03" => "02", "04" => "03", "05" => "04", "06" => "05", "07" => "06", "08" => "07", "09" => "08", "10" => "09", "11" => "10", "12" =>"11");
my %dateyearhash=("2015" => "2014", "2016" => "2015", "2017" => "2016");
my $datestring="$datecols[1]_$datecols[0]";
my $date_old="$datehash{$datecols[1]}_$datecols[0]";
my $year=$datecols[0];
if ($datecols[1] eq "01"){
	$year=$dateyearhash{$datecols[0]};
}
my $outfolder="/space1/databaseUpdates/";

#my $datestring="01_2014";
my $inmap="$datestring.ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt";
my $insql="$datestring.MYHPO.sql";
my $inobo="$datestring.hp.obo";
my $updatesql="$datestring.MYHPO_insert_new.sql";
my $changemeta="$datestring.HPO_changemeta.sql";

#download HPO files (need to change the sql file name)
unless (-e "$outfolder$inmap"){
	`wget http://compbio.charite.de/hudson/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt -O $inmap`;
	`mv $inmap $outfolder`;
}
unless (-e "$outfolder$insql"){
	`wget http://compbio.charite.de/hudson/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/MYHPO_$datestring.sql -O $insql`;
	`mv $insql $outfolder`;
}

unless (-e "$outfolder$inobo"){
	`wget http://compbio.charite.de/hudson/job/hpo/lastStableBuild/artifact/ontology/release/hp.obo -O $inobo`;
	`mv $inobo $outfolder`;
}

open(INMAP,"<$outfolder$inmap");
open(SQL,"<$outfolder$insql");
open(INOBO,"<$outfolder$inobo");
open(UPDATESQL,">$outfolder$updatesql");
open(CHANGEMETA,">$outfolder$changemeta");

my $md5sum = `md5sum $outfolder$inmap`;
my $file = `ls -lah $outfolder/$inmap`;
chomp($md5sum);
chomp($file);
chomp($date);


# make file for HPO table
<INMAP>;
print UPDATESQL "DROP TABLE IF EXISTS HPO_update;\n";
print UPDATESQL "DROP TABLE IF EXISTS HPOgraph_path_update;\n";
#print UPDATESQL "DROP TABLE IF EXISTS HPOphenotype_update;\n";
print UPDATESQL "DROP TABLE IF EXISTS HPOterm_update;\n";
#print UPDATESQL "DROP TABLE IF EXISTS altid;\n";
#print UPDATESQL "ALTER TABLE HPOterm add definition varchar(1500);\n";
#print UPDATESQL "ALTER TABLE HPOterm add synonyms varchar(1000);\n";


print UPDATESQL "CREATE TABLE HPOterm_update LIKE HPOterm;\n";
print UPDATESQL "CREATE TABLE HPO_update LIKE HPO;\n";
print UPDATESQL "CREATE TABLE HPOgraph_path_update like HPOgraph_path;\n";
#print UPDATESQL "CREATE TABLE HPOphenotype_update like HPOphenotype;\n";;
#print UPDATESQL "INSERT HPOphenotype_update select * from HPOphenotype;\n";
#print UPDATESQL "ALTER TABLE HPOphenotype_update add constraint HPOphenotype_ibfk2_$datestring FOREIGN KEY (patID) REFERENCES patient (id) on DELETE CASCADE;\n";
#print UPDATESQL "CREATE TABLE altid (oldid int(11), newid int (11), description varchar (254) )ENGINE=InnoDB;\n";

my @temp=();

# insert into  HPO_temp table
my $newhpo="INSERT INTO HPO_update VALUES ";
my @hpoarray=();
my %deschash;
while (my $line = <INMAP>){
    chomp($line);
    my @A = split/\t/,$line;
    my $omim = substr($A[0], 5);
    my $hp = substr($A[3], 3);
    my $hpi = int($hp);
    my $desc = $A[4];
    $desc =~s/'/\\'/g;
    $deschash{$hpi}=$desc;
    my $gene = $A[1];
    push(@hpoarray,"('NULL',$hpi,'$desc','$gene',0),");
}
@hpoarray=uniq(@hpoarray);
foreach my $hpo (@hpoarray){
	$newhpo.=$hpo;
}
$newhpo =~ s/,$/;\n/;
push(@temp, $newhpo);
#push(@temp,"UNLOCK TABLES;\n");


# parse hp.obo, find comment defenition, synonym
my %cmthash;
my %defhash;
my %synhash;
my %altidhash;
my @record;
while (my $line = <INOBO>){
	chomp($line);
	if ($line eq "[Term]"){
		my $id;
		my $altid;
		my $comment="NULL";
		my $def="NULL";
		if ($#record > 1){
			foreach my $rec (@record){
				$rec=~ s/\\*//g;
				if ($rec =~ /^id/){
					my @cols=split/:/,$rec;
					$id=$cols[2];
					$id=~ /(0*)(\d*)/;
					$id=$2;	
				}
				if ($rec=~ /^alt_id/){
					my @cols=split(/HP:/,$rec);
					$altid=$cols[1];
					$altid=~ /(0*)(\d*)/;
					$altid=$2;
					push(@{$altidhash{$id}},$altid);
				}
				if ($rec =~ /^def:/){
					my @cols=split(/"/,$rec);
					$def=$cols[1];
					$def=~ s/\'/\\'/g;
					$defhash{$id}=$def;
				}
				if ($rec =~ /^synonym:/){
					my @cols=split(/"/,$rec);
					push @{$synhash{$id}},$cols[1];
					
				}	
				if ($rec =~ /^comment:/){
					my @cols=split(/:/,$rec);
					$comment=$cols[1];
					$comment =~ s/^\s+|\s+$//g;
					$comment =~ s/\'/\\'/g;
					$cmthash{$id}=$comment;
				}
			}
		}
		@record=();
	}
	push (@record, $line);
}
=pod
my $altidline="INSERT INTO altid VALUES ";
foreach my $key (keys %altidhash){
	foreach my $altid (@{$altidhash{$key}}){
		$altidline.="($altid,$key,'$deschash{$key}'),";
	}
}
$altidline =~ s/,$/;\n/;
push(@temp,$altidline);
=cut
#make file for graph_path,term
while (my $line = <SQL>){
	chomp($line);
	if ($line =~ /^INSERT INTO `graph_path` VALUES*/){
		#push(@temp, "LOCK TABLES `HPOgraph_path_update` WRITE;\n");
		$line =~ s/graph_path/HPOgraph_path_update/g;
		push(@temp, "$line\n");
		#push(@temp, "UNLOCK TABLES;\n");
        }
	if ($line =~ /^INSERT INTO `term` VALUES*/){
                #push(@temp, "-- HPO term table\n");
		#push(@temp, "LOCK TABLES `term_temp` WRITE;\n");
		#push(@temp, "LOCK TABLES `HPOterm_update` WRITE;\n");
		procterm($line);
		#push(@temp,"UNLOCK TABLES;\n");
	}
}
#push (@temp,"update altid inner join HPO_update on altid.newid=HPO_update.hpoID set altid.descritpion=HPO_update.descritpion;\n");


push (@temp,"INSERT INTO HPOterm_update (select * from HPOterm where id not in (select id from HPOterm_update));\n");
push (@temp,"INSERT INTO HPOgraph_path_update (select * from HPOgraph_path where term1_id not in (select term1_id from HPOgraph_path_update));\n");
push (@temp,"INSERT INTO HPOgraph_path_update (select * from HPOgraph_path where term2_id not in (select term2_id from HPOgraph_path_update));\n");

#push(@temp, "UPDATE HPOphenotype_update set HPOID =(select altid.newid from altid where altid.oldid=HPOphenotype_update.HPOID) where HPOID =(select altid.oldid from altid where altid.oldid=HPOphenotype_update.HPOID);\n");
push(@temp, "INSERT INTO HPO_update (select * from HPO where not exists(select gene from HPO_update where HPO.gene=HPO_update.gene and HPO.hpoID=HPO_update.hpoID) and HPO.isCodified=(1));\n");
#push(@temp, "update HPO_update inner join altid on HPO_update.hpoID=altid.oldid set HPO_update.description=altid.description;\n");
#push(@temp, "update HPO_update inner join altid on HPO_update.hpoID=altid.newid set HPO_update.description=altid.description;\n");
#push(@temp, "update HPO_update set hpoID =(select altid.newid from altid where altid.oldid=HPO_update.hpoID) where HPOID =(select altid.oldid from altid where altid.oldid=HPO_update.hpoID);\n");

#push(@temp, "DROP table altid;\n");
foreach my $line(@temp){
	print UPDATESQL $line;
}

#print UPDATESQL "update HPO_update set isCodified=(select distinct(isCodified) from HPO where HPO_update.hpoID=HPO.hpoID and HPO_update.gene=HPO.gene);";

print UPDATESQL "ALTER TABLE HPOgraph_path_update add constraint fk_graph_path_term1_$datestring FOREIGN KEY (term1_id) references HPOterm_update(id);\n";
print UPDATESQL "ALTER TABLE HPOgraph_path_update add constraint fk_graph_path_term2_$datestring FOREIGN KEY (term2_id) references HPOterm_update(id);\n";
#print UPDATESQL "ALTER TABLE HPOphenotype_update add constraint HPOphenotype_ibfk1_$datestring FOREIGN KEY (HPOID) references HPOterm_update(id);\n";
print UPDATESQL "ALTER TABLE HPO_update add constraint HPO_ibfk_$datestring FOREIGN KEY (hpoID) references HPOterm_update(id);\n";

close UPDATESQL;

print CHANGEMETA "RENAME TABLE HPOgraph_path to HPOgraph_path_$date_old;\n";
print CHANGEMETA "RENAME TABLE HPOterm to HPOterm_$date_old;\n";
print CHANGEMETA "RENAME TABLE HPO to HPO_$date_old;\n";
#print CHANGEMETA "RENAME TABLE HPOphenotype to HPOphenotype_$date_old;\n";
print CHANGEMETA "RENAME TABLE HPOgraph_path_update to HPOgraph_path;\n";
print CHANGEMETA "RENAME TABLE HPOterm_update to HPOterm;\n";
print CHANGEMETA "RENAME TABLE HPO_update to HPO;\n";
#print CHANGEMETA "RENAME TABLE HPOphenotype_update to HPOphenotype;\n";
print CHANGEMETA "UPDATE METADATA set tablename='HPO_$date_old' where tablename='HPO';\n";
print CHANGEMETA "INSERT INTO METADATA VALUES ('HPO','$file','$md5sum','$date');\n";

close CHANGEMETA; 

sub uniq {
	return keys %{{ map { $_ => 1 } @_ }};
}



sub procterm(){
	 my $line=$_[0];
        $line =~ s/INSERT INTO `term` VALUES //g;
        my @cols=split(/\),\(/,$line);
	my $newterm="INSERT INTO HPOterm_update VALUES ";
        foreach my $record (@cols){
                $record =~ s/[\(|\);]//g;
                chomp($record);
                $record =~ /(\d+),{1}(.+),{1}(\d+),{1}(\d+),{1}(.+),{1}(.+),{1}(.+)/;
		my $synonym;
		my $cmt=$cmthash{$1};
                if (exists $synhash{$1}){
			 foreach (@{$synhash{$1}}){
				my $syn=$_;
				$syn =~ s/\'//g;			
                		$synonym.="$syn, ";
       			 }
			 $synonym =~ s/, $//g;
		}
		my $def;
		if (exists $defhash{$1}){
			$def=$defhash{$1};
			$def =~ s/\'//g;
		}
		$newterm.="($1,$2,$3,$4,$5,'$cmt',$7,'$def','$synonym'),";
	
        }
	$newterm =~ s/,$/;\n/;
	push (@temp,$newterm);
}
