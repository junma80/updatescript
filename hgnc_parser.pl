use strict;
use List::MoreUtils 'any';
use List::MoreUtils qw(firstidx);
my $inmap="hgnc_update.txt";
my $outfile="hgnc_hugo_update_fix.sql";
my $inhugo="hugoids.txt";

open(INMAP,"<$inmap");
open(OUTFILE, ">$outfile");
open(INHUGO,"<$inhugo");
my $header= <INMAP>;
chomp($header);
my @headercols=split(/\t/,$header);

my %hugonamehash;
while (my $line =<INHUGO>){
	chomp($line);
	$hugonamehash{$line}=1;
}
close(INHUGO);

my $newsymbolpos=firstidx {$_ eq 'Approved Symbol'} @headercols;
my $oldsymbolpos=firstidx {$_ eq 'Synonyms'} @headercols;
my $presymbolpos=firstidx {$_ eq 'Previous Symbols'} @headercols;
my $statuspos=firstidx {$_ eq 'Status'} @headercols;
my $hgncidpos=firstidx {$_ eq 'HGNC ID'} @headercols;

my %hgncnamehash;
while (my $line =<INMAP>){
	chomp($line);
	my @cols=split(/\t/,$line);
	my $newsymbol=$cols[$newsymbolpos];
	$hgncnamehash{$newsymbol}=1;
}
close(INMAP);

my %maphash;
open(INMAP,"<$inmap");
<INMAP>;
my %specialhash=("CCBP2" => "ACKR2", "DBC1" => "BRINP1","KAP" => "CDKN3", "CCRL1" => "ACKR4","RIG" => "DIRAS1",
		 "RAD26L" => "ERCC6L2","PAR4" => "PAWR", "F379" => "FAM138A" ,"C11orf2" => "VPS51", "IDDM1" => "HLA-DQB1",
		 "SMAP" => "KIFAP3", "MLL2" => "KMT2D", "MLL4" => "KMT2B", "PWCR"=>"SNRPN", "RITA" => "ZNF331", "U2" => "RNU2-1",
		 "U3" => "RNU3P1" ,"U4" =>"RNU4-1" ,"U6" => "RNU6-50P", "PAR1" => "PWAR1" ,"BOP" => "SMYD1",
		 "ODZ3"=>"TENM3","TCRBV12S2" => "TRBV12-2" ,"TCRBV2S1" => "TRBV20-1","TCRBV21S1" =>"TRBV21-1",
		 "TCRBV15S1"=>"TRBV24-1","TCRBV14S1" => "TRBV27");
my $newline ="INSERT INTO HGNC_fix (hugoID, hgncName, hgncID) VALUES";

while (my $line = <INMAP>){
	chomp($line);
	my @cols=split(/\t/,$line);
	my $newsymbol=$cols[$newsymbolpos];
	my $oldsymbol=$cols[$oldsymbolpos];
	my $presymbol=$cols[$presymbolpos];
	my $status=$cols[$statuspos];
	my $hgncID=$cols[$hgncidpos];
	if ($status eq "Approved"){
		if ($oldsymbol ne "" || $presymbol ne "" ){
			my @oldsymbolcols=split(/,\s/,$oldsymbol);
			foreach my $ele (@oldsymbolcols){
				my $uc=uc($ele);
				if (not exists $hgncnamehash{$ele} and not exists $specialhash{$ele} and $ele eq $uc and $hugonamehash{$ele}){					
					$ele =~ s/\'/\\'/g;
					#$newline.= "('$ele','$newsymbol',0),";
					$maphash{$ele}=$newsymbol;
					
				}elsif (exists $specialhash{$ele}){
					
					$maphash{$ele}=$specialhash{$ele};
				}
#				print "$ele\t$newsymbol\n";
			}
			my @presymbolcols=split(/,\s/,$presymbol);
			foreach my $ele (@presymbolcols){
				my $uc=uc($ele);
				if (not exists $hgncnamehash{$ele} and not exists $specialhash{$ele} and $ele eq $uc and $hugonamehash{$ele}){
					$ele =~ s/\'/\\'/g;
					#$newline.="('$ele','$newsymbol',0),";
					$maphash{$ele}=$newsymbol;
		
				}elsif (exists $specialhash{$ele}){
	
					$maphash{$ele}=$specialhash{$ele};
				}	
#				print "$ele\t$newsymbol\n";
			}
		}
		#$newline.="('$newsymbol','$newsymbol',0),";
		$maphash{$newsymbol}=$newsymbol;
	}	
}
foreach my $ele (keys %maphash){	
	$newline.="('$ele','$maphash{$ele}',0),";
}

$newline =~ s/,$/;\n/;
print OUTFILE $newline;

