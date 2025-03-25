## .pl <gff> <fasta>
use strict;


my ($gtf,$genome) = @ARGV;

my %genes; my %exon; my $geneID; 
my %gene_orient; my %chr_genes;
open GTF,"$gtf" or die $!;
while (my $line = <GTF>)
{
	if($line =~ /\A#/){next}

	my($chr,$flag,$start,$end,$orient,$info) = (split /\t/,$line)[0,2,3,4,6,-1];
	#print "$chr,$flag,$start,$end,$orient,$info\n";
	if($flag eq 'transcript')
	{
		#print "$line\n";
		my @info_words = split /[;= ]/,$info;
		for(my $i=0;$i<@info_words;$i+=2)
		{
			if($info_words[$i] eq 'ID')
			{
				my $j = $i +1;
				$geneID = $info_words[$j];
				#print "$geneID\n";
				last;
			}
		}
		#$genes{$geneID} = $start.'-'.$end;
		$gene_orient{$geneID} = $orient;
		$chr_genes{$geneID} = $chr;
		next;
	}
	if($flag eq 'CDS')
	{
		$exon{$geneID} .= $start.'-'.$end.'_';
	}
}
close GTF;

my %introns; 
foreach my $gene (keys %gene_orient)
{
	my $intro_num;
	my @exons = sort{$a<=>$b}split /_/,$exon{$gene};
	
	if(@exons >1)
	{
		print "@exons\t$gene_orient{$gene}\n";
		@exons = split /-/,(join "_",@exons);

		for(my $i=1;$i<@exons-1;$i++)
		{
			my($intro_s,$intro_e) = split /_/,$exons[$i];
			$intro_s++;
			$intro_e--;

			$intro_num++;
			my $intro_ID = $gene."-intron-$intro_num";
			$introns{$chr_genes{$gene}}{$intro_ID} = $gene_orient{$gene}."\t".$intro_s."\t".$intro_e;
		}		
	}
}

my $chr_ID = "IsItTheFirstLine";
my $fasta;
open GENOME,"$genome" or die $!;
while (my $line = <GENOME>)
{
	chomp($line);
	if($line =~ s/>//)
	{
		if($chr_ID eq "IsItTheFirstLine")
		{
			$chr_ID = (split /[\t ]/,$line)[0];
			next;
		}
		else
		{
			foreach my $intro (keys %{$introns{$chr_ID}})
			{
				my ($orient,$intro_s,$intro_e) = split /\t/,$introns{$chr_ID}{$intro};
				my $intro_length = $intro_e - $intro_s +1;
				print STDERR "$chr_ID\t$intro\t$introns{$chr_ID}{$intro}\t$intro_length\n";

				$intro_s--;
				my $intro_fa = substr($fasta,$intro_s,$intro_length);
				if($orient eq '-')
				{
					$intro_fa =~ tr/[AGCTNacgtn]/[TCGANtgcan]/;
					$intro_fa = reverse($intro_fa);
				}
				print ">$intro\n$intro_fa\n";
			}

			$chr_ID = (split /[\t ]/,$line)[0];
			$fasta = ();
		}
	}
	else
	{
		$fasta .= $line;
	}

}
close GENOME;

foreach my $intro (keys %{$introns{$chr_ID}})
{
	my ($orient,$intro_s,$intro_e) = split /\t/,$introns{$chr_ID}{$intro};
	my $intro_length = $intro_e - $intro_s +1;
	print STDERR "$chr_ID\t$intro\t$introns{$chr_ID}{$intro}\t$intro_length\n";

	$intro_s--;
	my $intro_fa = substr($fasta,$intro_s,$intro_length);
	if($orient eq '-')
	{
		$intro_fa =~ tr/[AGCTNacgtn]/[TCGANtgcan]/;
		$intro_fa = reverse($intro_fa);
	}
	print ">$intro\n$intro_fa\n";
}
