use strict;

my $Phage_Coords = $ARGV[0];

open(PC, $Phage_Coords) || die "Can't find coordinates file, $Phage_Coords $!";

my @Phage_Start = ();
my @Phage_End = ();

while(<PC>)  ## Read phage coordinates into the lists, @Phage_Start and @Phage_End
 {
  if($_ =~ /^[1-9]/)
  {
     my @line = split(/\s+/, $_);
     push @Phage_Start, $line[0];
     push @Phage_End, $line[1];  
  }
  
 }

my $Fasta = $ARGV[1];

## assume $accessions were parsed and written in the same order as HMMsearch rankings

open(FAS, $Fasta) || die "Can't find nucleotide fasta file, $Fasta $!";

my @Fasta = <FAS>;

my $head = '';

my @DNA = ();

my $l = 0;

foreach(@Fasta)
 {
      if($Fasta[$l] =~ /^>/)
        {
       	   $head = $Fasta[$l];
        }
      elsif($Fasta[$l] =~ /^(A|T|C|G|N)/i)
       {
         chomp $Fasta[$l];
         my @line = split('', $Fasta[$l]);
         push @DNA, @line; 
       }
     $l++;
  } 

my $ori = 0;


## Now, print each base pair of @DNA falling outside phage regions
##      otherwise, print 'N'

  my @Masked_Chromosome = ();
  
  my $dna = 0; ## Index for chromosome
  my $phage = 0; ## Index for prophage coordinates list
  
for(my $dna = 0; $dna < scalar @DNA; $dna++)
 {
     if($dna < $Phage_Start[0])
      {
	 # print "0..".$Phage_Start[0]."\n";
        push @Masked_Chromosome, @DNA[$dna]; 
      }
   
     my $limit = scalar @Phage_Start - 1;
  
    for(my $phage = 0; $phage < scalar @Phage_Start; $phage++)
     {
      if(($dna >= $Phage_Start[$phage]) && ($dna < $Phage_End[$phage]))
       {
         # print $Phage_Start[$phage]."..".$Phage_End[$phage]."\n";
	   push @Masked_Chromosome, 'N';
       }
      elsif(($dna >= $Phage_End[$phage - 1]) && ($dna < $Phage_Start[$phage]) )
       {
         # print $Phage_End[$phage - 1]."..".$Phage_Start[$phage]."\n";
           push @Masked_Chromosome, $DNA[$dna];   
       } 
     }
  
     if($dna >= $Phage_End[$limit])
     {
 	 # print $Phage_End[$limit]."..End\n";
         # print $DNA[$dna]."\n";
          push @Masked_Chromosome, $DNA[$dna];   
     }
       
      
 }

  
   my $pcon = 0;
   
   chomp $head;

   ## Now, print masked chromosome

  # foreach(@Masked_Chromosome)
  #  {
       print $head."\n";
   
       my $offset = 0;

       while( ($offset + 60) < scalar @Masked_Chromosome )
         {
	     print @Masked_Chromosome[$offset..($offset + 59)];
             print "\n";
           $offset = $offset + 60;
         }
 
    my $total_length = scalar @Masked_Chromosome - 1;

     print @Masked_Chromosome[$offset..$total_length];
     print "\n"; 
   

	$pcon++;
   # }
     
  
  print "\n";

  exit;
