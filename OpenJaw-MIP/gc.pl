my $seq = "CCCAAAATCTGTGATCTTGAC";
my @gc = $seq=~/G|C/ig;
print scalar(@gc)/length($seq),"\n";

