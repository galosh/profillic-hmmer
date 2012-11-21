#!/usr/bin/perl 
## Output the first eighteen lines untouched
for($i = 0;  ($i < 18); $i++) {
  $_ = <STDIN>;
  print;
}

## Read in all subsequent lines.
## Split into tokens separated by strings of spaces.
## If a token has a decimal point in it, convert it to exp(-token) and output
## Otherwise just output it.

while (<STDIN>) {
  foreach $token (split(/  */)) {
    if ( $token =~ m/\./g ) {$token = sprintf("%9.5f",exp(-$token));}
    print $token . "\t";
  }
  print "\n";
}
