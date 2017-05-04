#!/usr/bin/perl

while(<>){
    $_ =~ s/\r//;
    print $_;
}
