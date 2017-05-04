#!/usr/bin/perl

while(<>){
    if($_ =~/>/){
        print ">header\n";
    }else{
	$lestT .= $_;
    }
}

$lestT = uc $lestT;
$lestT =~ s/\n//g;
$lestT =~ s/\r//g;
$size = length($lestT);
$lest ="";
for($i=0;$i<$size;$i+=60){
    $lest .=substr($lestT,$i,60);
    $lest .="\n";
}
    $lest .=substr($lestT,$i);

print $lest;