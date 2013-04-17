#!/usr/bin/perl

use strict;
use DBI;

my $server='localhost';
my $db='transfac';
my $username='root';
my $password='chenyuan';

my $dbh=DBI->connect("dbi:mysql:$db:$server",$username,$password);

my $query="select * from site";
my $sth=$dbh->prepare($query);

$sth->execute();

while (my $row=$sth->fetchrow_arrayref) {
	print join("\t",@$row),"\n";
}

$dbh->disconnect;



