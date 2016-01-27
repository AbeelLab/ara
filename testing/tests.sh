#!/bin/sh

ARA_LOC='../build/ara.jar'

ARACMD="java -jar $ARA_LOC"

# Find unique markers for a particular gene
echo "TEST1"
$ARACMD unique-markers test1.fasta test1.fasta 21 test1.out

# Find unique markers that uniquely identify two orthologs
echo "TEST2"
$ARACMD unique-markers test1.fasta,test2.fasta test1.fasta,test2.fasta 21 test2.out

# Find unique markers for a particular gene that is unique also from anything in an ortholog, return redundant markers
echo "TEST3"
$ARACMD unique-markers test1.fasta test1.fasta,test2.fasta 21 test3.out R
