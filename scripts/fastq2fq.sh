#!/usr/bin/env sh

# Given a directory, convert all filenames containing 'fastq' with 'fq'
#
# usage: ./fastq2fq..sh dir

dir=$1

ls $dir/* | parallel 'mv {} {= s/fastq/fq/=}'
