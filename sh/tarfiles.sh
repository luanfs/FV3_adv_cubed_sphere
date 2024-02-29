#!/bin/bash
# Creates a tarball of the bacukp files
directory="FV3_adv_cubed_sphere"

# Source code files

srcdir="../"
sourcefiles="$srcdir/tools $srcdir/driver $srcdir/model $srcdir/sh $srcdir/run $srcdir/plot $srcdir/Makefile $srcdir/main.f90"
files="$sourcefiles"

# Output name
output=$directory".tar.bz2"

# Creates the tarball
tar cjfv $output $files

echo "File " $output " ready!"
echo
