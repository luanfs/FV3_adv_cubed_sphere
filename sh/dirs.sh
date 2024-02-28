#!/bin/bash
#Script to set directory structure

echo
echo "Setting directory structure..."
echo
if [ -d bin ] ; then
    echo "Binary files will be placed inside          bin/" 
else
    mkdir bin
    echo "Directory 'bin' created" 
    echo "Binary files will be placed inside          bin/" 
fi

if [ -d grid ] ; then
    echo "Output grid files will be put inside        grid/"
else
    mkdir grid
    echo "Directory 'grid' created" 
    echo "Grid files will be put inside               grid/"
fi

if [ -d data ] ; then
    echo "Output data files will be put inside        data/"
else
    mkdir data
    echo "Directory 'data' created" 
    echo "Output data files will be put inside        data/"
fi

if [ -d par ] ; then
    echo "Parameters files are, or should be, inside  par/"
else
    mkdir par
    echo "Directory 'par' created" 
    echo "Parameters files are, or should be, inside  par/"
fi

if [ -d graphs ] ; then
    echo "Graphs will be placed inside                graphs/"
else
    mkdir graphs
    echo "Directory 'graphs' created" 
    echo "Graphs will be places inside                graphs/"
fi
echo
echo "The mesh parameters file can be passed as argument, "
echo "or be named as 'mesh.par' in the directory par/"
echo

echo "The simulation and model parameters file must be passed as a file "
echo " named as 'par/*.par' "
echo
echo "End of directory structuring."
echo
echo " * Consider commenting unused modules from the Makefile to speed up compilation time"
echo

