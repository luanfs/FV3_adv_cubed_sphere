#!/bin/bash

#Script to create a link for the main executable on the main directory
version=` date +%y.%m.%d `
binname="bin/main$version"
echo $binname
cp bin/main $binname
rm -rf main
if [[ -f main ]] ; then
    ln -n -s -f $binname main
    echo "A link for the main executable (named 'main') was updated " $binname  
else
    if [[ -f bin/main ]] ; then
	ln -s $binname main
	echo "A link for the main executable (named 'main') was created " $binname 
    else
	echo "Could not create link, bin/main does not exist" $binname 
    fi
fi
echo
