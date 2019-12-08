#!/bin/bash
#prun -v -1 -np 1 ./apsp -c ../input/V5000-E500000 -n 5000 -e 500000

echo "Running APSP for all input files..."
threads=16
for dir in ../input/*
do
    	echo "------------------------NEXT FILE------------------------"
	tmp=$(basename "$dir")
	nodes=${tmp##*V}
	edges=${tmp##*E}
	nodes=${nodes%-*}

	prun -v -1 -np 1 ./apsp -c $dir -n $nodes -e $edges -p 16
	#echo "prun -v -1 -np 1 ./apsp -c $tmp -n $nodes -e $edges"
	
done

