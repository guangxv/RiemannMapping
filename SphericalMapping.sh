#!/bin/bash

<<COMMENT1
LevelSet_Cardiac
COMMENT1

echo "	Original Surface
	Spherical Conformal Map"

i=1
while [ $i -lt 2 ]
do
	inputPath="/home/liguangxu/Projects/Data/Meshes/Tobata/Chest/Case$i/"
	outputPath=$inputPath

	# do the levelset algorithm
	./Build/bin/SphericalMapping \
	"$inputPath"LevelSet_Manual_Cardiac.off\
	"$outputPath"SphericalMapping_Cardiac.off \
	0.0005 0.00001

i=`expr $i + 1`
done
