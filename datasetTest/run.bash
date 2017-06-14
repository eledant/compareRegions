#!/bin/bash

rm statsSummary.txt

path='test4_fixDatasetA_diffRegionB/20000000TotalBP_1%Overlap_Test3_B>A'
for id in '25' '50' '100' '150' '200' '250' '500' '750' '1000' '1250' '1500' '1750' '2000' '2500'
do
	../compareRegions.py -p -l all -m 100 -n 100 genome_test.bed $path/datasetA.bed $path/datasetB_$id.bed	
done	
