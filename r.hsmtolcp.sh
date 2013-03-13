#!/bin/sh
# Author: Micha Silver, August 2012; copyright: GPL >= 2
# Purpose: Calculate a Least Cost Path grid of Population Migration
# Input: 
# 	* A Habitat Suitability Model (from MaxEnt)
#	* A CSV file of sampling localities
#	* A list of pairs of localities with Fst = 0
#	* Localities will have a 3 letter code for convenience in naming raster

################################################################
# Part 1: Create friction raster
################################################################
# Import the HSM tiff into GRASS
r.in.gdal input=HSM_last.tif output=hsm title="Habitat Suitability Model"
# Set the GRASS region to match that raster
g.region -p rast=hsm
# Create friction raster as aritmetic inverse of HSM
r.mapcalc friction=1.0-hsm
r.support map=friction title="Friction map"
# DEBUG: info for friction map
# r.info friction

################################################################
# Part 2: Create cost raster for each locality
################################################################
# Change field separator to comma
orig_IFS=$IFS
IFS=,
while read name code lat lon; do 
	echo "Creating $code cost map"
	r.cost --verbose --overwrite input=friction out=cost_"$code" coord="$lon,$lat"; 
done < localities.csv 
IFS=$orig_IFS
# DEBUG: list results
# g.mlist type=rast pattern="cost*"

################################################################
# Part 3: Create LCP corridors for each non-structured pair
################################################################
# Create corridors for each locality pair by adding the two cost rasters 
# Loop thru all non structured pairs from file non_structured_pairs.txt
while read loc1 loc2; do 
	echo "Calculating LCP map for $loc1 and $loc2"
	# Round to integer values (for reclass step)
	r.mapcalc corridor_"$loc1"_"$loc2"="round(cost_"$loc1" + cost_"$loc2")"
	# DEBUG: Examine univariate statistics of corridor raster
	# r.univar -g corridor_"$loc1"_"$loc2"_int
	# Save minimum value
	min=`r.univar -g corridor_"$loc1"_"$loc2" | grep min | cut -d= -f2`

	# Create variables for reclass rules file
	two_percent=`echo $min*1.02/1 | bc`
	five_percent=`echo $min*1.05/1 | bc`
	#seven_percent=`echo $min*1.07/1 | bc`

	# Now put variables into reclass file
cat << EOF > "$loc1"_"$loc2"_rules.txt
$min = 3			Minimum
$min thru $two_percent = 2		One percent
$two_percent thru $five_percent = 1		Two percent
* = NULL			NULL
end
EOF
	# DEBUG: examine rules file
	# cat "$loc1"_"$loc2"_rules.txt

	# Make a reclassed raster based on that file
	r.reclass --overwrite input=corridor_"$loc1"_"$loc2" output=lcp_"$loc1"_"$loc2" rule="$loc1"_"$loc2"_rules.txt title="$loc1 - $loc2 Least Cost Path"
	# Set color ramp for lcp map
	r.colors --quiet -n map=lcp_"$loc1"_"$loc2" color=ryg

done < Fst_zero.txt
# End of operations on each pair of localities

#################################################################
# Part 4: Merge LCP rasters together
#################################################################
# Combine all reclass LCP maps into a new merged LCP map
lcp_maps=`g.mlist type=rast pattern="lcp*" separator=,`
r.series --overwrite input=$lcp_maps output=lcp_merged method=sum
# Set color ramp for lcp_merged map
r.colors --quiet -n map=lcp_merged color=ryg

