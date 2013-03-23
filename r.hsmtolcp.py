#!/usr/bin/env python
 
"""
MODULE:    r.hsmtolcp.py
AUTHOR(S): Micha Silver <micha at arava co il>
PURPOSE:   Create predicted migration corridors for insects
		based on a habitat suitability raster produced using Maxent "presence only" software and
		genetic distance from an Fst matrix produced using Arlequin 
COPYRIGHT: (C) 2013 Micha Silver, Pablo Fresia and the GRASS Development Team
This program is free software under the GNU General Public License
(>=v2). Read the file COPYING that comes with GRASS for details.
"""

#%module
#% description: Creates a predicted migration corridors raster from an HSM, output from MaxEnt, and an Fst matrix, output from Arlequin
#% keywords: raster
#% keywords: habitat suitability
#% keywords: genetic distance
#% keywords: population genetics
#%end

#%option
#% key: hsm
#% type: string
#% description: Name of the Habitat suitability raster (in GeoTiff format) from MaxEnt
#% required: yes
#%end
#%option
#% key: friction
#% type: string
#% description: Name for GRASS friction raster
#% required: no
#% answer: friction

#%option
#% key: fst_matrix
#% type: string
#% description: Text file of Fst matrix from arlequin output
#% required: yes
#%end 
#%option
#% key: fst_max
#% type: double
#% description: Maximum Fst value to consider a pair of localities as "close" populations genetically
#% required: no
#% answer: 0
#%end 
#%option
#% key: pvalue_matrix
#% type: string
#% description: Text file of p-values matrix from arlequin output
#% required: no
#%end 
#%option
#% key: p_max
#% type: double
#% description: Maximum p-value to consider a pair of Fst values as statistcally significant
#% required: no
#% answer: 0.05
#%end 
#%option
#% key: localities
#% type: string
#% description: Text file (csv) sampling localities, with locality_code, locality_num, and x,y coordinates
#% required: yes
#%end 
#%option
#% key: corridors
#% type: string
#% description: Name for corridors output raster
#% required: no
#% answer: corridors
#%end 
#%option
#% key: lcp_merged
#% type: string
#% description: Name for Least Cost Path merged output raster
#% required: no
#% answer: lcp_merged

import os, sys, csv
import grass.script as grass
import numpy as np


def create_friction(hsm, friction):
	"""
	Import hsm Geotiff into GRASS and create friction map
	"""

	# Import hsm.tif as GRASS raster
	exists = grass.find_file(name = 'hsm', element = 'raster', quiet = True)
	if exists:
		grass.message("Habitat raster: hsm already exists and will be overwritten")

	grass.run_command('r.in.gdal', input=hsm, output="hsm", overwrite=True, quiet=True)
	grass.run_command('g.region', raster=hsm, quiet=True)

	# Create friction map
	exists = grass.find_file(name = friction, element = 'raster', quiet = True)
	if exists:
		grass.message("friction raster: "+friction +" already exists and will be overwritten")
	
	friction_expr = friction+" = 1.0-hsm"
	grass.mapcalc(friction_expr)



def create_localities(loc):
	"""
	Read and parse the csv file of localities
	Create a python list, where each entry is a list of ('Locality Name','Locality Num', 'X coord', Y coord')
	"""
	try:
		csv_file = open(loc,"rb")
	except IOError:
		grass.fatal("Cannot open localities file: "+loc)

	locs=csv.reader(csv_file, delimiter=',')
	loc_list=[]
	for l in locs:
		# Add to the loc_list a tuple containing the code, X coord and Y coord
		loc_list.append([l[1],l[2],l[3]])

	csv_file.close()
	return loc_list



def create_cost_maps(loc_list, friction):
	"""
	Use the friction map to create a least cost cost raster for each of the localities
	based on the friction map and locality coordinates
	The cost rasters are all named "cost"_localitycode
	"""
	cnt=0
	for i in range(len(loc_list)):
		cnt +=1
		code, x, y = loc_list[i][0], loc_list[i][2], loc_list[i][3]
		outrast="cost_"+code
		coords=x+","+y
		grass.run_command('r.cost', input=friction, output=outrast, coordinates=coords, overwrite=True, quiet=True)

	return cnt


def create_fst_pairs(fst, f_max, pval, p_maxi, num_locs):
	"""
	Use NumPy to import Fst and p-value matrices
	Use boolean comparisons to find which Fst values are <= f_max
	and which p-value values are <= p_max
	Create a list of pairs (of localities) which match both conditions
	"""
	
	# Import arrays, skip header row, and slice off first column
	fst_mat = np.genfromtxt(fst, skip_headers=1, delimiter=',', usecols=range(1,num_locs+1))
	pval_mat = np.genfromtxt(pval, skip_headers=1, delimiter=',', usecols=range(1,num_locs+1))
	# Make boolean arrays that reflect the max conditions for Fst and p-value 
	fst_bool = fst_mat <= f_max
	pval_bool = pval_mat <= p_max
	# Create array where both conditions apply (logical AND)
	acc_bool = np.logical_and(fst_bool, pval_bool)
	# Get indices where accept matrix is TRUE
	accept = np.where(acc_bool)
	# And convert to list
	pairs = np.transpose(accept).tolist()

	return pairs


def create_lcp_corridors(pairs, locs, lcp_merged):
	"""
	Loop thru all pairs of localities in the pairs list
	FOr each pair, find the codes of that pair from the localities list,
	add together the two cost maps of those two localities to create a corridor
	Get the minimum for each corridor using r.univar, 
	Use that minimum to make a reclass file with values:
		min					= 3			Min value
		min + 2% 		= 2			Two percent above
		min + 5%		=	1			Five percent above
		*						= NULL	NULL
	Run r.reclass to create uniform reclass rasters for each pair
	"""
	for i in range(len(pairs)):
		id1, id2 = pairs[i][0], pairs[i][1]
		code1, code2  = locs[id1][0], locs[id2][0]
		cost1, cost2  = "cost_"+code1, "cost_"+code2
		corridor = "corr_"+code1+"_"+code2
		corridor_expr = corridor+"=round("+cost1+"+"+cost2+")"
		grass.mapcalc(corridor_expr)
		u = grass.read_command('r.univar', map=corridor, flags="g", quiet=True)
		udict = grass.parse_key_val(u)
		min = float(udict['min'])
		two_pc = min*1.02
		five_pc = min*1.05
		# Create reclass file
		tmp_reclass = grass.tempfile()
		trc = open(tmp_reclass, "w")
		trc.write(min + " = 3	\t Minimum\n")
		trc.write(two_pc + " = 2 \t Two percent\n")
		trc.write(five_pc + " = 1 \t Five percent\n")
		trc.write("* = NULL \t NULL\n")
		trc.close()
		# Now create reclass raster
		lcp = "lcp_"+code1+"_"+code2
		grass.run_command('r.reclass',input=corridor,output=lcp, rule=tmp_reclass, overwrite=True, quiet=True)
		grass.run_command('r.colors', map=lcp, color="ryg", quiet=True)
		# Get rid of tmp reclass file
		unlink(tmp_reclass)

	# Now merge all lcp_* maps
	lcp_maps = grass.read_command('g.mlist', type="rast", pattern="lcp_*", separator=",")
	grass.run_command('r.series',input = lcp_maps, output = lcp_merged, method=sum, overwrite=True, quiet=True)
	grass.run_command('r.colors', map=lcp_merged, color="ryg", quiet=True)
	return len(pairs)


def main():
	""" 
	Get command line parameters
	"""
	if "GISBASE" not in os.environ:
		print "You must be in GRASS GIS to run this program."
		sys.exit(1)

	hsm = options['hsm']
	loc = options['localities']
	fst = options['fst_matrix']
	f_max = options['f_max']
	pval = options['pvalue_matrix']
	p_max = options['p_max']
	corr = options['corridors']
	friction = options['friction']
	lcp_merged = options['lcp_merged']

	if not hsm:
		grass.fatal("Input habitat suitability raster is required")
	if not loc:
		grass.fatal("Input localities raster is required")
	if not fst:
		grass.fatal("Input Fst matrix is required")
	
	# Work starts here
	create_friction(hsm, friction)
	loc_list = create_localities(loc)
	cost_count = create_cost_maps(loc_list, friction)
	grass.message("Created: "+cost_count+" least cost maps")
	pairs=create_fst_pairs(fst, f_max, pval, p_max, len(loc_list))
	grass.message("Created list of: "+str(len(pairs))+" Fst pairs ")
	lcp_count = create_lcp_corridors(pairs, loc_list, lcp_merged)
	grass.message("Created: "+lcp_count+" Least Cost Path reclass maps")


if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())
