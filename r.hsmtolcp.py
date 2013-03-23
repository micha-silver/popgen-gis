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
	
	friction_expr = "friction = 1.0-hsm"
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
		loc_list.append([l[0],l[1],l[2],l[3]])

	csv_file.close()
	return loc_list


def create_cost_maps(loc_list, friction)
	"""
	Use the friction map to create a least cost cost raster for each of the localities
	based on the friction map and locality coordinates
	The cost rasters are all named "cost"_localitycode
	"""
	cnt=0
	for i in range(len(loc_list)):
		cnt +=1
		code, x, y = loc_list[i][0], loc_list[i][2], loc_list[i][3]
		outrast="cost"+code
		coords=x+","+y
		grass.run_command('r.cost', input=friction, output=outrast, coordinates=coords, overwrite=True, quiet=True)

	return cnt


def create_fst_pairs(fst, f_max, pval, p_maxi, num_locs)
	"""
	Use NumPy to import Fst and p-value matrices
	Use boolean comparisons to find which Fst values are <= f_max
	and which p-value values are <= p_max
	Create a list of pairs (of localities) which match the conditions
	"""
	
	# Import arrays, skip header row, and slice off first column
	fst_mat = np.genfromtxt(fst, skip_headers=1, delimiter=',', usecols=range(1,39))
	pval_mat = np.genfromtxt(pval, skip_headers=1, delimiter=',', usecols=range(1,39))
	fst_accept = fst_mat <= f_max
	pval_accept = pval_mat <= p_max
	# Create array where both conditions apply (logical AND)
	accept = np.logical_and(fst_accept, pval_accept)


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



if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())
