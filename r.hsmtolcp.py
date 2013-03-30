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

#%flag
#% key: e
#% description: Export the localities as shapefile and resulting Least Cost Network as GeoTiff
#% end
#%option
#% key: maxent
#% type: string
#% description: Name of the input Habitat Suitability Model raster (in ArcInfo ASCII format) from MaxEnt
#% required: yes
#%end
#%option
#% key: friction
#% type: string
#% description: Name for GRASS friction raster
#% required: no
#% answer: friction
#%end
#%option
#% key: fst_matrix
#% type: string
#% description: Input text file (csv) of Fst matrix from arlequin output
#% required: yes
#%end 
#%option
#% key: f_max
#% type: double
#% description: Maximum Fst value to consider a pair of localities as "close" populations genetically
#% required: no
#% answer: 0
#%end 
#%option
#% key: pvalue_matrix
#% type: string
#% description: Input text file (csv) of p-values matrix from arlequin output
#% required: yes
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
#% description: Input text file (csv) sampling localities, with locality_code, locality_num, and x,y coordinates
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
#% key: lcp_network
#% type: string
#% description: Name for Least Cost Path Network (merged) output raster
#% required: no
#% answer: lcp_network

import os, sys, csv
import grass.script as grass
import numpy as np


def create_friction(maxent, friction):
	"""
	Import hsm raster from MaxEnt output into GRASS and create friction map
	"""
	grass.message(" === Reading MaxEnt input raster to create friction map ===")
	# Import maxent as GRASS raster named "hsm"
	grass.run_command('r.in.gdal', input=maxent, output="hsm", overwrite=True, quiet=True, flags="o")
	r=grass.read_command('g.region', rast="hsm", flags="p")
	#grass.message("Current region settings: "+r)

	# Create friction map
	exists = grass.find_file(name = friction, element='cell')
	if exists:
		grass.message("friction raster: "+friction +" already exists and will be overwritten")
	
	friction_expr = friction+" = 1.0-hsm"
	#print "Running mapcalc: " + friction_expr
	grass.mapcalc(friction_expr, overwrite=True)



def create_localities(loc):
	"""
	Read and parse the csv file of localities
	The csv file must have exactly 5 folumns:
	"Locality name" (a string), index (integer), code (3-4 letters), X coord, Y coord (floats)
	Create a python list, where each entry is a list of:
	('Locality Name',Locality index,'Locality code, X coord, Y coord)
	Also Create a GRASS vector
	"""
	try:
		csv_file = open(loc,"rb")
	except IOError:
		grass.fatal("Cannot open localities file: "+loc)

	grass.message(" === Reading list of localities ===")
	locs=csv.reader(csv_file, delimiter=',')
	loc_list=[]
	for l in locs:
		# Add to the loc_list a tuple containing the code, X coord and Y coord
		loc_list.append([l[0],l[1],l[2],l[3],l[4]])

	csv_file.close()

	# Make a GRASS vector
	loc_vector = os.path.splitext(loc)[0]
	grass.run_command('v.in.ascii', input=loc, output=loc_vector, separator=",", x=4, y=5, 
			columns="loc_name varchar(32), id integer, code varchar(6), longitude double, latitude double", 
			quiet=True, overwrite=True)
	
	return loc_list, loc_vector




def create_fst_pairs(fst, f_max, pval, p_max, num_locs):
	"""
	Use NumPy to import Fst and p-value matrices
	Use boolean comparisons to find which Fst values are <= f_max
	and which p-value values are <= p_max
	Create a list of pairs (of localities) which match both conditions
	"""

	grass.message(" === Creating list of pairs from Fst and p-value matrices ===")
	# Import arrays, skip header row, and slice off first column
	fst_mat = np.genfromtxt(fst, skip_header=1, delimiter=',', usecols=range(1,num_locs+1))
	pval_mat = np.genfromtxt(pval, skip_header=1, delimiter=',', usecols=range(1,num_locs+1))
	# Make boolean arrays that reflect the max conditions for Fst and p-value 
	fst_bool = (fst_mat<=f_max)
	#print fst_bool
	pval_bool = (pval_mat<=p_max)
	#print pval_bool
	# Create array where both conditions apply (logical AND)
	acc_bool = np.logical_and(fst_bool, pval_bool)
	# Get indices where accept matrix is TRUE
	accept = np.where(acc_bool)
	# And convert to list
	pairs = np.transpose(accept).tolist()
	#print pairs

	return pairs


def create_lcp_corridors(friction, pairs, locs, lcp_network):
	"""
	Loop thru all pairs of localities in the pairs list
	FOr each pair, find the indexes, codes and X-Y coords of that pair from the localities list,
	Use r.cost to make corridor maps using the start_coordinate and end_coordinate paramters
	Get the minimum for each corridor using r.univar, 
	Use that minimum to make a reclass file with values:
		min					= 5			Min value
		min + 5% 		= 3			Two percent above
		min + 8%		=	1			Five percent above
		*						= NULL	NULL
	Run r.reclass to create uniform reclass rasters for each pair
	"""
	grass.message(" === Creating corridors ===")
	cnt = 0
	for i in range(len(pairs)):
		id1, id2 = pairs[i][0], pairs[i][1]
		this_pair = str(id1)+","+str(id2)
		for j in range(len(locs)):
			# Find the locality code that matches each index id1 and id2
			# Get X-Y coords for each locality
			if int(locs[j][1]) == id1:
				code1 = locs[j][2]
				x1, y1 = locs[j][3], locs[j][4]
			elif int(locs[j][1]) == id2:
				code2 = locs[j][2]
				x2, y2 = locs[j][3], locs[j][4]

		# Some pairs might not be in the sampling sites list at all.
		# Check if both codes are found before continuing
		try:
			code1 and code2
		except NameError:
			grass.message("The pair: "+this_pair+" is not in the localities list. Ignoring...")
		else:
			# Create both cost maps (one for each locality in the pair)
			cnt += 1
			cost1, cost2  = "cost_"+code1, "cost_"+code2
			start1=x1+","+y1
			start2=x2+","+y2
			grass.run_command('r.cost', input=friction, output=cost1, 
				start_coordinate=start1, overwrite=True, quiet=True)
			grass.run_command('r.cost', input=friction, output=cost2, 
				start_coordinate=start2,  overwrite=True, quiet=True)

			# Now combine the cost maps into a corridor
			corridor = "corr_"+code1+"_"+code2
			corridor_expr = corridor+"=round("+cost1+"+"+cost2+")"
			grass.mapcalc(corridor_expr, overwrite=True)
			# Get minimum value from the corridor map
			u = grass.read_command('r.univar', map=corridor, flags="g", quiet=True)
			udict = grass.parse_key_val(u)
			min = float(udict['min'])
			five_pc = min*1.05
			eight_pc = min*1.08
			# Create reclass file
			tmp_reclass = grass.tempfile()
			trc = open(tmp_reclass, "w")
			# Convert to int for output to the reclass file
			trc.write(str(int(min)) + " = 5	\t Minimum\n")
			trc.write(str(int(five_pc)) + " = 3 \t Five percent\n")
			trc.write(str(int(eight_pc)) + " = 1 \t Eight percent\n")
			trc.write("* = NULL \t NULL\n")
			trc.close()
			# Now create reclass raster
			lcp = "lcp_"+code1+"_"+code2
			grass.run_command('r.reclass',input=corridor,output=lcp, rule=tmp_reclass, overwrite=True, quiet=True)
			grass.run_command('r.colors', map=lcp, color="ryg", quiet=True)
			# Get rid of tmp reclass file and cost raster
			os.unlink(tmp_reclass)
			grass.run_command('g.mremove', rast="cost_*", quiet=True, flags="f")

	# Now merge all lcp_* maps
	lcp_maps = grass.read_command('g.mlist', type="rast", pattern="lcp_*", separator=",").rstrip()
	grass.run_command('r.series',input = lcp_maps, output = lcp_network, method="sum", overwrite=True, quiet=True)
	grass.run_command('r.colors', map=lcp_network, color="ryg", quiet=True)
	return cnt


def export_layers(vect, rast):
	"""
	Export the localities vector to a shapefile
	and the lcp_network to a GeoTiff
	"""
	grass.message(" === Exporting localities shapefile and lcp_network raster ===")
	shp = vect+".shp"
	gtiff=rast+".tif"
	if os.path.isfile(shp):
		os.unlink(shp)
		os.unlink(vect+".dbf")
		os.unlink(vect+".shx")
	if os.path.isfile(gtiff):
		os.unlink(gtiff)
	
	grass.run_command('v.out.ogr',input=vect, dsn=shp, overwrite=True, quiet=True, flags="e")
	grass.run_command('r.out.gdal', input=rast, output=gtiff, 
			overwrite=True, quiet=True, type="Float64", createopt="TFW=YES")
	
	return shp, gtiff

def main():
	""" 
	Get command line parameters
	"""
	if "GISBASE" not in os.environ:
		print "You must be in GRASS GIS to run this program."
		sys.exit(1)

	maxent = options['maxent']
	loc = options['localities']
	fst = options['fst_matrix']
	f_max = float(options['f_max'])
	pval = options['pvalue_matrix']
	p_max = float(options['p_max'])
	corr = options['corridors']
	friction = options['friction']
	lcp_network = options['lcp_network']
	export_bool = flags['e']

	if not maxent:
		grass.fatal("Input habitat suitability raster is required")
	if not loc:
		grass.fatal("Input localities raster is required")
	if not fst:
		grass.fatal("Input Fst matrix is required")
	
	# Work starts here
	create_friction(maxent, friction)
	loc_list, loc_vector = create_localities(loc)
	#cost_count = create_cost_maps(loc_list, friction)
	#grass.message("Created: "+str(cost_count)+" least cost maps")
	pairs=create_fst_pairs(fst, f_max, pval, p_max, len(loc_list))
	grass.message(" === Created list of: "+str(len(pairs))+" Fst pairs ===")
	lcp_count = create_lcp_corridors(friction, pairs, loc_list, lcp_network)
	grass.message(" === Created: "+str(lcp_count)+" Least Cost Path reclass maps, and merged into: "+lcp_network+" ===")
	if export_bool:
		shp, gtiff = export_layers(loc_vector, lcp_network)
		grass.message(" === Exported layers: " +shp+ " and " + gtiff+" ===")



if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())
