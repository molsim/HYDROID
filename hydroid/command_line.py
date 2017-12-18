#!/usr/bin/env python

import pkgutil
import os
from Bio.Seq import Seq
import tempfile


def test_exp():

	from hydroid.HYDROIDexp import assign_peaks_interactive

	temp = tempfile.NamedTemporaryFile()
	temp2 = tempfile.NamedTemporaryFile()
	temp.write(pkgutil.get_data('hydroid', 'pkgdata/test_data/lane_profiles.xls'))
	temp2.write(pkgutil.get_data('hydroid', 'pkgdata/test_data/lane_config.csv'))
	temp.seek(0)
	temp2.seek(0)
	lane_profile_file = temp.name
	lane_config_file = temp2.name

	lane_names=['scCSE4_601TA_BS']

	#Common code:
	# will iterate through the lanes opening an interactive window which will show
	# results of automatic peak identification.
	# Parameters of the automatic peak identification should be interactively adjusted
	# until locations of all peaks are identified correctly
	# This mainly applies to the OH-footprinting profiles that will be latter quantified,
	# sequencing profiles might be left as is.
	###################################
	for LN in lane_names:
		assign_peaks_interactive(lane_profile_file,lane_config_file,LN)
	temp.close()
	temp2.close()


def test_pred():

	hydroid.HYDROIDpred import get_DNA_H_SASA

	out_path="results"
	temp = tempfile.NamedTemporaryFile()
	temp.write(pkgutil.get_data('hydroid', 'pkgdata/test_data/test.pdb'))
	temp.seek(0)

	try:
		os.mkdir(out_path)
	except:
		pass

	#Set DNA sequence
	TS_seq=Seq('GGAGATACCCGGTGCTAAGGCCGCTTAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCTACCGCGTTTTAACCGCCAATAGGATTACTTACTAGTCTCTAGGCACGTGTAAGATATATACATCCT')
	BS_seq=TS_seq.reverse_complement()

	prof_data=[
		{'prof_name':'test','pdb_file':temp.name,'chain':'I','resids':range(-71,72),'seq':TS_seq,'vdw_set':'charmm36-rmin'}
		]


	#Common code:
	#will calculate H-SASA profiles from PDB structure with hydrogen atoms.
	#adjust n_threads to the number of CPU cores in your computer for optimal performance
	###################################
	for p in prof_data:
		get_DNA_H_SASA(p['pdb_file'],os.path.join(out_path,p['prof_name']+'_H-SASA.csv'),\
			chain=p['chain'],resids=p['resids'],seq=p['seq'],probe_radius=1.4,slicen=200,vdw_set=p['vdw_set'],\
			Hcontrib=[1.0]*7,n_threads=6,verbose=False)

	temp.close()



