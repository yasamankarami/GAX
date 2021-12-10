#!/usr/bin/env python3
# -*- coding: UTF8 -*-
# Author: Yasaman Karami -- yasaman.karami@pasteur.fr
# 2021-11-18

import os
import sys
import argparse
import numpy
import MDAnalysis
import os.path 
import time
import genetic_algorithm
import pickle

ADDR_RAM="/dev/shm"
"""
python GA_wrapper.py --addr /c7/scratch2/ykarami/MD/GA/code/data \
--pdb diubqt.pdb --trj diubqt_prot_fit.xtc --saxs SASDCG7.dat
"""
def split_trajectory(addr, traj_file_name, protein_file_name):
	u = MDAnalysis.Universe("%s/%s" %(addr, protein_file_name), "%s/%s" %(addr, traj_file_name))
	protein = u.select_atoms("protein")
	# get index first and last atom
	atm0 = protein[0].index  + 1
	atmN = protein[-1].index + 1 
	nb_atoms = u.atoms.n_atoms 
	nb_frames = len(u.trajectory)
	#########################################################
	counter = 0
	for ts in u.trajectory[0:nb_frames]:
		with MDAnalysis.Writer("%s/protein%d.pdb" %(ADDR_RAM,counter), protein.n_atoms) as W:
			W.write(protein)
		counter += 1
	#########################################################
	'''
	chuncks = 1000
	nchunks = int(float(nb_frames)/chuncks+0.9999999)
	for j in range(nchunks):
		start1 = j * chuncks
		end1 = (j+1) * chuncks
		if end1 > nb_frames:
			end1 = nb_frames
		with MDAnalysis.Writer("%s/protein%d.pdb" %(ADDR_RAM,j), protein.n_atoms) as W:
			for ts in u.trajectory[start1:end1]:
				W.write(protein)
	'''
	return nb_atoms, nb_frames#, nchunks

def read_saxs_file(addr, saxs_file_name):
	saxs_file = open("%s/%s" %(addr, saxs_file_name), "r").readlines()
	Qval = []; Ival = []
	for i in range(len(saxs_file)):
		if saxs_file[i].startswith("#"): continue
		if len(saxs_file[i].split()) != 3: continue
		Qval.append(float(saxs_file[i].split()[0]))
		Ival.append(float(saxs_file[i].split()[1]))
	return Qval, Ival

def run_foxs(addr, frames, saxs_id, saxs_file):
	for i in range(frames):
		os.system("foxs %s/protein%d.pdb %s/%s -h -m 1" %(ADDR_RAM, i, addr, saxs_file))
		os.system("rm %s/protein%d.pdb" %(ADDR_RAM, i))
		os.system("rm %s/protein%d_%s.dat" %(ADDR_RAM, i, saxs_id))
		os.system("rm %s/protein%d.pdb.dat" %(ADDR_RAM, i))
		#os.system("rm %s/protein%d_m*_%s.dat" %(ADDR_RAM, i, saxs_id))
		#os.system("rm %s/protein%d_m*.pdb.dat" %(ADDR_RAM, i))
	#file_name = "%s/protein%d_%s.fit" %(ADDR_RAM, (frames-1), saxs_id)
	#if not os.path.isfile(file_name):
	#	raise ValueError("FoXs failed!!! %s does not exist!" % file_name)
	return

def read_profile(addr, saxs_name, nbFrames):
	out_file = open("%s/Icalc.dat" %addr, "w")
	for i in range(nbFrames):
		in_file = open("%s/protein%d_%s.fit" %(ADDR_RAM, i, saxs_name), "r").readlines()
		for lines in in_file:
			if lines.startswith("#"): continue
			out_file.write("%s " %lines.split()[3])
		os.system("rm %s/protein%d_%s.fit" %(ADDR_RAM, i, saxs_name))
		out_file.write("\n")
	out_file.close()
	return

def print_results(addr, ensemble_size, nbRep):
	report_file = open("%s/GA_results.dat" %addr, "w")
	report_file.write("#ensemble_size replicate chi2 frames weights\n")
	for ii in range(len(ensemble_size)):
		size = ensemble_size[ii]
		for rep in range(nbRep):
			file_name = '%s/ga_saxs_1000_%d_%d.dat' %(addr, size, rep)
			if not os.path.exists(file_name): continue
			res_file = open(file_name, 'rb')
			ga = pickle.load(res_file)
			score = ga.score[0]
			frames = ga.component_ids[0]
			weights = ga.weights[0]
			chi2 = float(ga.score[0])
			report_file.write("%d %d %f " %(size, rep+1, chi2))
			for j in range(len(frames)):
				report_file.write("%d " %frames[j])
			for j in range(len(weights)):
				report_file.write("%f " %weights[j])
			report_file.write("\n")
	report_file.close()

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Extract ensemble of conformations matching best with the experimental SAXS data')
	parser.add_argument('--addr', type=str, help='Working directory',
	                    required=True)
	parser.add_argument('--pdb', type=str, help='Protein structure file',
	                    required=True)
	parser.add_argument('--trj', type=str, help='Trajectory file',
	                    required=True)
	parser.add_argument('--saxs', type=str, help='Experimental SAXS data',
	                    required=True)
	parser.add_argument('-s', help='Ensemble size', nargs='+', type=int,
	                    default=3)
	parser.add_argument('-r', help='Number of repeats', type=int,
	                    default=3)
	args = parser.parse_args()

	PATH = args.addr
	PDBFILENAME = args.pdb
	TRAJECTORYNAME = args.trj
	SAXSEXPNAME = args.saxs
	saxs_id = SAXSEXPNAME.split(".")[0]
	ENSEMBLE_SIZE = args.s
	NUMBER_OF_REPEATS = int(args.r)

	Qvalues, Ivalues = read_saxs_file(PATH, SAXSEXPNAME)
	### divide the trajectory into smaller chunks
	nb_atoms,nb_frames = split_trajectory(PATH, TRAJECTORYNAME, PDBFILENAME)
	print("pdb extracted!")
	### run FoXs to back calculate the profiles
	run_foxs(PATH, nb_frames, saxs_id, SAXSEXPNAME)
	print("FoXs done!")
	### reading the calculated profiles into a single file
	read_profile(PATH, saxs_id, nb_frames)
	### running the genetic algorithm
	for i in range(len(ENSEMBLE_SIZE)):
		ens_size = ENSEMBLE_SIZE[i]
		for j in range(NUMBER_OF_REPEATS):
			experiment_labels = ['saxs']
			n_generation = 1000
			n_ensemble = 1000
			window_size = 100
			rsd_stop = 0.0001
			experimental_weights = [1.0]
			score_types = ['chi2']
			exp_saxs = "%s/%s" %(PATH, SAXSEXPNAME)
			calc_saxs = "%s/Icalc.dat" %PATH
			#genetic_algorithm.run_ga_arg(experiment_labels, exp_saxs, calc_saxs, n_generation, n_ensemble, ens_size, experimental_weights, score_types, j, window_size, rsd_stop,PATH)
	### report the results
	print_results(PATH, ENSEMBLE_SIZE, NUMBER_OF_REPEATS)
	