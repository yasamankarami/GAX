#!/usr/bin/env python3
# -*- coding: UTF8 -*-

#############################################################################
# Author: Yasaman Karami -- yasaman.karami@pasteur.fr                       #
# https://research.pasteur.fr/en/member/fr-yasaman-karami/                  #
# Copyright (c) 2022 Institut Pasteur                                       #
#                                                                           #
#                                                                           #
#  Redistribution and use in source and binary forms, with or without       #
#  modification, are permitted provided that the following conditions       #
#  are met:                                                                 #
#                                                                           #
#  1. Redistributions of source code must retain the above copyright        #
#  notice, this list of conditions and the following disclaimer.            #
#  2. Redistributions in binary form must reproduce the above copyright     #
#  notice, this list of conditions and the following disclaimer in the      #
#  documentation and/or other materials provided with the distribution.     #
#  3. Neither the name of the copyright holder nor the names of its         #
#  contributors may be used to endorse or promote products derived from     #
#  this software without specific prior written permission.                 #
#                                                                           #
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS      #
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT        #
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR    #
#  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT     #
#  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   #
#  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT         #
#  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,    #
#  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY    #
#  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      #
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE    #
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.     #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#                                                                           #
#############################################################################

import os,sys
import argparse
import numpy as np
import os.path 
import time
import genetic_algorithm
import pickle
import matplotlib.pyplot as plt
from zipfile import ZipFile
import MDAnalysis

def check_ensemble(ens_list):
	if len(ens_list) == 0 or ens_list.startswith(",") or ens_list.endswith(","):
		print("Error in the list of ensemble size!!!")
	words = ens_list.split(",")
	if len(words) == 0:
		print("Error in the list of ensemble size!!!")
	ens = []
	for i in range(len(words)):
		ens.append(int(words[i]))
	return ens

def check_SAXS_file(work_dir, saxs_file_name):
	new_saxs_id = "SAXS_exp_modif"
	saxs_flag = "no"
	saxs_file = open(saxs_file_name, "r").readlines()
	modif_saxs_file = open("%s/%s.dat" %(work_dir, new_saxs_id), "w")
	nb_saxs_points = 0
	for i in range(len(saxs_file)):
		if saxs_file[i].startswith("#"): continue
		words = saxs_file[i].split()
		if len(words) != 3: 
			print("Error in the experimental SAXS data! Every row should contain 3 values (Q, I, err).")
		if float(words[1]) < 0:
			print("WARNING Negative intensity value: %s %s %s skipping remaining profile points" %(words[0], words[1], words[2]))
			saxs_flag = "yes"
			break
		modif_saxs_file.write(saxs_file[i])
		nb_saxs_points += 1
	return nb_saxs_points, new_saxs_id, saxs_flag

def check_profile_format(profile, nbSAXS):
	inp_file = open(profile,"r").readlines()
	for i in range(len(inp_file)):
		if inp_file[i].startswith("#"): continue
		words = inp_file[i].split()
		if len(words) != nbSAXS:
			print("Back-calculated profile format is not correct!\n The number of of columns should be equal to the number of experimental Intensities and the number of rows should be equal to the number of MD conformations!!!") 
	return

def split_trajectory(traj_file_name, protein_file_name, temp_addr):
	u = MDAnalysis.Universe(protein_file_name, traj_file_name)
	protein = u.select_atoms("protein")
	# get index first and last atom
	atm0 = protein[0].index  + 1
	atmN = protein[-1].index + 1 
	nb_atoms = u.atoms.n_atoms 
	nb_frames = len(u.trajectory)
	counter = 0
	for ts in u.trajectory[0:nb_frames]:
		with MDAnalysis.Writer("%s/protein%d.pdb" %(temp_addr,counter), protein.n_atoms) as W:
			W.write(protein)
		counter += 1
	return nb_atoms, nb_frames

def read_saxs_file(saxs_file_name):
	saxs_file = open(saxs_file_name, "r").readlines()
	Qval = []; Ival = []
	for i in range(len(saxs_file)):
		if saxs_file[i].startswith("#"): continue
		if len(saxs_file[i].split()) != 3: continue
		Qval.append(float(saxs_file[i].split()[0]))
		Ival.append(float(saxs_file[i].split()[1]))
	return Qval, Ival

def run_foxs(frames, saxs_id, saxs_file, temp_addr):
	for i in range(frames):
		os.system("foxs %s/protein%d.pdb %s -h -m 1" %(temp_addr, i, saxs_file))
		os.remove("%s/protein%d.pdb" %(temp_addr, i))
		os.remove("%s/protein%d_%s.dat" %(temp_addr, i, saxs_id))
		os.remove("%s/protein%d.pdb.dat" %(temp_addr, i))
	return

def read_profile(saxs_name, nbFrames, temp_addr):
	Icalc_file_name = "%s/Icalc.dat" %temp_addr 
	out_file = open(Icalc_file_name, "w")
	for i in range(nbFrames):
		in_file = open("%s/protein%d_%s.fit" %(temp_addr, i, saxs_name), "r").readlines()
		for lines in in_file:
			if lines.startswith("#"): continue
			out_file.write("%s " %lines.split()[3])
		os.remove("%s/protein%d_%s.fit" %(temp_addr, i, saxs_name))
		out_file.write("\n")
	out_file.close()
	return Icalc_file_name

def plot_results(addr, ensemble_size, nbRep, nb_generation):
	colors = ['olive', 'purple', 'cyan', 'salmon','darkkhaki', 'darkolivegreen', 'violet', 'green', 'blue', 'red']
	chi2 = np.zeros((nbRep,len(ensemble_size)), dtype=float)
	fig, axs = plt.subplots(len(ensemble_size), 1, sharex=True, squeeze=False)
	fig.suptitle('Ensemble-based back-calculated SAXS profile')
	for i in range(len(ensemble_size)):
		for j in range(nbRep):
			file_name = '%s/ga_saxs_%d_%d_%d.dat' %(addr, nb_generation, ensemble_size[i], j)
			ga1 = pickle.load(open(file_name, 'rb'))
			chi2[j][i] = float(ga1.score[0])
			if j==0:
				axs[i,0].semilogy(ga1.target[0][:, 0], ga1.target[0][:, 1], color='black',linewidth=2.5, label='exp')
			axs[i,0].semilogy(ga1.target[0][:, 0], ga1.get_models(0)[0], color=colors[j],linewidth=1.5, label='rep %d' %(j+1))
		axs[i,0].set_ylabel("I", rotation=90)
		axs[i,0].set_title('Ensemble size %d' %ensemble_size[i], position=(0.5, 0.75))
		legend = axs[i,0].legend(loc='lower left', shadow=False, framealpha=0)
		legend.get_frame()
	axs[i,0].set_xlabel("q")
	plt.savefig("%s/Intensities.jpg" %addr)
	###############################################
	fig, ax = plt.subplots()
	for i in range(nbRep):
		if len(ensemble_size) > 1:
			ax.semilogx(ensemble_size, chi2[i][:], color=colors[i])
		ax.plot(ensemble_size, chi2[i][:], marker="o", markersize=5, label='rep %d' %(i+1), color=colors[i])
	legend = ax.legend(loc='upper right', shadow=False, framealpha=0)
	ax.set_ylabel('chi2')
	ax.set_xlabel('Ensemble size')
	legend.get_frame()
	plt.savefig("%s/Chi2.jpg" %addr)

def print_results(addr, ensemble_size, nbRep, nb_generation, saxs_flag, new_saxs, profile_flag, bck_profile):
	zipObj = ZipFile('GAX_output.zip', 'w')
	report_file = open("%s/GAX_summary.txt" %addr, "w")
	report_file.write("#ensemble_size,replicate,chi2,frames,weights\n")
	for ii in range(len(ensemble_size)):
		size = ensemble_size[ii]
		for rep in range(nbRep):
			file_name = '%s/ga_saxs_%d_%d_%d.dat' %(addr, nb_generation, size, rep)
			if not os.path.exists(file_name): continue
			res_file = open(file_name, 'rb')
			ga = pickle.load(res_file)
			score = ga.score[0]
			frames = ga.component_ids[0]
			weights = ga.weights[0]
			chi2 = float(ga.score[0])
			report_file.write("%d,%d,%f," %(size, rep+1, chi2))
			for j in range(len(frames)):
				report_file.write("%d" %frames[j])
				if j < (len(frames)-1):
					report_file.write(";")
			report_file.write(",")
			for j in range(len(weights)):
				report_file.write("%f" %weights[j])
				if j < (len(frames)-1): 
					report_file.write(";")
			report_file.write("\n")
			zipObj.write(file_name)
			os.remove(file_name)
	if saxs_flag == "yes":
		zipObj.write(new_saxs)
		os.remove(new_saxs)
	if profile_flag == "no":
		zipObj.write(bck_profile)
		os.remove(bck_profile)
	zipObj.close()
	report_file.close()
	os.system("")

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Extract ensemble of conformations matching best with the experimental SAXS data')
	#parser.add_argument('--addr', type=str, help='Working directory',
	#                    required=True)
	parser.add_argument('--pdb', type=str, help='Protein structure file',
	                    required=True)
	parser.add_argument('--saxs', type=str, help='Experimental SAXS data',
	                    required=True)
	parser.add_argument('--mode', type=str, help='Do you have the back-calculated profile',
	                    required=True)
	####### Conditional inputs depending on the chosen mode ##########################
	parser.add_argument('--trj', type=str, help='Trajectory file')
	parser.add_argument('--calc_saxs', type=str, help='Back-calculated SAXS profile')
	##################################################################################
	parser.add_argument('-s', help='Comma separated list of ensemble size', type=str,
	                    default="3")
	parser.add_argument('-r', help='Number of repeats', type=int,
	                    default=3)
	parser.add_argument('-nb_gen', help='Number of generations', type=int,
	                    default=1000)
	parser.add_argument('-nb_ens', help='Number of ensembles', type=int,
	                    default=1000)
	parser.add_argument('-window', help='window size', type=int,
	                    default=100)
	args = parser.parse_args()

	###
	PDBFILENAME = args.pdb
	SAXSEXPNAME = args.saxs
	ENSEMBLE_SIZE_LIST = args.s
	NUMBER_OF_REPEATS = int(args.r)
	NUMBER_OF_GENERATIONS = args.nb_gen
	NUMBER_OF_ENSEMBLES = args.nb_ens
	WINDOW_SIZE = args.window
	PROFILE_OR_TRJ = args.mode
	
	###
	working_dir = "." #args.addr
	#saxs_id = SAXSEXPNAME.split("/")[-1].split(".")[0]
	
	### check the inputs
	ENSEMBLE_SIZE = check_ensemble(ENSEMBLE_SIZE_LIST)
	nbSAXSpoints, saxs_id, saxs_mod_flag = check_SAXS_file(working_dir, SAXSEXPNAME)
	new_saxs_file = "%s/%s.dat" %(working_dir, saxs_id)

	### check input choice
	if PROFILE_OR_TRJ == "yes":
		back_calc_profile = args.calc_saxs
		check_profile_format(back_calc_profile, nbSAXSpoints)
	if PROFILE_OR_TRJ == "no":
		TRAJECTORYNAME = args.trj
		temp_path = "./temp" #"/dev/shm"
		if not os.path.exists(temp_path):
			os.mkdir(temp_path)
		### read SAXS data
		Qvalues, Ivalues = read_saxs_file(new_saxs_file)
		### divide the trajectory into smaller chunks
		nb_atoms,nb_frames = split_trajectory(TRAJECTORYNAME, PDBFILENAME, temp_path)
		### run FoXs to back calculate the profiles
		run_foxs(nb_frames, saxs_id, new_saxs_file, temp_path)
		### reading the calculated profiles into a single file
		back_calc_profile = read_profile(saxs_id, nb_frames, temp_path)
	
	### running the genetic algorithm
	for i in range(len(ENSEMBLE_SIZE)):
		ens_size = ENSEMBLE_SIZE[i]
		for j in range(NUMBER_OF_REPEATS):
			experiment_labels = ['saxs']
			rsd_stop = 0.0001
			experimental_weights = [1.0]
			score_types = ['chi2']
			genetic_algorithm.run_ga_arg(experiment_labels, new_saxs_file, back_calc_profile, NUMBER_OF_GENERATIONS, NUMBER_OF_ENSEMBLES, ens_size, experimental_weights, score_types, j, WINDOW_SIZE, rsd_stop, working_dir)

	### remove temporary files
	if PROFILE_OR_TRJ == "no":
		new_profile = "%s/Icalc.dat" %working_dir
		os.system("cp %s %s" %(back_calc_profile, new_profile))
		os.remove(back_calc_profile)
		os.rmdir(temp_path)
	else:
		new_profile = back_calc_profile

	### report the results
	plot_results(working_dir, ENSEMBLE_SIZE, NUMBER_OF_REPEATS, NUMBER_OF_ENSEMBLES)
	print_results(working_dir, ENSEMBLE_SIZE, NUMBER_OF_REPEATS, NUMBER_OF_ENSEMBLES, saxs_mod_flag, new_saxs_file, PROFILE_OR_TRJ, new_profile)
	
