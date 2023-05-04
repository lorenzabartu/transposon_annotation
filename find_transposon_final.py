# find_transposon.py
# Written by Lorenza Bartu and Shahd ElNaggar
# Computational Genomics Final Project
# Input a bacterial DNA fasta file (required) to find predicted transposons/insertion sequences

import csv
import sys
import argparse
import os 
import pandas as pd
from glob import glob
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIWWW
from Bio import GenBank
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from dna_features_viewer import GraphicFeature, GraphicRecord
from dna_features_viewer import BiopythonTranslator
import subprocess
from matplotlib import pyplot as plt


# command-line arguments and input files
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--plasmid", type=str, required = True)
parser.add_argument("-o", "--outdir", metavar = "directory")
parser.add_argument("-i", "--identity", metavar = "%", type = float)
parser.add_argument("-c","--context",action="store_true")
parser.add_argument("-m", "--gbkmap",action="store_true")
parser.set_defaults(output = os.getcwd(),identity = 85)
plasname = (parser.parse_args()).plasmid
output = (parser.parse_args()).outdir
ident = (parser.parse_args()).identity

# blast sequence against database
def run_blastn(filename, outdir): 
	outfile = '.'.join(os.path.basename(filename).split('.')[:-1]) + ".blastn"
	if outfile in os.listdir(outdir + "/blastn"):
		print("\tPrevious analysis found for " + outfile + ". Passing.")
		pass

    # Run BLASTn if no .blastn file is found in the output directory.
	else:
		db = os.path.dirname(os.path.realpath(__file__)) + "/db/transposons.fna"
		print("\tBlasting file " + filename + '.')
		subprocess.call(["blastn",
						"-query", filename,
						"-subject", db,
						"-out", outdir + "/blastn/" + outfile,
                        "-outfmt", "6 qseqid qstart qend sseqid slen sstart send pident length evalue",
                        "-dust", "no",
                        "-word_size", "7",
                        "-evalue", "1e-5"])

	
	parse_results(outdir + "/blastn/" + outfile)


# parse blast results 
def parse_results(filename):

	plasmid = plasname.strip('.fasta')
	blast_file = pd.read_table(filename, header = None)
	blast_file.columns = ["Accession_ID", "Start_Query", "Stop_Query", "Subject_ID", 
							"Length_Subject", "Start_Sub", "Stop_Sub", "Percent_Identity", "Hit_Length", "E_value"]
	blast_file["Percent_Coverage"] = blast_file["Hit_Length"]*100/blast_file["Length_Subject"]

	# sort table by e-value
	blast_file = blast_file.sort_values(by = "E_value", ascending = True)

	# drop duplicate IDs for transposons at the same start site and filter by pct identity  
	blast_file = blast_file.drop(columns = ["Accession_ID"]).drop_duplicates(subset = ["Start_Query"]).drop_duplicates(subset = ["Stop_Query"])
	results_filtered = blast_file[blast_file["Percent_Identity"] >= ident]

	# update ID name 
	results_filtered = results_filtered.reset_index(drop = True)
	counter = 0
	for i in results_filtered["Subject_ID"]:
		i = i.split('|')
		results_filtered.at[counter, "Subject_ID"] = i[1] + "|" + i[2] 
		counter = counter + 1

	# sort results by start query 
	results_filtered = results_filtered.sort_values(by = ['Start_Query'])
	results_filtered = results_filtered.reset_index(drop = True)	

	# creating dictionary to store all MGE
	match = {}
	matches = []
	seq_match = 0
	match_number = 1
	match_name = "MGE_"

	for index, row in results_filtered.iterrows():
		if seq_match == 0: 
			seq_match = row["Stop_Query"]
			matches.append(row["Start_Query"])
			matches.append(row["Stop_Query"])


		else: 
			compare = row["Start_Query"]
			if (compare - seq_match) < 2500: 
				seq_match = row["Stop_Query"]
				matches.append(seq_match)

			else:
				key = match_name + str(match_number)
				match[key] = matches
				matches = []
				matches.append(compare)
				seq_match = row["Stop_Query"]
				matches.append(seq_match)
				match_number += 1

	key = match_name + str(match_number)
	match[key] = matches

	print('***************************************************************')
	print("Here is the number of MGEs we're isolating + their starts + stops")
	print(match)

	def get_group(row):
		for group, region in match.items():
			if row['Start_Query'] in region or row['Stop_Query'] in region:
				return group
		return 'Other'

	# apply the function to each row to get the group
	results_filtered['Group'] = results_filtered.apply(get_group, axis=1)

	final_transposons = pd.DataFrame(columns = ["Subject_ID", "Start", "Stop", "PI"])

	# group the dataframe by the groups
	entry = []
	grouped_df = results_filtered.sort_values(by = "Start_Query", ascending = True).groupby('Group')
	for key, item in grouped_df: 
		g_mge = item.groupby('Subject_ID')
		for name, group in g_mge: 
			entry.append(name)
			entry.append(group.get("Start_Query").values[0])
			entry.append(group.get("Stop_Query").values[len(group)-1])
			entry.append(group.get("Percent_Identity").values.mean())
			final_transposons.loc[len(final_transposons)] = entry
			entry =[]

	# new_res = results_filtered.groupby(by=["Subject_ID"])
	final_transposons = final_transposons.sort_values(["Start", "Stop"])
	removing = []

	current = final_transposons.iloc[0]
	for i in range(1,len(final_transposons)):
		next_row = final_transposons.iloc[i]
		if next_row["Start"] >= current["Stop"]:
			current = next_row
		elif next_row["PI"] < current["PI"]:
			removing.append(i)

	final_transposons = final_transposons.drop(removing)
	grouped = final_transposons.groupby("Subject_ID")
	missed_tn = pd.DataFrame(columns=["Subject_ID", "Start", "Stop", "PI"])
	for key,item in grouped:
		if len(item) == 2:
			entry.append(key)
			entry.append(item.get("Start").values[0])
			entry.append(item.get("Stop").values[len(item)-1])
			entry.append(item.get("PI").values.mean())
			missed_tn.loc[len(missed_tn)] = entry
			entry = []

	if not missed_tn.empty:
		for row_index, row in missed_tn.iterrows():
			name = row["Subject_ID"]
			final_transposons = final_transposons[final_transposons.Subject_ID != name]
			final_transposons.loc[len(final_transposons)] = row

	print(final_transposons)

	for i in match: # match or final_trasposons?
		rec = (match[i])
		rec = [i, rec[0], rec[-1]]
		if parser.parse_args().context:
			genetic_context_files(rec, plasmid)
			genetic_context(plasmid)
			genetic_context_df = pd.DataFrame(columns = ['Gene', 'Start', 'Stop']) # edit 

	final_transposons.to_csv(output + "/final.csv")
	if parser.parse_args().gbkmap:
		gbk_map(plasname,final_transposons)
	return final_transposons


# writes fasta files for regions upstream and downstream of transposon
def genetic_context_files(info, plasmid):
	name = info[0]
	start_mge = info[1]
	stop_mge = info[2]
	rec = SeqIO.parse(plasmid + ".fasta", "fasta")
	context = []
	for record in rec: 
		if start_mge < 1000:
			prior = record.seq[start_mge-1000:]+ record.seq[0:start_mge]
			prior = SeqRecord(seq = Seq(prior), id = "prior_"+str(start_mge))
			SeqIO.write(prior, "./"+plasmid+"_results/context/"+name+"_prior_context.fasta", "fasta")

		else: 
			prior = record.seq[start_mge-1000:start_mge]
			prior = SeqRecord(seq = Seq(prior), id = "prior_"+str(start_mge))
			SeqIO.write(prior, "./"+plasmid+"_results/context/"+name+"_prior_context.fasta", "fasta")

		if (len(record.seq) - stop_mge) < 1000: 
			extra = 1000 - len(record.seq) - stop_mge
			post = record.seq[stop_mge:] + record.seq[0:extra]
			post = SeqRecord(seq = Seq(post), id = "post_"+str(stop_mge))
			SeqIO.write(post, "./"+plasmid+"_results/context/"+name+"_post_context.fasta", "fasta")

		else: 
			post = record.seq[stop_mge: stop_mge+1000]
			post = SeqRecord(seq = Seq(pos,final_transposonst), id = "post_"+str(stop_mge))
			SeqIO.write(post, "./"+plasmid+"_results/context/"+name+"_post_context.fasta", "fasta")

	genetic_context(plasmid)


# protein blastx of upstream and downstream - is gene being interrupted? 
def genetic_context(plasmid):
	path = "./"+plasmid+"_results/context/"
	fastas = glob("%s/*.fasta" % path)
	for file in fastas: 
		print("Currently blasting " + file)
		name = (file[file.index("MGE"):-6])
		if path+name in os.listdir(output + "/blastn"):
			print("\tPrevious analysis found for " + outfile + ". Passing.")
			pass
		else:
			blastx_cline = NcbiblastxCommandline(query= file, db='refseq_protein', evalue=0.000001, outfmt=5,
	                                          out=path+name + ".xml", remote=True, max_target_seqs = 4)
			stdout, stderr = blastx_cline()

		result_handle = open(path+name +".xml")
		blast_records = NCBIXML.parse(result_handle)
		# to be continued - speed of blast must be improved in future 

# makes a genbank file of the features and then plots them as a map
def gbk_map(plasname,transposons):
	for record in SeqIO.parse(plasname,"fasta"):
		sequence = record.seq
		description = record.description


	# make genbank file 
	gbk_record = SeqRecord(Seq(sequence), 
                    description=description+" predicted composite transposons",
                    annotations={"accession":'.', "version":'.',
                                "organism":'.', "date":'today',
                                "data_file_division":"BCT",
                                "molecule_type": "DNA"})
	# add features
	map_feats = []
	for index, row in transposons.iterrows():
		if (row["Start"] < row["Stop"]):
			curr_strand = +1
		else:
			curr_strand = -1
		feature = SeqFeature(FeatureLocation(
	            row["Start"], row["Stop"], strand=curr_strand),
	            type="MGE",
	            qualifiers={"IS": row["Subject_ID"]})
		current = GraphicFeature(start = row["Start"], end = row["Stop"],
				strand = curr_strand, color = "#ffd700", label = row["Subject_ID"])
		map_feats.append(current)
		gbk_record.features.append(feature)
	outname = os.path.dirname(os.path.realpath(__file__)) + "/" + output + "/" + plasname.strip(".fasta") + ".gb"
	SeqIO.write(gbk_record, outname , "gb")
	positions = transposons["Stop"].append(transposons["Start"])
	seq_length = positions.max() + 500
	record = GraphicRecord(sequence_length=seq_length, features=map_feats)
	myplot = record.plot(figure_width = 5)
	plt.show()

def main():
	print("\tStarting analyses.")
	if not os.path.exists(output):
		os.mkdir(output)
	if not os.path.exists(output + "/blastn"):
		os.mkdir(output + "/blastn")
	run_blastn(plasname,output)

if __name__ == "__main__":
	main()