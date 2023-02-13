import numpy as np
from itertools import groupby
import json
import sys
import csv
import os
import pandas as pd
from tabulate import tabulate

#pd.set_option('display.float_format', lambda x: '%,g' % x)
with open(sys.argv[2]) as f:
	ff=f.readlines()

# estSizeFile = open()
#
# size_bp = estSizeFile.readlines()

estSize=int(float(sys.argv[4]))

if estSize == 0:
	size_bp = ff[0][:-1]
	print("No estimated Genome Size provided, using estimation from Genomescope2: ", size_bp)
else:
	print("Estimated Genome Size provided as:", estSize, " using this size for NG and LG stats")
	size_bp = estSize

def fasta_iter(fasta_file):
	"""
	Takes a FASTA file, and produces a generator of Header and Sequences.
	This is a memory-efficient way of analyzing a FASTA files -- without
	reading the entire file into memory.

	Parameters
	----------
	fasta_file : str
	The file location of the FASTA file

	Returns
	-------
	header: str
	The string contained in the header portion of the sequence record
	(everything after the '>')
	seq: str
	The sequence portion of the sequence record
	"""

	fh = open(fasta_file)
	fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in fa_iter:
		# drop the ">"
		header = next(header)[1:].strip()
		# join all sequence lines to one.
		seq = "".join(s.upper().strip() for s in next(fa_iter))
		yield header, seq


def read_genome(fasta_file):
	"""
	Takes a FASTA file, and produces 2 lists of sequence lengths. It also
	calculates the GC Content, since this is the only statistic that is not
	calculated based on sequence lengths.

	Parameters
	----------
	fasta_file : str
		The file location of the FASTA file

	Returns
	------
	contig_lens: list
		A list of lengths of all contigs in the genome.
	scaffold_lens: list
		A list of lengths of all scaffolds in the genome.
	gc_cont: float
	The percentage of total basepairs in the genome that are either G or C.
	"""

	gc = 0
	total_len = 0
	contig_lens = []
	scaffold_lens = []
	for _, seq in fasta_iter(fasta_file):
		scaffold_lens.append(len(seq))
		if "NN" in seq:
			contig_list = seq.split("NN")
		else:
			contig_list = [seq]
		for contig in contig_list:
			if len(contig):
				gc += contig.count('G') + contig.count('C')
				total_len += len(contig)
				contig_lens.append(len(contig))
				gc_cont = (gc / total_len) * 100
#	print("this is gc_content", gc_cont)
	return contig_lens, scaffold_lens, gc_cont


def calculate_stats(seq_lens, gc_cont):
	naming = sys.argv[3]
	stats = {}
	seq_array = np.array(seq_lens)
#    stats['Assembly:']=naming
	stats['sequence_count'] = seq_array.size
	testsize=stats['sequence_count']
	stats['number_of_gaps'] = 0
#    print("this is the count",naming," ", testsize)
	stats['gc_content (%)'] = gc_cont
	sorted_lens = seq_array[np.argsort(-seq_array)]
	stats['longest (bp)'] = int(sorted_lens[0])
	testlongest= stats['longest (bp)']
#    print("this is the longest", naming," ",testlongest)
	stats['shortest (bp)'] = int(sorted_lens[-1])
#    stats['median'] = np.median(sorted_lens)
#    stats['mean'] = np.mean(sorted_lens)
	stats['total_bps (bp)'] = int(np.sum(sorted_lens))
	testprint=stats['total_bps (bp)']
#    print("total_bp is", naming," ",testprint)
	stats['estimated_size (bp)'] = int(size_bp)
	csum = np.cumsum(sorted_lens)
	# if stats['total_bps (bp)'] < stats['estimated_size (bp)']:
    #     csum_ng = np.append(csum, stats['estimated_size (bp)'])
    # else:
    #     csum_ng=csum
	for level in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 96, 97, 98, 99, 100]:
		nx = int(stats['total_bps (bp)'] * (level / 100))
		csumn = min(csum[csum >= nx])
		l_level = int(np.where(csum == csumn)[0]) + 1
		stats['L' + str(level)] = l_level
	for level in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 96, 97, 98, 99, 100]:
#        print("the totalbps are:", stats['total_bps (bp)'])
		nx = int(stats['total_bps (bp)'] * (level / 100))
#        print("this is the nx", nx)
#        print("this is the csum", csum)
		csumn = min(csum[csum >= nx])
#       print("this is the csumn", csumn)
		l_level = int(np.where(csum == csumn)[0])
		n_level = int(sorted_lens[l_level])
		stats['N' + str(level)] = n_level
#        print(level, " ", n_level)
	for level in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 96, 97, 98, 99, 100]:
#        print("the estbps are:", stats['estimated_size (bp)'])
		ngx = int(stats['estimated_size (bp)'] * (level / 100))
#        print("this is the ngx", ngx)
#        print("this is the csum", csum_ng)
#        print("this is the csum", csum)
#        print("this is the [csum >= ngx]", np.array(csum >= ngx))
		new_array=np.array(csum >= ngx)
 #       print(np.any(new_array))
		if np.any(new_array) == False:
			csumng = csum[seq_array.size-1]
  #         print("this is the csumng", csumng)
			lg_level = int(np.where(csum == csumng)[0]) + 1
			stats['LG' + str(level)] = lg_level
		elif np.any(new_array) == True:
			csumng = min(csum[csum >= ngx])
  #         print("this is the csumng", csumng)
			lg_level = int(np.where(csum == csumng)[0]) + 1
			stats['LG' + str(level)] = lg_level
	for level in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 96, 97, 98, 99, 100]:
		ngx = int(stats['estimated_size (bp)'] * (level / 100))
		new_array=np.array(csum >= ngx)
 #       print(np.any(new_array))
		if np.any(new_array) == False:
			csumng = csum[seq_array.size-1]
  #          print("this is the csumng", csumng)
			lg_level = int(np.where(csum == csumng)[0])
			ng_level = int(sorted_lens[lg_level])
			stats['NG' + str(level)] = ng_level
		elif np.any(new_array) == True:
			csumng = min(csum[csum >= ngx])
   #         print("this is the csumng", csumng)
			lg_level = int(np.where(csum == csumng)[0])
			ng_level = int(sorted_lens[lg_level])
			stats['NG' + str(level)] = ng_level
	return stats


if __name__ == "__main__":
#    print(size_bp)
#    print(type(size_bp))
	naming = sys.argv[3]
	infilename = sys.argv[1]
	contig_lens, scaffold_lens, gc_cont = read_genome(infilename)
#    contig_stats = calculate_stats(contig_lens, gc_cont)
	scaffold_stats = calculate_stats(scaffold_lens, gc_cont)
	rounded_gc = round(gc_cont, 2)
#	print("this is the rounded_gc: ", rounded_gc)
	contig_stats = calculate_stats(contig_lens, gc_cont)
	gaps=contig_stats.get('sequence_count') - scaffold_stats.get('sequence_count')
	scaffold_stats['number_of_gaps'] = gaps
	contig_stats['number_of_gaps'] = gaps
#	print("this is the scaffold_stats: ", scaffold_stats)
#    print(scaffold_stats)
#    df_scaffold_all= pd.DataFrame.from_dict(scaffold_stats, orient= 'index')
#    print(df_scaffold_all)
#    df2 = pd.DataFrame([scaffold_stats]).T
#    print(df2)
#    s = pd.Series(scaffold_stats, name=naming)
#    s.index.name = 'Assembly:'
#    s.reset_index()
#    print(s)
	scaff_seq=pd.DataFrame(scaffold_stats.items(), columns=['Assembly:', naming])
#	print("this is the scaff_seq: ", scaff_seq)
	df_scaffold_top=scaff_seq.iloc[0:7,0:2]
#	print("this is GC CONTENT", gc_cont)
	df_scaffold_top[naming]=df_scaffold_top[naming].astype('int64').apply('{:,}'.format)
#	print("this is the top: ", df_scaffold_top)
	df_scaffold_top[naming][2] = rounded_gc
#	print("this is df_scaffold_top[naming][2]", df_scaffold_top[naming][2])
#	print("this is the top: ", df_scaffold_top)
    # df_scaffold_top=df_scaffold_top.style.hide_index()
#    df_scaffold_top[naming].round(decimals=0)
    # df_contig_all=pd.DataFrame(data=contig_stats)
	# df_contig_top=df_contig_all.iloc[0:6,0:2]
	df_scaffold_Nxx=pd.DataFrame(scaffold_stats.items(), columns=['Nxx Level (%)', 'Length of Nxx Scaffold (bp)'])
	df_scaffold_Nxx=df_scaffold_Nxx.iloc[31:55,0:2]
	df_scaffold_Nxx=df_scaffold_Nxx.reset_index()
#    print(df_scaffold_Nxx)
	df_scaffold_NGxx=pd.DataFrame(scaffold_stats.items(), columns=['NGxx Level (%)', 'Length of NGxx Scaffold (bp)'])
	df_scaffold_NGxx=df_scaffold_NGxx.iloc[79:104,0:2]
	df_scaffold_NGxx=df_scaffold_NGxx.reset_index()
#    print(df_scaffold_NGxx)
	df_scaffold_N_NG=pd.concat([df_scaffold_Nxx,df_scaffold_NGxx], axis=1)
	df_scaffold_N_NG=df_scaffold_N_NG.drop(df_scaffold_N_NG.columns[[0,3]], axis = 1)
	df_scaffold_N_NG['Length of Nxx Scaffold (bp)']=df_scaffold_N_NG['Length of Nxx Scaffold (bp)'].astype('int64').apply('{:,}'.format)
	df_scaffold_N_NG['Length of NGxx Scaffold (bp)']=df_scaffold_N_NG['Length of NGxx Scaffold (bp)'].astype('int64').apply('{:,}'.format)
    # df_scaffold_N_NG=df_scaffold_N_NG.style.hide_index()
	df_scaffold_Lxx=pd.DataFrame(scaffold_stats.items(), columns=['Lxx Level (%)', 'Count of scaffolds (for Lxx Level)'])
	df_scaffold_Lxx=df_scaffold_Lxx.iloc[7:31,0:2]
	df_scaffold_Lxx=df_scaffold_Lxx.reset_index()
#    print(df_scaffold_Nxx)
	df_scaffold_LGxx=pd.DataFrame(scaffold_stats.items(), columns=['LGxx Level (%)', 'Count of scaffolds (for LGxx Level)'])
	df_scaffold_LGxx=df_scaffold_LGxx.iloc[55:79,0:2]
	df_scaffold_LGxx=df_scaffold_LGxx.reset_index()
#    print(df_scaffold_NGxx)
	df_scaffold_L_LG=pd.concat([df_scaffold_Lxx,df_scaffold_LGxx], axis=1)
	df_scaffold_L_LG=df_scaffold_L_LG.drop(df_scaffold_L_LG.columns[[0,3]], axis = 1)
	df_scaffold_L_LG['Count of scaffolds (for Lxx Level)']=df_scaffold_L_LG['Count of scaffolds (for Lxx Level)'].astype('int64').apply('{:,}'.format)
	df_scaffold_L_LG['Count of scaffolds (for LGxx Level)']=df_scaffold_L_LG['Count of scaffolds (for LGxx Level)'].astype('int64').apply('{:,}'.format)
################################################################################################################

	contig_seq=pd.DataFrame(contig_stats.items(), columns=['Assembly:', naming])
	df_contig_top=contig_seq.iloc[0:7,0:2]
	df_contig_top[naming]=df_contig_top[naming].astype('int64').apply('{:,}'.format)
	df_contig_top[naming][2] = rounded_gc
    # df_scaffold_top=df_scaffold_top.style.hide_index()
#    df_scaffold_top[naming].round(decimals=0)
    # df_contig_all=pd.DataFrame(data=contig_stats)
	# df_contig_top=df_contig_all.iloc[0:6,0:2]
	df_contig_Nxx=pd.DataFrame(contig_stats.items(), columns=['Nxx Level (%)', 'Length of Nxx contig (bp)'])
	df_contig_Nxx=df_contig_Nxx.iloc[31:55,0:2]
	df_contig_Nxx=df_contig_Nxx.reset_index()
#    print(df_scaffold_Nxx)
	df_contig_NGxx=pd.DataFrame(contig_stats.items(), columns=['NGxx Level (%)', 'Length of NGxx contig (bp)'])
	df_contig_NGxx=df_contig_NGxx.iloc[79:104,0:2]
	df_contig_NGxx=df_contig_NGxx.reset_index()
#   print(df_scaffold_NGxx)
	df_contig_N_NG=pd.concat([df_contig_Nxx,df_contig_NGxx], axis=1)
	df_contig_N_NG=df_contig_N_NG.drop(df_contig_N_NG.columns[[0,3]], axis = 1)
	df_contig_N_NG['Length of Nxx contig (bp)']=df_contig_N_NG['Length of Nxx contig (bp)'].astype('int64').apply('{:,}'.format)
	df_contig_N_NG['Length of NGxx contig (bp)']=df_contig_N_NG['Length of NGxx contig (bp)'].astype('int64').apply('{:,}'.format)
	# df_scaffold_N_NG=df_scaffold_N_NG.style.hide_index()
	df_contig_Lxx=pd.DataFrame(contig_stats.items(), columns=['Lxx Level (%)', 'Count of contig (for Lxx Level)'])
	df_contig_Lxx=df_contig_Lxx.iloc[7:31,0:2]
	df_contig_Lxx=df_contig_Lxx.reset_index()
#    print(df_scaffold_Nxx)
	df_contig_LGxx=pd.DataFrame(contig_stats.items(), columns=['LGxx Level (%)', 'Count of contig (for LGxx Level)'])
	df_contig_LGxx=df_contig_LGxx.iloc[55:79,0:2]
	df_contig_LGxx=df_contig_LGxx.reset_index()
#    print(df_scaffold_NGxx)
	df_contig_L_LG=pd.concat([df_contig_Lxx,df_contig_LGxx], axis=1)
	df_contig_L_LG=df_contig_L_LG.drop(df_contig_L_LG.columns[[0,3]], axis = 1)
	df_contig_L_LG['Count of contig (for Lxx Level)']=df_contig_L_LG['Count of contig (for Lxx Level)'].astype('int64').apply('{:,}'.format)
	df_contig_L_LG['Count of contig (for LGxx Level)']=df_contig_L_LG['Count of contig (for LGxx Level)'].astype('int64').apply('{:,}'.format)
    # df_scaffold_L_LG=df_scaffold_L_LG.style.hide_index()
 #   print(df_contig_top)
#    print(scaffold_stats)
#    stat_output = {'Contig Stats': contig_stats,
#                   'Scaffold Stats': scaffold_stats}
#    print(json.dumps(stat_output, indent=2, sort_keys=False))
#    scaffold_out = {scaffold_stats}
#    contig_out = {contig_stats}


#    print(json.dumps(scaffold_stats, indent=2, sort_keys=False))
#    print(json.dumps(contig_stats, indent=2, sort_keys=False))
 #   dict_writer = csv.DictWriter(stat_output, fieldnames=None, delimiter='\t')
 #   print(dict_writer)
#    df_raw = pd.DataFrame.from_dict(data, orient='index')
 #   print(df_raw)
#    with open('statdata.txt', 'w') as outfile:
#        json.dump(scaffold_stats, outfile)
#    with open('contigdata.txt', 'w') as out2file:
#        json.dump(contig_stats, out2file)
    # scaff = csv.writer(open(naming + "_scaffold_stats.tsv", "w"), delimiter='\t')
    # for key, val in scaffold_stats.items():
    #     scaff.writerow([key, val])
	#
    # contig = csv.writer(open(naming + "_contig_stats.tsv", "w"), delimiter='\t')
    # for key, val in contig_stats.items():
    #
	with open(sys.argv[5], 'w') as outputfile:
# #    print('#' + libraryName, file=outputfile)
# #    print("Total Reads Processed (Paired):        " + total_processed + "   ( 100 %)", file=outputfile)
# #    print("Discarded reads (Paired):              " + discarded + "    ( "+str(discarded_perc)+"%)", file=outputfile)
# #    print("Successfully Processed reads (Paired): " + successfull + "   ( "+str(successfull_perc)+"%)", file=outputfile)
		print(df_scaffold_top.to_string(index=False), file=outputfile)
		print("", file=outputfile)
		print(df_scaffold_N_NG.to_string(index=False), file=outputfile)
		print("", file=outputfile)
		print(df_scaffold_L_LG.to_string(index=False), file=outputfile)

	with open(sys.argv[6], 'w') as outputfile2:
# #    print('#' + libraryName, file=outputfile)
# #    print("Total Reads Processed (Paired):        " + total_processed + "   ( 100 %)", file=outputfile)
# #    print("Discarded reads (Paired):              " + discarded + "    ( "+str(discarded_perc)+"%)", file=outputfile)
# #    print("Successfully Processed reads (Paired): " + successfull + "   ( "+str(successfull_perc)+"%)", file=outputfile)
		print(df_contig_top.to_string(index=False), file=outputfile2)
		print("", file=outputfile2)
		print(df_contig_N_NG.to_string(index=False), file=outputfile2)
		print("", file=outputfile2)
		print(df_contig_L_LG.to_string(index=False), file=outputfile2)

    #  with open(sys.argv[4], 'w') as outputRst:
    #      print(tabulate(df_scaffold_top, headers='keys',tablefmt="rst", showindex=False), file=outputRst)
    #      print("", file=outputRst)
    #      print(tabulate(df_scaffold_N_NG, headers='keys',tablefmt="rst", showindex=False), file=outputRst)
    #      print("", file=outputRst)
    #      print(tabulate(df_scaffold_L_LG, headers='keys',tablefmt="rst", showindex=False), file=outputRst)
    #      print("", file=outputRst)
	# #
    #  with open(sys.argv[5], 'w') as outputRst2:
    #      print(tabulate(df_contig_top, headers='keys',tablefmt="rst", showindex=False), file=outputRst2)
    #      print("", file=outputRst2)
    #      print(tabulate(df_contig_N_NG, headers='keys',tablefmt="rst", showindex=False), file=outputRst2)
    #      print("", file=outputRst2)
    #      print(tabulate(df_contig_L_LG, headers='keys',tablefmt="rst", showindex=False), file=outputRst2)
    #      print("", file=outputRst2)

	# with open(sys.argv[4], 'w') as outputRst:
	# 	print(tabulate(df_scaffold_top, headers='keys',tablefmt="pipe", showindex=False), file=outputRst)
	# 	print("", file=outputRst)
	# 	print(tabulate(df_scaffold_N_NG, headers='keys',tablefmt="pipe", showindex=False), file=outputRst)
	# 	print("", file=outputRst)
	# 	print(tabulate(df_scaffold_L_LG, headers='keys',tablefmt="pipe", showindex=False), file=outputRst)
	# 	print("", file=outputRst)
	# #
	# with open(sys.argv[5], 'w') as outputRst2:
	# 	print(tabulate(df_contig_top, headers='keys',tablefmt="pipe", showindex=False), file=outputRst2)
	# 	print("", file=outputRst2)
	# 	print(tabulate(df_contig_N_NG, headers='keys',tablefmt="pipe", showindex=False), file=outputRst2)
	# 	print("", file=outputRst2)
	# 	print(tabulate(df_contig_L_LG, headers='keys',tablefmt="pipe", showindex=False), file=outputRst2)
	# 	print("", file=outputRst2)
	# list_of_dfs=[df_scaffold_top,df_scaffold_N_NG,df_scaffold_L_LG]
    # for df in list_of_dfs:
    #     with open('all_dfs.tsv','a') as f:
    #         df.to_csv(f)
    #         f.write("\n")
#     with open("testok.tsv", 'w') as outputfile:
# #    print('#' + libraryName, file=outputfile)
# #    print("Total Reads Processed (Paired):        " + total_processed + "   ( 100 %)", file=outputfile)
# #    print("Discarded reads (Paired):              " + discarded + "    ( "+str(discarded_perc)+"%)", file=outputfile)
# #    print("Successfully Processed reads (Paired): " + successfull + "   ( "+str(successfull_perc)+"%)", file=outputfile)
#         print(tabulate(loadinTop, headers='keys',tablefmt="rst", showindex=False), file=outputfile)
# #    print('', file=outputfile)
#         print(tabulate(result, headers='keys', tablefmt="psql", showindex=False), file=outputfile)
#         print('', file=outputfile)
