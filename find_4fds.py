import os, sys, subprocess, shutil, re, argparse, time
from core import CAPTURE

parser = argparse.ArgumentParser(description='Find 4-Fold Degenerate Sites', prog=os.path.basename(__file__), usage='%(prog)s [options]', epilog='to find 4-fold degenerate sites in a fasta file.')
parser.add_argument('-g', metavar='</path/to/annotation.gff>', type=str, required=True, help='specify path to gff annotation')
parser.add_argument('-f', metavar='</path/to/reference.fasta>', type=str, required=True, help='specify path to reference genome')
parser.add_argument('-o', metavar='</path/to/output_directory>', type=str, required=True, help='specify path to output directory')
parser.add_argument('-l', metavar='<file_name>', nargs=3, required=True, help='specify ouput file names')

*paths, labels = vars(parser.parse_args()).values() # define user inputs

gff_file, fasta_file, out_dir = [ f'/{path.strip("/")}' for path in paths ] # ensure correct path format
assert(os.path.exists(gff_file)), f'Problem finding gff file.'# check path exists
assert(os.path.exists(fasta_file)), f'Problem finding reference genome file.'# check path exists
os.makedirs(out_dir, exist_ok=True) # create output directory as required
all_sites_file, consistent_sites_file, ignored_transcripts_file = [ f'{out_dir}/{label}' for label in labels] # specify output files

genes, mRNA = {},{} # establish region, gene, mRNA ditionaries
sequence, starts, ends, strand, transcripts = 'sequence', 'starts', 'ends', 'strand', 'transcripts' # define sub-dictionary categories
unassigned_mRNA_file, unassigned_CDS_file = [ f'{out_dir}/unassigned_{label}' for label in ['mRNA','CDS'] ]

degenerate_codons = { # 4-fold dengerate sites for forward strands
'CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG','CCT','CCC','CCA','CCG','CTT','CTC','CTA','CTG','TCT','TCC',
'TCA','TCG','ACT','ACC','ACA','ACG','GGT','GGC','GGA','GGG','GCT','GCC','GCA','GCG','GTT','GTC','GTA','GTG'}

print(f'\nPROCESSING GFF FILE...', flush=True)

with open(unassigned_mRNA_file,'w') as unassigned_mRNA, open(unassigned_CDS_file,'w') as unassigned_CDS:

    unassigned_features = 0

    for line in open(gff_file, 'r').readlines(): # organise gff enteries
        
        if not line.strip() or line.startswith('#'): continue  # ignore blank lines & those that are comments
        else:

            scaffold, source, SO, start, end, score, sign, phase, attributes = line = line.strip().split('\t') # extract gff info

            if not scaffold in genes: 
                genes[scaffold] = {} # establish region sub-dictionary as required
            
            attributes = sorted([ re.split('=|,',info) for info in attributes.rstrip('\n').split(';') if info.startswith('ID=' if SO == 'gene' else ('ID=','Parent='))]) # extract ids/parents for genes, mRNA & cds
                    
            if SO == 'gene': 
                (*_,ID), = attributes # extract info
                genes[scaffold][ID] = { # create entry for gene id
                    starts: int(start), 
                    ends: int(end), 
                    strand: sign, 
                    transcripts: [] } 

            if SO in ['mRNA','CDS']: 

                (*_,ID), (*_,parent) = attributes # extract info
                
                if SO == 'mRNA': 
                    mRNA[ID] = [] # create entry for mRNA id
                    if not parent in genes[scaffold]: 
                        unassigned_features+=1
                        print(ID, file=unassigned_mRNA)
                    else: genes[scaffold][parent][transcripts].append(ID) # store mRNA id with relevant gene

                if SO == 'CDS': 
                    if not parent in mRNA: 
                        unassigned_features+=1
                        print(ID, file=unassigned_CDS)
                    else: mRNA[parent].append( (int(start), int(end)) ) # store CDS coorindates with relevant mRNA

    print(f'\nGFF FILE PROCESSED', flush=True)
    
    if unassigned_features: print('\nWARNING: Some features could not be assigned to a parent; see files for further details.')


print('\nFINDING SITES...', flush=True)

with open(all_sites_file,'w') as all_output, open(consistent_sites_file,'w') as consistent_output, open(ignored_transcripts_file,'w') as ignored_transcripts:

    total_4fds, total_consistent_4fds = 0, 0
    
    for scaffold, entry in sorted(genes.items()): # cycle scaffolds

        all_scaffold_4fds, consistent_scaffold_4fds, scaffold_genes = 0, 0, len(entry)

        for i, (gene_id, gene_info) in enumerate( sorted(entry.items()),1 ): # cycle scaffold genes

            if i % 1000 == 0: print(f'{i} OF {scaffold_genes} GENES PROCESSED FOR SCAFFOLD {scaffold.upper()}\n', flush=True) # log progress

            gene_4fds = [] # establish entries for gene sites
                        
            if not gene_info[transcripts]: print(f'\n Skipping gene {gene_id}; no transcripts.', flush=True) # ignore genes with no transcripts
            else: # process transcripts

                sequence = ''.join(CAPTURE(f'samtools faidx {fasta_file} {scaffold}:{gene_info[starts]}-{gene_info[ends]}').upper().split('\n')[1:]) # extract gene sequence from fasta
                positions = range(gene_info[starts], gene_info[ends]+1) # calculate base positions        
                base_positions = list( zip(sequence,positions) ) # merge sequence & base positions

                for mRNA_id in sorted( gene_info[transcripts] ): # cycle scaffold gene transcripts

                    spliced_sequences, transcript_4fds = [],set() # establish current transcript sets

                    CDS_coordinates = sorted(mRNA[mRNA_id], key=lambda coordinates: coordinates[0] ) # sort coordinates in ascending order according to start & re-extract cds coorinates 
                    
                    for CDS_start, CDS_end in CDS_coordinates: # cycle transcript cds
                        offset = gene_info[starts]
                        splice = base_positions[ CDS_start-offset:CDS_end-offset+1 ] # adjust coordinates & slice sequence
                        spliced_sequences.extend(splice) # store cds base positions
                    
                    if gene_info[strand] != '+': # determine reverse compliment if required
                        spliced_sequences = [ (base.translate(str.maketrans({'G':'C','C':'G','T':'A','A':'T'})), site) for base,site in spliced_sequences ] # convert to forward strand
                        spliced_sequences.reverse() # reverse sequence
        
                    if len(spliced_sequences) % 3 != 0: # check transcript divisible by 3 (i.e. unambiguously translated)
                        print(mRNA_id, file=ignored_transcripts, flush=True) # log & output amibguously translated transcripts
                        print(f'\n Skipping transcript {mRNA_id}; reading frame ambiguous.', flush=True)
                    else:# process cds

                        for i in range(0, len(spliced_sequences), 3): # cycle triplicate CDS bases    
                        
                            (base1,*_), (base2,*_), (base3,position) = spliced_sequences[i:i+3] # extract base positions
                        
                            codon = base1+base2+base3 # re-form codon
                        
                            if codon in degenerate_codons: # compare codon against degenerate codons
                        
                                transcript_4fds.add(f'{scaffold}:{position}-{position}') # record transcript 4fds coordinate 
                    
                    gene_4fds.append(transcript_4fds) # record gene 4fds coorindates

                all_gene_4fds = set.union(*gene_4fds)  # extract all 4fds sites found in gene transcripts
                all_scaffold_4fds += len(all_gene_4fds) # update scaffold site totals

                consistent_gene_4fds = set.intersection(*gene_4fds) # extract 4fds sites found in all gene transcripts
                consistent_scaffold_4fds += len(consistent_gene_4fds) # update scaffold site totals
                
                outputs = [
                    (all_gene_4fds, all_output)
                    ,(consistent_gene_4fds, consistent_output)
                           ]

                [ [print(site,file=site_output, flush=True) for site in sorted(site_list, key=lambda info: int(info.split('-')[-1]))] for site_list,site_output in outputs ] # output coordinates
                
        total_4fds += all_scaffold_4fds # update genome site totals
        total_consistent_4fds += consistent_scaffold_4fds # update genome site totals
    
        print(f'\nALL {scaffold_genes} GENES PROCESSED FOR SCAFFOLD {scaffold.upper()}: {all_scaffold_4fds} SITES / {consistent_scaffold_4fds} CONSISTENT SITES', flush=True)
    
    print(f'\nFINDING SITES COMPLETE: {total_4fds} SITES / {total_consistent_4fds} CONSISTENT SITES\n', flush=True) # log total genome sites
