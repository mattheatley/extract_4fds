#!/usr/bin/python
# -*- coding: utf-8 -*-
# CONDA ENV SOFTWARE: python (3.7), samtools (v1.9), gatk4 (v4.1.4.0)

import os, sys, subprocess, shutil, re, argparse, math
from core import CAPTURE, SBATCH, ezSub, FindSupplementaryFile

pipe_script = os.path.realpath(__file__) # extract file path
find_script = f'{os.path.dirname(pipe_script)}/find_4fds.py'
script_name, *arguments = sys.argv # extract command line arguments

parser = argparse.ArgumentParser(description='Extract 4-fold degenerate sites', prog=script_name, usage='%(prog)s [options]', epilog='see the readme file for further details.')

site_types = ['all','consistent']
modes = parser.add_mutually_exclusive_group(required=True) # run modes
modes.add_argument('-setup', action='store_true', help='setup initial directories')
modes.add_argument('-index', action='store_true', help='index reference genome')
modes.add_argument('-find', action='store_true', help='find 4-fold degenerate sites')
modes.add_argument('-extract', metavar='<category>', choices=site_types, help='extract sites from vcfs')
modes.add_argument('-merge', metavar='<category>', choices=site_types, help='merge extracted site vcfs')

additional = parser.add_argument_group(description='additional modes:') # check modes
additional.add_argument('-pipe', action='store_true', help='submit pipe as task')
additional.add_argument('-test', action='store_true', help='test script locally')

inputs = parser.add_argument_group(description='user inputs:') # user inputs
inputs.add_argument('-u', metavar='<name>', type=str, help='specify user name')
inputs.add_argument('-p', metavar='</path/to/>', type=str, help='specify path')
inputs.add_argument('-d', metavar='<name>', type=str, default='pipe_4fds', help='specify working directory')
inputs.add_argument('-m', metavar='<name>', type=str, help='specify user email address')
inputs.add_argument('-e', metavar='<name>', type=str, default='ngs_env', help='specify conda environment')
inputs.add_argument('-l', metavar='<number>', type=int, default=100, help='specify parallel task limit')

hpcc = parser.add_argument_group(description='hpcc settings:') # overide hpcc settings
hpcc.add_argument('-pa', metavar='<name>', type=str, choices=['defq','shortq','devq','mmemq','hmemq','voltaq','visq'], help='specify partition')
hpcc.add_argument('-no', metavar='<number>', type=int, help='specify nodes')
hpcc.add_argument('-nt', metavar='<number>', type=int, help='specify ntasks')
hpcc.add_argument('-me', metavar='<number[units]>', type=str.lower, help='specify memory [k|m|g|t]')
hpcc.add_argument('-wt', metavar='<HH:MM:SS>', type=str, help='specify walltime')

setup, indexing, finding, extracting, merging, submitting_self, testing, user, wrk_path, wrk_dir, address, environment, limit, *hpcc_overides = vars(parser.parse_args()).values() # define user inputs

*_, memory, walltime = hpcc_overides # define hpcc overides
if memory and (sum(memory.count(unit) for unit in ['k','m','g','t']) != 1 or not memory.rstrip('kmgt').isdigit()): parser.error(f"argument -me: invalid str format: '{memory}'") # check memory format
if walltime and ( not 8 <= len(walltime) <= 9 or walltime[-6:-2:3].count(':') != 2 or not walltime.replace(':','').isdigit()): parser.error(f"argument -wt: invalid str format: '{walltime}'") # check walltime format

(pipe_flag, *_), *_ = [ info.option_strings for info in additional._group_actions ] # extract pipe flag

mode_flags = [ flag for sublist in [ mode.option_strings for mode in modes._group_actions ] for flag in sublist ] # extract mode flags

active_mode, = set(arguments).intersection(set(mode_flags)) # determine active mode
mode_description =  active_mode.strip("-") if setup else f'{active_mode.strip("-").rstrip("e")}ing' # specify active mode description
relevant_sites = extracting if extracting else merging if merging else None

resource_labels = 'partition','nodes','ntasks-per-node', 'memory', 'walltime' # specify subdictionary labels
hpcc_settings = { # default hpcc settings 
    'indexing': ['devq','1','1','4g', '01:00:00'], 
    'finding': ['devq','1','1','4g', '01:00:00'], 
    'extracting': ['shortq','1','1','4g' ,'12:00:00'],
    'merging': ['shortq','1','1','4g' ,'12:00:00'],
    '+': ['defq','1','1','4g' ,'168:00:00'] }
hpcc_settings = { key:  {label:info for label,info in zip(resource_labels,value)} for key,value in hpcc_settings.items() } # creatre subdictionaries

print(f'\n{mode_description.upper()}\n')

user, wrk_path = [ os.getenv(bash) if not argument else argument for argument,bash in [ (user,'USER'),(wrk_path,'HOME')] ] # find user & path
wrk_dir = f'/{wrk_path.strip("/")}/{wrk_dir.strip("/")}' # ensure correct path format

dir_labels = ['ref.gen','gff','vcf', 'coordinates', f'extracted_{relevant_sites}', 'merged']
batch_labels = [ f'x.slurm/{"pipe" if submitting_self  else mode_description}/{label}' for label in ['scripts','out.err','ids'] ]

dir_paths = [ f'{wrk_dir}/{label}' for label in [*dir_labels, *batch_labels] ] # specify directory paths
ref_dir, gff_dir, vcf_dir, coord_dir, extract_dir, merge_dir, sh_dir, oe_dir, id_dir = dir_paths # specify initial & stage directories
split_subdir = f'{coord_dir}/split_{relevant_sites}' # split directory
tmp = f'{split_subdir}/tmp' # temporary split directory

initial_dirs = [ref_dir, gff_dir, vcf_dir] # specify initial directories

if setup: 
    [ os.makedirs(path, exist_ok=True) for path in initial_dirs ] # make initial directories as required
    print('SETUP COMPLETE\n', flush=True)
    
else: # proceed with alternative pipeline mode 
    assert(os.path.exists(wrk_dir)), f'Problem finding the working directory "{wrk_dir}".'# check path exists
    assert(environment in CAPTURE('conda info --env')), f'Problem finding the conda environment "{environment}".' # check conda environment exists


    out_dir = coord_dir if finding else extract_dir if extracting else merge_dir
    stage_dirs = [out_dir, sh_dir, oe_dir, id_dir]
    if extracting: stage_dirs.append(tmp)
    [ os.makedirs(stage_dir, exist_ok=True) for stage_dir in stage_dirs ] # make stage-specific directories as required

    id_file = f'{id_dir}/task.ids' # specify slurm id file
    list_labels = ['all_sites.list','consistent_sites.list','ignored_transcripts']
    all_sites_output, consistent_sites_output, ignored_transcripts = [ f'{coord_dir}/{label}' for label in list_labels ] # specify initial directories

    hpcc_settings[mode_description].update({ label:setting for label,setting in zip(resource_labels, hpcc_overides) if setting }) # update hpcc settings as required
    partition, nodes, ntasks, memory, walltime = resources = [  hpcc_settings['+' if submitting_self else mode_description][label] for label in resource_labels ] # extract hpcc settings

    setting_info = zip(['user','working directory','conda environment',*resource_labels],[user,wrk_dir,environment,*resources])
    [ print(f'{label}: {info}') for label,info in setting_info ]

    print('\nFINDING SUPPLEMENTARY FILES...')
    list_suffix, vcf_suffix, extracted_infix, scaffold_prefix = '.list','vcf.gz', '4fds', 'scaffold_'
    if indexing or finding or extracting:
        fasta_file = FindSupplementaryFile(ref_dir, ('.fasta','.fa','.faa','.fas')) # find reference genome
        if indexing:
            if not fasta_file.endswith('.fasta'): # check suffix appropriate for GATK
                print('\tRENAMING REFERENCE GENOME')
                fasta_prefix, *_ = fasta_file.rsplit('.', 1) # extract reference prefix
                new_fasta_file = f'{fasta_prefix}.fasta' # specify new file name
                os.rename(fasta_file,new_fasta_file) # rename original file
                fasta_file = new_fasta_file # update reference genome variable
            chr_names = [ name.lstrip('>') for name in CAPTURE(f'grep \> {fasta_file}').split('\n') ] # extract chromosome names
            dodgy_names = any('|' in name for name in chr_names) # check if any problematic characters in chromosome names
            proceeding = False if dodgy_names else True
            while not proceeding: # allow user to manually change scaffold names
                response = input('\nSome scaffold names contain pipe characters that can cause problems with the software used here - would you like to quit & remove these? (y/n): ') # instruct task resubmission
                proceeding = response in ['y','n']
                if response == 'n': print('\nOk but you\'ve been warned.\n')
                if response == 'y': print('\nExiting.\n'); sys.exit(0)
    
    if finding:
        gff_file = FindSupplementaryFile(gff_dir, ('.gff','.gff3'))

    if extracting or merging: 
        vcf_file = FindSupplementaryFile(vcf_dir, (vcf_suffix))

    if extracting: 
        sites_file = all_sites_output if extracting == 'all' else consistent_sites_output 
        sites_file = FindSupplementaryFile(coord_dir, os.path.basename(sites_file))

    if merging:
        fai_file = FindSupplementaryFile(ref_dir, '.fai')


    print('\tSUPPLEMENTARY FILES FOUND!')

    if extracting or merging:
        
        vcf_base_name = os.path.basename(vcf_file).replace(vcf_suffix,'').strip('.')

        if extracting:

            print('\nSPLITTING SITE FILE...') # LONG WINDED APPROACH THAT IS FASTER THEN SPLITTING SINGLE FILE BY SCAFFOLD
            
            part_prefix = 'part-' # specify prefix for site list chunk            
            part_size = 100000 # specify lines per site list chunk
            total_sites,  *_ = CAPTURE(f'wc -l {sites_file}').split(' ') # calculate total 4fds site
            total_parts = math.ceil(int(total_sites)/part_size) # calculate total chunks

            subprocess.call(f'split -l {part_size} --numeric-suffixes=1 --additional-suffix={list_suffix} {sites_file} {tmp}/{part_prefix}', shell=True) # split site list into chunks

            split_files = sorted([ contents.path for contents in os.scandir(tmp) if contents.name.startswith(part_prefix) ]) # find site list chunks

            scaffolds_with_4fds = { scaffold:[] for scaffold in CAPTURE(f'cut -d ":" -f 1 {sites_file} | uniq').split('\n') } # create temporary store

            for part in split_files: [ scaffolds_with_4fds[scaffold].append(part) for scaffold in CAPTURE(f'cut -d ":" -f 1 {part} | uniq').split('\n') ] # assign site list chunks by scaffold

            for scaffold,relevant_split_files in scaffolds_with_4fds.items():
                if relevant_split_files: 
                    scaffold_sites_output = f'{split_subdir}/{scaffold_prefix}{scaffold}.list' # specify scaffold sites list

                    subprocess.call(f'> {scaffold_sites_output}', shell=True) # create scaffold sites list
                    [ subprocess.call(f'''awk -F: '$1=="{scaffold}"' {part} >> {scaffold_sites_output}''', shell=True) for part in relevant_split_files ] # extract sites & assign to scaffold site list

            shutil.rmtree(tmp)
            split_files = sorted([ contents.path for contents in os.scandir(split_subdir) if contents.name.endswith(list_suffix) ]) # find site list chunks
            
            print('\tSPLITTING COMPLETE!')

    scripts = []

    index_cmds = ['samtools faidx','gatk CreateSequenceDictionary -R'] # define index stages

    to_process = index_cmds if indexing else [1] if finding or merging else split_files # specify tasks to process

    for i,task in enumerate(to_process, 1):
        
        if indexing: tool, *_ = task.split(' ', 1) # extract software name from command
        
        task_id = mode_description if not indexing else tool # specify task id
        
        if extracting:
            part = os.path.basename(task).replace(list_suffix, '').strip('.') # extract scaffold id
            task_id += f'-{relevant_sites}_{part}' # include scaffold id

        if extracting or merging:
            output_vcf = f'{out_dir}/{vcf_base_name}.{extracted_infix}.{relevant_sites}.{{}}.{vcf_suffix}'.format(part if extracting else 'merged')

        sh_file = f'{sh_dir}/{task_id}.sh' # specify sh file

        scripts.append(sh_file) # record sh file

        with open(sh_file, 'w') as sh:
                    
            hpcc_directives = SBATCH(job_id=task_id, partition=partition, nodes=nodes, ntasks=ntasks, memory=memory, walltime=walltime, out_err=oe_dir, conda_env=environment) # specify hpcc directives (slurm)
            sh.write(hpcc_directives)                

            sh.write('echo TASK STARTED `date`\n') # log start time

            # INDEX
            
            if indexing: sh.write(f'{task} {fasta_file}\n') # index FASTA reference sequence

            # FIND

            if finding: sh.write(f'python {find_script} -g {gff_file} -f {fasta_file} -o {out_dir} -l {" ".join(list_labels)}\n') # run find script

            # EXRACT

            if extracting:

                sh.write('gatk SelectVariants \\\n'
                +f'-R {fasta_file} \\\n'
                +f'-V {vcf_file} \\\n'
                +f'-L {task} \\\n'
                +f'-O {output_vcf} \n')

            # MERGE

            if merging:

                all_scaffolds = CAPTURE(f'cut -f 1 {fai_file} | uniq').split('\n') # extract indexed scaffolds
                
                raw_prefix, raw_suffix = [re.escape(affix) for affix in [scaffold_prefix,vcf_suffix] ] # specify regex search pattern

                input_vcfs = { re.search(f'{raw_prefix}(.*)\.{raw_suffix}', contents.name).group(1): contents.path for contents in os.scandir(extract_dir) if contents.name.endswith(vcf_suffix) }

                scaffolds_with_4fds = sorted(input_vcfs, key=lambda scaffold: all_scaffolds.index(scaffold)) # sort scaffolds with 4fds by indexed position

                sh.write('gatk GatherVcfs \\\n')
                for scaffold in scaffolds_with_4fds:
                    sh.write(f'-I {input_vcfs[scaffold]} \\\n')
                sh.write(f'-O {output_vcf} \n'
                +'gatk IndexFeatureFile \\\n'
                +f'-F {output_vcf} \n') 

            sh.write('echo TASK COMPLETED `date` \n') # log end time 


    if testing:
        print('\nIndividual Task Scripts:\n')
        for sh_file in scripts: print(open(sh_file, 'r').read()) 
        print('\nTESTING COMPLETE\n')
    else: # submit tasks to hpcc

        if submitting_self: # submit pipeline to trickle tasks from hpcc
            pipe_id = f'pipe-{mode_description}'
            sh_script = f'{sh_dir}/{pipe_id}.sh'
            with open(sh_script, 'w') as sh:  
                hpcc_directives = SBATCH(job_id=pipe_id, partition=partition, nodes=nodes, ntasks=ntasks, memory=memory, walltime=walltime, out_err=oe_dir, conda_env=environment) # specify hpcc directives (slurm)
                sh.write(hpcc_directives)                
                arguments.remove(pipe_flag) # remove self submission flag
                print('python', pipe_script, *arguments, sep=' ', file=sh) # submit self as task
            pipe_id = CAPTURE(f'sbatch {sh_script}')
            print(f'\nPIPELINE SUBMITTED ({pipe_id})\n')
        
        else: # submit tasks
            print('\nSUBMITTING TASKS...')                    
            with open(id_file, 'a+') as id_output:
                for i,script in enumerate(scripts,1): 
                    ezSub(i=i, check=600, user=user, limit=limit) # maintain tasks below parellel task limit
                    sub_id = CAPTURE(f'sbatch -d singleton {script}') # submit task
                    print(sub_id, script, sep='\t', file=id_output, flush=True) # record task job id & shell script
            print('\tALL TASKS SUBMITTED\n')
            


