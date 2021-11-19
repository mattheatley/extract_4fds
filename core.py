import os, sys, subprocess, time

def CAPTURE(cmd): return subprocess.run(f'{cmd}', shell=True, capture_output=True).stdout.decode('utf-8').strip(' \n') # capture & format terminal output


def SBATCH(job_id, partition, nodes, ntasks, memory, walltime, out_err, task=None, email=None, conda_env=None):
    labels = ['job-name','partition','nodes','ntasks-per-node','mem','time','output','error', 'parsable'] # specify relevant directives
    inputs = [job_id, partition, nodes, ntasks, memory, walltime, *[ f'{out_err}/%x{f".{task}" if task else ""}.{suffix}' for suffix in ['out','err'] ], None] # organise settings
    if email: labels.extend(['mail-user','mail-type']); inputs.extend([email,'END']) # add optional user email address
    sbatch = ''.join([ f'#SBATCH --{option}{f"={info}" if info else ""} \n' for option,info in zip(labels,inputs) ]) # format directives & settings
    sbatch += '\nsource $HOME/.bash_profile\n'
    if conda_env: sbatch += f'\nconda activate {conda_env}'+'\necho ENVIRONMENT $CONDA_DEFAULT_ENV ACTIVE\n' # add optional conda environment
    return '#!/bin/bash\n'+f'{sbatch}\n'+'echo RUNNING ON `hostname`\n'


def ezSub(i, check, user, limit): # auto submission
    total = int(CAPTURE(f'squeue -u {user} -h | wc -l')) # find current tasks
    if total > limit: 
        print(f'\nSUBMISSION LIMIT REACHED: WAITING TO SUBMIT TASK {i}...')
        time.sleep(check)
        ezSub(i, check, user, limit) # check every 5 mins if task can be submitted
    else: 
        return # submit task (with dependancy) if possible


def REVIEW(id_file):
    
    print('\nREVIEWING TASKS...')
    with open(id_file, 'a+') as id_update:
        
        *non_fails, failed = states = ['PENDING', 'RUNNING', 'COMPLETED', 'CANCELLED', 'FAILED'] # specify slurm job state categories
        categories = {category:set() for category in states } # establish slurm job state categories
        
        id_update.seek(0) # reset file position to read from beginning
        sub_info = { sub_id:script for sub_id,script, *_ in [line.strip('\n').split('\t') for line in id_update.readlines()] } # extract task job ids & scripts info
        
        sub_ids = ",".join(sub_info.keys()) # specify task job id list to review
        headers, *sacct_info = [ line.split('|') for line in CAPTURE(f'sacct -p -j {sub_ids}').split('\n') ] # extract slurm accounting data for tasks
        
        for info in sacct_info:
            sub_id, step, state = [ info[headers.index(column)].strip('+') for column in ['JobID','JobName','State'] ] # extact specific slurm accounting data
            if not sub_id.endswith(('batch','extern')): categories[failed if not state in non_fails else state].add(sub_info[sub_id]) # categorise tasks by state
        pending, running, completed, cancelled, failed = categories.values() # extract categorised tasks
        problems = failed.difference( set().union(pending, running, completed) ) # extract failed tasks that are not running (i.e. re-submitted)
        if pending: print(f'\n - TASKS PENDING: {len(pending)}')
        if running: print(f'\n - TASKS RUNNING: {len(running)}')
        if (pending or running) and not problems: print('\nNo problems identified for current tasks.')
        if completed and not problems and not pending and not running: total = len(completed); print(f'\nAll {total} tasks have completed.')
        if problems: 
            print('\nIt seem\'s that there were problems with the following tasks that havn\'t been dealt with yet:\n')
            [ print(os.path.basename(task)) for task in sorted(problems) ]
        proceed = False
        if problems: 
            while proceed is False:
                response = input('\nwould you like to repeat these tasks now? (y/n): ') # instruct task resubmission
                proceed, repeat = response in ['y','n'], response == 'y' # interpret instructions
                if repeat is True:
                    for script_file,*_ in sorted(problems): 
                        resub_id = CAPTURE(f'sbatch {script_file}') # resubmit task
                        print(resub_id, script_file, '(RESUBMITTED)', sep='\t', file=id_update) # record resubmitted id
        print('\nREVIEW COMPLETE\n')


def FindSupplementaryFile(search_dir, search_suffix): # find required files
    assert(os.path.isdir(search_dir)), f'Problem finding the directory \"{search_dir}\" when searching for supplementary files.'
    found_files = [ found.path for found in os.scandir(search_dir) if found.name.endswith(search_suffix) ]
    assert(len(found_files) == 1), f'{"Multiple" if found_files else "No"} files found in \"{search_dir}\" ending with {JoinNeatly(search_suffix)}.'
    found_file, *_ = found_files
    return found_file


def JoinNeatly(to_join):
    *leading, trailing = [f'\"{info}\"' for info in to_join]
    if leading: leading, *_ = [(", ").join(leading)] if len(leading) > 1 else leading # join by commas if relevant
    return f'either {leading} or {trailing}' if leading else trailing # join by & if relevant
