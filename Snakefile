import os,time,shutil
configfile: 'configs/eg.yaml'

# argument
species = config['species']

# structure
indir=os.path.abspath(config['indir'])
outdir=os.path.abspath(config['outdir'])
Lsample=config['Lsample']
Plog=f'{outdir}/0log'
Pigblast=f'{outdir}/1igblast'

for p in [Pigblast]:
    os.makedirs(p, exist_ok=True)

onstart:
    for r in workflow.rules:
        os.makedirs(f'{Plog}/{r}', exist_ok=True)
    global time_start, smk_file_name
    time_start = time.strftime("%y%m%d_%H%M%S", time.localtime())
    smk_file_name = workflow.snakefile.split('/')[-1]
    shutil.copyfile(workflow.snakefile, f'{Plog}/all/{time_start}_{smk_file_name}')


rule all:
    input:
        igblast_tsv = expand(Pigblast + '/{sample}/igblast_airr.tsv', sample=Lsample),
    
    
rule igblast:
    input: fa = indir + '/{sample}.fa',
    output:
        igblast_tsv = Pigblast + '/{sample}/igblast_airr.tsv',
        igblast_txt = Pigblast + '/{sample}/igblast_blast.txt',
    log: e = Plog + '/igblast/{sample}.e', o = Plog + '/igblast/{sample}.o'
    benchmark: Plog + '/igblast/{sample}.bmk'
    resources: cpus=config['igblast_cpus']
    params: 
        igblast_ref_prefix = f'ref/igblast-ref/germline_db/{species}-VDJB/{species}_',
        igblast_igdata = "ref/igblast-ref",
        igblast_aux = f'ref/igblast-ref/optional_file/{species}_gl.aux'
    conda: 'envs/env.yaml'
    shell:"""
        export IGDATA={params.igblast_igdata}

        # airr format
        igblastn -organism {species} -domain_system imgt \\
            -auxiliary_data {params.igblast_aux} \\
            -germline_db_V {params.igblast_ref_prefix}V \\
            -germline_db_D {params.igblast_ref_prefix}D -germline_db_J {params.igblast_ref_prefix}J \\
            -show_translation -outfmt 19 -num_threads {resources.cpus} \\
            -query {input.fa} \\
            -out {output.igblast_tsv} 1>>{log.o} 2>>{log.e}
        
        # blast format
        igblastn -organism {species} -domain_system imgt \\
            -auxiliary_data {params.igblast_aux} \\
            -germline_db_V {params.igblast_ref_prefix}V \\
            -germline_db_D {params.igblast_ref_prefix}D -germline_db_J {params.igblast_ref_prefix}J \\
            -show_translation -outfmt '7 std qseq sseq btop' -num_threads {resources.cpus} \\
            -query {input.fa} \\
            -out {output.igblast_txt} 1>>{log.o} 2>>{log.e}
        """