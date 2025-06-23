import os,time,shutil
configfile: 'configs/eg.yaml'

# argument
species = config['species']

# structure
indir=os.path.abspath(config['indir'])
outdir=os.path.abspath(config['outdir'])
Lsample=config['Lsample']
Plog=f'{outdir}/0log'
Pparse=f'{outdir}/1parse'

for p in [Pparse]:
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
        igblast_tsv = expand(Pparse + '/{sample}/igblast_airr.tsv', sample=Lsample),
        changeo = expand(Pparse + '/{sample}/changeo_clone-pass_germ-pass.tsv', sample=Lsample),
        aa_fa = expand(Pparse + '/{sample}/seq_aa.fasta', sample=Lsample),
        anarci_H = expand(Pparse + '/{sample}/anarci_H.csv', sample=Lsample)
    
    
rule igblast:
    input: fa = indir + '/{sample}.fa',
    output:
        igblast_tsv = Pparse + '/{sample}/igblast_airr.tsv',
        igblast_txt = Pparse + '/{sample}/igblast_blast.txt',
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

rule changeo:
    input:
        fa = rules.igblast.input.fa,
        igblast_txt = rules.igblast.output.igblast_txt
    output: 
        changeo_db = Pparse + '/{sample}/changeo_db-pass.tsv',
        changeo_clone = Pparse + '/{sample}/changeo_clone-pass.tsv',
        changeo_fail = Pparse + '/{sample}/changeo_clone-fail.tsv',
        changeo  = Pparse + '/{sample}/changeo_clone-pass_germ-pass.tsv',
    log: e = Plog + '/changeo/{sample}.e', o = Plog + '/changeo/{sample}.o'
    benchmark: Plog + '/changeo/{sample}.bmk'
    resources: cpus=config['changeo_cpus']
    conda: 'envs/env.yaml'
    params: 
        changeo_VB_ref = f'ref/igblast-ref/germline_db/{species}-VDJB/IGV.fasta',
        changeo_DB_ref = f'ref/igblast-ref/germline_db/{species}-VDJB/IGD.fasta',
        changeo_JB_ref = f'ref/igblast-ref/germline_db/{species}-VDJB/IGJ.fasta',
    shell:"""
        MakeDb.py igblast -i {input.igblast_txt} \\
            -r {params.changeo_VB_ref} {params.changeo_DB_ref} {params.changeo_JB_ref} \\
            -s {input.fa} --outname changeo --regions default \\
            --failed --partial --format airr --extended 1>>{log.o} 2>>{log.e}
        
        DefineClones.py -d {output.changeo_db} --failed --act set --nproc {resources.cpus}\\
            --outname changeo --model ham --norm len --dist 0.15 1>>{log.o} 2>>{log.e}

        # in case all contigs are failed or all passed
        touch {output.changeo_clone}
        touch {output.changeo_fail}

        CreateGermlines.py -d {output.changeo_clone} -g dmask --cloned \\
            -r {params.changeo_VB_ref} {params.changeo_DB_ref} {params.changeo_JB_ref} \\
            1>>{log.o} 2>>{log.e}
        """

rule gen_aa:
    input: igblast_tsv = rules.igblast.output.igblast_tsv
    output: 
        aa_fa = Pparse + '/{sample}/seq_aa.fasta',
        gm_aa_fa = Pparse + '/{sample}/seq_gm_aa.fasta'
    log: notebook = Plog + '/gen_aa/{sample}.r.ipynb', 
            e = Plog + '/gen_aa/{sample}.e', o = Plog + '/gen_aa/{sample}.o'
    benchmark: Plog + '/gen_aa/{sample}.bmk'
    resources: cpus=1
    conda: 'envs/env.yaml'
    notebook: 'src/gen_aa.r.ipynb'

rule anarci:
    input: 
        aa_fa = rules.gen_aa.output.aa_fa,
        gm_aa_fa = rules.gen_aa.output.gm_aa_fa
    output:
        anarci_H = Pparse + '/{sample}/anarci_H.csv',
        anarci_KL = Pparse + '/{sample}/anarci_KL.csv',
        gm_anarci_H = Pparse + '/{sample}/gm_anarci_H.csv',
        gm_anarci_KL = Pparse + '/{sample}/gm_anarci_KL.csv',
        ext_anarci_H = Pparse + '/{sample}/ext_anarci_H.csv',
        ext_anarci_KL = Pparse + '/{sample}/ext_anarci_KL.csv',
    log: e = Plog + '/anarci/{sample}.e', o = Plog + '/anarci/{sample}.o'
    benchmark: Plog + '/anarci/{sample}.bmk'
    resources: cpus=config['anarci_cpus']
    params: 
        outdir = Pparse + '/{sample}', 
        ext_numbering = config['ext_numbering'],
        use_species_cmd = f'--use_species {species}' if config['species'] in ['human', 'mouse'] else ''
    conda: 'envs/env.yaml'
    shell:"""
        cd {params.outdir}
        # for mutation target profile 
        ANARCI -i {input.aa_fa} -o anarci -ht anarci_hittable.txt \\
            {params.use_species_cmd} --restrict ig -s imgt --csv --ncpu {resources.cpus} \\
            --assign_germline 1>>{log.o} 2>>{log.e}

        ANARCI -i {input.gm_aa_fa} -o gm_anarci -ht gm_anarci_hittable.txt \\
            {params.use_species_cmd} --restrict ig -s imgt --csv --ncpu {resources.cpus} \\
            --assign_germline 1>>{log.o} 2>>{log.e}

        # for extended numbering
        ANARCI -i {input.aa_fa} -o ext_anarci -ht ext_anarci_hittable.txt \\
            {params.use_species_cmd} --restrict ig -s {params.ext_numbering} --csv --ncpu {resources.cpus} \\
            --assign_germline 1>>{log.o} 2>>{log.e}
    """