samples = ["mutant"]
ref = "/home/marc/bwgs/refs/wildtype.fna"
gff = "/home/marc/bwgs/refs/wildtype.gff"
kraken2 = "/home/marc/Desktop/DBs/kraken2/pluspfp"
indir  = "/home/marc/bwgs/samples"

def get_final_output(wildcards):
    final_output = []
    
    final_output.extend(expand(
            "{dir}/{sample}_stats.txt",
            sample = samples, dir = ["01_map_to_ref", "05_assembly-mapping"]
    ))

    final_output.extend(expand(
            "04_contamination/{dir}/{sample}.report",
            sample = samples, dir = ["unmapped_reads", "filttered"]
    ))

    return final_output

rule QC:
    input:
        R1=f"{indir}/{{sample}}_R1.fastq",
        R2=f"{indir}/{{sample}}_R2.fastq"
    output:
        R1="00_filterred/{sample}_R1.fastq",
        R2="00_filterred/{sample}_R2.fastq",
        rep="00_filterred/{sample}_R2.html"
    threads: 16
    shell:
        """
        fastp --in1 {input.R1} --in2 {input.R2} \
            --out1 {output.R1} --out2 {output.R2} \
            --thread {threads} -h {output.rep}
        """

rule index_ref:
    input:
        ref
    output:
        directory('bowtie2_index/')
    threads: 20
    shell:
        """
        mkdir -p bowtie2_index
        bowtie2-build --threads {threads} {input} bowtie2_index/Org_ref
        """

rule map_to_ref:
    input:
        rules.index_ref.output,
        R1="00_filterred/{sample}_R1.fastq",
        R2="00_filterred/{sample}_R2.fastq",
    output:
        "01_map_to_ref/{sample}.sam"
    threads: 20
    shell:
        """
        bowtie2 --threads {threads} -x bowtie2_index/Org_ref \
                -1 {input.R1} \
                -2 {input.R2} \
                -S {output}
        """ 

rule fixmate:
    input: "01_map_to_ref/{sample}.sam"
    output: "01_map_to_ref/{sample}_fixmate.sam"
    shell: 
        """
        samtools sort -n \
                -O sam {input} | samtools fixmate \
                -m -O bam - {output}
        """
rule sort:
    input: "01_map_to_ref/{sample}_fixmate.sam"
    output: "01_map_to_ref/{sample}_fixmate.bam"
    shell: 
        """
        samtools sort -O bam {input} \
                -o {output}
        """

rule remove_duplicate:
    input: "01_map_to_ref/{sample}_fixmate.bam"
    output: "01_map_to_ref/{sample}_dedup.bam"
    shell: "samtools markdup -r -S {input} {output}"

rule mapping_stats:
    input: "01_map_to_ref/{sample}_dedup.bam"
    output: 
        sam="01_map_to_ref/{sample}_stats.txt",
        quali=directory("01_map_to_ref/{sample}_stats")
    threads: 2
    shell:
        """
        samtools flagstat {input} > {output.sam} &
        qualimap bamqc -bam {input} -outdir {output.quali}
        wait
        """

rule extract_mapped:
    input: "01_map_to_ref/{sample}_dedup.bam"
    output: 
        bam="01_map_to_ref/{sample}_concordant.bam",
        R1="02_Extract_reads/mapped_reads/{sample}_R1.fastq",
        R2="02_Extract_reads/mapped_reads/{sample}_R2.fastq",
        U="02_Extract_reads/mapped_reads/{sample}_U.fastq"
    shell: 
        """
        samtools view -b -f 3 {input} > 01_map_to_ref/{wildcards.sample}_concordant.bam
        samtools fastq -1 {output.R1} \
            -2 {output.R2} \
            -0 {output.U} 01_map_to_ref/{wildcards.sample}_concordant.bam
        """

rule extract_unmapped:
    input: "01_map_to_ref/{sample}_dedup.bam"
    output: 
        bam="01_map_to_ref/{sample}_unmapped.bam",
        R1="02_Extract_reads/unmapped_reads/{sample}_R1.fastq",
        R2="02_Extract_reads/unmapped_reads/{sample}_R2.fastq",
        U="02_Extract_reads/unmapped_reads/{sample}_U.fastq"
    shell: 
        """
        samtools view -b -f 4 {input} > 01_map_to_ref/{wildcards.sample}_unmapped.bam
        samtools fastq -1 {output.R1} \
            -2 {output.R2} \
            -0 {output.U} 01_map_to_ref/{wildcards.sample}_unmapped.bam
        """

rule assemble_RAW_reads:
    input:
        R1="00_filterred/{sample}_R1.fastq",
        R2="00_filterred/{sample}_R2.fastq"
    output:
        eldir=directory("03_Assembly/RAW/{sample}_filterd-Assembly"),
        asmpl="03_Assembly/RAW/{sample}_filterd-Assembly/assembly.fasta"
    threads: 32
    shell:
        """
    unicycler -1 {input.R1} -2 {input.R2} \
              -o {output.eldir} \
              --threads {threads} 
        """

rule assmble_qc_RAW_reads:
    input:
        "03_Assembly/RAW/{sample}_filterd-Assembly/assembly.fasta"
    output:
        directory("03_Assembly/RAW/{sample}_filterd-Assembly/QC")
    shell:
        "quast -o {output} {input}"

rule assemble_mapped_reads:
    input:
        R1="02_Extract_reads/mapped_reads/{sample}_R1.fastq",
        R2="02_Extract_reads/mapped_reads/{sample}_R2.fastq"
    output:
        eldir=directory("03_Assembly/mapped_reads/{sample}_mapped_reads-Assembly"),
        asmpl="03_Assembly/mapped_reads/{sample}_mapped_reads-Assembly/assembly.fasta"
    threads: 32
    shell:
        """
    unicycler -1 {input.R1} -2 {input.R2} \
              -o {output.eldir} \
              --threads {threads} 
        """

rule assmble_qc_mapped_reads:
    input:
        "03_Assembly/mapped_reads/{sample}_mapped_reads-Assembly/assembly.fasta"
    output:
        directory("03_Assembly/mapped_reads/{sample}_mapped_reads-Assembly/QC")
    shell:
        "quast -o {output} {input}"


rule kraken_RAW_reads:
    input:
        R1="00_filterred/{sample}_R1.fastq",
        R2="00_filterred/{sample}_R2.fastq"
    output:
        report="04_contamination/filttered/{sample}.report",
        log="04_contamination/filttered/{sample}.out"
    params:
        db=kraken2
    threads: 32
    shell:
        """
        kraken2 --db {params.db} --threads {threads} \
            --report {output.report} \
            --paired {input.R1} {input.R2} > {output.log}
        """

rule kraken_unmapped_reads:
    input:
        R1="02_Extract_reads/unmapped_reads/{sample}_R1.fastq",
        R2="02_Extract_reads/unmapped_reads/{sample}_R2.fastq"
    output:
        report="04_contamination/unmapped_reads/{sample}.report",
        log="04_contamination/unmapped_reads/{sample}.out"
    params:
        db=kraken2
    threads: 32
    shell:
        """
        kraken2 --db {params.db} --threads {threads} \
            --report {output.report} \
            --paired {input.R1} {input.R2} > {output.log}
        """

rule index_ref_assembly:
    input:
        "03_Assembly/RAW/{sample}_filterd-Assembly/assembly.fasta"
    output:
        directory("03_Assembly/RAW/{sample}_filterd-Assembly/bowtie2_index")
    threads: 32
    shell:
        """
        mkdir -p 03_Assembly/RAW/{wildcards.sample}_filterd-Assembly/bowtie2_index/
        bowtie2-build --threads {threads} {input} 03_Assembly/RAW/{wildcards.sample}_filterd-Assembly/bowtie2_index/assembly
        """

rule map_to_assembly:
    input:
        rules.index_ref_assembly.output,
        R1="00_filterred/{sample}_R1.fastq",
        R2="00_filterred/{sample}_R2.fastq",
    output:
        "05_assembly-mapping/{sample}.sam"
    threads: 32
    shell:
        """
        bowtie2 --threads {threads} -x 03_Assembly/RAW/{wildcards.sample}_filterd-Assembly/bowtie2_index/assembly \
                -1 {input.R1} \
                -2 {input.R2} \
                -S {output}
        """ 

rule fixmate_assembly:
    input: "05_assembly-mapping/{sample}.sam"
    output: "05_assembly-mapping/{sample}_fixmate.sam"
    shell: 
        """
        samtools sort -n \
                -O sam {input} | samtools fixmate \
                -m -O bam - {output}
        """

rule sort_assembly:
    input: "05_assembly-mapping/{sample}_fixmate.sam"
    output: "05_assembly-mapping/{sample}_fixmate.bam"
    shell: 
        """
        samtools sort -O bam {input} \
                -o {output}
        """

rule remove_duplicate_assembly:
    input: "05_assembly-mapping/{sample}_fixmate.bam"
    output: "05_assembly-mapping/{sample}_dedup.bam"
    shell: "samtools markdup -r -S {input} {output}"


rule mapping_stats_assembly:
    input: "05_assembly-mapping/{sample}_dedup.bam"
    output: 
        sam="05_assembly-mapping/{sample}_stats.txt",
        quali=directory("05_assembly-mapping/{sample}_stats")
    threads: 2
    shell:
        """
        samtools flagstat {input} > {output.sam} &
        qualimap bamqc -bam {input} -outdir {output.quali}
        wait
        """

        
rule multiqc:
    input:
        get_final_output
    
    conda: "../env/wes_gatk.yml"

    benchmark: "benchamrks/Multiqc/report.txt"

    output:
        "multiqc/multiqc_report.html"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        "multiqc . -o multiqc/"

