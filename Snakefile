import os
import sys
import config
import program
import reference

def coverage(filename):
    os.system("module load samtools && samtools idxstats " + filename + " > " + filename + ".idx")
    if os.path.isfile(filename + ".idx"):
        fl = open(filename + ".idx")
        ll = fl.readlines()
        fl.close()
        tb = 0
        tr = 0
        for ln in ll:
            tb = tb + int(ln.strip().split("\t")[1])
            tr = tr + int(ln.strip().split("\t")[2])
        return int(tr*150/tb)


file = open("contrasts.txt", 'r')
combine = []
lln = file.readlines()
file.close()

for line in lln:
	combine.append(line.split(",")[0] + "_" + line.split(",")[1].strip("\n"))

rule all:
    input:
        expand("{comb}_tumor_damage_estimator.png", comb=combine),
        expand("QC/{comb}_tumor/qualimapReport.html", comb=combine),
        "QC/indexcov/index.html",
        expand("{comb}_conpair_conctamination.txt", comb=combine),
        expand("{comb}_tumor.ancestry.out", comb=combine),
        expand("SV/{comb}/gridss/somatic.vcf.assembly.bam", comb=combine),
        expand("SV/{comb}/lumpy/lumpy.somatic.vcf", comb=combine),
        expand("SV/{comb}/svaba/svaba.svaba.somatic.sv.vcf", comb=combine),
        expand("SV/{comb}/whamg/whamg.somatic.vcf", comb=combine),
        expand("SV/{comb}/manta/manta.somatic.vcf", comb=combine),
        expand("SV/{comb}/annotsv.input.bed", comb=combine),
        expand("SV/{comb}/annotsv/annot.sv.txt", comb=combine),
        expand("SV/{comb}/melt/tumor/SVA.final_comp.vcf", comb=combine),
        expand("SV/{comb}/melt/normal/SVA.final_comp.vcf", comb=combine),
        expand("SV/{comb}/msisensor/{comb}", comb=combine),
        expand("SV/{comb}/svviz/svviz.svg", comb=combine),
        expand("SV/{comb}/mavis/output/summary/mavis_summary_all_{comb}.tab", comb=combine),
        expand("CNV/{comb}/freec/{comb}_tumor.bam_minipileup.pileup", comb=combine),
        expand("CNV/{comb}/sequenza/seqz.bin1000.gz", comb=combine),
        expand("CNV/{comb}/facets/facets.cnv.tsv", comb=combine)


rule bbtools:
    input:
        ts = "{comb}_tumor.bam",
        ns = "{comb}.bam"
    output:
        s1 = temp("{comb}_tumor_5million.bam"),
        s2 = temp("{comb}_normal_5million.bam")
    log:
        log1 = "{comb}_tumor_bbtools.log",
        log2 = "{comb}_normal_bbtools.log"
    threads: 8
    params:
        batch = "-l nodes=1:ppn=8"
    shell:
        r'''
            module load samtools && \
            {bbtools} -Xmx64g in={input.ts} out={output.s1} samplereads=5000000 ref={bbtool_fa} 2>{log.log1} && \
            {bbtools} -Xmx64g in={input.ns} out={output.s2} samplereads=5000000 ref={bbtool_fa} 2>{log.log2}
        '''

rule split_reads:
    input:
        ts = "{comb}_tumor_5million.bam",
        ns = "{comb}_normal_5million.bam"
    output:
        R1 = temp("{comb}_tumor_R1_5million.mpileup"),
        R2 = temp("{comb}_tumor_R2_5million.mpileup"),
        R3 = temp("{comb}_normal_R1_5million.mpileup"),
        R4 = temp("{comb}_normal_R2_5million.mpileup")
    log:
        log1 = "{comb}_tumor_splitreads.log",
        log2 = "{comb}_normal_splitreads.log"
    params:
        batch = "-l nodes=1:ppn=8"
    shell:
        r'''
            module load samtools && \
            source /installed_tools/Damage-estimator/.source_damage && \
            /installed_tools/Damage-estimator/split_mapped_reads1.pl \
            -bam {input.ts} -genome {bbtool_fa} -mpileup1 {output.R1} -mpileup2 {output.R2} 2>{log.log1} && \
            /installed_tools/Damage-estimator/split_mapped_reads1.pl \
            -bam {input.ns} -genome {bbtool_fa} -mpileup1 {output.R3} -mpileup2 {output.R4} 2>{log.log2}
        '''

rule damage_estimator:
    input:
        R1 = "{comb}_tumor_R1_5million.mpileup",
        R2 = "{comb}_tumor_R2_5million.mpileup",
        R3 = "{comb}_normal_R1_5million.mpileup",
        R4 = "{comb}_normal_R2_5million.mpileup"
    output:
        txt1 = "{comb}_tumor_estimate.txt",
        plot1 = "{comb}_tumor_damage_estimator.png",
        txt2 = "{comb}_normal_estimate.txt",
        plot2 = "{comb}_normal_damage_estimator.png"
    log:
        log1 = "{comb}_tumor_estimate.log",
        log2 = "{comb}_normal_estimate.log"
    params:
        batch = "-l nodes=1:ppn=8,mem=64gb",
        prefix1 = "{comb}_tumor",
        prefix2 = "{comb}_normal"
    shell:
        r'''
            source /installed_tools/Damage-estimator/.source_damage && \
            /installed_tools/Damage-estimator/estimate_damage.pl \
            --mpileup1 {input.R1} --mpileup2 {input.R2} --id {params.prefix1} >{output.txt1} && \
            /installed_tools/Damage-estimator/plot_damage.R {output.txt1} {output.plot1} && \
            /installed_tools/Damage-estimator/estimate_damage.pl \
            --mpileup1 {input.R3} --mpileup2 {input.R4} --id {params.prefix2} >{output.txt2} && \
            /installed_tools/Damage-estimator/plot_damage.R {output.txt2} {output.plot2}
        '''

rule qualimap:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out1 = "QC/{comb}_tumor/qualimapReport.html",
        out2 = "QC/{comb}_normal/qualimapReport.html"
    params:
        batch = "-l nodes=1:ppn=8,mem=24gb",
        prefix1 = "QC/{comb}_tumor",
        prefix2 = "QC/{comb}_normal"
    shell:
        r'''
            {program.qualimap} bamqc -bam {input.file1} -ip -c gd {org} -gff {target_bed} -outdir {params.prefix1} \
            -outformat HTML  -nt 8 --java-mem-size=24G  -nw 500 -p NON-STRAND-SPECIFIC && \
            /opt/nasapps/development/qualimap/2.2/qualimap bamqc -bam {input.file2} -ip -c gd {org} -gff {target_bed} -outdir {params.prefix2} \
            -outformat HTML  -nt 8 --java-mem-size=24G  -nw 500 -p NON-STRAND-SPECIFIC
        '''

rule conpair1:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out1 = "{comb}_tumor.pileup",
        out2 = "{comb}.pileup"
    params:
        batch = "-l nodes=1:ppn=8"
    shell:
        r'''
            source /installed_tools/Conpair/source_conpair  && \
            /installed_tools/Conpair/scripts/run_gatk_pileup_for_sample.py \
            -B {input.file1} -O {output.out1} -D /installed_tools/Conpair \
            -R  /installed_tools/Conpair/data/GRCh38/hg38.fa \
            -M /installed_tools/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed \
            -G {program.gatk} -t ./ --xmx_java=96g && \
            /installed_tools/Conpair/scripts/run_gatk_pileup_for_sample.py \
            -B {input.file2} -O {output.out2} -D /installed_tools/Conpair \
            -R  /installed_tools/Conpair/data/GRCh38/hg38.fa \
            -M /installed_tools/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed \
            -G {program.gatk} -t ./ --xmx_java=96g			
        '''

rule conpair2:
    input:
        file1 = "{comb}_tumor.pileup",
        file2 = "{comb}.pileup"
    output:
        out3 = "{comb}_conpair_conctamination.txt",
        out4 = "{comb}_concordance.txt"
    params:
        batch = "-l nodes=1:ppn=8"
    shell:
        r'''
            source /installed_tools/Conpair/source_conpair  && \
            python /installed_tools/Conpair/scripts/estimate_tumor_normal_contamination.py  \
            -T {input.file1} -N {input.file2} -D /installed_tools/Conpair \
            -M /installed_tools/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt \
            -O {output.out3} && \
            python /installed_tools/Conpair/scripts/estimate_tumor_normal_contamination.py  \
            -T {input.file1} -N {input.file2} -D /installed_tools/Conpair \
            -M /installed_tools/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt \
            -O {output.out4}			
        '''
		
rule ancestry:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out1 = "{comb}_tumor.ancestry.out",
        out2 = "{comb}.ancestry.out"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix1 = "{comb}_tumor",
        prefix2 = "{comb}"
    shell:
        r'''
            export PATH=/installed_tools/Anaconda/2.7/install/bin:$PATH && \
            python /installed_tools/ancestry/runancestry.py \
            -f {reference.fregs} \
            --bam {input.file1} -o {params.prefix1} \
            --path /installed_tools/ancestry --addchr -c 8 && \
            python /installed_tools/ancestry/runancestry.py \
            -f {reference.fregs} \
            --bam {input.file2} -o {params.prefix2} \
            --path /installed_tools/ancestry --addchr -c 8
        '''

rule indexcov:
    input:
        expand("{comb}_tumor.bam", comb = combine),
        expand("{comb}.bam", comb = combine)
    output:
        out = "QC/indexcov/index.html"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix = "QC/indexcov"
    shell:
        r'''
            export GOROOT=/installed_tools/go_binary/go1.9.1 && \
            export PATH==$PATH:$GOROOT/bin && \
            export PATH=$PATH:/installed_tools/go_left && \
            goleft  indexcov -d {params.prefix} --sex chrX,chrY --fai {reference.fai} {input}
        '''

rule gridss:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out = "SV/{comb}/gridss/gridss.somatic.vcf",
        out1 = "SV/{comb}/gridss/gridss.somatic.tumor.specific.vcf",
        outbam = "SV/{comb}/gridss/somatic.vcf.assembly.bam"
    log:
        log1 = "SV/{comb}/gridss/run.gridss.err",
        log2 = "SV/{comb}/gridss/run.gridss.log"
    params:
        batch = "-l nodes=1:ppn=16",
        prefix = "SV/{comb}/gridss"
    shell:
        r'''
            cd {params.prefix}
            module load bwa && \
            module load bedtools && \
            module load samtools && \
            java -Xmx128g -Dgridss.gridss.output_to_temp_file=true -Dsamjdk.create_index=true -Dsamjdk.use_async_io_read_samtools=true \
            -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=true -Dsamjdk.compression_level=1 \
            -cp /installed_tools/gridss/gridss-1.5.1-jar-with-dependencies.jar gridss.CallVariants \
            TMP_DIR=./ WORKING_DIR=./ ASSEMBLY=somatic.vcf.assembly.bam REFERENCE_SEQUENCE={reference.bwa_seqc} \
            INPUT= ../../../{input.file1} INPUT=../../../{input.file2} OUTPUT=gridss.somatic.vcf IGNORE_DUPLICATES=TRUE THREADS=16 2>run.gridss.err 1>run.gridss.log && \
            export PATH=/opt/nasapps/development/R/3.4.0/bin/:$PATH && \
            export R_LIBS=/installed_tools/R && \
            Rscript {gridss_script} gridss.somatic.vcf
        '''

rule lumpy:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out1 = "SV/{comb}/lumpy/lumpy.somatic.vcf",
        out2 = "SV/{comb}/lumpy/lumpy.somatic.tumor.only.vcf"
    log:
        log1 = "SV/{comb}/lumpy/run.lumpy.err",
        log2 = "SV/{comb}/lumpy/run.lumpy.log"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix = "SV/{comb}/lumpy"
    shell:
        r'''
            cd {params.prefix} && module load samtools && module load sambamba && \
            export PATH=/installed_tools/samblaster:$PATH && \
            /installed_tools/lumpy-sv/bin/lumpyexpress -B ../../../{input.file1},../../../{input.file2} \
            -R {reference.bwa_seqc} \
            -o lumpy.somatic.vcf -T ./ 2>run.lumpy.err 1>run.lumpy.log && \
            python /installed_tools/scripts/pick_vcf.py \
			lumpy.somatic.vcf lumpy.somatic.tumor.only.vcf lumpy
        '''

rule svaba:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out = "SV/{comb}/svaba/svaba.svaba.somatic.sv.vcf"
    log:
        log1 = "SV/{comb}/svaba/run.svaba.err",
        log2 = "SV/{comb}/svaba/run.svaba.log"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix = "SV/{comb}/svaba"
    shell:
        r'''
            cd {params.prefix}
            module load bwa && \
            module load bedtools && \
            module load samtools && \
            /installed_tools/svaba/bin/svaba run -t ../../../{input.file1} -n ../../../{input.file2} -a svaba \
            -p 8 -G {bwa_ref} -D {dbsnp_indel} 2>run.svaba.err 1>run.svaba.log
        '''

rule whamg:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out1 = "SV/{comb}/whamg/whamg.somatic.vcf",
        out2 = "SV/{comb}/whamg/whamg.somatic.tumor.only.vcf"
    log:
        log1 = "SV/{comb}/whamg/run.whamg.err",
        log2 = "SV/{comb}/whamg/run.whamg.log"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix = "SV/{comb}/whamg/list_{comb}"
    run:
        kee = open(params.prefix, 'w')
        kee.write(input.file1 + "\n" + input.file2)
        kee.close()
        if os.path.isfile(params.prefix):
            shell(r'''module load bwa && module load bedtools && module load samtools && \
            /installed_tools/wham/bin/whamg -x 8 -f {params.prefix} \
            -a {reference.bwa_seqc} > {output.out1} && \
            python /installed_tools/scripts/pick_vcf.py {output.out1} {output.out2} whamg''')

rule manta:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out = "SV/{comb}/manta/manta.somatic.vcf"
    log:
        log1 = "SV/{comb}/manta/run.manta.err",
        log2 = "SV/{comb}/manta/run.manta.log"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix = "SV/{comb}/manta"
    shell:
        r'''
            /installed_tools/manta/1.3.1/bin/configManta.py --referenceFasta={reference.bwa_seqc} \
            --bam={input.file2} --tumorBam={input.file1} --runDir={params.prefix} && cd {params.prefix} && \
            export PATH=/installed_tools/Anaconda/2.7/install/bin/:$PATH && \
            python runWorkflow.py -m local -j 8 --memGb=64 2>run.manta.err 1>run.manta.log && \
            gunzip -c results/variants/somaticSV.vcf.gz > manta.somatic.vcf
        '''

rule survivor:
    input:
        gridss_vcf = "SV/{comb}/gridss/gridss.somatic.tumor.specific.vcf",
        lumpy_vcf = "SV/{comb}/lumpy/lumpy.somatic.tumor.only.vcf",
        svaba_vcf = "SV/{comb}/svaba/svaba.svaba.somatic.sv.vcf",
        whamg_vcf = "SV/{comb}/whamg/whamg.somatic.tumor.only.vcf",
        manta_vcf = "SV/{comb}/manta/manta.somatic.vcf"
    output:
        outvcf = "SV/{comb}/survivor.vcf",
        outbed = "SV/{comb}/annotsv.input.bed"
    log:
        log1 = "SV/{comb}/run.survivor.err",
        log2 = "SV/{comb}/run.survivor.log"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix = "list_sv_{comb}"
    run:
        ke = open(params.prefix, 'w')
        ke.write(input.gridss_vcf + "\n" + lumpy_vcf + "\n" + input.svaba_vcf + "\n" + input.whamg_vcf + "\n" + input.manta_vcf)
        ke.close()
        if os.path.isfile(params.prefix):
            shell(r'''/installed_tools/SURVIVOR/SURVIVOR-1.0.3/Debug/SURVIVOR merge {params.prefix} 500 2 0 0 0 50 {output.outvcf} 2>{log.log1} 1>{log.log2} && \
            /installed_tools/SURVIVOR/SURVIVOR-1.0.3/Debug/SURVIVOR vcftobed {output.outvcf} 0 1000000000  {output.outbed}''')

rule annotsv:
    input:
        file1 = "SV/{comb}/annotsv.input.bed",
        #file2 = "{comb}.vcf.gz"
    output:
        out = "SV/{comb}/annotsv/annot.sv.txt"
    log:
        log1 = "SV/{comb}/annotsv/run.annotsv.err",
        log2 = "SV/{comb}/annotsv/run.annotsv.log"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix1 = "SV/{comb}/annotsv",
        prefix2 = "annot.sv.txt"
    shell:
        r'''
            export PATH=/installed_tools/Anaconda/3.6/install/bin/:$PATH && \
            export PYTHONPATH=/installed_tools/Anaconda/3.6/install/lib/python3.6/site-packages/:$PYTHONPATH && \
            export ANNOTSV=/installed_tools/annotsv/AnnotSV_1.1.1 && \
            /installed_tools/annotsv/AnnotSV_1.1.1/bin/AnnotSV \
            -SVinputFile {input.file1} -bedtools \
            /opt/nasapps/production/bedtools/2.26.0/bin/bedtools \
            -genomeBuild GRCh38 -outputDir {params.prefix1} -outputFile {params.prefix2} 2>{log.log1} 1>{log.log2}
        '''

rule melt_t:
    input:
        file1 = "{comb}_tumor.bam"
    output:
        out1 = "SV/{comb}/melt/tumor/SVA.final_comp.vcf"
    log:
        log1 = "SV/{comb}/melt/tumor/run.melt.err",
        log2 = "SV/{comb}/melt/tumor/run.melt.log"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix1 = "SV/{comb}/melt/tumor",
        prefix2 = "SV/{comb}/melt/tumor/{comb}_tumor.bam.idx"
    run:
        #ntumor = coverage (input.file1)
        ntumor = 0
        shell(r'''module load samtools && samtools idxstats {input.file1} > {params.prefix2}''')
        if os.path.isfile(params.prefix2):
            fl = open(params.prefix2)
            ll = fl.readlines()
            fl.close()
            tb = 0
            tr = 0
            for ln in ll:
                tb = tb + int(ln.strip().split("\t")[1])
                tr = tr + int(ln.strip().split("\t")[2])
            ntumor = int(tr*150/tb)
            if ntumor == 0:
                ntumor = 50
        shell(r'''java -Xmx31g -jar /installed_tools/MELT/MELTv2.1.5/MELT.jar Single \
        -bamfile {input.file1} \
        -a -e 450 -h {reference.bwa_seqc} \
        -bowtie /installed_tools/bowtie2/2.2.6/bin/bowtie2 \
        -t /installed_tools/MELT/MELTv2.1.5/me_refs/Hg38/list.txt \
        -n /installed_tools/MELT/MELTv2.1.5/add_bed_files/Hg38/Hg38.genes.bed \
        -w {params.prefix1} -c {ntumor} -z 10000 2>{log.log1} 1>{log.log2}''')

rule melt_n:
    input:
        file2 = "{comb}.bam"
    output:
        out2 = "SV/{comb}/melt/normal/SVA.final_comp.vcf"
    log:
        log3 = "SV/{comb}/melt/normal/run.melt.err",
        log4 = "SV/{comb}/melt/normal/run.melt.log"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix2 = "SV/{comb}/melt/normal",
        prefix3 = "SV/{comb}/melt/normal/{comb}.bam.idx"
    run:
        #ntumor = coverage (input.file1)
        nnormal = 0
        shell(r'''module load samtools && samtools idxstats {input.file2} > {params.prefix3}''')
        if os.path.isfile(params.prefix3):
            fl = open(params.prefix3)
            ll = fl.readlines()
            fl.close()
            tb = 0
            tr = 0
            for ln in ll:
                tb = tb + int(ln.strip().split("\t")[1])
                tr = tr + int(ln.strip().split("\t")[2])
            nnormal = int(tr*150/tb)
            if nnormal == 0:
                nnormal = 50
        shell(r'''java -Xmx31g -jar /installed_tools/MELT/MELTv2.1.5/MELT.jar Single \
        -bamfile {input.file2} \
        -a -e 450 -h {reference.bwa_seqc} \
        -bowtie /installed_tools/bowtie2/2.2.6/bin/bowtie2 \
        -t /installed_tools/MELT/MELTv2.1.5/me_refs/Hg38/list.txt \
        -n /installed_tools/MELT/MELTv2.1.5/add_bed_files/Hg38/Hg38.genes.bed \
        -w {params.prefix2} -c {nnormal} -z 10000 2>{log.log3} 1>{log.log4}''')

rule msisensor:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out = "SV/{comb}/msisensor/{comb}"
    log:
        log1 = "SV/{comb}/msisensor/run.msisensor.err",
        log2 = "SV/{comb}/msisensor/run.msisensor.log"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix = "SV/{comb}/msisensor"
    shell:
        r'''
            /installed_tools/msisensor/msisensor msi \
            -d /installed_tools/RefGenomes/hg38_SEQC/hg38.filter.msisensor.txt \
            -n {input.file2} -t {input.file1} -o {output.out} -c 15 -b 16 2>{log.log1} 1>{log.log2}
        '''
rule svviz:
    input:
        file1 = "SV/{comb}/survivor.vcf",
        file2 = "{comb}_tumor.bam"
    output:
        out = "SV/{comb}/svviz/svviz.svg"
    log:
        log1 = "SV/{comb}/svviz/run.svviz.err",
        log2 = "SV/{comb}/svviz/run.svviz.log"
    params:
        batch = "-l nodes=1:ppn=32",
        prefix1 = "SV/{comb}/svviz/svviz_tra_summary.tsv",
        prefix2 = "SV/{comb}/svviz/svviz_tra_input.vcf"
    shell:
        r'''
            grep -w -E '#|TRA' {input.file1}  | sed 's/STRANDS/STRAND/g' > {params.prefix2} && \
            source /installed_tools/svviz/svviz1/source_svviz  && \
            python  /installed_tools/svviz/svviz1/bin/svviz \
            --type batch --summary {params.prefix1} -b {input.file2} \
            {reference.bwa_seqc} {params.prefix2} \
            --processes 32 --skip-cigar --fast -A {reference.svviz_gtf} \
            -e {output.out} --format svg 2>{log.log1} 1>{log.log2}
        '''

rule mavis:
    input:
        file1 = "SV/{comb}/annotsv.input.bed",
        file2 = "{comb}_tumor.bam"
    output:
        out = "SV/{comb}/mavis/output/summary/mavis_summary_all_{comb}.tab"
    log:
        log1 = "SV/{comb}/mavis/run.mavis.err",
        log2 = "SV/{comb}/mavis/run.mavis.log"
    params:
        batch = "-l nodes=1:ppn=8",
        prefix1 = "SV/{comb}/mavis/tmp.in.mavis.tsv",
        prefix2 = "/installed_tools/script/header",
        prefix3 = "SV/{comb}/mavis/mavis.input.tsv",
        prefix4 = "SV/{comb}/mavis/mavis.cfg",
        prefix5 = "SV/{comb}/mavis",
        prefix6 = "{comb}"
    shell:
        r'''
            cut -f1-6 {input.file1} > {params.prefix1} && cat {params.prefix2} {params.prefix1} > {params.prefix3} && \
            module load blat && \
            export PYTHONPATH=/installed_tools/svviz/svviz2/lib/python3.6/site-packages/ && \
            export PATH=/installed_tools/Anaconda/3.6/install/bin/:$PATH && \
            export MAVIS_ALIGNER='blat' && \
            export MAVIS_ALIGNER_REFERENCE=/installed_tools/mavis/hg38/hg38.2bit && \
            export MAVIS_DGV_ANNOTATION=/installed_tools/mavis/hg38/dgv_hg38_variants.tab && \
            export MAVIS_ANNOTATIONS=/installed_tools/mavis/hg38/ensembl79_hg38_annotations.json && \
            export MAVIS_TEMPLATE_METADATA=/installed_tools/mavis/hg38/cytoBand.txt && \
            export MAVIS_REFERENCE_GENOME=/installed_tools/mavis/hg38/hg38.fa && \
            export MAVIS_MASKING=/installed_tools/mavis/hg38/GRCh38_masking.tab && \
            export min_clusters_per_file=100 && \
            export MAVIS_SCHEDULER=LOCAL && \
            export MAVIS_CONCURRENCY_LIMIT=8 && \
            python /installed_tools/mavis/2.1.1/bin/mavis config --library {params.prefix6} genome diseased FALSE {input.file2} \
            --input {params.prefix3} {params.prefix6} -w {params.prefix4} && \
            python /installed_tools/mavis/2.1.1/bin/mavis setup {params.prefix4} --skip_stage cluster --skip_stage validate -o {params.prefix5} && \
            python /installed_tools/mavis/2.1.1/bin/mavis schedule -o {params.prefix5} --submit 2>{log.log1} 1>{log.log2}
        '''

rule freec:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out1 = "CNV/{comb}/freec/{comb}_tumor.bam_minipileup.pileup",
        out2 = "CNV/{comb}/freec/{comb}.bam_minipileup.pileup"
    log:
        log1 = "CNV/{comb}/freec/run.freec.err",
        log2 = "CNV/{comb}/freec/run.freec.log"
    params:
        batch = "-l nodes=1:ppn=24",
        prefix1 = "CNV/{comb}/freec",
        prefix2 = "CNV/{comb}/freec/freec.config"
    run:
        lk = open (freec_path,'r')
        lkall = lk.read()
        lk.close()
        lkall = lkall.replace("$tumor", "../../../" + input.file1)
        lkall = lkall.replace("$normal", "../../../" + input.file2)
        lk1 = open(params.prefix2, 'w')
        lk1.write(lkall)
        lk1.close()
        if os.path.isfile(params.prefix2):
            shell(r''' module load samtools && module load bedtools && module load sambamba && \
			cd {params.prefix1} && /installed_tools/FREEC/FREEC-11.5/src/freec -conf freec.config 2>run.freec.err 1>run.freec.err''')

rule facets:
    input:
        file1 = "{comb}_tumor.bam",
        file2 = "{comb}.bam"
    output:
        out1 = "CNV/{comb}/facets/facets.cnv.gz",
        out2 = "CNV/{comb}/facets/facets.cnv.tsv"
    log:
        log1 = "CNV/{comb}/facets/run.facets.err",
        log2 = "CNV/{comb}/facets/run.facets.log"
    params:
        batch = "-l nodes=1:ppn=24",
        prefix1 = "CNV/{comb}/facets"
    shell:
        r'''
            cd {params.prefix1} && export LD_LIBRARY_PATH=/installed_tools/htslib/lib:$PATH && \
            /installed_tools/facets/inst/extcode/snp-pileup \
            -g -q 20 -Q 20 -r 15,0 /is2/projects/CCR-SF/active/RefGenomes/hg38_SEQC/facets_common_all_20170710.vcf \
            facets.cnv ../../../{input.file2} ../../../{input.file1} 2>run.facets.err 1>run.facets.log && \
            export PATH=/installed_tools/R/3.4.0/bin/:$PATH && \
            export R_LIBS=/installed_tools/R && \
            Rscript {facets_script} facets.cnv.gz
        '''

rule sequenza:
    input:
        file1 = "CNV/{comb}/freec/{comb}_tumor.bam_minipileup.pileup",
        file2 = "CNV/{comb}/freec/{comb}.bam_minipileup.pileup"
    output:
        out1 = "CNV/{comb}/sequenza/seqz.gz",
        out2 = "CNV/{comb}/sequenza/seqz.bin1000.gz"
    params:
        batch = "-l nodes=1:ppn=32",
        prefix1 = "CNV/{comb}/sequenza"
    shell:
        r'''
            export PATH=/installed_tools/samtools/1.7/bin/:/installed_tools/Anaconda/3.6/install/bin/:$PATH && \
            module load samtools && cd {params.prefix1} && \
            sequenza-utils bam2seqz -p -n ../../../{input.file2} -t ../../../{input.file1} \
            -gc /installed_tools/RefGenomes/hg38_SEQC/sequenza_gc_hg38.txt \
            -F /installed_tools/RefGenomes/hg38_SEQC/hg38.filter.fa | gzip > seqz.gz && \
            sequenza-utils seqz_binning -w 1000 -s seqz.gz | gzip > seqz.bin1000.gz && \
            export PATH=/installed_tools/R/3.4.0/bin/:$PATH && \
            export R_LIBS=/installed_tools/R && \
            Rscript {sequenza_script} seqz.bin1000.gz
        '''
