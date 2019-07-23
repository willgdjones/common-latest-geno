#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/beagle5
========================================================================================
 lifebit-ai/beagle5 Imputation Pipeline.
 #### Homepage / Documentation
 https://github.com/lifebit-ai/beagle5
----------------------------------------------------------------------------------------
*/

params.panel_path = "s3://lifebit-featured-datasets/modules/beagle5/ref"
params.fasta_path = "s3://lifebit-featured-datasets/modules/beagle5/hs37d5.fa.gz"
params.genotypes_path = "s3://lifebit-featured-datasets/modules/genotypes"
params.chromosome_regions = "s3://lifebit-featured-datasets/beagle5/modules/hs37d5.txt"
params.regex = ~/chr(.*).1kg.phase3.v5a.b37.bref3/

/*--------------------------------------------------
  Genotype input files from kindof 23AndMe style
---------------------------------------------------*/

if(params.genotypes_path) {
  if(params.genotypes_path.endsWith(".txt")) {

     Channel
        .fromPath( "${params.genotypes_path}" )
        .map { row -> [ file(row).baseName, [ file( row ) ] ] }
        .ifEmpty { exit 1, "${params.genotypes_path} not found"}
        .set { genoChannel }

  } else {
    Channel
        .fromFilePairs("${params.genotypes_path}/*", size: 1)
        .ifEmpty { exit 1, "${params.genotypes_path}/*.txt not found"}
        .set{genoChannel}
  }
} else {
  exit 1, "please specify --genotypes OR --genotypes_path"
}

/*--------------------------------------------------
  Genome panel reference files from 1k Genomes
---------------------------------------------------*/

if(params.panel_path) {
    Channel
      .fromFilePairs("${params.panel_path}/chr*", size: 1){ file ->
                def chr_match = file.getName() =~ params.regex
                chr_match[0][1]
        }
      .filter{ it ->
        def match = it[0] =~ ~/^[0-9]*/
        match[0]
      }
      .ifEmpty { exit 1, "${params.panel_path} not found"}
      .set{referencePanel}
} else {
  exit 1, "please specify --panel_path"
}

/*--------------------------------------------------
  Genome panel reference regions to impute
---------------------------------------------------*/

if(params.chromosome_regions) {
    Channel
      .fromPath(params.chromosome_regions)
      .splitText(by: 1)
      .map{ line ->
        def chr = line.split(':')
        [chr[0], line]
      }
      .ifEmpty { exit 1, "${params.chromosome_regions} not found"}
      .phase(referencePanel)
      .map{ chr_regions ->
        [chr_regions[0][0], chr_regions[0][1].strip(), chr_regions[1][1]]
      }
      .set{regionsReferencePanel}
} else {
  exit 1, "please specify --chromosome_regions"
}


/*--------------------------------------------------
  Validate if FASTA file is to be uncompressed
---------------------------------------------------*/

if("${params.fasta_path}".endsWith(".gz")){

  Channel
    .fromPath("${params.fasta_path}")
    .ifEmpty { exit 1, "Cannot find genome fasta assemby: ${params.fasta_path}"}
    .set {genome2unzip}

  process unzip {
    input:
    file(fasta_file) from genome2unzip

    output:
    file('*.fa') into (genomeToBeIndexed, genomeAssembly)

    script:
    """
    gunzip -d -f ${fasta_file}
    """

  }

} else {
  if(params.fasta_path) {
    Channel
      .fromPath(params.fasta_path)
      .ifEmpty { exit 1, "Genome fasta assembly file: ${params.fasta_path} not found"}
      .into{ genomeToBeIndexed; genomeAssembly }
  } else {
    exit 1, "please specify --fasta_path"
  }
}

/*--------------------------------------------------
  Index FASTA genome
---------------------------------------------------*/

process index {
  tag "${fasta_file}"
  container 'lifebitai/samtools'

  input:
  file(fasta_file) from genomeToBeIndexed

  output:
  file("${fasta_file}.fai") into genomeIndex
  file("${fasta_file}.vcf") into reheaderIndex

  script:
  """
  samtools faidx ${fasta_file}
  awk '{printf("##contig=<ID=%s,length=%d>\\n",\$1,\$2);}' ${fasta_file}.fai > ${fasta_file}.vcf
  """
}

/*--------------------------------------------------
  Convert from 23AndMe format to VCF
---------------------------------------------------*/

process bcftools {
  tag "${name}"
  container 'lifebitai/bcftools'

  input:
  set val(name), file(genotype_file) from genoChannel
  file fasta_file from genomeAssembly.collect()
  file fai from genomeIndex.collect()

  output:
  set val(name), file("${name}.vcf.gz") into vcfGenotypes

  script:
  """
  # convert dos files to unix
  sed -i 's/\r\$//' ${genotype_file}
  # convert ancestry files
  sed 's/23/X/g; s/24/Y/g; s/25/MT/g; s/26/X/g' ${genotype_file} > tmp.txt

  # for males where genotype contains one allele (non-autosomal chrs), duplicate the allele so that it's homozygous
  grep -v '^#' tmp.txt | awk '{print \$2}' | sort -u > chrs.txt
  if grep -q "Y" chrs.txt; then
    awk '{OFS="\t"; if (((length(\$4) == 1 && \$2 == "X"))) \$4=\$4\$4; print \$0}' tmp.txt > ${genotype_file}
  fi

  bcftools convert --tsv2vcf ${genotype_file}  -f ${fasta_file} -s $name -Oz -o ${name}.tmp.vcf.gz
  bcftools filter --set-GTs . -e 'FMT/GT="."' -Oz -o ${name}.filt.vcf.gz ${name}.tmp.vcf.gz
  bcftools view -t "^MT" -f PASS -Oz -o ${name}.vcf.gz ${name}.filt.vcf.gz
  """
}

regionsReferencePanel
  .combine(vcfGenotypes)
  .dump(tag: 'vcfGenotypes_here')
  .set{chrs2Genotype}


/*--------------------------------------------------
  Run imputation with Beagle5
---------------------------------------------------*/

process beagle {

  tag "${region},${chr_genotype_file}"
  container 'lifebitai/beagle5'

  input:
  set val(chr), val(region), file(ref_file), val(name), file(chr_genotype_file) from chrs2Genotype

  output:
  set val(name), file('*.imputed.vcf.gz') into imputedChrs

  script:
  if( !task.memory ){
    log.info "[Beagle] Available memory not known - defaulting to 2GB. Specify process memory requirements to change this."
    avail_mem = 2
  } else {
    avail_mem = task.memory.toGiga()
  }
  """
  beagle -Xmx${avail_mem}g ref=${ref_file} gt=${chr_genotype_file} out=${chr_genotype_file}.imputed chrom=${region} gp=true
  """
}

imputedChrs
  .groupTuple()
  .combine(reheaderIndex)
  .dump(tag: 'sampleImputed')
  .set{sampleImputedChrs}


/*--------------------------------------------------
  Concatenate sample-chromosome imputed calls
---------------------------------------------------*/

process mergeChromosomes {
  container 'vandhanak/bcftools:1.3.1'
  tag "${name}"
  publishDir 'results', mode: 'copy'

  input:
  set val(name), file('imputed_chrs_*.vcf.gz'), file(reheader) from sampleImputedChrs

  output:
  file("${name}.vcf.bgz") into sampleImputed
  file("${name}.vcf.tbi") into sampleImputedIndex

  script:
  """
  echo "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  $name" >> $reheader
  ls  *.vcf.gz | while read vcf; do tabix -fp vcf \$vcf; done
  bcftools concat -Oz -o ${name}.tmp.vcf.gz imputed_chrs_*.vcf.gz
  bcftools reheader -h $reheader -o ${name}.vcf.gz  ${name}.tmp.vcf.gz

  # sort the output VCF file
  gunzip ${name}.vcf.gz
  grep "^#" ${name}.vcf | uniq > output.vcf
  grep -v "^#" ${name}.vcf | sort -k1,1V -k2,2g >> output.vcf
  cat output.vcf > ${name}.vcf
  bgzip -c ${name}.vcf > ${name}.vcf.bgz
  tabix -p vcf ${name}.vcf.bgz.tbi
  """
}
