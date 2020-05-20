#!/usr/bin/dev nextflow

params.ref = '/fs1/resources/ref/hg19/fasta/human_g1k_v37_decoy.fasta'
params.candidates_bcf = '/data/bnf/dev/sima/varloci/candidates.vcf'
params.purity = 0.75
//params.smpl_id ='KTC-HD829-C2'
smpl_id ='KTC-HD829-C2'
params.normal_bam='/fs1/results/myeloid/bam/KTC-HD829-C2.normal.dedup.bam'
params.tumor_bam ='/fs1/results/myeloid/bam/KTC-HD829-C2.tumor.dedup.bam'
params.outdir ='/data/bnf/dev/sima/varloci'
OUTDIR = params.outdir
params.vcfs_path = '/fs1/results/myeloid/vcf/'

process merge_vcfs{
	publishDir "$OUTDIR/vcf", mode :'copy'
	when:
	params.merge_vcfs
	
	output:
	file('CandidateVariants.vcf') into  candiate_vcf
	script:
	"""
	python /data/bnf/dev/sima/varloci/bin/merge_for_varlociraptor.py --callers 'vardict,tnscope,freebayes' --dir ${params.vcfs_path} --output CandidateVariants.vcf
	"""
	}

process  varloci_pryepro{
	publishDir "$OUTDIR/vcf", mode :'copy'	
	input:
		file(candidates) from candiate_vcf
		 //set val smpl_id,file bam_normal,file bam_tumor from bams_ch.view()
		//set val(group), smpl_id , val(type) from meta_manta.groupTuple()
	
	output:
		set val(smpl_id), file("${smpl_id}.tumor.observations.bcf"), file("${smpl_id}.normal.observations.bcf") into tum_normal_obs
		//set(smpl_id), ${sampl_id}.tumor.observations.bcf to normal_obs
		
	script:
		/*
		Tumor_index = type.findIndexOf{ it == 'tumor' }
		ID_Tumor = smpl_id[Tumor_index]
		tumor_index= id.findIndexOf{it == "$ID_Tumor" }
		bam_tumor = bam[tumor_index]
  		Normal_index = type.findIndexOf{ it == 'normal' }
		ID_normal = smpl_id[Normal_index]
		normal_index = id.findIndexOf{it == "$ID_normal" }
		bam_normal = bam[normal_index]
		*/
	//smpl_id = ${params.smpl_id}
	"""
	 varlociraptor preprocess variants ${params.ref} --bam ${params.normal_bam} --output ${smpl_id}.normal.observations.bcf < ${candidates}
	 varlociraptor preprocess variants ${params.ref} --bam ${params.tumor_bam} --output ${smpl_id}.tumor.observations.bcf < ${candidates}
	"""
}

process varloci_calling {
	publishDir "$OUTDIR/vcf", mode :'copy'
	input:
		set val(smpl_id), file(tumor_bcf), file(normal_bcf) from tum_normal_obs
	output:
		file "${smpl_id}.calls.bcf"  into calls_bcf
	script:
	
	"""
	varlociraptor call variants tumor-normal --purity ${params.purity} --tumor ${tumor_bcf} --normal ${normal_bcf} > ${smpl_id}.calls.bcf
	"""
}


process  FDR_filtering{

	 when:
		params.FDR
	 input:
		set val(smpl_id), file(calls_file) from calls_bcf
	 output:
		set val(smpl_id), file ("${smpl_id}.calls.filter.bcf") into callsfiltered
	 script:
	 """
		varlociraptor filter-calls control-fdr ${calls_file} --events SOMATIC_TUMOR --fdr 0.05 --var SNV > ${smpl_id}.calls.filter.bcf
	 """
	 }

process posteriorOdds_filtering{

	when:
		params.posterior
	input:
                set val(smpl_id), file(calls_file) from calls_bcf
	output:
                set val(smpl_id), file ("${smpl_id}.calls.posteriorFilter.bcf") into callsfiltered_post
        script:
	"""
	varlociraptor filter-calls posterior-odds --events SOMATIC_TUMOR --odds strong < ${calls_file} > ${smpl_id}.calls.posteriorFilter.bcf
	"""
	 }



