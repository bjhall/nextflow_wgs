
// Input channels for various meta information //
Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, file(row.sv_vcf), file(row.snv_vcf) ) }
	.set{ input }

Channel
    .fromPath(file(params.ped))
    .into{ strip_ped; ped_compound}
OUTDIR = params.outdir+'/'+params.subdir

ped = file(params.ped)
mode = ped.countLines() > 1 ? "family" : "single"
println(mode)

process strip_score {
    cpus 1
    time '30m'
    memory '5 GB'

    input:
        set group, file(sv_vcf), file(snv_vcf) from input
        file(ped) from strip_ped
    
    output:
        set group, file("${group}.sv_stripped.vcf"), file(snv_vcf) into stripped_vcf

    """
    rescore_vcf.pl --vcf $sv_vcf --ped $ped > ${group}.sv_stripped.vcf
    """
}

process score_sv {
	cpus 5
	tag "$group $mode"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	memory '10 GB'
	time '2h'

	input:
		set group, file(vcf), file(snv) from stripped_vcf

	output:
		set group, file("${group}.sv.scored.sorted.vcf.gz"), file("${group}.sv.scored.sorted.vcf.gz.tbi"), file(snv) into sv_rescore
				
	script:
	
		if (mode == "family") {
			"""
			genmod score -i $group -c $params.svrank_model -r $vcf -o ${group}.sv.scored_tmp.vcf
			bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored_tmp.vcf 
			bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
			tabix ${group}.sv.scored.sorted.vcf.gz -f
			"""
		}
		else {
			"""
			genmod score -i $group -c $params.svrank_model_s -r $vcf -o ${group}.sv.scored.vcf
			bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored.vcf
			bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
			tabix ${group}.sv.scored.sorted.vcf.gz -f
			"""
		}
}

process compound_finder {
	cpus 5
	tag "$group $mode"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	memory '10 GB'
	time '2h'

	when:
		mode == "family"

	input:
		set group, file(vcf), file(tbi), file(snv) from sv_rescore
		file(ped) from ped_compound

	output:
		set group, file("${group}.snv.rescored.sorted.vcf.gz"), file("${group}.snv.rescored.sorted.vcf.gz.tbi"), \
			file("${group}.sv.rescored.sorted.vcf.gz"), file("${group}.sv.rescored.sorted.vcf.gz.tbi") into vcf_yaml
				
	script:
		"""
		compound_finder.pl \\
			--sv $vcf --ped $ped --snv $snv \\
			--osv ${group}.sv.rescored.sorted.vcf \\
			--osnv ${group}.snv.rescored.sorted.vcf 
		bgzip -@ ${task.cpus} ${group}.sv.rescored.sorted.vcf -f
		bgzip -@ ${task.cpus} ${group}.snv.rescored.sorted.vcf -f
		tabix ${group}.sv.rescored.sorted.vcf.gz -f
		tabix ${group}.snv.rescored.sorted.vcf.gz -f
		"""

}