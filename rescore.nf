csv = file(params.csv)
mode = csv.countLines() > 2 ? "family" : "single"
trio = csv.countLines() > 3 ? true : false

process score_sv {
	cpus 5
	tag "$group $mode"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	memory '10 GB'
	time '2h'

	input:
		set group, file(vcf) from annotatedSV

	output:
		set group, file("${group}.sv.scored.sorted.vcf.gz"), file("${group}.sv.scored.sorted.vcf.gz.tbi") into sv_rescore
		file("${group}.INFO") into sv_INFO
		set group, file("${group}.sv.scored.sorted.vcf.gz") into svvcf_bed, svvcf_pod
				
	script:
	
		if (mode == "family") {
			"""
			genmod score -i $group -c $params.svrank_model -r $vcf -o ${group}.sv.scored_tmp.vcf
			bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored_tmp.vcf 
			bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
			tabix ${group}.sv.scored.sorted.vcf.gz -f
			echo "SV	${OUTDIR}/vcf/${group}.sv.scored.sorted.vcf.gz" > ${group}.INFO
			"""
		}
		else {
			"""
			genmod score -i $group -c $params.svrank_model_s -r $vcf -o ${group}.sv.scored.vcf
			bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored.vcf
			bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
			tabix ${group}.sv.scored.sorted.vcf.gz -f
			echo "SV	${OUTDIR}/vcf/${group}.sv.scored.sorted.vcf.gz" > ${group}.INFO
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
		set group, file(vcf), file(tbi) from sv_rescore
		file(ped) from ped_compound
		set group, file(snv), file(tbi) from snv_sv_vcf

	output:
		set group, file("${group}.snv.rescored.sorted.vcf.gz"), file("${group}.snv.rescored.sorted.vcf.gz.tbi"), \
			file("${group}.sv.rescored.sorted.vcf.gz"), file("${group}.sv.rescored.sorted.vcf.gz.tbi") into vcf_yaml
		file("${group}.INFO") into svcompound_INFO
				
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
		echo "SVc	${OUTDIR}/vcf/${group}.sv.rescored.sorted.vcf.gz,${OUTDIR}/vcf/${group}.snv.rescored.sorted.vcf.gz" > ${group}.INFO
		"""

}