echo Set up!

SEEDFILE=data.cell_level_qc.txt

bam=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
nam=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2}')
bc=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $3}')

bash /stanley/levin_dr/kwanho/git_repos/CellLevel_QC/scripts/PrepSTARSolo.sh \
	${bam} \
	/stanley/levin_xinjin/kwanho/BD_Rhapsody/ref/pre-built_BD_ref/BD_Rhapsody_Reference_Files/gencode.vM31.primary_assembly.annotation-phix_genome-processed.gtf \
	${nam}
	
java -jar /stanley/levin_dr/kwanho/git_repos/CellLevel_QC/Jar/SingleCellQC.jar \
	-i ${nam}.label.bam \
	-c ${bc} \
	-z \
	-o qc_${nam}.txt \
	-q STARSolo

echo DONE

