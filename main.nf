$HOSTNAME = ""
params.outdir = 'results'  

// global param
params.mate = "pair"
params.mate2 = "single"

// Process Parameters:
// Process Parameters for Assemble_pairs_assemble_pairs:
params.Assemble_pairs_assemble_pairs.method = "align" 
params.Assemble_pairs_assemble_pairs.coord = "illumina" 
params.Assemble_pairs_assemble_pairs.rc = "tail"
params.Assemble_pairs_assemble_pairs.head_fields_R1 = "UMI TRIM"
params.Assemble_pairs_assemble_pairs.head_fields_R2 = ""
params.Assemble_pairs_assemble_pairs.failed = "false"
params.Assemble_pairs_assemble_pairs.fasta = "false"
params.Assemble_pairs_assemble_pairs.nproc = "${params.nproc}"
params.Assemble_pairs_assemble_pairs.alpha = 0.00001
params.Assemble_pairs_assemble_pairs.maxerror = 0.3
params.Assemble_pairs_assemble_pairs.minlen = 8
params.Assemble_pairs_assemble_pairs.maxlen = 1000
params.Assemble_pairs_assemble_pairs.scanrev = "false"
params.Assemble_pairs_assemble_pairs.minident = 0.5
params.Assemble_pairs_assemble_pairs.evalue = 0.00001
params.Assemble_pairs_assemble_pairs.maxhits = 100
params.Assemble_pairs_assemble_pairs.fill = "false"
params.Assemble_pairs_assemble_pairs.aligner = "blastn"
params.Assemble_pairs_assemble_pairs.// align_exec =  ""
params.Assemble_pairs_assemble_pairs.// dbexec =  ""   
params.Assemble_pairs_assemble_pairs.gap = 0
params.Assemble_pairs_assemble_pairs.usearch_version = "11.0.667"
params.Assemble_pairs_assemble_pairs.assemble_reference =  ''

// Process Parameters for Assemble_pairs_parse_log_AP:
params.Assemble_pairs_parse_log_AP.field_to_parse = "ID REFID LENGTH OVERLAP GAP ERROR IDENTITY PVALUE EVALUE1 EVALUE2"

// Process Parameters for Filter_Sequence_Quality_filter_seq_quality:
params.Filter_Sequence_Quality_filter_seq_quality.method = "quality"
params.Filter_Sequence_Quality_filter_seq_quality.nproc = "${params.nproc}"
params.Filter_Sequence_Quality_filter_seq_quality.q = "20"
params.Filter_Sequence_Quality_filter_seq_quality.n_length = "35"
params.Filter_Sequence_Quality_filter_seq_quality.n_missing = "10"

// Process Parameters for Parse_header_parse_headers:
params.Parse_header_parse_headers.method = "add"
params.Parse_header_parse_headers.act = "min"
params.Parse_header_parse_headers.args = "-f SAMPLE -u ${params.sample_name}"

// Process Parameters for collapse_sequences_collapse_seq:
params.collapse_sequences_collapse_seq.max_missing = 20
params.collapse_sequences_collapse_seq.inner = "true"
params.collapse_sequences_collapse_seq.fasta = "false"
params.collapse_sequences_collapse_seq.act = "set" 
params.collapse_sequences_collapse_seq.uf = "UMI" 
params.collapse_sequences_collapse_seq.cf = "SAMPLE"
params.collapse_sequences_collapse_seq.nproc = "${params.nproc}"
params.collapse_sequences_collapse_seq.failed = "false"

// Process Parameters for split_sequences_split_seq:
params.split_sequences_split_seq.field = "DUPCOUNT"
params.split_sequences_split_seq.num = 2
params.split_sequences_split_seq.fasta = "true"

if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.mate2){params.mate2 = ""} 

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_6_reads_g_11}
 } else {  
	g_6_reads_g_11 = Channel.empty()
 }

Channel.value(params.mate).into{g_7_mate_g_11;g_7_mate_g1_15;g_7_mate_g1_19;g_7_mate_g1_12}
Channel.value(params.mate2).into{g_10_mate_g2_5;g_10_mate_g2_0}


process unizp {

input:
 set val(name),file(reads) from g_6_reads_g_11
 val mate from g_7_mate_g_11

output:
 set val(name),file("*.fastq")  into g_11_reads0_g1_12

script:

readArray = reads.toString().split(' ')	
R1 = readArray.grep(~/.*R1.*/)[0]
R2 = readArray.grep(~/.*R2.*/)[0]

R1_zip = "readlink ${$R1}"
R2_zip = "readlink ${$R2}"

"""
case "$R1" in
*.gz | *.tgz ) 
        gunzip -k \$(${R1_zip}) 
        ;;
*)
        cp $R1 ./R2.fastq
        echo "$R1 not gzipped"
        ;;
esac

case "$R2" in
*.gz | *.tgz ) 
        gunzip -k \$(${R2_zip}) 
        ;;
*)
        cp $R2 ./R2.fastq
        echo "$R2 not gzipped"
        ;;
esac
"""
}


process Assemble_pairs_assemble_pairs {

input:
 set val(name),file(reads) from g_11_reads0_g1_12
 val mate from g_7_mate_g1_12

output:
 set val(name),file("*_assemble-pass.f*")  into g1_12_reads0_g2_0
 set val(name),file("AP_*")  into g1_12_logFile1_g1_15
 set val(name),file("*_assemble-fail.f*") optional true  into g1_12_reads_failed22
 set val(name),file("out*")  into g1_12_logFile33

script:
method = params.Assemble_pairs_assemble_pairs.method
coord = params.Assemble_pairs_assemble_pairs.coord
rc = params.Assemble_pairs_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_assemble_pairs.failed
fasta = params.Assemble_pairs_assemble_pairs.fasta
nproc = params.Assemble_pairs_assemble_pairs.nproc
alpha = params.Assemble_pairs_assemble_pairs.alpha
maxerror = params.Assemble_pairs_assemble_pairs.maxerror
minlen = params.Assemble_pairs_assemble_pairs.minlen
maxlen = params.Assemble_pairs_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_assemble_pairs.scanrev
minident = params.Assemble_pairs_assemble_pairs.minident
evalue = params.Assemble_pairs_assemble_pairs.evalue
maxhits = params.Assemble_pairs_assemble_pairs.maxhits
fill = params.Assemble_pairs_assemble_pairs.fill
aligner = params.Assemble_pairs_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_assemble_pairs.// dbexec
gap = params.Assemble_pairs_assemble_pairs.gap
usearch_version = params.Assemble_pairs_assemble_pairs.usearch_version
assemble_reference = params.Assemble_pairs_assemble_pairs.assemble_reference
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi
	
	echo \$align_exec
	echo \$dbexec
	
	
	
	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc} >> out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process Assemble_pairs_parse_log_AP {

input:
 set val(name),file(log_file) from g1_12_logFile1_g1_15
 val mate from g_7_mate_g1_15

output:
 set val(name),file("*.tab")  into g1_15_logFile0_g1_19

script:
field_to_parse = params.Assemble_pairs_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Assemble_pairs_report_assemble_pairs {

input:
 set val(name),file(log_files) from g1_15_logFile0_g1_19
 val matee from g_7_mate_g1_19

output:
 file "*.rmd"  into g1_19_rMarkdown0_g1_22



shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')
	assemble = readArray[0]
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("assemble_length", "Histogram showing the distribution assembled sequence lengths in 
	                            nucleotides for the Align step (top) and Reference step (bottom).")
	figures("assemble_overlap", "Histogram showing the distribution of overlapping nucleotides between 
	                             mate-pairs for the Align step (top) and Reference step (bottom).
	                             Negative values for overlap indicate non-overlapping mate-pairs
	                             with the negative value being the number of gap characters between
	                             the ends of the two mate-pairs.")
	figures("assemble_error", "Histograms showing the distribution of paired-end assembly error 
	                           rates for the Align step (top) and identity to the reference germline 
	                           for the Reference step (bottom).")
	figures("assemble_pvalue", "Histograms showing the distribution of significance scores for 
	                            paired-end assemblies. P-values for the Align mode are shown in the top
	                            panel. E-values from the Reference step's alignment against the 
	                            germline sequences are shown in the bottom panel for both input files
	                            separately.")
	```
	
	```{r, echo=FALSE, warning=FALSE}
	assemble_log <- loadLogTable(file.path(".", "!{assemble}"))
	
	# Subset to align and reference logs
	align_fields <- c("ERROR", "PVALUE")
	ref_fields <- c("REFID", "GAP", "EVALUE1", "EVALUE2", "IDENTITY")
	align_log <- assemble_log[!is.na(assemble_log$ERROR), !(names(assemble_log) %in% ref_fields)]
	ref_log <- assemble_log[!is.na(assemble_log$REFID), !(names(assemble_log) %in% align_fields)]
	
	# Build log set
	assemble_list <- list()
	if (nrow(align_log) > 0) { assemble_list[["Align"]] <- align_log }
	if (nrow(ref_log) > 0) { assemble_list[["Reference"]] <- ref_log }
	plot_titles <- names(assemble_list)
	```
	
	# Paired-End Assembly
	
	Assembly of paired-end reads is performed using the AssemblePairs tool which 
	determines the read overlap in two steps. First, de novo assembly is attempted 
	using an exhaustive approach to identify all possible overlaps between the 
	two reads with alignment error rates and p-values below user-defined thresholds. 
	This method is denoted as the `Align` method in the following figures. 
	Second, those reads failing the first stage of de novo assembly are then 
	mapped to the V-region reference sequences to create a full length sequence, 
	padding with Ns, for any amplicons that have insufficient overlap for 
	de novo assembly. This second stage is referred to as the `Reference` step in the
	figures below.
	
	## Assembled sequence lengths
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="length", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_length")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="overlap", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_overlap")`
	
	## Alignment error rates and significance
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="error", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="pvalue", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```

	`r figures("assemble_pvalue")`

	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}
}


process Assemble_pairs_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "assemble_pairs_report/$filename"}
input:
 file rmk from g1_19_rMarkdown0_g1_22

output:
 file "*.html"  into g1_22_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Filter_Sequence_Quality_filter_seq_quality {

input:
 set val(name),file(reads) from g1_12_reads0_g2_0
 val mate from g_10_mate_g2_0

output:
 set val(name), file("*_${method}-pass.fastq")  into g2_0_reads0_g3_15
 file "FS_*"  into g2_0_logFile1_g2_5
 set val(name), file("*_${method}-fail.fastq") optional true  into g2_0_reads22

script:
method = params.Filter_Sequence_Quality_filter_seq_quality.method
nproc = params.Filter_Sequence_Quality_filter_seq_quality.nproc
q = params.Filter_Sequence_Quality_filter_seq_quality.q
n_length = params.Filter_Sequence_Quality_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Quality_filter_seq_quality.n_missing
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="quality"){
	q = "-q ${q}"
	n_length = ""
	n_missing = ""
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = ""
		n_length = ""
		n_missing = "-n ${n_missing}"
	}
}

readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --outname ${name}_R1 --nproc ${nproc} --log FS_R1_${name}.log --failed
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --outname ${name}_R2 --nproc ${nproc} --log FS_R2_${name}.log --failed
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --outname ${name} --nproc ${nproc} --log FS_${name}.log --failed
	"""
}


}


process Parse_header_parse_headers {

input:
 set val(name), file(reads) from g2_0_reads0_g3_15

output:
 set val(name),file("*${out}")  into g3_15_reads0_g4_16

script:
method = params.Parse_header_parse_headers.method
act = params.Parse_header_parse_headers.act
args = params.Parse_header_parse_headers.args


if(method=="collapse" || method=="copy" || method=="rename" || method=="merge"){
	out="_reheader.fastq"
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} --act ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} ${args}
			"""	
	}else{
		out="_reheader.fastq"
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}


}


process collapse_sequences_collapse_seq {

input:
 set val(name), file(reads) from g3_15_reads0_g4_16

output:
 set val(name),  file("*_collapse-unique.fast*")  into g4_16_reads0_g5_20
 set val(name),  file("*_collapse-duplicate.fast*") optional true  into g4_16_reads_duplicate11
 set val(name),  file("*_collapse-undetermined.fast*") optional true  into g4_16_reads_undetermined22
 file "CS_*"  into g4_16_logFile33

script:
max_missing = params.collapse_sequences_collapse_seq.max_missing
inner = params.collapse_sequences_collapse_seq.inner
fasta = params.collapse_sequences_collapse_seq.fasta
act = params.collapse_sequences_collapse_seq.act
uf = params.collapse_sequences_collapse_seq.uf
cf = params.collapse_sequences_collapse_seq.cf
nproc = params.collapse_sequences_collapse_seq.nproc
failed = params.collapse_sequences_collapse_seq.failed

inner = (inner=="true") ? "--inner" : ""
fasta = (fasta=="true") ? "--fasta" : ""
act = (act=="none") ? "" : "--act ${act}"
cf = (cf=="") ? "" : "--cf ${cf}"
uf = (uf=="") ? "" : "--uf ${uf}"
failed = (failed=="false") ? "" : "--failed"

"""
CollapseSeq.py -s ${reads} -n ${max_missing} ${fasta} ${inner} ${uf} ${cf} ${act} --log CS_${name}.log ${failed}
"""

}


process split_sequences_split_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_atleast-.*.fast.*$/) "reads/$filename"}
input:
 set val(name),file(reads) from g4_16_reads0_g5_20

output:
 set val(name), file("*_atleast-*.fast*")  into g5_20_reads00

script:
field = params.split_sequences_split_seq.field
num = params.split_sequences_split_seq.num
fasta = params.split_sequences_split_seq.fasta

readArray = reads.toString()

if(num!=0){
	num = " --num ${num}"
}else{
	num = ""
}

fasta = (fasta=="false") ? "" : "--fasta"

"""
SplitSeq.py group -s ${readArray} -f ${field} ${num} ${fasta}
"""

}


process Filter_Sequence_Quality_parse_log_FS {

input:
 set val(name), file(log_file) from g2_0_logFile1_g2_5
 val mate from g_10_mate_g2_5

output:
 set val(name), file("*.tab")  into g2_5_logFile00

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
