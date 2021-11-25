version 1.0

# ================================================================================
# Workflow 
# ================================================================================
workflow estimate {
input {
 Array[Pair[File, File]]+ inputData
 String outputFileNamePrefix = "ESTIMATE"
}

scatter(i in inputData) {
  File rsemFiles = i.left
  File starFiles = i.right
}

call preProcessRsem { input: rsemData = rsemFiles, starData = starFiles, outputPrefix = outputFileNamePrefix }
call launchEstimate { input: inRSEM = preProcessRsem.gFpkm, outputFileNamePrefix = outputFileNamePrefix }

parameter_meta {
  parameter1: "Input file with the first mate reads."
  parameter2: " Input file with the second mate reads."
  outputFileNamePrefix: "Output prefix, customizable. Default is the first file's basename."
}

meta {
    author: "Peter Ruzanov"
    email: "peter.ruzanov@oicr.on.ca"
    description: "ESTIMATE provides researchers with scores for tumor purity, the level of stromal cells present, and the infiltration level of immune cells in tumor tissues based on expression data"
    dependencies: [
      {
        name: "estimate/1.0.13",
        url:  "http://R-Forge.R-project.org"
      },
      {
        name: "rstats/4.0",
        url: "https://www.r-project.org/"
      }
    ]
    output_meta: {
       gRcounts: "File with RAW counts",
       gCounts: "File with estimated counts",
       gFpkm: "FPKMS from RSEM",
       gTpm: "TPMS from RSEM",
       estimateFile: "File with results from ESTIMATE"
    }
}

output {
  File gRcounts = preProcessRsem.gRcounts
  File gCounts  = preProcessRsem.gCounts
  File gFpkm = preProcessRsem.gFpkm
  File gTpm = preProcessRsem.gTpm
  File estimateFile = launchEstimate.estimateFile
}
}


# ====================
#    PRE-PROCESS RSEM
# ====================
task preProcessRsem {
input {
  Int  jobMemory = 8
  Int  timeout   = 20
  String outputPrefix
  Array[File]+ rsemData
  Array[File]+ starData
  String tmpDir = "tmp"
  String dataDir = "data"
}

command <<<

TMP='~{tmpDir}'
DATA='~{dataDir}'
mkdir $TMP
mkdir $DATA

cp ~{sep=' ' rsemData} $DATA/
cp ~{sep=' ' starData} $DATA/

STARG=$(ls $DATA/*.tab | head -1);
if [ ! -z $STARG ]; then
  awk 'NR>3 {print $1}' $STARG | sed 's/N_ambiguous/gene_id/' > $TMP/sgene;
fi;

RSEMG=$(ls $DATA/*.genes.results | head -1);
if [ ! -z $RSEMG ]; then
  cut -f 1 $RSEMG  > $TMP/genes;
fi;

# We will use basename as a sample ID here
for t in $DATA/*results;do
  BASE=$(basename $t | sed s/.results$//);
  NAME=$(echo $BASE | sed 's/\..*//');
  echo $t;
  echo $NAME > $TMP/$NAME.fpkm;
  cut -f 7 $t | awk 'NR>1' >> $TMP/$NAME.fpkm;
  echo $NAME > $TMP/$NAME.tpm;
  cut -f 6 $t | awk 'NR>1' >> $TMP/$NAME.tpm;
  echo $NAME > $TMP/$NAME.count;
  cut -f 5 $t | awk 'NR>1' >> $TMP/$NAME.count;
  echo $NAME > $TMP/$NAME.rcount;
  awk 'NR>4 {if ($4 >= $3) print $4; else print $3}' $DATA/$NAME.ReadsPerGene.out.tab >> $TMP/$NAME.rcount;
done

# Merging
paste $TMP/sgene $TMP/*.rcount > ~{outputPrefix}_genes_all_samples_RCOUNT.txt;
paste $TMP/genes $TMP/*.count > ~{outputPrefix}_genes_all_samples_COUNT.txt;
paste $TMP/genes $TMP/*.fpkm > ~{outputPrefix}_genes_all_samples_FPKM.txt;
paste $TMP/genes $TMP/*.tpm > ~{outputPrefix}_genes_all_samples_TPM.txt;
>>>

parameter_meta {
 jobMemory: "Memory allocated to the task."
 timeout: "Timeout in hours, needed to override imposed limits."
 rsemData: "RSEM files ordered the same way as STAR files, passed as array"
 starData: "STAR files ordered the same way as RSEM files, passed as array"
 tmpDir: "temporary directory"
 dataDir: "data directory"
}

runtime {
  memory:  "~{jobMemory} GB"
  timeout: "~{timeout}"
}

output {
  File gRcounts = "~{outputPrefix}_genes_all_samples_RCOUNT.txt"
  File gCounts  = "~{outputPrefix}_genes_all_samples_COUNT.txt"
  File gFpkm = "~{outputPrefix}_genes_all_samples_FPKM.txt"
  File gTpm = "~{outputPrefix}_genes_all_samples_TPM.txt"
}
}

# =========================
#   LAUNCH ESTIMATE
# =========================
task launchEstimate {
input {
  Int  jobMemory = 8
  Int  timeout   = 20
  String estimateScript
  String rsemZscoreRScript
  File inRSEM
  String ensFile
  # String gmtFile
  String dataDir = "."
  String modules = "estimate/1.0.13"
  String outputFileNamePrefix = "ESTIMATE"
}

command <<<
 set -euo pipefail
 Rscript ~{estimateScript} ~{inRSEM} ~{dataDir} ~{ensFile} ~{rsemZscoreRScript} ~{outputFileNamePrefix}
>>>

parameter_meta {
 jobMemory: "Memory allocated to the task."
 modules: "Names and versions of required modules. This needs to be customized by shesmu"
 timeout: "Timeout in hours, needed to override imposed limits."
 estimateScript: "Script to run ESTIMATE"
 rsemZscoreRScript: "calculation of zScore for ESTIMATE results"
 inRSEM: "RSEM inputs, pre-processed"
 ensFile: "file for converting Ensembl gene_id to HUGO symbol"
 # gmtFile: "GMT file"
 dataDir: "data directory"
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File estimateFile = "~{outputFileNamePrefix}.estimate.gct"
}
}

