# estimate

ESTIMATE provides researchers with scores for tumor purity, the level of stromal cells present, and the infiltration level of immune cells in tumor tissues based on expression data

## Dependencies

* [estimate 1.0.13](http://R-Forge.R-project.org)
* [rstats 4.0](https://www.r-project.org/)


## Usage

### Cromwell
```
java -jar cromwell.jar run estimate.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputData`|Array[Pair[File,File]]+|Input files from RSEM and STAR.
`launchEstimate.estimateScript`|String|Script to run ESTIMATE
`launchEstimate.rsemZscoreRScript`|String|calculation of zScore for ESTIMATE results
`launchEstimate.ensFile`|String|path of a file for converting Ensembl gene_id to HUGO symbol


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|"ESTIMATE"|Output prefix, customizable. Default is the first file's basename.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`preProcessRsem.jobMemory`|Int|8|Memory allocated to the task.
`preProcessRsem.timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`preProcessRsem.tmpDir`|String|"tmp"|temporary directory
`preProcessRsem.dataDir`|String|"data"|data directory
`launchEstimate.jobMemory`|Int|8|Memory allocated to the task.
`launchEstimate.timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`launchEstimate.dataDir`|String|"."|data directory
`launchEstimate.modules`|String|"estimate/1.0.13"|Names and versions of required modules.


### Outputs

Output | Type | Description
---|---|---
`gRcounts`|File|File with RAW counts
`gCounts`|File|File with estimated counts
`gFpkm`|File|FPKMS from RSEM
`gTpm`|File|TPMS from RSEM
`estimateFile`|File|File with results from ESTIMATE


## Commands
 This section lists command(s) run by estimate workflow
 
 * Running ESTIMATE
 
 ### Merge data
 
 Bash code is used to extract data from RSEM and STAR inputs into
 separate tables for TPMs, FPKMs and counts.
 
 '''
 
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
 '''
 
 ### Run ESTIMATE using FPKM values
 
 '''
  set -euo pipefail
  Rscript ~{estimateScript} ~{inRSEM} ~{dataDir} ~{ensFile} ~{rsemZscoreRScript} ~{outputFileNamePrefix}
 '''
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
