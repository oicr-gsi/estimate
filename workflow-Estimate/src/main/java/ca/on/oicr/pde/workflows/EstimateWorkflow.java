package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import org.apache.commons.io.FilenameUtils;
import org.mortbay.log.Log;
import java.text.SimpleDateFormat;  

/**
 * <p>
 * For more information on developing workflows, see the documentation at
 * <a href="http://seqware.github.io/docs/6-pipeline/java-workflows/">SeqWare
 * Java Workflows</a>.</p>
 *
 * Quick reference for the order of methods called: 1. setupDirectory 2.
 * setupFiles 3. setupWorkflow 4. setupEnvironment 5. buildWorkflow
 *
 * See the SeqWare API for
 * <a href="http://seqware.github.io/javadoc/stable/apidocs/net/sourceforge/seqware/pipeline/workflowV2/AbstractWorkflowDataModel.html#setupDirectory%28%29">AbstractWorkflowDataModel</a>
 * for more information.
 */
/**
 *
 * @author prath@oicr.on.ca
 */
public class EstimateWorkflow extends OicrWorkflow {

    //dir
    private String dataDir, tmpDir;
    private String outDir;

    // Input Data
    private String outputFilenamePrefix;
    private String ensFile;
    private String gmtFile;
    
    private String inputRSEMFiles;
    private String inputSTARFiles;
     

    //Tools
    private String estimateScript;
    private String rpath;
    private String rScript;
    private String rLib;
    


    //Memory allocation
    private Integer estimateMem;


    private boolean manualOutput;
    private String queue;

    
    // metatypes
    private String TXT_METATYPE="txt/plain";
    
    //
    Date date = new Date();
    SimpleDateFormat formatter = new SimpleDateFormat("MM_dd_yyy");
    String currDateStr = formatter.format(date);

    private void init() {
        try {
            //dir
            dataDir = "data/";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            inputRSEMFiles = getProperty("rsem_inputs");
            inputSTARFiles = getProperty("star_inputs");
            gmtFile = getWorkflowBaseDir() + "/dependencies/ensemble_conversion.txt";
            ensFile = getWorkflowBaseDir() + "/dependencies/dahaner2017_liu2015_immune_genesets.gmt";
            
            outputFilenamePrefix = this.currDateStr + "_" + getProperty("study_title");
            
            //tools
            estimateScript = getWorkflowBaseDir() + "/dependencies/estimate.R";
            rpath = getProperty("rpath");
            rScript = getProperty("rpath") + "/bin/Rscript";
            rLib=getProperty("rLib");

            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            estimateMem = Integer.parseInt(getProperty("estimate_mem"));

        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void setupDirectory() {
        init();
        this.addDirectory(dataDir);
        this.addDirectory(tmpDir);
        if (!dataDir.endsWith("/")) {
            dataDir += "/";
        }
        if (!tmpDir.endsWith("/")) {
            tmpDir += "/";
        }
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        /**
         * Provisioning multiple RSEM files
         */
        Map<String,List<String>> inputFileMap = this.getRsemStarMap(this.inputRSEMFiles, this.inputSTARFiles);
        for (String key : inputFileMap.keySet()){
            SqwFile file0 = this.createFile("RSEM_"+key);
            SqwFile file1 = this.createFile("STAR_"+key);
            String rsemFile = inputFileMap.get(key).get(0);
            String starFile = inputFileMap.get(key).get(1);
            file0.setSourcePath(rsemFile);
            file0.setType(TXT_METATYPE);
            file0.setIsInput(true);
            file1.setSourcePath(starFile);
            file1.setType(TXT_METATYPE);
            file1.setIsInput(true);
        }  
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {
        Job parentJob = null;
        this.outDir = this.outputFilenamePrefix + "_output/";
        String postProcessedRSEM = this.dataDir + this.outputFilenamePrefix + "_genes_all_samples_RCOUNT.txt";
        
        String estimateGCT = this.dataDir + this.outputFilenamePrefix + ".txt.estimate.gct";
        String ssGSEA = this.dataDir + this.outputFilenamePrefix + ".txt.ssGSEA.txt";
        
        Job preProcess = postProcessRSEM(this.inputRSEMFiles, this.inputSTARFiles, postProcessedRSEM);
        parentJob = preProcess;
        
        Job runEstimate = launchEstimate(postProcessedRSEM);
        runEstimate.addParent(parentJob);
        parentJob = runEstimate;
        
        // Provision out estimae and ssgsea outputs
        SqwFile estimateOutput = createOutputFile(estimateGCT, TXT_METATYPE, this.manualOutput);
        estimateOutput.getAnnotations().put("estimate_output", "estimate");
        parentJob.addFile(estimateOutput);
        
        SqwFile ssGSEAOutput = createOutputFile(ssGSEA, TXT_METATYPE, this.manualOutput);
        ssGSEAOutput.getAnnotations().put("ssgsea_output", "ssgsea");
        parentJob.addFile(ssGSEAOutput);
//        
    }
    
    private Job launchEstimate(String inRSEM) {
        Job runEst = getWorkflow().createBashJob("run_estimate_ssgsea");
        Command cmd = runEst.getCommand();
        cmd.addArgument("export R_LIBS=" + rLib + ";");
        cmd.addArgument(this.rScript);
        cmd.addArgument(this.estimateScript);
        cmd.addArgument(inRSEM);
        cmd.addArgument(this.dataDir);
        cmd.addArgument(this.gmtFile);
        cmd.addArgument(this.ensFile);
        cmd.addArgument(getWorkflowBaseDir() + "/dependencies/convert_rsem_results_zscore.r");
        runEst.setMaxMemory(Integer.toString(this.estimateMem * 1024));
        runEst.setQueue(getOptionalProperty("queue", ""));
        return runEst;
    }   
    
    private Job postProcessRSEM(String inRSEMs, String inSTARs, String postProcessedRSEM) {
        Job postProcessRSEMGeneCounts = getWorkflow().createBashJob("post_process_RSEM");
        Command cmd = postProcessRSEMGeneCounts.getCommand();
        Map<String, List<String>> map = this.getRsemStarMap(inRSEMs, inSTARs);
        for (String key: map.keySet()){
//            Log.debug(key + " ... " + map.get(key));
            String geneCount = this.tmpDir + key + ".count";
            String geneRcount = this.tmpDir + key + ".rcount";
            String gene = getFiles().get("RSEM_"+key).getProvisionedPath();
            String rtab = getFiles().get("STAR_"+key).getProvisionedPath();
            cmd.addArgument("echo \"" + key + "\" > " + geneCount + ";");
            cmd.addArgument("cut -f5 " + gene + " | awk 'NR>1' >> " + geneCount + ";");
            cmd.addArgument("echo \"" + key + "\" > " 
                    + geneRcount 
                    + ";");
            cmd.addArgument("awk 'NR>4 {if ($4 >= $3) print $4; else print $3}' " 
                    + rtab + " >> " + geneRcount + ";");
        }
        cmd.addArgument("paste " + this.tmpDir + "*.rcount > " + postProcessedRSEM);
        postProcessRSEMGeneCounts.setMaxMemory(Integer.toString(this.estimateMem * 1024));
        postProcessRSEMGeneCounts.setQueue(getOptionalProperty("queue", ""));
        return postProcessRSEMGeneCounts; 
    }
    
    private String getSampleName(String fileBaseName, String extn){
        String[] sampleBaseName = fileBaseName.split(extn);
        List<String> sampleTokens = new ArrayList<String>(Arrays.asList(sampleBaseName[0].split("_")));
        List<String> sampleDesc = new ArrayList<String>();
        for (int i = 2; i < sampleTokens.size(); i++){
            String token = sampleTokens.get(i);
            sampleDesc.add(token);
        }
        String sampleName = String.join("_", sampleDesc);
        return sampleName;
    }
    
    private Map<String, List<String>> getRsemStarMap(String commaSeparatedRSEM, String commaSeparatedSTAR){
        /**
         * Given a list of comma Separated RSEM files;
         * get matching STAR
         */
        String[] rsemFilePaths = commaSeparatedRSEM.split(",");
        String[] starFilePaths = commaSeparatedSTAR.split(",");
        Map<String,List<String>> rsemStarMap = new HashMap<String, List<String>>();
        for (String rsemFile : rsemFilePaths){
            String rsemBaseName = FilenameUtils.getBaseName(rsemFile);
            String rsemSampleName = this.getSampleName(rsemBaseName, ".genes");
            List<String> vls = new ArrayList<String> ();
            vls.add(rsemFile);
            for (String starFile : starFilePaths){
                String starBaseName = FilenameUtils.getBaseName(starFile);
                String starSampleName = this.getSampleName(starBaseName, ".ReadsPerGene.out");
                if (starSampleName.equals(rsemSampleName)){
                    vls.add(starFile);
                }
            }
            if (vls.size() != 2 ){
                Log.debug("Matching files not present");
                continue;
            }
            rsemStarMap.put(rsemSampleName, vls);
        }
        return rsemStarMap;
    }
}