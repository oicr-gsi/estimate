package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import org.apache.commons.io.FilenameUtils;
import org.mortbay.log.Log;

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
public class EstimateWorkflow extends OicrWorkflow {

    //dir
    private String dataDir, tmpDir;
    private String outDir;

    // Input Data
    private String inputRCOUNT;
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

    private void init() {
        try {
            //dir
            dataDir = "data/";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            inputRCOUNT = getProperty("input_rsem_file");
            inputRSEMFiles = getProperty("rsem_inputs");
            inputSTARFiles = getProperty("star_inputs");
            gmtFile = getWorkflowBaseDir() + "/dependencies/ensemble_conversion.txt";
            ensFile = getWorkflowBaseDir() + "/dependencies/dahaner2017_liu2015_immune_genesets.gmt";
            
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
        SqwFile file0 = this.createFile("RSEMinput");
        file0.setSourcePath(inputRCOUNT);
        file0.setType(TXT_METATYPE);
        file0.setIsInput(true);     
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {
        Job parentJob = null;
        this.outDir = this.outputFilenamePrefix + "_output/";
        String inRSEM = getFiles().get("RSEMinput").getProvisionedPath();
        
//        String[] tokens = inRSEM.split("\\.(?=[^\\.]+$)");
//        this.outputFilenamePrefix = tokens[0];
//        
        this.outputFilenamePrefix = FilenameUtils.getBaseName(inRSEM);

        String estimateGCT = this.dataDir + this.outputFilenamePrefix + ".txt.estimate.gct";
        String ssGSEA = this.dataDir + this.outputFilenamePrefix + ".txt.ssGSEA.txt";
        
        Job runEstimate = launchEstimate(inRSEM);
        parentJob = runEstimate;
        
        Job preProcess = postProcessRSEM(this.inputRSEMFiles, this.inputSTARFiles);
        preProcess.addParent(parentJob);
        parentJob = preProcess;
        
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
    
    private Job postProcessRSEM(String inRSEMs, String inSTARs) {
        Job postProcessRSEMGeneCounts = getWorkflow().createBashJob("post_process_RSEM");
        Command cmd = postProcessRSEMGeneCounts.getCommand();
        String[] rsems = inRSEMs.split(",");
        String[] stars = inSTARs.split(",");
        for (int i = 0; i < rsems.length; i++) {
            String rsembasename = FilenameUtils.getBaseName(rsems[i]).split(".")[0];
            String starbasename = FilenameUtils.getBaseName(stars[i]).split(".")[0];
            if (rsembasename != starbasename) {
                continue;
            }
            String gene = rsems[i];
            String rtab = stars[i];
            String geneCount = this.tmpDir + rsembasename + ".count";
            String geneRcount = this.tmpDir + rsembasename + ".rcount";
            cmd.addArgument("GENE=" + rsems[i] + ";");
            cmd.addArgument("RTAB=" + stars[i] + ";");
            cmd.addArgument("echo \"" + rsembasename + "\" > " + geneCount + "; cut -f5 " + gene + " | awk 'NR>1' >> " + geneCount + ";");
            cmd.addArgument("echo \"$EXTT\" > " + geneCount + "; awk 'NR>4 {if ($4 >= $3) print $4; else print $3}'" + rtab + " >> " + geneRcount + ";");

        }
        if (rsems.length != stars.length){
            Log.debug("PROCESSING only for "+ Integer.toString(rsems.length) + " RSEM files");
        }
        postProcessRSEMGeneCounts.setMaxMemory(Integer.toString(this.estimateMem * 1024));
        postProcessRSEMGeneCounts.setQueue(getOptionalProperty("queue", ""));
        return postProcessRSEMGeneCounts; 
    }
}