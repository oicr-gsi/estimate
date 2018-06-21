package ca.on.oicr.pde.deciders;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;

/**
 *
 * @author prath@oicr.on.ca
 */
public class EstimateDecider extends OicrDecider {

    private final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;

    private String templateType = "WT";
    private String queue = "";
    private String externalID;

    private final static String TXT_METATYPE = "text/plain";
//    private String tumorType;
//    private List<String> results;

    public EstimateDecider() {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("template-type", "Required. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type").withRequiredArg();
        parser.accepts("queue", "Optional: Set the queue (Default: not set)").withRequiredArg();
        parser.accepts("tumor-type", "Optional: Set tumor tissue type to something other than primary tumor (P), i.e. X . Default: Not set (All)").withRequiredArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(Arrays.asList(TXT_METATYPE));
        this.setHeadersToGroupBy(Arrays.asList(Header.FILE_SWA));

        ReturnValue rv = super.init();
        rv.setExitStatus(ReturnValue.SUCCESS);

        //Group by sample if no other grouping selected
        if (this.options.has("group-by")) {
            Log.error("group-by parameter passed, but this decider does not allow overriding the default grouping (by Donor + Library Type)");
        }

        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        }

        if (this.options.has("template-type")) {
            if (!options.hasArgument("template-type")) {
                Log.error("--template-type requires an argument, WT");
                rv.setExitStatus(ReturnValue.INVALIDARGUMENT);
                return rv;
            } else {
                this.templateType = options.valueOf("template-type").toString();
                if (!this.templateType.equals("WT")) {
                    Log.stderr("NOTE THAT ONLY WT template-type SUPPORTED, WE CANNOT GUARANTEE MEANINGFUL RESULTS WITH OTHER TEMPLATE TYPES");
                    rv.setExitStatus(ReturnValue.INVALIDARGUMENT);
                }
            }
        }

        return rv;
    }

    /**
     * Final check
     *
     * @param commaSeparatedFilePaths
     * @param commaSeparatedParentAccessions
     *
     * @return
     */
    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String[] filePaths = commaSeparatedFilePaths.split(",");
        boolean haveFiles = false;

        for (String p : filePaths) {
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p)) {
                    continue;
                }
                String tt = bs.getTissueType();

                if (!tt.isEmpty()) {
                    haveFiles = true;
                }
            }
        }
        if (haveFiles) {
            return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        }
       Log.error("Data not available, WON'T RUN");
        return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);
        String currentTtype = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        String currentTissueType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_tissue_type");

        if (null == currentTissueType) {
            return false; // we need only those which have their tissue type set
        }

        // Filter the data of a different template type if filter is specified
        if (!this.templateType.equalsIgnoreCase(currentTtype)) {
            Log.warn("Excluding file with SWID = [" + returnValue.getAttribute(Header.FILE_SWA.getTitle())
                    + "] due to template type/geo_library_source_template_type = [" + currentTtype + "]");
            return false;
        }


        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        Log.debug("Number of files from file provenance = " + vals.size());

        // get files from study
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        // Override the supplied group-by value
        for (ReturnValue currentRV : vals) {
            boolean metatypeOK = false;
            boolean fileExtnOK = false;

            for (int f = 0; f < currentRV.getFiles().size(); f++) {
                try {
                    if (currentRV.getFiles().get(f).getMetaType().equals(TXT_METATYPE)) {
                        metatypeOK = true;
                    }
                    if (currentRV.getFiles().get(f).getFilePath().endsWith(".ReadsPerGene.out.tab") 
                            || currentRV.getFiles().get(f).getFilePath().endsWith(".genes.results")){
                        fileExtnOK = true;
                    } else {
                        Log.debug("Undesired file type "+currentRV.getFiles(). get(f).getFilePath());
                    }
                } catch (Exception e) {
                    Log.stderr("Error checking a file");
                }
            }
            
//            Log.debug(fileExtnOK);
           

            if (!metatypeOK && !fileExtnOK) {
                continue; // Go to the next value
            }

            BeSmall currentSmall = new BeSmall(currentRV);
            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);
            //make sure you only have the most recent single file for each
            //sequencer run + lane + barcode + workflow-name
            String fileDeets = currentSmall.getIusDetails() + "_" + currentSmall.getWorkflowDetails();
            Date currentDate = currentSmall.getDate();

            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                iusDeetsToRV.put(fileDeets, currentRV);
            } //if there is an entry, compare the current value to the 'old' one in
            //the map. if the current date is newer than the 'old' date, replace
            //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);
                BeSmall oldSmall = fileSwaToSmall.get(oldRV.getAttribute(Header.FILE_SWA.getTitle()));
                Date oldDate = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    iusDeetsToRV.put(fileDeets, currentRV);
                }
            }
            
        }
        
//        for (ReturnValue rv: iusDeetsToRV.values()){
//            BeSmall currSmall = new BeSmall(rv);
//            String sampleDeets = currSmall.getSampleNameDetails();
//            String fileSuffix = currSmall.getPath();
//            String workflowName = currSmall.getWorkflowDetails();
//            if (workflowName == "RSEM"){
//                
//            }
//            
//        }

        //only use those files that entered into the iusDeetsToRV
        //since it's a map, only the most recent values
        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();

        //group files according to the designated header (e.g. sample SWID)
        for (ReturnValue r : newValues) {
            String currVal = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).getGroupByAttribute();
            Log.debug(currVal);
            Log.debug(r.getFiles().get(0).getFilePath());
            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }

        return map;
    }

    @Override
    protected String handleGroupByAttribute(String attribute) {
        String a = super.handleGroupByAttribute(attribute);
        BeSmall small = fileSwaToSmall.get(a);
        if (small != null) {
            return small.getGroupByAttribute();
        }
        return attribute;
    }

    

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {

        String[] filePaths = commaSeparatedFilePaths.split(",");
        List<String> rsemGeneCounts = new ArrayList<String>();
        List<String> starGeneCounts = new ArrayList<String>();
//        Map<String,Integer> rsemSTARMap = new HashMap<String,Integer>();

        for (String p : filePaths) {

               for (BeSmall bs : fileSwaToSmall.values()) {
                String tt = bs.getTissueType();
//                String sampleName = bs.getSampleNameDetails();
                if (!tt.isEmpty()) {
//                        Log.stdout("WRITING TO INI FILE ... " + bs.getPath());
                    if (!bs.getPath().equals(p)) {
                        continue;
                    }
                    if (p.endsWith(".genes.results")) {
                        rsemGeneCounts.add(p);
//                        rsemSTARMap.put(sampleName, );
                    } else if (p.endsWith(".ReadsPerGene.out.tab")){
                        
                        starGeneCounts.add(p);
                    }

                } else {
                    Log.error("THE DONOR does not have data to run the workflow");
                    abortSchedulingOfCurrentWorkflowRun();
                }
            }

        }
        Map<String, String> iniFileMap = super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);

        iniFileMap.put("data_dir", "data");
        iniFileMap.put("rsem_inputs", String.join(",", rsemGeneCounts));
        iniFileMap.put("star_inputs", String.join(",", starGeneCounts));
        iniFileMap.put("template_type", this.templateType);
//        iniFileMap.put("external_name", this.externalID);
        if (!this.queue.isEmpty()) {
            iniFileMap.put("queue", this.queue);
        }

        return iniFileMap;

//        return super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);
    }

    public static void main(String args[]) {

        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(EstimateDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }
    private class BeSmall {

        private Date date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String tissueType = null;
        private String path = null;
        private String extName = null;
        private String groupID = null;
        private String groupDescription = null;
        private String workflowDetails = null;
        private String sampleNameDetails = null;

        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            workflowDetails = rv.getAttribute(Header.WORKFLOW_NAME.getTitle());
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode() + "_"+ workflowDetails;
            sampleNameDetails =iusDetails = fa.getLibrarySample() + "_" + fa.getSequencerRun() + fa.getLane() + fa.getBarcode();
            tissueType = fa.getLimsValue(Lims.TISSUE_TYPE);
            extName = rv.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_external_name");
            //fa.getLimsValue(Lims.TUBE_ID);
            if (null == extName || extName.isEmpty()) {
                extName = "NA";
            }
            groupID = fa.getLimsValue(Lims.GROUP_ID);
            if (null == groupID || groupID.isEmpty()) {
                groupID = "NA";
            }
            groupDescription = fa.getLimsValue(Lims.GROUP_DESC);
            if (null == groupDescription || groupDescription.isEmpty()) {
                groupDescription = "NA";
            }
            groupByAttribute = fa.getStudy() + ":" + fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE) + fa.getMetatype();
            path = rv.getFiles().get(0).getFilePath() + "";
            
        }

        public String getSampleNameDetails() {
            return sampleNameDetails;
        }

        public String getWorkflowDetails() {
            return workflowDetails;
        }

        public Date getDate() {
            return date;
        }

        public void setDate(Date date) {
            this.date = date;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }

        public String getTissueType() {
            return tissueType;
        }

        public String getIusDetails() {
            return this.iusDetails;
        }

        public void setIusDetails(String iusDetails) {
            this.iusDetails = iusDetails;
        }

        public String getPath() {
            return path;
        }

        public String getExtName() {
            return extName;
        }

        public String getGroupID() {
            return groupID;
        }

        public String getGroupDescription() {
            return groupDescription;
        }

        public void setPath(String path) {
            this.path = path;
        }
    }


}
