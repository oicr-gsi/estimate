package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.testing.workflow.DryRun;
import ca.on.oicr.pde.testing.workflow.TestDefinition;
import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Map.Entry;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import org.apache.commons.io.FileUtils;
import static org.testng.Assert.assertTrue;
import org.testng.annotations.Test;

/**
 *
 * @author prath
 */
public class EstimateWorkflowTest {

    public EstimateWorkflowTest() {
    }

    @Test
    public void validateRegressionTestDefinition() throws IllegalAccessException, InstantiationException, IOException, Exception {
        TestDefinition td = TestDefinition.buildFromJson(FileUtils.readFileToString(new File("src/test/resources/tests.json")));
        for (TestDefinition.Test t : td.getTests()) {
            DryRun d = new DryRun(System.getProperty("bundleDirectory"), t.getParameters(), EstimateWorkflow.class);
            AbstractWorkflowDataModel workflow = d.buildWorkflowModel();
            d.validateWorkflow();

            //validate input files
            Map<String, SqwFile> files = workflow.getFiles();
            for (Entry<String, SqwFile> f : files.entrySet()) {
                String sqwfileId = f.getKey();
                String provisionInFilePath = f.getValue().getProvisionedPath();
                assertTrue(provisionInFilePath.endsWith(".out.tab") || provisionInFilePath.endsWith(".genes.results"),
                        sqwfileId + " has invalid provision in file path: [" + provisionInFilePath + "]");
            }
        }
    }

}
