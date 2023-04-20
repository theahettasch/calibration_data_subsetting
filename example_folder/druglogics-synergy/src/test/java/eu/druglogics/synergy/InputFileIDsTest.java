package eu.druglogics.synergy;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertEquals;

class InputFileIDsTest {

    @Test
    void test_contains() {
        assertNull(InputFileIDs.contains("a_file_with_no_associated_file_ids"));
        assertNull(InputFileIDs.contains("topology"));

        assertEquals(InputFileIDs.contains("toy_ags_training_data.tab"), "training");
        assertEquals(InputFileIDs.contains("toy_ags_perturbations.tab"), "perturbations");
        assertEquals(InputFileIDs.contains("toy_ags_modeloutputs.tab"), "modeloutputs");
        assertEquals(InputFileIDs.contains("toy_ags_config.tab"), "config");
        assertEquals(InputFileIDs.contains("toy_ags_drugpanel.tab"), "drugpanel");
        assertEquals(InputFileIDs.contains("/home/user/toy_ags_network.sif"), "network");

        assertEquals("training", InputFileIDs.TRAINING_DATA_FILE.getFileID());
        assertEquals("perturbations", InputFileIDs.PERTURBATIONS_FILE.getFileID());
        assertEquals("modeloutputs", InputFileIDs.MODEL_OUTPUTS_FILE.getFileID());
        assertEquals("config", InputFileIDs.CONFIGURATION_FILE.getFileID());
        assertEquals("drugpanel", InputFileIDs.DRUG_PANEL_FILE.getFileID());
        assertEquals("network", InputFileIDs.NETWORK_FILE.getFileID());
    }
}
