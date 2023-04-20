package eu.druglogics.synergy;

public enum InputFileIDs {

    TRAINING_DATA_FILE("training"),
    PERTURBATIONS_FILE("perturbations"),
    MODEL_OUTPUTS_FILE("modeloutputs"),
    CONFIGURATION_FILE("config"),
    DRUG_PANEL_FILE("drugpanel"),
    NETWORK_FILE("network");

    private String fileID;

    InputFileIDs(String fileID) {
        this.fileID = fileID;
    }

    public String getFileID() {
        return this.fileID;
    }

    public static String contains(String filename){
        for (InputFileIDs InputFileID : InputFileIDs.values()) {
            String fileID = InputFileID.getFileID();
            if (filename.contains(fileID)) {
                return fileID;
            }
        }
        return null;
    }
}
