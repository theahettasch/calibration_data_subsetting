package eu.druglogics.synergy;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import java.util.ArrayList;
import java.util.List;

@Parameters(separators = "=")
class CommandLineArgs {

    @Parameter
    private List<String> parameters = new ArrayList<>();

    @Parameter(names = { "--project", "-p" },
            description = "Name of the project", order = 0)
    private String projectName;

    @Parameter(names = { "--inputDir", "-i" }, required = true,
            description = "Directory with all input files", order = 1)
    private String directoryInput;

    String getProjectName() {
        return projectName;
    }

    String getDirectoryInput() {
        return directoryInput;
    }

}
