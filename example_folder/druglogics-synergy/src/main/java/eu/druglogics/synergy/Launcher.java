package eu.druglogics.synergy;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import eu.druglogics.drabme.Drabme;
import eu.druglogics.gitsbe.Gitsbe;

import javax.naming.ConfigurationException;
import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import static eu.druglogics.gitsbe.util.Util.abort;
import static eu.druglogics.gitsbe.util.Util.getFileName;

public class Launcher {

	// Global input
	private String projectName;
	private String directoryOutput;

	// Gitsbe input
	private String filenameNetwork = "";
	private String filenameTrainingData = "";
	private String filenameModelOutputs = "";
	private String filenameConfig = "";
	private String directoryTmpGitsbe;

	// Drabme input
	private String filenameDrugs = "";
	private String filenamePerturbations = "";
	private String directoryModels;
	private String directoryTmpDrabme;

	protected Launcher() {}

	public static void main(String[] args) {
		Launcher drugLogicsLauncher = new Launcher();
		drugLogicsLauncher.start(args);
	}

	private void start(String[] args) {
		setupAndValidateInput(args);

		runGitsbe();
		runDrabme();
	}

	private void setupAndValidateInput(String[] args) {
		try {
			CommandLineArgs arguments = new CommandLineArgs();
			JCommander.newBuilder().addObject(arguments).build().parse(args);

			projectName = arguments.getProjectName();
			String directoryInput = arguments.getDirectoryInput();

			// projectName is not required, but we have to set it
			if (projectName == null) {
				projectName = getFileName(directoryInput);
			}

			DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmmss");

			directoryOutput = new File(directoryInput, projectName + "_" +
				dateFormat.format(Calendar.getInstance().getTime())).getAbsolutePath();
			directoryTmpGitsbe = new File(directoryOutput, "gitsbe_tmp").getAbsolutePath();
			directoryTmpDrabme = new File(directoryOutput, "drabme_tmp").getAbsolutePath();
			directoryModels = new File(directoryOutput, "models").getAbsolutePath();

			System.out.println("Finding input files in directory: " + directoryInput);
			File[] files = new File(directoryInput).listFiles();

			for (File file : files) {
				if (!file.isDirectory()) {
					String filename = file.getName();
					loadFileFromDirectory(filename, directoryInput);
				}
			}

			printInputFilesMessage();

		} catch (ParameterException parEx) {
			System.out.println("\nOptions preceded by an asterisk are required.");
			parEx.getJCommander().setProgramName("eu.druglogics.synergy.Launcher");
			parEx.usage();
			abort();
		} catch (ConfigurationException configEx) {
			configEx.printStackTrace();
			abort();
		}
	}

	/**
	 * Loads the given file from the input directory and checks for duplicates based on the
	 * existence of specified strings in the name of the file. Returns true if there
	 * is no duplicate file with the same string ID and the file was loaded correctly.
	 *
	 * @param filename
	 * @param directoryInput
	 */
	void loadFileFromDirectory(String filename, String directoryInput) throws ConfigurationException {
		String fileID = InputFileIDs.contains(filename);
		String filenameWithFullPath = new File(directoryInput, filename).getAbsolutePath();

		if (fileID != null) {
			switch (fileID) {
				case "training":
					if (filenameTrainingData.length() > 0) {
						throw new ConfigurationException("Aborting, multiple training data files detected: "
								+ filenameWithFullPath + ", " + filenameTrainingData);
					} else {
						filenameTrainingData = filenameWithFullPath;
					}
					break;
				case "perturbations":
					if (filenamePerturbations.length() > 0) {
						throw new ConfigurationException("Aborting, multiple perturbation files detected: "
                                + filenameWithFullPath + ", " + filenamePerturbations);
					} else {
						filenamePerturbations = filenameWithFullPath;
					}
					break;
				case "modeloutputs":
					if (filenameModelOutputs.length() > 0) {
						throw new ConfigurationException("Aborting, multiple model output files detected: "
                                + filenameWithFullPath + ", " + filenameModelOutputs);
					} else {
						filenameModelOutputs = filenameWithFullPath;
					}
					break;
				case "config":
					if (filenameConfig.length() > 0) {
						throw new ConfigurationException("Aborting, multiple config files detected: "
                                + filenameWithFullPath + ", " + filenameConfig);
					} else {
						filenameConfig = filenameWithFullPath;
					}
					break;
				case "drugpanel":
					if (filenameDrugs.length() > 0) {
						throw new ConfigurationException("Aborting, multiple drug definition "
                                + "files detected: " + filenameWithFullPath + ", " + filenameDrugs);
					} else {
						filenameDrugs = filenameWithFullPath;
					}
					break;
				case "network":
					if (filenameNetwork.length() > 0) {
						throw new ConfigurationException("Aborting, multiple network files "
                                + "detected: " + filenameWithFullPath + ", " + filenameNetwork);
					} else {
						filenameNetwork = filenameWithFullPath;
					}
					break;
			}
		}
	}

	private void printInputFilesMessage() {
		System.out.println("\nGitsbe input files:" + "\nConfig:       \t" + filenameConfig
				+ "\nNetwork:      \t" + filenameNetwork +  "\nTraining data:\t"
                + filenameTrainingData + "\nModel outputs:\t" + filenameModelOutputs
                + "\n\nDrabme input files:" + "\nDrugs:        \t" + filenameDrugs
                + "\nPerturbations:\t" + filenamePerturbations + "\nModel outputs:\t"
                + filenameModelOutputs);
		System.out.println("\nOutput directory: " + directoryOutput + "\n");
	}

	private void runDrabme() {
		Thread thread = new Thread(new Drabme(
				projectName,
				filenameDrugs,
				filenamePerturbations,
				filenameModelOutputs,
				filenameConfig,
				directoryModels,
				directoryOutput,
				directoryTmpDrabme
		));
		execute(thread);
	}

	private void runGitsbe() {
		Thread thread = new Thread(new Gitsbe(
				projectName,
				filenameNetwork,
				filenameTrainingData,
				filenameModelOutputs,
				filenameConfig,
				filenameDrugs,
				directoryOutput,
				directoryTmpGitsbe
		));
		execute(thread);
	}

	private void execute(Thread thread) {
		thread.start();
		try {
			thread.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
			abort();
		}
	}

}