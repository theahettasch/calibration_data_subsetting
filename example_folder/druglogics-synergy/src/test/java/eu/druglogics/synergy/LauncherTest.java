package eu.druglogics.synergy;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import javax.naming.ConfigurationException;
import java.io.File;

import static org.junit.jupiter.api.Assertions.assertThrows;

public class LauncherTest {
	@TempDir
	File tempDir;

	@Test
	void test_load_duplicate_file_from_dir() throws ConfigurationException {
		Launcher launcher = new Launcher();
		launcher.loadFileFromDirectory("training", tempDir.getAbsolutePath());

		// load file that has the same id (here: "training")
		assertThrows(ConfigurationException.class, () ->
			launcher.loadFileFromDirectory("training_2", tempDir.getAbsolutePath()));
	}
}
