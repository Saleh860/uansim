package uansim;

import java.io.FileNotFoundException;
import java.io.PrintStream;

public class LogFile {

	static final PrintStream logFile = LogFile.initLogFile();

	static PrintStream initLogFile() {
		try {
			return new PrintStream("simulator.log");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return null;
	}

}
