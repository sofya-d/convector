package chimeric_solver;

import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.concurrent.*;

final class SmithWaterman implements Callable<Boolean> {
    private Integer maxPercOfErrors = 20;
    private Double gapCost = -1.0;
    private Double matchCost = 1.0;
    private Double misMatchCost = -1.0;
    private String chimera;
    private String reference;

    public SmithWaterman(String chim, String ref) {
        this.chimera = chim;
        this.reference = ref;
    }

    public Boolean call() {
        if (chimera.equals(reference)) {
            return false;
        }

        int maxNumOfErrors = (int) Math.ceil(0.01 * chimera.length() * maxPercOfErrors) + 1;
        int minScoreToAlign = (chimera.length() - maxNumOfErrors);

        Double maxScore = 0.0;
        Double[][] matrixOfDp = new Double[reference.length() + 1][chimera.length() + 1];
        
        for (int i = 0; i < reference.length() + 1; i++) {
            for (int j = 0; j < chimera.length() + 1; j++) {
                matrixOfDp[i][j] = 0.0;
            }
        }
        
        for (int i = 1; i < 1 + reference.length(); i++) {
            for (int j = 1; j < chimera.length() + 1; j++) {
                Double leftGap = gapCost;
                Double rightGap = gapCost;
                Double match = 0.0;
                if (reference.charAt(i - 1) == chimera.charAt(j - 1)) {
                    match = matchCost;
                } else if (reference.charAt(i - 1) == '*') {
                    match = 0.0;
                    leftGap = 0.0;
                    rightGap = 0.0;
                } else match = misMatchCost;
                
                matrixOfDp[i][j] = Math.max(Math.max(matrixOfDp[i][j - 1] + leftGap, matrixOfDp[i - 1][j]) + rightGap, Math.max(matrixOfDp[i - 1][j - 1] + match, 0));
                
                if (matrixOfDp[i][j] > maxScore) {
                    maxScore = matrixOfDp[i][j];
                }
                if (maxScore > minScoreToAlign) {
                    return true;
                }
            }
        }
        return false;
    }
}

class Solver {
    private ArrayList<Boolean> canBeAligned = new ArrayList<>();
    private ArrayList<Integer> coverageIncrease = new ArrayList<>();
    private ArrayList<String> references = new ArrayList<>();
    private ArrayList<String> chimeras = new ArrayList<>();

    // Variables are now set via constructor
    private String fileWithChimeras;
    private String fileWithReferences;
    private String tmpOutput;

    private void parseFiles() {
        try (BufferedReader br = new BufferedReader(new FileReader(fileWithChimeras))) {
            for (String line; (line = br.readLine()) != null; ) {
                String[] tmpListForMap = line.split("\\t");
                chimeras.add(tmpListForMap[0]);
                canBeAligned.add(false);
            }
        } catch (IOException e) {
            System.err.println("\nError reading chimeras file: " + fileWithChimeras);
            System.err.println(e.getMessage());
            System.exit(1);
        }

        try (BufferedReader br = new BufferedReader(new FileReader(fileWithReferences))) {
            for (String line; (line = br.readLine()) != null; ) {
                String[] tmpListForMap = line.split("\\t");
                references.add(tmpListForMap[0]);
                coverageIncrease.add(0);
            }
        } catch (IOException e) {
            System.err.println("\nError reading references file: " + fileWithReferences);
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    private void writeOutput() {
        Path file = Paths.get(tmpOutput);
        try {
            Files.deleteIfExists(file);
            Files.createFile(file);
        } catch (IOException e) {
            System.err.format("createFile error: %s%n", e);
        }

        try (BufferedWriter bw = new BufferedWriter(new FileWriter(tmpOutput))) {
            for (Integer elem : coverageIncrease) {
                bw.write(elem + "\n");
            }
        } catch (IOException e) {
            System.err.format("Something went wrong with output file: %s%n", tmpOutput);
        }
    }

    public Solver(String chimeraPath, String referencePath, String outputPath) {
        this.fileWithChimeras = chimeraPath;
        this.fileWithReferences = referencePath;
        this.tmpOutput = outputPath;
        
        parseFiles();
        
        ExecutorService service = Executors.newFixedThreadPool(16);
        try {
            for (int i = 0; i < chimeras.size(); i++) {
                CompletionService<Boolean> pool = new ExecutorCompletionService<>(service);
                for (int j = 0; j < references.size(); j++) {
                    pool.submit(new SmithWaterman(chimeras.get(i), references.get(j)));
                }

                for (int j = 0; j < references.size(); j++) {
                    try {
                        Boolean result = pool.take().get(); // Efficiently wait for the first available result
                        if (result) {
                            canBeAligned.set(i, true);
                            coverageIncrease.set(j, coverageIncrease.get(j) + 1);
                            break; 
                        }
                    } catch (InterruptedException | ExecutionException e) {
                        e.printStackTrace();
                    }
                }
            }
        } finally {
            service.shutdown();
        }
        writeOutput();
    }
}

public class Main {
    public static void main(String[] args) {
        if (args.length < 3) {
            System.out.println("Usage: java chimeric_solver.Main <chimeras_path> <references_path> <output_path>");
            System.exit(1);
        }

        String chimerasPath = args[0];
        String referencesPath = args[1];
        String outputPath = args[2];

        new Solver(chimerasPath, referencesPath, outputPath);
    }
}