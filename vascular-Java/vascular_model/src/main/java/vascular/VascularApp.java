package vascular;

import vascular.analysis.PatternAnalyzer;
import vascular.config.ParameterPresets;
import vascular.config.SimulationConfig;
import vascular.io.StateExporter;
import vascular.model.ModelParameters;
import vascular.model.Species;
import vascular.simulation.GillespieSimulator;
import vascular.simulation.SimulationObserver;
import vascular.visualization.VascularViewer;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Main entry point for the Vascular Patterning Model.
 *
 * <p>Usage:
 *   java VascularApp --gui                         Launch the GUI viewer
 *   java VascularApp --gui --DH=0.1 --DB=10.0      GUI with custom parameters
 *   java VascularApp --headless                    Run default 2D simulation
 *   java VascularApp --headless --k1=0.002 --DH=0.1  Headless with overrides
 *   java VascularApp --preset 0                    Run preset #0 headless
 *   java VascularApp --preset 2 --k9=0.01 --DM=3.0  Preset with overrides
 *   java VascularApp --list                        List available presets
 */
public class VascularApp {

    public static void main(String[] args) {
        // Step 1: Separate overrides (--KEY=VALUE) from commands
        Map<String, Double> overrides = new LinkedHashMap<>();
        List<String> commands = new ArrayList<>();

        for (String arg : args) {
            if (arg.startsWith("--") && arg.contains("=")) {
                // Parse --KEY=VALUE
                String kv = arg.substring(2); // remove leading --
                int eq = kv.indexOf('=');
                String key = kv.substring(0, eq);
                try {
                    double value = Double.parseDouble(kv.substring(eq + 1));
                    overrides.put(key, value);
                } catch (NumberFormatException e) {
                    System.err.println("Warning: invalid number in " + arg + " — ignored");
                }
            } else {
                commands.add(arg);
            }
        }

        // Step 2: Dispatch command
        String cmd = commands.isEmpty() ? "--gui" : commands.get(0);

        if ("--gui".equals(cmd)) {
            launchGUI(overrides);
        } else if ("--headless".equals(cmd)) {
            SimulationConfig config = ParameterPresets.hbDefault2D();
            config = applyOverrides(config, overrides);
            runHeadless(config);
        } else if ("--preset".equals(cmd) && commands.size() > 1) {
            try {
                int idx = Integer.parseInt(commands.get(1));
                SimulationConfig[] presets = ParameterPresets.allPresets();
                if (idx >= 0 && idx < presets.length) {
                    SimulationConfig config = applyOverrides(presets[idx], overrides);
                    runHeadless(config);
                } else {
                    System.err.println("Invalid preset index. Available: 0-" + (presets.length - 1));
                    listPresets(presets);
                }
            } catch (NumberFormatException e) {
                System.err.println("Invalid preset number: " + commands.get(1));
            }
        } else if ("--list".equals(cmd)) {
            listPresets(ParameterPresets.allPresets());
        } else {
            printUsage();
        }
    }

    /** Applies CLI overrides to a config's parameters. */
    private static SimulationConfig applyOverrides(SimulationConfig config, Map<String, Double> overrides) {
        if (overrides.isEmpty()) return config;

        ModelParameters newParams = config.parameters().applyOverrides(overrides);
        System.out.println("Parameter overrides: " + overrides);
        System.out.printf("  DB/DH ratio: %.1f → %.1f%n",
                config.parameters().DBoverDH(), newParams.DBoverDH());

        return new SimulationConfig(
                config.name() + "-custom",
                config.model(),
                newParams,
                config.domain(),
                config.maxTime(),
                config.maxEvents(),
                config.observerInterval(),
                config.seed()
        );
    }

    private static void launchGUI(Map<String, Double> overrides) {
        System.setProperty("awt.useSystemAAFontSettings", "on");
        System.setProperty("swing.aatext", "true");

        SwingUtilities.invokeLater(() -> {
            VascularViewer viewer = new VascularViewer();
            if (!overrides.isEmpty()) {
                viewer.setParameterOverrides(overrides);
                System.out.println("GUI launching with overrides: " + overrides);
            }
            viewer.setVisible(true);
        });
    }

    private static void runHeadless(SimulationConfig config) {
        System.out.println("=== Vascular Patterning Simulation ===");
        System.out.println("Configuration: " + config);
        System.out.println();

        GillespieSimulator sim = new GillespieSimulator(config.domain(), config.parameters(), config.seed());
        sim.setMaxEvents(config.maxEvents());
        sim.setMaxTime(config.maxTime());
        sim.setObserverInterval(config.observerInterval());

        // Progress observer
        sim.setObserver(new SimulationObserver() {
            private long lastReportTime = System.currentTimeMillis();
            private long lastEvents = 0;

            @Override
            public void onStep(long eventCount, double simTime, vascular.model.Lattice lattice) {
                long now = System.currentTimeMillis();
                long dt = now - lastReportTime;
                if (dt >= 2000) { // report every 2 seconds
                    long de = eventCount - lastEvents;
                    double rate = de * 1000.0 / dt;
                    System.out.printf("  Events: %,d | Time: %.2e | Rate: %,.0f/s | Total H: %,d | Max H: %,d%n",
                            eventCount, simTime, rate,
                            lattice.totalCount(Species.H),
                            lattice.maxCount(Species.H));
                    lastReportTime = now;
                    lastEvents = eventCount;
                }
            }

            @Override
            public void onComplete(long eventCount, double simTime) {
                System.out.println("  Simulation complete.");
            }
        });

        System.out.println("Starting simulation...");
        long startTime = System.currentTimeMillis();
        sim.runSync();
        long elapsed = System.currentTimeMillis() - startTime;

        // Final statistics
        System.out.println();
        System.out.printf("Completed in %.1f seconds%n", elapsed / 1000.0);
        System.out.printf("Total events: %,d%n", sim.eventCount());
        System.out.printf("Simulation time: %.4e%n", sim.simTime());
        System.out.println(StateExporter.summarize(sim.lattice()));

        // Pattern analysis
        int[] dims = config.domain().dimensions();
        PatternAnalyzer analyzer = new PatternAnalyzer(50);
        List<PatternAnalyzer.BundleInfo> bundles = analyzer.detectBundles2D(sim.lattice(), dims[2] / 2);

        System.out.println("=== Pattern Analysis ===");
        System.out.printf("Detected %d vascular bundles (threshold H >= 50)%n", bundles.size());
        double wavelength = analyzer.estimateWavelength(bundles);
        if (!Double.isNaN(wavelength)) {
            System.out.printf("Estimated wavelength: %.2f voxels%n", wavelength);
        }
        System.out.printf("Average bundle size: %.1f voxels%n",
                bundles.stream().mapToInt(b -> b.area()).average().orElse(0));

        // Export results
        try {
            File outputDir = new File("output");
            outputDir.mkdirs();
            File csvFile = new File(outputDir, config.name().replace(' ', '_') + "_final.csv");
            StateExporter.exportCSV(sim.lattice(), dims[2] / 2, csvFile);
            System.out.println("Results exported to: " + csvFile.getAbsolutePath());

            File ckptFile = new File(outputDir, config.name().replace(' ', '_') + ".ckpt");
            StateExporter.saveCheckpoint(sim.lattice(), sim.simTime(), sim.eventCount(), ckptFile);
            System.out.println("Checkpoint saved to: " + ckptFile.getAbsolutePath());
        } catch (IOException e) {
            System.err.println("Failed to export results: " + e.getMessage());
        }
    }

    private static void listPresets(SimulationConfig[] presets) {
        System.out.println("Available presets:");
        for (int i = 0; i < presets.length; i++) {
            System.out.printf("  %d: %s%n", i, presets[i]);
        }
    }

    private static void printUsage() {
        System.out.println("Vascular Patterning Model — Turing-like Stochastic Reaction-Diffusion");
        System.out.println("Based on Hearn (2019) PLoS ONE 14(7): e0219055");
        System.out.println();
        System.out.println("Usage:");
        System.out.println("  java VascularApp --gui                       Launch GUI viewer");
        System.out.println("  java VascularApp --gui --DH=0.1 --DB=10.0    GUI with custom params");
        System.out.println("  java VascularApp --headless                  Run default headless");
        System.out.println("  java VascularApp --headless --k1=0.002       Headless with overrides");
        System.out.println("  java VascularApp --preset N                  Run preset #N");
        System.out.println("  java VascularApp --preset 2 --k9=0.01        Preset with overrides");
        System.out.println("  java VascularApp --list                      List all presets");
        System.out.println();
        System.out.println("Override any parameter with --KEY=VALUE:");
        System.out.println("  k1-k12, DH, DB, DM, a, V");
        System.out.println("  Example: --DH=0.2 --DB=20.0 --k2=1.0");
    }
}
