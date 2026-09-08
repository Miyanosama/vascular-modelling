package vascular.visualization;

import vascular.config.ParameterPresets;
import vascular.config.SimulationConfig;
import vascular.model.*;
import vascular.simulation.GillespieSimulator;
import vascular.simulation.SimulationObserver;

import javax.swing.*;
import java.awt.*;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Main Swing GUI for the vascular patterning simulation.
 *
 * <p>Layout:
 * - Left: 2D heatmap
 * - Right (scrollable): Controls + Parameter Editor + Stats
 * - Bottom: Status bar
 */
public class VascularViewer extends JFrame {

    private GillespieSimulator simulator;
    private SimulationConfig config;
    private Thread simulationThread;
    private Map<String, Double> pendingOverrides;

    // Components
    private LatticePanel2D latticePanel;
    private JLabel timeLabel, eventLabel, fpsLabel, statusLabel;
    private JButton startButton, pauseButton, stopButton;
    private JComboBox<String> speciesCombo, colorMapCombo;
    private JCheckBox logScaleCheck;

    // Parameter editor fields
    private final Map<String, JTextField> paramFields = new LinkedHashMap<>();

    // Stats
    private long lastObserverEvent;
    private long lastObserverTime;

    public VascularViewer() {
        super("Vascular Patterning Model — Turing-like Stochastic Reaction-Diffusion");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        buildUI();
        setSize(950, 750);
        setLocationRelativeTo(null);
    }

    /** Called from CLI when --gui is used with --KEY=VALUE overrides. */
    public void setParameterOverrides(Map<String, Double> overrides) {
        this.pendingOverrides = overrides;
    }

    private void buildUI() {
        setLayout(new BorderLayout());

        // ---- Center: Lattice view ----
        latticePanel = new LatticePanel2D();
        latticePanel.setPreferredSize(new Dimension(500, 500));
        add(latticePanel, BorderLayout.CENTER);

        // ---- Right: Controls + Parameters (scrollable) ----
        JPanel rightPanel = new JPanel();
        rightPanel.setLayout(new BoxLayout(rightPanel, BoxLayout.Y_AXIS));
        rightPanel.setBorder(BorderFactory.createEmptyBorder(5, 8, 5, 8));

        buildControlSection(rightPanel);
        rightPanel.add(Box.createVerticalStrut(10));
        buildParameterSection(rightPanel);
        rightPanel.add(Box.createVerticalStrut(10));
        buildDisplaySection(rightPanel);
        rightPanel.add(Box.createVerticalStrut(10));
        buildStatsSection(rightPanel);

        JScrollPane scrollPane = new JScrollPane(rightPanel,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setPreferredSize(new Dimension(280, 700));
        scrollPane.setBorder(BorderFactory.createEmptyBorder());
        add(scrollPane, BorderLayout.EAST);

        // ---- Bottom: Status bar ----
        JPanel bottomPanel = new JPanel(new BorderLayout());
        bottomPanel.setBorder(BorderFactory.createLoweredBevelBorder());
        JLabel infoLabel = new JLabel("  Stochastic Reaction-Diffusion SSA | Hearn (2019) PLoS ONE | Vascular Patterning Model");
        infoLabel.setFont(new Font("SansSerif", Font.PLAIN, 11));
        bottomPanel.add(infoLabel, BorderLayout.WEST);
        add(bottomPanel, BorderLayout.SOUTH);
    }

    // ========== Control Section ==========

    private void buildControlSection(JPanel parent) {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        panel.setBorder(BorderFactory.createTitledBorder("Controls"));
        panel.setMaximumSize(new Dimension(260, 200));

        JPanel buttonPanel = new JPanel(new GridLayout(3, 1, 4, 4));
        startButton = new JButton("▶ Start");
        pauseButton = new JButton("⏸ Pause");
        stopButton = new JButton("⏹ Stop");
        pauseButton.setEnabled(false);
        stopButton.setEnabled(false);

        startButton.addActionListener(e -> onStart());
        pauseButton.addActionListener(e -> onPause());
        stopButton.addActionListener(e -> onStop());

        buttonPanel.add(startButton);
        buttonPanel.add(pauseButton);
        buttonPanel.add(stopButton);
        buttonPanel.setMaximumSize(new Dimension(240, 90));
        panel.add(buttonPanel);
        panel.add(Box.createVerticalStrut(8));

        // Preset selector
        panel.add(new JLabel("Load Preset:"));
        String[] presets = {"HB Default 2D", "HBP Default 3D", "HBPM Default 3D",
                "HBPM High M", "Large Stem", "Medium Stem", "Small Stem"};
        JComboBox<String> presetCombo = new JComboBox<>(presets);
        presetCombo.setMaximumSize(new Dimension(240, 25));
        presetCombo.addActionListener(e -> loadPreset(presetCombo.getSelectedIndex()));
        panel.add(presetCombo);

        parent.add(panel);
    }

    // ========== Parameter Editor Section ==========

    private void buildParameterSection(JPanel parent) {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        panel.setBorder(BorderFactory.createTitledBorder("Parameters"));

        // Define parameters to show: key, label, tooltip
        String[][] paramDefs = {
                {"k1", "k1 (autocatalysis)", "2H+B→3H rate. Higher = stronger autocatalysis"},
                {"k2", "k2 (B production)", "B production rate. Higher = more substrate"},
                {"k3", "k3 (B degradation)", "B degradation rate"},
                {"k4", "k4 (H degradation)", "H degradation rate. Higher = less activator"},
                {"k5", "k5 (H production)", "Basal H production rate"},
                {"DH", "DH (H diffusion)", "Activator diffusion. Higher = larger spots"},
                {"DB", "DB (B diffusion)", "Substrate diffusion. DB/DH controls spot density"},
        };

        for (String[] def : paramDefs) {
            JPanel row = new JPanel(new BorderLayout(4, 0));
            row.setMaximumSize(new Dimension(240, 24));
            JLabel label = new JLabel(def[1]);
            label.setFont(new Font("SansSerif", Font.PLAIN, 11));
            label.setToolTipText(def[2]);
            label.setPreferredSize(new Dimension(130, 20));
            row.add(label, BorderLayout.WEST);

            JTextField field = new JTextField(6);
            field.setFont(new Font("Monospaced", Font.PLAIN, 11));
            field.setToolTipText(def[2]);
            field.setMaximumSize(new Dimension(100, 22));
            row.add(field, BorderLayout.EAST);

            panel.add(row);
            panel.add(Box.createVerticalStrut(2));
            paramFields.put(def[0], field);
        }

        panel.add(Box.createVerticalStrut(5));

        // Apply button
        JButton applyBtn = new JButton("Apply & Reset");
        applyBtn.setToolTipText("Stop simulation, apply new parameters, reset lattice");
        applyBtn.setMaximumSize(new Dimension(240, 28));
        applyBtn.addActionListener(e -> applyParameters());
        panel.add(applyBtn);

        parent.add(panel);
    }

    private void applyParameters() {
        if (config == null) return;

        // Stop current simulation
        if (simulator != null) {
            simulator.stop();
        }

        // Read parameters from text fields
        Map<String, Double> overrides = new LinkedHashMap<>();
        for (var entry : paramFields.entrySet()) {
            String text = entry.getValue().getText().trim();
            if (!text.isEmpty()) {
                try {
                    overrides.put(entry.getKey(), Double.parseDouble(text));
                } catch (NumberFormatException ex) {
                    statusLabel.setText("Invalid number: " + entry.getKey() + "=" + text);
                    statusLabel.setForeground(Color.RED);
                    return;
                }
            }
        }

        if (overrides.isEmpty()) {
            statusLabel.setText("No parameters changed");
            statusLabel.setForeground(Color.BLACK);
            return;
        }

        // Build new config with overridden parameters
        ModelParameters newParams = config.parameters().applyOverrides(overrides);
        this.config = new SimulationConfig(
                config.name() + "-edited", config.model(), newParams,
                config.domain(), config.maxTime(), config.maxEvents(),
                config.observerInterval(), config.seed()
        );

        // Recreate simulator
        this.simulator = new GillespieSimulator(config.domain(), config.parameters(), config.seed());
        this.simulator.setMaxEvents(config.maxEvents());
        this.simulator.setMaxTime(config.maxTime());
        this.simulator.setObserverInterval(config.observerInterval());
        this.simulator.setObserver(new UIObserver());

        latticePanel.setLattice(simulator.lattice());
        updateStats(0, 0.0);
        statusLabel.setText("Ready — " + config.name() + " | DB/DH=" +
                String.format("%.0f", config.parameters().DBoverDH()));
        statusLabel.setForeground(new Color(0, 100, 0));

        startButton.setEnabled(true);
        pauseButton.setEnabled(false);
        stopButton.setEnabled(false);

        System.out.println("Parameters applied: " + overrides);
    }

    /** Refresh parameter text fields from current config. */
    private void refreshParamFields() {
        if (config == null) return;
        ModelParameters p = config.parameters();
        setField("k1", p.k1());
        setField("k2", p.k2());
        setField("k3", p.k3());
        setField("k4", p.k4());
        setField("k5", p.k5());
        setField("DH", p.DH());
        setField("DB", p.DB());
    }

    private void setField(String key, double value) {
        JTextField f = paramFields.get(key);
        if (f != null) f.setText(String.valueOf(value));
    }

    // ========== Display Section ==========

    private void buildDisplaySection(JPanel parent) {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        panel.setBorder(BorderFactory.createTitledBorder("Display"));
        panel.setMaximumSize(new Dimension(260, 120));

        panel.add(new JLabel("Species:"));
        speciesCombo = new JComboBox<>(new String[]{"H", "B", "P", "M"});
        speciesCombo.setMaximumSize(new Dimension(240, 25));
        speciesCombo.addActionListener(e -> {
            String sel = (String) speciesCombo.getSelectedItem();
            if (sel != null) latticePanel.setDisplaySpecies(Species.valueOf(sel));
        });
        panel.add(speciesCombo);

        panel.add(Box.createVerticalStrut(5));
        panel.add(new JLabel("Color Map:"));
        colorMapCombo = new JComboBox<>(new String[]{"Viridis", "Magma", "Inferno", "Plasma", "Blue-Red"});
        colorMapCombo.setMaximumSize(new Dimension(240, 25));
        colorMapCombo.addActionListener(e -> {
            String sel = (String) colorMapCombo.getSelectedItem();
            if (sel != null) {
                latticePanel.setColorMap(ColorMap.valueOf(sel.toUpperCase().replace('-', '_')));
            }
        });
        panel.add(colorMapCombo);

        logScaleCheck = new JCheckBox("Log Scale", true);
        logScaleCheck.addActionListener(e -> latticePanel.setLogScale(logScaleCheck.isSelected()));
        panel.add(logScaleCheck);

        parent.add(panel);
    }

    // ========== Stats Section ==========

    private void buildStatsSection(JPanel parent) {
        JPanel panel = new JPanel(new GridLayout(5, 1, 2, 2));
        panel.setBorder(BorderFactory.createTitledBorder("Statistics"));
        panel.setMaximumSize(new Dimension(260, 150));

        timeLabel = new JLabel("Time: 0.00");
        eventLabel = new JLabel("Events: 0");
        fpsLabel = new JLabel("Rate: --");
        statusLabel = new JLabel("Status: Ready");
        JLabel dbdhLabel = new JLabel("DB/DH: --");

        panel.add(timeLabel);
        panel.add(eventLabel);
        panel.add(fpsLabel);
        panel.add(dbdhLabel);
        panel.add(statusLabel);

        // Store reference to DB/DH label for updates
        panel.putClientProperty("dbdhLabel", dbdhLabel);

        parent.add(panel);
    }

    // ========== Preset Loading ==========

    private void loadPreset(int index) {
        SimulationConfig cfg = switch (index) {
            case 0 -> ParameterPresets.hbDefault2D();
            case 1 -> ParameterPresets.hbpDefault3D();
            case 2 -> ParameterPresets.hbpmDefault3D();
            case 3 -> ParameterPresets.hbpmHighM();
            case 4 -> ParameterPresets.largeStem2D();
            case 5 -> ParameterPresets.mediumStem2D();
            case 6 -> ParameterPresets.smallStem2D();
            default -> ParameterPresets.hbDefault2D();
        };

        // Apply any pending CLI overrides
        if (pendingOverrides != null && !pendingOverrides.isEmpty()) {
            ModelParameters newParams = cfg.parameters().applyOverrides(pendingOverrides);
            cfg = new SimulationConfig(cfg.name() + "-cli", cfg.model(), newParams,
                    cfg.domain(), cfg.maxTime(), cfg.maxEvents(),
                    cfg.observerInterval(), cfg.seed());
            pendingOverrides = null;
        }

        this.config = cfg;
        this.simulator = new GillespieSimulator(cfg.domain(), cfg.parameters(), cfg.seed());
        this.simulator.setMaxEvents(cfg.maxEvents());
        this.simulator.setMaxTime(cfg.maxTime());
        this.simulator.setObserverInterval(cfg.observerInterval());
        this.simulator.setObserver(new UIObserver());

        latticePanel.setLattice(simulator.lattice());
        refreshParamFields();
        updateStats(0, 0.0);
        statusLabel.setText("Ready — " + cfg.name() + " | DB/DH=" +
                String.format("%.0f", cfg.parameters().DBoverDH()));
        statusLabel.setForeground(Color.BLACK);

        startButton.setEnabled(true);
        pauseButton.setEnabled(false);
        stopButton.setEnabled(false);
    }

    // ========== Simulation Controls ==========

    private void onStart() {
        if (simulator == null) return;

        startButton.setEnabled(false);
        pauseButton.setEnabled(true);
        stopButton.setEnabled(true);
        statusLabel.setText("Running...");
        statusLabel.setForeground(new Color(0, 120, 0));

        lastObserverEvent = 0;
        lastObserverTime = System.nanoTime();

        simulationThread = Thread.ofVirtual().start(() -> {
            simulator.runSync();
            SwingUtilities.invokeLater(() -> {
                statusLabel.setText("Complete (" + String.format("%,d", simulator.eventCount()) + " events)");
                statusLabel.setForeground(new Color(0, 0, 180));
                startButton.setEnabled(true);
                pauseButton.setEnabled(false);
                stopButton.setEnabled(false);
            });
        });
    }

    private void onPause() {
        if (simulator != null) {
            simulator.pause();
            statusLabel.setText("Paused");
            statusLabel.setForeground(new Color(180, 120, 0));
            startButton.setEnabled(true);
            pauseButton.setEnabled(false);
        }
    }

    private void onStop() {
        if (simulator != null) {
            simulator.stop();
            statusLabel.setText("Stopped");
            statusLabel.setForeground(new Color(180, 0, 0));
            startButton.setEnabled(true);
            pauseButton.setEnabled(false);
            stopButton.setEnabled(false);
        }
    }

    private void updateStats(long events, double time) {
        timeLabel.setText(String.format("Time: %.2e", time));
        eventLabel.setText(String.format("Events: %,d", events));
    }

    // ========== Simulation Observer ==========

    private class UIObserver implements SimulationObserver {
        @Override
        public void onStep(long eventCount, double simTime, Lattice lattice) {
            SwingUtilities.invokeLater(() -> {
                updateStats(eventCount, simTime);

                long now = System.nanoTime();
                long dt = now - lastObserverTime;
                if (dt > 0) {
                    long de = eventCount - lastObserverEvent;
                    fpsLabel.setText(String.format("Rate: %,.0f evt/s", de * 1e9 / dt));
                }
                lastObserverTime = now;
                lastObserverEvent = eventCount;

                latticePanel.repaint();
            });
        }

        @Override
        public void onComplete(long eventCount, double simTime) {
            SwingUtilities.invokeLater(() -> {
                updateStats(eventCount, simTime);
                statusLabel.setText("Complete (" + String.format("%,d", eventCount) + " events)");
                statusLabel.setForeground(new Color(0, 0, 180));
                startButton.setEnabled(true);
                pauseButton.setEnabled(false);
                stopButton.setEnabled(false);
                latticePanel.repaint();
            });
        }
    }
}
