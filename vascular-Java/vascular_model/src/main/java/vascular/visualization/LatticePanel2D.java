package vascular.visualization;

import vascular.model.*;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;

/**
 * Renders a 2D heatmap of molecular concentrations for a single z-slice.
 * Uses a BufferedImage for fast pixel-level rendering.
 */
public class LatticePanel2D extends JPanel {

    private transient Lattice lattice;
    private DomainMask domain;
    private Species displaySpecies = Species.H;
    private ColorMap colorMap = ColorMap.VIRIDIS;
    private boolean logScale = true;

    private int pixelSize = 4; // pixels per voxel
    private BufferedImage buffer;

    public LatticePanel2D() {
        setPreferredSize(new Dimension(400, 400));
        setBackground(Color.BLACK);
    }

    public void setLattice(Lattice lattice) {
        this.lattice = lattice;
        this.domain = lattice != null ? lattice.domain() : null;
        repaint();
    }

    public void setDisplaySpecies(Species species) { this.displaySpecies = species; repaint(); }
    public void setColorMap(ColorMap cm) { this.colorMap = cm; repaint(); }
    public void setLogScale(boolean log) { this.logScale = log; repaint(); }
    public Species displaySpecies() { return displaySpecies; }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        if (lattice == null || domain == null) {
            g.setColor(Color.DARK_GRAY);
            g.drawString("No simulation data", getWidth() / 2 - 50, getHeight() / 2);
            return;
        }

        int[] dims = domain.dimensions();
        int width = dims[0];
        int height = dims[1];

        // Ensure buffer is the right size
        int imgW = width * pixelSize;
        int imgH = height * pixelSize;
        if (buffer == null || buffer.getWidth() != imgW || buffer.getHeight() != imgH) {
            buffer = new BufferedImage(imgW, imgH, BufferedImage.TYPE_INT_RGB);
        }

        int maxCount = lattice.maxCount(displaySpecies);
        if (maxCount <= 0) maxCount = 1;

        Graphics2D bg = buffer.createGraphics();
        bg.setColor(Color.BLACK);
        bg.fillRect(0, 0, imgW, imgH);

        double cx = domain.center()[0];
        double cy = domain.center()[1];
        double radius = domain.radius();

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                // Check if inside domain
                double dx = x - cx;
                double dy = y - cy;
                if (dx * dx + dy * dy > radius * radius) {
                    continue; // outside domain, leave black
                }

                int count = lattice.get(displaySpecies, x, y, 0);
                Color color;
                if (count <= 0) {
                    color = new Color(10, 10, 20); // very dark blue for empty voxels
                } else if (logScale) {
                    color = colorMap.mapLog(count, maxCount);
                } else {
                    color = colorMap.mapLinear(count, maxCount);
                }

                bg.setColor(color);
                bg.fillRect(x * pixelSize, y * pixelSize, pixelSize, pixelSize);
            }
        }
        bg.dispose();

        // Scale and draw
        int panelW = getWidth();
        int panelH = getHeight();
        double scale = Math.min((double) panelW / imgW, (double) panelH / imgH);
        int drawW = (int) (imgW * scale);
        int drawH = (int) (imgH * scale);
        int drawX = (panelW - drawW) / 2;
        int drawY = (panelH - drawH) / 2;

        g.drawImage(buffer, drawX, drawY, drawW, drawH, null);

        // Draw domain boundary
        g.setColor(new Color(80, 80, 80));
        g.drawOval(
                drawX + (int)((cx - radius) * pixelSize * scale),
                drawY + (int)((cy - radius) * pixelSize * scale),
                (int)(2 * radius * pixelSize * scale),
                (int)(2 * radius * pixelSize * scale)
        );

        // Color scale legend
        drawColorBar(g, panelW - 30, 20, 16, panelH - 60, maxCount);
    }

    private void drawColorBar(Graphics g, int x, int y, int w, int h, int maxCount) {
        for (int i = 0; i < h; i++) {
            double t = 1.0 - (double) i / h;
            Color c;
            if (logScale) {
                c = colorMap.map(t);
            } else {
                c = colorMap.map(t);
            }
            g.setColor(c);
            g.fillRect(x, y + i, w, 1);
        }
        g.setColor(Color.WHITE);
        g.drawRect(x, y, w, h);
        g.setFont(new Font("SansSerif", Font.PLAIN, 10));
        g.drawString(String.valueOf(maxCount), x + w + 4, y + 10);
        g.drawString("0", x + w + 4, y + h);
        g.drawString("[" + displaySpecies + "]", x + w + 4, y + h / 2);
    }
}
