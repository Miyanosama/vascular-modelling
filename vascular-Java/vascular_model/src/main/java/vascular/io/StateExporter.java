package vascular.io;

import vascular.model.*;

import java.io.*;
import java.nio.charset.StandardCharsets;

/**
 * Exports simulation state to various formats for external visualization and analysis.
 */
public final class StateExporter {

    private StateExporter() {}

    /**
     * Exports a 2D slice of the lattice as a CSV file.
     * Format: x,y,count_H,count_B,count_P,count_M
     */
    public static void exportCSV(Lattice lattice, int zSlice, File outputFile) throws IOException {
        DomainMask domain = lattice.domain();
        int[] dims = domain.dimensions();

        try (PrintWriter writer = new PrintWriter(
                new OutputStreamWriter(new FileOutputStream(outputFile), StandardCharsets.UTF_8))) {
            writer.println("x,y,H,B,P,M");

            for (int y = 0; y < dims[1]; y++) {
                for (int x = 0; x < dims[0]; x++) {
                    if (!domain.isInside(x, y, zSlice)) continue;
                    writer.printf("%d,%d,%d,%d,%d,%d%n",
                            x, y,
                            lattice.get(Species.H, x, y, zSlice),
                            lattice.get(Species.B, x, y, zSlice),
                            lattice.get(Species.P, x, y, zSlice),
                            lattice.get(Species.M, x, y, zSlice)
                    );
                }
            }
        }
    }

    /**
     * Exports a summary of the lattice state (total counts per species).
     */
    public static String summarize(Lattice lattice) {
        StringBuilder sb = new StringBuilder();
        sb.append("=== Lattice Summary ===\n");
        sb.append(String.format("Total voxels: %d\n", lattice.totalVoxels()));
        for (Species s : Species.values()) {
            int total = lattice.totalCount(s);
            int max = lattice.maxCount(s);
            sb.append(String.format("  %s: total=%d, max=%d\n", s, total, max));
        }
        return sb.toString();
    }

    /**
     * Saves the full 3D lattice state as a compressed binary file.
     */
    public static void saveCheckpoint(Lattice lattice, double simTime, long eventCount, File file) throws IOException {
        try (DataOutputStream dos = new DataOutputStream(
                new BufferedOutputStream(new FileOutputStream(file)))) {
            // Header
            dos.writeLong(0x56415343434B5054L); // "VASCICKPT" magic
            dos.writeInt(1); // version
            dos.writeDouble(simTime);
            dos.writeLong(eventCount);

            int[] dims = lattice.domain().dimensions();
            dos.writeInt(dims[0]);
            dos.writeInt(dims[1]);
            dos.writeInt(dims[2]);

            int n = lattice.totalVoxels();
            dos.writeInt(n);

            for (Species s : Species.values()) {
                int[] counts = lattice.getCounts(s);
                for (int i = 0; i < n; i++) {
                    dos.writeInt(counts[i]);
                }
            }
        }
    }

    /**
     * Loads a checkpoint into a lattice.
     */
    public static CheckpointData loadCheckpoint(DomainMask domain, File file) throws IOException {
        try (DataInputStream dis = new DataInputStream(
                new BufferedInputStream(new FileInputStream(file)))) {

            long magic = dis.readLong();
            if (magic != 0x56415343434B5054L) {
                throw new IOException("Invalid checkpoint file (bad magic)");
            }
            int version = dis.readInt();
            double simTime = dis.readDouble();
            long eventCount = dis.readLong();

            int w = dis.readInt(), h = dis.readInt(), d = dis.readInt();
            int n = dis.readInt();

            Lattice lattice = new Lattice(domain);

            for (Species s : Species.values()) {
                for (int i = 0; i < n; i++) {
                    lattice.set(s, i, dis.readInt());
                }
            }

            return new CheckpointData(lattice, simTime, eventCount, w, h, d);
        }
    }

    public record CheckpointData(Lattice lattice, double simTime, long eventCount,
                                  int width, int height, int depth) {}
}
