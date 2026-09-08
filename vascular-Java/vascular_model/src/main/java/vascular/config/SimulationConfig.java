package vascular.config;

import vascular.model.DomainMask;
import vascular.model.ModelParameters;

/**
 * Complete simulation configuration bundling model parameters, domain, and run settings.
 */
public record SimulationConfig(
        String name,
        String model,  // "HB", "HBP", "HBPM"
        ModelParameters parameters,
        DomainMask domain,
        double maxTime,
        long maxEvents,
        long observerInterval,
        long seed
) {

    /** Builder for programmatic construction. */
    public static final class Builder {
        private String name = "unnamed";
        private String model = "HB";
        private ModelParameters parameters = ModelParameters.defaultHB();
        private DomainMask domain;
        private double maxTime;
        private long maxEvents = 100_000_000;
        private long observerInterval = 1_000_000;
        private long seed = 42;

        public Builder name(String v) { name = v; return this; }
        public Builder model(String v) { model = v; return this; }
        public Builder parameters(ModelParameters v) { parameters = v; return this; }
        public Builder domain(DomainMask v) { domain = v; return this; }
        public Builder maxTime(double v) { maxTime = v; return this; }
        public Builder maxEvents(long v) { maxEvents = v; return this; }
        public Builder observerInterval(long v) { observerInterval = v; return this; }
        public Builder seed(long v) { seed = v; return this; }

        public SimulationConfig build() {
            if (domain == null) {
                throw new IllegalStateException("Domain mask is required");
            }
            return new SimulationConfig(name, model, parameters, domain, maxTime, maxEvents, observerInterval, seed);
        }
    }

    @Override
    public String toString() {
        int[] dims = domain.dimensions();
        return String.format("%s [%s] %dx%dx%d, maxEvents=%d, seed=%d",
                name, model, dims[0], dims[1], dims[2], maxEvents, seed);
    }
}
