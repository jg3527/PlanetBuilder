package pb.g6;

public class CandidatePush implements Comparable<CandidatePush> {
    public int ast_idx;
    public long time;
    public double energy;
    public double angle;

    public CandidatePush() {
        this.ast_idx = -1;
        this.time = 0;
        this.energy = 0.0;
        this.angle = 0.0;
    }

    public CandidatePush(int idx, long t, double e, double a) {
        this.ast_idx = idx;
        this.time = t;
        this.energy = e;
        this.angle = a;
    }

    @Override
    public int compareTo(CandidatePush o) {
        return (int) Math.ceil(this.energy - o.energy);
    }
}
