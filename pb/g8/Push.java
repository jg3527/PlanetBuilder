package pb.g8;

public class Push {

    public long asteroid_id;
    public double energy;
    public double direction;
    public long time_of_push;
    public long time_of_collision;

    public Push(long asteroid_id, double energy, double direction, long time_of_push) {
        this.asteroid_id = asteroid_id;
        this.energy = energy;
        this.direction = direction;
        this.time_of_push = time_of_push;
        // we don't know it yet
        this.time_of_collision = -1;
    }

    public Push(long asteroid_id, double energy, double direction, long time_of_push, long time_of_collision) {
        this.asteroid_id = asteroid_id;
        this.energy = energy;
        this.direction = direction;
        this.time_of_push = time_of_push;
        this.time_of_collision = time_of_collision;
    }

    @Override
    public String toString() {
        return "Push{" +
                "asteroid_id=" + asteroid_id +
                ", energy=" + energy +
                ", direction=" + direction +
                ", time_of_push=" + time_of_push +
                ", time_of_collision=" + time_of_collision +
                '}';
    }
}