package pb.g8;

/**
 * Created by naman on 10/11/15.
 */
public class Push {

    public int asteroid_id;
    public double energy;
    public double direction;
    public long time_of_push;

    public Push(int asteroid_id, double energy, double direction, long time_of_push) {
        this.asteroid_id = asteroid_id;
        this.energy = energy;
        this.direction = direction;
        this.time_of_push = time_of_push;
    }
}
