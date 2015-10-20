package pb.g8;

/**
 * Created by naman on 10/20/15.
 */

public class Tuple {
    public final Long x;
    public final Long y;
    public Tuple(Long x, Long y) {
        this.x = x;
        this.y = y;
    }

    @Override
    public String toString() {
        return "(" + x + "," + y + ")";
    }

    public boolean containsId(long id) {
        return x == id || y == id;
    }

}