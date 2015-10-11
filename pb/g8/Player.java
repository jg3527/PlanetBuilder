package pb.g8;


import java.util.ArrayList;
import java.util.Random;

import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;
import pb.sim.Orbit;
import pb.sim.Point;

public class Player implements pb.sim.Player {

    //iteration number
    int iteration =1;

    // used to pick asteroid and velocity boost randomly
    private Random random = new Random();

    // current time, time limit
    private long time = -1;
    private long time_limit = -1;


    private Point origin = new Point(0,0);

    // time until next push
    private long time_of_push = 0;

    // number of retries
    private int retries_per_turn = 1;
    private int turns_per_retry = 3;

    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit) {
        if (Orbit.dt() != 24 * 60 * 60)
            throw new IllegalStateException("Time quantum is not a day");
        this.time_limit = time_limit;
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids, double[] energy, double[] direction) {
        // if not yet time to push do nothing
        if (++time <= time_of_push) return;
        //System.out.println("Year: " + (1 + time / 365));
        //System.out.println("Day: "  + (1 + time % 365));
        for (int retry = 1 ; retry <= retries_per_turn ; ++retry) {
            // pick a random asteroid and get its velocity
            int i = 0;
            if(iteration == 1) {
                i = getFarthestAsteroid(asteroids);
            } else {
                i = getLightestAsteroid(asteroids);
            }

//            i = getClosestPairAsteroid(asteroids);
            Point v = asteroids[i].orbit.velocityAt(time);
            // add 5-50% of current velocity in magnitude
            //System.out.println("Try: " + retry + " / " + retries_per_turn);
            double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
            double v2 = 0.0;
            for (double k=0; k< 5; k=k+0.1) {
                v2 = v1 * (k * 0.45 + 0.05);
                //System.out.println("  Speed: " + v1 + " +/- " + v2);
                // apply push at -π/8 to π/8 of current angle
                double d1 = Math.atan2(v.y, v.x);
                double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
                //System.out.println("  Angle: " + d1 + " -> " + d2);
                // compute energy
                double E = 0.5 * asteroids[i].mass * v2 * v2;
                // try to push asteroid
                Asteroid a1 = null;
                try {
                    a1 = Asteroid.push(asteroids[i], time, E, d2);
                } catch (InvalidOrbitException e) {
                    System.out.println("  Invalid orbit: " + e.getMessage());
                    continue;
                }
                // avoid allocating a new Point object for every position
                Point p1 = v, p2 = new Point();
                // search for collision with other asteroids
                if (iteration == 1) {
                    for (int j = 0; j != asteroids.length; ++j) {
                        if (i == j) continue;
                        Asteroid a2 = asteroids[j];
                        double r = a1.radius() + a2.radius();
                        // look 10 years in the future for collision
                        boolean willCollide = willCollide(a1, a2, 3650);
                        boolean willCollideOrigin = willCollideOrigin(a1, a2, energy, direction, i, E, d2);
                        if(willCollide != willCollideOrigin){
                            debug("Wrong !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                        }
                    }
                    System.out.println("  No collision ...");
                } else {
                    int j = getHeaviestAsteroid(asteroids);
                    // search for collision with other asteroids
                    Asteroid a2 = asteroids[j];
                    double r = a1.radius() + a2.radius();

                    // look 10 years in the future for collision
                    boolean willCollide = willCollide(a1, a2, 3650);

                    boolean willCollideOrigin = willCollideOrigin(a1, a2, energy, direction, i, E, d2);
                    if(willCollide != willCollideOrigin){
                        debug("Wrong !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    }
                }
            }
        }
        time_of_push = time + turns_per_retry;
    }

    public boolean willCollideOrigin(Asteroid a1, Asteroid a2,  double[] energy, double[] direction, int i, double E, double d2){
        double r = a1.radius() + a2.radius();
        Point p1 = new Point();
        Point p2 = new Point();
        for (long ft = 0; ft != 3650; ++ft) {
            long t = time + ft;
            if (t >= time_limit) break;
            a1.orbit.positionAt(t - a1.epoch, p1);
            a2.orbit.positionAt(t - a2.epoch, p2);
            // if collision, return push to the simulator

            if (Point.distance(p1, p2) < r) {
                energy[i] = E;
                direction[i] = d2;
                // do not push again until collision happens
                time_of_push = t + 1;
                System.out.println("  Collision prediction !");
                System.out.println("  Year: " + (1 + t / 365));
                System.out.println("  Day: " + (1 + t % 365));
                debug("current time: " + time);
                debug("origin time: " + t + " || ");
                //System.out.println("will overlap: " + willOverlap(p1, a1.radius(), p2, a2.radius()));                         
                return true;
            }

        }
        return false;
    }

    private int getClosestPairAsteroid(Asteroid[] asteroids) {
        int index1 = 0;
        int index2 = 0;
        double minDistance = Integer.MAX_VALUE;
        Point p1 = new Point();
        Point p2 = new Point();
        for(int aa = 0; aa < asteroids.length; aa++) {
            Asteroid a1 = asteroids[aa];
            a1.orbit.positionAt(time - a1.epoch, p1);
            for(int bb = 0; bb < asteroids.length; bb++) {
                if(aa == bb) {
                    continue;
                }
                Asteroid a2 = asteroids[bb];
                a2.orbit.positionAt(time - a2.epoch, p2);
                double distance = Point.distance(p1, p2);
                if (distance < minDistance) {
                    minDistance = distance;
                    index1 = aa;
                    index2 = aa;
                }
            }
        }
        asteroids[index1].orbit.positionAt(time, p1);
        asteroids[index2].orbit.positionAt(time, p1);
        if(Point.distance(p1, origin) > Point.distance(p2, origin)) {
            return index1;
        }
        return index2;
    }

    public int getFarthestAsteroid(Asteroid[] asteroids) {
        int index = 0;
        double maxDistance = 0;
        Point point = new Point();
        for(int aa = 0; aa < asteroids.length; aa++) {
            Asteroid a1 = asteroids[aa];
            a1.orbit.positionAt(time - a1.epoch, point);
            double distance = Point.distance(point, origin);
            if(distance > maxDistance) {
                maxDistance = distance;
                index = aa;
            }
        }
        return index;
    }

    public boolean willCollide(Asteroid a1, Asteroid a2, long timeInterval){
        //Make sure a1 is always the one has bigger period to short the loop time
        if(a1.orbit.period() > a2.orbit.period()){
            Asteroid tmp = a1;
            a1 = a2;
            a2 = tmp;
        }
        long threshold = 10;
        ArrayList<Long> cts = getCollisonPoints(a1, a2);
        for(long l : cts){
            //System.out.println(l + " meeting time final result");
        }
        Point p1 = new Point();
        Point p2 = new Point();
        a1.orbit.positionAt(time - a1.epoch, p1);
        a2.orbit.positionAt(time - a2.epoch, p2);
        if(willOverlap(p1, a1.radius(), p2, a2.radius())){
            debug("overlap ***********");
            return true;
        }

        for(int i = 0; i < cts.size(); i++){
            long startTime = cts.get(i) - threshold;
            long endTime = cts.get(i) + threshold;

            while(startTime < timeInterval + time){
                //debug("checked Time: " + startTime + " to " + endTime);
                endTime = endTime > timeInterval + time? timeInterval + time: endTime;
                for(int t = (int)startTime; t <= endTime; t++){
                    a1.orbit.positionAt(t - a1.epoch, p1);
                    a2.orbit.positionAt(t - a2.epoch, p2);
                    if(willOverlap(p1, a1.radius(), p2, a2.radius()))
                        return true;
                }
                startTime = startTime + a1.orbit.period();
                endTime = endTime + a1.orbit.period();
            }
        }
        return false;

    }

    public boolean willOverlap(Point p1, double r1, Point p2, double r2){
        double distance = Point.distance(p1, p2);
        if(distance < (r1 + r2))
            return true;
        //debug("will not overlap");
        return false;
    }

    public ArrayList<Long> getCollisonPoints(Asteroid a1, Asteroid a2){

        double period = a1.orbit.period();
        Point tmp = new Point();
        Point center = a2.orbit.center();
        double a = a2.orbit.a;
        double threshold = 2 * a1.radius();
        Point foci1 = new Point(0, 0);
        Point foci2 = new Point(2 * center.x, 2 * center.y);
        ArrayList<Long> cts = new ArrayList<Long>();
        for(int i = 0; i < period; i++){
            a1.orbit.positionAt(time + i - a1.epoch, tmp);
            double distance =Point.distance(foci1, tmp) + Point.distance(foci2, tmp);
            if(Math.abs(2 * a - distance) <= threshold){
                Long ct = time + i;
                cts.add(ct);
            }
        }
        for(long l : cts){
            //System.out.println(l + " meeting time raw result");
        }
        ArrayList<Long> ret = new ArrayList<Long>();
        if(cts.size() > 0){
            ret.add(cts.get(0));
            for(int i = 1; i < cts.size(); i++){
                if(cts.get(i) - cts.get(i - 1) > 1){
                    ret.add(cts.get(i));
                }
            }
        }

        return ret;
    }

    public void debug(String str){
        System.out.print("debug: " + str + "\n");
    }

    private int getLightestAsteroid(Asteroid asteroids[]) {
        int min = 1;
        for (int i = 0; i < asteroids.length; i++) {
            if (asteroids[i].mass < asteroids[min].mass)
                min = i;
        }
        return min;
    }

    private int getHeaviestAsteroid(Asteroid asteroids[]) {
        int max = 1;
        for (int i = 0; i < asteroids.length; i++) {
            if (asteroids[i].mass > asteroids[max].mass)
                max = i;
        }
        return max;
    }

}

