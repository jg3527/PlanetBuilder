package pb.g8;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Arrays;
import java.util.Random;

public class Player implements pb.sim.Player {

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
    public void init(Asteroid[] asteroids, long time_limit)
    {
        if (Orbit.dt() != 24 * 60 * 60)
            throw new IllegalStateException("Time quantum is not a day");
        this.time_limit = time_limit;
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids,
                     double[] energy, double[] direction)
    {
        // if not yet time to push do nothing
//        System.out.println(time + Arrays.toString(asteroids));
        if (++time <= time_of_push) return;
        System.out.println("Year: " + (1 + time / 365));
        System.out.println("Day: "  + (1 + time % 365));
        for (int retry = 1 ; retry <= retries_per_turn ; ++retry) {
            // pick a random asteroid and get its velocity
            int i = 0;
            i = getFarthestAsteroid(asteroids);
//            i = getClosestPairAsteroid(asteroids);
            Point v = asteroids[i].orbit.velocityAt(time);
            // add 5-50% of current velocity in magnitude
            System.out.println("Try: " + retry + " / " + retries_per_turn);
            double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
            double v2 = 0.0;
            for(int aa = 1; aa < 10; aa++) {
                v2 = v1 * (0.05 * aa);
                System.out.println("  Speed: " + v1 + " +/- " + v2);
                // apply push at -π/8 to π/8 of current angle
                double d1 = Math.atan2(v.y, v.x);
                double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
                System.out.println("  Angle: " + d1 + " -> " + d2);
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
                for (int j = 0; j != asteroids.length; ++j) {
                    if (i == j) continue;
                    Asteroid a2 = asteroids[j];
                    double r = a1.radius() + a2.radius();
                    // look 10 years in the future for collision
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
                            return;
                        }
                    }
                }
            }
            System.out.println("  No collision ...");
        }
        time_of_push = time + turns_per_retry;
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
}
