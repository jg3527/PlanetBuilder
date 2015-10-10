package pb.g6;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class Player implements pb.sim.Player {

    // used to pick asteroid and velocity boost randomly
    private Random random = new Random();

    // current time, time limit
    private long time = -1;
    private long time_limit = -1;

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
//		for (int i=0; i < asteroids.length; i++){
//			System.out.println(asteroids[i].mass);
//			System.out.println(asteroids[i].radius());
//		}
//		System.exit(0);
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids, double[] energy, double[] direction) {
        // if not yet time to push do nothing
        if (++time <= time_of_push) return;
//        System.out.println("Year: " + (1 + time / 365));
//        System.out.println("Day: " + (1 + time % 365));
        for (int retry = 1; retry <= retries_per_turn; ++retry) {

            // pick an asteroid with the smallest mass
//			double minMass = asteroids[0].mass;
//			int minMassIndex = 0;
//			for (int i = 0; i < asteroids.length; i++){
//				if (asteroids[i].mass < minMass){
//					minMass = asteroids[i].mass;
//					minMassIndex = i;
//				}
//			}
//			int i = minMassIndex;

            ArrayList<CandidatePush> cands = new ArrayList<CandidatePush>();
            for(int i = 0; i < asteroids.length; ++i) {


                // pick a random asteroid and get its velocity
//                int i = random.nextInt(asteroids.length);

                Point v = asteroids[i].orbit.velocityAt(time);
                // add 5-50% of current velocity in magnitude
//                System.out.println("Try: " + retry + " / " + retries_per_turn);
                double v1 = Math.sqrt(v.x * v.x + v.y * v.y);

                //Changing velocity addition from 0.05 ~ 0.5 to 0.05 ~ 0.2
                double v2 = v1 * (random.nextDouble() * 0.15 + 0.05);
                //Random deviation from energy
//                System.out.println("  Speed: " + v1 + " +/- " + v2);
                // apply push at -π/8 to π/8 of current angle
                double d1 = Math.atan2(v.y, v.x);
                //changing angle select from 0.25 to 0.35
                double d2 = d1 + (random.nextDouble() - 0.7) * Math.PI * 0.35;
//                System.out.println("  Angle: " + d1 + " -> " + d2);
                // compute energy
                double E = 0.5 * asteroids[i].mass * v2 * v2;
                // try to push asteroid
                Asteroid a1;
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
                            cands.add(new CandidatePush(i, t, E, d2));
//                            energy[i] = E;
//                            direction[i] = d2;
                            // do not push again until collision happens
//                            time_of_push = t + 1;
//                            System.out.println("  Collision prediction !");
//                            System.out.println("  Year: " + (1 + t / 365));
//                            System.out.println("  Day: " + (1 + t % 365));
                        }
                    }
                }
            }

            if (cands.size() != 0 ) {
                Collections.sort(cands);
                CandidatePush chosen = cands.get(0);
                energy[chosen.ast_idx] = chosen.energy;
                direction[chosen.ast_idx] = chosen.angle;
                time_of_push = chosen.time + 1;

//                System.out.println("  Collision prediction !");
//                System.out.println("  Year: " + (1 + chosen.time / 365));
//                System.out.println("  Day: " + (1 + chosen.time % 365));

                return;
            }

//            System.out.println("  No collision ...");
        }
        time_of_push = time + turns_per_retry;
    }
}
