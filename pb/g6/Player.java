package pb.g6;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.ArrayList;
import java.util.Random;

public class Player implements pb.sim.Player {

    // used to pick asteroid and velocity boost randomly
    private Random random = new Random();

    // current time, time limit
    private long time = -1;
    private long time_limit = -1;

    // time until next push
    private long time_of_push = 0;
    private int magnet_index = -1;

    private double point_to_distance(Point p) {return Math.sqrt(p.x*p.x+p.y*p.y);}


    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit) {
        if (Orbit.dt() != 24 * 60 * 60)
            throw new IllegalStateException("Time quantum is not a day");
        this.time_limit = time_limit;
        magnet_index = asteroids.length/2;
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids,
                     double[] energy, double[] direction) {
        

        // if not yet time to push do nothing
        if (++time <= time_of_push) return;
        System.out.println("Year: " + (1 + time / 365));
        System.out.println("Day: " + (1 + time % 365)); 

        double min_E = Double.MAX_VALUE;
        double min_dir = Double.MAX_VALUE;
        int min_index = -1;
        long min_time = 0;

        Asteroid a2 = asteroids[magnet_index];
        Point magnet_position = new Point();
        a2.orbit.positionAt(time-a2.epoch, magnet_position);
        double r2 = point_to_distance(magnet_position);
        double omega2 = Math.sqrt(Orbit.GM / Math.pow(r2,3));


        Point v2 = a2.orbit.velocityAt(time-a2.epoch);
        double d2 = Math.atan2(v2.y,v2.x); 
        
        for (int i=asteroids.length-1; i>=0; i--){

            if (i==magnet_index) continue;

            Asteroid a1 = asteroids[i];
            double r1 = a1.radius();
            
            for (long ft = 0; ft < 365 * 10; ++ft) {
                long t = time + ft;
                if (t >= time_limit) break;

                Point v1 = a1.orbit.velocityAt(t-a1.epoch);
                double d1 = Math.atan2(v1.y,v1.x);
                double tH = Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM));
                double thresh = r1 + r2;
                v2 = a2.orbit.velocityAt(t-a2.epoch);
                d2 = Math.atan2(v2.y, v2.x);

                if ( Math.abs(d1 + Math.PI - d2 - tH*omega2) < thresh / r2) {
                    double v_new = Math.sqrt(Orbit.GM / r1) * (Math.sqrt(2*r2/(r1+r2)) - 1);
                    double calculated_E = 0.5*a1.mass * v_new * v_new; 

                    try {
                        Asteroid a1_pushed = Asteroid.push(a1, t, calculated_E, d1);

                        if(calculated_E < min_E){
                            System.out.println("min_index: " + min_index);
                            System.out.println("min_E: " + min_E);
                            System.out.println(" ===========Min_E updated ============= ");
                            min_E = calculated_E;
                            min_dir = d1;
                            min_index = i;
                            min_time = ft;
                        }
                    } catch (InvalidOrbitException e) {
                        System.out.println("  Invalid orbit: " + e.getMessage());
                    }
                 }
             }
             
            
        }
        if(min_index != -1){
            energy[min_index] = min_E;
            direction[min_index] = min_dir;
            // do not push again until collision happens
            time_of_push = min_time;
            System.out.println("  Collision prediction !");
            System.out.println("  Year: " + (1 + time_of_push / 365));
            System.out.println("  Day: " + (1 + time_of_push % 365));

        }

    }
}
