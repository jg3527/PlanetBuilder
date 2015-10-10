package pb.g1;

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
		if (++time <= time_of_push) return;
		System.out.println("Year: " + (1 + time / 365));
		System.out.println("Day: "  + (1 + time % 365));
		for (int retry = 1 ; retry <= retries_per_turn ; ++retry) {
			// pick a random asteroid and get its velocity
			//int i = random.nextInt(asteroids.length);
			
			int i= getLightestAsteroid(asteroids);
			Point v = asteroids[i].orbit.velocityAt(time);
			// add 5-50% of current velocity in magnitude
			System.out.println("Try: " + retry + " / " + retries_per_turn);
			double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
			double v2 = v1 * (random.nextDouble() * 0.45 + 0.05);
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
			for (int j = 0 ; j != asteroids.length ; ++j) {
				if (i == j) continue;
				Asteroid a2 = asteroids[j];
				double r = a1.radius() + a2.radius();
				// look 10 years in the future for collision
				for (long ft = 0 ; ft != 3650 ; ++ft) {
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
						System.out.println("  Day: "  + (1 + t % 365));
						return;
					}
				}
			}
			System.out.println("  No collision ...");
		}
		time_of_push = time + turns_per_retry;
	}
	
	public static int getLightestAsteroid(Asteroid[] asteroids){
		double min_weight = 0;
		int min_weight_index = 0;
		for(int i=0;i<asteroids.length;i++){
			if (asteroids[i].mass< min_weight){
				min_weight = asteroids[i].mass;
				min_weight_index = i;
			}
		}
		
		return min_weight_index;
		
	}
	
	public static int getHeaviestAsteroid(Asteroid[] asteroids){
		double max_weight = 0;
		int max_weight_index = 0;
		for(int i=0;i<asteroids.length;i++){
			if (asteroids[i].mass> max_weight){
				max_weight = asteroids[i].mass;
				max_weight_index = i;
			}
		}
		
		return max_weight_index;
		
	}
	
	
	public static int getFarthestAsteroid(Asteroid[] asteroids){
		System.out.println("Method Called");
		double max_major_radius = 0;
		int max_major_radius_index = 0;
		for(int i=0;i<asteroids.length;i++){
			if (asteroids[i].orbit.a > max_major_radius){
				System.out.println("Asteroid "+i+" asteroids[i].orbit.a ");
				max_major_radius = asteroids[i].orbit.a ;
				max_major_radius_index = i;
			}
		}
		
		return max_major_radius_index;
		
	}
}
