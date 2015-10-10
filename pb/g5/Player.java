package pb.g5;

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
	
	private boolean collisionImminent;
	private long nextCollision;
	private int n_asteroids;

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit)
	{
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
		collisionImminent = false;
		nextCollision = Long.MIN_VALUE;
		n_asteroids = asteroids.length;
	}
	
	public Push findCollision(Asteroid a1, Asteroid a2, Point collisionPoint, long collisionTime) {
		Point v1 = a1.orbit.velocityAt(time - a1.epoch);
		Point currentLocation = a1.orbit.positionAt(time - a1.epoch);
		double speed = Math.hypot(v1.x, v1.y);
		double energy = 0.5 * a1.mass * speed * speed;
		double[] energies = new double[50];
		double[] directions = new double[60];
		double d1 = Math.atan2(currentLocation.y - collisionPoint.y, currentLocation.x - collisionPoint.x);
		for(int i = 0; i < directions.length; ++i) {
			directions[i] = d1 + Math.toRadians(30-i);
		}
		for(int i = 0; i < energies.length; ++i) {
			double pushSpeed = speed*(i+1)/100;

			energies[i] = pushSpeed*pushSpeed*0.5*a1.mass;
			if(energies[i] <0){
			    System.out.println("pushedSpeed : " + pushSpeed + " , mass " + a1.mass);
			}
			for(Double direction : directions) {
				Asteroid pushed;
				try {
					pushed = Asteroid.push(a1, time, energies[i], direction);
				} catch (NumberFormatException  | InvalidOrbitException e) {
					continue;
				}
				Point locationAtCollisionTime = pushed.orbit.positionAt(collisionTime+time - pushed.epoch);
				double sumRadii = a1.radius() + a2.radius();
				if(Point.distance(collisionPoint, locationAtCollisionTime) <= sumRadii) {
                    System.out.println(currentLocation);
                    System.out.println(sumRadii);
                    System.out.println(Point.distance(collisionPoint, locationAtCollisionTime));
                    return new Push(energies[i], direction);
				}
			}
		}
		return null;
	}

	// try to push asteroid
	public void play(Asteroid[] asteroids,
	                 double[] energy, double[] direction)
	{
		++time;
		if(time % 365 == 0){
        System.out.println("Year: " + (1 + time / 365));
		System.out.println("Day: "  + (1 + time % 365));}
		if(time < 365) {
			return;
		}
		// if not yet time to push do nothing
		if (time < nextCollision) {
			return;
		} else if(time == nextCollision) {
			collisionImminent = false;
		}
		
		if(asteroids.length <= 1)
			return;


		Asteroid a1 = asteroids[0];
		Asteroid a2 = asteroids[1];

        
        Push optimalPush = null;
        double lowestEnergy = Double.MAX_VALUE;
        Point bestCollisionPoint = null;
        long bestCollisionTime = 0;
        long[] collisionTimes = new long[100];
        for(int i = 0; i < collisionTimes.length; ++i) {
        	collisionTimes[i] = 500+i*10;
        	Point collisionPoint = a2.orbit.positionAt(time + collisionTimes[i] - a2.epoch);
            Push push = findCollision(a1, a2, collisionPoint, collisionTimes[i]);
            if(push != null && push.energy < lowestEnergy) {
            	optimalPush = push;
            	lowestEnergy = push.energy;
            	bestCollisionPoint = collisionPoint;
            	bestCollisionTime = collisionTimes[i];
            }
        }
        

        if(optimalPush == null) {
			//System.out.println("failed");
		} else {
			Push push = optimalPush;
			System.out.println("curLoc: " + a1.orbit.positionAt(time));
			System.out.println("collisionPoint: " + bestCollisionPoint);
			System.out.println("Push: "+push);
            System.out.println("Year: " + (1 + time / 365));
            System.out.println("Day: "  + (1 + time % 365));
			energy[0] = push.energy;
			direction[0] = push.direction;
			nextCollision = bestCollisionTime + time;
			collisionImminent = true;
		}

        //System.out.println(Arrays.toString(energy));
/*		for (int retry = 1 ; retry <= retries_per_turn ; ++retry) {
			// pick a random asteroid and get its velocity
			int i = random.nextInt(asteroids.length);
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
			a1 = null;
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
				a2 = asteroids[j];
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
		}*/
		time_of_push = time + turns_per_retry;
	}
}
