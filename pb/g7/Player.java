package pb.g7;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

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
	    time++;
	    // if not yet time to push do nothing
	    /*	    System.out.println("Year: " + (1 + time / 365));
		    System.out.println("Day: "  + (1 + time % 365));*/
	    int inner = -1;
	    int outer = -1;
	    double innerR = 0;
	    //boolean stabilized = false;
	    if (inner == -1) {
		if (asteroids[0].orbit.a < asteroids[1].orbit.a) {
		    inner = 0;
		    outer = 1;
		    innerR = asteroids[inner].orbit.a;
		}
		else {
		    inner = 1;
		    outer = 0;
		    innerR = asteroids[inner].orbit.a;
		}
	    }
	    double r1 = asteroids[inner].orbit.a;	    
	    double r2 = asteroids[outer].orbit.a;	    
	    double omega_outer = Math.sqrt(Orbit.GM / Math.pow(r2,3));
	    Point v_inner = asteroids[inner].orbit.velocityAt(time);
	    double normv = Math.sqrt(v_inner.x*v_inner.x+v_inner.y*v_inner.y);
	    double theta_inner = Math.atan2(v_inner.y,v_inner.x);
	    Point v_outer = asteroids[outer].orbit.velocityAt(time);
	    double theta_outer = Math.atan2(v_outer.y,v_outer.x);
	    double tH = Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM));
	    double thresh = asteroids[inner].radius() + asteroids[outer].radius();
	    /*Point asteroidLocation = new Point();
	      asteroids[inner].orbit.positionAt(time, asteroidLocation);*/
	    if ( Math.abs(theta_inner + Math.PI - theta_outer - tH*omega_outer) < thresh / r2) {
		double vnew = Math.sqrt(Orbit.GM / r1) * (Math.sqrt( 2*r2 / (r1+r2)) - 1);
		energy[inner] = 0.5*asteroids[inner].mass * vnew * vnew;
		direction[inner] = theta_inner;
	    }
	    /*	    else if ( !stabilized && Math.abs(Math.sqrt(asteroidLocation.x*asteroidLocation.x + asteroidLocation.y*asteroidLocation.y) - r2) < 0.001 * (asteroids[inner].radius() + asteroids[outer].radius())) {		
		double r1 = innerR;
		double dist = Math.sqrt(asteroidLocation.x*asteroidLocation.x + asteroidLocation.y*asteroidLocation.y);
		Point next = new Point();
		asteroids[inner].orbit.positionAt(time+1,next);
		if (Math.abs(Math.sqrt(next.x*next.x-next.y*next.y)-r2) > Math.abs(dist-r2)) {
		    double vnew = Math.sqrt(Orbit.GM / r2) * (1 - Math.sqrt(2*r1/(r1+r2)));
		    energy[inner] = 0.5*asteroids[inner].mass * vnew * vnew;
		    direction[inner] = theta;
		    stabilized = true;
		}
		}*/
	}
}
