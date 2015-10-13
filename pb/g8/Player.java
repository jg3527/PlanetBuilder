package pb.g8;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;
import pb.sim.Orbit;
import pb.sim.Point;

public class Player implements pb.sim.Player {

    //iteration number
    int iteration = 1;

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
    private int total_number;
    private HashMap<Integer, List<Long>> asteroidClusters;
    private int cluster_number = 0;
    //key is the id of ast, value is the index of ast
    private HashMap<Long, Asteroid> indexMap;

    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit) 
    {
        if (Orbit.dt() != 24 * 60 * 60)
            throw new IllegalStateException("Time quantum is not a day");
        this.time_limit = time_limit;
        this.total_number = asteroids.length;
        refreshIndexMap(asteroids);
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids, double[] energy, double[] direction) 
    {
    	refreshIndexMap(asteroids);
        // if not yet time to push do nothing
        if (++time <= time_of_push) return;
        for(Asteroid a: asteroids) {
            System.out.println(a.id);
        }
        long time_left_per_asteroid = (time_limit - time)/asteroids.length;
        time_left_per_asteroid = Math.max(time_left_per_asteroid, 3650);
        //System.out.println("Year: " + (1 + time / 365));
        //System.out.println("Day: "  + (1 + time % 365));
        List<Push> pushes = new ArrayList<Push>();
        for (int retry = 1 ; retry <= retries_per_turn ; ++retry) 
        {
            // pick a random asteroid and get its velocity
            int i = 0;
            if(iteration == 1) 
            {
                i = getFarthestAsteroid(asteroids);
            }
            else
            {
                i = getHeaviestAsteroid(asteroids);
            }

            // i = getClosestPairAsteroid(asteroids);
            
            Point v = asteroids[i].orbit.velocityAt(time);
            
            // add 5-50% of current velocity in magnitude
            //System.out.println("Try: " + retry + " / " + retries_per_turn);
            
            double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
            double v2 = 0.0;
            for (double k=0; k< 1; k=k+0.001) 
            {
                v2 = v1 * (k * 0.45 + 0.05); //Can't change this
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
                    //System.out.println("  Invalid orbit: " + e.getMessage());
                    continue;
                }
                // avoid allocating a new Point object for every position
                Point p1 = v, p2 = new Point();
            
                // search for collision with other asteroids
                if (iteration == 1) 
                {
                    for (int j = 0; j != asteroids.length; ++j) 
                    {
                        if (i == j) continue;
                        Asteroid a2 = asteroids[j];
                        
                        // look in the future for collision         
                        boolean willCollide = willCollide(a1, a2, time_left_per_asteroid, p1, p2);
                        
                        if (willCollide) 
                        {
                        	debug("COllide !!" + willCollide);
                           pushes.add(new Push(i, E, d2, time_of_push));
//                           energy[i] = E;
//                           direction[i] = d2;
                           // do not push again until collision happens
                           long t = time_of_push - 1;
                           System.out.println("  Collision prediction !");
                           System.out.println("  Year: " + (1 + t / 365));
                           System.out.println("  Day: " + (1 + t % 365));
                           debug("current time: " + time);
                           debug("origin time: " + t + " || ");
                            //iteration number breaks here
                           iteration++;
                           //System.out.println("will overlap: " + willOverlap(p1, a1.radius(), p2, a2.radius()));                         
                           break;
                       }
                        
                    }
                    //System.out.println("  No collision ...");
                } 
                
                else 
                {
                    // search for collision with other asteroids
                    Asteroid a2 = asteroids[i];

                    // look 10 years in the future for collision
                    //boolean willCollideOrigin = willCollideOrigin(a1, a2, energy, direction, i, E, d2);
                    boolean willCollide = willCollide(a1, a2, time_left_per_asteroid, p1, p2);
                  
                    if(willCollide)
                    {
                    	debug("COllide " + willCollide);
                        energy[i] = E;
                        direction[i] = d2;
                        // do not push again until collision happens
                        long t = time_of_push - 1;
                        System.out.println("  Collision prediction !");
                        System.out.println("  Year: " + (1 + t / 365));
                        System.out.println("  Day: " + (1 + t % 365));
                        //iteration number breaks here
                        iteration++;
                        //System.out.println("will overlap: " + willOverlap(p1, a1.radius(), p2, a2.radius()));                         
                        return;
                    
                    }
                    	
                }
            }
        }
        if(pushes.size() > 0) 
        {
            Push min_push = pushes.get(0);
            for(Push push: pushes) 
            {
                System.out.println("Energy:" + push);
                if(push.energy < min_push.energy) 
                {
                    min_push = push;
                }
            }
            System.out.println("Min Energy:" + min_push);
            energy[min_push.asteroid_id] = min_push.energy;
            direction[min_push.asteroid_id] = min_push.direction;
            time_of_push = min_push.time_of_push;
            return;
        }
        time_of_push = time + turns_per_retry;
    }

    public int getFarthestAsteroid(Asteroid[] asteroids) 
    {
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

    public boolean willCollide(Asteroid a11, Asteroid a22, long timeInterval,Point p1, Point p2)
    {    	
    	Asteroid a1, a2;
        //Make sure a1 is always the one has bigger period to short the loop time
        if(a11.orbit.period() > a22.orbit.period())
        {
            a1 = a22;
            a2 = a11;
        }
        else
        {
        	a1 = a11;
        	a2 = a22;
        }
        long threshold = 20;
        ArrayList<Long> cts = getCollisonPoints(a1, a2);
       
        a1.orbit.positionAt(time - a1.epoch, p1);
        a2.orbit.positionAt(time - a2.epoch, p2);
        if(willOverlap(p1, a1.radius(), p2, a2.radius()))
        {
            debug("overlap ***********");
            time_of_push = time + 2;
            return true;
        }

        for(int i = 0; i < cts.size(); i++)
        {
            long startTime = cts.get(i) - threshold;
            long endTime = cts.get(i) + threshold;

            while(startTime < timeInterval + time){
                //debug("checked Time: " + startTime + " to " + endTime);
                endTime = endTime > timeInterval + time? timeInterval + time: endTime;
                for(int t = (int)startTime; t <= endTime; t++){
                    a1.orbit.positionAt(t - a1.epoch, p1);
                    a2.orbit.positionAt(t - a2.epoch, p2);
                    if(willOverlap(p1, a1.radius(), p2, a2.radius()))
                    {
                    	debug("will overlap");
	                    	time_of_push = t + 1;
	                        return true;
                    }
                }
                startTime = startTime + a1.orbit.period();
                endTime = endTime + a1.orbit.period();
            }
        }
        return false;

    }

    public boolean willOverlap(Point p1, double r1, Point p2, double r2)
    {
        double distance = Point.distance(p1, p2);
        if(distance < (r1 + r2))
            return true;
        //debug("will not overlap");
        return false;
    }

    public ArrayList<Long> getCollisonPoints(Asteroid a1, Asteroid a2)
    {
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
            if(Math.abs(2 * a - distance) <= threshold)
            {
                Long ct = time + i;
                cts.add(ct);
            }
        }
        ArrayList<Long> ret = new ArrayList<Long>();
        if(cts.size() > 0){
            ret.add(cts.get(0));
            for(int i = 1; i < cts.size(); i++)
            {
                if(cts.get(i) - cts.get(i - 1) > 1)
                {
                    ret.add(cts.get(i));
                }
            }
        }

        return ret;
    }

    public void debug(String str)
    {
        //System.out.print("debug: " + str + "\n");
    }

    private int getLightestAsteroid(Asteroid asteroids[]) 
    {
        int min = 1;
        for (int i = 0; i < asteroids.length; i++)
        {
            if (asteroids[i].mass < asteroids[min].mass)
                min = i;
        }
        return min;
    }

    private int getHeaviestAsteroid(Asteroid asteroids[])
    {
        int max = 1;
        for (int i = 0; i < asteroids.length; i++)
        {
            if (asteroids[i].mass > asteroids[max].mass)
                max = i;
        }
        return max;
    }
    private List<Long> getOutersidePart(Asteroid[] asteroids){
    	ArrayList<Long> ret = new ArrayList<Long>();
    	Arrays.sort(asteroids, new Comparator<Asteroid>() {
    	    public int compare(Asteroid a1, Asteroid a2) {
    	        return (int)(a1.orbit.a - a2.orbit.a);
    	    }
    	});
    	double sum = 0;
    	for(int i = 0; i < asteroids.length; i++){
    		sum += asteroids[i].mass;
    	}
    	double half = 0;
    	for(int i = asteroids.length - 1; i >= 0; i++){
    		if(half >= sum / 2){
    			return ret;
    		}
    		half += asteroids[i].mass;
    		ret.add(asteroids[i].id);		
    	}
    	return ret;	
    }
    private double numberOfClusters(Asteroid[] asteroids)
    {
        double clusters = asteroids.length;
        for (Asteroid a1: asteroids) {
            double nearbyAsteroids = 0;
            for(Asteroid a2: asteroids) {
                if(a1.id != a2.id) {
                    //distance between them less than some threshold
                    //if(Math.abs(a1.orbit.a - a2.orbit.a) < ) {
                    	nearbyAsteroids++;
                    //}
                }
            }
            clusters -= nearbyAsteroids / (nearbyAsteroids + 1);
        }
        return clusters;
    }
    private double tryToCollideOutside(Asteroid[] asteroids){
    	
    	List<Long> ids = new ArrayList<Long>();

    	Point origin = new Point(0, 0);
    	for(int i = 0; i < cluster_number; i++){
    		ids = asteroidClusters.get(i);
    		Collections.sort(ids, new Comparator<Long>() {
				@Override
				public int compare(Long l1, Long l2) {
					Asteroid a1 = indexMap.get(l1);
					Asteroid a2 = indexMap.get(l2);
					double d1 = Point.distance(origin, a1.orbit.positionAt(time - a1.epoch));
					double d2 = Point.distance(origin, a2.orbit.positionAt(time - a2.epoch));
					return (int)(d1 - d2);
				}
			});
    		if(ids.size() == 0)
    			continue;
    		int size = ids.size();
    		//loop within this cluster
    		if(size == 1){
    			//TODO push it to the next cluster
    		}else{
	    		for(int j = 1; j < size; j++){
	    			Asteroid a1 = indexMap.get(ids.get(0));
	    		}
    		}
    		
    	}
    }
    private void refreshIndexMap(Asteroid[] asteroids){
    	indexMap = new HashMap<Long, Asteroid>();
    	for(int i = 0; i < asteroids.length; i++){
    		indexMap.put(asteroids[i].id, asteroids[i]);
    	}
    }
    private void checkCollisons(Asteroid a1, List<Asteroid> asteroids){
    	Point v = a1.orbit.velocityAt(time);
        
    	double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
        double v2 = 0.0;
    }

}

