package pb.g8;


import java.util.*;


import net.sf.javaml.clustering.Clusterer;
import net.sf.javaml.clustering.KMeans;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.DenseInstance;
import net.sf.javaml.core.Instance;
import pb.sim.Asteroid;
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
    private HashMap<Integer, Long> time_of_push;

    // number of retries
    private int retries_per_turn = 1;
    private int turns_per_retry = 3;
    private int total_number;

    private HashMap<Integer, List<Long>> asteroidClusters;
    private int cluster_number = 0;
    
    private HashMap<Long, Asteroid> indexMap;
  //key is the id of ast, value is the index of ast
    private HashMap<Long, Integer> indexHashMap; 
    private double clusterThreshold;
    private HashMap<Integer, Boolean> pushAgain;

    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit) 
    {
        if (Orbit.dt() != 24 * 60 * 60)
            throw new IllegalStateException("Time quantum is not a day");
        this.time_limit = time_limit;
        this.total_number = asteroids.length;

        refreshIndexMap(asteroids);
        reorderCluster(asteroids);
        time_of_push = new HashMap<Integer, Long>();
        for(int i=0; i < cluster_number; i++) {
        	time_of_push.put(i, 0l);
        }
        pushAgain = new HashMap<Integer, Boolean>();
        System.out.println("Initialization done!");
    }

    private void reorderCluster(Asteroid[] asteroids) {
        Set<Long> relevantAsteroidIds = getOutersidePart(asteroids);
        int clusterCount = numberOfClusters(relevantAsteroidIds);
        System.err.println(clusterCount);
        asteroidClusters = getAsteroidClusters(relevantAsteroidIds, clusterCount);
        cluster_number = asteroidClusters.size()-1;
    }

   

    private HashMap<Integer, List<Long>> getAsteroidClusters(Set<Long> relevantAsteroidIds, int clusterCount) {
        // gets the clusters for the provided dataset
        Dataset[] clusters = computeClusters(relevantAsteroidIds, clusterCount);
        HashMap<Integer, List<Long>> result = new HashMap<Integer, List<Long>>();
        for(int ii = 0; ii < clusters.length; ii++) {
            List<Long> asteroidIds = new ArrayList<Long>();
            for (Instance instance : clusters[ii]) {
                // gets asteroid associated with the instance
                for(Long asteroidId: relevantAsteroidIds) {
                    Asteroid asteroid = indexMap.get(asteroidId);
                    if (asteroid.orbit.a == instance.value(0)) {
                        asteroidIds.add(asteroid.id);
                    }
                }
            }
            result.put(ii, asteroidIds);
        }
        // add all asteroids that we don't want to consider in cluster -1
        Set<Long> remainingAsteroidIds = indexMap.keySet();
        remainingAsteroidIds.removeAll(relevantAsteroidIds);
        result.put(-1, new ArrayList<Long>(remainingAsteroidIds));
        return result;
    }

    private Set<Long> getAsteroidIds(Asteroid[] asteroids) {
        Set<Long> ids = new HashSet<Long>();
        for(Asteroid asteroid: asteroids) {
            ids.add(asteroid.id);
        }
        return ids;
    }

    private Dataset[] computeClusters(Set<Long> relevantAsteroidIds, int clusterCount) {
        // create default data set
        Dataset data = new DefaultDataset();
        // create instances and populate the data set
        for(Long asteroidId: relevantAsteroidIds) {
            Asteroid asteroid = indexMap.get(asteroidId);
            Instance instance = new DenseInstance(new double[] {asteroid.orbit.a});
            data.add(instance);
        }
        /* Create a new instance of the KMeans algorithm, with no options
         * specified. By default this will generate 4 clusters.
         */
        Clusterer km = new KMeans(clusterCount);
        /*  Cluster the data, it will be returned as an array of data sets,
         *  with each dataset representing a cluster.
         */
        return km.cluster(data);
    }


    private int numberOfClusters(Set<Long> relevantAsteroidIds)
    {
        double clusters = relevantAsteroidIds.size();
        for (Long id1: relevantAsteroidIds) {
            double nearbyAsteroids = 0;
            Asteroid a1 = indexMap.get(id1);
            for (Long id2: relevantAsteroidIds) {
                if (id1.longValue() != id2.longValue()) {
                    Asteroid a2 = indexMap.get(id2);
                    if(Math.abs(a1.orbit.a - a2.orbit.a) < clusterThreshold) {
                        nearbyAsteroids++;
                    }
                }
            }
            clusters -= nearbyAsteroids / (nearbyAsteroids + 1);
        }
        return (int) clusters;
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids, double[] energy, double[] direction) 
    {
    	time++;
    	refreshIndexMap(asteroids);
        // if not yet time to push do nothing
        updateClusters(asteroids);
        int count = 0;
        for(int i = 0; i < cluster_number; i++){
        	if(time <= time_of_push.get(i)){
        		pushAgain.put(i, false);
        		count++;
        	}else{
        		pushAgain.put(i, true);
        	}
        }
        if(count == cluster_number){
        	return;
        }
        long time_left_per_asteroid = (time_limit - time)/asteroids.length;
        time_left_per_asteroid = Math.max(time_left_per_asteroid, 3650);
        //System.out.println("Year: " + (1 + time / 365));
        //System.out.println("Day: "  + (1 + time % 365));
        List<Push> pushes = new ArrayList<Push>();
        //============================================
        tryToCollideOutside(asteroids, direction, energy);
        
        //============================================
        /*for (int retry = 1 ; retry <= retries_per_turn ; ++retry) 
        {
        	for(Map.Entry<Integer, List<Long>> entry: asteroidClusters.entrySet()) {
        		int clusterId = entry.getKey();
        		List<Long> asteroidIds = entry.getValue();
        		if(clusterId == -1) {
        			continue;
        		}
            		long i = getClosestAsteroid(asteroidIds);
               		Point v = indexMap.get(i).orbit.velocityAt(time);
            		double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
                    double v2 = 0.0;
                    for (double k=0; k< 1; k=k+0.001) 
                    {
                    	v2 = v1 * (k * 0.45 + 0.05); //Can't change this
                        double d1 = Math.atan2(v.y, v.x);
                        double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
                        double E = 0.5 * indexMap.get(i).mass * v2 * v2;
                        Asteroid a1 = null;
                        try {
                            a1 = Asteroid.push(indexMap.get(i), time, E, d2);
                        } catch (InvalidOrbitException e) {
                            continue;
                        }
                        Point p1 = v, p2 = new Point();
                		
                        for(Long j: asteroidIds) {
                         {
                             if (i == j) continue;
                             Asteroid a2 = indexMap.get(j);
                             boolean willCollide = willCollide(a1, a2, time_left_per_asteroid, p1, p2);
                             if (willCollide) 
                             {
                             	debug("COllide !!" + willCollide);
                                pushes.add(new Push((int)i, E, d2, time_of_push));
                                long t = time_of_push - 1;
                                System.out.println("  Collision prediction !");
                                System.out.println("  Year: " + (1 + t / 365));
                                System.out.println("  Day: " + (1 + t % 365));
                                debug("current time: " + time);
                                debug("origin time: " + t + " || ");
                                 break;
                    }
                }
                    }
                    if(pushes.size() > 0) 
                    {
                        Push min_push = pushes.get(0);
                        for(Push push: pushes) 
                        {
                        	System.out.println("Energy: " + push);
                            if(push.energy < min_push.energy) 
                            {
                                min_push = push;
                            }
                        }
                        System.out.println("Min Energy: " + min_push);
                        energy[min_push.asteroid_id] = min_push.energy;
                        direction[min_push.asteroid_id] = min_push.direction;
                        time_of_push = min_push.time_of_push;
                        return;
                    }
                    time_of_push = time + turns_per_retry;
    }
        }
        	}
        }*/
    }

    private void updateClusters(Asteroid[] asteroids) {
        int count = 0;
        for(Map.Entry<Integer, List<Long>> entry: asteroidClusters.entrySet()) {
            count += entry.getValue().size();
        }
        if(count == asteroids.length) {
            return;
        }
        Set<Long> newAsteroidIds = new HashSet<Long>();
        for(Asteroid asteroid: asteroids) {
            newAsteroidIds.add(asteroid.id);
        }
        int clusterId = -1;
        Set<Long> idsToRemove = new HashSet<Long>();
        int removedAsteroidCount = 0;
        long newId = -1;
        // go through all the asteroids present in previous cluster config
        for(Map.Entry<Integer, List<Long>> entry: asteroidClusters.entrySet()) {
            List<Long> ids = entry.getValue();
            // looping on each asteroid for a cluster
            for(Long id: ids) {
                // if the cluster asteroid still exist, remove it from new asteroids
                // else need to remove it from cluster
                if(newAsteroidIds.contains(id)) {
                    newAsteroidIds.remove(id);
                } else {
                    removedAsteroidCount++;
                    // If removed asteroids are greater than 2, time to reorder
                    if(removedAsteroidCount > 2) {
                        reorderCluster(asteroids);
                        return;
                    }
                    // if everything is ok, remove this from cluster
                    clusterId = entry.getKey();
                    idsToRemove.add(id);
                }
            }
        }
        // if more than 1 new asteroid, reorder and return
        if(newAsteroidIds.size() != 1) {
            reorderCluster(asteroids);
            return;
        }
        for(Long id: newAsteroidIds) {
            newId = id;
        }
        // remove the asteroids from the cluster
        List<Long> newAsteroidIdsForCluster = new ArrayList<Long>();
        for(Long id: asteroidClusters.get(clusterId)) {
            if(!idsToRemove.contains(id)) {
                newAsteroidIdsForCluster.add(id);
            }
        }
        asteroidClusters.put(clusterId, newAsteroidIdsForCluster);
        // add new id to the cluster
        asteroidClusters.get(clusterId).add(indexMap.get(newId).id);
        }

    private long getClosestAsteroid(List<Long> asteroidIds) {
    	long min = asteroidIds.get(0);
    	for(Long asteroidId: asteroidIds)
    	{
    		Asteroid a = indexMap.get(asteroidId);
    		if(a.orbit.a < indexMap.get(min).orbit.a)
    			min = a.id;
    	}
    	return min;
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
/*
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
*/
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

    private Set<Long> getOutersidePart(Asteroid[] asteroids){
    	Set<Long> ret = new HashSet<Long>();
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
    	long firstId = asteroids[asteroids.length - 1].id;
    	long secondId = asteroids[asteroids.length - 1].id;
    	for(int i = asteroids.length - 1; i >= 0; i--){
    		if(half >= sum / 2){
    			break;
    		}
    		half += asteroids[i].mass;
    		ret.add(asteroids[i].id);
    		secondId = asteroids[i].id;
    	}
    	clusterThreshold = Math.abs(indexMap.get(firstId).orbit.a - indexMap.get(secondId).orbit.a) / (ret.size() - 1);
    	return ret;	
    }
    
    private void tryToCollideOutside(Asteroid[] asteroids, double[] energy, double[] direction){	
    	List<Long> ids = new ArrayList<Long>();

    	Point origin = new Point(0, 0);
    	for(int i = 0; i < cluster_number; i++){
    		if(!pushAgain.get(i))
    			continue;
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
    		int size = ids.size();
    		if(ids.size() == 0) {
    			System.out.println("continuing because size is zero");
    			continue;
    		}
    		
    		//loop within this cluster
    		if(size == 1){
    			//TODO push it to the next cluster
    		}else{
	    		//TODO 
    			Asteroid a1 = indexMap.get(ids.get(0));
    			for(int j = 1; j < ids.size(); j++){
    				Asteroid a2 = indexMap.get(ids.get(j));
    				Push push = calculateFirstPush(a1, a2, 356 * 40);
	    			if(push != null){
	    				System.out.println("Real push");
	    				direction[push.asteroid_id] = push.direction;
	    				energy[push.asteroid_id] = push.energy;
	    				System.out.println(push.energy);
	    				time_of_push.put(i, push.time_of_push);
	    				continue;
	    			}
	    			//push it to the near outside one
	    		}
    		}
    	}
    }
    private void refreshIndexMap(Asteroid[] asteroids){
    	indexMap = new HashMap<Long, Asteroid>();
    	indexHashMap = new HashMap<Long, Integer>();
    	for(int i = 0; i < asteroids.length; i++){
    		indexMap.put(asteroids[i].id, asteroids[i]);
    		indexHashMap.put(asteroids[i].id, i);
    	}
    }
    

    //Push a1 to a2
    private Push calculateFirstPush(Asteroid a1, Asteroid a2, long t){
    	double r1 = a1.orbit.a;
    	double r2 = a2.orbit.a;
    	double vnew = Math.sqrt(Orbit.GM / r1) * (Math.sqrt(2 * r2 / (r1 + r2)) - 1);
    	double direction = Double.MAX_VALUE;
    	double energy = 0.5*a1.mass * vnew * vnew;
//    	System.out.println("Energy:" + energy);
    	long time_push = Long.MAX_VALUE;
    	
    	for(long ft = 0; ft < t; ft++) {
	    	Point v1 = a1.orbit.velocityAt(time + ft - a1.epoch);
	    	
	    	direction = Math.atan(v1.y / v1.x);
	    	double theta1 = Math.atan2(v1.y, v1.x);
	    	double timeH = Math.PI* Math.sqrt(Math.pow((r1 + r2), 3)/(8*Orbit.GM));
	    	double thresh = a1.radius() + a2.radius();
	    	Point v2 = a2.orbit.velocityAt(time + t - a2.epoch);
	    	double theta2 = Math.atan2(v2.y, v2.x);
	    	double omega2 = Math.sqrt(Orbit.GM / Math.pow(r2, 3));
	    	
	    	if(Math.abs(theta1 + Math.PI - theta2 - timeH*omega2) < thresh/2) {
	    		System.out.println("Energy:" + energy);
	    		System.out.println("Above energy tried to be pushed");
	    		direction = theta1;
	    		time_push = time + ft;
	    		//the first parameter is the index
	    		return new Push(indexHashMap.get(a1.id), energy, direction, time_push);
	    	}
    	}
    	return null;
    }
  
    private ArrayList<Double> calculatePush(double speed, double angle, double targetSpeed, double targetAngle) {
    	double vx = speed * Math.cos(angle);
    	double vy = speed * Math.sin(angle);
    	double target_vx = targetSpeed * Math.cos(targetAngle);
    	double target_vy = targetSpeed * Math.sin(targetAngle);
    	double push_vx = target_vx - vx;
    	double push_vy = target_vy - vy;
    	double push_angle = Math.atan2(push_vy, push_vx);
    	double push_speed = l2norm(push_vx, push_vy);
    	ArrayList<Double> parameters = new ArrayList<Double>();
    	parameters.add(new Double(push_speed));
    	parameters.add(new Double(push_angle));
    	return parameters;
        }
    private double l2norm(Point p) {return Math.sqrt(p.x*p.x+p.y*p.y);}
    private double l2norm(double x, double y) {return Math.sqrt(x*x+y*y);}
}

