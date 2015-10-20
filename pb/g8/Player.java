package pb.g8;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

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
	private long number_of_ast = 0;
	private Point origin = new Point(0,0);

	// time until next push
	private HashMap<Integer, Push> time_of_push;

	// number of retries
	private int retries_per_turn = 1;
	private int turns_per_retry = 3;
	private int total_number;

	private HashMap<Integer, List<Long>> asteroidClusters;
	private int cluster_number = 0;

	private HashMap<Long, Asteroid> asteroidMap;
	//key is the id of ast, value is the index of ast
	private HashMap<Long, Integer> asteroidIndexMap;
	private Set<Long> asteroidsForCircularPush;
	private double clusterThreshold;


	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit) 
	{
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
		this.total_number = asteroids.length;
		this.number_of_ast = asteroids.length;
		refreshIndexMap(asteroids);
		reorderCluster(asteroids);
		time_of_push = new HashMap<Integer, Push>();
		for(int i=0; i < cluster_number; i++) {
			time_of_push.put(i, null);
		}
		//        pushAgain = new HashMap<Integer, Boolean>();
		System.out.println("Initialization done!");
	}

	private void reorderCluster(Asteroid[] asteroids) {
		Set<Long> relevantAsteroidIds = getOutersidePart(asteroids);
		int clusterCount = numberOfClusters(relevantAsteroidIds);
		asteroidClusters = getAsteroidClusters(relevantAsteroidIds, clusterCount);
		cluster_number = asteroidClusters.size() - 1;
		mergeOnlyOneAstCluster();
	}

	private HashMap<Integer, List<Long>> getAsteroidClusters(Set<Long> relevantAsteroidIds, int clusterCount) {
		// gets the clusters for the provided dataset
		Dataset[] clusters = computeClusters(relevantAsteroidIds, clusterCount);
		HashMap<Integer, List<Long>> result = new HashMap<Integer, List<Long>>();
		Arrays.sort(clusters, new Comparator<Dataset>() {
			@Override
			public int compare(Dataset o1, Dataset o2) {
				return (int) (o1.get(0).value(0) - o2.get(0).value(0));
			}
		});
		for(int ii = 0; ii < clusters.length; ii++) {
			List<Long> asteroidIds = new ArrayList<Long>();
			for (Instance instance : clusters[ii]) {
				// gets asteroid associated with the instance
				for(Long asteroidId: relevantAsteroidIds) {
					Asteroid asteroid = asteroidMap.get(asteroidId);
					if (asteroid.orbit.a == instance.value(0)) {
						asteroidIds.add(asteroid.id);
					}
				}
			}
			result.put(ii, asteroidIds);
		}
		// add all asteroids that we don't want to consider in cluster -1
		Set<Long> remainingAsteroidIds = new HashSet<Long>(asteroidMap.keySet());
		remainingAsteroidIds.removeAll(relevantAsteroidIds);
		result.put(-1, new ArrayList<Long>(remainingAsteroidIds));

		return result;
	}

	

	private Dataset[] computeClusters(Set<Long> relevantAsteroidIds, int clusterCount) {
		// create default data set
		Dataset data = new DefaultDataset();
		// create instances and populate the data set
		for(Long asteroidId: relevantAsteroidIds) {
			Asteroid asteroid = asteroidMap.get(asteroidId);
			Instance instance = new DenseInstance(new double[] {asteroid.orbit.a});
			data.add(instance);
		}
		/* Create a new instance of the KMeans algorithm, with no options
		 * specified. By default this will generate 4 clusters.
		 */
		Clusterer km = new KMeans(clusterCount, 10000);
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
			Asteroid a1 = asteroidMap.get(id1);
			for (Long id2: relevantAsteroidIds) {
				if (id1.longValue() != id2.longValue()) {
					Asteroid a2 = asteroidMap.get(id2);
					if(Math.abs(a1.orbit.a - a2.orbit.a) < clusterThreshold) {
						nearbyAsteroids++;
					}
				}
			}
			clusters -= nearbyAsteroids / (nearbyAsteroids + 1);
		}
		int divider = (int)(number_of_ast / 100);
		divider = (int)Math.pow(2, divider);
		return (int) clusters / 2 / divider == 0 ? 1 :(int) clusters / 2 / divider;
	}

	// try to push asteroid
	public void play(Asteroid[] asteroids, double[] energy, double[] direction) 
	{
		time++;
		refreshIndexMap(asteroids);
		updateClusters(asteroids);

		int count = 0;
		Set<Integer> keys = time_of_push.keySet();

		for(Integer key: keys) {
			Push push = time_of_push.get(key);
			if(push == null) {
				continue;
			}
			count++;
			if(asteroidIndexMap.containsKey(push.asteroid_id)) {
				//                int index = asteroidIndexMap.get(push.asteroid_id);
				////            System.out.println("time of push is: " + push.time_of_push);
				//                if (push.time_of_push == time) {
				//                    System.out.println("aha pushing");
				//                    energy[index] = push.energy;
				//                    direction[index] = push.direction;
				//
				//                }
			} else {
				time_of_push.put(key, null);
			}
		}
		for(Long asteroidId: asteroidsForCircularPush) {
			circularPush(asteroidMap.get(asteroidId), energy, direction);
		}
		if(asteroidsForCircularPush.size() > 0) {
			return;
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
		tryToCollideOutside(asteroids, energy, direction);
		//        System.out.println(time_of_push);

		//============================================

		//========================================*/
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
		for(Long id: newAsteroidIds) {
			Asteroid newAsteroid = asteroidMap.get(id);
			if(newAsteroid.orbit.a != newAsteroid.orbit.b) {
				asteroidsForCircularPush.add(newAsteroid.id);
			}
		}
		// if more than 1 new asteroid, reorder and return
		// make new asteroids circular
		if(newAsteroidIds.size() != 1) {
			reorderCluster(asteroids);
		}else{
			for(Long id: newAsteroidIds) {
				newId = id;
			}
			// remove the asteroids from the cluster
			HashMap<Integer, List<Long>> newAsteroidCluster = new HashMap<Integer, List<Long>>();
			for(Map.Entry<Integer, List<Long>> entry: asteroidClusters.entrySet()) {
				List<Long> newAsteroidIdsForCluster = new ArrayList<Long>();
				for (Long id : entry.getValue()) {
					if (!idsToRemove.contains(id)) {
						newAsteroidIdsForCluster.add(id);
					}
				}
				newAsteroidCluster.put(entry.getKey(), newAsteroidIdsForCluster);
			}
			asteroidClusters = newAsteroidCluster;
			// add new id to the cluster
			asteroidClusters.get(clusterId).add(newId);
		}
		//merge clusters which only have one ast left with nearest clusters
		mergeOnlyOneAstCluster();
	}

	public void mergeOnlyOneAstCluster(){
		boolean changed = false;
		if(asteroidClusters.size() == 2)
			return;
		for(int i = 0; i < cluster_number; i++){
			List<Long> list = asteroidClusters.get(i);
			if(list.size() == 1){
				changed = true;
				int index = i == cluster_number - 1 ? cluster_number - 2 : i + 1;
				asteroidClusters.get(index).addAll(asteroidClusters.get(i));
				asteroidClusters.put(i, null);
			}/*else{
				Asteroid a1 = asteroidMap.get(list.get(0));
				Asteroid a2 = asteroidMap.get(list.get(list.size() - 1));
				double dis = Math.abs(a1.orbit.a - a2.orbit.a);
				double radius = a1.radius() + a2.radius();
				if(dis < radius){
					int index = i == cluster_number - 1 ? cluster_number - 2 : i + 1;
					asteroidClusters.get(index).addAll(asteroidClusters.get(i));
					asteroidClusters.put(i, null);
				}
			}*/
		}
		if(changed){
			int newClusterNumber = 0;
			HashMap<Integer, List<Long>> map = new HashMap<Integer, List<Long>>();
			for(int i = -1; i < cluster_number - 1; i++){
				if(asteroidClusters.get(i) != null){
					map.put(newClusterNumber, asteroidClusters.get(i));
					newClusterNumber++;
				}
			}
			cluster_number = newClusterNumber;
			asteroidClusters = map;
		}    
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
//		Asteroid middleAsteroid = asteroids[asteroids.length / 2];
//		Point v = middleAsteroid.orbit.velocityAt(time - middleAsteroid.epoch);
//		double vRandom = Math.sqrt(v.x * v.x + v.y * v.y) * 0.25;
//		boolean alreadySkip = false;
//		double abandonedMass = 0;
//		int indexOfMedian = 0;
		for(int i = asteroids.length - 1; i >= 0; i--){
			/*double a = asteroids[i].orbit.a;
			if(!alreadySkip){
				double vHomman = Math.sqrt(Orbit.GM / a) * (Math.sqrt(2 * a / (a + middleAsteroid.orbit.a) - 1));
				if(vHomman > vRandom && (abandonedMass + asteroids[i].mass < half)){
					System.out.println("abondan this ast: " + a);
					continue;
				}
			}
			*/
			if(half >= sum * 0.55){
//				indexOfMedian = i;
				break;
			}
//			alreadySkip = true;
			half += asteroids[i].mass;
			ret.add(asteroids[i].id);
			secondId = asteroids[i].id;
		}
		//ret.remove(o)
		clusterThreshold = Math.abs(asteroidMap.get(firstId).orbit.a - asteroidMap.get(secondId).orbit.a) / ret.size();
		return ret;	
	}

	private void tryToCollideOutside(Asteroid[] asteroids, double[] energy, double[] direction){	
		List<Long> ids = new ArrayList<Long>();
		Point origin = new Point(0, 0);
		//        System.out.println("clusters: " + asteroidClusters);
		//        System.out.println("cluster number: " + cluster_number);
		for(int i = 0; i < asteroidClusters.size() - 1; i++){
			debug("i: " + i);
			if(time_of_push.get(i) != null) {
				debug("i skip: " + i);
				continue;
			}
			ids = asteroidClusters.get(i);
			//List<Long> list = asteroidClusters.get(cluster_number - 1);
			//Asteroid a2 = asteroidMap.get(list.get(list.size() - 1));
			Collections.sort(ids, new Comparator<Long>() {
				@Override
				public int compare(Long l1, Long l2) {
					Asteroid a1 = asteroidMap.get(l1);
					Asteroid a2 = asteroidMap.get(l2);
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
			//    		System.out.println("cluster id: " + i + " size: " + size);
			//loop within this cluster
			if(size == 1){
				//TODO push it to the next cluster
				//                System.out.println("There is only 1");
				System.out.println("bug!!! there should not be any cluster which only have one ast");
			}else{
				//TODO
				//                System.out.println("There are multiple too!");
				Asteroid a2 = asteroidMap.get(ids.get(size - 1));
				for(int j = 0; j < ids.size() - 1; j++){
					Asteroid a1 = asteroidMap.get(ids.get(j));
//                    Push push = calculateFirstPush(a1, a2, 356 * 40, energy, direction);
                    Push push = calculateFirstPushReverse(a1, a2, 356 * 40, energy, direction);
					if(push != null){
						System.out.println("Real push");
						// do this at the time of pushnot,  immdiately
						System.out.println("time to push: " + push.time_of_push);
						System.out.println("energy: " + push.energy);
						System.out.println("collision time: " + push.time_of_collision);
						time_of_push.put(i, push);
					}
					//push it to the near outside one
				}
			}
		}


	}
	private void refreshIndexMap(Asteroid[] asteroids){
		asteroidMap = new HashMap<Long, Asteroid>();
		asteroidIndexMap = new HashMap<Long, Integer>();
		for(int i = 0; i < asteroids.length; i++){
			asteroidMap.put(asteroids[i].id, asteroids[i]);
			asteroidIndexMap.put(asteroids[i].id, i);
		}
		asteroidsForCircularPush = new HashSet<Long>();
	}


	//Push a1 to a2
	private Push calculateFirstPush(Asteroid a1, Asteroid a2, long t, double[] energy, double[] direction){
		
		double r1 = a1.orbit.a;
		double r2 = a2.orbit.a;
		Point v1 = a1.orbit.velocityAt(time - a1.epoch);
		Point v2 = a2.orbit.velocityAt(time - a2.epoch);
		double theta1 = Math.atan2(v1.y, v1.x);
		double theta2 = Math.atan2(v2.y, v2.x);
        double vnew;
        vnew = Math.sqrt(Orbit.GM / r1) * (Math.sqrt(2 * r2 / (r1 + r2)) - 1);
        double E = 0.5 * a1.mass * vnew * vnew;

		double timeH = Math.PI* Math.sqrt(Math.pow((r1 + r2), 3)/(8 * Orbit.GM));
		double thresh = a1.radius() + a2.radius();
		double omega2 = Math.sqrt(Orbit.GM / Math.pow(r2, 3));

		if(Math.abs(theta1 + Math.PI - theta2 - timeH * omega2) < thresh / r2) {
			System.out.println("Energy:" + E);
			System.out.println("Above energy tried to be pushed");
			long time_push = time;
			System.out.println("time_push: " + time_push);

			//the first parameter is the index
			System.out.println("adding new push");
			long time_of_collision = (new Double(timeH)).longValue() + time;
			//long time_of_collision_2 = calCollisionTime(a1, a2, 0, t, p1, p2);
			System.out.println(time_of_collision);
			System.out.println(timeH);
			//System.out.println(time_of_collision_2);
			//            if(time_of_collision < t) {
			int index = asteroidIndexMap.get(a1.id);
			System.out.println("index:" + index);
			energy[index] = E;
            direction[index] = theta1;
			return new Push(a1.id, E, theta1, time_push, time_of_collision);
			//            }
		}
		return null;
	}

    //Push a1 to a2
    private Push calculateFirstPushReverse(Asteroid a1, Asteroid a2, long t, double[] energy, double[] direction) {
        if(a1.orbit.a < a2.orbit.a) {
            Asteroid tmp = a2;
            a2 = a1;
            a1 = tmp;
        }
        double r1 = a2.orbit.a;
        double r2 = a1.orbit.a;
        double omegaOuter = Math.sqrt(Orbit.GM/ Math.pow(r1, 3));
        Point velocityInner = a2.orbit.velocityAt(time - a2.epoch);
        double angleInner = Math.atan2(velocityInner.y, velocityInner.x);

        Point velocityOuter = a1.orbit.velocityAt(time - a1.epoch);
        double angleOuter = Math.atan2(velocityOuter.y, velocityOuter.x);

        double timeTransfer = Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM));

        double sumRadii = a1.radius() + a2.radius();

        double alignThreshold = sumRadii/r1;

        double alignment = Math.abs(Math.PI - timeTransfer*omegaOuter - (angleInner - angleOuter));

        if (alignment < alignThreshold)
        {
            double anglePush = 0.0;
            if (angleOuter > 0)
            {
                anglePush = angleOuter - Math.PI;
            }
            else
            {
                anglePush = Math.PI + angleOuter;
            }
            //double velocityNew = Math.sqrt(Orbit.GM / r1) * (Math.sqrt( 2*r2 / (r1+r2)) - 1);
            double velocityNew = Math.sqrt(Orbit.GM / r2) * (1 - Math.sqrt( 2*r1 / (r1+r2)) );
            int index = asteroidIndexMap.get(a1.id);
            System.out.println("index:" + index);
            long time_push = time;
            long time_of_collision = (new Double(timeTransfer)).longValue() + time;
            energy[index] = 0.5*a1.mass * velocityNew * velocityNew;
            direction[index] = anglePush;
            return new Push(a1.id, energy[index], direction[index], time_push, time_of_collision);
            //store mass
            //asteroidCombinedMass = (asteroids[outerIndex].mass + asteroids[innerIndex].mass);
        }
        return null;
    }

	//make the asteroid circular again
	private void circularPush(Asteroid a, double[] energy, double[] direction){
		Point location = a.orbit.positionAt(time - a.epoch);
		Point velocity = a.orbit.velocityAt(time - a.epoch);
		double distance = Point.distance(location, origin);
		double orbit_speed = Math.sqrt(Orbit.GM / distance);
		double tangent_theta = Math.PI/2 + Math.atan2(location.y, location.x);
		double normv2 = l2norm(velocity);
		double theta2 = Math.atan2(velocity.y, velocity.x);
		ArrayList<Double> parameters = calculatePush(normv2, theta2, orbit_speed, tangent_theta);
		int index = asteroidIndexMap.get(a.id);
		energy[index] = 0.5 * a.mass * Math.pow(parameters.get(0),2);
		direction[index] = parameters.get(1);
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


	public long calCollisionTime(Asteroid a11, Asteroid a22, long startTime, long timeInterval, Point p1, Point p2)
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
			debug("overlap now");
			return time + 1;
		}

		for(int i = 0; i < cts.size(); i++)
		{
			startTime = cts.get(i) - threshold;
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
						return t;
					}
				}
				startTime = startTime + a1.orbit.period();
				endTime = endTime + a1.orbit.period();
			}
		}
		return -1;
	}
}

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
	private void debugCluster(){
		//debug the cluster
		for(int i = 0; i < cluster_number; i++){
			System.out.println("cluster id: " + i);
			for(Long id: asteroidClusters.get(i)){
				System.out.println("asteroid " + id + " a: " + asteroidMap.get(id).orbit.a);
			}
		}
	}
	private Set<Long> getAsteroidIds(Asteroid[] asteroids) {
		Set<Long> ids = new HashSet<Long>();
		for(Asteroid asteroid: asteroids) {
			ids.add(asteroid.id);
		}
		return ids;
	}
}*/