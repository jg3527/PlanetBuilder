package pb.g6;

// import pb.sim.Point;
// import pb.sim.Orbit;
// import pb.sim.Asteroid;
// import pb.sim.InvalidOrbitException;

import java.util.Comparator;

/**
 * Created by llxtt on 10/16/15.
 */

public class AsteroidSort {
	private int index;
	private double mass;
	private double radius;

    public AsteroidSort(int index, double mass, double radius) {
    	this.index = index;
    	this.mass = mass;
    	this.radius = radius;
    }
    public int getIndex() {
        return index;
    }
    public double getMass() {
        return mass;
    }
    public double getRadius() {
        return radius;
    }
    // Comparator for sorting the list by mass
    public static Comparator<AsteroidSort> AsteroidMassComparator = new Comparator<AsteroidSort>() {
    public int compare(AsteroidSort a1, AsteroidSort a2) {
        double a1_mass = a1.getMass();
        double a2_mass = a2.getMass();
        // descending order
        return Double.compare(a2_mass, a1_mass);
    }};
    // Comparator for sorting the list by radius
    public static Comparator<AsteroidSort> AsteroidRadiusComparator = new Comparator<AsteroidSort>() {
    public int compare(AsteroidSort a1, AsteroidSort a2) {
        double a1_radius = a1.getRadius();
        double a2_radius = a2.getRadius();
        // ascending order
        return Double.compare(a1_radius, a2_radius);
    }};
    @Override
    public String toString() {
        return "[ index=" + index + ", mass=" + mass + ", radius=" + radius + "]";
    }
}


















