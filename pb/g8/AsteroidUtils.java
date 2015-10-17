package pb.g8;

import pb.sim.Asteroid;
import pb.sim.Point;

import java.util.List;

/**
 * Created by naman on 10/17/15.
 */
public class AsteroidUtils {

    private static Point origin = new Point(0,0);

    public static int getFarthestAsteroid(Asteroid[] asteroids, long time)
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

    public static int getNearestAsteroid(Asteroid[] asteroids, long time)
    {
        int index = 0;
        double minDistance = Double.MAX_VALUE;
        Point point = new Point();
        for(int aa = 0; aa < asteroids.length; aa++) {
            Asteroid a1 = asteroids[aa];
            a1.orbit.positionAt(time - a1.epoch, point);
            double distance = Point.distance(point, origin);
            if(distance < minDistance) {
                minDistance = distance;
                index = aa;
            }
        }
        return index;
    }

    public static int getLightestAsteroid(Asteroid asteroids[])
    {
        int min = 0;
        for (int i = 0; i < asteroids.length; i++)
        {
            if (asteroids[i].mass < asteroids[min].mass)
                min = i;
        }
        return min;
    }

    public static int getHeaviestAsteroid(Asteroid asteroids[])
    {
        int max = 0;
        for (int i = 0; i < asteroids.length; i++)
        {
            if (asteroids[i].mass > asteroids[max].mass)
                max = i;
        }
        return max;
    }
}
