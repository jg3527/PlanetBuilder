package pb.g6;

import pb.sim.Point;
//import pb.sim.Orbit;
import pb.sim.Asteroid;
//import pb.sim.InvalidOrbitException;

import java.util.ArrayList;

/**
 * Example usage:
 * Collision collision = new Collision(asteroids);
 * for allKindsOfPush:
 * collision.updateAsteroid(0, pushedAsteroid);
 * nearestTime = collision.findCollisionTime(asteroidIndex1, asteroidIndex2);
 */
public class Collision {
    public final long TIME_LIMIT = 1000 * 24 * 60 * 60;
    public Asteroid[] asteroids;
    public int curPos;  // position of the current push
    public Asteroid oriAsteroid;
    public final double DELTA = 10;

    public Collision(Asteroid[] asteroids) {
        this.asteroids = asteroids;  // is shallow copy good enough?
        oriAsteroid = null;
    }

    public void updateAsteroid(int pushedAsteroidPos, Asteroid pushedAsteroid) {
        if (oriAsteroid != null) {
            // reverse the previous push
            asteroids[curPos] = oriAsteroid;
        }
        assert pushedAsteroidPos < asteroids.length : "Failure message";
        curPos = pushedAsteroidPos;
        oriAsteroid = asteroids[curPos];
        asteroids[curPos] = pushedAsteroid;
    }

    public Point[] getFoci(double a, double b, double A, double h, double k) {
        double c = 0.0;
        if (a == b) {
            return new Point[]{new Point(c, c)};
        } else if (a > b) {
            c = Math.sqrt(Math.pow(a, 2) - Math.pow(b, 2));
            Point f1 = new Point(c, 0.0);
            Point f2 = new Point(-c, 0.0);

            f1 = new Point(f1.x * Math.cos(A) - f1.y * Math.sin(A) + h, f1.x * Math.sin(A) + f1.y * Math.cos(A) + k);
            f2 = new Point(f2.x * Math.cos(A) - f2.y * Math.sin(A) + h, f2.x * Math.sin(A) + f2.y * Math.cos(A) + k);

            return new Point[]{f1, f2};
        } else {
            c = Math.sqrt(Math.pow(b, 2) - Math.pow(a, 2));
            Point f1 = new Point(0.0, c);
            Point f2 = new Point(0.0, -c);

            f1 = new Point(f1.x * Math.cos(A) - f1.y * Math.sin(A) + h, f1.x * Math.sin(A) + f1.y * Math.cos(A) + k);
            f2 = new Point(f2.x * Math.cos(A) - f2.y * Math.sin(A) + h, f2.x * Math.sin(A) + f2.y * Math.cos(A) + k);

            return new Point[]{f1, f2};
        }
    }

    public Point findPoint(double a, double b, double A, double h, double k, double t) {
        t = (t * Math.PI) / 180;
        return new Point(a * Math.cos(t) * Math.cos(A) - b * Math.sin(t) * Math.sin(A) + h, a * Math.cos(t) * Math.sin(A) + b * Math.sin(t) * Math.cos(A) + k);
    }

    public ArrayList<ArrayList<Point>> findIntersection(Asteroid a, int ignore_idx) {
        ArrayList<ArrayList<Point>> intersect = new ArrayList<ArrayList<Point>>(asteroids.length);

        Ellipse[] ellipses = new Ellipse[asteroids.length];

        double a1 = a.orbit.a;
        double b1 = a.orbit.b;
        double A1 = a.orbit.A;
        Point c1 = a.orbit.center();

        for (int i = 0; i < asteroids.length; i++) {
            ArrayList<Point> inner = new ArrayList<Point>();
            intersect.add(inner);

            double a2 = asteroids[i].orbit.a;
            double b2 = asteroids[i].orbit.b;
            double A2 = asteroids[i].orbit.A;
            Point c2 = asteroids[i].orbit.center();
            Point[] f = getFoci(a2, b2, A2, c2.x, c2.y);

            ellipses[i] = new Ellipse(a2, b2, A2, c2.x, c2.y, f);
        }


        Point p = null;
        for (double angle = 0; angle < 360; angle += 0.2d) {
            p = findPoint(a1, b1, A1, c1.x, c1.y, angle);
            for (int i = 0; i < asteroids.length; i++) {
                if (i == ignore_idx) {
                    continue;
                }

                double sum = 0;
                for (Point f : ellipses[i].foci) {
                    sum += Point.distance(f, p);
                }

                double radius = asteroids[i].radius();
                if (ellipses[i].foci.length == 0) {
                    throw new Error();
                } else if (ellipses[i].foci.length == 1) {
                    if (sum > ellipses[i].a - radius && sum < ellipses[i].a + radius) {
                        intersect.get(i).add(p);
                    }
                } else {
                    if (sum > 2 * ellipses[i].a - radius && sum < 2 * ellipses[i].a + radius) {
                        intersect.get(i).add(p);
                    }
                }
            }
        }

        return intersect;
    }

    /**
     * Find time of asteroid at position p relative to the creation time
     */
    public long timeAt(Asteroid a, Point p) {
        long res;
        long period = a.orbit.period();
        Point curPos = new Point();
        a.orbit.positionAt(0, curPos);
        double curAngle = Math.atan2(curPos.y, curPos.x);
        double lastAngle;
        double pAngle = Math.atan2(p.y, p.x);

        long dt = period / ((long) 2.0);
        long dt_test = dt;
        long t = 0;
        // System.out.println("==================================================");
        while (Point.distance(curPos, p) > a.radius()) {
            lastAngle = curAngle;
            // if the jump step is 0, t + 1 should be the solution
            if (dt == 0) {
                return t + 1;
            }
            // jump
            t += dt;
            a.orbit.positionAt(t, curPos);
            curAngle = Math.atan2(curPos.y, curPos.x);
            // System.out.println("try " + dt);
            // System.out.println("lastAngle " + lastAngle);
            // System.out.println("pAngle " + pAngle);
            // System.out.println("curAngle " + curAngle);
            // System.out.println("time " + t);


            if (lastAngle > 0 && curAngle < 0) {
                if ((curAngle-pAngle)*(lastAngle-pAngle) <= 0) {
                    continue;
                }
            }
            else if ((curAngle-pAngle)*(lastAngle-pAngle) >= 0) {
                continue;
            }  // not passing the point

            /* passing the point */
            // come back to before jumping
            curAngle = lastAngle;
            t -= dt; 

            // if (dt > 10) {
                dt /= 2;
            // } else {
            //     dt = 4 * dt / 5;
            // }
        }
        System.err.println("timeAt found");
        return t;
    }

    /**
     * timeAt p relative to the beginning of time
     * return -1 if cannot find any value
     */
    public long absTimeAt(Asteroid a, Point p) {
        Point p0 = new Point();

        // check for the exact point
        for (int i = -5; i < 15; i++) {
            a.orbit.positionAt(timeAt(a, p) + i, p0);
            if (Point.distance(p, p0) < a.radius()) {
                // System.out.println("timeAt: found i = " + i);
                return timeAt(a, p) + i + a.epoch;
            }
        }
        return -1;
    }

    /**
     * Find collision time when push asteroids[ignore_idx] to the state of asteroid a
     *
     * @return list of all collision time under TIME_LIMIT.
     */
    public ArrayList<Long> findCollisionTime(Asteroid a, int ignore_idx) {
        ArrayList<ArrayList<Point>> intersect = findIntersection(a, ignore_idx);
        ArrayList<Long> collision = new ArrayList<Long>();
        long t1, t2, period1, period2, colTime;
        double k1, k2;


        for (int i = 0; i < intersect.size(); i++) {
            if (i == ignore_idx) {
                continue;
            }

            // find all collisions between a and asteroids[i]
            period1 = a.orbit.period();
            period2 = asteroids[i].orbit.period();
            for (int j = 0; j < intersect.get(i).size(); j++) {
                Point p = intersect.get(i).get(j);

                t1 = absTimeAt(a, p);
                t2 = absTimeAt(asteroids[i], p);
                if (t1 == -1 || t2 == -1) {
                    // if can't find t1 or t2
                    continue;
                }

                // find int values of k1 and k2 such that t1 + k1*period1 = t2 + k2*period2
                k2 = 1;
                do {

                    k1 = (double)(k2*period2 + t2 - t1) / (double)period1;
                    ++k2;
                } while (k2 < (TIME_LIMIT / period2) && Math.round(k1) != k1);

                colTime = t1 + (long)k1*period1;  // collision time

                // check
                Point p1 = new Point();
                Point p2 = new Point();
                a.orbit.positionAt(colTime - a.epoch, p1);
                asteroids[i].orbit.positionAt(colTime - asteroids[i].epoch, p2);
                double r = a.radius() + asteroids[i].radius();
                if (Point.distance(p1, p2) >= r) {
                    System.out.println("No collision");
                    System.out.println("k1 " + k1 + " k2 " + (k2-1));
                    System.out.println("t1 " + colTime + " t2 " + (t2 + (k2-1)*period2));
                    continue;
                }

                System.out.println("Collision between " + ignore_idx + " & " + i);
                System.out.println("  Year: " + (1 + colTime / 365));
                System.out.println("  Day: " + (1 + colTime % 365));
                collision.add(colTime);
            }
        }

        return collision;
    }
}

