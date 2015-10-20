package pb.g5;


import pb.sim.Asteroid;
import pb.sim.Orbit;
import pb.sim.Point;
import pb.g5.Push;

public class Hohmann {

    public static double transfer(Asteroid asteroid, double r2){

        double r1 = asteroid.orbit.a;
        double deltaV = Math.sqrt(Orbit.GM/r1) * (Math.sqrt((2 * r2) / (r1 + r2)) - 1);

        double hohmannEnergy = 0.5 * asteroid.mass * deltaV * deltaV;

        return hohmannEnergy;
    }

    public static double reverseTransfer(Asteroid asteroid, double r2){

        double r1 = asteroid.orbit.a;
        double deltaV = Math.sqrt(Orbit.GM/r2) * (1 - Math.sqrt((2 * r2) / (r1 + r2)));

        double hohmannEnergy = 0.5 * asteroid.mass * deltaV * deltaV;

        return hohmannEnergy;
    }

    public static double direction(Asteroid asteroid, long time){
        Point p1 = asteroid.orbit.velocityAt(time - asteroid.epoch);
        return Math.atan2(p1.y, p1.x);
    }

    public static int getLowestAverageHohmanTransfer(Asteroid[] asteroids){

        int bestIndex = -1;
        double bestEnergy = Double.MAX_VALUE;

        for(int i = 0; i < asteroids.length; ++i){
            double energySum = 0;
            for(int j = 0; j < asteroids.length; ++j){
                if(i == j) continue;
                energySum += transfer(asteroids[j], asteroids[i].orbit.a);
            }

            energySum /= asteroids[i].mass;

            if(energySum < bestEnergy){
                bestEnergy = energySum;
                bestIndex = i;
            }

        }

        return bestIndex;

    }

    public static Push findTransferByTime(Asteroid asteroid, Asteroid collideWith, long currentTime, long collideWithin){

        Push bestPush = null;

        try {

            // total overlap radius
            double rTotal = asteroid.radius() + collideWith.radius();

            // direction of pushed asteroid
            double direction = direction(asteroid, currentTime);

            // creating points up front to save memory
            Point asteroidPosition = new Point(), collideWithPosition = new Point();

            for (long time = currentTime; time < currentTime + collideWithin; ++time) {

                // get the radius of the asteroid being collided with at the future time
                collideWith.orbit.positionAt(time - collideWith.epoch, collideWithPosition);
                double collideRadius =  collideWith.orbit.a;

                // calculate the energy required to move the second asteroid to that radius
                double energy = transfer(asteroid, collideRadius);

                // simulate the push
                Asteroid pushed = Asteroid.push(asteroid, currentTime, energy, direction);

                // get the point at which the new asteroid will be at the future time
                pushed.orbit.positionAt(time - pushed.epoch, asteroidPosition);

                long days = Utils.asteroidsCollide(pushed, collideWith, currentTime, 3650);
                if(days > currentTime){

                    if(bestPush == null){
                        bestPush = new Push(energy, direction);
                    } else if( energy < bestPush.energy ) {
                        bestPush = new Push(energy, direction);
                    }
                }
            }
        } catch (Exception e ) {
            //e.printStackTrace();
        }

        return bestPush;
    }

}
