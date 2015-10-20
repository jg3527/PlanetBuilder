package pb.g5;

import pb.g5.GradientDescent;
import pb.g5.Push;
//import pb.g5.PushWithTime;
import pb.sim.Asteroid;
import pb.sim.Point;

public class GD_Response {
	Asteroid toPush, target;
	
	GradientDescent gd;
	int pushedIndex;
	Push push;
	long pushTime;
	long predictedTime;
	double maxEnergy;
	
	public GD_Response(Asteroid toPush, Asteroid target, long pushTime, double maxEnergy, int pushedIndex) {
		super();
		this.toPush = toPush;
		this.target = target;
		this.maxEnergy = maxEnergy;
		this.gd = new GradientDescent(toPush, target, pushTime, maxEnergy);
		this.pushedIndex = pushedIndex;
		this.pushTime = pushTime;
	}
	
	public GradientDescent getGd() {
		return gd;
	}
	public int getPushedIndex() {
		return pushedIndex;
	}
	
	public Push tuneGd() {
		push = gd.tune();
		pushTime = gd.timeOfStart;
		predictedTime = gd.predictedTime;
		return push;
	}
	
	public void tunePushTime(long minTime) {
		long delta = 300;
		int numInDirection = 0;
		while(delta > 1) {
			if(Math.abs(numInDirection) > 4) {
				delta += delta;
				numInDirection = 0;
			}
			
			GradientDescent l = new GradientDescent(toPush, target, pushTime+delta, maxEnergy);
			GradientDescent e = new GradientDescent(toPush, target, pushTime-delta, maxEnergy);
			
			Push later = l.tune();
			Push earlier = e.tune();
			
			if(later != null && later.energy < push.energy) {
				pushTime += delta;
				push = later;
				predictedTime = l.predictedTime;
				if(numInDirection < 0) numInDirection = 1;
				else ++numInDirection;
			} else if(earlier != null && earlier.energy < push.energy) {
				if(pushTime - delta < minTime)
					return;
				pushTime -= delta;
				push = earlier;
				predictedTime = e.predictedTime;
				if(numInDirection < 0) numInDirection = -1;
				else --numInDirection;
			} else {
				delta /= 2;
			}
		}
		
	}
	
}
