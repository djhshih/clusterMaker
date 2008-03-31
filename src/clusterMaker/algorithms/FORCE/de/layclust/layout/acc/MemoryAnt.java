package de.layclust.layout.acc;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;

public class MemoryAnt extends JumpingAntWithIncreasingViewField {

	protected int sizeOfMemory;
	protected LinkedList memory;
	protected double minSimilarityForJump = 0.5;
	protected boolean performedMemoryJump;
	
	public MemoryAnt(Playground playground) {
		super(playground);
		this.sizeOfMemory = 20;
		this.memory = new LinkedList();
	}
	
	public void makeStep(){
		makeStep(false);	
	}
	
	
	public void makeStep(boolean last) {

		//when ant carries no item, check for new item to pick up:
		//it should always carry an item...
		if(item == 0) {
			performPickUpJump();
		//otherwise compute probability to unload the item
		} else {
			//System.out.print("Ant carries item "+item+" and moves from "+Arrays.toString(position));
			performMovement();
			//System.out.println(" to "+Arrays.toString(position));
			if(playground.getLocation(position) == 0) {
				double prob = computeDropProbability();
				if(Math.random() < prob) {
					playground.setLocation(item, position);
					playground.setItemPosition(item, position);
					playground.addUncarriedItem(item);
					//System.out.println("Dropping "+item+"at "+position[0]+"/"+position[1]);
					memory.add(item);
					item = 0;
					if(memory.size() > sizeOfMemory) {
						memory.removeFirst();
					}
					if(!last) {
						performPickUpJump();
					}
				} else {
					updateAlpha(false);
				}
			} else {
				//System.out.println("Place already used.");
			}
		}
		steps++;
		//increase view-size linear:
		if(steps % 10000 == 0) {
			//System.out.println("Steps: "+steps);
			//System.out.println(kp);
		}
		if(view != maxViewSize && steps > view * maxSteps / maxViewSize) {
			view++;
			//System.out.println("View: "+ view);
		}
		//decreasing pick-up-probability:
		kp = basic_kp * (double) ( maxSteps-steps) / (maxSteps); 
	}
	
	public boolean performMemoryJump(){
		int bestItem = 0;
		double bestSim = minSimilarityForJump;
		for (Iterator iter = memory.iterator(); iter.hasNext();) {
			int element =  ((Integer)iter.next()).intValue();
			double sim = playground.getSimilarity(item-1, element-1);
			if(sim > bestSim) {
				bestSim = sim;
				bestItem = element;
			}
		}
		if(bestItem != 0) {
			//System.out.print("Jump with "+item+" from " +position[0]+"/"+position[1]);
			int[]pos = playground.getItemLocation(bestItem);
			for (int i = 0; i < pos.length; i++) {
				this.position[i] = pos[i];
			}
			//System.out.println(" to "+position[0]+"/"+position[1]+" with "+bestItem);
			return true;
		}
		return false;
	}
	
	protected void performPickUpJump() {
		//while (this.item == 0) {
			int it = playground.getUncarriedItem();
			//System.out.println("Try to pick up "+it);
			int[] pos = playground.getItemLocation(it);
			this.position = pos;
			double prob = computePickUpProbability();
			item = playground.getLocation(position);
			if(Math.random() < prob) {
				//System.out.println("Pickup "+item+"at "+position[0]+"/"+position[1]);
				playground.setLocation(0, position);
				//System.out.println("Item before memoryjump: "+item);
				if(this.memory != null && performMemoryJump()) {
//					//move until ant is on a place without an item
					while(playground.getLocation(position) == 0) {
						performMovement(1);
					}
				}
			} else {
				playground.addUncarriedItem(item);
				item = 0;
			}
		//}
	}
	
}
