package uansim;

import java.util.ArrayList;
import java.util.function.Supplier;

import utils.Geometry;
import utils.SequentialIdGenerator;

public class Node {
	int id;
	
	public double[] pos;
	
	static Supplier<Integer> defaultIdGenerator = new SequentialIdGenerator();
	
	/**Construct a physical node and connect it to the physical channel
	 * @param pos initial coordinates vector (x, y, z), where z is the non-negative depth underwater
	 */
	public Node(double[] pos) {
		this(defaultIdGenerator, pos);
	}

	public Node(Supplier<Integer> idGenenerator, double[] pos) {
		this.id = idGenenerator.get();
		this.pos = pos.clone();
		this.adapters = new ArrayList<Adapter>(1);
	}
	
	public int getId() {
		return id;
	}
	
	public double[] getPosition() {
		return pos;
	}
	
	public ArrayList<Adapter> adapters;	
	
	public Adapter getDefaultAdapter() {
		if(adapters!=null) {
			if(adapters.size()>0) {
				return adapters.get(0);
			}
		}
		return null;
	}

	public void installAdapter(Adapter adapter) {
		adapters.add(adapter);		
	}
//
//	public double distanceTo(Node node) {
//		return Geometry.distance(this.getPosition(), node.getPosition()); 
//	}
}
