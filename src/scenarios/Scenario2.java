package scenarios;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import uansim.Channel;
import uansim.DBR;
import uansim.DPR;
import uansim.Flood;
import uansim.Physical;
import uansim.DPR.NodeLayer;
import uansim.Application;
import uansim.Scheduler;
import utils.Geometry;
import uansim.Application;
import uansim.Scheduler;
import uansim.Node;

public class Scenario2 {

	double m_NetworkDensity=1;
	int m_numAttackers=0;
	double m_boundary=1500;
	double m_depth=100;
	long topologySeed=12341;
	int numberOfMessages=2;
	double generationRate=1000;
	int messageLength=100;
	boolean debug=true;
	double COMMUNICATION_RANGE=150;

	public void run() {
		Channel channel=new Channel();
		Physical physical = new Physical(channel);
		DPR dbr = new DPR(physical);
		Physical.setFramingBitOverhead(0);
		Flood.MAX_TTL=5;
		Physical.BIT_RATE=30000;
		DBR.MAXIMUM_FORWARDING_DELAY=100;
		DBR.FORWARDING_THRESHOLD=0;
		Arrays.fill(DPR.FORWARDING_PROBABILITY,0.0);
		

		int m_numNodes=(int) Math.ceil((m_NetworkDensity) * m_boundary * m_boundary  * m_depth /
                (4/3 * Math.PI * Math.pow(COMMUNICATION_RANGE,3)));
		
		log("Number of underwater nodes = " + m_numNodes);
		

//#ifdef UAN_PROP_BH_INSTALLED
//Ptr<UanPropModelBh> prop = CreateObjectWithAttributes<UanPropModelBh> ("ConfigFile", StringValue ("exbhconfig.cfg"));
//#else 
//Ptr<UanPropModelIdeal> prop = CreateObjectWithAttributes<UanPropModelIdeal> ("CommunicationRange", DoubleValue(m_commRange));
//#endif //UAN_PROP_BH_INSTALLED
//Ptr<UanChannel> channel = CreateObjectWithAttributes<UanChannel> ("PropagationModel", PointerValue (prop));


//Good nodes locations
		Random rand = new Random(topologySeed); 
		class Vector{
			double[] pos;
			public Vector(double x, double y, double z) {
				pos=new double[] {x, y, z};
			}
			public Vector() {
				 pos=new double[]{rand.nextDouble()*m_boundary,
						 rand.nextDouble()*m_boundary,
						 rand.nextDouble()*m_depth};
			}
			
			double distanceTo(Vector v) {
				return Geometry.distance(pos, v.pos); 
			}
			
			double x() {return pos[0];}
			double y() {return pos[1];}
			double z() {return pos[2];}
		}
		
		Vector sinkPos = new Vector(m_boundary/2.0, m_boundary/2.0, 0);
		ArrayList<Vector> connectedTopology = new ArrayList<Vector>();
		connectedTopology.add(new Vector(m_boundary/2.0, m_boundary/2.0, 0.0));
		for(int i=0; i<m_numNodes; i++) {
			Vector p;
			boolean isConnected;
			do {
				p = new Vector();
				isConnected = false;
				for(int j=0; j<connectedTopology.size(); j++){
					Vector v = connectedTopology.get(j);
					double distance = p.distanceTo(v);
					if(distance <= COMMUNICATION_RANGE && p.z() > v.z()) {
						log("Node #" + (i+1) + " is connected to Node #" + (j));
						isConnected = true;
						break;
					}
				}
			} while(!isConnected);
			
			connectedTopology.add(p);
		}
		
		
		int sourceNode=-1;
		double maxDistance=0;

		//Choose the furthest node from the sink to be the source
		for(int i=0; i<m_numNodes; i++) {
			Vector v = connectedTopology.get(i+1);
			//pos->Add ( v );
		
			double distance = v.distanceTo(sinkPos);
			log("Node #" + (i+1) + " at\t" + v.x() + ", \t" + v.y() + ", \t" + v.z() + "  \td(sink) = " + distance);
			// Choose the node furthest from the sink to be source node
			// Should be out of the communication range of the sink
			if(distance>maxDistance) {
			   log("New max d(sink)");
			   maxDistance = distance;
			   sourceNode = i;
			 }
		}
		log("Selected node #" + (sourceNode+1) + " as source.");
		
		//Choose the attacker node location to affect as many neighbors of the source node as possible
		Vector sourcePos = connectedTopology.get(sourceNode);
		Vector attackerPos = new Vector(sourcePos.x(), sourcePos.y(), sourcePos.z() +1.0);

		//Verify the network density by calculating the number of forwarders within the communication
		//range of the source node
		//Verify that the attacker can reach all neighbors of source but reaches no other nodes
		ArrayList<Integer> sourceNeighbors=new ArrayList<Integer>();
		ArrayList<Integer> vulnerableNeighbors=new ArrayList<Integer>();
		for(int i=0; i<m_numNodes; i++) {
			boolean neighbor=false, vulnerable=false; 
			Vector nodePos = connectedTopology.get(i);
			if(i!=sourceNode){
				if(nodePos.distanceTo(sourcePos) <= COMMUNICATION_RANGE) {
					sourceNeighbors.add(i);
					neighbor = true;
				}
				if(nodePos.distanceTo(attackerPos) <= COMMUNICATION_RANGE) {
					 vulnerableNeighbors.add(i);
					 vulnerable = true;
				}
				if(neighbor!=vulnerable) {
					log("Failed to position attacker in an optimal location");
				}    
			}
		}

		log("Number of neighbors of source nodes = " + sourceNeighbors.size());


		//Create nodes
		NodeLayer sink = dbr.new NodeLayer(new Node(sinkPos.pos));
		
		NodeLayer[] nodes = new NodeLayer[m_numNodes];
		
		for(int i=0; i<m_numNodes; i++) {
			Vector v = connectedTopology.get(i);
			nodes[i] = dbr.new NodeLayer(new Node(v.pos));
		}
		NodeLayer source = nodes[sourceNode];
		
		if(m_numAttackers>0) {
			DPR.AttackNodeLayer attacker = dbr.new AttackNodeLayer(new Node(attackerPos.pos));
		}
		
		Application app = new Application(0, source, sink, 0.25, 12345)
				.setNumberOfMessages(numberOfMessages)
				.setGenerationBitRate(generationRate)
				.setMessageBitLength(messageLength); 
		

		if(debug) {
			Scheduler.it.enableDebug();
			dbr.enableDebug();
		}

		Scheduler.it.handleAll();

		if(debug) {
			physical.printStatistics();
			app.printStatistics();
		}
		
		int forwardedCount = 0;
		for(int i=0; i<nodes.length; i++) {
			if(i!=sourceNode)
				forwardedCount += nodes[i].getPhysical().statistics.sentCount;
		}


		log("total packets sent by source = " + source.getPhysical().statistics.sentCount);
		log("total packets received = " + app.receivedCount);
		log("total packets forwarded = " + forwardedCount);

		System.out.println("p\tattacks\tnodes\tdegree\tsent\trecved\tforwd");
		System.out.println(DPR.FORWARDING_PROBABILITY[0] +"\t" + m_numAttackers + "\t" +
				m_numNodes + "\t" + sourceNeighbors.size() + "\t" + source.getPhysical().statistics.sentCount+"\t" +
				app.receivedCount + "\t" + forwardedCount);
		
	}

	private void log(String string) {
		if(debug)
			System.out.println("INFO:"+string);
	}

	public static void main(String[] args ) {
		new Scenario2().run();
	}
}
