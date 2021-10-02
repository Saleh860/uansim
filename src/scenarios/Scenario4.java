package scenarios;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import uansim.DPR;
import uansim.DPR.NodeLayer;
import uansim.Flood;
import uansim.Physical;
import uansim.Application;
import uansim.Channel;
import uansim.DBR;
import uansim.Scheduler;
import uansim.Node;

public class Scenario4 {

	double COMMUNICATION_RANGE=150;	// (m)
	
	public class Topology {
		public double[] sinkPos;
		public double[][] nodePositions;
		public double[] attackerPos;
		public int sourceNode;
		public int nodeCount;

		public Topology() {
			this.nodeCount = 0;
			this.sinkPos = null;
			this.nodePositions = null;
			this.attackerPos = null;
			this.sourceNode = -1;
		}

		public void updateNodeCount() {
			double networkDensity = topologyParameters.networkDensity;
			nodeCount=(int) Math.ceil((networkDensity) * 
					Math.pow(topologyParameters.deploymentSideLength,3)  /
		            (4/3 * Math.PI * Math.pow(topologyParameters.communicationRange,3)));
			
			log("Number of underwater nodes = " + nodeCount);
			
			sinkPos = new double[] {topologyParameters.deploymentSideLength/2.0, 
					topologyParameters.deploymentSideLength/2.0, 0};
			
			this.nodePositions = null;
			this.attackerPos = null;
			this.sourceNode = -1;			
		}
		
		public void create() throws Exception {
			updateNodeCount();
			
			generateConnectedTopology();
			
			chooseFurthestSource();		
			
			//Choose the attacker node location to affect as many neighbors of the source node as possible
			placeAttacker();
			
			//Count source node neighbors
			countSourceNeighbors();
		}

		public void generateConnectedTopology() {
			//Good nodes locations, including the sink
			nodePositions = new double[nodeCount+1][];
			rand = new Random(topologyParameters.topologySeed); 
		
			nodePositions[0] = sinkPos;
			for(int i=0; i<nodeCount; i++) {
				double[] p;
				boolean isConnected;
				do {
					p = randVector();
					isConnected = false;
					for(int j=0; j<=i; j++){
						double[] v = nodePositions[j];
						double distance = distance(p,v);
						if(distance <= COMMUNICATION_RANGE && p[2] > v[2]) {
							log("Node #" + (i+1) + " is connected to Node #" + (j));
							isConnected = true;
							break;
						}
					}
				} while(!isConnected);
				
				nodePositions[i+1]=p;
			}
			sourceNode = -1;
			attackerPos = null;
		}

		public int countSourceNeighbors() throws Exception {
			if(nodePositions==null || sourceNode<0) {
				throw new Exception("Must generate topology and select source first");
			}
			
			double[] sourcePos = nodePositions[sourceNode];

			//Verify the network density by calculating the number of forwarders within the communication
			//range of the source node
			ArrayList<Integer> sourceNeighbors=new ArrayList<Integer>();
			for(int i=0; i<nodePositions.length; i++) {
				double[] nodePos = nodePositions[i];
				if(i!=sourceNode){
					if(distance(nodePos,sourcePos) <= COMMUNICATION_RANGE) {
						sourceNeighbors.add(i);
					}
				}
			}
			log("Number of neighbors of source nodes = " + sourceNeighbors.size());		
			return sourceNeighbors.size();
		}
		
		public void chooseFurthestSource() throws Exception {
			if(nodePositions==null) {
				throw new Exception("Must generate topology first");
			}
			sourceNode=-1;
			double maxDistance=0;
			double[] sinkPos = nodePositions[0];
			//Choose the furthest node from the sink to be the source
			for(int i=1; i<nodePositions.length; i++) {
				double[] v = nodePositions[i];
				//pos->Add ( v );
			
				double distance = distance(v, sinkPos);
				log("Node #" + (i) + " at\t" + v[0] + ", \t" + v[1] + ", \t" + v[2] + "  \td(sink) = " + distance);
				// Choose the node furthest from the sink to be source node
				// Should be out of the communication range of the sink
				if(distance>maxDistance) {
				   log("New max d(sink)");
				   maxDistance = distance;
				   sourceNode = i;
				 }
			}
			log("Selected node #" + sourceNode + " as source.");
			attackerPos = null;
		}

		public void placeAttacker() throws Exception {
			if(nodePositions==null || sourceNode<0) {
				throw new Exception("Must generate topology and select source first");
			}

			double[] sourcePos = nodePositions[sourceNode];
			double[] attackerPos = new double[] 
					{sourcePos[0], sourcePos[1], sourcePos[2] +1.0};

			//Verify that the attacker can reach all neighbors of source but reaches no other nodes
			ArrayList<Integer> vulnerableNeighbors=new ArrayList<Integer>();
			for(int i=0; i<nodePositions.length; i++) {
				boolean neighbor=false, vulnerable=false; 
				double[] nodePos = nodePositions[i];
				if(i!=sourceNode){
					if(distance(nodePos, sourcePos) <= COMMUNICATION_RANGE) {
						neighbor = true;
					}
					if(distance(nodePos, attackerPos) <= COMMUNICATION_RANGE) {
						 vulnerableNeighbors.add(i);
						 vulnerable = true;
					}
					if(neighbor!=vulnerable) {
						log("Failed to position attacker in an optimal location");
					}    
				}
			}
			log("Placed attacker at " + attackerPos[0] + "\t" + attackerPos[1] + "\t" +attackerPos[2]);
		}
	}

	public static class TopologyParameters {
		public double networkDensity;
		public int attackersCount;
		public double deploymentSideLength;
		public double deploymentDepth;
		public long topologySeed;
		public double communicationRange;

		public TopologyParameters() {
			this.networkDensity = 1;
			this.attackersCount = 0;
			this.deploymentSideLength = 500;
			this.deploymentDepth = 100;
			this.topologySeed = 12341;
			this.communicationRange = 150;
		}
	}

	public static class RoutingParameters {
		public int maximumTTL;
		public double maximumForwardingDelay;
		public double forwardingThreshold;
		public double forwardingProbability;
		public double K;

		public RoutingParameters(int maximumTTL, double maximumForwardingDelay, double forwardingThreshold,
				double forwardingProbability, double k) {
			this.maximumTTL = maximumTTL;
			this.maximumForwardingDelay = maximumForwardingDelay;
			this.forwardingThreshold = forwardingThreshold;
			this.forwardingProbability = forwardingProbability;
			K = k;
		}
	}

	public static class TrafficParameters {
		public int messageCount;
		public double generationRate;
		public int messageLength;
		public int framingOverhead;
		public double channelBitRate;

		public TrafficParameters(int numberOfMessages, double generationRate, int messageLength, int framingOverhead,
				double channelBitRate) {
			this.messageCount = numberOfMessages;
			this.generationRate = generationRate;
			this.messageLength = messageLength;
			this.framingOverhead = framingOverhead;
			this.channelBitRate = channelBitRate;
		}
	}
	
	public static class Results {
		public double deliveryRate;
		public int generatedCount;
		public int sentCount;
		public int receivedCount;
		public int forwardedCount;

		public Results() {
			this.deliveryRate = 0;
			this.generatedCount = 0;
			this.sentCount = 0;
			this.receivedCount = 0;
			this.forwardedCount = 0;
		}
	}

	public TopologyParameters topologyParameters = new TopologyParameters();

	public RoutingParameters routingParameters = new RoutingParameters(5, 10, 0.0, 0.0, 2.0/1.0);

	public TrafficParameters trafficParameters = new TrafficParameters(100, 100, 100, 0, 30000/8);
	
	public Results results = new Results();

	public boolean debug=false;
	
	public Topology topology = new Topology();

	Random rand =  new Random();
	
	double[] randVector() {
		double[] v = new double[3];
		v[0] = rand.nextDouble()*topologyParameters.deploymentSideLength;
		v[1] = rand.nextDouble()*topologyParameters.deploymentSideLength;
		v[2] = rand.nextDouble()*topologyParameters.deploymentDepth;	
		return v;
	}
	
	double distance(double[] u, double[] v) {
			return Math.sqrt(Math.pow(u[0]-v[0],2)+Math.pow(u[1]-v[1],2)+Math.pow(u[2]-v[2],2));
	}

	public void run() throws Exception {
		//Set model parameters in corresponding classes
		Physical.setFramingBitOverhead(trafficParameters.framingOverhead);
		COMMUNICATION_RANGE=topologyParameters.communicationRange;	
		DBR.MAX_TTL=routingParameters.maximumTTL;
		Physical.BIT_RATE=trafficParameters.channelBitRate;			
		DBR.MAXIMUM_FORWARDING_DELAY=routingParameters.maximumForwardingDelay;	
		DBR.FORWARDING_THRESHOLD=routingParameters.forwardingThreshold;
		Arrays.fill(DPR.FORWARDING_PROBABILITY,0.0);

		if(topology.nodePositions==null)
			throw new Exception("Must create network topology first");

		//Run experiment
		Channel channel=new Channel();
		Physical physical = new Physical(channel);
		DPR dbr = new DPR(physical);
		

		//#ifdef UAN_PROP_BH_INSTALLED
		//Ptr<UanPropModelBh> prop = CreateObjectWithAttributes<UanPropModelBh> ("ConfigFile", StringValue ("exbhconfig.cfg"));
		//#else 
		//Ptr<UanPropModelIdeal> prop = CreateObjectWithAttributes<UanPropModelIdeal> ("CommunicationRange", DoubleValue(m_commRange));
		//#endif //UAN_PROP_BH_INSTALLED
		//Ptr<UanChannel> channel = CreateObjectWithAttributes<UanChannel> ("PropagationModel", PointerValue (prop));

		
		//Create nodes
		NodeLayer sink = dbr.new NodeLayer(new Node(topology.sinkPos));
		
		NodeLayer[] nodes = new NodeLayer[topology.nodePositions.length-1];
		
		for(int i=1; i<topology.nodePositions.length; i++) {
			double[] v = topology.nodePositions[i];
			nodes[i-1] = dbr.new NodeLayer(new Node(v));
		}
		NodeLayer source = nodes[topology.sourceNode-1];
		
		if(topologyParameters.attackersCount>0) {
			DPR.AttackNodeLayer attacker = dbr.new AttackNodeLayer(new Node(topology.attackerPos));
		}
		
		Application app = new Application( 0, source, sink, 0.25, 12345)
				.setNumberOfMessages(trafficParameters.messageCount)
				.setGenerationBitRate(trafficParameters.generationRate)
				.setMessageBitLength(trafficParameters.messageLength); 
		

		if(debug) {
			Scheduler.it.enableDebug();
			dbr.enableDebug();
		}

		Scheduler.it.handleAll();

		if(debug) {
			physical.printStatistics();
			app.printStatistics();
		}
		
		results.forwardedCount = 0;
		for(int i=0; i<nodes.length; i++) {
			if(i!=topology.sourceNode)
				results.forwardedCount += nodes[i].getPhysical().statistics.sentCount;
		}
		results.generatedCount=app.generatedCount; 
		results.sentCount=source.getPhysical().statistics.sentCount;
		results.receivedCount=app.receivedCount;


		log("total packets sent by source = " + results.sentCount);
		log("total packets received = " + results.receivedCount);
		log("total packets forwarded = " + results.forwardedCount);

		results.deliveryRate = ((double)app.receivedCount)/(double)source.getPhysical().statistics.sentCount;
	}

	public void updateForwardingProbability() {
		for(int i=0; i<DPR.FORWARDING_PROBABILITY.length; i++)
			DPR.FORWARDING_PROBABILITY[i] = (
					DPR.FORWARDING_PROBABILITY[i] - results.deliveryRate + 1
					) / routingParameters.K;		
	}


	private void log(String string) {
		if(debug)
			System.out.println("INFO:"+string);
	}

	public static void main(String[] args ) {
		Scenario4 sim = new Scenario4();
		try {
			sim.topology.create();
			sim.run();
			System.out.println("p\tattacks\tnodes\tdegree\tgenrtd\tsent\trecved\tforwd");
			System.out.println(String.format("%4.2f", DPR.FORWARDING_PROBABILITY[0]) +"\t" + 
					sim.topologyParameters.attackersCount + "\t" +
					sim.topology.nodeCount + "\t" + 
					sim.topology.countSourceNeighbors()+ "\t" + 
					sim.results.generatedCount + "\t" + sim.results.sentCount+"\t" + 
					sim.results.receivedCount + "\t" + sim.results.forwardedCount);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
