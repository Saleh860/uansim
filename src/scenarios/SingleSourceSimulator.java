package scenarios;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import scenarios.Scenario4.Results;
import scenarios.Scenario4.RoutingParameters;
import scenarios.Scenario4.TopologyParameters;
import scenarios.Scenario4.TrafficParameters;
import uansim.DBR.AttackNodeLayer;
import uansim.Application;
import uansim.Channel;
import uansim.DBR;
import uansim.DPR;
import uansim.Flood;
import uansim.Network;
import uansim.Physical;
import uansim.Scheduler;
import uansim.DPR.NodeLayer;
import uansim.Application;
import uansim.Scheduler;
import uansim.Node;

public class SingleSourceSimulator {
	
	double COMMUNICATION_RANGE=150;	// (m)

	public class Topology {
		public double[] sinkPos;
		public double[][] nodePositions;
		public double[] attackerPos;
		public int sourceNode;
		public int nodeCount;
		private Random rand;

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
						
			this.nodePositions = null;
			this.attackerPos = null;
			this.sourceNode = -1;			
		}
		
		public void prepare() throws Exception {
			updateNodeCount();			
			
			generateConnectedTopology();
			
			chooseFurthestSource();
			
			//Choose the attacker node location to affect as many neighbors of the source node as possible
			placeAttacker();			
		}
		
		public void createNodes() throws Exception {
			if(nodePositions==null)
				prepare();
			
			makeChannel();
			
			//Create nodes
			makeSink(sinkPos);
			
			for(int i=1; i<nodePositions.length; i++) {
				double[] v = nodePositions[i];
				makeNode(v);
			}

			setSource(sourceNode);
			
			if(topologyParameters.attackersCount>0) {
				makeAttacker(attackerPos);
			}
			
		}
		


		public void generateConnectedTopology() {
			//Good nodes locations, including the sink
			sinkPos = new double[] {topologyParameters.deploymentSideLength/2.0, 
					topologyParameters.deploymentSideLength/2.0, 0};

			nodePositions = new double[nodeCount+1][];
			seedRandom(topologyParameters.topologySeed); 
		
			nodePositions[0] = sinkPos;
			for(int i=0; i<nodeCount; i++) {
				double[] p;
				boolean isConnected;
				do {
					p = randVector();
					isConnected = false;
					for(int j=0; j<=i; j++){
						double[] v = nodePositions[j];
						if(inCommunicationRange(p,v) && p[2] > v[2]) {
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
		
		public boolean inCommunicationRange(double[] p, double[] v) {
			return distance(p,v) <= COMMUNICATION_RANGE;
		}

		public void seedRandom(long topologySeed) {
			rand = new Random(topologySeed);
		}

		public double[] randVector() {
			double[] v = new double[3];
			v[0] = rand.nextDouble()*topologyParameters.deploymentSideLength;
			v[1] = rand.nextDouble()*topologyParameters.deploymentSideLength;
			v[2] = rand.nextDouble()*topologyParameters.deploymentDepth;	
			return v;
		}
		
		public double distance(double[] u, double[] v) {
				return Math.sqrt(Math.pow(u[0]-v[0],2)+Math.pow(u[1]-v[1],2)+Math.pow(u[2]-v[2],2));
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

				// Choose the node furthest from the sink to be source node
				// Should be out of the communication range of the sink
				if(distance>maxDistance) {
				   maxDistance = distance;
				   sourceNode = i;
				 }
			}
			attackerPos = null;
		}

		public void placeAttacker() throws Exception {
			if(nodePositions==null || sourceNode<0) {
				throw new Exception("Must generate topology and select source first");
			}

			double[] sourcePos = nodePositions[sourceNode];
			attackerPos = new double[] 
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
						System.err.println("Failed to position attacker in an optimal location");
					}    
				}
			}
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
	
	private Physical physical;
	private DPR network;
	private NodeLayer sink;
	private NodeLayer source;
	private AttackNodeLayer attacker;

	public DPR makeChannel() {
		Channel channel = new Channel();
		physical = new Physical(channel);
		return network = new DPR(physical);
	}
	
	public NodeLayer makeSink(double[] pos) {
		return sink = network.new NodeLayer(new Node(pos));
	}
	
	public NodeLayer makeNode(double[] pos) {
		return network.new NodeLayer(new Node(pos));
	}

	public NodeLayer setSource(int nodeId) {
		source = (NodeLayer) network.getNodeById(nodeId);
		return source;
	}
	
	public void makeAttacker(double[] pos) {
		attacker = network.new AttackNodeLayer(new Node(pos));
	}
	
	public void run() throws Exception {
		if(source==null)
			throw new Exception("Must set source first");
		if(sink==null)
			throw new Exception("Must create sink first");
		
		//Set model parameters in corresponding classes
		Physical.setFramingBitOverhead(trafficParameters.framingOverhead);
		COMMUNICATION_RANGE=topologyParameters.communicationRange;	
		Flood.MAX_TTL=routingParameters.maximumTTL;
		Physical.BIT_RATE=trafficParameters.channelBitRate;			
		DBR.MAXIMUM_FORWARDING_DELAY=routingParameters.maximumForwardingDelay;	
		DBR.FORWARDING_THRESHOLD=routingParameters.forwardingThreshold;
		DPR.FORWARDING_PROBABILITY[topology.sourceNode]=0.5;

		Application app = new Application(0, source, sink, 0.25, 12345)
				.setNumberOfMessages(trafficParameters.messageCount)
				.setGenerationBitRate(trafficParameters.generationRate)
				.setMessageBitLength(trafficParameters.messageLength); 
		
		if(debug) {
			Scheduler.it.enableDebug();
			network.enableDebug();
		}

		Scheduler.it.handleAll();

		if(debug) {
			physical.printStatistics();
			app.printStatistics();
		}
		
		results.forwardedCount = 0;
		for(uansim.Network.Layer node: network.getNodes()) {
			if(node!=source)
				results.forwardedCount += node.getPhysical().statistics.sentCount;
		}
	
		results.generatedCount=app.generatedCount; 
		results.sentCount=source.getPhysical().statistics.sentCount;
		results.receivedCount=app.receivedCount;
		results.deliveryRate = ((double)app.receivedCount)/(double)source.getPhysical().statistics.sentCount;		
	}
	
	public void updateForwardingProbability() {
		routingParameters.forwardingProbability = 
				routingParameters.K * (
						routingParameters.forwardingProbability 
						- results.deliveryRate + 1) ;		
	}
	
	public static void main(String[] args ) {
		SingleSourceSimulator sim = new SingleSourceSimulator();
		try {
			sim.topology.createNodes();
			sim.run();
			System.out.println("Testing " + SingleSourceSimulator.class);
			System.out.println("p\tattacks\tnodes\tdegree\tgenrtd\tsent\trecved\tforwd");
			System.out.println(String.format("%4.2f", DPR.FORWARDING_PROBABILITY[sim.topology.sourceNode]) +"\t" + 
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
