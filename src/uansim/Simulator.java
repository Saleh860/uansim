package uansim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.RejectedExecutionException;
import java.util.function.Supplier;

import uansim.DBR.AttackNodeLayer;
import uansim.DPR.NodeLayer;
import utils.SequentialIdGenerator;

public class Simulator {
	
	
	public class Topology {
		public double[] sinkPos;
		public double[][] nodePositions;
		public double[][] attackerPositions;
		public int[] sourceNodes;
		public int nodeCount;
		private Random rand=new Random();

		public Topology() {
			this.nodeCount = 0;
			this.sinkPos = null;
			this.nodePositions = null;
			this.attackerPositions = new double[0][];
			this.sourceNodes = null;
		}

		public void updateNodeCount() {
			double networkDensity = topologyParameters.networkDensity;
			nodeCount=(int) Math.ceil((networkDensity) * 
					Math.pow(topologyParameters.deploymentSideLength,3)  /
		            (4/3 * Math.PI * Math.pow(topologyParameters.communicationRange,3)));
						
			this.nodePositions = null;
			this.attackerPositions = null;
			this.sourceNodes = null;			
		}
		
		/**Update the topology to reflect topology parameters and 
		 * adjust attacker positions at the default positions
		 * @throws Exception
		 */
		public void prepare() throws Exception {
			updateNodeCount();			
			
			generateConnectedTopology();
			
			for(int i=0; i<topologyParameters.sourcesCount; i++)
				chooseNextFurthestSource();
			
			//Choose the attacker node location to affect as many neighbors of the source node as possible
			for(int i=0; i<sourceNodes.length && i<topologyParameters.attackersCount; i++)
				placeAttacker(sourceNodes[i]);
		}
		
		public void createModel() throws Exception {
			if(nodePositions==null)
				prepare();

			Physical.setFramingBitOverhead(trafficParameters.framingBitOverhead);
//			Physical.COMMUNICATION_RANGE=topologyParameters.communicationRange;	
//			Physical.INTERFERENCE_RANGE=topologyParameters.interferenceRange;
			
			makeNetwork();
						
			//Create nodes
			makeSink(sinkPos);
			nodes = new NodeLayer[nodePositions.length-1];
			for(int i=1; i<nodePositions.length; i++) {
				double[] v = nodePositions[i];
				nodes[i-1] = makeNode(v);
			}

			setSources(sourceNodes);

			makeAttackers(attackerPositions);
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
			sourceNodes = null;
			attackerPositions = null;
		}
		
		public boolean inCommunicationRange(double[] p, double[] v) {
			return distance(p,v) <= topologyParameters.communicationRange;
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


		public int countSourceNeighbors(int sourceId) throws Exception {
			
			if(nodePositions==null || sourceId<0) {
				throw new Exception("Must generate topology and select source first");
			}
			
			double[] sourcePos = nodePositions[sourceId];

			//Verify the network density by calculating the number of forwarders within the communication
			//range of the source node
			ArrayList<Integer> sourceNeighbors=new ArrayList<Integer>();
			for(int i=0; i<nodePositions.length; i++) {
				double[] nodePos = nodePositions[i];
				if(i!=sourceId){
					if(inCommunicationRange(nodePos,sourcePos)) {
						sourceNeighbors.add(i);
					}
				}
			}
			return sourceNeighbors.size();
		}
		
		public void chooseNextFurthestSource() throws Exception {
			if(nodePositions==null) {
				throw new Exception("Must generate topology first");
			}
			if(sourceNodes==null)
				sourceNodes=new int[0];
			
			int sourceNode=-1;
			double maxDistance=0;
			double[] sinkPos = nodePositions[0];
			//Choose the furthest node from the sink to be the source
			for(int i=1; i<nodePositions.length; i++) {
				boolean excluded=false;
				for(int j:sourceNodes) {
					if(i==j) {
						excluded=true;
						break;
					}
				}
				if(!excluded) {
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
			}
			sourceNodes = Arrays.copyOf(sourceNodes, sourceNodes.length+1);
			sourceNodes[sourceNodes.length-1] = sourceNode;
			attackerPositions = new double[0][];
		}

		
		public void clearAttackers() {
			attackerPositions = new double[0][];
		}

		public void placeAttacker(double[] attackerPosition) throws Exception {
			attackerPositions = Arrays.copyOf(attackerPositions, attackerPositions.length+1);
			attackerPositions [attackerPositions.length-1] = attackerPosition; 
		}

		public void placeAttacker(int sourceNode) throws Exception {
			if(nodePositions==null || sourceNode<0) {
				throw new Exception("Must generate topology and select source first");
			}

			double[] sourcePos = nodePositions[sourceNode];
			double[] attackerPosition = new double[] 
					{sourcePos[0], sourcePos[1], sourcePos[2] +1.0};

			placeAttacker(attackerPosition);

			//Verify that the attacker can reach all neighbors of source but reaches no other nodes
			ArrayList<Integer> vulnerableNeighbors=new ArrayList<Integer>();
			for(int i=0; i<nodePositions.length; i++) {
				boolean neighbor=false, vulnerable=false; 
				double[] nodePos = nodePositions[i];
				if(i!=sourceNode){
					if(inCommunicationRange(nodePos, sourcePos)) {
						neighbor = true;
					}
					if(inCommunicationRange(nodePos, attackerPosition)) {
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

	/**The simulated topology parameters, specified before running the simulation and used for generating the simulated network.
	 * Accessed through the {@link Simulator#topologyParameters} 
	 * @author Saleh
	 *
	 */
	public static class TopologyParameters {
		/**
		 * The average number of neighbors
		 */
		public double networkDensity;
		/**
		 * The number of nodes generating traffic
		 */
		public int sourcesCount;
		/**
		 * The number of malicious nodes performing the depth-spoofing attack
		 */
		public int attackersCount;
		/**
		 * Deployment area is square-shaped and this is the side length in meters
		 */
		public double deploymentSideLength;
		/**
		 * The deployment volume depth, i.e. the maximum depth of any node
		 */
		public double deploymentDepth;
		/**
		 * The seed of the random number generator controlling the generation of nodes and other topology-related configuration
		 */
		public long topologySeed;
		/**
		 * The nominal communication range (at which the received signal to noise ratio is expected to be well above the detection threshold)
		 */
		public double communicationRange;
		/**
		 * The nominal interference range (at which the non-intended received signal would causes a noise powerful enough to 
		 * make the signal to noise ration below the detection threshold)
		 */
		public double interferenceRange;

		protected TopologyParameters() {
			this.networkDensity = 1;
			this.sourcesCount = 1;
			this.attackersCount = 0;
			this.deploymentSideLength = 500;
			this.deploymentDepth = 100;
			this.topologySeed = 8792766196518812341L;
			this.communicationRange = 150;
			this.interferenceRange = 2*communicationRange;
		}
	}

	/**
	 * The simulated topology parameters, specified before running the simulation and used for generating the simulated network.
	 */
	public final TopologyParameters topologyParameters = new TopologyParameters();

	
	/**The simulated routing protocol parameters, specified before running the simulation and used for controlling routing decisions during simulation.
	 * @author Saleh
	 *
	 */
	public static class RoutingParameters {
		/**
		 * The maximum number of hops before dropping the packet. Useful for avoiding routing loops.
		 */
		public int maximumTTL=5;
		/**
		 * The DBR maximum forwarding delay
		 */
		public double maximumForwardingDelay=10;
		public double forwardingThreshold=0;
		public double[] forwardingProbability=new double[DPR.FORWARDING_PROBABILITY.length];
		public double K=0;
		public int decisionSeed=12342;
	}

	public static class TrafficParameters {
		public int messageCount=100;
		public double generationBitRate=100;
		public int messageBitLength=100;
		public int framingBitOverhead=0;
		public double channelBitRate=30000;
		public int trafficSeed=101010;
	}
	
	public static class Results {
		public double[] deliveryRatio; //per source
		public double[] e2eDelay;	//per source
		public double finishTime;
		public int generatedCount;
		public int sentCount;
		public int receivedCount;
		public int forwardedCount;
		public double averageE2EDelay;
		public int[] numberOfTransmissions;

		public Results() {
			this.generatedCount = 0;
			this.sentCount = 0;
			this.receivedCount = 0;
			this.forwardedCount = 0;
		}
	}

	public RoutingParameters routingParameters = new RoutingParameters();

	public TrafficParameters trafficParameters = new TrafficParameters();
	
	public Results results = new Results();

	public boolean debug=false;
	
	public Topology topology = new Topology();
	
	public Physical physical;
	public DPR network;
	public NodeLayer sink;
	public NodeLayer[] sources;
	public NodeLayer[] nodes;
	public AttackNodeLayer[] attackers;
	Supplier<Integer> idGenerator;


	public DPR makeNetwork() {
		idGenerator = new SequentialIdGenerator();
		Channel channel=new Channel();
		physical = new Physical(channel);
		return network = new DPR(physical, new Random(routingParameters.decisionSeed));
	}
	
	public NodeLayer makeSink(double[] pos) {
		return sink = network.new NodeLayer(new Node(idGenerator, pos));
	}
	
	public NodeLayer makeNode(double[] pos) {
		return network.new NodeLayer(new Node(idGenerator, pos)); 
	}

	public NodeLayer[] setSources(int[] sourceNodes) {
		
		sources = new NodeLayer[sourceNodes.length];
		for(int i=0; i<sourceNodes.length;i++)
			sources[i] = (NodeLayer) network. getNodeById(sourceNodes[i]);
		return sources;
	}
	
	public void makeAttackers(double[][] attackerPositions) {
		attackers = new AttackNodeLayer[attackerPositions.length];
		for(int i=0; i<attackerPositions.length; i++) {
			attackers[i] = network.new AttackNodeLayer(new Node(idGenerator, attackerPositions[i]));
		}
	}
	
	public void run() throws Exception {
		if(sources==null)
			throw new Exception("Must set source first");
		if(sink==null)
			throw new Exception("Must create sink first");

		//Set model parameters in corresponding classes
		Flood.MAX_TTL=routingParameters.maximumTTL;
		Physical.BIT_RATE=trafficParameters.channelBitRate;			
		DBR.MAXIMUM_FORWARDING_DELAY=routingParameters.maximumForwardingDelay;	
		DBR.FORWARDING_THRESHOLD=routingParameters.forwardingThreshold;
		DPR.FORWARDING_PROBABILITY = routingParameters.forwardingProbability;

		Random appRand = new Random(trafficParameters.trafficSeed);
		Application[] apps = new Application[sources.length];
		for(int i=0; i<apps.length; i++) {
			apps[i]= new Application( 
					0,
//					apps.length*appRand.nextDouble(), 
					sources[i], sink, 0.25, appRand.nextInt())
					.setNumberOfMessages(trafficParameters.messageCount)
					.setGenerationBitRate(trafficParameters.generationBitRate)
					.setMessageBitLength(trafficParameters.messageBitLength); 
		}
		
		if(debug) {
			Scheduler.it.enableDebug();
			network.enableDebug();
		}

		Scheduler.it.handleAll();

		if(debug) {
			physical.printStatistics();
			for(int i=0; i<apps.length; i++) {
				apps[i].printStatistics();
			}
		}
		
		results.numberOfTransmissions=new int[nodes.length];
		for(int i=0;i<nodes.length;i++) {
			results.numberOfTransmissions[i]=nodes[i].statistics.numberOfTransmissions;
		}
		
		results.forwardedCount = 0;
		for(NodeLayer node: nodes) {
			results.forwardedCount += node.getPhysical().statistics.sentCount;
		}
		results.sentCount=0;
		for(NodeLayer node: sources) {
			results.sentCount += node.getPhysical().statistics.sentCount;
		}

		results.generatedCount=0;
		results.receivedCount=0;
		for(int i=0; i<apps.length; i++) {
			results.generatedCount+=apps[i].generatedCount;
			results.receivedCount+=apps[i].receivedCount;
			results.averageE2EDelay +=apps[i].getAverageEndToEndDelay();
		}
		results.averageE2EDelay /= apps.length;
		results.forwardedCount-=results.generatedCount;
		
		results.deliveryRatio = new double[sources.length];
		results.e2eDelay = new double[sources.length];
		for(int i=0; i<sources.length; i++) {
			results.deliveryRatio[i] = (double)(apps[i].receivedCount)
									/(double)(apps[i].generatedCount);
			results.e2eDelay[i] = apps[i].getAverageEndToEndDelay();
		}
		
		results.finishTime = Scheduler.it.now();
	}
	
	public void updateForwardingProbability() {
		for(int i=0; i<topology.sourceNodes.length; i++) {
			routingParameters.forwardingProbability[topology.sourceNodes[i]] = 
				routingParameters.K * (
						routingParameters.forwardingProbability[topology.sourceNodes[i]] - results.deliveryRatio[i] + 1);
		}
	}
	
	/**Example using Simulator object
	 * @param args
	 * @throws InterruptedException 
	 * @throws IllegalStateException 
	 * @throws IllegalArgumentException 
	 * @throws ExecutionException 
	 * @throws RejectedExecutionException 
	 */
	public static void main(String[] args ) throws IllegalArgumentException, IllegalStateException, InterruptedException, RejectedExecutionException, ExecutionException {

		Simulator sim = new Simulator();
		sim.topologyParameters.topologySeed = 17943;
		sim.topologyParameters.networkDensity = 3;
		sim.topologyParameters.deploymentSideLength=500;
		sim.topologyParameters.deploymentDepth=250;
		sim.topologyParameters.sourcesCount = 1;
		sim.topologyParameters.attackersCount =1;

		sim.trafficParameters.messageCount=10;
		sim.trafficParameters.generationBitRate=1;
		sim.trafficParameters.messageBitLength=100;
		sim.trafficParameters.channelBitRate=100000;
		sim.trafficParameters.framingBitOverhead=160;

		sim.routingParameters.K = 0.8;
		sim.routingParameters.forwardingProbability = new double[256];
		Arrays.fill(sim.routingParameters.forwardingProbability, 0.8);
		sim.routingParameters.forwardingThreshold = 0;
		sim.routingParameters.maximumForwardingDelay=1000;
		
		try {
			sim.debug=true;
			sim.topology.sinkPos=new double[]
				{250, 250, 0};
			sim.topology.nodePositions = new double[][]{
				{250,250,0},
				{144.93236342234,197.108960356535,16.5855567737351},
				{203.295819690989,317.76941185012,89.0652249256978},
				{261.063293152552,274.007180403145,105.556680348683},
				{158.23655278988,233.837341794375,179.58811647348},
				{168.306871948773,176.2345376969,223.331846382113},
				{363.406238089536,204.65417837108,137.295722962193},
				{185.060935307332,91.1040989783746,84.0628171552127},
				{37.0073150843427,132.961361804754,34.8180079326821},
				{362.355873292609,300.708853501805,236.607883105458},
				{352.615563245572,152.311323356518,142.460088090444},
				{269.771444384837,112.087630558154,240.38536108136},
				{134.221010431667,400.773586838063,172.541141718101},
				{134.372446398112,60.6148830838347,123.94769064553},
				{47.3263778530943,162.688480366709,120.10065221854},
				{307.773267705585,137.341562680081,195.10604839079},
				{20.1759502619379,435.986770826845,224.335796357047},
				{203.744650880151,10.6967135988879,170.210177401343},
				{131.621366650576,11.5634295939092,243.544389092804},
				{50.1781187217191,434.084366716015,183.088805216717},
				{35.4502827542157,416.140366156516,221.275988477404},
				{363.432013117747,345.969359908596,6.79980721325205},
				{400.088339158508,371.09813557361,66.1812618449161},
				{393.550915854219,460.402981944853,148.115013878514},
				{455.154859585989,118.771558520916,230.742792069311},
				{477.108255213613,487.385051936932,212.944844689929},
				{7.59838895216775,321.339120309239,195.538835999059},
				{494.836021659927,294.818209603396,18.0586900829593},
				{430.001653012442,353.956910358624,192.243524433513},
				{468.150398767779,162.530004520884,31.6822393558034},
				{451.478485327843,31.6549316713198,38.6887703946467}
				};
				
				sim.topology.attackerPositions=new double[][] {
			   {131.6213666505755,   11.5634295939092,    244.5443890928039}
			};
			sim.topology.sourceNodes = new int[] {18};

			sim.topology.createModel();
			sim.run();
			
			System.out.println("Testing " + Simulator.class);
			System.out.println("p\tattacks\tnodes\tdegree\tgenrtd\tsent\trecved\tforwd\te2eDly\tfinTime");
			System.out.println(String.format("%4.2f", DPR.FORWARDING_PROBABILITY[sim.topology.sourceNodes[0]]) +"\t" + 
					sim.topologyParameters.attackersCount + "\t" +
					sim.topology.nodeCount + "\t" + 
					sim.topology.countSourceNeighbors(sim.topology.sourceNodes[0])+ "\t" + 
					sim.results.generatedCount + "\t" + sim.results.sentCount+"\t" + 
					sim.results.receivedCount + "\t" + sim.results.forwardedCount +"\t" + 
					String.format("%6.2f",sim.results.averageE2EDelay) + "\t" + 
					String.format("%6.2f",sim.results.finishTime) );
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}		
}
