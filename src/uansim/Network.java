package uansim;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Supplier;

import uansim.DPR.NodeLayer;


/**Routing base class
 * includes basic broadcasting and unicasting using TTL-limited flooding
 * @author saleh
 *
 */
public class Network {
	
	/**
	 * The duration in seconds after which a sequence number maybe reused by a source
	 */
	public static double CACHE_TIMEOUT = 10000; // (s)
	/**
	 * The delay of making the forwarding decision
	 */
	public static double PROCESSING_DELAY=0.001;
	/**
	 * Special node ID identifying broadcast packets intended to be delivered application layer in all nodes
	 */
	public static final int BROADCAST_ID=255;
	
	/**Basic routing layer packet format. Containing basic routing fields and application message
	 * @author saleh
	 *
	 */
	public static interface Packet{
		
		/**Initialize this packet content from the given PDU
		 * @param pdu protocol data unit
		 * @return packet with all fields set with values extracted from the given PDU
		 */
		public Packet fromByteArray(byte[] pdu); 
		
		/**convert packet to a byte array, to serve as an SDU to the lower layer
		 * @return a pdu with the packet fields embedded
		 */
		public byte[] toByteArray(); 	
		
		/**
		 * @return application message
		 */
		public byte[] getSDU(); 
		
		/**Get packet destination address
		 * @return destination address
		 */
		public int getDestination();

		/**Get packet source address
		 * @return source address
		 */
		public int getSource(); 
		
		/**Check whether a packet received for the first time should be forwarded
		 * @return
		 */
		public boolean isRoutable();
		
		/**Initialize a unicast packet fields given a source, destination and application message
		 * @param source unicast packet source node
		 * @param destination unicast packet destination node
		 * @param message application message
		 */
		public void init(Network.Layer source, Network.Layer destination, byte[] message); 
		
		/**Initialize a broadcast packet fields given source and application message
		 * @param source
		 * @param message
		 */
		public void init(Network.Layer source, byte[] message); 
		
		/**Update packet while being forwarded. Basically decreases TTL.
		 * @param relay the node currently forwarding the packet
		 */
		public void update(Network.Layer relay); 

		/**
		 *A packet is summarized by the hash of its source ID and sequence number
		 */
		@Override
		public int hashCode(); 

		/**
		 *Two packets are considered equal if they at least have the same 
		 *source ID and sequence numbers
		 */
		@Override
		public boolean equals(Object obj); 
		
		
		@Override
		public String toString();
	}
	
	
	/**
	 * The underlying physical layer
	 */
	Physical physical;
	public Physical getPhysical() {
		return physical;
	}
	
	ArrayList<Layer> nodes = new ArrayList<Layer>();
	
	/**
	 * Default constructor
	 */
	public Network(Physical physical) {
		this.physical = physical;
	}

	/**
	 * @param nodeId
	 * @return the node with the given nodeId or <b>null</b> if no such node exists.
	 */
	public Layer getNodeById(int nodeId) {
		for(Layer node: nodes)
			if(node.getId()==nodeId)
				return node;
		return null;
	}

	public ArrayList<Layer> getNodes() {
		return nodes;
	}
	
	/**Keeps track of packets seen so far, to avoid duplicate forwarding
	 * @author saleh
	 *
	 */
	static class PacketCache extends HashMap<Packet,Double> {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1414194002947732771L;

		/**check whether the given packet has been seen before
		 * @param packet a received packet to be checked for duplication
		 * @return
		 */
		boolean contains(Packet packet) {
			return super.containsKey(packet);
		}
		
		void put(Packet packet,double seenTime) {
			super.put(packet, seenTime);
		}
		
		/**Get the time at which a matching packet was first seen
		 * @param header
		 * @return
		 */
		Double get(Packet packet) {
			return super.get(packet);
		}
	}
	
	/**The network layer of a node
	 * @author saleh
	 *
	 */
	public class Layer {
		int seqNum=0;
		class Statistics {
			int numberOfTransmissions=0;
		}
		
		public Statistics statistics=new Statistics();
		
		public int nextSequenceNumber() {
			int next = seqNum;
			seqNum++;
			if(seqNum>=BROADCAST_ID)
				seqNum=0;
			return next;
		}
		
		/**
		 * supplier of new packets
		 */
		public Supplier<Packet> newPacket = ()-> null; 

		/**
		 * The underlying physical layer of the same node
		 */
		private Physical.NodeLayer phyLayer; 
		public Physical.NodeLayer getPhysical() {
			return phyLayer;
		}
		
		/**
		 *Packets already processed 
		 */
		PacketCache packetCache=new PacketCache();
		
		/**construct a network layer node and the underlying physical layer with default behavior
		 * @param node contains initial x, y and z coordinates of the node
		 */
		public Layer(Node node) {
			initNodeLayer(physical.new NodeLayer(node));
		}
		
		public Layer(Physical.NodeLayer node) {
			initNodeLayer(node);
		}
		
		void initNodeLayer(uansim.Physical.NodeLayer phyLayer) {
			nodes.add(this);

			this.phyLayer = phyLayer;
			this.phyLayer.actionReceiveAborted=  (n)-> log("collision detected");
			this.phyLayer.actionReceiveCompleted = (data)->actionReceiveCompletedDefault(data); 
		}

		/**
		 * convert packet to byte array
		 */
		public Function<Packet, byte[]> packer = (packet)-> packet.toByteArray(); //PDU.pack(packet.getHeader(), packet.getSDU()).toByteArray();
		public Function<byte[], Packet> unpacker = (data)->newPacket.get().fromByteArray(data);
		
		/**Delivers the received packet to the application
		 */
		public Consumer<byte[]> deliver = (message) -> 
			log("delivered message " + Arrays.toString(message));			

		/** Drops/discards/ignores the given received packet 
		 * @param header the packet header
		 * @param payload the application message
		 */
		public Consumer<Packet> drop = (packet) -> log("packet dropped");
		

		/** Forwards the given received packet on to neighbors
		 * @param header the DBR packet header
		 * @param payload the application message
		 */
		public Consumer<Packet> forward = (packet) -> {
			packet.update(this);
			if(phyLayer.send(packet.toByteArray())) {
				log("packet forwarded");
				statistics.numberOfTransmissions++;
			}else {
				log("packet forwarding failed");
			};
		};
		
		
		/** Called when a decision is made to forward the packet
		 * @param header the packet header
		 * @param payload the application message
		 */			
		public Consumer<Packet> forwardDecision = (packet) -> {
			Scheduler.it.schedule(PROCESSING_DELAY, 
					(s)->forward.accept(packet), 
					()->"Node#" + String.format("%3d", getId()) +"::ForwardPacket("+ packet +")");				
		};

		/** Called when a decision is made to discarded because its TTL expired
		 * @param header the packet header
		 * @param payload the application message
		 */			
		public Consumer<Packet> notRoutableDecision = (packet) -> {
			Scheduler.it.schedule(PROCESSING_DELAY, 
					(s)->drop.accept(packet), 
					()->"Node#" + String.format("%3d", getId()) +"::TTLExpired(" + packet + ")");					
		};

		public Consumer<Packet> deliverDecision = (packet) -> {
			Scheduler.it.schedule(PROCESSING_DELAY, 
					(s)->deliver.accept(packet.getSDU()), 
					()->"Node#" + String.format("%3d", getId()) +"::Deliver(" + packet + ")");				
		};

		/** Called when a decision is made to discarded because it is a duplicate
		 * @param header the packet header
		 * @param payload the application message
		 */			
		public Consumer<Packet> duplicateDecision = (packet) -> {
			//Duplicate, discard
			Scheduler.it.schedule(PROCESSING_DELAY, 
					(s)->drop.accept(packet), 
					()->"Node#" + String.format("%3d", getId()) +"::Duplicate(" + packet + ")");				
		};

		/**Transmit a unicast message
		 * @param message	the application message
		 * @param destination	the destination node (? extends Node)
		 */
		public boolean sendTo(byte[] message, Layer destination) {
			Packet packet = newPacket.get();
			packet.init(this, destination, message);
			packetCache.put(packet,now());

			byte[] pdu = packer.apply(packet);
			return phyLayer.send(pdu);
		}

		
		/**Transmit a broadcast message
		 * @param message the application message to be received by all nodes 
		 */
		public boolean broadcast(byte[] message) {
			Packet packet = newPacket.get();
			packet.init(this, message);
			packetCache.put(packet,now());
			byte[] pdu = packer.apply(packet);
			return phyLayer.send(pdu);
		}

		/**
		 * Four possible outcomes:
		 * <ol>
		 * <li>Deliver</li>
		 * <li>DiscardDuplicate</li>
		 * <li>DiscardExpiredTTL</li>
		 * <li>Forward</li>
		 * </ol> 
		 */				
		private void actionReceiveCompletedDefault(byte[] data) {
			if(data.length==0)
				log("received noise");
			
			else {
				Packet packet = unpacker.apply(data);
				log("received packet " + packet);

				//Drop packets that already exist in the cache
				if(packetCache.contains(packet) && 
						(now() - packetCache.get(packet))<CACHE_TIMEOUT) {
					log("should drop (duplicate)");
					/*********************************/
					duplicateDecision.accept(packet);
					/*********************************/						
				}
				//Handle fresh packets
				else {
					// Deliver packets destined to this node and broadcast messages
					if(packet.getDestination()==this.getId() //Unicast sent to this node
							|| (packet.getDestination()==BROADCAST_ID && //Broadcast
								packet.getSource() !=this.getId()		//  not sent from this node
								)) {
						log("should deliver");
						/*********************************/
						deliverDecision.accept(packet);
						/*********************************/
					}

					//Forward packets that are not destined specifically to this node 
					if(packet.getDestination()!=this.getId()) {
						//Do not forward packets whose TTL expired 
						if(packet.isRoutable()) {
							log("should drop (packet not routable)");
							/*********************************/
							notRoutableDecision.accept(packet);
							/*********************************/
						}
						else {
							log("should forward");
							/*********************************/
							forwardDecision.accept(packet);
							/*********************************/
						}
					}
				}
				//Refresh cache
				packetCache.put(packet,now());
			}
		};
		

		public int getId() {
			return phyLayer.getId();
		}

		public void log(String str) {
			if(debug) {
				System.out.println(String.format("t=%8.4f, Node#%3d %s\t"+str,
					now(), getId(), this.getClass().getName()));
				LogFile.logFile.println(String.format("t=%8.4f, Node#%3d %s\t"+str,
						now(), getId(), this.getClass().getName()));
			}
		}

		public double getDepth() {
			return getPhysical().getPos(2);
		}
	}
	
	public static Double now() {
		return Scheduler.it.now();
	}

	static boolean debug=false;
	public void enableDebug() {
		debug=true;
	}	
}
