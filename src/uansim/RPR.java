package uansim;

import java.util.Arrays;
import java.util.Collection;
import java.util.Dictionary;
import java.util.HashMap;
import java.util.Objects;
import java.util.Random;
import java.util.function.Consumer;

import uansim.DBR.NodeLayer;
import uansim.Network.Packet;

/** Resilient pressure-based routing protocol 
 * Zuba 2014
 * @author saleh
 *
 */
public class RPR extends Network {
	
	public static double COMMUNICATION_RANGE=150, INTERFERENCE_RANGE=250; // (m)
	public static final double PROPAGATION_SPEED=1500;	//(m/s)

	/**
	 * Depth measurement is quantized to be encoded in packet headers. 
	 */
	public static double DEPTH_UNIT = 4; // (m)
	
	/**
	 * Maximum depth of the deployment
	 */
	public static double MAX_DEPTH=250;
	
	/**
	 * Packet lifetime and end-to-end delay time unit
	 */
	public static double DELAY_UNIT = 0.1; // (s)
	
	/**
	 * Maximum number of holding time units (r2)
	 */
	public static int MAX_HOLDING = 10;

	/**
	 * The maximum number of forwarding hops before dropping a packet to avoids loops.
	 */
	public static int MAX_TTL = 5;
	
	/**
	 * Maximum number of retransmissions until an I-ACK is received
	 */
	public static int MAX_RETRANSMIT_COUNT=3;
	
	/**
	 * Number of transmission rounds of the discovery phase
	 */
	public static int DISCOVERY_ROUNDS=5;

	public RPR(Physical physical) {
		super(physical);
	}
	
	public RPR(Physical physical, Random random) {
		this(physical);
		rand = random;
	}

	Random rand = new Random();	

	/**RPR Packet format
	 * Source ID   	- 2 bytes
	 * Forwarder ID - 2 bytes
	 * Seq. No.		- 2 bytes
	 * Depth        - 1 byte
	 * Tmin			- 1 byte
	 * Tmax 		- 1 byte
	 * TTL			- 1 byte
	 * 
	 */
	static class Packet implements Network.Packet {
		int sourceId;
		int forwarderId;
		int seqNum;
		double forwarderDepth;
		double tMin;
		double tMax;
		double pT;		//Packet lifetime
		int TTL;
		byte[] sdu;
		double Hn;
		int transmissionCount=0;
				
		@Override
		public uansim.Network.Packet fromByteArray(byte[] pdu) {
			this.sourceId = ((pdu[0]&0xFF) | ((pdu[1]&0xFF)<<8));
			this.forwarderId = ((pdu[2]&0xFF) | ((pdu[3]&0xFF)<<8));
			this.seqNum= ((pdu[4]&0xFF) | ((pdu[5]&0xFF)<<8));
			this.forwarderDepth = (pdu[6]&0xFF) * DEPTH_UNIT;
			this.tMin = (pdu[7]&0xFF) * DEPTH_UNIT;
			this.tMax = (pdu[8]&0xFF) * DEPTH_UNIT;
			this.pT = (pdu[9]&0xFF) * DELAY_UNIT;
			this.TTL = pdu[10] & 0xFF;
			this.sdu = Arrays.copyOfRange(pdu, 11, pdu.length);
			return this;
		}
		
		@Override
		public byte[] toByteArray() {
			byte[] pdu = new byte[11+sdu.length];
			pdu[0] = (byte) (this.sourceId & 0xFF);
			pdu[1] = (byte) ((this.sourceId & 0xFFFF) >> 8);
			pdu[2] = (byte) (this.forwarderId & 0xFF);
			pdu[3] = (byte) ((this.forwarderId & 0xFFFF) >> 8);
			pdu[4] = (byte) (this.seqNum & 0xFF);
			pdu[5] = (byte) ((this.seqNum & 0xFFFF) >> 8);
			pdu[6] = (byte) (this.forwarderDepth/DEPTH_UNIT);
			pdu[7] = (byte) (this.tMin/DEPTH_UNIT);
			pdu[8] = (byte) (this.tMax/DEPTH_UNIT);
			pdu[9] = (byte) (this.pT/DELAY_UNIT);
			pdu[10] = (byte) this.TTL; 
			System.arraycopy(sdu, 0, pdu, 11, sdu.length); 
			return pdu;
		}
		
		@Override
		public byte[] getSDU() {
			return sdu;
		}
		
		@Override
		public int getDestination() {
			return 0;
		}
		
		@Override
		public int getSource() {
			return this.sourceId;
		}
		
		@Override
		public boolean isRoutable() {
			return this.TTL>0;
		}
		
		@Override
		public void init(Network.Layer source, Network.Layer destination, byte[] message) {
			init(source, message);
		}
		
		@Override
		public void init(Network.Layer source, byte[] message) {
			this.sourceId = source.getId();
			this.forwarderId = source.getId();
			this.seqNum = source.nextSequenceNumber();
			this.forwarderDepth = ((Layer) source).getDepth();
			this.tMin=((Layer)source).tMin;
			this.tMax=((Layer)source).tMax;
			this.TTL=MAX_TTL;
			this.sdu=message.clone();
		}
		
		@Override
		public void update(Network.Layer relay) {
			if(transmissionCount==0) {
				this.TTL--;
				this.forwarderDepth=relay.getDepth();
			}
			this.pT+=Hn;
			this.forwarderId=relay.getId();
			transmissionCount++;
		}
		
		@Override
		public String toString() {
			return "RPR(SourceID="+sourceId+", ForwarderID="+forwarderId+
					", SeqNo="+seqNum+", senderDepth="+forwarderDepth+", Tmin="+tMin
					+", Tmax="+tMax+", Pt="+pT+", TTL="+TTL+", Hn="+Hn+
					", transmissionCount="+transmissionCount+", sdu="+Arrays.toString(sdu)+")";
		}
		/**
		 *A packet is summarized by the hash of its source ID and sequence number
		 */
		
		/**
		 *A packet is summarized by the hash of its source ID and sequence number
		 */
		@Override
		public int hashCode() {
			return Objects.hash(sourceId,seqNum);
		}

		/**
		 *Two packets are considered equal if they at least have the same 
		 *source ID and sequence numbers
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (!(obj instanceof Packet))
				return false;
			Packet other = (Packet) obj;
			return seqNum == other.seqNum  && sourceId == other.sourceId;
		}
		
		
	}
	
	
	public class Layer extends Network.Layer {
		HashMap<Integer, Double> neighbors = new HashMap<Integer,Double>();
		PacketCache holdQueue = new PacketCache();
		
		public double tMin;
		public double tMax;
		public int transmissionRounds=0;
				
		public boolean inDiscovery() {
			return transmissionRounds<DISCOVERY_ROUNDS;
		}
		boolean isSink=false;
		boolean isRoutingVoid=false;
		
		public Layer(Node node) {
			super(node);
		}
		
		public Layer(Physical.NodeLayer phy) {
			super(phy);
		}
				
		public void initNodeLayer(Physical.NodeLayer phy) {
			super.initNodeLayer(phy);
			tMin=this.getDepth();
			tMax=Math.max(0,this.getDepth()-COMMUNICATION_RANGE);
			this.getPhysical().actionReceiveCompleted = (data)->listenerReceiveCompleted(data); 
			this.newPacket = ()-> new Packet();
			this.forwardDecision = (packet)->listenerHold((Packet)packet);
		}
		
		public void setSink() {
			isSink=true;
		}

		// Actions
		Consumer<Void> actionRoutingVoidDetected = (V) -> {
			log("no neighbors disovered - routing void detected ");
			isRoutingVoid=true;
			clearHoldQueue();				
		};
		
		Consumer<Packet> actionForwardHeldPacket = (packet)->forwardHeldPacket(packet);
		
		void forwardHeldPacket(Packet packet){
			//Check if the packet is still in the hold queue
			if(! holdQueue.contains(packet) || holdQueue.get(packet)<0) {
				if(packet.transmissionCount==0)
					log("packet dropped");
				else
					log("packet acknowledged successfully.");
			}
			else {
				if(inDiscovery()) { 
					//Packet being transmitted for the first time 
					if((packet).transmissionCount==0) {
						packet.tMin=tMin;
						packet.tMax=tMax;
						forward.accept(packet);		// will call update which increments packet.transmissionCount 
						transmissionRounds++;
					}
					else {
						log("no I-ACK received - routing void detected");
						isRoutingVoid=true;
						clearHoldQueue();
					}
				}
				else if(!inDiscovery()) {
					if(isRoutingVoid) {
						clearHoldQueue();
					}
					else {
						if(neighbors.size()==0) {
							log("no neighbors discovered - routing void detected - flushing hold queue");
							isRoutingVoid=true;
							clearHoldQueue();						
						}
						else {
							Object[] neighborDepths = neighbors.values().toArray();
							//Set threshold window
							tMin=(double) neighborDepths[rand.nextInt(neighborDepths.length)];
							tMax=(double) neighborDepths[rand.nextInt(neighborDepths.length)];
							if(tMax>tMin) {
								double t=tMin;
								tMin=tMax;
								tMax=t;
							}
							((RPR.Packet)packet).tMin=tMin;
							((Packet)packet).tMax=tMax;
							forward.accept(packet);
							if(((Packet)packet).transmissionCount<=MAX_RETRANSMIT_COUNT) {
								log("updating threshold window, tMin="+tMin+", tMax="+tMax+", chosen from "+
										neighbors.keySet()+" neighbors");
								double waitForACKtime = 2*COMMUNICATION_RANGE/PROPAGATION_SPEED
										+ MAX_HOLDING*DELAY_UNIT;
								((Packet)packet).Hn = waitForACKtime;
								Scheduler.it.schedule(waitForACKtime, 
										(s)->forwardHeldPacket(packet), 
										()->"Node#" + String.format("%3d", getId()) +"::RetransmitPacket("+ packet +")");				
							}
							else {
								holdQueue.remove(packet);
							}
						}
					}
				}
			}
		}
		
		Consumer<Packet> actionHold = (packet) -> listenerHold(packet);
		Consumer<Packet> actionDuplicateReceived = (packet)->listenerDuplicateReceived(packet); 
		
		void listenerHold(Packet packet) {
			holdQueue.put(packet, now());
			double Hn = rand.nextInt(MAX_HOLDING+1)*DELAY_UNIT;
			double AVGe2e = MAX_DEPTH/PROPAGATION_SPEED
					+(0+MAX_HOLDING)/2.0*MAX_DEPTH/COMMUNICATION_RANGE;
			if(packet.pT>AVGe2e) {
				if(this.getDepth()>2*COMMUNICATION_RANGE)
					Hn/=2.0;
				else
					Hn/=3.0;
			}
			packet.Hn = Hn;
			log("Node#" + String.format("%3d", getId()) +" hold time = " + Hn);
			Scheduler.it.schedule(Hn, 
					(s)->forwardHeldPacket(packet), 
					()->"Node#" + String.format("%3d", getId()) +"::ForwardHeldPacket("+ packet +")");				
		};
		
		private void listenerDuplicateReceived(Packet packet) {
			log("should drop (duplicate).");
			/*********************************/
			duplicateDecision.accept(packet);
			/*********************************/

			if(holdQueue.contains(packet)) {
				if(packet.forwarderDepth<this.getDepth()) {
					log("found better candidate, held packet dropped.");
					holdQueue.remove(packet);
				}
				else {
					log("no better candidate found, held packet retained.");
				}
			}
			else {
				log("not found in held queue.");
			}
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
		void listenerReceiveCompleted(byte[] data) {
			if(data.length==0)
				log("received noise");
			else {
				if(!isSink && !isRoutingVoid && !inDiscovery() && neighbors.size()==0) {
					actionRoutingVoidDetected.accept(null);
				}
				else if(isRoutingVoid) {
					log("should drop (routing void).");					
				}
				else {
					Packet packet = (Packet) unpacker.apply(data);
					boolean isFresh = !packetCache.contains(packet) || 
							(now() - packetCache.get(packet))>CACHE_TIMEOUT;
							
					log((isFresh?"fresh":"duplicate")+" packet received ("+(inDiscovery()?"discovery ":"")+
							(isRoutingVoid?"void":"")+",neighbors="+neighbors.size()+") " + packet);

					//Update neighbors tables
					if(!isSink 
							&& (this.getDepth() >= packet.forwarderDepth) // neighbor closer to the surface
							&& this.getDepth()-packet.forwarderDepth<=COMMUNICATION_RANGE //neigbor not further than communication range
							) {
						neighbors.put(packet.forwarderId,packet.forwarderDepth);
					}
					
					if(!isFresh) {
						actionDuplicateReceived.accept(packet);
					}
					//Handle fresh packets
					else {
						if(isSink) {
							log("should deliver");
							/*********************************/
							deliverDecision.accept(packet);
							/*********************************/						
							//Send acknowledgement
							packet.tMin=0;
							packet.tMax=0;						
							Scheduler.it.schedule(PROCESSING_DELAY, 
									(event)->{this.forward.accept(packet);}
									, ()->{return "Sink acknowledgement";});
							//Update cache
							packetCache.put(packet,now());
						}
						else {
							//Check constraints
							if(this.getDepth()-packet.tMin<=DEPTH_UNIT && packet.tMax-this.getDepth()<=DEPTH_UNIT) {
	
								//Do not forward packets whose TTL expired 
								if(!packet.isRoutable()) {
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
									//Update cache
									packetCache.put(packet,now());
								}
							}	
							else {
								log("should drop (depth outside threshold: " + Math.floor(this.getDepth()) + ">" + packet.tMin +
										" || " + Math.ceil(this.getDepth()) + "<" + packet.tMax + " )");								
							}
						}
					}
				}
			}
		};

		private void clearHoldQueue() {
			log("flushing hold queue");
			holdQueue.clear();			
		};

	}

	public static class Test	{
		
		public static void main(String[] args) {
			Channel channel = new Channel();
			Physical physical = new Physical(channel);
			RPR network = new RPR(physical);

			Layer source = network.new Layer(physical.new NodeLayer(new Node(new double[]{200,200,100})));
			network.new Layer(physical.new NodeLayer(new Node(new double[]{200,100,50})));
			network.new Layer(physical.new NodeLayer(new Node(new double[]{100,220,40})));
			Layer sink = network.new Layer(physical.new NodeLayer(new Node(new double[]{100,100,0})));
			network.new Layer(physical.new NodeLayer(new Node(new double[]{100,0,70})));
			sink.setSink();
			
			if(true){
				byte[] message = new byte[] {0,1,2,3};
				Scheduler.it.schedule(0, 
						(s)->source.sendTo(message, sink), 
						()->"Node#" + source.getId() +"::SendMessage("+ Arrays.toString(message) +")");
			}
			
			
			if(false){
				byte[] message = new byte[] {4,5,6,7};
				Scheduler.it.schedule(1, 
					(s)->source.sendTo(message, sink), 
					()->"Node#" + source.getId() +"::SendMessage("+ Arrays.toString(message) +")");

			}
			
			
			if(false){
				byte[] message = new byte[] {8,9,10,3};
				Scheduler.it.schedule(10, 
					(s)->source.broadcast(message), 
					()->"Node#" + source.getId() +"::BroadcastMessage("+ Arrays.toString(message) +")");
			}
			
			for(int i=0; i<10; i++) {
				byte[] message = new byte[] {(byte)i,(byte)(i+1),(byte)(i+2),(byte)(i+3)};
				Scheduler.it.schedule(2*i, 
						(s)->source.sendTo(message, sink), 
						()->"Node#" + source.getId() +"::SendMessage("+ Arrays.toString(message) +")");				
			}
			Scheduler.it.enableDebug();
			network.debug=true;
			Scheduler.it.handleAll();
		}
	}
	
	public class AttackNodeLayer extends Layer {
		public AttackNodeLayer(Node node) {
			super(node);
		}

		public AttackNodeLayer(Physical.NodeLayer node) {
			super(node);

			forwardDecision = (packet) -> {
				Packet p = ((Packet)packet);
				p.update(this);
				p.forwarderDepth=0;
				p.tMin=0;
				p.tMax=0;
				Scheduler.it.schedule(Network.PROCESSING_DELAY, 
						(s)->forward.accept(packet), 
						()->"Node#" + String.format("%3d", getId()) +"::SpoofDepth("+ packet +")");				
			};
		}
		@Override
		public double getDepth() {
			return (getPhysical().getPos(2) > COMMUNICATION_RANGE)?
					(getPhysical().getPos(2) - COMMUNICATION_RANGE):0;
		}
	}


}
