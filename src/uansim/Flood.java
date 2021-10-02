package uansim;

import java.util.Arrays;
import java.util.Objects;

import uansim.Physical.NodeLayer;


public class Flood extends Network {
	/**
	 * The maximum number of forwarding hops before dropping a packet to avoids loops.
	 */
	public static int MAX_TTL = 5;
	
	public Flood(Physical physical) {
		super(physical);
	}
	
	public static class Packet implements Network.Packet{
		private int ttl;
		private int sourceId;
		private int destinationId;
		private int sequenceNumber;
		/**
		 * Application message
		 */
		private byte[] sdu;

		/**Fill in the packet contents
		 * @param ttl time-to-live
		 * @param sourceId source ID
		 * @param destinationId destination ID
		 * @param sequenceNumber source-unique sequence number
		 * @param sdu application message
		 * @return
		 */
		private Packet set(byte ttl, byte sourceId, byte destinationId, byte sequenceNumber, byte[] sdu) {
			this.ttl = ttl&0xFF;
			this.sourceId = sourceId&0xFF;
			this.destinationId = destinationId&0xFF;
			this.sequenceNumber = sequenceNumber&0xFF;
			this.sdu = sdu;
			return this;
		}
		
		/**
		 * Default constructor. <br/>
		 * set(.) should be called to initialize the packet before usage.
		 */
		public Packet() {	
		}
		
		@Override
		public Packet fromByteArray(byte[] pdu) {
			return set(pdu[0], pdu[1], pdu[2], pdu[3], Arrays.copyOfRange(pdu, getHeaderLength(), pdu.length));
		}
		
		@Override
		public byte[] toByteArray() {
			int l1 = getHeaderLength();
			int l2 = sdu.length;
			byte[] pdu = new byte[l1+l2];
			writeHeader(pdu);
			for(int i=0; i<l2; i++) {
				pdu[l1+i] = sdu[i];
			}
			return pdu;
		}
		
		/**
		 * @return application message
		 */
		@Override
		public byte[] getSDU() {
			return sdu;
		}
		
		@Override
		public int getDestination() {
			return destinationId;
		}
		
		@Override
		public int getSource() {
			return sourceId;
		}
		
		@Override
		public boolean isRoutable() {
			return ttl==0;
		}
		
		/**
		 * @return the packet header size in bytes
		 */
		int getHeaderLength() {
			return 4;
		}

		/**
		 * @param header buffer to store the returned packet header bytes
		 * @return packet header stored in byte array
		 */
		byte[] writeHeader(byte[] header) {
			header[0]=(byte)ttl;
			header[1]=(byte)sourceId;
			header[2]=(byte)destinationId;
			header[3]=(byte)sequenceNumber;
			return header;
		}
		
		/**Initialize a unicast packet fields given a source, destination and application message
		 * @param source unicast packet source node
		 * @param destination unicast packet destination node
		 * @param message application message
		 */
		@Override
		public void init(Network.Layer source, Network.Layer destination, byte[] message) {
			this.sourceId=source.getId();
			this.destinationId=destination.getId();
			this.sequenceNumber=source.nextSequenceNumber();
			this.ttl=MAX_TTL;			
			this.sdu=message.clone();
		}
		
		/**Initialize a broadcast packet fields given source and application message
		 * @param source
		 * @param message
		 */
		@Override
		public void init(Network.Layer source, byte[] message) {
			this.sourceId=source.getId();
			this.destinationId=BROADCAST_ID;
			this.sequenceNumber=source.nextSequenceNumber();
			this.ttl=MAX_TTL;		
			this.sdu=message.clone();
		}
		
		/**Update packet while being forwarded. Basically decreases TTL.
		 * @param relay the node currently forwarding the packet
		 */
		@Override
		public void update(Network.Layer relay) {
			this.ttl--;			
		}


		/**
		 *A packet is summarized by the hash of its source ID and sequence number
		 */
		@Override
		public int hashCode() {
			return Objects.hash(sequenceNumber, sourceId);
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
			return sequenceNumber == other.sequenceNumber && sourceId == other.sourceId;
		}
		
		@Override
		public String toString() {
			return "src="+sourceId+",seq="+sequenceNumber+
					",ttl="+ttl+",dst="+destinationId;
		}
	}
	
	class NodeLayer extends Network.Layer {

		public NodeLayer(Node node) {
			super(physical.new NodeLayer(node));
			initNodeLayer();
		}

		public NodeLayer(Physical.NodeLayer node) {
			super(node);
			initNodeLayer();
		}

		private void initNodeLayer() {
			newPacket = ()->new Packet();
		}
	}
	
	public static class Test{
		
		public static void main(String[] args) {
			Channel channel = new Channel();
			Physical physical = new Physical(channel);
			Flood network = new Flood(physical);

			Flood.NodeLayer source = network.new NodeLayer(new Node(new double[]{0,0,0}));
			NodeLayer relay = network.new NodeLayer(new Node(new double[]{100,0,0}));
			Flood.NodeLayer sink = network.new NodeLayer(new Node(new double[]{200,0,0}));
			
			source.getPhysical().getAdapter().enableDebug();
			relay.getPhysical().getAdapter().enableDebug();
			sink.getPhysical().getAdapter().enableDebug();
			
			if(true){
				byte[] message = new byte[] {0,1,2,3};
				Scheduler.it.schedule(0, 
						(s)->source.sendTo(message, sink), 
						()->"Node#" + source.getId() +"::SendMessage("+ Arrays.toString(message) +")");
			}
			
			
			if(true){
				byte[] message = new byte[] {4,5,6,7};
				Scheduler.it.schedule(1, 
					(s)->source.sendTo(message, sink), 
					()->"Node#" + source.getId() +"::SendMessage("+ Arrays.toString(message) +")");

			}
			
			
			if(true){
				byte[] message = new byte[] {8,9,10,3};
				Scheduler.it.schedule(10, 
					(s)->source.broadcast(message), 
					()->"Node#" + source.getId() +"::BroadcastMessage("+ Arrays.toString(message) +")");
			}
			
			for(int i=0; i<10; i++) {
				byte[] message = new byte[] {(byte)i,(byte)(i+1),(byte)(i+2),(byte)(i+3)};
				Scheduler.it.schedule(20+2*i, 
						(s)->source.sendTo(message, sink), 
						()->"Node#" + source.getId() +"::SendMessage("+ Arrays.toString(message) +")");				
			}
//			Scheduler.it.enableDebug();
			network.enableDebug();
 			network.physical.enableDebug();
			network.physical.channel.enableDebug();
			Scheduler.it.handleAll();
		}
	}
}
