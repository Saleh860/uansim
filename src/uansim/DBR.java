package uansim;

import java.util.Arrays;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.function.BiFunction;
import java.util.function.Consumer;

import uansim.Physical.NodeLayer;

public class DBR extends Flood{
	public static double COMMUNICATION_RANGE=150, INTERFERENCE_RANGE=250; // (m)
	public static final double PROPAGATION_SPEED=1500;	//(m/s)
	/**
	 * Depth measurement is quantized to be encoded in packet headers. 
	 */
	public static double DEPTH_UNIT = 4; // (m)
	
	/**
	 * The forwarding delay when the depth gain is minimum
	 */
	public static double MAXIMUM_FORWARDING_DELAY = 1; // (s)
	
	/**
	 * The minimum depth gain. Usually a non-negative number.
	 */
	public static double FORWARDING_THRESHOLD = 0; // (m)

	public DBR(Physical physical) {
		super(physical);
		
	}
	
	public static class Packet extends Flood.Packet{
		
		double senderDepth;
		
		public double getSenderDepth() {
			return senderDepth;
		}

		@Override
		public void init(Network.Layer source, Network.Layer destination, byte[] message) {
			super.init(source, destination, message);
			assert(source instanceof NodeLayer);
			this.senderDepth = ((NodeLayer) source).getDepth();
		}
		
		@Override
		public void update(Network.Layer relay) {
			super.update(relay);
			this.senderDepth = ((NodeLayer)relay).getDepth();
		}
		@Override
		public void init(Network.Layer source, byte[] message) {
			super.init(source, message);
			this.senderDepth = ((NodeLayer) source).getDepth();
		}

		public Packet() {
			
		}
		
		@Override
		int getHeaderLength() {
			return 5;
		}
		
		@Override
		public Packet fromByteArray(byte[] pdu) {
			super.fromByteArray(pdu);
			this.senderDepth = (pdu[4]&0xFF)*DEPTH_UNIT;
			return this;
		}
		
		@Override
		public String toString() {
			return super.toString() +",sndrDpth="+senderDepth;
		}

		@Override
		byte[] writeHeader(byte[] header) {
			// TODO Auto-generated method stub
			header[4]=(byte) (this.senderDepth/DEPTH_UNIT);
			return super.writeHeader(header);
		}
		

	}

	public class NodeLayer extends Flood.NodeLayer {
		
		public Consumer<Packet> delayedForwardDecision;
		public BiFunction<Packet, Double, Void> updateDelayedDecision;
		public Consumer<Packet> cancelDelayedDecision;
		
		public NodeLayer(Node node) {
			super(physical.new NodeLayer(node));
			initNodeLayer();
		}

		public NodeLayer(Physical.NodeLayer node) {
			super(node);
			initNodeLayer();
		}

		private void initNodeLayer() {
			this.newPacket = ()->new Packet();
			
			this.forwardDecision = (netpacket) -> {
				Packet packet = (Packet) netpacket;

				double forwardingDelay = calculateForwardingDelay(packet.getSenderDepth(), this.getDepth());

				// if forwardingDelay>=0 the packet may be forwarded if fresh
				if(forwardingDelay>=0) {
					delayedForwardingQueue.add(packet, now()+forwardingDelay);
					log("--->     Delayed forwarding time = " + String.format("%8.4f", now()+forwardingDelay));
				}
				else {
					log("--->     Dropped because depth gain < threshold");
				}
			};

			this.delayedForwardDecision = (packet) -> {
				this.forward.accept(packet);
			};

			
			this.updateDelayedDecision = (packet, newForwardingTime) -> {				
				delayedForwardingQueue.update(packet, newForwardingTime);	
				log("Update forwarding time to " + newForwardingTime);
				return null;
			};
			
			this.cancelDelayedDecision = (packet)-> {
				delayedForwardingQueue.remove(packet);
				log("Delayed forwarding canceled");
			};
						
			Consumer<uansim.Network.Packet> superDuplicateDecision = super.duplicateDecision;
			this.duplicateDecision = (netpacket) -> {
				Packet packet = (Packet) netpacket;
				superDuplicateDecision.accept(packet);

				//Duplicate packets can change forwarding time
				if(delayedForwardingQueue.contains(packet)) {
					double newForwardingDelay = calculateForwardingDelay(packet.getSenderDepth(), this.getDepth());

					// if the sender depth is smaller than this node, 
					// the delayed packet should be dropped 
					if(newForwardingDelay<0) {
						cancelDelayedDecision.accept(packet);
					}
					// Attempt to update the forwarding time of the delayed packet
					else {
						double oldForwardingTime = delayedForwardingQueue.getForwardingTime(packet);

						if(oldForwardingTime>now()+newForwardingDelay) {
							updateDelayedDecision.apply(packet, now()+newForwardingDelay);
						}
						else {
							
						}
					}
				}
			};
			
			
		}

		public double calculateForwardingDelay(double senderDepth, double nodeDepth) {
//			double tau = Physical.COMMUNICATION_RANGE / Physical.PROPAGATION_SPEED;
			double d = senderDepth-nodeDepth;
			
			if(d>COMMUNICATION_RANGE) {
				log("    ---> This shouldn't happen!");
				return -1;
			}
			
			//double delta = 2*tau / MAXIMUM_FORWARDING_DELAY * Routing.COMMUNICATION_RANGE;
//			double delta = 2*Physical.COMMUNICATION_RANGE * Physical.COMMUNICATION_RANGE
//					/ Physical.PROPAGATION_SPEED / MAXIMUM_FORWARDING_DELAY;
			if(d>FORWARDING_THRESHOLD) {
//				return (2*tau/delta)*(Routing.COMMUNICATION_RANGE-d);				
				return MAXIMUM_FORWARDING_DELAY *(COMMUNICATION_RANGE-d)/ COMMUNICATION_RANGE;
			}
			else {
				return -1;
			}
		}

		/**Buffer for packets kept for delayed forwarding 
		 * 
		 *
		 */
		class DelayedForwardingQueue {
			PriorityQueue<Record> queue= new PriorityQueue<Record>();
			HashMap<Packet,Record> map = new HashMap<Packet,Record>();				
			Scheduler.Event nextEvent = null;
			
			class Record implements Comparable<Record>{
				public Packet packet;
				public byte[] sdu;
				public double forwardingTime;
				public Record(Packet packet, double forwardingTime) {
					this.packet=packet;
					this.forwardingTime=forwardingTime;
				}
				@Override
				public int compareTo(Record o) {
					return Double.compare(forwardingTime, o.forwardingTime);
				}
			}
			
			Record makeRecord(Packet packet, double forwardingTime) {
				return new Record(packet,forwardingTime);
			}
			
			public void add(Packet packet, double forwardingTime) {
				if(!map.containsKey(packet)) {
					Record record = makeRecord(packet, forwardingTime);
					map.put(packet, record);
					queue.add(record);
					if(queue.peek()==record) {
						updateNextEvent();
					}
				}
				else {
					log("     ---> packet already in the delayed forwarding queue");
				}
			}

			/**
			 * Schedule an event for the queue head, canceling any prior events
			 */
			private void updateNextEvent() {
				//Remove old event (if any)
				if(nextEvent!=null)
					nextEvent.cancel();
				
				nextEvent=null;

				//Schedule new event (if any)
				Record record = queue.peek();
				if(record!=null) {
					nextEvent = Scheduler.it.schedule(
							record.forwardingTime - now(), 
							(s)->forwardQueueHead(), ()->"Delayed-forwarding timeout");
				}
			}

			/**
			 * Forward the packet at the top of the queue and schedule the next event
			 */
			private void forwardQueueHead() {
				Record record = queue.peek();
				if(record!=null) {
					NodeLayer.this.delayedForwardDecision.accept(record.packet);
					queue.remove(record);
					map.remove(record.packet);
				}
				updateNextEvent();
			}

			/**Delete the packet with the given header from the queue
			 * @param header
			 */
			public void remove(Packet packet) {
				Record record = map.get(packet);
				queue.remove(record);
				map.remove(packet);
				updateNextEvent();
			}

			public boolean contains(Packet packet) {
				return map.containsKey(packet);
			}

			public void update(Packet packet, double newForwardingTime) {
				Record record = map.get(packet);
				queue.remove(record);
				record.forwardingTime = newForwardingTime;
				queue.add(record);
				updateNextEvent();
			}

			public double getForwardingTime(Packet packet) {
				Record record = map.get(packet);
				return record.forwardingTime;
			}
		}
		
		DelayedForwardingQueue delayedForwardingQueue = new DelayedForwardingQueue();
		
	}
	
	public class AttackNodeLayer extends NodeLayer {
		private Consumer<uansim.Network.Packet> attackForewardDecision =  (packet) -> {
			packet.update(this);
			Scheduler.it.schedule(Network.PROCESSING_DELAY, 
					(s)->forward.accept(packet), 
					()->"Node#" + String.format("%3d", getId()) +"::SpoofDepth("+ packet +")");	
			log("--->     Spoofed depth ");
		};

		public AttackNodeLayer(Node node) {
			super(node);
			forwardDecision = attackForewardDecision ;
		}

		public AttackNodeLayer(Physical.NodeLayer node) {
			super(node);
			forwardDecision = attackForewardDecision ;
		}

		@Override
		public double getDepth() {
			return (getPhysical().getPos(2) > COMMUNICATION_RANGE)?
					(getPhysical().getPos(2) - COMMUNICATION_RANGE):0;
		}
	}

	public static class Test	{
		
		public static void main(String[] args) {
			Channel channel = new Channel();
			Physical physical = new Physical(channel);
			DBR network = new DBR(physical);

			NodeLayer source = network.new NodeLayer(physical.new NodeLayer(new Node(new double[]{200,200,100})));
			network.new NodeLayer(physical.new NodeLayer(new Node(new double[]{200,100,50})));
			network.new NodeLayer(physical.new NodeLayer(new Node(new double[]{100,220,40})));
			network.new NodeLayer(physical.new NodeLayer(new Node(new double[]{100,100,30})));
			NodeLayer sink = network.new NodeLayer(physical.new NodeLayer(new Node(new double[]{100,0,70})));
			
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
			
			for(int i=0; i<0; i++) {
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
}
