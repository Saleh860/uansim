package uansim;

import java.util.Arrays;
import java.util.function.Consumer;
import java.util.function.Function;

public class Physical {
	
	public Channel channel;

	public static double BIT_RATE=1000;  // (bps)

	//TODO move to a MAC layer
	public static void setFramingBitOverhead(int framingBitOverhead) {
		Physical.framingBitOverhead = new byte[(framingBitOverhead+7)/8];
	}
//	private static int FRAMING_OVERHEAD=0;	//(Bytes)

	public static class Statistics {
		public int sentCount=0;
		public int receivedCount=0;
		public int receivedNoiseCount=0;
		public int collisionCount=0;
	}
	public Statistics statistics = new Statistics();
	
	public static byte[] framingBitOverhead=new byte[0];
	
	public static class PDU{
		int headerLength;
		int sduLength;
		private byte[] data;
		private byte[] header;
		private byte[] sdu;
		private static byte[] emptyHeader = new byte[0];

		private PDU() {
			
		}
		public PDU(byte[] data) {
			
		}
		public static PDU pack(byte[] sdu) {
			return pack(emptyHeader,sdu);
		}
		public static PDU pack(byte[] header,byte[] sdu) {
			PDU pdu = new PDU();
			if(header==null) header = emptyHeader;
			if(sdu==null) sdu = emptyHeader;
			pdu.header=header.clone();
			pdu.sdu=sdu.clone();
			pdu.headerLength=header.length;
			pdu.sduLength=sdu.length;
			pdu.data = new byte[header.length+sdu.length];
			for(int i=0; i<header.length; i++)
				pdu.data[i] = header[i];
			for(int i=0; i<sdu.length; i++)
				pdu.data[header.length+i] = sdu[i];
			return pdu;
		}
		public static PDU unpack(byte[] data, int headerLength) {
			PDU pdu = new PDU();
			if(data==null) {
				pdu.data = emptyHeader;
				pdu.header=emptyHeader;
				pdu.sdu=emptyHeader;
				pdu.headerLength=0;
				pdu.sduLength=0;
			}
			else {
				pdu.data=data;
				pdu.headerLength=headerLength;
				pdu.sduLength=data.length - headerLength;
				pdu.header = new byte[headerLength];
				pdu.sdu = new byte[pdu.sduLength];
			}			
			for(int i=0; i<headerLength; i++)
				pdu.header[i] = data[i];
			for(int i=0; i<pdu.sduLength; i++)
				pdu.sdu[i] = data[headerLength+i];
			return pdu;
		}
		public final byte[] toByteArray() {
			return data.clone();
		}
		public int length() {
			return data.length;
		}
		public byte[] getHeader() {
			return header.clone();
		}
		public byte[] getSDU() {
			return sdu.clone();
		}
	}

	public class NodeLayer {		

		public Statistics statistics = new Statistics();
		
		Adapter adapter;


		/**Construct a physical node and connect it to the physical channel
		 * @param node2 coordinates (x, y, z) where the depth is non-negative underwater
		 */
		public NodeLayer(Node node2) {
			adapter = new Adapter(node2);
			channel.connect(adapter);
			adapter.actionReceiveAborted = (n)->{
				statistics.collisionCount++;
				Physical.this.statistics.collisionCount++;
				log("collision detected");
				this.actionReceiveAborted.accept(this);
			};
			adapter.actionReceiveCompleted = (data)->{
				if(data==null) {
					statistics.receivedNoiseCount++;
					Physical.this.statistics.receivedNoiseCount++;
					log("received noise");
				}
				else {
					statistics.receivedCount++;
					Physical.this.statistics.receivedCount++;
					log("received data" + Arrays.toString(data));
					this.actionReceiveCompleted.accept(unpacker.apply(data).getSDU());
				}
			};
			adapter.actionReceiveStarted = (n)->{
				this.actionReceiveStarted.accept(this);
			};
			adapter.actionSendCompleted = (data)->{
				this.statistics.sentCount++;
				Physical.this.statistics.sentCount++;
				log("data sent" + Arrays.toString(data));
				this.actionSendCompleted.accept(data);
			};
		}
		
		/**
		 * Converts sdu byte array to PDU object
		 */
		public Function<byte[],PDU> packer = (sdu)->PDU.pack(Physical.framingBitOverhead, sdu);
		
		/**
		 * Converts pdu byte array to a PDU object
		 */
		public Function<byte[], PDU> unpacker = (data)->PDU.unpack(data,Physical.framingBitOverhead.length);

		/**Start transmitting a physical frame immediately
		 * @param data link layer frame
		 */
		public boolean send(byte[] data) {
			statistics.sentCount++;
			return adapter.send(packer.apply(data).data, BIT_RATE);
		}
		

		/**
		 * Called when the first bit of a frame is received
		 */
		public Consumer<NodeLayer> actionReceiveStarted = (n)->{}; 

		/**
		 * Called when a collision occurs during reception
		 */
		public Consumer<NodeLayer> actionReceiveAborted = (n)->{}; 
		
		/**Called when the last bit of a frame is received.
		 * @param data the link layer frame, if successfully received, 
		 * or <b>null</b> if a collision occurred.
		 */
		public Consumer<byte[]> actionReceiveCompleted = (pdu)->{};
		
		/**
		 * Called when the transmission is complete 
		 */
		public Consumer<byte[]> actionSendCompleted=(data)->{};

//		public double distanceTo(NodeLayer receiver) {
//			return adapter.distanceTo(receiver.adapter); 
//		}
//		
		/**
		 * @return node identified (serial number, unique over the same channel)
		 */
		public int getId() {
			return adapter.getId();
		}
		
		public boolean isIdle() {
			return adapter.isIdle(); 
		}
		
		public double whenSendingIsDone() {
			return adapter.getSendingDoneTime();
		}
		
		
		public void log(String str) {
			if(debug)
				System.out.println(String.format("t=%8.4f, Node#%3d %s\t"+str,
					now(), getId(), this.getClass().getName()));
		}

		public double getPos(int i) {
			return adapter.getPosition(i);
		}

		public Adapter getAdapter() {
			return this.adapter;
		}
	}

	boolean debug=false;
	public void enableDebug() {
		debug=true;
	}

	public void printStatistics() {
		System.out.println("Channel transmitted " + statistics.sentCount + " frames.");
		
	}
	
	public Physical(Channel channel) {
		this.channel = channel;
	}
	
	public double now() {
		return Scheduler.it.now();
	}
	
	static class Test{
		static class MyPhy extends Physical {

			public MyPhy(Channel channel) {
				super(channel);
			}

			class MyNode extends NodeLayer {
				
				public MyNode(double x, double y, double z) {
					super(new Node(new double[] {x, y, z}));
					packer=(sdu)->PDU.pack(new byte[] {(byte) getId()}, sdu);
					unpacker=(data)->PDU.unpack(data,1);
					actionReceiveStarted = (n)->n.log("receive started");
					actionSendCompleted = (data)->log("send completed, " + Arrays.toString(data));
					actionReceiveAborted = (n)->n.log("collision detected");
					actionReceiveCompleted = (pdu)-> {
						if(pdu.length==0)
							log("received noise");
						else
							log("received data " + Arrays.toString(pdu));
					};
				}
			}
		}
		
		public static void main(String[] args) {
			Channel channel = new Channel();
			MyPhy phy = new MyPhy(channel);
			MyPhy.MyNode[] nodes = {
					phy.new MyNode(0,0,0),	
					phy.new MyNode(0,100,0),
					phy.new MyNode(0,200,0),
					phy.new MyNode(0,300,0),
					phy.new MyNode(0,400,0),
					phy.new MyNode(0,500,0),
					phy.new MyNode(0,600,0),
					phy.new MyNode(0,700,0)
			};
			for(MyPhy.MyNode n:nodes) {
				n.getAdapter().enableDebug();
			}

			phy.enableDebug();
			phy.channel.enableDebug();
			
			nodes[0].send(new byte[] {1,2,3,4,5,6,7,8,9,10});
			Scheduler.it.schedule(0.1, 
					(s)->nodes[7].send(new byte[] {10,9,8,7,6,5,4,3,2,1}), ()->"2nd send");
			Scheduler.it.schedule(0.279, 
					(s)->nodes[3].send(new byte[] {1,2,3,4,5,5,4,3,2,1}), ()->"3rd send");
	
			Scheduler.it.handleAll();
		}
	}
}
