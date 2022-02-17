package uansim;

import java.util.Random;
import java.util.function.Consumer;

import uansim.Node;

public class DPR extends DBR {
	public static double[] FORWARDING_PROBABILITY=new double[256];
	public double getForwardingProbability(int i) {
		return FORWARDING_PROBABILITY[i];
	}
	public void setForwardingProbability(int i, double p) {
		FORWARDING_PROBABILITY[i]=p;
	}
	public DPR(Physical physical) {
		super(physical);
		rand = new Random();
	}

	public DPR(Physical physical, Random decisionRand) {
		super(physical);
		this.rand = decisionRand;
	}
	
	public Random rand;

	public class NodeLayer extends DBR.NodeLayer {
		
		public class DelayedForwardingQueue2 extends DelayedForwardingQueue{
			
			class Record extends DelayedForwardingQueue.Record {
				boolean markedForDropping=false;

				public Record(Packet packet, double forwardingTime) {
					super(packet, forwardingTime);
				}
				
			}

			@Override
			DelayedForwardingQueue.Record
				makeRecord(Packet packet, double forwardingTime) {
				return new Record(packet, forwardingTime);
			}
			
			public void markForDropping(Packet packet) {
				Record record = (Record) map.get(packet);
				record.markedForDropping=true;
			}
			
			public boolean isMarkedForDropping(Packet packet) {
				Record record = (Record) map.get(packet);
				if(record==null) {
					log("packet not found in DF queue!!!");
					return false;
				}
				return record.markedForDropping;

			}
			
		}
		
		public NodeLayer(Node node) {
			super(physical.new NodeLayer(node));
			initNodeLayer();
		}

		public NodeLayer(Physical.NodeLayer node) {
			super(node);
			initNodeLayer();
		}

		private void initNodeLayer() {
			delayedForwardingQueue = new DelayedForwardingQueue2();

			cancelDelayedDecision = (packet) -> {
				((DelayedForwardingQueue2)this.delayedForwardingQueue).markForDropping(packet);
				log("packet marked for dropping");
			};
			
			Consumer<Packet> superDelayedForwardDecision = delayedForwardDecision;
			delayedForwardDecision = (packet) -> {
				if(packet.getSource()<0 ||packet.getSource()>=FORWARDING_PROBABILITY.length)
					System.err.println("src="+packet.getSource()+", FP.length="+FORWARDING_PROBABILITY.length);
				if(!((DelayedForwardingQueue2)this.delayedForwardingQueue)
						.isMarkedForDropping(packet)
						|| rand.nextDouble()<FORWARDING_PROBABILITY[packet.getSource()]) 
					superDelayedForwardDecision.accept(packet);
				else {
					cancelDelayedDecision.accept(packet);	
				}
			};			
		}
	}
	
	static class Test{
		public static void main(String[] args) {
			Channel channel = new Channel();
			Physical physical = new Physical(channel);
			DPR network = new DPR(physical);

			DBR.MAXIMUM_FORWARDING_DELAY=100;
			DPR.FORWARDING_PROBABILITY[0]=0.0;

			NodeLayer source = network.new NodeLayer(physical.new NodeLayer(new Node(new double[]{200,200,100})));
			AttackNodeLayer attacker  = network.new AttackNodeLayer(physical.new NodeLayer(new Node(new double[]{200,200,101})));
			NodeLayer relay1 = network.new NodeLayer(physical.new NodeLayer(new Node(new double[]{200,100,50})));
//			NodeLayer relay2 = network.new Node(physical.new Node(new NodeLayer(new double[]{100,220,40})));
			NodeLayer relay3 = network.new NodeLayer(physical.new NodeLayer(new Node(new double[]{100,100,30})));
			NodeLayer sink = network.new NodeLayer(physical.new NodeLayer(new Node(new double[]{100,100,0})));
			
			Application app = new Application(0, source, sink, 0.25, 12345)
					.setNumberOfMessages(100)
					.setGenerationBitRate(10)
					.setMessageBitLength(100); 
			
			Scheduler.it.enableDebug();
			network.enableDebug();

			Scheduler.it.handleAll();
			physical.printStatistics();
			app.printStatistics();
						
		}
	}
}
