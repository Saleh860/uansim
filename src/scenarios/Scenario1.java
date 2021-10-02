package scenarios;

import java.util.Arrays;

import uansim.DBR;
import uansim.DPR;
import uansim.Physical;
import uansim.DPR.NodeLayer;
import uansim.Application;
import uansim.Channel;
import uansim.Scheduler;
import uansim.Node;

public class Scenario1 {

	public static void main(String[] args ) {
		Channel channel=new Channel();
		Physical physical = new Physical(channel);
		DPR network = new DPR(physical);
		DBR.MAXIMUM_FORWARDING_DELAY=100;
		Arrays.fill(DPR.FORWARDING_PROBABILITY,0.0);
		NodeLayer source = network.new NodeLayer(new Node(new double[]{200,200,100}));
//		Channel.AttackNode attacker  = channel.new AttackNode(new double[]{200,200,101});
		NodeLayer relay1 = network.new NodeLayer(new Node(new double[]{200,100,50}));
//		Channel.Node relay2 = channel.new Node(new double[]{100,220,40});
		NodeLayer relay3 = network.new NodeLayer(new Node(new double[]{100,100,30}));
		NodeLayer sink = network.new NodeLayer(new Node(new double[]{100,100,0}));
		
		Application app = new Application(0, source, sink, 0, 12345)
				.setNumberOfMessages(100)
				.setGenerationBitRate(10)
				.setMessageBitLength(100); 
		
//		channel.getScheduler().setDebug();
//		channel.enableDebug();

		Scheduler.it.handleAll();
		physical.printStatistics();
		app.printStatistics();
		
	}
}
