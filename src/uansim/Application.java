package uansim;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.util.Random;
import java.util.function.Consumer;

import uansim.Network.Layer;
import uansim.Scheduler.Event;

public class Application {

	public Application(double startDelay, Layer source, Layer sink, double randFactor, int trafficSeed) {
		sourceNodeId = (byte) source.getId();
		this.sendFunction=(message)->source.sendTo(message, sink);
		sink.deliver = sink.deliver.andThen((message)->deliver((byte[])message));
		
		this.startDelay=startDelay;
		intervalSeconds = messageBitLength/generationBitRate;
		setStartDelay(startDelay);
		rand = new Random(trafficSeed);
		this.randFactor = randFactor;
	}

	private double randFactor;
	private Random rand;
	private Event startEvent;
	//Connectors
	private Consumer<byte[]> sendFunction;

	byte sourceNodeId;
	
	//Data generation parameters
	double startDelay=0;	// (s)
	int messageBitLength=800;	// (Bits)
	double generationBitRate=100;	// (Bps)
	int numberOfMessages=100;	// (message)
	double intervalSeconds;
	
	//Statistics
	public int generatedCount=0;
	public int receivedCount=0;
	private double endToEndDelay=0;

	public double getDeliveryRatio() {
		return (100.0* receivedCount)/generatedCount;
	}
	public double getAverageEndToEndDelay() {
		
		return (receivedCount==0?-1:endToEndDelay/receivedCount);
	}
	public void printStatistics() {
		System.out.println("Application sent "+generatedCount+
				" messages and received "+ receivedCount + " message.");
		System.out.println("Delivery ratio = " + 
				String.format("%4.1f%%", getDeliveryRatio()) );
		System.out.println("Average end-to-end delay = " + 
				String.format("%8.4f s", getAverageEndToEndDelay()) );
	}
	
	private void deliver(byte[] message) {
		if(message[1]==sourceNodeId) {
			// Extract timestamp from application message to calculate end-to-end delay
			if(messageBitLength>=80) {
				ByteBuffer buffer = ByteBuffer.wrap(message,2,8);
				DoubleBuffer doubleBuffer = buffer.asDoubleBuffer();
				double sendingTime = doubleBuffer.get();
				endToEndDelay+=Scheduler.it.now()-sendingTime;
			}
			receivedCount++;
		}
	}
	
	byte messageId;
	
	private void run() {
		messageId=0;
		for(int i=0; i<numberOfMessages; i++) {
			Scheduler.it.schedule(
					i*intervalSeconds  + (randFactor)*intervalSeconds*rand.nextDouble(),
//					i*(1-randFactor)*intervalSeconds+(2*randFactor) *intervalSeconds*rand.nextDouble(), 
					(s)->send(), 
					()->"Send application message");
		}
	}
	
	private void send() {
		byte[] message = new byte[(messageBitLength+7)/8];
		message[0] = (byte) messageId++;
		message[1] = sourceNodeId;
		
		// Embed timestamp in application message to calculate end-to-end delay at the sink
		if(messageBitLength>=80) {
			ByteBuffer buffer = ByteBuffer.wrap(message,2,8);
			DoubleBuffer doubleBuffer = buffer.asDoubleBuffer();
			doubleBuffer.put(Scheduler.it.now());
		}
		sendFunction.accept(message);
		generatedCount++;
	}

	public double getStartDelay() {
		return startDelay;
	}

	public Application setStartDelay(double startDelay) {
		this.startDelay = startDelay;
		if(startEvent!=null) 
			this.startEvent.cancel();
		startEvent = Scheduler.it.schedule(startDelay, (s)->run(), ()->"Application started.");

		return this;
	}

	public int getMessageLength() {
		return messageBitLength;
	}

	public Application setMessageBitLength(int messageBitLength) {
		if(messageBitLength<80)
			System.err.println("Can't calculate end-to-end delays for messges shorter than 10 bytes long");
		this.messageBitLength = messageBitLength;
		intervalSeconds = messageBitLength/generationBitRate;
		return this;
	}

	public double getGenerationBitRate() {
		return generationBitRate;
	}

	public Application setGenerationBitRate(double generationBitRate) {
		this.generationBitRate = generationBitRate;
		intervalSeconds = messageBitLength/generationBitRate;
		return this;
	}

	public int getNumberOfMessages() {
		return numberOfMessages;
	}

	public Application setNumberOfMessages(int numberOfMessages) {
		this.numberOfMessages = numberOfMessages;
		return this;
	}
}
