package uansim;

import java.util.function.Consumer;

import uansim.Channel.PowerDelayProfile;
import utils.Geometry;

/**A network interface 
 * @author saleh
 *
 */
public class Adapter {
	public static double TRANSMISSION_POWER = 70; // dBm
	/**Receiver sensitivity, i.e., minimum detectable SNR
	 * 
	 */
	public static double RECEIVER_SENSITIVITY=25;  //dB

	Node node;

	public Adapter(Node node) {
		this.node = node;
		node.installAdapter(this);
	}

	public static class State {
		public Adapter receivingFrom=null;
		public double signalPower=0.0;
		public double interferencePower=0;
		
		public boolean collision=false;
		public boolean busySending=false;
		public double sendingFinishTime=-1;
	}
	
	protected State state = new State();
	
	public State getState() {
		return state;
	}
	
	public static class Statistics {
		public int sentCount=0;
		public int receivedCount=0;
		public int collisionCount=0;
		
		public void reset() {
			sentCount=0;
			receivedCount=0;
		}
	}
	
	/**
	 * Adapter monitoring statistics
	 */
	protected Statistics statistics = new Statistics();

	/**
	 * Called when the first bit of a frame is arrives
	 */
	public Consumer<Adapter> actionReceiveStarted = (v)->{}; 

	/**
	 * Called when a collision occurs during reception
	 */
	public Consumer<Adapter> actionReceiveAborted = (v)->{}; 
	
	/**Called when the last bit of a frame is received.
	 * @param data the physical layer frame, if successfully received, 
	 * or <b>null</b> if a collision occurred.
	 */
	public Consumer<byte[]> actionReceiveCompleted=(data)->{};
	
	/**Called when the last bit of a frame is sent out
	 * @param data the physical layer frame
	 */
	public Consumer<byte[]> actionSendCompleted=(data)->{};
	
	/**
	 * Called when the first bit of a link layer frame arrives at the node 
	 * 
	 * @param sender adapter sending the received signal
	 * @param power of the received signal power in dBm
	 */
	protected void beginReceive(Adapter sender, double power) {
		// if not already sending nor receiving
		if(! state.busySending && (state.receivingFrom==null)) {
			//TODO use the modulation frequency band
			double snr=power-state.interferencePower-channel.getAmbientNoisePower(new double[]{0,0});
			if(snr>RECEIVER_SENSITIVITY) {
				state.signalPower=power;
				state.receivingFrom=sender;
				actionReceiveStarted.accept(this);				
			}
			else {
				state.interferencePower+=power;
			}
		}
		//if already sending or receiving
		else {
			//consider the received as noise
			state.interferencePower+=power;
			
			// if already receiving and no collision has been detected yet
			if(state.receivingFrom!=null && !state.collision) {
				assert(state.receivingFrom!=sender);
				
				//TODO use the modulation frequency band
				double snr=state.signalPower-state.interferencePower-channel.getAmbientNoisePower(new double[]{0,0});
				if(snr<RECEIVER_SENSITIVITY) {
					state.collision=true;
					actionReceiveAborted.accept(this);
				}
			}
		}
	}
	
	/**Called when the last bit of a link layer frame arrives at the node
	 * @param power received signal power in dBm
	 * @param data the link layer frame, if no collision occurred, or 
	 * <b>null</b> if a collision occurred
	 */
	protected void endReceive(Adapter sender, double power, byte[] data) {
		if(state.receivingFrom==sender) {
			if(!state.collision){
				statistics.receivedCount++;
				actionReceiveCompleted.accept(data);
			}
			else {
				this.statistics.collisionCount++;				
			}
			state.receivingFrom=null;
			state.collision=false;
		}
		else {
			state.interferencePower-=power;
			if(Math.abs(state.interferencePower)<1e-6)
				state.interferencePower=0;
		}
	}
	
	/**Sending time is over for the given data frame
	 * @param data the data that was being sent <b>null</b> if sending failed.
	 * Update state and signal completion to observers. 
	 */
	protected void endSend(byte[] data) {
		state.busySending=false;
		state.collision = false;
		actionSendCompleted.accept(data);
	}
	
	protected double[] getPosition() {
		return node.getPosition();
	}

	/**
	 * @return node identified (serial number, unique over the same channel)
	 */
	public int getId() {
		return node.getId();
	}
	
	public boolean isIdle() {
		return !this.state.busySending && this.state.receivingFrom==null;
	}
	
	public double getSendingDoneTime() {
		if(!state.busySending)
			return Scheduler.it.now();
		else
			return state.sendingFinishTime;
	}
	
	boolean debug=false;
	
	public void enableDebug() {
		debug=true;
	}
	
	protected void log(String str) {
		if(debug) {
			System.out.println(String.format("t=%8.4f, Node#%3d %s\t"+str,
				Scheduler.it.now(), getId(), this.getClass().getName()));				
			LogFile.logFile.println(String.format("t=%8.4f, Node#%3d %s\t"+str,
					Scheduler.it.now(), getId(), this.getClass().getName()));				

		}
	}

	public double getPosition(int i) {
		return getPosition()[i];
	}

	
	Channel channel=null;
	
	public void connect(Channel channel) {
		this.channel=channel;
	}
	
	public boolean isConnected() {
		return channel!=null;
	}

	/**Start transmitting a link layer frame immediately
	 * @param data link layer frame
	 * @param bitRate the transmission rate
	 */
	public boolean send(byte[] data, double bitRate) {
		if(!isConnected())
			return false;
		
		if(state.busySending) {
			return false;
		}
		if(state.receivingFrom!=null) {
			state.collision=true;
			actionReceiveAborted.accept(this);
		}
		//Keep statistics
		channel.statistics.sentCount++;
		statistics.sentCount++;
		
		state.busySending=true;
		log("BeginSend");
		double transmissionTime = (data.length)*8/bitRate;
		state.sendingFinishTime = transmissionTime+Scheduler.it.now();
		Scheduler.it.schedule(transmissionTime, 
				(s)->endSend(data), 
				()->"EndSend("+getId()+")");
		
		for(Adapter receiver: channel.getAdapters()) {
			if(this!=receiver) {
//				double distance=this.node.distanceTo(receiver.node);
//				if(distance<=interferenceRange) {
					//TODO use modulation bandwidth
					PowerDelayProfile pdp = channel.getPowerDelayProfile(getPosition(), receiver.getPosition(), new double[]{0,0});
					double delay=pdp.getDelay();
					double power=TRANSMISSION_POWER;		
					double gain=pdp.getPower();
					double receivePower=power+gain;
					Scheduler.it.schedule(delay, 
							(s)->receiver.beginReceive(this,receivePower),
							()->"BeginReceive("+this.getId()+"->"+receiver.getId()+")");

					Scheduler.it.schedule(delay+transmissionTime, 
							(s)->receiver.endReceive(this,receivePower,data), 
							()->"EndReceive("+this.getId()+"->"+receiver.getId()+")::data");
					
//					if(distance<commRange)
//						Scheduler.it.schedule(delay+transmissionTime, 
//								(s)->receiver.endReceive(this,receivePower,data), 
//								()->"EndReceive("+this.getId()+"->"+receiver.getId()+")::data");
//					else
//						Scheduler.it.schedule(delay+transmissionTime, 
//								(s)->receiver.endReceive(this,receivePower,null), 
//								()->"EndReceive("+this.getId()+"->"+receiver.getId()+")::noise");						
				}
			}
//		}
		return true;
	}

	public void send(byte[] data) {
		this.send(data,1000);
	}

//	public double distanceTo(Adapter adapter) {
//		if(this.channel==adapter.channel)
//			return this.node.distanceTo(adapter.node);
//		else
//			return Double.POSITIVE_INFINITY;
//	}
//
}
