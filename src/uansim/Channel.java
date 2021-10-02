package uansim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Supplier;
import java.util.function.ToDoubleBiFunction;

import uansim.Channel.PowerDelayProfile;
import uansim.Network.Layer;
import utils.Geometry;

/**
 * @author Saleh
 *An under-water channel carries acoustic signals from sources to destinations adding ambient noise
 */
public class Channel {
	
	/**
	 * Power-Delay profile is a representation of the received signal power and delay over multiple paths.<br>
	 * A <i>tap</i> represents a single copy of the signal over one path, which specifies the 
	 * delay of the path and its loss (always negative).
	 */
	public static class PowerDelayProfile {
		double[][] taps;
		
		/**Create a PowerDelayProfile with the specified number of taps. The taps should then be set individually using {@link #setTap(int, double, double) setTap}
		 * @param numTaps number of taps. For single path propagation the number of taps = 1.
		 */
		public PowerDelayProfile(int numTaps) {
			taps=new double[numTaps][2];
		}
		
		/**Get the number of taps in the profile
		 * @return the number of taps.
		 */
		public int getNumTaps() {
			return taps.length;
		}
		
		/**Set the specified tap to the given power gain and delay
		 * @param i tap index
		 * @param power path loss (should be negative)
		 * @param delay path delay
		 */
		public PowerDelayProfile setTap(int i, double power, double delay) {
			taps[i][0]=power;
			taps[i][1]=delay;
			return this;
		}
		
		/**Get the path loss of the specified tap
		 * @param i tap index
		 * @return path loss (should be negative)
		 */
		public double getTapPower(int i) {
			return taps[i][0];
		}
		
		/**Get the path delay of the specified tap
		 * @param i path index
		 * @return the delay of the specified path
		 */
		public double getTapDelay(int i) {
			return taps[i][1];
		}
		
		/**Return the delay of the main path, i.e. the path with the lowest loss 
		 * (usually the tap with the lowest delay).
		 * @return the delay of the main path
		 */
		public double getDelay() {
			int j=0;
			double p=getTapPower(j);
			for(int i=1; i<taps.length;i++) {
				if(getTapPower(i)>p) {
					j=i;
					p=getTapPower(j);
				}
			}
			return getTapDelay(j);
		}
		
		/**Get the loss of the main path, i.e. the path with the lowest loss.
		 * @return the loss of the main path (should be negative)
		 */
		public double getPower() {
			int j=0;
			double p=getTapPower(j);
			for(int i=1; i<taps.length;i++) {
				if(getTapPower(i)>p) {
					j=i;
					p=getTapPower(j);
				}
			}
			return p;
		}
	}

	/**
	 * The propagation model which determines the received signal 
	 */
	public interface PropagationModel {

		/**Get the Power-Delay profile corresponding to the multipath between the given source and destination
		 * @param srcPos the coordinates of the source
		 * @param dstPos the coordinates of the destination
		 * @param f frequency band, i.e. f[0] low frequency, f[1] high frequency
		 * @return the power-delay profile for the transmission from the given source to the given destination.
		 */
		PowerDelayProfile getPowerDelayProfile(double[] srcPos, double[] dstPos, double f[]);

	}
	
	/**
	 * The ideal propagation model assumes single path propagation with no path loss.
	 * 
	 */
	public static class IdealPropagationModel implements PropagationModel {
		/**
		 * Propagation speed in (m/s)
		 */
		double propagationSpeed=1500;	//(m/s)

		/**Get the Power-Delay profile corresponding to the multipath between the given source and destination
		 * @param srcPos the coordinates of the source
		 * @param dstPos the coordinates of the destination
		 * @param f frequency band, i.e. f[0] low frequency, f[1] high frequency
		 * @return the power-delay profile for the transmission from the given source to the given destination.
		 */
		@Override
		public PowerDelayProfile getPowerDelayProfile(double[] srcPos, double[] dstPos, double[] f) {
			double d = Geometry.distance(srcPos,dstPos);
			return new PowerDelayProfile(1).setTap(0, 0, d/propagationSpeed);
		}
	}

	/**
	 * Thorp's single path propagation model. The path loss has two components<ol>
	 * <li>Absorptive attenuation
	 * <li>Spreading</ol>
	 * Therefore, two parameters determine the path loss according to Thorp's model<ol>
	 * <li> frequency-dependent absorption coefficient (see {@link ThorpPropagationModel#getAttenuationRate(double) getAttenuationRate}),
	 * <li> spreading coefficient, which depends on the geometry of the antenna and the depth of the sea floor.</ol>
	 */
	public static class ThorpPropagationModel implements PropagationModel {
		/**
		 * Propagation speed in (m/s)
		 */
		double propagationSpeed;
		
		/**
		 * spreading coefficient in dB
		 */
		double spreadingCoefficient;
		
		/**Initialize Thorp's propagation model with the given spreading coefficient
		 * @param propagationSpeed the constant propagation speed of acoustic waves in water
		 * @param spreadingCoefficient the rate of signal power loss due to geometric spreading
		 */
		public ThorpPropagationModel(double propagationSpeed, double spreadingCoefficient) {
			this.spreadingCoefficient=spreadingCoefficient;
			this.propagationSpeed=propagationSpeed;
		}
		
		/**
		 * Get the attenuation rate in dB per meter
		 */
		double getAttenuationRate(double f0) {
			f0 *= 1e-3;
			double f2 = f0 * f0;
			 
			if (f0 >= 0.4)
				return 1e-3*(0.11 * f2 / (1 + f2) + 44 * f2 / (4100 + f2) + 2.75e-4 * f2 + 3e-3);
			else
				return 1e-3*(0.002 + 0.11 * (f0 / (1 + f0)) + 0.011 * f0);
		}
		
		/**Get the Power-Delay profile corresponding to the multipath between the given source and destination
		 * @param srcPos the coordinates of the source
		 * @param dstPos the coordinates of the destination
		 * @param f frequency band, i.e. f[0] low frequency, f[1] high frequency
		 * @return the power-delay profile for the transmission from the given source to the given destination.
		 */
		@Override
		public PowerDelayProfile getPowerDelayProfile(double[] srcPos, double[] dstPos, double f[]) {
			double d = Geometry.distance(srcPos,dstPos);
			double p = 10*spreadingCoefficient*Math.log10(d)+getAttenuationRate((f[0]+f[1])/2)*d;
			
			PowerDelayProfile result=new PowerDelayProfile(1);
			result.setTap(0, -p, d/propagationSpeed);
			return result;
		}				
		
	}

	
	/**Represents one source of ambient noise which can be queried to 
	 * calculate the ambient noise power within a given frequency band.
	 * 
	 */
	public static abstract class AmbientNoiseSource{
		
		/**Calculate the ambient noise power within the given frequency band
		 * @param f frequency band, i.e. f[0] low frequency, f[1] high frequency
		 * @return ambient noise power 
		 */
		public double getPower(double f[]) {
			double mean=0;
			double df=(f[1]-f[0])/100 ;
			for(double fi=f[0]+df/2; fi<f[1]; fi+=df) {
				mean+=0.01*getPowerDensity(fi);
			}
			return (f[1]-f[0])*mean;
		}
		/**Calculate the ambient noise power density at the given frequency
		 * @param f frequency
		 * @return ambient power density
		 */
		public abstract double getPowerDensity(double f);
	}
	
	public static class MultisourceAmbientNoise extends AmbientNoiseSource{
		ArrayList<AmbientNoiseSource> sources=new ArrayList<AmbientNoiseSource>();		
		
		public MultisourceAmbientNoise addSource(AmbientNoiseSource src) {
			sources.add(src);
			return this;
		}

		@Override
		public double getPowerDensity(double f) {
			double maxP=0;
			for(AmbientNoiseSource ans: sources) {
				maxP=Math.max(maxP, ans.getPowerDensity(f));
			}
			return maxP;
		}
	}
	
	public static class TurbulenceNoise extends AmbientNoiseSource{

		@Override
		public double getPowerDensity(double f) {
			return 17.0 - 30*Math.log10(f);
		}
		
	}

	public static class ShippingNoise extends AmbientNoiseSource{
		/**
		 * Shipping activity factor: 0=low, 1=high
		 */
		final double s;
		/**
		 * @param s Shipping activity factor: 0=low, 1=high
		 */
		public ShippingNoise(double s) {
			this.s=s;
		}
		@Override
		public double getPowerDensity(double f) {
			return 40.0 + 20*(s-0.5)+26*Math.log10(f)-60*Math.log10(f+0.03);
		}
	}
	
	public static class WaveNoise extends AmbientNoiseSource{
		/**
		 * Wind speed
		 */
		final double w;
		/**
		 * @param w the wind speed (m/s)
		 */
		public WaveNoise(double w) {
			this.w=w;
		}
		
		@Override
		public double getPowerDensity(double f) {
			return 50.0 + 7.5*Math.sqrt(w)+20*Math.log10(f)-40*Math.log10(f+0.4);
		}
	}
	
	public static class ThermalNoise extends AmbientNoiseSource{
		@Override
		public double getPowerDensity(double f) {
			return -15+20*Math.log10(f);
		}
	}
	
	PropagationModel propagationModel = new IdealPropagationModel();
	
	/**Calculates the received signal power-delay profile given the source position, destination position and frequency band. 
	 * By default the {@link IdealPropagationModel} is used for the calculation. 
	 * To change the default behavior use {@link #setPropagationModel(PropagationModel) setPropagationModel}
	 * @param srcPos coordinates of the transmitted
	 * @param dstPos coordinates of the receiver
	 * @param f frequency band, i.e. f[0] low frequency, f[1] high frequency
	 * @return power-delay profile of the received signal
	 */
	public PowerDelayProfile getPowerDelayProfile(double[] srcPos, double[] dstPos, double[] f) {
		return propagationModel.getPowerDelayProfile(srcPos, dstPos, f);
	}

	AmbientNoiseSource ambientNoiseSource=null;
	
	/**Calculates the ambient noise power in the specified frequency band. 
	 * By default the ambient noise is ignored, and this function returns 0.
	 * To change the default behavior use {@link #setAmbientNoiseSource(AmbientNoiseSource) setAmbientNoiseSource}.
	 * @param f frequency band, i.e. f[0] low frequency, f[1] high frequency
	 * @return ambient noise power
	 */
	public double getAmbientNoisePower(double[] f) {
		if(ambientNoiseSource!=null)
			return ambientNoiseSource.getPower(f); //(dBm)
		else
			return 0;
	}
		
	/**Specify the propagation model object used for calculating the power-delay profile of received signals
	 * @param propagationModel PropagationModel object 
	 * @return this Channel object for chaining
	 */
	public Channel setPropagationModel(PropagationModel propagationModel) {
		this.propagationModel=propagationModel;
		return this;
	}
	
	/**Specify the ambient noise power calculation function. 
	 * @param ambientNoiseSource an object which implements the getAmbientNoisePower function. 
	 * See {@link MultisourceAmbientNoise}, {@link TurbulenceNoise}, {@link ShippingNoise}, {@link WaveNoise}, and {@link ThermalNoise}
	 * @return this Channel object for chaining
	 */
	public Channel setAmbientNoiseSource(AmbientNoiseSource ambientNoiseSource) {
		this.ambientNoiseSource=ambientNoiseSource;
		return this;
	}
	
	public static class Statistics {
		public int sentCount=0;
		public int collisionCount=0;
	}
	
	public Statistics statistics = new Statistics();
	
	public static class NodeStatistics {
		public int sentCount=0;
		public int receivedCount=0;
		public int collisionCount=0;
		
		public void reset() {
			sentCount=0;
			receivedCount=0;
		}
	}

	public static class NodeState {
		public int busyReceiving=0;
		public boolean collision=false;
		public boolean busySending=false;
		public double sendingFinishTime=-1;
	}
	
	public void addNodes(Node[] nodes) {
		for(Node node: nodes) {
			for(Adapter adapter: node.adapters) {
				if(adapter instanceof uansim.Adapter)
					adapters.add((uansim.Adapter) adapter);
			}
		}
	}
	boolean debug=false;

	public void enableDebug() {
		debug=true;
	}

	public void printStatistics() {
		System.out.println("Channel transmitted " + statistics.sentCount + " frames.");
		
	}


	/**
	 * List of all nodes connected to the channel
	 */
	private ArrayList<uansim.Adapter> adapters;
	
	/** Connect an adapter to the channel
	 * @param adapter the adapter to be connected
	 */
	public void connect(uansim.Adapter adapter) {
		if(!adapters.contains(adapter)) {
			adapters.add(adapter);
			adapter.connect(this);
		}
	}

	/**Get list of adapters
	 * @return list of all nodes connected to the channel
	 */
	public ArrayList<uansim.Adapter> getAdapters() {
		return adapters;
	}
	
	/**Get the (first) connected adapter with the given nodeId
	 * @param nodeId the node for which the adapter is requested
	 * @return the first connected adapter for with the given nodeId or <b>null</b> if no such node exists.
	 */
	public Adapter getAdapterById(int adapterId) {
		for(Adapter adapter: adapters)
			if(adapter.getId()==adapterId)
				return adapter;
		return null;
	}
	
	/**
	 * Construct a new Channel object. 
	 */
	public Channel() {
		this.adapters= new ArrayList<uansim.Adapter>();
	}
	
	static class Test{

		public static void main(String[] args) {
			
			Node[] nodes = {
				new Node(new double[] {0,0,0}),	
				new Node(new double[] {0,100,0}),
				new Node(new double[] {0,200,0}),
				new Node(new double[] {0,300,0}),
				new Node(new double[] {0,400,0}),
				new Node(new double[] {0,500,0}),
				new Node(new double[] {0,600,0}),
				new Node(new double[] {0,700,0})
			};
			
			Channel channel = new Channel()
					.setAmbientNoiseSource(
							new MultisourceAmbientNoise()
							.addSource(new TurbulenceNoise())
							.addSource(new WaveNoise(0))	//wind speed = 0
							.addSource(new ShippingNoise(0)) //low shipping activity
							.addSource(new ThermalNoise()))
					.setPropagationModel(
							new ThorpPropagationModel(1500, 2.0));
			channel.enableDebug();
			
			for(Node node: nodes) {
				uansim.Adapter adapter = new uansim.Adapter(node);
				adapter.enableDebug();
				channel.connect(adapter);

				adapter.actionReceiveAborted = (n)->n.log("collision detected");
				adapter.actionReceiveCompleted = (data)-> {
					if(data==null) {
						adapter.log("received noise");			
					}
					else {
						adapter.log("received data " + Arrays.toString(data));
					}
				};
			}
			
			nodes[0].getDefaultAdapter().send(new byte[] {1,2,3,4,5,6,7,8,9,10},1000);
			Scheduler.it.schedule(0.1, 
					(s)->nodes[7].getDefaultAdapter().send(new byte[] {10,9,8,7,6,5,4,3,2,1},1000), ()->"2nd send");
			Scheduler.it.schedule(0.279, 
					(s)->nodes[3].getDefaultAdapter().send(new byte[] {1,2,3,4,5,5,4,3,2,1},1000), ()->"3rd send");
			
			Scheduler.it.schedule(0.34, 
					(s)->System.out.println(
							"Node #4: isReceiving=" + (nodes[4].getDefaultAdapter().state.receivingFrom!=null) + 
							", signalPower=" +nodes[4].getDefaultAdapter().state.signalPower +
							", collision="+nodes[4].getDefaultAdapter().state.collision
							), ()->"checking node #4 state");

			Scheduler.it.handleAll();
			
			for(Adapter a: channel.adapters) {
				System.out.println("Adapted #"+a.getId()+" interference=" +a.state.interferencePower);
			}
		}
	}

}
