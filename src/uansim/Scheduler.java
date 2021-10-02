package uansim;
import java.io.PrintStream;
import java.util.PriorityQueue;
import java.util.function.*;

/**
 * An event priority queue
 * @author saleh
 *
 */
public class Scheduler {
	/**
	 * The one and only Scheduler object
	 */
	public static Scheduler it = new Scheduler();
	
	/**
	 * Make sure there can only be one Scheduler object
	 */
	private Scheduler() {
		
	}
	
	/**
	 * current time
	 */
	private double time = 0;
	
	/**
	 * flag indicating that debug is enabled 
	 */
	private boolean debug = false;
	
	/**
	 * event queue, events are sorted in ascending order of time 
	 */
	private PriorityQueue<Event> queue = new PriorityQueue<Event>();
	
	
	/**Cancels the given event, by removing it from the event queue
	 * @param event to be canceled
	 * @return <b>true</b> if event was removed from the queue 
	 * and <b>false</b> if event wasn't found in the queue.
	 */
	private boolean cancel(Event event) {
		return queue.remove(event);
	}
	
	/** 
	 * A scheduled action with a corresponding description
	 *
	 */
	public class Event implements Comparable<Event>{
		/**
		 * Event time
		 */
		private double time;
		/**
		 * Action function
		 */
		private Consumer<Scheduler> actionFn;
		/**
		 * Description function
		 */
		private Supplier<String> descriptionFn;
		
		/**
		 * @param time scheduled time
		 * @param actionFn action function
		 * @param descriptionFn description of event
		 */
		Event(double time, Consumer<Scheduler> actionFn, Supplier<String> descriptionFn) {
			this.time=time;
			this.actionFn=actionFn;
			this.descriptionFn=descriptionFn;
		}
		
		/**
		 * Perform the scheduled action
		 */
		private void action() {
			this.actionFn.accept(Scheduler.this);
		}
		
		/**
		 * order events by ascending time
		 */
		@Override
		public int compareTo(Event o) {
			return Double.compare(this.time, o.time);
		}
		
		/**
		 *Return a string containing the time and description of the event.
		 */
		@Override
		public String toString() {
			return "t=" + String.format("%8.4f",time) + ", " + descriptionFn.get();
		}
		
		/**
		 * Remove the event from the scheduler event queue
		 * @return <b>true</b> if event was removed from the queue 
		 * and <b>false</b> if event wasn't found in the queue, for instance, 
		 * if the event has already been removed earlier.
		 */
		public boolean cancel() {
			return it.cancel(this);
		}
	}
	
	
	/**
	 * Display debugging information during Scheduler operation
	 */
	public void enableDebug() {
		debug=true;
	}
	
	/**
	 * Get the current time
	 * @return the current scheduler time
	 */
	public double now() {
		return time;
	}
	
	/**Get the number of queued future events
	 * @return number of scheduled events waiting in the queue
	 */
	public int getQueueLength() {
		return this.queue.size();
	}

	/**Schedule an event at now()+delay that is handled by the given handler and has the given description
	 * @param delay event time relative to now(), negative values are treated as zero
	 * @param actionFn a function to handle the event
	 * @param descriptionFn function that returns the description. this is good when the description is dynamic
	 * @return the event that has just been scheduled
	 */
	public Event schedule(double delay, Consumer<Scheduler> actionFn, Supplier<String> descriptionFn) {
		Event event = new Event(time+delay, actionFn, descriptionFn);
		queue.add(event);
		return event;
	}

		
	/** Execute the event action at the top of the event queue
	 * @return true if an action was executed, false if the queue has been found empty
	 */
	public boolean handleNext() {
		if(queue.isEmpty())
			return false;
		Event e = queue.poll();
		time = e.time;
		if(debug) {
			System.out.println(e);
			LogFile.logFile.println(e);
		}
		e.action();
		return true;
	}
	
	/**
	 * Execute all event actions in the event queue until the queue becomes empty. 
	 * The execution of an event action can potentially enqueue further events which must also be handled before handleAll returns.
	 */
	public void handleAll() {
		while(handleNext());
	}
	
	/**Print the times and descriptions of all events on the event queue
	 * @param out the output stream to which the event queue should be dumped
	 */
	public void printQueue(PrintStream out) {
		for(Event i : queue) {
			out.println(i);
		}
	}
	
	/**
	 * Test the operation of the Scheduler class
	 * @author Saleh
	 *
	 */
	public static class Test{
		private static Event event2;
		
		private static void event1(int i, Scheduler S) {
			System.out.println("Event #" + i + " handled at time " +S.now() + 
					" remaining= " + S.getQueueLength() + " events in queue.");
			S.cancel(event2);
		}
		/**Run the Scheduler Tests
		 * @param args command line arguments
		 */
		public static void main(String[] args) {
			Scheduler S = new Scheduler();
			S.schedule(5, (s)->Test.event1(1,s), ()->"Event1");
			event2=S.schedule(20, (s)->Test.event1(2,s), ()->"Event2");
			S.schedule(3, (s)
					->s.schedule(10, (t)
							->Test.event1(4, t), ()->"Event4"), ()->"Event3");
			S.schedule(4, (s)->s.schedule(20, 
					(t)->System.out.println("Event #4 handled at time "+t.now()+ 
							" remaining= " + S.getQueueLength() + " events in queue."), ()->"Event6"), ()->"Event5");
			S.enableDebug();
			S.printQueue(System.out);
			S.handleAll();
		}
	}
}
