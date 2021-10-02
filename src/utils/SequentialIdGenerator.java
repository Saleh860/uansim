package utils;

import java.util.function.Supplier;

/**
 * Generates automatically incremented IDs for nodes 
 */
public class SequentialIdGenerator implements Supplier<Integer> {
	Integer nextId;

	@Override
	public Integer get() {
		return nextId++;
	}
	
	public SequentialIdGenerator() {
		nextId=0;
	}
	
	public SequentialIdGenerator(int firstId) {
		nextId=firstId;
	}
}
