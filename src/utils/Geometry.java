package utils;

public class Geometry {
	public static double distance(double[] u, double[] v) {
		return Math.sqrt(Math.pow(u[0]-v[0],2)+Math.pow(u[1]-v[1],2)+
				Math.pow(u[2]-v[2],2));
	}
}
