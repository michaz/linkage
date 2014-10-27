package linkage;/*
 * Created on Nov 23, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

/**
 * @author zilske
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class MyDoubleVector implements Cloneable {
	private double[] v;
	public int first;
	public int lastVertex() {
		return( ((v.length+first-2) % v.length)+1);
	}
	public MyDoubleVector(int n) {
		v = new double[n];
		first = 1;
	}
	public double getAt(int n) {
		return v[(v.length+n+first-2) % v.length];
	}
	public boolean equalIx(int m, int n) {
			return (m % v.length == n % v.length);
	}
	public void setAt(int n,double d) {
		v[(v.length+n+first-2) % v.length] = d;
	}
	public void addTo(int n,double d) {
		v[(v.length+n+first-2) % v.length] += d;
	}
	public void setFirst(int n) {
		first = first+n-1;
	}
	public int length() {
		return v.length;
	}
	public void empty() {
		for(int i=0; i<v.length; i++) {
			v[i] = 0.0;
		}
	}
	public double max() {
		double m=0.0;
		for(int i=0; i<v.length; i++) {
			m=Math.max(m,v[i]);
		}
		return m;
	}
	
	public void angleNormalize() {
		double ngrada = 0.0;
		for(int i=0; i<v.length; i++) {
			ngrada += Math.abs(v[i]);
		}
		for(int i=0; i<v.length; i++) {
			v[i] = v[i]/ngrada;
		}
	}
	
	public void normalize() {
		double norm= 0.0;
		for(int i=0; i<v.length; i++) {
			norm += Math.pow(v[i],2);
		}
		for(int i=0; i<v.length; i++) {
			v[i] = v[i]/norm;
		}
	}
	
	public void resizeAndEmpty(int nn) {
		if(nn != v.length)
			v = new double[nn];
	}
	public double[] getArray() {
		return (double[]) v.clone();
	}
}
