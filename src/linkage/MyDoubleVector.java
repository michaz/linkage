package linkage;/*
 * Created on Nov 23, 2004
 */


class MyDoubleVector implements Cloneable {
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

    public void angleNormalize() {
		double ngrada = 0.0;
        for (double aV : v) {
            ngrada += Math.abs(aV);
        }
		for(int i=0; i<v.length; i++) {
			v[i] = v[i]/ngrada;
		}
	}

    public double[] getArray() {
		return v.clone();
	}
}
