package linkage;

import java.lang.Math;

public class CIDOUnfolder {
	
	public final int suchSchritte = 10;
	
	private double x0, y0;
	private MyDoubleVector x;
	private MyDoubleVector y;
	private MyDoubleVector lengths;
	private MyDoubleVector angles;
	public double lastVertexOrientation;
	private boolean isClosed;
	private int n;
	public int iterations = 0;

    // KLUDGE:
	public int first=0;
	
	public static double signum(double x){
		int c = Double.compare(x,0.0);
		if(c<0) return -1;
		if(c>0) return 1;
		return 0;
	}
	
	public static double hypot(double y, double x) {
		return Math.sqrt(y*y + x*x);
	}
	
	public int lastVertex() { // zwischen 1 und n
		return( x.lastVertex()   );
	}
	
	public static double minEucDistance(	MyDoubleVector sx, 
											MyDoubleVector sy,
											int n,
											boolean closed) {
		double mindist = Double.POSITIVE_INFINITY;
		for(int v=1; v<=(closed ? n :(n-1)); v++) { // Kanten
			for(int u=1; u<=n; u++) { // Ecken
				if(sx.equalIx(u,v) || sx.equalIx(u,(v+1))) continue;
				double d = 
					eucDistanceEdgeVertex(sx.getAt(v),sy.getAt(v),sx.getAt(v+1),sy.getAt(v+1),sx.getAt(u),sy.getAt(u));
				mindist = Math.min(mindist,d);
			}
		}
		return mindist;
	}
	
	public static double stepSizeBound(		MyDoubleVector sx,
											MyDoubleVector sy,
											MyDoubleVector lengths,
											MyDoubleVector grada,
											int n,
											boolean closed) {
		double virtualLength = 0.0;
		for(int i=1; i<=(n-1); i++) {
			virtualLength += Math.abs(grada.getAt(i)) * lengths.getAt(i);
		}
		double intersectionBound = minEucDistance(sx,sy,n,closed) / virtualLength;
		if(closed) {
			double overlap = overlap(sx,sy,lengths,n);
			double closureBound = overlap / virtualLength;
			return Math.min(closureBound,intersectionBound);
		} else {
			return intersectionBound;
		}
	}
	
	public static double overlap( 	MyDoubleVector sx,
									MyDoubleVector sy,
									MyDoubleVector lengths,
									int i) {
		double dx = sx.getAt(i-1)-sx.getAt(i+1);
		double dy = sy.getAt(i-1)-sy.getAt(i+1);
		double dsq = dx*dx+dy*dy;
		double d = Math.sqrt(dsq);
        double r = Math.max(lengths.getAt(i),lengths.getAt(i-1));
		double s = Math.min(lengths.getAt(i),lengths.getAt(i-1));

		double overlap = Math.min(r+s- d, d +s-r);
		return Math.max(overlap,0.0);
	}
	
	public static double ellEnergy(	MyDoubleVector sx, 
									MyDoubleVector sy,
									int n,
									boolean closed) {
		double energy = 0.0;
		for(int v=1; v<=(closed ? n :(n-1)); v++) { // Kanten
			for(int u=1; u<=n; u++) { // Ecken
				if(sx.equalIx(u,v) || sx.equalIx(u,(v+1))) continue;
				energy+=1/ellDistanceEdgeVertex(	sx.getAt(v),sy.getAt(v),
													sx.getAt(v+1),sy.getAt(v+1),
													sx.getAt(u),sy.getAt(u));
			}
		}
		return energy;
	}
	
	public static void ellEnergyGradient(	MyDoubleVector sx, 
											MyDoubleVector sy,
											MyDoubleVector angles,
											MyDoubleVector lengths,
											int n,
											boolean closed,
											MyDoubleVector grada/* out */) {
		MyDoubleVector gradpx = new MyDoubleVector(n);
		MyDoubleVector gradpy = new MyDoubleVector(n);
		for(int v=1; v<=(closed ? n :(n-1)); v++) { // Kanten
			for(int u=1; u<=n; u++) { // Ecken
				if(sx.equalIx(u,v) || sx.equalIx(u,(v+1))) continue;
				double duxvx = sx.getAt(u)-sx.getAt(v);
				double duyvy = sy.getAt(u)-sy.getAt(v);
				double duxwx = sx.getAt(u)-sx.getAt(v+1);
				double duywy = sy.getAt(u)-sy.getAt(v+1);
				double dvxwx = sx.getAt(v)-sx.getAt(v+1);
				double dvywy = sy.getAt(v)-sy.getAt(v+1);
				double duv = hypot(duxvx,duyvy);
				double duw = hypot(duxwx,duywy);
				double dvw = hypot(dvxwx,dvywy);
				double q = duv + duw - dvw;
				double p = -1 / (q*q);
				double dtux = p * (duxvx / duv + duxwx / duw);
				double dtuy = p * (duyvy / duv + duywy / duw);
				double dtvx = p * ((-duxvx) / duv - dvxwx / dvw);
				double dtvy = p * ((-duyvy) / duv - dvywy / dvw);
				double dtwx = p * ((-duxwx) / duw - (-dvxwx) / dvw);
				double dtwy = p * ((-duywy) / duw - (-dvywy) / dvw);
				gradpx.addTo(v,dtvx);
				gradpy.addTo(v,dtvy);
				gradpx.addTo(v+1,dtwx);
				gradpy.addTo(v+1,dtwy);
				gradpx.addTo(u,dtux);
				gradpy.addTo(u,dtuy);
			}
		}
		for(int i=(n-1); i>=1; i--) {
			gradpx.addTo(i,gradpx.getAt(i+1));
			gradpy.addTo(i,gradpy.getAt(i+1));
		}			
		for(int i=1; i<=(n-1); i++) {
			grada.setAt(i, 
				(-Math.sin(angles.getAt(i))*gradpx.getAt(i+1)
				 +Math.cos(angles.getAt(i))*gradpy.getAt(i+1)) * lengths.getAt(i)
			);
		}
		if(closed) projectForClosure(sx,sy,angles,lengths,n,grada);
	}
	
	
	// Euclidean distance
	public static double eucDistanceEdgeVertex(double vx,double vy,double wx,double wy,double ux,double uy) {
		double vvx = wx-vx; double vvy = wy-vy;
		double wwx = ux-vx; double wwy = uy-vy;
		double c1 = vvx*wwx + vvy*wwy;
		double c2 = vvx*vvx + vvy*vvy;
		if(c1<=0)
			return hypot(ux-vx,uy-vy);
		if(c2<=c1)
			return hypot(ux-wx,uy-wy);
		double pbx = vx + (c1/c2) * vvx;
		double pby = vy + (c1/c2) * vvy;
		return hypot(ux-pbx,uy-pby);
	}
	
	// Elliptic distance
	public static double ellDistanceEdgeVertex(double vx,double vy,double wx,double wy,double ux,double uy) {
		return hypot(uy-vy,ux-vx) + hypot(uy-wy,ux-wx) - hypot(vy-wy,vx-wx);
	}

	// Convert angle/length parametrization to point representation
	public static void anglesToPoints(	MyDoubleVector angles,
								MyDoubleVector lengths,
								int n,
								boolean closed,
								double x0,
								double y0,
								double lastVertexOrientation,
								MyDoubleVector nx,
								MyDoubleVector ny) {
		nx.setAt(1,x0);
		ny.setAt(1,y0);
		for(int i=2; i<=(closed ? (n-1) : n); i++) {
			nx.setAt(i,nx.getAt(i-1) + lengths.getAt(i-1)*Math.cos(angles.getAt(i-1)));
			ny.setAt(i,ny.getAt(i-1) + lengths.getAt(i-1)*Math.sin(angles.getAt(i-1)));
		}
		if(closed) { // Calculate position of nth point
			double dx = nx.getAt(n-1)-nx.getAt(1);
			double dy = ny.getAt(n-1)-ny.getAt(1);
			double dsq = dx*dx+dy*dy;
            double cp = Math.pow(lengths.getAt(n),2)-Math.pow(lengths.getAt(n-1),2)+dsq;
			double c = cp / (2*dsq);
			double ax = nx.getAt(1) + c*dx;
			double ay = ny.getAt(1) + c*dy;
			double er = Math.pow(lengths.getAt(n),2)/dsq-Math.pow(c,2);
			double e = Math.sqrt(er);
			// DEBUG:
			if(Double.isNaN(e)) {
				System.out.println("Last vertex panic!");
			}
			double vnx = ax + lastVertexOrientation * (-dy) * e;
			double vny = ay + lastVertexOrientation * dx * e;	
			nx.setAt(n,vnx);
			ny.setAt(n,vny);
		}
	}

	public static void projectForClosure(MyDoubleVector x,
								MyDoubleVector y,
								MyDoubleVector angles,
								MyDoubleVector lengths,
								int n,
								MyDoubleVector grada) {
		MyDoubleVector p = new MyDoubleVector(n-1);
		double dx = x.getAt(n) - x.getAt(1);
		double dy = y.getAt(n) - y.getAt(1);
		double d = 0.0;
		double pvalsq = 0.0;
		// normal of projection plane
		for(int i=1;i<=(n-1);i++) {
			double di = lengths.getAt(i)*( -Math.sin(angles.getAt(i))*dx 
					                       +Math.cos(angles.getAt(i))*dy );
			p.setAt(i,di);
			d += grada.getAt(i) * di;
			pvalsq += Math.pow(p.getAt(i),2);
		}
		d=d/pvalsq;
		for(int i=1;i<=(n-1);i++) {
			grada.setAt(i,grada.getAt(i)- d*p.getAt(i));
		}
	}
	
	
	public CIDOUnfolder(float[] points, boolean b) {
		n = points.length/2;
		x = new MyDoubleVector(n);
		y = new MyDoubleVector(n);
		isClosed = b;
		
		for(int i=1; i<=n; i++) {
			x.setAt(i,points[2*(i-1)]);
			y.setAt(i,points[2*(i-1)+1]);
		}
		double xx; double yy;
		x0 = x.getAt(1); 
		y0 = y.getAt(1);
		if(isClosed) {
			angles = new MyDoubleVector(n);
			lengths = new MyDoubleVector(n);
		} else {
			angles = new MyDoubleVector(n-1);
			lengths = new MyDoubleVector(n-1);
		}
		for(int i=1;i<=n-1;i++) {
			xx = x.getAt(i+1) - x.getAt(i);
			yy = y.getAt(i+1) - y.getAt(i);
			lengths.setAt(i,hypot(yy,xx));
			angles.setAt(i,Math.atan2(yy,xx));
		}
		if(isClosed) {
			xx = x.getAt(1) - x.getAt(n);
			yy = y.getAt(1) - y.getAt(n);
			lengths.setAt(n,hypot(yy,xx));
			angles.setAt(n,Math.atan2(yy,xx));
		}
		if(isClosed) {
			lastVertexOrientation=signum(
					(x.getAt(n)-x.getAt(n-1)) * (y.getAt(1)-y.getAt(n-1))
					+(y.getAt(n)-y.getAt(n-1)) * (-x.getAt(1)+x.getAt(n-1)));
		}
		anglesToPoints(angles,lengths,n,isClosed,x0,y0,lastVertexOrientation,x,y);
	}

	public double step(int n) {
		double deltae=0.0;
		for(int i=0;i<n && !finished();i++) {
			deltae+=step();
		}
		return deltae;
	}
	
	public double step() {
		double energy=ellEnergy(x,y,n,isClosed);
		
		MyDoubleVector grada = new MyDoubleVector(n-1);
		if(isClosed) {
			fixNewPoint();
		} 
		ellEnergyGradient(x,y,angles,lengths,n,isClosed,grada);
		grada.angleNormalize();
		double deltat=stepSizeBound(x,y,lengths,grada,n,isClosed);
		MyDoubleVector newangles = new MyDoubleVector(isClosed?n:n-1);
		MyDoubleVector nx = new MyDoubleVector(n);
		MyDoubleVector ny = new MyDoubleVector(n);
		for(int i=1; i<=(isClosed ? (n-2): (n-1)); i++) {
			double deltai = grada.getAt(i) * deltat;
			if(Double.isNaN(deltai)) deltai = 0; 
			newangles.setAt(i,angles.getAt(i)-deltai);
			if(newangles.getAt(i) < (-Math.PI)) newangles.addTo(i,2*Math.PI);
			if(newangles.getAt(i) > Math.PI) newangles.addTo(i,-2*Math.PI);
		}
		anglesToPoints(newangles,lengths,n,isClosed,x0,y0,lastVertexOrientation,nx,ny);			
		
		double ea = 0.0;
		double eb = energy - ellEnergy(nx,ny,n,isClosed);
		double a = 0.0;
		double b = deltat;
		MyDoubleVector testangles = new MyDoubleVector(isClosed?n:n-1);
		MyDoubleVector testx = new MyDoubleVector(n);
		MyDoubleVector testy = new MyDoubleVector(n);
		double teste;
		double d = deltat * 0.5;
		for(int h=0; h<suchSchritte;h++) {
			for(int i=1; i<=(isClosed ? (n-2): (n-1)); i++) {
				double deltai = grada.getAt(i) * (a+d);
				if(Double.isNaN(deltai)) deltai = 0; 
				testangles.setAt(i,angles.getAt(i)-deltai);
				if(testangles.getAt(i) < (-Math.PI)) testangles.addTo(i,2*Math.PI);
				if(testangles.getAt(i) > Math.PI) testangles.addTo(i,-2*Math.PI);
			}
			anglesToPoints(testangles,lengths,n,isClosed,x0,y0,lastVertexOrientation,testx,testy);
			teste=energy - ellEnergy(testx,testy,n,isClosed);
			if(eb>ea) {
				if(teste<ea) {
					break;
				} else {
					a+=d;
					ea=teste;
				}
			} else{
				if(teste<eb) {
					break;
				} else {
					b-=d;
					eb=teste;
				}
			}
			d*=0.5;
		}
		double deltae;
		if(ea>eb) { 
			deltat=a; 
			deltae=ea;
		}
		else {
			deltat=b;
			deltae=eb;
		}
		for(int i=1; i<=(isClosed ? (n-2): (n-1)); i++) {
			double deltai = grada.getAt(i) * deltat;
			if(Double.isNaN(deltai)) deltai = 0; 
			newangles.setAt(i,angles.getAt(i)-deltai);
			if(newangles.getAt(i) < (-Math.PI)) newangles.addTo(i,2*Math.PI);
			if(newangles.getAt(i) > Math.PI) newangles.addTo(i,-2*Math.PI);
		}
		anglesToPoints(newangles,lengths,n,isClosed,x0,y0,lastVertexOrientation,nx,ny);
		
		if(isClosed) { // calculate last two angles
			double xx = nx.getAt(1) - nx.getAt(n);
			double yy = ny.getAt(1) - ny.getAt(n);
			newangles.setAt(n,Math.atan2(yy,xx));
			xx = nx.getAt(n) - nx.getAt(n-1);
			yy = ny.getAt(n) - ny.getAt(n-1);
			newangles.setAt(n-1,Math.atan2(yy,xx));
		}
		angles = newangles;
		x = nx;
		y = ny;
		iterations++;
		return deltae;
	}
	
	public void fixNewPoint() {
		// absolute turn angle nearest to 90 degrees
		int besti = 0; double mindiff = Math.PI/2;
		for(int i=1; i<=n; i++) {
			double turnangle = Math.abs(angles.getAt(i-1)-angles.getAt(i));
			if(turnangle > Math.PI) turnangle = 2*Math.PI - turnangle;
			double anglediff = Math.abs(Math.PI/2 - turnangle);
			if(anglediff < mindiff) {
				mindiff = anglediff;
				besti = i;
			}
		}
		setSpecialPoint(besti);
	}
	
	public void setSpecialPoint(int i) {
		x0 = x.getAt(i+1);
		y0 = y.getAt(i+1);
		lastVertexOrientation=signum(
				 (x.getAt(i)-x.getAt(i-1)) * (y.getAt(i+1)-y.getAt(i-1))
				+(y.getAt(i)-y.getAt(i-1)) * (-x.getAt(i+1)+x.getAt(i-1)));
		x.setFirst(i+1);
		y.setFirst(i+1);
		// KLUDGE:
		first = (first - i + n) % n;
		
		angles.setFirst(i+1);
		lengths.setFirst(i+1);
	}
	
	public boolean finished() {
		for(int i=1; i<=(n-1); i++) {
			if(Math.abs(angles.getAt(1)-angles.getAt(i)) > 0.1)
				return false;
		}
		return true;
	}
	
	public boolean isClosed() {
		return isClosed;
	}
	
	public double[] getX() {
		return x.getArray();
	}
	
	public double[] getY() {
		return y.getArray();
	}

}
	
