package linkage;/*
 * Created on 24.11.2004
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
public class VectorTest {

	public static void main(String argv[]) {
		MyDoubleVector dv = new MyDoubleVector(6);
		dv.setAt(1,1.0);
		dv.setAt(2,2.0);
		dv.setAt(3,3.0);
		dv.setAt(4,4.0);
		dv.setAt(5,5.0);
		dv.setAt(6,6.0);
		
		dv.addTo(1,4.0);
		dv.addTo(2,4.0);
		dv.addTo(3,4.0);
		dv.addTo(4,4.0);
		dv.addTo(5,4.0);
		dv.addTo(6,4.0);
		
		
		
	}
	
}
