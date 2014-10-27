package linkage; /*
 * LinkageTest.java
 *
 * Created on September 11, 2004, 4:53 PM
 */

import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.*;
import java.io.*;

public class LinkageTest {
	public static final double ENOUGH = 0.00000001;
	public static final boolean KEYTIMES = false;
	public static float[] readLinkage(String filename) throws IOException {
		
		
		Vector nums = new Vector();
		BufferedReader r = new BufferedReader(new FileReader(filename));
		StreamTokenizer st = new StreamTokenizer(r);
		st.parseNumbers();
		st.eolIsSignificant(false);
		while(st.nextToken()!=StreamTokenizer.TT_EOF) {
			if(st.ttype == StreamTokenizer.TT_NUMBER) {
				nums.add(new Double(st.nval));
			}
		}
		
		int s = nums.size();
		int n = s / 2;
		float[] points = new float[s];
		for(int i=0;i<n;i++) {
			points[2*i] = ((Double) nums.elementAt(i)).floatValue();
		}
		for(int i=0;i<n;i++) {
			points[2*i+1] = ((Double) nums.elementAt(i+n)).floatValue();
		}
		return points;
	}
	
	public static void main(String argv[]) {
		float[] points = null;
		try { 
			// points = readLinkage("spider-SLOW.lsf");
			// points = readLinkage("triangleTree.lsf");
		    points = readLinkage("doubleSpiral.lsf");
			// points = readLinkage("letterW.lsf");
			// points = readLinkage("teeth.lsf");
			// points = readLinkage("fastQuadrat.lsf");
			// points = readLinkage("pineTree.lsf");
			// points = readLinkage("hookedW.lsf");
			// points = readLinkage("zickzack.lsf");
//			points = readLinkage(argv[0]);
		} catch(IOException e) {
			System.out.println("Lesefehler\n");
			System.exit(1);
		}
		CIDOUnfolder luf = new CIDOUnfolder(points,true);


		LinkageDrawer ldr = new LinkageDrawer(luf);
		Frame f = new Frame("CIDO Linkage Unfolding");
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {System.exit(0);}
			//            public void windowDeiconified(WindowEvent e) { demo.start(); }
			//            public void windowIconified(WindowEvent e) { demo.stop(); }
		});
		f.add(ldr);
		f.pack();
		f.setSize(new Dimension(400,300));
		f.setVisible(true);
		
		Thread me = Thread.currentThread();

		
		StringBuffer bbox = new StringBuffer();
		StringBuffer times = new StringBuffer();
		Vector energies = new Vector();
		boolean first = true;
		
		System.out.println("<?xml version=\"1.0\" standalone=\"no\"?>");
		System.out.println("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">");
		System.out.println("<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" preserveAspectRatio=\"xMidYMid meet\">");
//		System.out.println("<desc>Entfaltung von "+argv[0]+"</desc>");
		System.out.println("<polygon fill=\"lime\" stroke=\"blue\" stroke-width=\"0.3%\">");
		System.out.println("<animate attributeName=\"points\" attributeType=\"XML\" values=\"");
	
		double deltae;
		do {  
			
			double[] x = luf.getX();
			double[] y = luf.getY();
			// int maxx = (int) x[0];
			// int maxy = (int) y[0];
			double maxx = x[0];
			double maxy = y[0];
			double minx = maxx;
			double miny = maxy;
			int n = x.length;
			int k = luf.first;
			
			if(!first) {
				System.out.println(";");
				bbox.append(";\n");
			}
			first = false;
			
			for(int i=0;i<n;i++) {
				int j = (i+k)%n;
				System.out.print(x[j]);
				System.out.print(",");
				System.out.print(y[j]);
				System.out.print(" ");
				maxx=Math.max(maxx,x[j]);
				maxy=Math.max(maxy,y[j]);
				minx=Math.min(minx,x[j]);
				miny=Math.min(miny,y[j]);
			}
			bbox.append(minx).append(" ").append(miny).append(" ");
			bbox.append(maxx-minx).append(" ").append(maxy-miny);
			deltae = luf.step(10);
			energies.add(new Double(deltae));
			// ldr.repaint();
			// System.out.println(deltae);
			// System.out.println(luf.ellEnergy());
		} while (deltae > ENOUGH);
		double sume=0.0;
		for(Iterator i=energies.listIterator(); i.hasNext();) {
			sume+=((Double) i.next()).doubleValue();
		}
		first = true;
		for(Iterator i=energies.listIterator(); i.hasNext();) {
			if(!first) times.append("\n;");
			first=false;
			times.append(((Double) i.next()).doubleValue() / sume);
		}
		if(KEYTIMES) {
			System.out.println("\" keytimes=\"");
			System.out.print(times);
		}
		System.out.println("\" dur=\"10s\" fill=\"freeze\"/>");
		System.out.println("</polygon>");
		System.out.println("<animate attributeName=\"viewBox\" attributeType=\"XML\" values=\"");
		System.out.print(bbox);
		if(KEYTIMES) {
			System.out.println("\" keytimes=\"");
			System.out.print(times);
		}
		System.out.println("\" dur=\"10s\" fill=\"freeze\"/>");
		
		System.out.println("</svg>");
		// System.out.print("Iterationen: ");
		// System.out.println(luf.iterations);
	}
}  