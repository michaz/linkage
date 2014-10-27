package linkage; /*
 * LinkageTest.java
 *
 * Created on September 11, 2004, 4:53 PM
 */

import java.io.IOException;
import java.util.Vector;

public class Linkage2SVG {

    public static final double ENOUGH = 0.00000001;

    public static final boolean KEY_TIMES = false;

    public static void main(String argv[]) throws IOException {
        float[] points = Linkages.readLinkage("doubleSpiral.lsf");

        CIDOUnfolder luf = new CIDOUnfolder(points, true);

        StringBuffer bbox = new StringBuffer();
        StringBuffer times = new StringBuffer();
        Vector<Double> energies = new Vector<>();
        boolean first = true;

        System.out.println("<?xml version=\"1.0\" standalone=\"no\"?>");
        System.out.println("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">");
        System.out.println("<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" preserveAspectRatio=\"xMidYMid meet\">");
        System.out.println("<polygon fill=\"lime\" stroke=\"blue\" stroke-width=\"0.3%\">");
        System.out.println("<animate attributeName=\"points\" attributeType=\"XML\" values=\"");

        double deltae;
        do {
            double[] x = luf.getX();
            double[] y = luf.getY();
            double maxx = x[0];
            double maxy = y[0];
            double minx = maxx;
            double miny = maxy;
            int n = x.length;
            int k = luf.first;

            if (!first) {
                System.out.println(";");
                bbox.append(";\n");
            }
            first = false;

            for (int i = 0; i < n; i++) {
                int j = (i + k) % n;
                System.out.print(x[j]);
                System.out.print(",");
                System.out.print(y[j]);
                System.out.print(" ");
                maxx = Math.max(maxx, x[j]);
                maxy = Math.max(maxy, y[j]);
                minx = Math.min(minx, x[j]);
                miny = Math.min(miny, y[j]);
            }
            bbox.append(minx).append(" ").append(miny).append(" ");
            bbox.append(maxx - minx).append(" ").append(maxy - miny);
            deltae = luf.step(1);
            energies.add(deltae);
        } while (deltae > ENOUGH);

        double sume = 0.0;
        for (Double energy : energies) {
            sume += energy;
        }
        first = true;
        for (Double energy : energies) {
            if (!first) times.append("\n;");
            first = false;
            times.append(energy / sume);
        }
        if (KEY_TIMES) {
            System.out.println("\" keytimes=\"");
            System.out.print(times);
        }
        System.out.println("\" dur=\"10s\" fill=\"freeze\"/>");
        System.out.println("</polygon>");
        System.out.println("<animate attributeName=\"viewBox\" attributeType=\"XML\" values=\"");
        System.out.print(bbox);
        if (KEY_TIMES) {
            System.out.println("\" keytimes=\"");
            System.out.print(times);
        }
        System.out.println("\" dur=\"10s\" fill=\"freeze\"/>");

        System.out.println("</svg>");
    }
}  