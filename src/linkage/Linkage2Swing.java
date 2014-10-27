package linkage; /*
 * LinkageTest.java
 *
 * Created on September 11, 2004, 4:53 PM
 */

import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;

public class Linkage2Swing {

    public static final double ENOUGH = 0.00000001;

    public static void main(String argv[]) throws IOException {
        float[] points = Linkages.readLinkage("doubleSpiral.lsf");

        CIDOUnfolder luf = new CIDOUnfolder(points, true);

        LinkageDrawer ldr = new LinkageDrawer(luf);
        Frame f = new Frame("CIDO Linkage Unfolding");
        f.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }
        });
        f.add(ldr);
        f.pack();
        f.setSize(new Dimension(400, 300));
        f.setVisible(true);

        double deltae;
        do {
            deltae = luf.step(1);
            ldr.repaint();
        } while (deltae > ENOUGH);

    }
}  