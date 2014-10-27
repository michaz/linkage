package linkage;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.Vector;

class Linkages {
    public static float[] readLinkage(String filename) throws IOException {
        Vector<Double> nums = new Vector<>();
        BufferedReader r = new BufferedReader(new FileReader(filename));
        StreamTokenizer st = new StreamTokenizer(r);
        st.parseNumbers();
        st.eolIsSignificant(false);
        while (st.nextToken() != StreamTokenizer.TT_EOF) {
            if (st.ttype == StreamTokenizer.TT_NUMBER) {
                nums.add(st.nval);
            }
        }
        int s = nums.size();
        int n = s / 2;
        float[] points = new float[s];
        for (int i = 0; i < n; i++) {
            points[2 * i] = nums.elementAt(i).floatValue();
        }
        for (int i = 0; i < n; i++) {
            points[2 * i + 1] = nums.elementAt(i + n).floatValue();
        }
        return points;
    }
}
