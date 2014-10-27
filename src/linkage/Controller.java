package linkage;

import javafx.fxml.FXML;
import javafx.scene.layout.Pane;
import javafx.scene.shape.Polyline;

import java.io.IOException;

public class Controller {

    @FXML
    Pane pane;

    @FXML
    void initialize() throws IOException {
        float[] points = Linkages.readLinkage("doubleSpiral.lsf");
        CIDOUnfolder luf = new CIDOUnfolder(points, true);
        Polyline path = new Polyline();
        double[] x=luf.getX();
        double[] y=luf.getY();
        int n = x.length;
        path.getPoints().addAll(x[0], y[0]);
        for (int i = 1; i < n; i++) {
            path.getPoints().addAll(x[i], y[i]);
        }
        if(luf.isClosed()) {
            path.getPoints().addAll(x[0], y[0]);
        }
        path.setStrokeWidth(0.001);
        path.setScaleX(100);
        path.setScaleY(100);
        pane.getChildren().add(path);
    }

}
