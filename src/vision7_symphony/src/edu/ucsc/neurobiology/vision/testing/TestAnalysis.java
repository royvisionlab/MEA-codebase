package edu.ucsc.neurobiology.vision.testing;

import java.io.IOException;

import edu.ucsc.neurobiology.vision.analysis.MovingBar;
import edu.ucsc.neurobiology.vision.io.ParametersFile;


/**
 * @author nobody, anyone can change
 */
public class TestAnalysis {




    public static void movingBarTest() throws IOException {
        String[] names = {"data017", "data016", "data015", "data014", "data013"};

        MovingBar[] b = new MovingBar[names.length];
        for (int i = 0; i < names.length; i++) {
            b[i] = new MovingBar("f:\\Good Data\\2005-04-26-0", names[i]);
        }

        ParametersFile pf = new ParametersFile(
            "f:\\Good Data\\2005-04-26-0\\data009\\data009.params");
        int[] id = pf.getNeuronsInClass("All/LED");

        for (int i = 0; i < id.length; i++) {
            MovingBar.showResponse1(id[i], b);
        }

//        MovingBar.showTunning(new int[] {id}, b, "");
        /*
                ParametersFile pf = new ParametersFile(
                    "f:\\data\\2005-04-26-0\\data009\\data009.params");
                String[] c = {
                             "All/OFF/Amacrine",
                             "All/ON/Parasol",
                             "All/OFF/Parasol",
                             "All/OFF/YFast",
                };

                for (int i = 0; i < c.length; i++) {
//            MovingBar.showTunning(pf.getNeuronsInClass(c[i]), b, c[i]);
                    MovingBar.showAverageTunning(pf.getNeuronsInClass(c[i]), b, c[i]);
                }
         */
    }


    public static void main(String[] args) throws IOException {
        movingBarTest();
    }
}
