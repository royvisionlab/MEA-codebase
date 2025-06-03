


package edu.ucsc.neurobiology.vision.test;

import java.io.*;

import java.awt.*;

import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.math.MathUtil;
import java.util.ArrayList;


/**
 *
 * @author Matthew Grivich
 */
public class SpectrumGenerator {
    int minWavelength, maxWavelength, interval;
    static double[][] spectrum;
    static double[] wavelengths;

    public SpectrumGenerator(int minWavelength, int maxWavelength, int interval) {
        this.minWavelength = minWavelength;
        this.maxWavelength = maxWavelength;
        this.interval = interval;

        this.spectrum = new double[3][ (maxWavelength - minWavelength) / interval + 1];
        wavelengths = new double[spectrum[0].length];
        for (int i = 0; i < spectrum[0].length; i++) {
            wavelengths[i] = minWavelength + i * interval;
        }

    }


    public void calculateGovardovskii(double[] peaks, boolean addBeta) {
        double A = 69.7;
        double B = 28;
        double b = .922;
        double C = -14.9;
        double c = 1.104;
        double D = .674;
        double ABeta = .26;

        for (int col = 0; col < 3; col++) {
            if (peaks[col] != 0.0) {
                double a = .8795 +
                           .0459 *
                           Math.exp( - (peaks[col] - 300.0) * (peaks[col] - 300.0) /
                                    11940.0);
                double lambdaMBeta = 189.0 + .315 * peaks[col];
                double bBeta = -40.5 + .195 * peaks[col];
                for (int i = 0; i < spectrum[0].length; i++) {
                    double x = peaks[col] / wavelengths[i];
                    spectrum[col][i] = 1.0 /
                                       (Math.exp(A * (a - x)) + Math.exp(B * (b - x)) +
                                        Math.exp(C * (c - x)) + D);
                    if (addBeta) {
                        spectrum[col][i] += ABeta *
                            Math.exp( - (wavelengths[i] - lambdaMBeta) *
                                     (wavelengths[i] - lambdaMBeta) / bBeta / bBeta);
                    }
                }
            }

        }

    }


    public void calculateBaylor(double peaks[]) {
        double[] a = { -5.2734, -87.403, 1228.4, -3346.3, -5070.3, 30881, -31607};
        for (int col = 0; col < 3; col++) {
            for (int i = 0; i < spectrum[0].length; i++) {
                spectrum[col][i] = 0;
                double x = 1000.0 * peaks[col] / wavelengths[i] / 561.0;

                for (int j = 0; j < 7; j++) {
                    spectrum[col][i] += a[j] * Math.pow(Math.log(x) / Math.log(10.0), j);
                }
                spectrum[col][i] = Math.pow(10, spectrum[col][i]);
            }

            //convert from photon based to energy based
            //Before, spectrum has units of signal/photon
            //After, has units of signal/watt
            //Constants normalized out.
//            for (int i = 0; i < spectrum[0].length; i++) {
//                spectrum[col][i] *= wavelengths[i];
//            }
            MathUtil.divide(spectrum[col], MathUtil.max(spectrum[col]));
        }

//        for(int i=0; i<spectrum[0].length; i++) {
//            System.out.println(spectrum[0][i] + "   " + spectrum[1][i] + "   " + spectrum[2][i]);
//        }


    }


    public static PlotPanel makePlot(double[][] spectrum, Color firstColor) {
        PlotPanel plotPanel = new PlotPanel();
        Color[] rgb = {firstColor, Color.green, Color.blue};
        double[] errors = new double[spectrum[0].length];

        for (int cIndex = 0; cIndex < 3; cIndex++) {
            ScatterPlot scatter = new ScatterPlot(wavelengths, spectrum[cIndex], errors);
            plotPanel.addData(scatter, new ScatterPlotStyle(
                SymbolType.NONE, 3, rgb[cIndex], true, rgb[cIndex], 1));
        }
        plotPanel.autoscale();
        PlotUtil.showData("", plotPanel);
        return plotPanel;
    }


    public void write(String fileName) {
        String fieldSeparator = "\t";
        String newline = System.getProperty("line.separator");
        try {
//            File f = new File(fileName);

            RandomAccessFile file = new RandomAccessFile(fileName, "rw");
            //If you overwrite a previous file with the same name, this line is necessary.
            file.setLength(0);
            for (int i = 0; i < spectrum[0].length; i++) {
                file.writeBytes(spectrum[0][i] + fieldSeparator + spectrum[1][i] +
                                fieldSeparator + spectrum[2][i] + newline);
            }
            file.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    public static double[][] readSpectralFile(String fileName) {
        DataInputStream dis = null;
        String record = null;
        int recCount = 0;
        double[][] data = null;
        try {

            File f = new File(fileName);
            FileInputStream fis = new FileInputStream(f);
            BufferedReader br
                = new BufferedReader(new InputStreamReader(fis));
            ArrayList<String> lines = new ArrayList<String> ();
            while ( (record = br.readLine()) != null) {
                if (record.startsWith(";"))continue;
                lines.add(record);
            }
            data = new double[3][lines.size()];
            for (int i = 0; i < lines.size(); i++) {
                String[] vals = lines.get(i).split("\t");
                for (int j = 0; j < 3; j++) {
                    data[j][i] = new Double(vals[j]).doubleValue();
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
//        show("", data);
        return data;
    }


    public static void main(String args[]) throws Exception {

        /*     SpectrumGenerator generator = new SpectrumGenerator(370, 730, 1);
             double[] peaks = {0, 529.9, 400.4};
             generator.calculateGovardovskii(peaks, true);
             generator.show("GPGreenBlue.txt");
             generator.write("c:\\GPGreenBlue.txt");
                     double[] peaks2 = {499.3, 529.9, 400.4};

                     generator.calculateGovardovskii(peaks2, true);
             generator.show("GPRodGreenBlue.txt");
             generator.write("c:\\GPRodGreenBlue.txt");
         */
//        generator.calculateBaylor(peaks);
//        generator.show("Baylor");
//         generator.write("c:\\baylor.txt");

        SpectrumGenerator generator = new SpectrumGenerator(370, 730, 1);
        double[][] emissionSpectrum = generator.readSpectralFile(
            "C:\\Documents and Settings\\mgrivich\\Desktop\\calibration\\2000-09-04-1\\RGB.1");

//        double[] peaks = {561, 531, 430}; //Macaque
        double[] peaks = {500, 530, 400};  //Guinea Pig
        generator.calculateGovardovskii(peaks, true);
        PlotPanel p = generator.makePlot(spectrum, Color.black);

        p.setLabelFont(GrivichThesis.labelFont);
        p.getXAxis().setLabel("Wavelength (nm)");
        p.getYAxis().setLabel("Normalized Photon Absorption Probability");
        p.axesBorder.setTopPadding(4);
//        p.saveAsEPS(new File(
//      "g:\\Publications\\thesis\\figures\\appendix\\AbsorptionSpectrum.eps"), 4, 2, true);
/////

//        generator.calculateBaylor(peaks);
//        spectrum = generator.readSpectralFile("C:\\Documents and Settings\\mgrivich\\Desktop\\calibration\\mf-lms.1");
//        generator.makePlot(spectrum, Color.red);

//        double[] emissionIntensities = {2.99E1/4.0, 3.88E1/4.0, 4.18E1/4.0};
//        double radius = 1342; //in microns

        //Divide by one two for white-> grey and one two for ndf
        double[] emissionIntensities = {1.58E1 / 2.0 /2.0, 1.95E1 / 2.0/2.0, 2.21E1 / 2.0/2.0}; //in watts
        double area = 1.93E6; //in square microns

        //Change emission spectrum to Watts/wavelength
        for (int i = 0; i < 3; i++) {

            MathUtil.multiply(emissionSpectrum[i],
                              emissionIntensities[i] / MathUtil.sum(emissionSpectrum[i]));
            System.out.println(MathUtil.sum(emissionSpectrum[i]));
        }

        //Change emission spectrum to Photons/wavelength

        for(int i=0; i<3; i++) {
            for(int j=0; j<emissionSpectrum[0].length; j++) {
                emissionSpectrum[i][j] *= wavelengths[j]*1E-9/1000/6.625E-34/3E14/area;

            }
        }

         p = generator.makePlot(emissionSpectrum, Color.red);

                  p.setLabelFont(GrivichThesis.labelFont);
             p.getXAxis().setLabel("Wavelength (nm)");
             p.getYAxis().setLabel("Photons/\u03bcm\u00b2");
             p.axesBorder.setTopPadding(4);

//             p.saveAsEPS(new File(
//           "g:\\Publications\\thesis\\figures\\appendix\\EmissionSpectrum.eps"), 4, 2, true);
//////
        for (int color = 0; color < spectrum.length; color++) {
            double absorbedPhotons = 0.0;
            for (int i = 0; i < spectrum[color].length; i++) {
                absorbedPhotons += spectrum[color][i] * (emissionSpectrum[0][i] +
                    emissionSpectrum[1][i] + emissionSpectrum[2][i]);
            }
            //Convert to photons at peak
            System.out.println(color + ": " +
                               absorbedPhotons /* 1E-9 * peaks[color] / 1000 /
                               6.625E-34 / 3.0E14 */);
        }
    }


}
