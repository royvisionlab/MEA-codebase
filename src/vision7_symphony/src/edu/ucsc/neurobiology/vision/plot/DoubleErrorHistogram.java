package edu.ucsc.neurobiology.vision.plot;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DoubleErrorHistogram /*extends Tag*/
    implements HistogramData, HistogramErrors, Cloneable {

    private HistogramData data;
    private HistogramData error;


    public DoubleErrorHistogram(HistogramData data, HistogramData error) {
//        super(NeuroTagSet.ERROR_HISTOGRAM_TAG, 0);
        this.data = data;
        this.error = error;
    }


    public double getBinInterval() {
        return data.getBinInterval();
    }


    public String getDescription() {
        return data.getDescription();
    }


    public int getBinCount() {
        return data.getBinCount();
    }


    public double getMin() {
        return data.getMin();
    }


    public double getMax() {
        return data.getMax();
    }


    public double getBin(int i) {
        return data.getBin(i);
    }


    public double getBinError(int bin) {
        return error.getBin(bin);
    }


    /*
        public Tag read(int tagId, TaggedInputStream input, int len) throws IOException {
            ErrorHistogram h = (ErrorHistogram)this.clone();

            h.data = (HistogramData)input.readTag();
            h.error = (HistogramData)input.readTag();

            return h;
        }


        public void write(int tagId, TaggedOutputStream output) throws IOException {
            output.writeTag((Tag)data);
            output.writeTag((Tag)error);
        }
     */

    public Object clone() {
        try {
            return super.clone();
        } catch (CloneNotSupportedException e) {
            return null;
        }
    }


    public String toString() {
        return getDescription();
    }
}
