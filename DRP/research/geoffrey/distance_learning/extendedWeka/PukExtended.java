import weka.classifiers.functions.supportVector.Puk;
import weka.core.Instances;


public class PukExtended
    extends Puk {

    /**
    * default constructor - does nothing.
    */
    public PukExtended() {
        super();
    }
  
    /**
    * Constructor. Initializes m_kernelPrecalc[].
    * 
    * @param data	the data to use
    * @param cacheSize	the size of the cache
    * @param omega	the exponent
    * @param sigma	the bandwidth
    * @throws Exception	if something goes wrong
    */
    public PukExtended(Instances data, int cacheSize, double omega, double sigma)
        throws Exception {
    
    super(data, cacheSize, omega, sigma);
    
    }

    
    /**
    * Calculates a dot product between two instances in the new metric space
    * 
    * @param inst1	the first instance
    * @param inst2	the second instance
    * @return 		the dot product of the two instances.
    * @throws Exception	if an error occurs
    */
    protected final double newDotProd(Instance inst1, Instance inst2)
        throws Exception {
    
        double result = 0;
    
        // we can do a fast dot product
        int n1 = inst1.numValues();
        int n2 = inst2.numValues();
        int classIndex = m_data.classIndex();
        for (int p1 = 0, p2 = 0; p1 < n1 && p2 < n2;) {
            int ind1 = inst1.index(p1);
            int ind2 = inst2.index(p2);
            if (ind1 == ind2) {
                if (ind1 != classIndex) {
                    result += inst1.valueSparse(p1) * inst2.valueSparse(p2);
                }
                p1++;
                p2++;
            } else if (ind1 > ind2) {
                p2++;
            } else {
                p1++;
            }
        }
        return (result);
    }

    protected double evaluate(int id1, int id2, Instance inst1)
    throws Exception {

        if (id1 == id2) {
            return 1.0;
        } else {
            double precalc1;
            if (id1 == -1){
                precalc1 = newDotProd(inst1, inst1);
            } else{
                precalc1 = m_kernelPrecalc[id1];
            }
            Instance inst2 = m_data.instance(id2);
            double squaredDifference = -2.0 * newDotProd(inst1, inst2) + precalc1 + m_kernelPrecalc[id2];
            double intermediate = m_factor * Math.sqrt(squaredDifference);
            double result = 1.0 / Math.pow(1.0 + intermediate * intermediate, getOmega());
            return result;
        }
    }
}
    
