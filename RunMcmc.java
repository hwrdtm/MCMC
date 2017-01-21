import java.io.*;
import java.text.*;

/**
 *  This class contains a <tt>main</tt> for running MCMC and printing
 *  the results periodically.  You are welcome to write your own
 *  <tt>main</tt>, but your MCMC code must also work correctly with
 *  the <tt>main</tt> provided here.
 */
public class RunMcmc {

    /**
     *  This is a simple <tt>main</tt> that you can use for testing.
     *  Once MCMC has been implemented, this <tt>main</tt> simply runs
     *  MCMC for some number of rounds, printing partial results with
     *  a specified frequency, on a given Bayes net and query.  The
     *  command-line arguments are: name of the file describing the
     *  Bayes net; name of the file describing the query; maximum
     *  number of iterations; frequency with which results should be
     *  printed.
     */
    public static void main(String argv[])
	throws FileNotFoundException, IOException {

	int max_iters = 0;
	int freq = 0;
	String bn_filename = null;
	String query_filename = null;

	try {
	    bn_filename = argv[0];
	    query_filename = argv[1];
	    max_iters = Integer.parseInt(argv[2]);
	    freq = Integer.parseInt(argv[3]);
	} catch(Exception e) {
	    System.err.println("Arguments: <bn_filename> <query_filename> <max_iterations> <printing_freq>");
	    return;
	}

	BayesNet bn = new BayesNet(bn_filename);

	Query q = new Query(bn, query_filename);


    // System.out.println(q.queryVar);

    // for (int evidence : q.evidence)
    //     System.out.println(evidence);

    System.out.println("Loaded BayesNet & Query.");

	Mcmc mcmc = new Mcmc(bn, q);

    System.out.println("MCMC Initialization was OK.");

	NumberFormat nf = new DecimalFormat("0.00000000");

	int max_len = Integer.toString(max_iters).length();

	for(int t = freq; t <= max_iters; t += freq) {
	    double[] p = mcmc.runMoreIterations(freq);
	    for(int i = Integer.toString(t).length(); i < max_len; i++)
		System.out.print(" ");
	    System.out.print(t + "   ");
	    for (int i = 0; i < bn.numValues[q.queryVar]; i++)
		System.out.print(nf.format(p[i]) + "  ");
	    System.out.println();
	}

    }

}
