/**
 *  This class is a template for implementing the MCMC algorithm.  You
 *  need to fill in the constructor and the
 *  <tt>runMoreIterations()</tt> method.
 */

import java.util.*;

public class Mcmc {

    public int currentState[];
    public double counts[];
    public ArrayList<Integer> hiddenVarIndices = new ArrayList<Integer>();
    public Random r = new Random();

    public BayesNet bayesNet;
    public Query query;

    /**
     *  This is the constructor for the class that you need to fill
     *  in.  Any initialization of the MCMC algorithm should go here.
     *  The parameters to the constructor specify the Bayes net and
     *  query on which MCMC is to be run.
     */
    public Mcmc(BayesNet bn, Query q) {

        // Store bn & q
        bayesNet = bn;
        query = q;

        currentState = new int[bn.numVariables];
        counts = new double[bn.numValues[q.queryVar]];

        // Setting evidence variables, initialize random hidden variables
        for (int i = 0; i < q.evidence.length; i++) {
            if (q.evidence[i] != -1) {
                currentState[i] = q.evidence[i];
            } else {
                // System.out.println("vvvvvv hiddenVarIndices vvvvvv");
                // System.out.println(i);

                hiddenVarIndices.add(i);
                currentState[i] = r.nextInt(bn.numValues[i]);
            }
        }

        // System.out.println("CURRENT STATE:");
        // for (int s : currentState)
        //     System.out.println(s);
        // System.out.println("^^^^^^^^^^^^^^^^^^^^^^^^^^^");
    }

    /**
     *  This method, which must be filled in, runs <tt>n</tt>
     *  <i>additional</i> iterations of the MCMC algorithm on the
     *  Bayes net and query that were specified when this object was
     *  constructed.  It is important to remember that the method must
     *  <i>continue</i> a previous execution of MCMC.  It should
     *  <i>not</i> restart from scratch each time it is called.  The
     *  method returns an array with the estimated probability of each
     *  value of the query variable as estimated by MCMC following the
     *  <tt>n</tt> additional iterations.
     */
    public double[] runMoreIterations(int n) {

        for (int iter = 0; iter < n; iter++) {
            // Choose random non-evidence variable
            int chosenHiddenVar = hiddenVarIndices.get(r.nextInt(hiddenVarIndices.size()));

            // System.out.println("vvvvvv chosenHiddenVar vvvvvv");
            // System.out.println(chosenHiddenVar);

            // Store distribution for si = 0, si = 1, si = 2.... into ArrayList
            ArrayList<Double> hiddenVarProbDist = new ArrayList<Double>();

            // Copy currentState
            int currentStateCopy[] = new int[currentState.length];
            System.arraycopy(currentState, 0, currentStateCopy, 0, currentState.length);

            // Find children
            ArrayList<Integer> chosenHiddenVarChildren = new ArrayList<Integer>();
            for (int i = 0; i < bayesNet.numVariables; i++) {
                // Check if this variable has parents which is the chosenHiddenVar
                for (int p = 0; p < bayesNet.numParents[i]; p++) {
                    if (bayesNet.parents[i][p] == chosenHiddenVar) {
                        chosenHiddenVarChildren.add(i);
                    }
                }
            }

            // Compute P(si | mb(si)). Let condProb = P(si | mb(si))
            for (int i = 0; i < bayesNet.numValues[chosenHiddenVar]; i++) {
                // Change chosenHiddenVar to each new value i
                currentStateCopy[chosenHiddenVar] = i;

                // Calculate prob for this value and add to ArrayList
                double condProb = bayesNet.getCondProb(chosenHiddenVar, currentStateCopy);

                // Refresh childrenProbProduct each time
                double childrenProbProduct = 1.0;

                for (int c : chosenHiddenVarChildren) {
                    childrenProbProduct *= bayesNet.getCondProb(c, currentStateCopy);
                }

                // Update hiddenVarProbDist values to be actual values given mb(si)
                double totalProb = condProb * childrenProbProduct;
                hiddenVarProbDist.add(totalProb);
            }

            // Calculate normalization factor alpha
            double normalizationFactor = 0;
            for (double s : hiddenVarProbDist) {
                normalizationFactor += s;
            }

            normalizationFactor = 1.0 / normalizationFactor;

            // Normalize probabilities
            for (int i = 0; i < hiddenVarProbDist.size(); i++) {
                hiddenVarProbDist.set(i, hiddenVarProbDist.get(i) * normalizationFactor);
            }

            // System.out.println("hiddenVarProbDist FINAL:");
            // for (double s : hiddenVarProbDist) {
            //     System.out.println(s);
            // }
            // System.out.println("^^^^^^^^^^^^^^^^^^^^^^^^^^^");

            // Now that we have the actual probDist, sample from it!
            double randomDouble = r.nextDouble();
            double boundary = 0.0;

            // System.out.println("hiddenVarProbDist is this big:");
            // System.out.println(hiddenVarProbDist.size());

            for (int i = 0; i < hiddenVarProbDist.size(); i++) {
                boundary += hiddenVarProbDist.get(i);

                if (randomDouble < boundary) {
                    // Need to only update count if the chosen hidden variable is equal to the query variable
                    if (chosenHiddenVar == query.queryVar) {
                        counts[i] += 1.0;
                    }
                    currentState[chosenHiddenVar] = i;
                    break;
                }
            }
        }

        // Normalize counts
        double normalization = 0;
        for (double c : counts) {
            normalization += c;
        }

        normalization = 1.0 / normalization;

        for (int i = 0; i < counts.length; i++) {
            counts[i] *= normalization;
        }

        return counts;
    }
}
