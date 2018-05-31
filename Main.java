import java.util.*;
import java.io.*;
import java.lang.*;

class Tuple {
  double prob; 
  boolean label;
  int weight;
  public Tuple(double prob, boolean label, int weight) {
    this.prob = prob; this.label = label; this.weight = weight;
  }
}

public class Main {
  public static int N, E, K = 4;
  public static double sgdStepSize = 0.01;    // may be different for different data
  public static final int negDictSize = 100000000;  // 100M
  public static String posFilename;

  public static double paramA = 1.0, paramB = 2.0, paramC = -1.0;
  public static double rAlpha = 1.5;
  public static double rLimit = 10;
  public static boolean WEIGHTED = false;
  public static double thetaReg = 1e-6;
  public static double tolerance = 1e-6;
  public static int MAX_EDGES = (int)1e8;
  public static boolean verbose = true;
  public static double tp = 0.9;
  public static final int shuffleSeed = 42;
  public static String outPrefix = "./";

  public static double nsw;
  public static Random rand;

  public static String outFilenameR, outFilenameTheta, outFilename1DTheta;
  public static List<Integer> edgeSources, edgeTargets, weights;
  public static List<Boolean> isRecip;
  public static Map<String, Integer> map;
  public static Map<Integer, String> invMap;
  public static int[] negDict = new int[negDictSize];
  public static Map<Integer, Integer> outDegree;

  public static void init() {
    rand = new Random(200);
    map = new HashMap<String, Integer>();
    invMap = new HashMap<Integer, String>();
    edgeSources = new ArrayList<Integer>();
    edgeTargets = new ArrayList<Integer>();
    weights = new ArrayList<Integer>();
    outDegree = new HashMap<Integer, Integer>();
    isRecip = new ArrayList<Boolean>();
    outFilenameR = outPrefix + "/res_r.txt";
    outFilenameTheta = outPrefix + "/res_z.txt";
    outFilename1DTheta = outPrefix + "/res_theta.txt";
    System.out.printf("[Info]: results will be saved to %s and %s\n", outFilenameR, outFilenameTheta);

    try {
      N = FileParser.readCSVDict(posFilename, "", map, invMap);

      E = FileParser.readCSVGraph(posFilename, map, edgeSources, edgeTargets, weights, isRecip, N, negDict, outDegree, WEIGHTED);

      // use readCSVGraphApprox when there are too many edges (i.e. E > 100M)
      //E = FileParser.readCSVGraphApprox(posFilename, map, edgeSources, edgeTargets, weights, isRecip, N, negDict, outDegree, WEIGHTED);

      nsw = 5;
    } catch (IOException e) {
      e.printStackTrace();
      System.exit(0);
    } 

    if (verbose) {
      System.out.printf("[Info] Number of nodes = %d\n", N);
      System.out.printf("[Info] Number of edges (pos) = %d, edges (neg) = %d\n", E, edgeSources.size()-E);
      System.out.printf("[Info] Density = %g, Negative sample weight = %g\n", 1.0*E/N/(N-1), nsw);
    }
  }

  /** 
   * calculate objective function (average probability)
   */
  public static double calcObj(double[] R, double[][] theta, List<Integer> allEdges) {
    double res = 0.0;
    double pos_p = 0.0, neg_p = 0.0;
    int pos_count = 0, neg_count = 0;
    double Z = Math.pow(rLimit, rAlpha+1) / (rAlpha+1);

    for (int _e = 0; _e < 2*E*tp; _e++) {
      int e = allEdges.get(_e);
      int s = edgeSources.get(e), t = edgeTargets.get(e);   // s cites t
      double lw = 1.0;
      double dTheta = l2dist(theta[s], theta[t]);
      double dij = paramA * (R[s]-R[t]) * (1-1/(1+dTheta)) - paramB * dTheta + paramC;
      double pij = logis(dij);
      if (e < E) {
	double log_pij = Math.log(pij) * lw;
	log_pij *= weights.get(e);
	pos_count += weights.get(e);
	pos_p += log_pij;
	res += log_pij;
      } else {
	double log_pij = Math.log(1-pij) * lw;
	log_pij *= weights.get(e);
	neg_count += weights.get(e);
	neg_p += log_pij;
	res += nsw * log_pij;
      }
    }
    if (pos_count+neg_count != 0 && verbose) {
      System.out.printf("Aver of pos = %f\n", pos_p/pos_count);
      System.out.printf("Aver of neg = %f, (without nsw) = %f\n", neg_p/neg_count*nsw, neg_p/neg_count);
      System.out.printf("Before regularization: %f\n", res/(pos_count+nsw*neg_count));
    }

    for (int n = 0; n < N; n++) {
      res += Math.pow(R[n], rAlpha) / Z;
      for (int k = 0; k < K; k++) {
	res += 0.5 * thetaReg * theta[n][k] * theta[n][k];
      }
    }

    if (pos_count+neg_count != 0) return res/(pos_count+nsw*neg_count);
    else return -1;
  }

  public static long nextLong(Random rng, long n) {
   // error checking and 2^x checking removed for simplicity.
   long bits, val;
   do {
      bits = (rng.nextLong() << 1) >>> 1;
      val = bits % n;
   } while (bits-val+(n-1) < 0L);
   return val;
  }

  /* return sigmoid(x) */
  public static double logis(double x) {
    if (x < -20) return 1e-9;
    else {
      double v = 1.0 / (1.0 + Math.exp(-x));
      if (v > 1-1e-9) return 1-1e-9;
      else return v;
    }
  }

  
  /* update parameters for a single link */
  public static void update(double[] R, double[][] theta, int i, int j, boolean label, double[][] grad) {
    /*
     * [output]
     *	  grad[0][0~(K-1)]: gradient w.r.t. theta_i
     *	  grad[0][K]: gradient w.r.t. R_i 
     *	  grad[1][0~(K-1)]: gradient w.r.t. theta_j
     *	  grad[1][K]: gradient w.r.t. R_j 
     */ 
    double lw = 1.0;
    double dTheta = l2dist(theta[i], theta[j]);
    double dij = paramA * (R[i]-R[j]) * (1-1/(1+dTheta)) - paramB * dTheta + paramC;
    double pij = logis(dij);
    int yij = label ? 1 : 0;

    double useNsw = label ? 1 : nsw;

    for (int k = 0; k < K+1; k++) {
      if (k < K) {
	// dTheta
	double g = (yij-pij) * ( 2 * paramA * (R[i]-R[j]) * (theta[i][k]-theta[j][k]) / ((1+dTheta)*(1+dTheta))
	  - 2 * paramB * (theta[i][k]-theta[j][k]) );
	grad[0][k] += g * lw * useNsw;
	grad[1][k] -= g * lw * useNsw;

	if (label)
	  grad[0][k] -= thetaReg * theta[i][k] / outDegree.get(i) * nsw / (1+nsw);	// Gaussian prior
	else
	  grad[0][k] -= thetaReg * theta[i][k] / outDegree.get(i) * 1.0 / (1+nsw);	// Gaussian prior
      } else {
	// dR
	if (!label) continue;
	double g = (yij-pij) * paramA * (1-1/(1+dTheta));
	grad[0][k] += g * lw * useNsw / outDegree.get(i);
	grad[1][k] -= g * lw * useNsw / outDegree.get(i);

	grad[0][k] += (R[i] != 0) ? (rAlpha / R[i]) / outDegree.get(i) : 100;	    // power law prior
      }
    }
  }

  public static double l2dist(double[] vec1, double[] vec2) {
    double res = 0;
    int len = vec1.length;
    for (int k = 0; k < len; k++) res += (vec1[k]-vec2[k]) * (vec1[k]-vec2[k]);
    return res;
  }

   public static double drCalcObj(List<Double> sim, double[] t, List<Integer> allEdges) {
    double res = 0;
    for (int _e = 0; _e < 2*E*tp; _e++) {
      int e = allEdges.get(_e);
      int i = edgeSources.get(e), j = edgeTargets.get(e);
      if (e < E) {
	double w = sim.get(e);
	res += w * Math.log(logis(CosSineTable.getCos(t[i]-t[j])));	  // likelihood
      } else {
	double w = sim.get(e);
	res += w * Math.log(1-logis(CosSineTable.getCos(t[i]-t[j])));
      }
    }
    res /= (2*E*tp);
    return res;
  }

  public static void drUpdate(double w, double[] t, int i, int j, boolean label, double[] grad) {
    int y = label ? 1 : 0;
    double s = logis(CosSineTable.getCos(t[i] - t[j]));
    double g = -w * (y-s) * CosSineTable.getSine(t[i] - t[j]);    // weighted
    grad[0] += g; grad[1] -= g;
  }

  public static void dimReduce(double[][] theta, double[] t) {
    List<Double> sim = new ArrayList<Double>();
    for (int e = 0; e < 2*E; e++) {
      int i = edgeSources.get(e), j = edgeTargets.get(e);
      double w = l2dist(theta[i], theta[j]);
      if (w < 0.5) 
	sim.add(1.0);
      else
	sim.add(0.0);
    }

    List<Integer> allEdges = new ArrayList<Integer>(2*E);     // index of all edges
    for (int i = 0; i < 2*E; i++) allEdges.add(i);
    Collections.shuffle(allEdges, new Random(shuffleSeed));
    double oldRes = drCalcObj(sim, t, allEdges), newRes = 0.0;

    long numEdges = 0;
    for (int _e = 0; _e < 2*E*tp; _e++) {
      int e = allEdges.get(_e);
      numEdges += weights.get(e);
    }
    int[][] edgeTable = new int[4][1<<30];
    long part = 0; int cur = 0;
    for (long i = 0; i < numEdges; i++) {
      if (i+1 > part) {
	part += weights.get(allEdges.get(cur));
	cur++;
      }
      int row = (int) (i >>> 30);
      int col = (int) (i & ((1 << 30) -1));
      edgeTable[row][col] = allEdges.get(cur-1);
    }

    for (int count = 0; count < 2*MAX_EDGES; count++) {
      if (count%(MAX_EDGES/10) == 0 && verbose) {
	System.out.printf("\n%d\n", count);
	newRes = drCalcObj(sim, t, allEdges);
	System.out.printf("[Train] obj = %f\n", newRes);
	test_performance(t);
      }
      long randl = nextLong(rand, numEdges);
      int row = (int) (randl >>> 30);
      int col = (int) (randl & ((1 << 30) - 1));
      int e = edgeTable[row][col];
      int i = edgeSources.get(e), j = edgeTargets.get(e);
      double w = sim.get(e);
      if (w == 0) continue;

      double[] grad = new double[2];
      if (e < E) 
	drUpdate(w, t, i, j, true, grad);
      else 
	drUpdate(0.5, t, i, j, false, grad);

      double thisStepSize = sgdStepSize * (1.0 - count / MAX_EDGES);
      t[i] += thisStepSize * grad[0];
      t[j] += thisStepSize * grad[1];
    }
  }

  /**
   * stochastic gradient descent
   */
  public static void runSGD(double[] R, double[][] theta, double convergenceTol) {

    List<Integer> allEdges = new ArrayList<Integer>(2*E);     // index of all edges
    for (int i = 0; i < 2*E; i++) allEdges.add(i);
    Collections.shuffle(allEdges, new Random(shuffleSeed));
    double oldRes = calcObj(R, theta, allEdges), newRes = 0.0;
    //System.out.printf("[Train] obj (init) = %f\n", oldRes);

    long numEdges = 0;
    for (int _e = 0; _e < 2*E*tp; _e++) {
      int e = allEdges.get(_e);
      numEdges += weights.get(e);
    }
    if (verbose) System.out.printf("[Info] Number of edges in training, including multiplicity = %d\n", numEdges);
    int[][] edgeTable = new int[4][1<<30];
    long part = 0; int cur = 0;
    for (long i = 0; i < numEdges; i++) {
      if (i+1 > part) {
	part += weights.get(allEdges.get(cur));
	cur++;
      }
      int row = (int) (i >>> 30);
      int col = (int) (i & ((1 << 30) -1));
      edgeTable[row][col] = allEdges.get(cur-1);
    }

    for (int count = 0; count < MAX_EDGES; count++) {
	if (count%(MAX_EDGES/10) == 0 && verbose) {
	  System.out.printf("\n%d\n", count);
	  newRes = calcObj(R, theta, allEdges);
	  System.out.printf("[Train] obj = %f\n", newRes);
	  test_performance(R, theta);
	  if (count/(MAX_EDGES/10) >= 3 && -(newRes-oldRes)/oldRes < convergenceTol) break;
	  oldRes = newRes;
	}

	long randl = nextLong(rand, numEdges);
	int row = (int) (randl >>> 30);
	int col = (int) (randl & ((1 << 30) - 1));
	int e = edgeTable[row][col];
	int s = edgeSources.get(e), t = edgeTargets.get(e);   // s -> t

	double[][] grad = new double[2][K+1];
	if (e < E) {
	  update(R, theta, s, t, true, grad);
	} else {
	  update(R, theta, s, t, false, grad);
	}

	// stepsize
	double thisStepSize = sgdStepSize * (1.0 - count / MAX_EDGES);
	for (int k = 0; k < K+1; k++) {
	  if (k < K) {
	    theta[s][k] += thisStepSize * grad[0][k];
	    theta[t][k] += thisStepSize * grad[1][k];
	  } else {
	    R[s] += thisStepSize * grad[0][k];
	    R[t] += thisStepSize * grad[1][k];
	  }
	}

	// projected gd
	if (R[s] < 0) R[s] = 0;
	if (R[t] < 0) R[t] = 0;
	if (R[s] > rLimit) R[s] = rLimit;
	if (R[t] > rLimit) R[t] = rLimit;
    }
  }

  // 1d theta only
  public static void test_performance(double[] theta) {
    List<Integer> allEdges = new ArrayList<Integer>(2*E);	// index of all edges
    for (int i = 0; i < 2*E; i++) allEdges.add(i);
    Collections.shuffle(allEdges, new Random(shuffleSeed));

    double auc = 0.0, posCount = 0, negCount = 0;
    List<Tuple> res = new ArrayList<Tuple>();
    for (int _e = (int)(2*E*tp); _e < 2*E; _e++) {
      int e = allEdges.get(_e);
      int i = edgeSources.get(e), j = edgeTargets.get(e);	// s -> t
      int w = weights.get(e);
      double sij = CosSineTable.getCos(theta[i]-theta[j]);
      if (e < E) {
	res.add(new Tuple(sij, true, w));
	posCount += w;
      } else {
	res.add(new Tuple(sij, false, w));
	negCount += w;
      }
    }

    Collections.sort(res, new Comparator<Tuple>() {
      @Override
      public int compare(Tuple t1, Tuple t2) {
	if (t1.prob < t2.prob) return 1;
	else if (t1.prob > t2.prob) return -1;
	else return 0;
      }
    });

    double nx = 0.0, ox = 0.0, ny = 0.0;
    for (Tuple t: res) {
      if (t.label) {
	ny += 1.0 * t.weight / posCount;
      } else {
	nx += 1.0 * t.weight / negCount;
      }
      auc += ny * (nx-ox);
      ox = nx;
    }
    System.out.printf("[Test] AUC = %f, nx = %f, ny = %f\n", auc, nx, ny);
  }


  // 1-dim R and K-dim theta (aka z)
  public static void test_performance(double[] R, double[][] theta) {
    List<Integer> allEdges = new ArrayList<Integer>(2*E);	// index of all edges
    for (int i = 0; i < 2*E; i++) allEdges.add(i);
    Collections.shuffle(allEdges, new Random(shuffleSeed));

    double auc = 0.0, posCount = 0, negCount = 0;
    List<Tuple> res = new ArrayList<Tuple>();
    for (int _e = (int)(2*E*tp); _e < 2*E; _e++) {
      int e = allEdges.get(_e);
      int s = edgeSources.get(e), t = edgeTargets.get(e);	// s -> t
      int w = weights.get(e);
      double dTheta = l2dist(theta[s], theta[t]);
      double dij = paramA * (R[s]-R[t]) * (1-1/(1+dTheta)) - paramB * dTheta + paramC;
      if (e < E) {
	res.add(new Tuple(dij, true, w));
	posCount += w;
      } else {
	res.add(new Tuple(dij, false, w));
	negCount += w;
      }
    }

    Collections.sort(res, new Comparator<Tuple>() {
      @Override
      public int compare(Tuple t1, Tuple t2) {
	if (t1.prob < t2.prob) return 1;
	else if (t1.prob > t2.prob) return -1;
	else return 0;
      }
    });

    double nx = 0.0, ox = 0.0, ny = 0.0;
    int lg = 0;
    for (Tuple t: res) {
      lg++;
      if (t.label) {
	ny += 1.0 * t.weight / posCount;
      } else {
	nx += 1.0 * t.weight / negCount;
      }
      auc += ny * (nx-ox);
      ox = nx;
    }
    System.out.printf("[Test] AUC = %f, nx = %f, ny = %f\n", auc, nx, ny);
  }


  public static void start(String[] args) {
    init();
    if (verbose) System.out.println("[Info] Init done.\n");

    double[] R = new double[N];
    double[][] theta = new double[N][K];
    for (int n = 0; n < N; n++) {
      R[n] = 1;
      for (int k = 0; k < K; k++) {
	theta[n][k] = 1.0 * (rand.nextDouble() - 0.5);
      }
    }

    runSGD(R, theta, tolerance);
    if (verbose) System.out.println("");

    try {
      FileParser.output_2d_1(R, outFilenameR, invMap, false);
      FileParser.output_2d_2(theta, outFilenameTheta, invMap);
    } catch (IOException e) {
      System.out.println("[I/O] File output error.");
    }

    /* 
     * the above scripts are enough to learn the embedding (r and z) in an
     * information network
     *
     * the scripts below are used for visualization (polarized coordinates)
     * uncomment if you want to see the polarized representation of data points
     */

    /*
    double[] d1t = new double[N];   // 1 dim theta
    for (int n = 0; n < N; n++) d1t[n] = rand.nextDouble() * 2 * Math.PI;
    dimReduce(theta, d1t);

    if (!verbose) test_performance(R, theta);
    if (!verbose) test_performance(d1t);

    try {
      FileParser.output_2d_1(d1t, outFilename1DTheta, invMap, true);
    } catch (IOException e) {
      System.out.println("[I/O] File output error.");
    }
    */
  }


  public static void main(String[] args) {
  /*
    switch (args.length) {
      case 1:	    // java Main <config_file>
	FileParser.readConfig(args[0]);
	break;
      case 4:	    // java Main <K> <dirR> <dirZ> <dirTheta>
	K = Integer.parseInt(args[0]);
	outFilenameR = args[1];
	outFilenameTheta = args[2];
	outFilename1DTheta = args[3];
	System.out.printf("[Info] Output to %s and %s\n", outFilenameR, outFilenameTheta);
	break;
      default:
	outFilenameR = "./res_r.txt";
	outFilenameTheta = "./res_z.txt";
	outFilename1DTheta = "./res_theta.txt";
	break;
    }
    */

    try {
      ArgumentParser.parse(args);
    } catch (NumberFormatException e) {
      //e.printStackTrace();
      System.out.println("\nIllegal arguments.");
      ArgumentParser.help();
      System.exit(0);
    }

    //System.out.printf("Parameters: %s\n%s\n%s\n", outFilenameR, outFilenameTheta, outFilename1DTheta);
    //System.out.printf("%f, %f, %f\n%f\n%d\n", paramA, paramB, paramC, thetaReg, K);

    long _start = System.currentTimeMillis();
    start(args);
    long _finish = System.currentTimeMillis();
    System.out.printf("Total time: %d seconds\n", (_finish-_start)/1000);
  }

}
