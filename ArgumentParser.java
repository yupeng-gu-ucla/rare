import java.io.*;
import java.util.*;

/* options:
 *   -k: dimension of embedding; 2 by default
 *   -data: the location of data file
 *   -weighted: whether the network is weighted (0: unweighted; other: weighted); 0 by default
 *   -out: the directory of output files; "." (current directory) by default
 *   -lr: initial stepsize; 0.05 by default
 *   -lambda_r: coefficient of dr (first component), must be positive; 1.0 by default
 *   -lambda_z: coefficient of dz (second component), must be positive; 1.0 by default
 *   -lambda_0: bias term (third component), should be negative for sparse networks; -1.0 by default
 *   -alpha_r: parameter for power law prior (on r), recommended value: 1.0~2.0; 1.5 by default
 *   -reg_z: l2 regularization coefficient (on z); 1e-6 by default
 *   -tp: fraction of nodes for training, must be between 0 and 1; 0.9 by default
 *   -tolerance: stopping criterion on log-likelihood improvement (negative number for no tolerance); 1e-6 by default
 *   -nsw: negative sample weight, must be positive; 5 by default
 *   -iter: value of maximum edges to be sampled (in thousands); 10000 thousand by default
 *   -verbose: whether see verbose output or not (0: show limited outputs; other: show all outputs); 1 by default
 */

public class ArgumentParser {
  public static void parse(String[] args) {
    int k = 0, K = args.length;
    while (k < K) {
      switch (args[k]) {
	case "-k":
	  Main.K = Integer.parseInt(args[++k]); break;
	case "-data":
	  Main.posFilename = args[++k]; break;
	case "-weighted":
	  int w = Integer.parseInt(args[++k]);
	  Main.WEIGHTED = (w != 0); break;
	case "-out":
	  Main.outPrefix = args[++k]; break;
	case "-lr":
	  Main.sgdStepSize = Double.parseDouble(args[++k]); break;
	case "-lambda_r":
	  Main.paramA = Double.parseDouble(args[++k]); break;
	case "-lambda_z":
	  Main.paramB = Double.parseDouble(args[++k]); break;
	case "-lambda_0":
	  Main.paramC = Double.parseDouble(args[++k]); break;
	case "-alpha_r":
	  Main.rAlpha = Double.parseDouble(args[++k]); break;
	case "-reg_z":
	  Main.thetaReg = Double.parseDouble(args[++k]); break;
	case "-tp":
	  Main.tp = Double.parseDouble(args[++k]); break;
	case "-tolerance":
	  Main.tolerance = Double.parseDouble(args[++k]); break;
	case "-nsw":
	  Main.nsw = Double.parseDouble(args[++k]); break;
	case "-iter":
	  Main.MAX_EDGES = 1000 * Integer.parseInt(args[++k]); break;
	case "-verbose":
	  int v = Integer.parseInt(args[++k]);
	  Main.verbose = (v != 0); break;
	case "-help":
	case "-h":
	  help();
	  System.exit(0);
      }
      k++;
    }
  }

  public static void help() {
    String h = "\nUsage: java Main <options>\n" 
      + "options:\n" 
      + "\t-k: dimension of embedding; 2 by default\n"
      + "\t-data: the location of data file\n"
      + "\t-weighted: whether the network is weighted (0: unweighted; other: weighted); 0 by default\n"
      + "\t-out: the directory of output files; \".\" (current directory) by default\n"
      + "\t-lr: initial stepsize; 0.05 by default\n"
      + "\t-lambda_r: coefficient of dr (first component), must be positive; 1.0 by default\n"
      + "\t-lambda_z: coefficient of dz (second component), must be positive; 1.0 by default\n"
      + "\t-lambda_0: bias term (third component), should be negative for sparse networks; -1.0 by default\n"
      + "\t-alpha_r: parameter for power law prior (on r), recommended value: 1.0~2.0; 1.5 by default\n"
      + "\t-reg_z: l2 regularization coefficient (on z); 1e-6 by default\n"
      + "\t-tp: fraction of nodes for training, must be between 0 and 1; 0.9 by default\n"
      + "\t-tolerance: stopping criterion on log-likelihood improvement (negative number for no tolerance); 1e-6 by default\n"
      + "\t-nsw: negative sample weight, must be positive; 5 by default\n"
      + "\t-iter: value of maximum edges to be sampled (in thousands); 10000 thousand by default\n"
      + "\t-verbose: whether see verbose output or not (0: show limited outputs; other: show all outputs); 0 by default\n"
      + "\t-h or -help: show this help\n";
    System.out.println(h);
  }
}

