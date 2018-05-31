import java.io.*; 
import java.util.*; 
import java.lang.*;

public class FileParser { 
  public static int countLines(String filename) throws IOException {
    InputStream is = new BufferedInputStream(new FileInputStream(filename)); 
    try {
      byte[] c = new byte[1024]; 
      boolean empty = true; 
      int count = 0; 
      int readChars = 0; 
      while ((readChars = is.read(c)) != -1) { 
	empty = false; 
	for (int i = 0; i < readChars; ++i) { 
	  if (c[i] == '\n') {
	    count++;
	  }
	}
      }
      return (count == 0 && !empty) ? 1 : count;
    } 
    finally {is.close();}
  }

  public static void readConfig(String configFile) {
    try {
      BufferedReader br = new BufferedReader(new FileReader(configFile));
      String currentLine;
      while ((currentLine = br.readLine()) != null) {
	currentLine = currentLine.trim();
	if (currentLine.equals("")) continue;		// empty line
	if (currentLine.charAt(0) == '#') continue;	// comment
	String[] tokens = currentLine.split("=");
	String x = tokens[0].trim(), y = tokens[1].trim();
	switch (x) {
	  case "dirR":
	    Main.outFilenameR = y;
	    break;
	  case "dirTheta":
	    Main.outFilenameTheta = y;
	    break;
	  case "dir1DTheta":
	    Main.outFilename1DTheta = y;
	    break;
	  case "A":
	    double _a = Double.parseDouble(y);
	    if (_a <= 0) throw new NumberFormatException();
	    Main.paramA = _a;
	    break;
	  case "B":
	    double _b = Double.parseDouble(y);
	    if (_b <= 0) throw new NumberFormatException();
	    Main.paramB = _b;
	    break;
	  case "C":
	    double _c = Double.parseDouble(y);
	    if (_c >= 0) throw new NumberFormatException();
	    Main.paramC = _c;
	    break;
	  case "K":
	    int _k = Integer.parseInt(y);
	    if (_k <= 0) throw new NumberFormatException();
	    Main.K = _k;
	    break;
	  /*
	  case "gammaK":
	    double _gk = Double.parseDouble(y);
	    if (_gk < 0) throw new NumberFormatException();
	    Main.gammaK = _gk;
	    break;
	  case "gammaTheta":
	    double _gt = Double.parseDouble(y);
	    if (_gt < 0) throw new NumberFormatException();
	    Main.gammaTheta = _gt;
	    break;
	  */
	  case "thetaReg":
	    double _tr = Double.parseDouble(y);
	    if (_tr < 0 || _tr >= 1) throw new NumberFormatException();
	    Main.thetaReg = _tr;
	    break;
	}
      }
    } catch (IOException e) {
      System.out.println("I/O Error");
      e.printStackTrace();
      System.exit(0);
    } catch (NumberFormatException e) {
      System.out.println("Parameter Error");
      e.printStackTrace();
      System.exit(0);
    }
  }


  /**
   * readCSVDict: 
   *	build a map (original ID: new ID in memory) and an inverted map
   */
  public static int readCSVDict(String filename1, String filename2, Map<String, Integer> map, Map<Integer, String> invMap) 
  throws FileNotFoundException, IOException {
    BufferedReader br = new BufferedReader(new FileReader(filename1));
    try {
      int newID = 0;
      String currentLine;
      int lineCount = 0;
      while ((currentLine = br.readLine()) != null) {
	if (Main.verbose && lineCount % 1000000 == 0)	// print every 1M records
	  System.out.println(lineCount);
	lineCount += 1;

	// each line: sampleID, OTU_ID, weigth
	String[] tokens = currentLine.split("\t");
	String x = tokens[0];
	String y = tokens[1];
	if (x == "" || y == "") {System.out.println("hi"); continue;}
	if (!map.containsKey(x)) {
	  map.put(x, newID); newID += 1;
	}
	if (!map.containsKey(y)) {
	  map.put(y, newID); newID += 1;
	}
      }

      // set inverse map
      for (Map.Entry<String, Integer> e: map.entrySet()) {
	invMap.put(e.getValue(), e.getKey());
      }
      return map.size();
    } finally { br.close(); }
  }


  /*
   * readCSVGraphApprox: read csv graph into a list of (x,y)'s 
   * WITH a dictionary in memory; weight indicated by WEIGHTED.
   *
   * Negative samples are approximate, in other words, non-existing links 
   * may actually overlap with existing links.
   *
   */
  public static int readCSVGraphApprox(
      String filename, Map<String, Integer> map, 
      List<Integer> edgeSources, List<Integer> edgeTargets, List<Integer> weights, List<Boolean> isRecip,
      int N, int[] negDict, Map<Integer, Integer> outDegree, boolean WEIGHTED
  ) throws FileNotFoundException, IOException {
    Map<Integer, Integer> inDegree = new HashMap<Integer, Integer>();	    // in-degree of every node
    for (int i = 0; i < N; i++) outDegree.put(i, 0);
    for (int i = 0; i < N; i++) inDegree.put(i, 0);
    Random rand = new Random(0);

    BufferedReader br = new BufferedReader(new FileReader(filename));
    int count = 0;
    try {
      String currentLine = "";
      int lineCount = 0;
      while ((currentLine = br.readLine()) != null) {
	// print every 1M records
	if (Main.verbose && lineCount % 1000000 == 0) System.out.println(lineCount);
	lineCount += 1;

	String[] tokens = currentLine.split("\t");
	if (tokens[0] == "" || tokens[1] == "") continue;
	int x = map.get(tokens[0]), y = map.get(tokens[1]);
	int w = 1;
	if (WEIGHTED) w = Integer.parseInt(tokens[2]);
	edgeSources.add(x); edgeTargets.add(y); weights.add(w);

	// update out-degree, in-degree and #edges
	outDegree.put(x, outDegree.get(x)+w);
	inDegree.put(y, inDegree.get(y)+w);
	count++;
      }
    } finally { br.close(); }

    // build negative samples dictionary
    if (Main.verbose) System.out.println("[Info] Building negative sample dictionary");
    int negDictSize = negDict.length;
    double cur_sum = 0, tot_sum = 0, part = 0;
    int n = 0;
    double power = 0.75;
    for (int i = 0; i < N; i++) 
      //tot_sum += Math.pow(inDegree.get(i), power);
      tot_sum += Math.pow(outDegree.get(i), power);
      //tot_sum += 1;
    for (int i = 0; i < negDictSize; i++) {
      if (1.0 * (i+1) / negDictSize > part) {
	//cur_sum += Math.pow(inDegree.get(n), power);
	cur_sum += Math.pow(outDegree.get(n), power);
	//cur_sum += 1;
	part = cur_sum / tot_sum;
	n++;
      }
      negDict[i] = n-1;
    }
    if (N != n) System.out.printf("N = %d, n = %d\n", N, n);	// should be equal to N

    // node with no out links (out_degree = 0)
    /*
    int nol = 0;
    for (int i = 0; i < N; i++) if (outDegree.get(i) == 0) {
      int repeat = 1;
      for (int _z = 0; _z < repeat; _z++) {
	int j = negDict[rand.nextInt(negDictSize)];
	//int j = rand.nextInt(N);
	edgeSources.add(i); edgeTargets.add(j); weights.add(1);
	count++;
      }
      nol++;
      outDegree.put(i, repeat);
    }
    System.out.println("nol = " + nol);
    */

    // build negative links
    for (int e = 0; e < count; e++) {
      int i = edgeSources.get(e);
      if (outDegree.get(i) == 0) {
	e--; continue;
      }
      int w = weights.get(e);
      int j = negDict[rand.nextInt(negDictSize)];
      edgeSources.add(i);
      edgeTargets.add(j);
      weights.add(w);
    }

    return count;
  }



  /**
   * readCSVGraph: read csv graph into a list of (x,y)'s 
   * WITH a dictionary in memory; weight indicated by WEIGHTED
   */
  public static int readCSVGraph(
      String filename, Map<String, Integer> map, 
      List<Integer> edgeSources, List<Integer> edgeTargets, List<Integer> weights, List<Boolean> isRecip,
      int N, int[] negDict, Map<Integer, Integer> outDegree, boolean WEIGHTED
  ) throws FileNotFoundException, IOException {
    Map<Integer, Integer> inDegree = new HashMap<Integer, Integer>();	    // in-degree of every node
    for (int i = 0; i < N; i++) outDegree.put(i, 0);
    for (int i = 0; i < N; i++) inDegree.put(i, 0);
    Map<Long, List<Integer>> pairMap = new HashMap<Long, List<Integer>>();  // (x,y) -> indexes
    Set<Long> secondPairMap = new HashSet<Long>();
    Random rand = new Random(0);

    BufferedReader br = new BufferedReader(new FileReader(filename));
    int count = 0;
    try {
      String currentLine = "";
      int lineCount = 0;
      while ((currentLine = br.readLine()) != null) {
	// print every 1M records
	if (Main.verbose && lineCount % 1000000 == 0) System.out.println(lineCount);
	lineCount += 1;

	String[] tokens = currentLine.split("\t");
	if (tokens[0] == "" || tokens[1] == "") continue;
	if (tokens[0] == tokens[1]) continue;
	int x = map.get(tokens[0]), y = map.get(tokens[1]);
	int w = 1;
	if (WEIGHTED) w = Integer.parseInt(tokens[2]);
	edgeSources.add(x); edgeTargets.add(y); weights.add(w);

	// update out-degree and in-degree
	outDegree.put(x, outDegree.get(x)+w);
	inDegree.put(y, inDegree.get(y)+w);

	/*
	// detect repcpricated links: delete for efficiency?
	long key = x * N + y;
	if (pairMap.containsKey(key)) {
	  List<Integer> tmp = pairMap.get(key);
	  tmp.add(count); pairMap.put(key, tmp);
	} else {
	  List<Integer> tmp = new ArrayList<Integer>();
	  tmp.add(count); pairMap.put(key, tmp);
	}

	long revKey = y * N + x;
	if (secondPairMap.contains(revKey)) {
	  isRecip.add(true);
	  count++;
	  continue;
	}

	if (pairMap.containsKey(revKey)) {
	  secondPairMap.add(key);
	  secondPairMap.add(revKey);

	  isRecip.add(true);
	  List<Integer> revIndexes = pairMap.get(revKey);
	  for (int revIndex: revIndexes) 
	    isRecip.set(revIndex, true);
	} else {
	  isRecip.add(false);
	}
	*/

	count++;
      }
    } finally { br.close(); }

    /*
    // write reciprocated links to file
    String recFileName = "./recip.txt";
    if (Main.verbose) System.out.println("[Info] Output recriprocated info to " + recFileName);
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(recFileName)))) {
      for (boolean b: isRecip) {
	writer.printf("%b\n", b);
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
    */

    // build negative samples dictionary
    if (Main.verbose) System.out.println("[Info] Building negative sample dictionary");
    int negDictSize = negDict.length;
    double cur_sum = 0, tot_sum = 0, part = 0;
    int n = 0;
    double power = 0.75;
    for (int i = 0; i < N; i++) 
      //tot_sum += Math.pow(inDegree.get(i), power);
      tot_sum += Math.pow(outDegree.get(i), power);
      //tot_sum += 1;
    for (int i = 0; i < negDictSize; i++) {
      if (1.0 * (i+1) / negDictSize > part) {
	//cur_sum += Math.pow(inDegree.get(n), power);
	cur_sum += Math.pow(outDegree.get(n), power);
	//cur_sum += 1;
	part = cur_sum / tot_sum;
	n++;
      }
      negDict[i] = n-1;
    }
    if (N != n) System.out.printf("N = %d, n = %d\n", N, n);	// should be equal to N

    // node with no out links (out_degree = 0)
    ///*
    int nol = 0;
    for (int i = 0; i < N; i++) if (outDegree.get(i) == 0) {
      int repeat = 1;
      for (int _z = 0; _z < repeat; _z++) {
	int j = negDict[rand.nextInt(negDictSize)];
	//int j = rand.nextInt(N);
	edgeSources.add(i);
	edgeTargets.add(j);
	weights.add(1);
	count++;
      }
      nol++;
      outDegree.put(i, repeat);
    }
    //System.out.println("nol = " + nol);
    //*/

    // build negative links
    for (int e = 0; e < count; e++) {
      //int i = rand.nextInt(N);
      int i = edgeSources.get(e);
      int w = weights.get(e);
      while (true) {
	int j = negDict[rand.nextInt(negDictSize)];
	long xyid = i * N + j;
	if (!pairMap.containsKey(xyid)) {
	  edgeSources.add(i); edgeTargets.add(j); weights.add(w);
	  break;
	}
      }
    }

    return count;
  }



  /** 
   * output_2d
   *	output to file
   *	arr_s.get(t) is a n*1 array 
   */
  public static void output_2d_1(double[] arr, String filename, Map<Integer, String> invMap, boolean constraint) 
  throws IOException {
    PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
    try {
      for (int i = 0; i < arr.length; i++) {
	double v = arr[i];
	if (constraint) {
	  while (v < 0) v += 2 * Math.PI;
	  while (v > 2 * Math.PI) v -= 2 * Math.PI;
	}

	String newline = invMap.get(i) + "\t" + Double.toString(v);
	writer.printf("%s\n", newline);
      }
    } finally {
      writer.close();
    }
    return;
  }


  public static void output_2d_2(double[][] arr, String filename, Map<Integer, String> invMap) 
  throws IOException {
    PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
    try {
      for (int i = 0; i < arr.length; i++) {
	String newline = invMap.get(i);
	for (int k = 0; k < Main.K; k++) {
	  newline = newline + "\t" + Double.toString(arr[i][k]);
	}
	writer.printf("%s\n", newline);
      }
    } finally {
      writer.close();
    }
    return;
  }

  /* arr_s.get(t) is a n*1 array */
  public static void
  output_1d(List<double[]> arr_s, String fileDir) {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
      for (int i = 0; i < arr_s.size(); i++) {
	double[] arr = arr_s.get(i);
        writer.printf("%d ", i);
	for (int j = 0; j < arr.length; j++) {
	  writer.printf("%f ", arr[j]);
	}
        writer.printf("\n");
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /* arr_s.get(t) is a n*1 array */
  public static void
  output(Map<Integer, Integer> id_map, String fileDir) {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
      for (Map.Entry<Integer, Integer> e: id_map.entrySet()) {
	int globalID = e.getKey();
	int localID = e.getValue();
	writer.printf("%d %d\n", globalID, localID);
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

}
