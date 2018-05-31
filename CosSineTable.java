public class CosSineTable {
  private static int N = 36000;
  private static int M = N/360;
  private static double[] cos = new double[N+1];
  private static double[] sin = new double[N+1];
  private static CosSineTable table = new CosSineTable();

  private CosSineTable() {
    for (int i = 0; i <= N; i++) {
      cos[i] = Math.cos(Math.toRadians(1.0*i/M));
      sin[i] = Math.sin(Math.toRadians(1.0*i/M));
    }
  }

  public static double getSine(double angle) {
    double deg_angle = Math.toDegrees(angle);
    while (deg_angle < 0) deg_angle += 360;
    while (deg_angle > 360) deg_angle -= 360;
    int angleCircle = (int)(deg_angle*M);
    return sin[angleCircle];
  }

  public static double getCos(double angle) {
    double deg_angle = Math.toDegrees(angle);
    while (deg_angle < 0) deg_angle += 360;
    while (deg_angle > 360) deg_angle -= 360;
    int angleCircle = (int)(deg_angle*M);
    return cos[angleCircle];
  }

  public static CosSineTable getTable() {
    return table;
  }
}

