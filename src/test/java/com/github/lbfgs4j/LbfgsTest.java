package com.github.lbfgs4j;

import com.github.lbfgs4j.liblbfgs.Function;
import org.junit.Test;

import static org.junit.Assert.*;
import static java.lang.Math.pow;

public class LbfgsTest {

  /**
   * Example of using L-BFGS to find global minimum of Rosenbrock function
   * (http://en.wikipedia.org/wiki/Rosenbrock_function)
   *
   * f(x,y) = (1-x)^2 + 100*(y-x^2)^2
   * f(x,y) is non-convex and has global minimum at (x,y) = (1,1) where f(x,y) = 0
   *
   * f_x = -2*(1-x) - 400*(y-x^2)*x
   * f_y = 200*(y-x^2)
   */
  @Test
  public void rosenbrock() {
    
    Function f = new Function() {
      public int getDimension() {
        return 2;
      }
      public double valueAt(double[] x) {
        return pow(1-x[0], 2) + 100 * pow(x[1] - pow(x[0], 2), 2); 
      }
      public double[] gradientAt(double[] x) {
        double[] g = new double[2];
        g[0] = -2*(1-x[0]) - 400 * (x[1] - pow(x[0], 2)) * x[0];
        g[1] = 200 * (x[1] - pow(x[0], 2));
        return g;
      }
    };
   
    boolean verbose = true;
    LbfgsMinimizer minimizer = new LbfgsMinimizer(verbose);
    double[] x = minimizer.minimize(f);
    double min = f.valueAt(x);

    System.out.printf("The function achieves its minimum value = %.5f at: ", min);
    printOut(x);
    
    assertEquals(min, 0.0, 1e-3);
    assertArrayEquals(x, new double[] {1.0,  1.0}, 1e-3);
  }

  @Test
  public void simple() {

    // f(x) = (x-5)^2 + 1
    Function f = new Function() {
      public int getDimension() {
        return 1;
      }
      public double valueAt(double[] x) {
        return Math.pow(x[0]-5, 2) + 1;
      }
      public double[] gradientAt(double[] x) {
        return new double[] { 2*(x[0]-5) };
      }
    };

    boolean verbose = false;
    LbfgsMinimizer minimizer = new LbfgsMinimizer(verbose);
    double[] x = minimizer.minimize(f); // x should be {5}
    double min = f.valueAt(x);          // min should be 1

    System.out.printf("The function achieves its minimum value = %.5f at: ", min);
    printOut(x);

    assertEquals(min, 1.0, 1e-3);
    assertArrayEquals(x, new double[] {5.0}, 1e-3);
  }
  
  // Helper function
  public void printOut(double[] x) {
    System.out.printf("[");
    for (double v: x)
      System.out.printf(" %.2f", v);
    System.out.printf(" ]\n");	  
  }
}
