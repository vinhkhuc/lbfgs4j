package com.github.lbfgs4j.liblbfgs;

/**
 * Function interface
 */
public interface Function {

  /**
   * Get function's dimension 
   */
  int getDimension();
  
  /**
   * Compute function's value 
   * @param x
   * @return function's value at x
   */
  double valueAt(double[] x);
  
  /**
   * Compute function's gradient 
   * @param x
   * @return function's gradient at x
   */
  double[] gradientAt(double[] x);
}
