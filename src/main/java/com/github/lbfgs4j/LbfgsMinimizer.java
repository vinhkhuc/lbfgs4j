package com.github.lbfgs4j;

import com.github.lbfgs4j.liblbfgs.Function;
import com.github.lbfgs4j.liblbfgs.Lbfgs;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Param;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Evaluate;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Progress;
import com.github.lbfgs4j.liblbfgs.MutableDouble;

import static com.github.lbfgs4j.liblbfgs.LbfgsConstant.ReturnValue.LBFGS_SUCCESS;
import static java.lang.System.out;

/**
 * Optimization solver using LibLBFGS
 */
public class LbfgsMinimizer {

  // Parameters for the L-BFGS optimization
  private LBFGS_Param param;
  private boolean verbose;

  public LbfgsMinimizer() {
    this(Lbfgs.defaultParams(), true);
  }

  public LbfgsMinimizer(LBFGS_Param param) {
    this(param, true);
  }

  public LbfgsMinimizer(boolean verbose) {
    this(Lbfgs.defaultParams(), verbose);
  }

  public LbfgsMinimizer(LBFGS_Param param, boolean verbose) {
    this.param   = param;
    this.verbose = verbose;
  }

  public double[] minimize(Function f) {
    int dim = f.getDimension();
    double[] x = new double[dim];

    MutableDouble fx  = new MutableDouble();
    Evaluate eval     = new Evaluate(f);
    Progress progress = verbose? new Progress() : null;

    // Start the L-BFGS optimization; this will invoke 
    // the callback functions evaluate() and progress() when necessary.
    Lbfgs.lbfgs(dim, x, fx, eval, progress, null, param);
    
    return x;
  }
  
  /**
   * Function evaluator
   */
  class Evaluate implements LBFGS_Evaluate {
    
    private Function func;
    
    public Evaluate(Function func) {
      this.func = func;
    }

    public double eval(
        Object   instance,
        double[] x,
        double[] g,
        int      n,
        double   step
        )
    {
      // Get function's gradient at x 
      double[] grad = func.gradientAt(x);
      System.arraycopy(grad, 0, g, 0, grad.length);
      
      // Return function's value at x 
      return func.valueAt(x);
    }
  }
  
  /**
   * Solving progress report
   */
  class Progress implements LBFGS_Progress {
    public LbfgsConstant.ReturnValue eval(
        Object   instance,
        double[] x,
        double[] g,
        double   fx,
        double   xnorm,
        double   gnorm,
        double   step,
        int      n,
        int      k,
        int      ls
        )
    {
      out.printf("Iteration %d:\n", k);
      out.printf("\tfx = %f, xnorm = %f, gnorm = %f, step = %f\n", fx, xnorm, gnorm, step);
      out.printf("\tn = %d, k = %d, ls = %d\n\n", n, k, ls);
      return LBFGS_SUCCESS;
    }
  }
}
