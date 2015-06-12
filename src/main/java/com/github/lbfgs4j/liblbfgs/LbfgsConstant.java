package com.github.lbfgs4j.liblbfgs;

/*
 * Copyright (c) 1990, Jorge Nocedal
 * Copyright (c) 2007-2010 Naoaki Okazaki
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

public class LbfgsConstant {

	/**
	 * Return value of lbfgs() 
	 */
	public static enum ReturnValue {
	    LBFGS_SUCCESS(0),
	    /** L-BFGS reaches convergence. */
	    LBFGS_CONVERGENCE(0),
	    LBFGS_STOP(1),
	    /** The initial variables already minimize the objective function. */
	    LBFGS_ALREADY_MINIMIZED(2),

	    /** Unknown error. */
	    LBFGSERR_UNKNOWNERROR(-1024),
	    /** Logic error. */
	    LBFGSERR_LOGICERROR(-1023),
	    /** Insufficient memory. */
	    LBFGSERR_OUTOFMEMORY(-1022),
	    /** The minimization process has been canceled. */
	    LBFGSERR_CANCELED(-1021),
	    /** Invalid number of variables specified. */
	    LBFGSERR_INVALID_N(-1020),
	    /** Invalid number of variables (for SSE) specified. */
	    LBFGSERR_INVALID_N_SSE(-1019),
	    /** The array x must be aligned to 16 (for SSE). */
	    LBFGSERR_INVALID_X_SSE(-1018),
	    /** Invalid parameter LBFGS_Param::epsilon specified. */
	    LBFGSERR_INVALID_EPSILON(-1017),
	    /** Invalid parameter LBFGS_Param::past specified. */
	    LBFGSERR_INVALID_TESTPERIOD(-1016),
	    /** Invalid parameter LBFGS_Param::delta specified. */
	    LBFGSERR_INVALID_DELTA(-1015),
	    /** Invalid parameter LBFGS_Param::linesearch specified. */
	    LBFGSERR_INVALID_LINESEARCH(-1014),
	    /** Invalid parameter LBFGS_Param::max_step specified. */
	    LBFGSERR_INVALID_MINSTEP(-1013),
	    /** Invalid parameter LBFGS_Param::max_step specified. */
	    LBFGSERR_INVALID_MAXSTEP(-1012),
	    /** Invalid parameter LBFGS_Param::ftol specified. */
	    LBFGSERR_INVALID_FTOL(-1011),
	    /** Invalid parameter LBFGS_Param::wolfe specified. */
	    LBFGSERR_INVALID_WOLFE(-1010),
	    /** Invalid parameter LBFGS_Param::gtol specified. */
	    LBFGSERR_INVALID_GTOL(-1009),
	    /** Invalid parameter LBFGS_Param::xtol specified. */
	    LBFGSERR_INVALID_XTOL(-1008),
	    /** Invalid parameter LBFGS_Param::max_linesearch specified. */
	    LBFGSERR_INVALID_MAXLINESEARCH(-1007),
	    /** Invalid parameter LBFGS_Param::orthantwise_c specified. */
	    LBFGSERR_INVALID_ORTHANTWISE(-1006),
	    /** Invalid parameter LBFGS_Param::orthantwise_start specified. */
	    LBFGSERR_INVALID_ORTHANTWISE_START(-1005),
	    /** Invalid parameter LBFGS_Param::orthantwise_end specified. */
	    LBFGSERR_INVALID_ORTHANTWISE_END(-1004),
	    /** The line-search step went out of the interval of uncertainty. */
	    LBFGSERR_OUTOFINTERVAL(-1003),
	    /** A logic error occurred; alternatively, the interval of uncertainty
	        became too small. */
	    LBFGSERR_INCORRECT_TMINMAX(-1002),
	    /** A rounding error occurred; alternatively, no line-search step
	        satisfies the sufficient decrease and curvature conditions. */
	    LBFGSERR_ROUNDING_ERROR(-1001),
	    /** The line-search step became smaller than LBFGS_Param::min_step. */
	    LBFGSERR_MINIMUMSTEP(-1000),
	    /** The line-search step became larger than LBFGS_Param::max_step. */
	    LBFGSERR_MAXIMUMSTEP(-999),
	    /** The line-search routine reaches the maximum number of evaluations. */
	    LBFGSERR_MAXIMUMLINESEARCH(-998),
	    /** The algorithm routine reaches the maximum number of iterations. */
	    LBFGSERR_MAXIMUMITERATION(-997),
	    /** Relative width of the interval of uncertainty is at most 
	        LBFGS_Param::xtol. */
	    LBFGSERR_WIDTHTOOSMALL(-996),
	    /** A logic error (negative line-search step) occurred. */
	    LBFGSERR_INVALIDPARAMETERS(-995),
	    /** The current search direction increases the objective function value. */
	    LBFGSERR_INCREASEGRADIENT(-994);
	    
	    public int val;
	    
	    private ReturnValue() {}
	    private ReturnValue(int val) { this.val = val; }
	}
	
	/**
	 * Line search algorithms.
	 */
	public static enum LineSearch {
		
	    /** The default algorithm (MoreThuente method). */
	    LBFGS_LINESEARCH_DEFAULT(0),
	    /** MoreThuente method proposed by More and Thuente. */
	    LBFGS_LINESEARCH_MORETHUENTE(0),
	    /**
	     * Backtracking method with the Armijo condition.
	     *  The backtracking method finds the step length such that it satisfies
	     *  the sufficient decrease (Armijo) condition,
	     *    - f(x + a * d) <= f(x) + LBFGS_Param::ftol * a * g(x)^T d,
	     *
	     *  where x is the current point, d is the current search direction, and
	     *  a is the step length.
	     */
	    LBFGS_LINESEARCH_BACKTRACKING_ARMIJO(1),
	    /** The backtracking method with the default (regular Wolfe) condition. */
	    LBFGS_LINESEARCH_BACKTRACKING(2),
	    /**
	     * Backtracking method with regular Wolfe condition.
	     *  The backtracking method finds the step length such that it satisfies
	     *  both the Armijo condition (LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
	     *  and the curvature condition,
	     *    - g(x + a * d)^T d >= LBFGS_Param::wolfe * g(x)^T d,
	     *
	     *  where x is the current point, d is the current search direction, and
	     *  a is the step length.
	     */
	    LBFGS_LINESEARCH_BACKTRACKING_WOLFE(2),
	    /**
	     * Backtracking method with strong Wolfe condition.
	     *  The backtracking method finds the step length such that it satisfies
	     *  both the Armijo condition (LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
	     *  and the following condition,
	     *    - |g(x + a * d)^T d| <= LBFGS_Param::wolfe * |g(x)^T d|,
	     *
	     *  where x is the current point, d is the current search direction, and
	     *  a is the step length.
	     */
	    LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE(3);
		
	    public int val;
	    
	    private LineSearch() {}
	    private LineSearch(int val) { this.val = val; }
	}
	
	/**
	 * L-BFGS optimization parameters.
	 *  Call lbfgs_parameter_init() function to initialize parameters to the
	 *  default values.
	 */
	public static class LBFGS_Param {

		/**
	     * The number of corrections to approximate the inverse Hessian matrix.
	     *  The L-BFGS routine stores the computation results of previous \ref m
	     *  iterations to approximate the inverse Hessian matrix of the current
	     *  iteration. This parameter controls the size of the limited memories
	     *  (corrections). The default value is \c 6. Values less than \c 3 are
	     *  not recommended. Large values will result in excessive computing time.
	     */
	    public int             m;

	    /**
	     * Epsilon for convergence test.
	     *  This parameter determines the accuracy with which the solution is to
	     *  be found. A minimization terminates when
	     *      ||g|| < \ref epsilon * max(1, ||x||),
	     *  where ||.|| denotes the Euclidean (L2) norm. The default value is
	     *  \c 1e-5.
	     */
	    public double          epsilon;

	    /**
	     * Distance for delta-based convergence test.
	     *  This parameter determines the distance, in iterations, to compute
	     *  the rate of decrease of the objective function. If the value of this
	     *  parameter is zero, the library does not perform the delta-based
	     *  convergence test. The default value is \c 0.
	     */
	    public int             past;

	    /**
	     * Delta for convergence test.
	     *  This parameter determines the minimum rate of decrease of the
	     *  objective function. The library stops iterations when the
	     *  following condition is met:
	     *      (f' - f) / f < \ref delta,
	     *  where f' is the objective value of \ref past iterations ago, and f is
	     *  the objective value of the current iteration.
	     *  The default value is \c 0.
	     */
	    public double          delta;

	    /**
	     * The maximum number of iterations.
	     *  The lbfgs() function terminates an optimization process with
	     *  ::LBFGSERR_MAXIMUMITERATION status code when the iteration count
	     *  exceeds this parameter. Setting this parameter to zero continues an
	     *  optimization process until a convergence or error. The default value
	     *  is \c 0.
	     */
	    public int             max_iterations;

	    /**
	     * The line search algorithm.
	     *  This parameter specifies a line search algorithm to be used by the
	     *  L-BFGS routine.
	     */
	    public LineSearch     lineSearch;

	    /**
	     * The maximum number of trials for the line search.
	     *  This parameter controls the number of function and gradients evaluations
	     *  per iteration for the line search routine. The default value is \c 20.
	     */
	    public int             max_linesearch;

	    /**
	     * The minimum step of the line search routine.
	     *  The default value is \c 1e-20. This value need not be modified unless
	     *  the exponents are too large for the machine being used, or unless the
	     *  problem is extremely badly scaled (in which case the exponents should
	     *  be increased).
	     */
	    public double          min_step;

	    /**
	     * The maximum step of the line search.
	     *  The default value is \c 1e+20. This value need not be modified unless
	     *  the exponents are too large for the machine being used, or unless the
	     *  problem is extremely badly scaled (in which case the exponents should
	     *  be increased).
	     */
	    public double          max_step;

	    /**
	     * A parameter to control the accuracy of the line search routine.
	     *  The default value is \c 1e-4. This parameter should be greater
	     *  than zero and smaller than \c 0.5.
	     */
	    public double          ftol;

	    /**
	     * A coefficient for the Wolfe condition.
	     *  This parameter is valid only when the backtracking line-search
	     *  algorithm is used with the Wolfe condition,
	     *  ::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE or
	     *  ::LBFGS_LINESEARCH_BACKTRACKING_WOLFE .
	     *  The default value is \c 0.9. This parameter should be greater
	     *  the \ref ftol parameter and smaller than \c 1.0.
	     */
	    public double          wolfe;

	    /**
	     * A parameter to control the accuracy of the line search routine.
	     *  The default value is \c 0.9. If the function and gradient
	     *  evaluations are inexpensive with respect to the cost of the
	     *  iteration (which is sometimes the case when solving very large
	     *  problems) it may be advantageous to set this parameter to a small
	     *  value. A typical small value is \c 0.1. This parameter should be
	     *  greater than the \ref ftol parameter (\c 1e-4) and smaller than
	     *  \c 1.0.
	     */
	    public double          gtol;

	    /**
	     * The machine precision for floating-point values.
	     *  This parameter must be a positive value set by a client program to
	     *  estimate the machine precision. The line search routine will terminate
	     *  with the status code (::LBFGSERR_ROUNDING_ERROR) if the relative width
	     *  of the interval of uncertainty is less than this parameter.
	     */
	    public double          xtol;

	    /**
	     * Coefficient for the L1 norm of variables.
	     *  This parameter should be set to zero for standard minimization
	     *  problems. Setting this parameter to a positive value activates
	     *  Orthant-Wise Limited-memory Quasi-Newton (OWL-QN) method, which
	     *  minimizes the objective function F(x) combined with the L1 norm |x|
	     *  of the variables, {F(x) + C |x|}. This parameter is the coefficient
	     *  for the |x|, i.e., C. As the L1 norm |x| is not differentiable at
	     *  zero, the library modifies function and gradient evaluations from
	     *  a client program suitably; a client program thus have only to return
	     *  the function value F(x) and gradients G(x) as usual. The default value
	     *  is zero.
	     */
	    public double          orthantwise_c;

	    /**
	     * Start index for computing L1 norm of the variables.
	     *  This parameter is valid only for OWL-QN method
	     *  (i.e., \ref orthantwise_c != 0). This parameter b (0 <= b < N)
	     *  specifies the index number from which the library computes the
	     *  L1 norm of the variables x,
	     *      |x| := |x_{b}| + |x_{b+1}| + ... + |x_{N}| .
	     *  In other words, variables x_1, ..., x_{b-1} are not used for
	     *  computing the L1 norm. Setting b (0 < b < N), one can protect
	     *  variables, x_1, ..., x_{b-1} (e.g., a bias term of logistic
	     *  regression) from being regularized. The default value is zero.
	     */
	    public int             orthantwise_start;

	    /**
	     * End index for computing L1 norm of the variables.
	     *  This parameter is valid only for OWL-QN method
	     *  (i.e., \ref orthantwise_c != 0). This parameter e (0 < e <= N)
	     *  specifies the index number at which the library stops computing the
	     *  L1 norm of the variables x,
	     */
	    public int             orthantwise_end;

	    /**
	     * Constructor
	     */
		public LBFGS_Param(int m, double epsilon, int past, double delta,
				int max_iterations, LineSearch lineSearch, int max_linesearch,
				double min_step, double max_step, double ftol, double wolfe,
				double gtol, double xtol, double orthantwise_c,
				int orthantwise_start, int orthantwise_end) 
		{
			this.m                 = m;
			this.epsilon           = epsilon;
			this.past              = past;
			this.delta             = delta;
			this.max_iterations    = max_iterations;
			this.lineSearch        = lineSearch;
			this.max_linesearch    = max_linesearch;
			this.min_step          = min_step;
			this.max_step          = max_step;
			this.ftol              = ftol;
			this.wolfe             = wolfe;
			this.gtol              = gtol;
			this.xtol              = xtol;
			this.orthantwise_c     = orthantwise_c;
			this.orthantwise_start = orthantwise_start;
			this.orthantwise_end   = orthantwise_end;
		}
		
		/**
		 * Copy constructor
		 */
		public LBFGS_Param(LBFGS_Param p) {
			this(
				p.m, p.epsilon, p.past, p.delta,
				p.max_iterations, p.lineSearch, p.max_linesearch,
				p.min_step, p.max_step, p.ftol, p.wolfe,
				p.gtol, p.xtol, p.orthantwise_c,
				p.orthantwise_start, p.orthantwise_end
			);
		}
	}
	
	/**
	 * Callback interface to provide objective function and gradient evaluations.
	 *
	 *  The lbfgs() function call this function to obtain the values of objective
	 *  function and its gradients when needed. A client program must implement
	 *  this function to evaluate the values of the objective function and its
	 *  gradients, given current values of variables.
	 *  
	 *  @param  instance    The user data sent for lbfgs() function by the client.
	 *  @param  x           The current values of variables.
	 *  @param  g           The gradient vector. The callback function must compute
	 *                      the gradient values for the current variables.
	 *  @param  n           The number of variables.
	 *  @param  step        The current step of the line search routine.
	 *  @return double      The value of the objective function for the current
	 *                          variables.
	 */
	public static interface LBFGS_Evaluate {
		public double eval(
                Object instance,
                double[] x,
                double[] g,
                int n,
                double step
        );
	}
	
	/**
	 * Callback interface to receive the progress of the optimization process.
	 *
	 *  The lbfgs() function call this function for each iteration. Implementing
	 *  this function, a client program can store or display the current progress
	 *  of the optimization process.
	 *
	 *  @param  instance    The user data sent for lbfgs() function by the client.
	 *  @param  x           The current values of variables.
	 *  @param  g           The current gradient values of variables.
	 *  @param  fx          The current value of the objective function.
	 *  @param  xnorm       The Euclidean norm of the variables.
	 *  @param  gnorm       The Euclidean norm of the gradients.
	 *  @param  step        The line-search step used for this iteration.
	 *  @param  n           The number of variables.
	 *  @param  k           The iteration count.
	 *  @param  ls          The number of evaluations called for this iteration.
	 *  @return int         Zero to continue the optimization process. Returning a
	 *                      non-zero value will cancel the optimization process.
	 */
	public static interface LBFGS_Progress {
		public ReturnValue eval(
                Object instance,
                double[] x,
                double[] g,
                double fx,
                double xnorm,
                double gnorm,
                double step,
                int n,
                int k,
                int ls
        );
	}
}
