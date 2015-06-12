package com.github.lbfgs4j.liblbfgs;

import static java.lang.Math.*;
import static com.github.lbfgs4j.liblbfgs.Arithmetic.*;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.*;
import static com.github.lbfgs4j.liblbfgs.LbfgsConstant.LineSearch.*;
import static com.github.lbfgs4j.liblbfgs.LbfgsConstant.ReturnValue.*;

/*
 *      Limited memory BFGS (L-BFGS).
 *
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
public class Lbfgs {

	static class CallbackData {
		int n;
		Object instance;
		LBFGS_Evaluate proc_evaluate;
		LBFGS_Progress proc_progress;
	}
	
	static class IterationData {
		double   alpha;
		double[] s;      /* [n] */
		double[] y;      /* [n] */
		double   ys;     /* vecdot(y, s) */
	}
	
	static LBFGS_Param _defparam = new LBFGS_Param(
				6, 1e-5, 0, 1e-5,
			    0, LBFGS_LINESEARCH_DEFAULT, 40,
			    1e-20, 1e20, 1e-4, 0.9, 0.9, 1.0e-16,
			    0.0, 0, -1
	    );
	
	static interface LineSearchProc {
		public ReturnValue eval(
                int n,
                double[] x,
                MutableDouble f,
                double[] g,
                double[] s,
                MutableDouble stp,
                double[] xp,
                double[] gp,
                double[] wa,
                CallbackData cd,
                LBFGS_Param param
        );
	}
	
	/**
	 * Get L-BFGS's parameter structure with the default values.
	 *
	 *  Call this function to get a parameter structure with the default values
	 *  and then overwrite parameter values if necessary.
	 *
	 *  @return  The parameter structure.
	 */
	public static LBFGS_Param defaultParams() {
		return new LBFGS_Param(_defparam);
	}
	
	/*
	A user must implement a function compatible with ::LBFGS_Evaluate (evaluation
	callback) and pass the pointer to the callback function to lbfgs() arguments.
	Similarly, a user can implement a function compatible with ::LBFGS_Progress
	(progress callback) to obtain the current progress (e.g., variables, function
	value, ||G||, etc) and to cancel the iteration process if necessary.
	Implementation of a progress callback is optional: a user can pass \c NULL if
	progress notification is not necessary.

	This algorithm terminates an optimization
	when:

	    ||G|| < \epsilon \cdot \max(1, ||x||) .

	In this formula, ||.|| denotes the Euclidean norm.
	*/

	/**
	 * Start a L-BFGS optimization.
	 *
	 *  @param  n               The number of variables.
	 *  
	 *  @param  x               The array of variables. A client program can set
	 *                          default values for the optimization and receive the
	 *                          optimization result through this array. This array
	 *                          must be allocated by ::lbfgs_malloc function
	 *                          for libLBFGS built with SSE/SSE2 optimization routine
	 *                          enabled. The library built without SSE/SSE2
	 *                          optimization does not have such a requirement
	 *                          (Note: lbfgs_malloc() is removed in the Java port). 
	 *                      
	 *  @param  ptr_fx          The pointer to the variable that receives the final
	 *                          value of the objective function for the variables.
	 *                          This argument can be set to \c NULL if the final
	 *                          value of the objective function is unnecessary.
	 *                      
	 *  @param  proc_evaluate   The callback function to provide function and
	 *                          gradient evaluations given a current values of
	 *                          variables. A client program must implement a
	 *                          callback function compatible with \ref
	 *                          LBFGS_Evaluate and pass the pointer to the
	 *                          callback function.
	 *                          
	 *  @param  proc_progress   The callback function to receive the progress
	 *                          (the number of iterations, the current value of
	 *                          the objective function) of the minimization
	 *                          process. This argument can be set to \c NULL if
	 *                          a progress report is unnecessary.
	 *                          
	 *  @param  instance        A user data for the client program. The callback
	 *                          functions will receive the value of this argument.
	 *                      
	 *  @param  _param          The pointer to a structure representing parameters for
	 *                          L-BFGS optimization. A client program can set this
	 *                          parameter to \c NULL to use the default parameters.
	 *                          Call lbfgs_parameter_init() function to fill a
	 *                          structure with the default values.
	 *                      
	 *  @return int             The status code. This function returns zero if the
	 *                          minimization process terminates without an error. A
	 *                          non-zero value indicates an error.
	 */
	public static ReturnValue lbfgs(
		    int n,
		    double[] x,
		    MutableDouble ptr_fx,
		    LBFGS_Evaluate proc_evaluate,
		    LBFGS_Progress proc_progress,
		    Object instance,
		    LBFGS_Param _param
		    )
	{
		ReturnValue ret;
	    int i, j, k, end, bound;
	    ReturnValue ls;
	    MutableDouble step = new MutableDouble();
	    
	    /* Constant parameters and their default values. */
	    LBFGS_Param param = _param != null? _param : _defparam;
	    int m = param.m;
	    
	    double[] xp;
	    double[] g, gp, pg;
	    double[] d, w, pf;
	    IterationData[] lm;
	    IterationData   it;
	    double ys, yy;
	    double xnorm, gnorm, beta;
	    MutableDouble fx = new MutableDouble(0);
	    double rate = 0.;
	    LineSearchProc linesearch = new LineSearchMorethuente();
	    
	    /* Construct a callback data. */
	    CallbackData cd = new CallbackData();
	    cd.n = n;
	    cd.instance = instance;
	    cd.proc_evaluate = proc_evaluate;
	    cd.proc_progress = proc_progress;
	    
	    /* Check the input parameters for errors. */
	    if (n <= 0) {
	        return LBFGSERR_INVALID_N;
	    }
	    if (param.epsilon < 0.) {
	        return LBFGSERR_INVALID_EPSILON;
	    }
	    if (param.past < 0) {
	        return LBFGSERR_INVALID_TESTPERIOD;
	    }
	    if (param.delta < 0.) {
	        return LBFGSERR_INVALID_DELTA;
	    }
	    if (param.min_step < 0.) {
	        return LBFGSERR_INVALID_MINSTEP;
	    }
	    if (param.max_step < param.min_step) {
	        return LBFGSERR_INVALID_MAXSTEP;
	    }
	    if (param.ftol < 0.) {
	        return LBFGSERR_INVALID_FTOL;
	    }
	    if (param.lineSearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE ||
	        param.lineSearch == LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE) {
	        if (param.wolfe <= param.ftol || 1. <= param.wolfe) {
	            return LBFGSERR_INVALID_WOLFE;
	        }
	    }
	    if (param.gtol < 0) {
	        return LBFGSERR_INVALID_GTOL;
	    }
	    if (param.xtol < 0) {
	        return LBFGSERR_INVALID_XTOL;
	    }
	    if (param.max_linesearch <= 0) {
	        return LBFGSERR_INVALID_MAXLINESEARCH;
	    }
	    if (param.orthantwise_c < 0) {
	        return LBFGSERR_INVALID_ORTHANTWISE;
	    }
	    if (param.orthantwise_start < 0 || n < param.orthantwise_start) {
	        return LBFGSERR_INVALID_ORTHANTWISE_START;
	    }
	    if (param.orthantwise_end < 0) {
	        param.orthantwise_end = n;
	    }
	    if (n < param.orthantwise_end) {
	        return LBFGSERR_INVALID_ORTHANTWISE_END;
	    }
	    if (param.orthantwise_c != 0) {
	        switch (param.lineSearch) {
	        case LBFGS_LINESEARCH_BACKTRACKING:
	            linesearch = new LineSearchBacktrackingOWLQN();
	            break;
	        default:
	            /* Only the backtracking method is available. */
	            return LBFGSERR_INVALID_LINESEARCH;
	        }
	    } else {
	        switch (param.lineSearch) {
	        case LBFGS_LINESEARCH_DEFAULT:
	        case LBFGS_LINESEARCH_MORETHUENTE:
	            linesearch = new LineSearchMorethuente();
	            break;
	        case LBFGS_LINESEARCH_BACKTRACKING_ARMIJO:
	        case LBFGS_LINESEARCH_BACKTRACKING_WOLFE:
	        case LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE:
	            linesearch = new LineSearchBacktracking();
	            break;
	        default:
	            return LBFGSERR_INVALID_LINESEARCH;
	        }
	    }
	    
	    xp = new double[n];
	    g  = new double[n];
	    gp = new double[n];
	    d  = new double[n];
	    w  = new double[n];
	    
	    pg = null;
	    if (param.orthantwise_c != 0) {
	        /* Allocate working space for OWL-QN. */
	    	pg = new double[n];
	    }
	    
	    /* Allocate limited memory storage. */
	    lm = new IterationData[m];
	    
	    /* Initialize the limited memory. */
	    for (i = 0;i < m; ++i) {
	    	it = new IterationData();
	    	it.alpha = 0;
	    	it.ys = 0;
	    	it.s = new double[n];
	    	it.y = new double[n];
	    	lm[i] = it;
	    }

	    /* Allocate an array for storing previous values of the objective function. */
	    pf = null;
	    if (0 < param.past) {
	    	pf = new double[param.past];
	    }
	    
	    /* Evaluate the function value and its gradient. */
	    fx.val = cd.proc_evaluate.eval(cd.instance, x, g, cd.n, 0);
	    if (0 != param.orthantwise_c) {
	        /* Compute the L1 norm of the variable and add it to the object value. */
	        xnorm = owlqnX1Norm(x, param.orthantwise_start, param.orthantwise_end);
	        fx.val += xnorm * param.orthantwise_c;
	        owlqnPseudoGradient(
	            pg, x, g, n,
	            param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
	            );
	    }
	    
	    /* Store the initial value of the objective function. */
	    if (pf != null) {
	        pf[0] = fx.val;
	    }
	    
	    /**
	     *  Compute the direction;
	     *  we assume the initial hessian matrix H_0 as the identity matrix. 
	     */
	    if (param.orthantwise_c == 0.) {
	        vecncpy(d, g, n);
	    } else {
	        vecncpy(d, pg, n);
	    }
	    
	    /**
	     *  Make sure that the initial variables are not a minimizer.
	     */
	    xnorm = vec2norm(x, n);
	    if (param.orthantwise_c == 0.) {
	    	gnorm = vec2norm(g, n);
	    } else {
	    	gnorm = vec2norm(pg, n);
	    }
	    if (xnorm < 1.0) xnorm = 1.0;
	    if (gnorm / xnorm <= param.epsilon) {
	        return LBFGS_ALREADY_MINIMIZED;
	    }
	    
	    /** 
	     * Compute the initial step:
	     *  step = 1.0 / sqrt(vecdot(d, d, n))
	     */
	    step.val = vec2norminv(d, n);
	    
	    k = 1;
	    end = 0;
	    for (;;) {
	        /* Store the current position and gradient vectors. */
	        veccpy(xp, x, n);
	        veccpy(gp, g, n);

	        /* Search for an optimal step. */
	        if (param.orthantwise_c == 0.) {
	            ls = linesearch.eval(n, x, fx, g, d, step, xp, gp, w, cd, param);
	        } else {
	            ls = linesearch.eval(n, x, fx, g, d, step, xp, pg, w, cd, param);
	            owlqnPseudoGradient(
	                pg, x, g, n,
	                param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
	                );
	        }
	        if (ls.val < 0) {
	            /* Revert to the previous point. */
	            veccpy(x, xp, n);
	            veccpy(g, gp, n);
	            return ls;
	        }

	        /* Compute x and g norms. */
	        xnorm = vec2norm(x, n);
	        if (param.orthantwise_c == 0.) {
	        	gnorm = vec2norm(g, n);
	        } else {
	        	gnorm = vec2norm(pg, n);
	        }

	        /* Report the progress. */
	        if (cd.proc_progress != null) {
	            ret = cd.proc_progress.eval(cd.instance, x, g, fx.val, xnorm, gnorm, step.val, cd.n, k, ls.val);
	            if (ret.val > 0) {
	            	return ret;
	            }
	        }

	        /**
	         *  Convergence test.
	         *  The criterion is given by the following formula:
	         *      |g(x)| / \max(1, |x|) < \epsilon
	         */
	        if (xnorm < 1.0) xnorm = 1.0;
	        if (gnorm / xnorm <= param.epsilon) {
	            /* Convergence. */
	            ret = LBFGS_SUCCESS;
	            break;
	        }

	        /**
	         *  Test for stopping criterion.
	         *  The criterion is given by the following formula:
	         *      (f(past_x) - f(x)) / f(x) < \delta
	         */
	        if (pf != null) {
	            /* We don't test the stopping criterion while k < past. */
	            if (param.past <= k) {
	                /* Compute the relative improvement from the past. */
	                rate = (pf[k % param.past] - fx.val) / fx.val;

	                /* The stopping criterion. */
	                if (rate < param.delta) {
	                    ret = LBFGS_STOP;
	                    break;
	                }
	            }

	            /* Store the current value of the objective function. */
	            pf[k % param.past] = fx.val;
	        }

	        if (param.max_iterations != 0 && param.max_iterations < k+1) {
	            /* Maximum number of iterations. */
	            ret = LBFGSERR_MAXIMUMITERATION;
	            break;
	        }

	        /**
	         *  Update vectors s and y:
	         *      s_{k+1} = x_{k+1} - x_{k} = \step * d_{k}.
	         *      y_{k+1} = g_{k+1} - g_{k}.
	         */
	        it = lm[end];
	        vecdiff(it.s, x, xp, n);
	        vecdiff(it.y, g, gp, n);

	        /**
	         *  Compute scalars ys and yy:
	         *      ys = y^t \cdot s = 1 / \rho.
	         *      yy = y^t \cdot y.
	         *  Notice that yy is used for scaling the hessian matrix H_0 (Cholesky factor).
	         */
	        ys = vecdot(it.y, it.s, n);
	        yy = vecdot(it.y, it.y, n);
	        it.ys = ys;

	        /**
	         *  Recursive formula to compute dir = -(H \cdot g).
	         *      This is described in page 779 of:
	         *      Jorge Nocedal.
	         *      Updating Quasi-Newton Matrices with Limited Storage.
	         *      Mathematics of Computation, Vol. 35, No. 151,
	         *      pp. 773--782, 1980.
	         */
	        bound = (m <= k) ? m : k;
	        ++k;
	        end = (end + 1) % m;

	        /* Compute the steepest direction. */
	        if (param.orthantwise_c == 0.) {
	            /* Compute the negative of gradients. */
	            vecncpy(d, g, n);
	        } else {
	            vecncpy(d, pg, n);
	        }

	        j = end;
	        for (i = 0;i < bound;++i) {
	            j = (j + m - 1) % m;    /* if (--j == -1) j = m-1; */
	            it = lm[j];
	            /* \alpha_{j} = \rho_{j} s^{t}_{j} \cdot q_{k+1}. */
	            it.alpha = vecdot(it.s, d, n);
	            it.alpha /= it.ys;
	            /* q_{i} = q_{i+1} - \alpha_{i} y_{i}. */
	            vecadd(d, it.y, -it.alpha, n);
	        }

	        vecscale(d, ys / yy, n);

	        for (i = 0;i < bound;++i) {
	            it = lm[j];
	            /* \beta_{j} = \rho_{j} y^t_{j} \cdot \gamma_{i}. */
	            beta = vecdot(it.y, d, n);
	            beta /= it.ys;
	            /* \gamma_{i+1} = \gamma_{i} + (\alpha_{j} - \beta_{j}) s_{j}. */
	            vecadd(d, it.s, it.alpha - beta, n);
	            j = (j + 1) % m;        /* if (++j == m) j = 0; */
	        }

	        /**
	         *  Constrain the search direction for orthant-wise updates.
	         */
	        if (param.orthantwise_c != 0.) {
	            for (i = param.orthantwise_start;i < param.orthantwise_end;++i) {
	                if (d[i] * pg[i] >= 0) {
	                    d[i] = 0;
	                }
	            }
	        }

	        /**
	         *  Now the search direction d is ready. We try step = 1 first.
	         */
	        step.val = 1.0;
	    }

	    /* Return the final value of the objective function. */
	    if (ptr_fx != null) {
	        ptr_fx.val = fx.val;
	    }
	    
	    return ret;
	}
	
	static class LineSearchBacktracking implements LineSearchProc {
		public ReturnValue eval(
			    int n,
			    double[] x,
			    MutableDouble f,
			    double[] g,
			    double[] s,
			    MutableDouble stp,
			    double[] xp,
			    double[] gp,
			    double[] wp,
			    CallbackData cd,
			    LBFGS_Param param
			    )
		{
		    int count = 0;
		    double width, dg;
		    double finit, dginit = 0, dgtest;
		    double dec = 0.5, inc = 2.1;
		    
		    if (stp.val <= 0) {
		    	return LBFGSERR_INVALIDPARAMETERS;
		    }
		    
		    /* Compute the initial gradient in the search direction. */
		    dginit = vecdot(g, s, n);
	
		    /* Make sure that s points to a descent direction. */
		    if (0 < dginit) {
		        return LBFGSERR_INCREASEGRADIENT;
		    }
		    
		    /* The initial value of the objective function. */
		    finit = f.val;
		    dgtest = param.ftol * dginit;
	
		    for (;;) {
		        veccpy(x, xp, n);
		        vecadd(x, s, stp.val, n);
	
		        /* Evaluate the function and gradient values. */
		        f.val = cd.proc_evaluate.eval(cd.instance, x, g, cd.n, stp.val);
	
		        ++count;
	
		        if (f.val > finit + stp.val * dgtest) {
		            width = dec;
		        } else {
		            /* The sufficient decrease condition (Armijo condition). */
		            if (param.lineSearch == LBFGS_LINESEARCH_BACKTRACKING_ARMIJO) {
		                /* Exit with the Armijo condition. */
		                return LBFGS_SUCCESS;
			        }
	
			        /* Check the Wolfe condition. */
		          dg = vecdot(g, s, n);
			        if (dg < param.wolfe * dginit) {
		    		    width = inc;
			        } else {
				        if(param.lineSearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE) {
				            /* Exit with the regular Wolfe condition. */
				            return LBFGS_SUCCESS;
				        }
	
				        /* Check the strong Wolfe condition. */
				        if(dg > -param.wolfe * dginit) {
				            width = dec;
				        } else {
				            /* Exit with the strong Wolfe condition. */
				            return LBFGS_SUCCESS;
				        }
		          }
		        }
	
		        if (stp.val < param.min_step) {
		            /* The step is the minimum value. */
		            return LBFGSERR_MINIMUMSTEP;
		        }
		        if (stp.val > param.max_step) {
		            /* The step is the maximum value. */
		            return LBFGSERR_MAXIMUMSTEP;
		        }
		        if (param.max_linesearch <= count) {
		            /* Maximum number of iteration. */
		            return LBFGSERR_MAXIMUMLINESEARCH;
		        }
	
		        stp.val *= width;
		    }
		}
	}
	
	static class LineSearchBacktrackingOWLQN implements LineSearchProc {
		public ReturnValue eval(
			    int n,
			    double[] x,
			    MutableDouble f,
			    double[] g,
			    double[] s,
			    MutableDouble stp,
			    double[] xp,
			    double[] gp,
			    double[] wp,
			    CallbackData cd,
			    LBFGS_Param param
				)
		{
		    int i;
		    int count = 0;
		    double width = 0.5, norm = 0.;
		    double finit = f.val, dgtest;
	
		    /* Check the input parameters for errors. */
		    if (stp.val <= 0) {
		        return LBFGSERR_INVALIDPARAMETERS;
		    }
	
		    /* Choose the orthant for the new point. */
		    for (i = 0;i < n;++i) {
		        wp[i] = (xp[i] == 0.) ? -gp[i] : xp[i];
		    }
	
		    for (;;) {
		        /* Update the current point. */
		        veccpy(x, xp, n);
		        vecadd(x, s, stp.val, n);
	
		        /* The current point is projected onto the orthant. */
		        owlqnProject(x, wp, param.orthantwise_start, param.orthantwise_end);
	
		        /* Evaluate the function and gradient values. */
		        f.val = cd.proc_evaluate.eval(cd.instance, x, g, cd.n, stp.val);
	
		        /* Compute the L1 norm of the variables and add it to the object value. */
		        norm = owlqnX1Norm(x, param.orthantwise_start, param.orthantwise_end);
		        f.val += norm * param.orthantwise_c;
	
		        ++count;
	
		        dgtest = 0.;
		        for (i = 0;i < n;++i) {
		            dgtest += (x[i] - xp[i]) * gp[i];
		        }
	
		        if (f.val <= finit + param.ftol * dgtest) {
		            /* The sufficient decrease condition. */
		            return LBFGS_SUCCESS;
		        }
	
		        if (stp.val < param.min_step) {
		            /* The step is the minimum value. */
		            return LBFGSERR_MINIMUMSTEP;
		        }
		        if (stp.val > param.max_step) {
		            /* The step is the maximum value. */
		            return LBFGSERR_MAXIMUMSTEP;
		        }
		        if (param.max_linesearch <= count) {
		            /* Maximum number of iteration. */
		            return LBFGSERR_MAXIMUMLINESEARCH;
		        }
	
		        stp.val *= width;
		    }
		}
	}
	
	static class LineSearchMorethuente implements LineSearchProc {
		public ReturnValue eval(
			    int n,
			    double[] x,
			    MutableDouble f,
			    double[] g,
			    double[] s,
			    MutableDouble stp,
			    double[] xp,
			    double[] gp,
			    double[] wa,
			    CallbackData cd,
			    LBFGS_Param param
			    )
		{
		    int count = 0;
		    MutableDouble brackt;
		    int stage1;
		    ReturnValue uinfo = LBFGS_SUCCESS;
		    MutableDouble dg;
		    MutableDouble stx, fx, dgx;
		    MutableDouble sty, fy, dgy;
		    MutableDouble fxm, dgxm, fym, dgym;
		    MutableDouble fm, dgm;
		    double finit, ftest1, dginit, dgtest;
		    double width, prev_width;
		    double stmin, stmax;

		    // Initialize objects
		    dg   = new MutableDouble();
		    stx  = new MutableDouble();  fx   = new MutableDouble(); dgx = new MutableDouble();
		    sty  = new MutableDouble();  fy   = new MutableDouble(); dgy = new MutableDouble();
		    fxm  = new MutableDouble();  dgxm = new MutableDouble(); fym = new MutableDouble(); 
		    dgym = new MutableDouble();  fm   = new MutableDouble(); dgm = new MutableDouble(); 
		    
		    /* Check the input parameters for errors. */
		    if (stp.val <= 0.) {
		        return LBFGSERR_INVALIDPARAMETERS;
		    }
	
		    /* Compute the initial gradient in the search direction. */
		    dginit = vecdot(g, s, n);
	
		    /* Make sure that s points to a descent direction. */
		    if (0 < dginit) {
		        return LBFGSERR_INCREASEGRADIENT;
		    }
	
		    /* Initialize local variables. */
		    brackt = new MutableDouble(0);
		    stage1 = 1;
		    finit = f.val;
		    dgtest = param.ftol * dginit;
		    width = param.max_step - param.min_step;
		    prev_width = 2.0 * width;
	
		    /**
		     *  The variables stx, fx, dgx contain the values of the step,
		     *  function, and directional derivative at the best step.
		     *  The variables sty, fy, dgy contain the value of the step,
		     *  function, and derivative at the other endpoint of
		     *  the interval of uncertainty.
		     *  The variables stp, f, dg contain the values of the step,
		     *  function, and derivative at the current step.
		     */
		    stx.val = sty.val = 0;
		    fx.val = fy.val = finit;
		    dgx.val = dgy.val = dginit;
	
		    for (;;) {
		        /**
		         *  Set the minimum and maximum steps to correspond to the
		         *  present interval of uncertainty.
		         */
		        if (brackt.val > 0) {
		            stmin = min(stx.val, sty.val);
		            stmax = max(stx.val, sty.val);
		        } else {
		            stmin = stx.val;
		            stmax = stp.val + 4.0 * (stp.val - stx.val);
		        }
	
		        /* Clip the step in the range of [stpmin, stpmax]. */
		        if (stp.val < param.min_step) stp.val = param.min_step;
		        if (param.max_step < stp.val) stp.val = param.max_step;
	
		        /**
		         *  If an unusual termination is to occur then let
		         *  stp be the lowest point obtained so far.
		         */
		        if ((brackt.val > 0 && ((stp.val <= stmin || stmax <= stp.val) || param.max_linesearch <= count + 1 || uinfo.val != 0)) || (brackt.val > 0 && (stmax - stmin <= param.xtol * stmax))) {
		            stp.val = stx.val;
		        }
	
		        /**
		         *  Compute the current value of x:
		         *      x <- x + (*stp) * s.
		         */
		        veccpy(x, xp, n);
		        vecadd(x, s, stp.val, n);
	
		        /* Evaluate the function and gradient values. */
		        f.val  = cd.proc_evaluate.eval(cd.instance, x, g, cd.n, stp.val);
		        dg.val = vecdot(g, s, n);
	
		        ftest1 = finit + stp.val * dgtest;
		        ++count;
	
		        /* Test for errors and convergence. */
		        if (brackt.val > 0 && ((stp.val <= stmin || stmax <= stp.val) || uinfo.val != 0)) {
		            /* Rounding errors prevent further progress. */
		            return LBFGSERR_ROUNDING_ERROR;
		        }
		        if (stp.val == param.max_step && f.val <= ftest1 && dg.val <= dgtest) {
		            /* The step is the maximum value. */
		            return LBFGSERR_MAXIMUMSTEP;
		        }
		        if (stp.val == param.min_step && (ftest1 < f.val || dgtest <= dg.val)) {
		            /* The step is the minimum value. */
		            return LBFGSERR_MINIMUMSTEP;
		        }
		        if (brackt.val > 0 && (stmax - stmin) <= param.xtol * stmax) {
		            /* Relative width of the interval of uncertainty is at most xtol. */
		            return LBFGSERR_WIDTHTOOSMALL;
		        }
		        if (param.max_linesearch <= count) {
		            /* Maximum number of iteration. */
		            return LBFGSERR_MAXIMUMLINESEARCH;
		        }
		        if (f.val <= ftest1 && abs(dg.val) <= param.gtol * (-dginit)) {
		            /* The sufficient decrease condition and the directional derivative condition hold. */
		            return LBFGS_SUCCESS;
		        }
	
		        /**
		         *  In the first stage we seek a step for which the modified
		         *  function has a nonpositive value and nonnegative derivative.
		         */
		        if (stage1 > 0 && f.val <= ftest1 && min(param.ftol, param.gtol) * dginit <= dg.val) {
		            stage1 = 0;
		        }
	
		        /**
		         *  A modified function is used to predict the step only if
		         *  we have not obtained a step for which the modified
		         *  function has a nonpositive function value and nonnegative
		         *  derivative, and if a lower function value has been
		         *  obtained but the decrease is not sufficient.
		         */
		        if (stage1 > 0 && ftest1 < f.val && f.val <= fx.val) {
		            /* Define the modified function and derivative values. */
		            fm.val   = f.val - stp.val * dgtest;
		            fxm.val  = fx.val - stx.val * dgtest;
		            fym.val  = fy.val - sty.val * dgtest;
		            dgm.val  = dg.val - dgtest;
		            dgxm.val = dgx.val - dgtest;
		            dgym.val = dgy.val - dgtest;
	
		            /**
		             *  Call update_trial_interval() to update the interval of
		             *  uncertainty and to compute the new step.
		             */
		            uinfo = updateTrialInterval(
		                stx, fxm, dgxm,
		                sty, fym, dgym,
		                stp, fm, dgm,
		                stmin, stmax, brackt
		                );
	
		            /* Reset the function and gradient values for f. */
		            fx.val = fxm.val + stx.val * dgtest;
		            fy.val = fym.val + sty.val * dgtest;
		            dgx.val = dgxm.val + dgtest;
		            dgy.val = dgym.val + dgtest;
		        } else {
		            /**
		             *  Call update_trial_interval() to update the interval of
		             *  uncertainty and to compute the new step.
		             */
		            uinfo = updateTrialInterval(
		                stx, fx, dgx,
		                sty, fy, dgy,
		                stp, f, dg,
		                stmin, stmax, brackt
		                );
		        }
	
		        /**
		         *  Force a sufficient decrease in the interval of uncertainty.
		         */
		        if (brackt.val > 0) {
		            if (0.66 * prev_width <= abs(sty.val - stx.val)) {
		                stp.val = stx.val + 0.5 * (sty.val - stx.val);
		            }
		            prev_width = width;
		            width = abs(sty.val - stx.val);
		        }
		    }
	
		    // return LBFGSERR_LOGICERROR; // Unreachable code
		}
	}
	
	/**
	 * Find a minimizer of an interpolated cubic function.
	 *  @param  u       The value of one point, u.
	 *  @param  fu      The value of f(u).
	 *  @param  du      The value of f'(u).
	 *  @param  v       The value of another point, v.
	 *  @param  fv      The value of f(v).
	 *  @param  du      The value of f'(v).
	 *  @return double	The minimizer of the interpolated cubic.
	 */
	static double cubicMinimizer(
			double u, 
			double fu, 
			double du, 
			double v, 
			double fv, 
			double dv
			) 
	{
		double d, theta, p, q, r, s, a, gamma;
		double cm;
	    d = v - u;
	    theta = (fu - fv) * 3 / d + du + dv;
	    p = abs(theta);
	    q = abs(du);
	    r = abs(dv);
	    s = max(max(p, q), r);
	    /* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */
	    a = theta / s;
	    gamma = s * sqrt(a * a - (du / s) * (dv / s));
	    if ((v) < (u)) gamma = -gamma;
	    p = gamma - du + theta;
	    q = gamma - du + gamma + dv;
	    r = p / q;
	    cm = u + r * d;
	    return cm;
	}
	
	/**
	 * Find a minimizer of an interpolated cubic function.
	 *  @param  u       The value of one point, u.
	 *  @param  fu      The value of f(u).
	 *  @param  du      The value of f'(u).
	 *  @param  v       The value of another point, v.
	 *  @param  fv      The value of f(v).
	 *  @param  du      The value of f'(v).
	 *  @param  xmin    The minimum value.
	 *  @param  xmax    The maximum value.
	 *  @return double	The minimizer of the interpolated cubic.
	 */
	static double cubicMinimizer2(
			double u, 
			double fu, 
			double du, 
			double v, 
			double fv, 
			double dv, 
			double xmin, 
			double xmax
			)
	{
		double d, theta, p, q, r, s, a, gamma;
		double cm;
	    d = v - u;
	    theta = (fu - fv) * 3 / d + du + dv;
	    p = abs(theta);
	    q = abs(du);
	    r = abs(dv);
	    s = max(max(p, q), r);
	    /* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */
	    a = theta / s;
	    gamma = s * sqrt(max(0, a * a - (du / s) * (dv / s)));
	    if (u < v) gamma = -gamma;
	    p = gamma - dv + theta;
	    q = gamma - dv + gamma + du;
	    r = p / q;
	    if (r < 0. && gamma != 0.) {
	        cm = v - r * d;
	    } else if (a < 0) {
	        cm = xmax;
	    } else {
	        cm = xmin;
	    }
	    return cm;
	}
	
	/**
	 * Find a minimizer of an interpolated quadratic function.
	 *  @param  u       The value of one point, u.
	 *  @param  fu      The value of f(u).
	 *  @param  du      The value of f'(u).
	 *  @param  v       The value of another point, v.
	 *  @param  fv      The value of f(v).
	 *  @return	double  The minimizer of the interpolated quadratic.
	 */
	static double quardMinimizer(
			double u, 
			double fu, 
			double du, 
			double v, 
			double fv
			) 
	{
		double a = v - u;
	    double qm = u + du / ((fu - fv) / a + du) / 2 * a;
	    return qm;
	}
	
	/**
	 * Find a minimizer of an interpolated quadratic function.
	 *  @param  u       The value of one point, u.
	 *  @param  du      The value of f'(u).
	 *  @param  v       The value of another point, v.
	 *  @param  dv      The value of f'(v).
	 *  @return double  The minimizer of the interpolated quadratic.
	 */
	static double quardMinimizer2(
			double u, 
			double du, 
			double v, 
			double dv
			) 
	{
		double a = u - v;
		double qm = v + dv / (dv - du) * a;
		return qm;
	}
	
	/**
	 * Update a safeguarded trial value and interval for line search.
	 *
	 *  The parameter x represents the step with the least function value.
	 *  The parameter t represents the current step. This function assumes
	 *  that the derivative at the point of x in the direction of the step.
	 *  If the bracket is set to true, the minimizer has been bracketed in
	 *  an interval of uncertainty with endpoints between x and y.
	 *
	 *  @param  x       The pointer to the value of one endpoint.
	 *  @param  fx      The pointer to the value of f(x).
	 *  @param  dx      The pointer to the value of f'(x).
	 *  @param  y       The pointer to the value of another endpoint.
	 *  @param  fy      The pointer to the value of f(y).
	 *  @param  dy      The pointer to the value of f'(y).
	 *  @param  t       The pointer to the value of the trial value, t.
	 *  @param  ft      The pointer to the value of f(t).
	 *  @param  dt      The pointer to the value of f'(t).
	 *  @param  tmin    The minimum value for the trial value, t.
	 *  @param  tmax    The maximum value for the trial value, t.
	 *  @param  brackt  The pointer to the predicate if the trial value is
	 *                  bracketed.
	 *  @retrun int     Status value. Zero indicates a normal termination.
	 *  
	 *  @see
	 *      Jorge J. More and David J. Thuente. Line search algorithm with
	 *      guaranteed sufficient decrease. ACM Transactions on Mathematical
	 *      Software (TOMS), Vol 20, No 3, pp. 286-307, 1994.
	 */
	static ReturnValue updateTrialInterval(
  			MutableDouble x,
  			MutableDouble fx,
  			MutableDouble dx,
		    MutableDouble y,
		    MutableDouble fy,
		    MutableDouble dy,
		    MutableDouble t,
		    MutableDouble ft,
		    MutableDouble dt,
		    double tmin,
		    double tmax,
		    MutableDouble brackt
			)
	{
	    int bound;
	    boolean dsign = fsigndiff(dt.val, dx.val);
	    double mc;              /* minimizer of an interpolated cubic. */
	    double mq;              /* minimizer of an interpolated quadratic. */
	    double newt;            /* new trial value. */
	    // double a, d, gamma, theta, p, q, r, s;    /* for CUBIC_MINIMIZER and QUARD_MINIMIZER. */

	    /* Check the input parameters for errors. */
	    if (brackt.val > 0) {
	        if (t.val <= min(x.val, y.val) || max(x.val, y.val) <= t.val) {
	            /* The trival value t is out of the interval. */
	            return LBFGSERR_OUTOFINTERVAL;
	        }
	        if (0 <= dx.val * (t.val - x.val)) {
	            /* The function must decrease from x. */
	            return LBFGSERR_INCREASEGRADIENT;
	        }
	        if (tmax < tmin) {
	            /* Incorrect tmin and tmax specified. */
	            return LBFGSERR_INCORRECT_TMINMAX;
	        }
	    }

	    /**
	     * Trial value selection.
	     */
	    if (fx.val < ft.val) {
	        /**
	         *  Case 1: a higher function value.
	         *  The minimum is brackt. If the cubic minimizer is closer
	         *  to x than the quadratic one, the cubic one is taken, else
	         *  the average of the minimizers is taken.
	         */
	        brackt.val = 1;
	        bound = 1;
	        mc = cubicMinimizer(x.val, fx.val, dx.val, t.val, ft.val, dt.val);
	        mq = quardMinimizer(x.val, fx.val, dx.val, t.val, ft.val);
	        if (abs(mc - x.val) < abs(mq - x.val)) {
	            newt = mc;
	        } else {
	            newt = mc + 0.5 * (mq - mc);
	        }
	    } else if (dsign) {
	        /**
	         *  Case 2: a lower function value and derivatives of
	         *  opposite sign. The minimum is brackt. If the cubic
	         *  minimizer is closer to x than the quadratic (secant) one,
	         *  the cubic one is taken, else the quadratic one is taken.
	         */
	        brackt.val = 1;
	        bound = 0;
	        mc = cubicMinimizer(x.val, fx.val, dx.val, t.val, ft.val, dt.val);
	        mq = quardMinimizer2(x.val, dx.val, t.val, dt.val);
	        if (abs(mc - t.val) > abs(mq - t.val)) {
	            newt = mc;
	        } else {
	            newt = mq;
	        }
	    } else if (abs(dt.val) < abs(dx.val)) {
	        /**
	         *  Case 3: a lower function value, derivatives of the
	         *  same sign, and the magnitude of the derivative decreases.
	         *  The cubic minimizer is only used if the cubic tends to
	         *  infinity in the direction of the minimizer or if the minimum
	         *  of the cubic is beyond t. Otherwise the cubic minimizer is
	         *  defined to be either tmin or tmax. The quadratic (secant)
	         *  minimizer is also computed and if the minimum is brackt
	         *  then the the minimizer closest to x is taken, else the one
	         *  farthest away is taken.
	         */
	        bound = 1;
	        mc = cubicMinimizer2(x.val, fx.val, dx.val, t.val, ft.val, dt.val, tmin, tmax);
	        mq = quardMinimizer2(x.val, dx.val, t.val, dt.val);
	        if (brackt.val > 0) {
	            if (abs(t.val - mc) < abs(t.val - mq)) {
	                newt = mc;
	            } else {
	                newt = mq;
	            }
	        } else {
	            if (abs(t.val - mc) > abs(t.val - mq)) {
	                newt = mc;
	            } else {
	                newt = mq;
	            }
	        }
	    } else {
	        /**
	         *  Case 4: a lower function value, derivatives of the
	         *  same sign, and the magnitude of the derivative does
	         *  not decrease. If the minimum is not brackt, the step
	         *  is either tmin or tmax, else the cubic minimizer is taken.
	         */
	        bound = 0;
	        if (brackt.val > 0) {
	        	newt = cubicMinimizer(t.val, ft.val, dt.val, y.val, fy.val, dy.val);
	        } else if (x.val < t.val) {
	            newt = tmax;
	        } else {
	            newt = tmin;
	        }
	    }

	    /**
	     *  Update the interval of uncertainty. This update does not
	     *  depend on the new step or the case analysis above.
		 *  
	     *  - Case a: if f(x) < f(t),
	     *      x <- x, y <- t.
	     *  - Case b: if f(t) <= f(x) && f'(t)*f'(x) > 0,
	     *      x <- t, y <- y.
	     *  - Case c: if f(t) <= f(x) && f'(t)*f'(x) < 0, 
	     *      x <- t, y <- x.
	     */
	    if (fx.val < ft.val) {
	        /* Case a */
	        y.val  = t.val;
	        fy.val = ft.val;
	        dy.val = dt.val;
	    } else {
	        /* Case c */
	        if (dsign) {
	            y.val  = x.val;
	            fy.val = fx.val;
	            dy.val = dx.val;
	        }
	        /* Cases b and c */
	        x.val  = t.val;
	        fx.val = ft.val;
	        dx.val = dt.val;
	    }

	    /* Clip the new trial value in [tmin, tmax]. */
	    if (tmax < newt) newt = tmax;
	    if (newt < tmin) newt = tmin;

	    /**
	     *  Redefine the new trial value if it is close to the upper bound
	     *  of the interval.
	     */
	    if (brackt.val > 0 && bound > 0) {
	        mq = x.val + 0.66 * (y.val - x.val);
	        if (x.val < y.val) {
	            if (mq < newt) newt = mq;
	        } else {
	            if (newt < mq) newt = mq;
	        }
	    }

	    /* Return the new trial value. */
	    t.val = newt;
	    return LBFGS_SUCCESS; // or maybe LBFGS_CONVERGENCE since both have value 0
	}
	
	static double owlqnX1Norm(
		    double[] x,
		    int start,
		    int n
		    )
	{
	    int i;
	    double norm = 0.;

	    for (i = start;i < n;++i) {
	        norm += abs(x[i]);
	    }

	    return norm;
	}

	static void owlqnPseudoGradient(
	    double[] pg,
	    double[] x,
	    double[] g,
	    int n,
	    double c,
	    int start,
	    int end
	    )
	{
	    int i;

	    /* Compute the negative of gradients. */
	    for (i = 0;i < start;++i) {
	        pg[i] = g[i];
	    }

	    /* Compute the pseudo-gradients. */
	    for (i = start;i < end;++i) {
	        if (x[i] < 0.) {
	            /* Differentiable. */
	            pg[i] = g[i] - c;
	        } else if (0. < x[i]) {
	            /* Differentiable. */
	            pg[i] = g[i] + c;
	        } else {
	            if (g[i] < -c) {
	                /* Take the right partial derivative. */
	                pg[i] = g[i] + c;
	            } else if (c < g[i]) {
	                /* Take the left partial derivative. */
	                pg[i] = g[i] - c;
	            } else {
	                pg[i] = 0.;
	            }
	        }
	    }

	    for (i = end;i < n;++i) {
	        pg[i] = g[i];
	    }
	}

	static void owlqnProject(
	    double[] d,
	    double[] sign,
	    int start,
	    int end
	    )
	{
	    int i;

	    for (i = start;i < end;++i) {
	        if (d[i] * sign[i] <= 0) {
	            d[i] = 0;
	        }
	    }
	}
}
