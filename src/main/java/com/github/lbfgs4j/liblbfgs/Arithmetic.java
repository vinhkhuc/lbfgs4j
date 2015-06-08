package com.github.lbfgs4j.liblbfgs;

/*
 *      Implementation of vector operations.
 *
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

/* $Id$ */

public class Arithmetic {

	public static boolean fsigndiff(double x, double y) {
		return (x * (y / Math.abs(y)) < 0);
	}
	
	public static void vecset(double[] x, double c, int n) {
		for (int i = 0; i < n; i++) {
			x[i] = c;
		}
	}
	
	public static void veccpy(double[] y, double[] x, int n) {
		for (int i = 0; i < n; i++) {
			y[i] = x[i];
		}
	}
	
	public static void vecncpy(double[] y, double[] x, int n) {
		for (int i = 0; i < n; i++) {
			y[i] = -x[i];
		}
	}
	
	public static void vecadd(double[] y, double[] x, double c, int n) {
		for (int i = 0; i < n; i++) {
			y[i] += c * x[i];
		}
	}
	
	public static void vecdiff(double[] z, double[] x, double[] y, int n) {
		for (int i = 0; i < n; i++) {
			z[i] = x[i] - y[i];
		}
	}
	
	public static void vecscale(double[] y, double c, int n) {
		for (int i = 0; i < n; i++) {
			y[i] *= c;
		}
	}
	
	public static void vecmul(double[] y, double[] x, int n) {
		for (int i = 0; i < n; i++) {
			y[i] *= x[i];
		}
	}
	
	public static double vecdot(double[] x, double[] y, int n) {
		double s = 0;
		for (int i = 0; i < n; i++) {
			s += x[i] * y[i];
		}
		return s;
	}
	
	public static double vec2norm(double[] x, int n) {
		return Math.sqrt(vecdot(x, x, n));
	}
	
	public static double vec2norminv(double[] x, int n) {
		return 1 / vec2norm(x, n);
	}
}
