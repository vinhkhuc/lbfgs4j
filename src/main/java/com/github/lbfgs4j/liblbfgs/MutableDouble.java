package com.github.lbfgs4j.liblbfgs;

public class MutableDouble {
  
	public double val;
	
	public MutableDouble() { 
		this(0); 
	}
	
	public MutableDouble(double val) { 
		this.val = val;
	}
}