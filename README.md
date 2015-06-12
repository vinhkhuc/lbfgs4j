This is a direct Java port of [liblbfgs](http://www.chokkan.org/software/liblbfgs/), a library of Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS).

The L-BFGS optimization method is used to solve unconstrained minimization problem:
<p>
minimize f(x), where x = (x1, x2, ... xn), f(x) and its gradient g(x) are computable.
</p>

This library also includes the implementation of the Orthant-Wise Limited-memory Quasi-Newton (OWL-QN) method which is used to minimize the function f(x) + C|x| where C is a positive constant number.

### Maven dependency
```xml
<dependency>
  <groupId>com.github.vinhkhuc</groupId>
  <artifactId>lbfgs4j</artifactId>
  <version>0.1</version>
</dependency>
```

### Usage example
Find the minimum value of the function f(x) = (x-5)^2 + 1.
```java
import com.github.lbfgs4j.liblbfgs.Function;
import com.github.lbfgs4j.LbfgsMinimizer;
...

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

LbfgsMinimizer minimizer = new LbfgsMinimizer();
double[] x = minimizer.minimize(f); // x should be [5]
double min = f.valueAt(x);          // min should be 1
```

The OWL-QN method will be used when initializing the minimizer with, for example, ```new LbfgsMinimizer(1.0)```, here 1.0 is the coefficient for |x|.

Other parameters such as memory size, maximum number of iterations, etc. can be set using ```LBFGS_Param```. 

### License

lbfgs4j is released under the MIT License (see the LICENSE file for details).