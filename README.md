This is a direct Java port of [liblbfgs](http://www.chokkan.org/software/liblbfgs/), a library of Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS).

The L-BFGS optimization method is used to solve unconstrained minimization problem:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;minimize f(x), where x = (x_1, x_2, ... x_n)

with the condition that f(x) and its gradient g(x) can be computed.

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

### License

lbfgs4j is released under the MIT License (see the LICENSE file for details).