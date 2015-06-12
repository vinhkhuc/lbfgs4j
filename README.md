This is a direct Java port of [liblbfgs](http://www.chokkan.org/software/liblbfgs/), a library of Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS).

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
double[] x = minimizer.minimize(f); // x should be {5}
double min = f.valueAt(x);          // min should be 1
```

### License

lbfgs4j is released under the MIT License (see the LICENSE file for details).