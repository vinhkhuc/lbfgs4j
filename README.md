This is a direct Java port of [liblbfgs](http://www.chokkan.org/software/liblbfgs/), a library of <strong>Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS)</strong>.

The L-BFGS optimization method is used to solve unconstrained minimization problem:
<p>
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
minimize f(x), where x = (x1, x2, ... xn) <br>
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
f(x) and its gradient are computable.
</p>

This library also includes the implementation of the <strong>Orthant-Wise Limited-memory Quasi-Newton (OWL-QN) </strong> method which is used to minimize the function f(x) + C|x|, where C is a positive number.

### Maven dependency
```xml
<dependency>
  <groupId>com.github.vinhkhuc</groupId>
  <artifactId>lbfgs4j</artifactId>
  <version>0.2.1</version>
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

Other parameters such as memory size, maximum number of iterations, etc. can be set using the class ```LBFGS_Param```. 

### License

lbfgs4j is released under the MIT License (see the LICENSE file for details).

### References

This library is used in several projects:
1. Moment-Based Quantile Sketches for Efficient High Cardinality Aggregation Queries. Edward Gan, Jialin Ding, Kai Sheng Tai, Vatsal Sharan, Peter Bailis. VLDB 2018. ([paper](https://dl.acm.org/doi/10.14778/3236187.3236212))
2. Humio, Log Management at Scale. ([acknowledgement](https://docs.humio.com/third-party-licenses/))
