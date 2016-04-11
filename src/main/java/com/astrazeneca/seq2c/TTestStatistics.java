package com.astrazeneca.seq2c;

import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;

public class TTestStatistics {

    public static double getT(double[] x, double[] y) {

        int significance = 95;
        double alpha = (double)(100-significance)/2;
        alpha /= 100;

        double meanDifference = StatUtils.mean(x) - StatUtils.mean(y);
        int df1= x.length - 1;
        int df2 = y.length - 1;

        double variance1 = StatUtils.variance(x);
        double variance2 = StatUtils.variance(y);
        double variance = (df1 * variance1 + df2 * variance2) / (df1 + df2);

        double standartDeviation = FastMath.sqrt(variance);

        StandardDeviation stDeviation = new StandardDeviation();
        double standartError1 = stDeviation.evaluate(x) / FastMath.sqrt(x.length);
        double standartError2 = stDeviation.evaluate(y) / FastMath.sqrt(y.length);

        double standartErrorEqual = FastMath.sqrt(1.0/x.length + 1.0/y.length) * standartDeviation;
        double standartErrorUnEqual = FastMath.sqrt(FastMath.pow(standartError1, 2) + FastMath.pow(standartError2, 2));

        double dfEqual = df1 + df2;
        double dfUnequal = (FastMath.pow(standartError1, 4) / df1 + FastMath.pow(standartError2, 4) / df2) != 0 ?
                (FastMath.pow(FastMath.pow(standartError1, 2) + FastMath.pow(standartError2, 2), 2) /
                        (FastMath.pow(standartError1, 4) / df1 + FastMath.pow(standartError2, 4) / df2))
                : dfEqual;

        double tStatEqual = standartErrorEqual != 0 ? FastMath.abs(meanDifference/standartErrorEqual) : 99;
        double tStatUnequal = standartErrorUnEqual != 0 ? FastMath.abs(meanDifference/standartErrorUnEqual) : 99;

        double fStat;

        if (Double.compare(variance1, variance2) >= 0) {
            fStat = variance2 != 0 ? variance1 / variance2 : 99;
        } else {
            int tmp = df1;
            df1 = df2;
            df2 = tmp;
            fStat = variance1 != 0 ? variance2 / variance1 : 99;
        }

        double fCutOff = FCriticalValue(alpha, df1, df2);
        double t;
        double degreesOfFreedom;
        if (Double.compare(fStat, fCutOff) <= 0) {
            t = Precision.round(tStatEqual, 12);
            degreesOfFreedom = dfEqual;
        } else {
            t = Precision.round(tStatUnequal, 12);
            degreesOfFreedom = dfUnequal;
        }

        TDistribution dist = new TDistribution(FastMath.floor(degreesOfFreedom));
        return 1 - FastMath.abs(Precision.round(dist.cumulativeProbability(t), 5) - dist.cumulativeProbability(-t));
    }

    public static double FCriticalValue(double p, int df1, int df2) {

        double  fval;
        double  maxf = 99999.0;     /* maximum possible F ratio */
        double  minf = .000001;     /* minimum possible F ratio */

        if (p <= 0.0 || p >= 1.0)
            return (0.0);

        fval = 1.0 / p; /* the smaller the p, the larger the F */

        while (FastMath.abs(maxf - minf) > .000001) {
            if (FProbability(fval, df1, df2) < p) /* F too large */
                maxf = fval;
            else /* F too small */
                minf = fval;
            fval = (maxf + minf) * 0.5;
        }

        return (fval);
    }

    public static double FProbability(double F, int df1, int df2) {

        int     i, j;
        int     a, b;
        double  w, y, z, d, p;

        if ((FastMath.abs(F) < 10e-10) || df1 <= 0 || df2 <= 0)
            return (1.0);
        a = (df1%2 == 1) ? 1 : 2;
        b = (df2%2 == 1) ? 1 : 2;
        w = (F * df1) / df2;
        z = 1.0 / (1.0 + w);
        if (a == 1)
            if (b == 1) {
                p = FastMath.sqrt (w);
                y = 1/Math.PI; /* 1 / 3.14159 */
                d = y * z / p;
                p = 2.0 * y * Math.atan (p);
            } else {
                p = FastMath.sqrt (w * z);
                d = 0.5 * p * z / w;
            } else if (b == 1) {
            p = FastMath.sqrt (z);
            d = 0.5 * z * p;
            p = 1.0 - p;
        } else {
            d = z * z;
            p = w * z;
        }
        y = 2.0 * w / z;
        for (j = b + 2; j <= df2; j += 2) {
            d *= (1.0 + a / (j - 2.0)) * z;
            p = (a == 1 ? p + d * y / (j - 1.0) : (p + w) * z);
        }
        y = w * z;
        z = 2.0 / z;
        b = df2 - 2;
        for (i = a + 2; i <= df1; i += 2) {
            j = i + b;
            d *= y * j / (i - 2.0);
            p -= z * d / j;
        }

        // correction for approximation errors suggested in certification
        if (p < 0.0)
            p = 0.0;
        else if (p > 1.0)
            p = 1.0;
        return (1.0-p);
    }
}
