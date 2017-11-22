using System;

namespace EigenfacesForRecognition
{
    class EigenVV
    {
        private byte[,,] phiRaster;
        private double[,] L;
        private int M, h, m, w;

        public EigenVV(byte[,,] phiRaster, int M, int h, int m, int w)
        {
            this.phiRaster = phiRaster;
            this.M = M;
            this.h = h;
            this.m = m;
            this.w = w;
            L = new double[M, M];

            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    double sum = 0;

                    for (int k = 0; k < h; k++)
                    {
                        for (int l = 0; l < w; l++)
                        {
                            sum += phiRaster[i, k, l] * phiRaster[j, k, l];
                        }
                    }

                    L[i, j] = sum;
                }
            }
        }

        // functions translated from C code found in
        // "A Numerical Library in C for Scientists and Engineers"
        // by H.T. Lau, PhD Chapters 1, 2, and 3

        private void ichcol(int l, int u, int i, int j, double[,] a)
        {
            double r;

            for (; l <= u; l++)
            {
                r = a[l, i];
                a[l, i] = a[l, j];
                a[l, j] = r;
            }
        }

        private void ichrow(int l, int u, int i, int j, double[,] a)
        {
            double r;

            for (; l <= u; l++)
            {
                r = a[i, l];
                a[i, l] = a[j, l];
                a[j, l] = r;
            }
        }

        private double matmat(int l, int u, int i, int j, double[,] a, double[,] b)
        {
            int k;
            double s;

            s = 0.0;
            for (k = l; k <= u; k++) s += a[i, k] * b[k, j];
            return (s);
        }

        private double mattam(int l, int u, int i, int j, double[,] a, double[,] b)
        {
            int k;
            double s = 0;

            for (k = l; k <= u; k++) s += a[i, k] * b[j, k];
            return (s);
        }

        private double matvec(int l, int u, int i, double[,] a, double[] b)
        {
            int k;
            double s;

            s = 0.0;
            for (k = l; k <= u; k++) s += a[i, k] * b[k];
            return (s);
        }

        private void rotcol(int l, int u, int i, int j, double[,] a, double c, double s)
        {
            double x, y;

            for (; l <= u; l++)
            {
                x = a[l, i];
                y = a[l, j];
                a[l, i] = x * c + y * s;
                a[l, j] = y * c - x * s;
            }
        }

        private void rotrow(int l, int u, int i, int j, double[,] a, double c, double s)
        {
            double x, y;

            for (; l <= u; l++)
            {
                x = a[i, l];
                y = a[j, l];
                a[i, l] = x * c + y * s;
                a[j, l] = y * c - x * s;
            }
        }

        private double tammat(int l, int u, int i, int j, double[,] a, double[,] b)
        {
            int k;
            double s = 0.0;

            for (k = l; k <= u; k++) s += a[k, i] * b[k, j];
            return (s);
        }

        private double tamvec(int l, int u, int i, double[,] a, double[] b)
        {
            int k;
            double s;

            s = 0.0;
            for (k = l; k <= u; k++) s += a[k, i] * b[k];
            return (s);
        }

        private void tfmreahes(double[,] a, int n, double[] em, int[] index)
        {
            int i, j, j1, k, l;
            double s, t, machtol, macheps, norm;

            double[] b = new double[n];
            macheps = em[0];
            norm = 0.0;
            for (i = 1; i <= n; i++)
            {
                s = 0.0;
                for (j = 1; j <= n; j++) s += Math.Abs(a[i, j]);
                if (s > norm) norm = s;
            }
            em[1] = norm;
            machtol = norm * macheps;
            index[1] = 0;
            for (j = 2; j <= n; j++)
            {
                j1 = j - 1;
                l = 0;
                s = machtol;
                for (k = j + 1; k <= n; k++)
                {
                    t = Math.Abs(a[k, j1]);
                    if (t > s)
                    {
                        l = k;
                        s = t;
                    }
                }
                if (l != 0)
                {
                    if (Math.Abs(a[j, j1]) < s)
                    {
                        ichrow(1, n, j, l, a);
                        ichcol(1, n, j, l, a);
                    }
                    else
                        l = j;
                    t = a[j, j1];
                    for (k = j + 1; k <= n; k++) a[k, j1] /= t;
                }
                else
                    for (k = j + 1; k <= n; k++) a[k, j1] = 0.0;
                for (i = 1; i <= n; i++)
                    b[i - 1] = a[i, j] +=
                        ((l == 0) ? 0.0 : matmat(j + 1, n, i, j1, a, a)) -
                        matvec(1, (j1 < i - 2) ? j1 : i - 2, i, a, b);
                index[j] = l;
            }
        }

        private void bakreahes2(double[,] a, int n, int n1, int n2, int[] index,
                       double[,] vec)
        {
            int i, l, k;
            double[] u = new double[n + 1];

            for (i = n; i >= 2; i--)
            {
                for (k = i - 2; k >= 1; k--) u[k + 1] = a[i, k];
                for (k = n1; k <= n2; k++) vec[i, k] += tamvec(2, i - 1, k, vec, u);
                l = index[i];
                if (l > i) ichrow(n1, n2, i, l, vec);
            }
        }

        private void eqilbr(double[,] a, int n, double[] em, double[] d, int[] inter)
        {
            int i, im, i1, p, q, j, t, count, exponent, ni;
            double c, r, eps, omega, factor, di;

            factor = 1.0 / (2.0 * Math.Log(2.0));
            eps = em[0];
            omega = 1.0 / eps;
            t = p = 1;
            q = ni = i = n;
            count = ((n + 1) * n) / 2;
            for (j = 1; j <= n; j++)
            {
                d[j] = 1.0;
                inter[j] = 0;
            }
            i = (i < q) ? i + 1 : p;
            while (count > 0 && ni > 0)
            {
                count--;
                im = i - 1;
                i1 = i + 1;
                c = Math.Sqrt(tammat(p, im, i, i, a, a) + tammat(i1, q, i, i, a, a));
                r = Math.Sqrt(mattam(p, im, i, i, a, a) + mattam(i1, q, i, i, a, a));
                if (c * omega <= r * eps)
                {
                    inter[t] = i;
                    ni = q - p;
                    t++;
                    if (p != i)
                    {
                        ichcol(1, n, p, i, a);
                        ichrow(1, n, p, i, a);
                        di = d[i];
                        d[i] = d[p];
                        d[p] = di;
                    }
                    p++;
                }
                else
                    if (r * omega <= c * eps)
                    {
                        inter[t] = -i;
                        ni = q - p;
                        t++;
                        if (q != i)
                        {
                            ichcol(1, n, q, i, a);
                            ichrow(1, n, q, i, a);
                            di = d[i];
                            d[i] = d[q];
                            d[q] = di;
                        }
                        q--;
                    }
                    else
                    {
                        exponent = (int)(Math.Log(r / c) * factor);
                        if (Math.Abs(exponent) > 1.0)
                        {
                            ni = q - p;
                            c = Math.Pow(2.0, exponent);
                            r = 1.0 / c;
                            d[i] *= c;
                            for (j = 1; j <= im; j++)
                            {
                                a[j, i] *= c;
                                a[i, j] *= r;
                            }
                            for (j = i1; j <= n; j++)
                            {
                                a[j, i] *= c;
                                a[i, j] *= r;
                            }
                        }
                        else
                            ni--;
                    }
                i = (i < q) ? i + 1 : p;
            }
        }

        private void baklbr(int n, int n1, int n2, double[] d, int[] inter, double[,] vec)
        {
            int i, j, k, p, q;
            double di;

            p = 1;
            q = n;
            for (i = 1; i <= n; i++)
            {
                di = d[i];
                if (di != 1)
                    for (j = n1; j <= n2; j++) vec[i, j] *= di;
                k = inter[i];
                if (k > 0)
                    p++;
                else
                    if (k < 0) q--;
            }
            for (i = p - 1 + n - q; i >= 1; i--)
            {
                k = inter[i];
                if (k > 0)
                {
                    p--;
                    if (k != p) ichrow(n1, n2, k, p, vec);
                }
                else
                {
                    q++;
                    if (-k != q) ichrow(n1, n2, -k, q, vec);
                }
            }
        }

        private int reaqri(double[,] a, int n, double[] em, double[] val, double[,] vec)
        {
            int m1, i, i1, m, j, q, max, count;
            double w, shift, kappa, nu=0, mu=0, r, tol, s, machtol, elmax, t, delta, det;
            double[] tf = new double[n + 1];

            machtol = em[0] * em[1];
            tol = em[1] * em[2];
            max = (int)em[4];
            count = 0;
            elmax = 0.0;
            m = n;
            for (i = 1; i <= n; i++)
            {
                vec[i, i] = 1.0;
                for (j = i + 1; j <= n; j++) vec[i, j] = vec[j, i] = 0.0;
            }
            do
            {
                m1 = m - 1;
                i = m;
                do
                {
                    q = i;
                    i--;
                } while ((i >= 1) ? (Math.Abs(a[i + 1, i]) > tol ? true : false) : false);
                if (q > 1)
                    if (Math.Abs(a[q, q - 1]) > elmax) elmax = Math.Abs(a[q, q - 1]);
                if (q == m)
                {
                    val[m] = a[m, m];
                    m = m1;
                }
                else
                {
                    delta = a[m, m] - a[m1, m1];
                    det = a[m, m1] * a[m1, m];
                    if (Math.Abs(delta) < machtol)
                        s = Math.Sqrt(det);
                    else
                    {
                        w = 2.0 / delta;
                        s = w * w * det + 1.0;
                        s = (s <= 0.0) ? -delta * 0.5 : w * det / (Math.Sqrt(s) + 1.0);
                    }
                    if (q == m1)
                    {
                        val[m] = a[m, m] += s;
                        val[q] = a[q, q] -= s;
                        t = (Math.Abs(s) < machtol) ? (s + delta) / a[m, q] : a[q, m] / s;
                        r = Math.Sqrt(t * t + 1.0);
                        nu = 1.0 / r;
                        mu = -t * nu;
                        a[q, m] -= a[m, q];
                        rotrow(q + 2, n, q, m, a, mu, nu);
                        rotcol(1, q - 1, q, m, a, mu, nu);
                        rotcol(1, n, q, m, vec, mu, nu);
                        m -= 2;
                    }
                    else
                    {
                        count++;
                        if (count > max)
                        {
                            em[3] = elmax;
                            em[5] = count;
                            return m;
                        }
                        shift = a[m, m] + s;
                        if (Math.Abs(delta) < tol)
                        {
                            w = a[m1, m1] - s;
                            if (Math.Abs(w) < Math.Abs(shift)) shift = w;
                        }
                        a[q, q] -= shift;
                        for (i = q; i <= m1; i++)
                        {
                            i1 = i + 1;
                            a[i1, i1] -= shift;
                            kappa = Math.Sqrt(a[i, i] * a[i, i] + a[i1, i] * a[i1, i]);
                            if (i > q)
                            {
                                a[i, i - 1] = kappa * nu;
                                w = kappa * mu;
                            }
                            else
                                w = kappa;
                            mu = a[i, i] / kappa;
                            nu = a[i1, i] / kappa;
                            a[i, i] = w;
                            rotrow(i1, n, i, i1, a, mu, nu);
                            rotcol(1, i, i, i1, a, mu, nu);
                            a[i, i] += shift;
                            rotcol(1, n, i, i1, vec, mu, nu);
                        }
                        a[m, m1] = a[m, m] * nu;
                        a[m, m] = a[m, m] * mu + shift;
                    }
                }
            } while (m > 0);
            for (j = n; j >= 2; j--)
            {
                tf[j] = 1.0;
                t = a[j, j];
                for (i = j - 1; i >= 1; i--)
                {
                    delta = t - a[i, i];
                    tf[i] = matvec(i + 1, j, i, a, tf) /
                                ((Math.Abs(delta) < machtol) ? machtol : delta);
                }
                for (i = 1; i <= n; i++) vec[i, j] = matvec(1, j, i, vec, tf);
            }
            em[3] = elmax;
            em[5] = count;
            return m;
        }

        private void reascl(double[,] a, int n, int n1, int n2)
        {
            int i, j;
            double s;

            for (j = n1; j <= n2; j++)
            {
                s = 0.0;
                for (i = 1; i <= n; i++)
                    if (Math.Abs(a[i, j]) > Math.Abs(s)) s = a[i, j];
                if (s != 0.0)
                    for (i = 1; i <= n; i++) a[i, j] /= s;
            }
        }

        private int reaeig3(double[,] a, int n, double[] em, double[] val, double[,] vec)
        {
            int i;

            int[] ind = new int[n + 1];
            int[] ind0 = new int[n + 1];
            double[] d = new double[n + 1];

            eqilbr(a, n, em, d, ind0);
            tfmreahes(a, n, em, ind);
            i = reaqri(a, n, em, val, vec);

            if (i == 0)
            {
                bakreahes2(a, n, 1, n, ind, vec);
                baklbr(n, 1, n, d, ind0, vec);
                reascl(vec, n, 1, n);
            }

            return i;
        }

        public double[] FindEigenValuesVectors(out double[,] vec, out int count)
        {
            double[] em = new double[6];
            double[] val = new double[M + 1];
            double[,] a = new double[M + 1, M + 1];

            em[0] = 1.0e-8;
            em[2] = 1.0e-7;
            em[4] = 256.0;
            vec = new double[M + 1, M + 1];

            for (int i = 0; i < M; i++)
                for (int j = 0; j < M; j++)
                    a[i + 1, j + 1] = L[i, j];

            count = reaeig3(a, M, em, val, vec);
            return val;
        }
    }
}