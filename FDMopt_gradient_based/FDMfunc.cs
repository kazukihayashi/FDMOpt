using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Special;
using Rhino.Geometry;
using ln = MathNet.Numerics.LinearAlgebra;
using Grasshopper.Kernel.Types;
using NLoptNet;

namespace ForceDensityMethod
{
    public class FDMfunc
    {
        public const double ee = 1;

        /// <summary>
        /// free nodal coordinates as a function of force density method 
        /// </summary>

        public static void Sens(ref double[] q, ref List<Point3d> gp, ref List<Curve> tp, ref List<int> lp, ref double[] cs, ref double compliance, ref double[] rfload, ref int[] fix, ref int[] free, ref ln.Vector<double> Pfix, ref ln.Matrix<double> C, ref int[] istart, ref int[] iend, ref double[] scompliance, ref List<double[]> srfload)
        {
            int nk = gp.Count;
            int nm = tp.Count;
            int nfree = free.Length;
            int nfix = fix.Length;

            ln.Matrix<double> Ct = C.Transpose();

            double[] x = new double[nk];
            double[] y = new double[nk];
            double[] z = new double[nk];

            for (int i = 0; i < nk; i++)
            {
                x[i] = gp[i].X;
                y[i] = gp[i].Y;
                z[i] = gp[i].Z;
            }

            ln.Matrix<double> qd = ln.Matrix<double>.Build.DenseOfDiagonalArray(q);
            List<ln.Matrix<double>> sqd = new List<ln.Matrix<double>>(nm);
            for (int i = 0; i < nm; i++)
            {
                var sq = ln.Vector<double>.Build.Dense(nm);
                sq[i] = 1.0;
                sqd.Add(ln.Matrix<double>.Build.DenseOfDiagonalVector(sq));
            }

            ln.Matrix<double> Q0 = Ct * qd * C;
            var sQ0 = new List<ln.Matrix<double>>(nm);
            for (int i = 0; i < nm; i++)
            {
                sQ0.Add(Ct * sqd[i] * C);
            }

            ///Qfree & sQfree
            var Qfree0 = Q0.Clone();
            for (int i = nfix - 1; i >= 0; i--)
            {
                Qfree0 = Qfree0.RemoveColumn(fix[i]).RemoveRow(fix[i]);
            }
            var sQfree0 = new List<ln.Matrix<double>>(nm);
            foreach (ln.Matrix<double> m in sQ0)
            {
                sQfree0.Add(m);
            }
            for (int j = 0; j < nm; j++)
            {
                for (int i = nfix - 1; i >= 0; i--)
                {
                    sQfree0[j] = sQfree0[j].RemoveColumn(fix[i]).RemoveRow(fix[i]);
                }
            }

            ///Qfix & sQfix
            var Qfix0 = Q0.Clone();
            for (int i = nfree - 1; i >= 0; i--)
            {
                Qfix0 = Qfix0.RemoveColumn(free[i]).RemoveRow(free[i]);
            }
            var sQfix0 = new List<ln.Matrix<double>>(nm);
            foreach (ln.Matrix<double> m in sQ0)
            {
                sQfix0.Add(m);
            }
            for (int j = 0; j < nm; j++)
            {
                for (int i = nfree - 1; i >= 0; i--)
                {
                    sQfix0[j] = sQfix0[j].RemoveColumn(free[i]).RemoveRow(free[i]);
                }
            }

            ///Qlink & sQlink
            var Qlink0 = Q0.Clone();
            for (int i = nfree - 1; i >= 0; i--)
            {
                Qlink0 = Qlink0.RemoveColumn(free[i]);
            }
            for (int i = nfix - 1; i >= 0; i--)
            {
                Qlink0 = Qlink0.RemoveRow(fix[i]);
            }
            var sQlink0 = new List<ln.Matrix<double>>(nm);
            foreach (ln.Matrix<double> m in sQ0)
            {
                sQlink0.Add(m);
            }
            for (int j = 0; j < nm; j++)
            {
                for (int i = nfree - 1; i >= 0; i--)
                {
                    sQlink0[j] = sQlink0[j].RemoveColumn(free[i]);
                }
                for (int i = nfix - 1; i >= 0; i--)
                {
                    sQlink0[j] = sQlink0[j].RemoveRow(fix[i]);
                }
            }


            var Qfree = ln.Matrix<double>.Build.Dense(nfree * 3, nfree * 3);
            var Qfix = ln.Matrix<double>.Build.Dense(nfix * 3, nfix * 3);
            var Qlink = ln.Matrix<double>.Build.Dense(nfree * 3, nfix * 3);

            for (int i = 0; i < 3; i++)
            {
                Qfree.SetSubMatrix(nfree * i, nfree * i, Qfree0);
                Qfix.SetSubMatrix(nfix * i, nfix * i, Qfix0);
                Qlink.SetSubMatrix(nfree * i, nfix * i, Qlink0);
            }

            var sQfree = new List<ln.Matrix<double>>(nm);
            var sQfix = new List<ln.Matrix<double>>(nm);
            var sQlink = new List<ln.Matrix<double>>(nm);
            for (int i = 0; i < nm; i++)
            {
                sQfree.Add(ln.Matrix<double>.Build.Dense(nfree * 3, nfree * 3));
                sQfix.Add(ln.Matrix<double>.Build.Dense(nfix * 3, nfix * 3));
                sQlink.Add(ln.Matrix<double>.Build.Dense(nfree * 3, nfix * 3));
            }

            for (int j = 0; j < nm; j++)
            {
                for (int i = 0; i < 3; i++)
                {
                    sQfree[j].SetSubMatrix(nfree * i, nfree * i, sQfree0[j]);
                    sQfix[j].SetSubMatrix(nfix * i, nfix * i, sQfix0[j]);
                    sQlink[j].SetSubMatrix(nfree * i, nfix * i, sQlink0[j]);
                }
            }

            var Xfree = ln.Vector<double>.Build.Dense(nfree * 3);
            var Xfix = ln.Vector<double>.Build.Dense(nfix * 3);
            for (int i = 0; i < nfix; i++)
            {
                Xfix[i] = gp[fix[i]].X;
                Xfix[nfix + i] = gp[fix[i]].Y;
                Xfix[nfix * 2 + i] = gp[fix[i]].Z;
            }

            Xfree = -(Qfree).Inverse() * Qlink * Xfix;

            var sQfreeXfree = ln.Matrix<double>.Build.Dense(nfree * 3, nm);
            for (int i = 0; i < nm; i++)
            {
                sQfreeXfree.SetColumn(i, sQfree[i] * Xfree);
            }

            var sQlinkXfix = ln.Matrix<double>.Build.Dense(nfree * 3, nm);
            for (int i = 0; i < nm; i++)
            {
                sQlinkXfix.SetColumn(i, sQlink[i] * Xfix);
            }

            var sXfree = ln.Matrix<double>.Build.Dense(nfree * 3, nm);
            sXfree = -Qfree.Inverse() * (sQfreeXfree + sQlinkXfix);

            ///return nodal coordinates
            for (int i = 0; i < nfree; i++)
            {
                x[free[i]] = Xfree[i];
                y[free[i]] = Xfree[nfree + i];
                z[free[i]] = Xfree[nfree * 2 + i];
            }
            double[,] sx = new double[nk,nm];
            double[,] sy = new double[nk,nm];
            double[,] sz = new double[nk,nm];
            for (int j = 0; j < nm; j++)
            {
                for (int i = 0; i < nfree; i++)
                {
                    sx[free[i], j] = sXfree[i, j];
                    sy[free[i], j] = sXfree[nfree + i, j];
                    sz[free[i], j] = sXfree[nfree * 2 + i, j];
                }
            }

            for (int i = 0; i < nfree; i++)
            {
                gp[free[i]] = new Point3d(x[free[i]], y[free[i]], z[free[i]]);
            }

            for (int i = 0; i < nm; i++)
            {
                tp[i] = new LineCurve(gp[istart[i]], gp[iend[i]]);
            }

            double[] length = new double[nm];
            double[] length2 = new double[nm];
            double[,] slength2 = new double[nm,nm];

            for (int i = 0; i < nm; i++)
            {
                length[i] = Math.Sqrt(Math.Pow(x[iend[i]] - x[istart[i]], 2) + Math.Pow(y[iend[i]] - y[istart[i]], 2) + Math.Pow(z[iend[i]] - z[istart[i]], 2));
                length2[i] = Math.Pow(length[i], 2);
            }

            for (int j = 0; j < nm; j++)
            {
                for (int i = 0; i < nm; i++)
                {
                    slength2[i, j] = 2 * ((x[iend[i]] - x[istart[i]]) * (sx[iend[i],j] - sx[istart[i],j])+ (y[iend[i]] - y[istart[i]]) * (sy[iend[i], j] - sy[istart[i], j])+ (z[iend[i]] - z[istart[i]]) * (sz[iend[i], j] - sz[istart[i], j]));
                }
            }

            var sQlinkXfree = ln.Matrix<double>.Build.Dense(nfix * 3, nm);
            var QlinksXfree = ln.Matrix<double>.Build.Dense(nfix * 3, nm);
            var sQfixXfix = ln.Matrix<double>.Build.Dense(nfix * 3, nm);
            for (int i = 0; i < nm; i++)
            {
                sQlinkXfree.SetColumn(i, sQlink[i].Transpose() * Xfree);
                QlinksXfree.SetColumn(i, Qlink.Transpose() * sXfree.Column(i));
                sQfixXfix.SetColumn(i, sQfix[i] * Xfix);
            }

            ///reaction force
            var rf = ln.Vector<double>.Build.Dense(nfix * 3);
            var srf = ln.Matrix<double>.Build.Dense(nfix * 3, nm);
            rf = Qlink.Transpose() * Xfree + (Qfix) * Xfix;
            srf = sQlinkXfree + QlinksXfree + sQfixXfix;

            int iload = 0;
            int nl = lp.Count;

            for (int i = 0; i < nfix; i++)
            {
                if (iload == nl)
                {
                    break;
                }
                if (fix[i] == lp[iload])
                {
                    rfload[0 + 3 * iload] = rf[i] - Pfix[i];
                    rfload[1 + 3 * iload] = rf[nfix + i] - Pfix[nfix + i];
                    rfload[2 + 3 * iload] = rf[nfix * 2 + i] - Pfix[nfix * 2 + i];
                    iload++;
                }
            }

            iload = 0;
            srfload = new List<double[]>(nl*3);
            for (int i = 0; i < nl*3; i++)
            {
                srfload.Add(new double[nm]);
            }
            for (int j = 0; j < nm; j++)
            {
                for (int i = 0; i < nfix; i++)
                {
                    if (iload == nl)
                    {
                        break;
                    }
                    if (fix[i] == lp[iload])
                    {
                        srfload[0 + 3 * iload][j] = srf[i, j];
                        srfload[1 + 3 * iload][j] = srf[nfix + i, j];
                        srfload[2 + 3 * iload][j] = srf[nfix * 2 + i, j];
                        iload++;
                    }
                }
            }

            ///cross-sectional area
            cs = new double[nm];
            for (int i = 0; i < nm; i++)
            {
                cs[i] = Math.Sqrt(Math.Pow(q[i], 2) + .1e-5) * length[i];
            }

            ///volume
            double vol = 0;
            for (int i = 0; i < nm; i++)
            {
                vol += Math.Sqrt(Math.Pow(q[i], 2) + .1e-5) * length2[i];
            }

            ///strain energy
            double se = 0;
            double[] sse = new double[nm];
            for (int i = 0; i < nm; i++)
            {
                se += .5e0 * length2[i] * Math.Sqrt(Math.Pow(q[i], 2) + .1e-5) / ee;
            }
            for (int j = 0; j < nm; j++)
            {
                for (int i = 0; i < nm; i++)
                {
                    sse[j] += .5e0 * slength2[i,j] * Math.Sqrt(Math.Pow(q[i], 2) + .1e-5) / ee;
                    sse[j] += .5e0 * length2[i] * q[i] / Math.Sqrt(Math.Pow(q[i], 2) + .1e-5) * ee;
                }
            }

            ///compliance
            compliance = se * .2e1;
            scompliance = new double[nm];
            for (int j = 0; j < nm; j++)
            {
                scompliance[j] = sse[j] * .2e1;
            }

        }
    }
}