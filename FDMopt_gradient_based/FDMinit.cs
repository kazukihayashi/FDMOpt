using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Special;
using Rhino.Geometry;
using ln = MathNet.Numerics.LinearAlgebra;
using Grasshopper.Kernel.Types;

namespace ForceDensityMethod
{
    public class FDMinit
    {
        public const double ee= 1;

        protected List<Point3d> gpinit;
        protected List<Curve> tpinit;
        protected List<Curve> spinit;
        protected List<int> lpinit;
        protected List<Vector3d> lvinit;
        protected double[] q;

        /// <summary>
        /// initial force density
        /// </summary>
        public static void Init(ref double[] q, ref List<Point3d> gpinit, ref List<Curve> tpinit, ref List<Curve> spinit, ref List<int> lpinit, ref List<Vector3d> lvinit)
        {
            int nk = gpinit.Count;//number of nodes
            int nm = tpinit.Count;//number of members

            double[,,] tt = new double[6, 6, nm];

            int[] istart = new int[nm];
            int[] iend = new int[nm];
            double[] dll = new double[nm];
            for (int i = 0; i < nm; i++)
            {
                istart[i] = Rhino.Collections.Point3dList.ClosestIndexInList(gpinit, tpinit[i].PointAtStart);
                iend[i] = Rhino.Collections.Point3dList.ClosestIndexInList(gpinit, tpinit[i].PointAtEnd);
                double dx = gpinit[iend[i]].X - gpinit[istart[i]].X;
                double dy = gpinit[iend[i]].Y - gpinit[istart[i]].Y;
                double dz = gpinit[iend[i]].Z - gpinit[istart[i]].Z;
                dll[i] = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                tt[0, 0, i] = dx / dll[i];
                tt[1, 0, i] = dy / dll[i];
                tt[2, 0, i] = dz / dll[i];
                if (Math.Abs(tt[0, 0, i]) >= 0.9) { tt[1, 1, i] = .1e1; }
                else { tt[0, 1, i] = .1e1; }
                tt[0, 2, i] = tt[1, 0, i] * tt[2, 1, i] - tt[2, 0, i] * tt[1, 1, i];
                tt[1, 2, i] = tt[2, 0, i] * tt[0, 1, i] - tt[0, 0, i] * tt[2, 1, i];
                tt[2, 2, i] = tt[0, 0, i] * tt[1, 1, i] - tt[1, 0, i] * tt[0, 1, i];
                double ddd = Math.Sqrt(tt[0, 2, i] * tt[0, 2, i] + tt[1, 2, i] * tt[1, 2, i] + tt[2, 2, i] * tt[2, 2, i]);
                for (int j = 0; j < 3; j++){ tt[j, 2, i] /= ddd; }
                tt[0, 1, i] = tt[1, 2, i] * tt[2, 0, i] - tt[2, 2, i] * tt[1, 0, i];
                tt[1, 1, i] = tt[2, 2, i] * tt[0, 0, i] - tt[0, 2, i] * tt[2, 0, i];
                tt[2, 1, i] = tt[0, 2, i] * tt[1, 0, i] - tt[1, 2, i] * tt[0, 0, i];
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++) { tt[j + 3, k + 3, i] = tt[j, k, i]; }
                }
            }


            int nsp = spinit.Count;

            ///fix0 is an index of supported & loaded points
            System.Collections.ArrayList fix0 = new System.Collections.ArrayList();
            //add supported points' indices into fix0
            for (int i = 0; i < nk; i++)
            {
                for (int j = 0; j < nsp; j++)
                {
                    double t = new double();
                    bool scparam = spinit[j].ClosestPoint(gpinit[i], out t);
                    Point3d scp = spinit[j].PointAt(t);
                    if ((gpinit[i] - scp).IsTiny())
                    {
                        fix0.Add(i);
                    }
                }
            }

            fix0.Sort();//maybe unnecessary??
            //convert to integer array "fix1"
            int[] fix1 = (int[])fix0.ToArray(typeof(int));
            //convert to HashSet(HashSet can remove overlapping values and sort automatically and quickly)
            HashSet<int> fix2 = new HashSet<int>(fix1);
            //convert to integer array "fix"
            int[] fix = new int[fix2.Count];
            //int[] fix = new int[fix2.Count];
            fix2.CopyTo(fix, 0);

            int nfix = fix.Length;
            int ifix = 0;
            int[] free = new int[nk - nfix];
            int nfree = 0;
            for (int i = 0; i < nk; i++)
            {
                if (ifix != nfix)
                {
                    if (i == fix[ifix])
                    {
                        ifix++;
                    }
                    else
                    {
                        free[nfree] = i;
                        nfree++;
                    }
                }

                else
                {
                    free[nfree] = i;
                    nfree++;
                }
            }

            int ifree = 0;
            int[,] ird = new int[3,nk];

            for(int i=0;i<nfix;i++)
            {
                ird[0, fix[i]] = -1;
                ird[1, fix[i]] = -1;
                ird[2, fix[i]] = -1;
            }

            for (int i = 0; i < nk; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if (ird[j,i] == 0)
                    {
                        ird[j, i] = ifree;
                        ifree++;
                    }
                    else
                    {
                        ird[j, i] = nfree * 3;
                    }
                }
            }

            int[,] ir = new int[6, nm];
            for (int i = 0; i < nm; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    ir[j, i] = ird[j, istart[i]];
                    ir[j + 3, i] = ird[j, iend[i]];
                }
            }

            double[,] tp = new double[nm, nfree * 3];
            for (int j = 0; j < 6; j++)
            {
                for (int i = 0; i < nm; i++)
                {
                    if(ir[j,i]<nfree*3)
                    {
                        tp[i, ir[j, i]] -= tt[j, 0, i];
                        tp[i, ir[j, i]] += tt[j, 3, i];
                    }
                }
            }
            for (int j = 0; j < nfree*3; j++)
            {
                for (int i = 0; i < nm; i++)
                {
                    { tp[i,j] /= dll[i]; }
                }
            }

            double[,] dk = new double[6, 6];
            dk[0, 0] = .1e1;
            dk[0, 3] = -.1e1;
            dk[3, 0] = -.1e1;
            dk[3, 3] = .1e1;

            double[,,] ddk = new double[6, 6, nm];
            for (int ii = 0; ii < nm; ii++)
            {
                double[,] ff1 = new double[6, 6];
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        for (int k = 0; k < 6; k++)
                        {
                            ff1[i, k] += tt[i, j, ii] * dk[j, k];
                        }
                    }
                }

                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        for (int k = 0; k < 6; k++)
                        {
                            ddk[i, k, ii] += (ff1[i, j] * tt[k, j, ii]);
                        }
                    }
                }
            }

            ///ddk=ddk*E*A/L, E=A=1.0
            for (int i = 0; i < nm; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    for (int k = 0; k < 6; k++)
                    {
                        ddk[j, k, i] /= dll[i];
                    }
                }
            }

            double[,] ak = new double[nfree * 3, nfree * 3];
            for (int j = 0; j < 6; j++)
            {
                for (int k = 0; k < 6; k++)
                {
                    for (int i = 0; i < nm; i++)
                    {
                        if (ir[k, i] < nfree * 3)
                        {
                            if (ir[j, i] < nfree * 3)
                            { ak[ir[j, i], ir[k, i]] += ddk[j, k, i]; }
                        }
                    }
                }
            }

            double[] pp = new double[nfree * 3];
            for (int i = 0; i < lpinit.Count; i++)
            {
                pp[ird[0, lpinit[i]]] = lvinit[i].X;
                pp[ird[1, lpinit[i]]] = lvinit[i].Y;
                pp[ird[2, lpinit[i]]] = lvinit[i].Z;
            }

            var ak2 = ln.Matrix<double>.Build.DenseOfArray(ak);
            var ak3 = ak2.SubMatrix(0, nfree * 3, 0, nfree * 3);
            var ak4 = ln.Matrix<double>.Build.Dense(nfree, nfree);
            var pp2 = ln.Vector<double>.Build.DenseOfArray(pp);
            var ws = ln.Vector<double>.Build.Dense(nfree * 3);

            ak4 = ak3.Inverse();
            ws = ak3.Inverse() * pp2;

            double[] st = new double[nm];
            for (int i = 0; i < nm; i++)
            {
                for (int j = 0; j < nfree * 3; j++)
                {
                    st[i] += tp[i, j] * ws[j];
                }
            }

            double[] af = new double[nm];
            ///af=st*E*A, E=1.0, A=1.0
            for (int i = 0; i < nm; i++)
            {
                af[i] = st[i]*1.0*1.0;
                q[i] = af[i] / dll[i];
            }
        }

        public static void Init2D(ref double[] q, ref List<Point3d> gpinit, ref List<Curve> tpinit, ref List<Curve> spinit, ref List<int> lpinit, ref List<Vector3d> lvinit)
        {
            int nk = gpinit.Count;//number of nodes
            int nm = tpinit.Count;//number of members

            double[,,] tt = new double[6, 6, nm];

            int[] istart = new int[nm];
            int[] iend = new int[nm];
            double[] dll = new double[nm];
            for (int i = 0; i < nm; i++)
            {
                istart[i] = Rhino.Collections.Point3dList.ClosestIndexInList(gpinit, tpinit[i].PointAtStart);
                iend[i] = Rhino.Collections.Point3dList.ClosestIndexInList(gpinit, tpinit[i].PointAtEnd);
                double dx = gpinit[iend[i]].X - gpinit[istart[i]].X;
                double dy = gpinit[iend[i]].Y - gpinit[istart[i]].Y;
                double dz = gpinit[iend[i]].Z - gpinit[istart[i]].Z;
                dll[i] = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                tt[0, 0, i] = dx / dll[i];
                tt[1, 0, i] = dy / dll[i];
                tt[2, 0, i] = dz / dll[i];
                if (Math.Abs(tt[0, 0, i]) >= 0.9) { tt[1, 1, i] = .1e1; }
                else { tt[0, 1, i] = .1e1; }
                tt[0, 2, i] = tt[1, 0, i] * tt[2, 1, i] - tt[2, 0, i] * tt[1, 1, i];
                tt[1, 2, i] = tt[2, 0, i] * tt[0, 1, i] - tt[0, 0, i] * tt[2, 1, i];
                tt[2, 2, i] = tt[0, 0, i] * tt[1, 1, i] - tt[1, 0, i] * tt[0, 1, i];
                double ddd = Math.Sqrt(tt[0, 2, i] * tt[0, 2, i] + tt[1, 2, i] * tt[1, 2, i] + tt[2, 2, i] * tt[2, 2, i]);
                for (int j = 0; j < 3; j++) { tt[j, 2, i] /= ddd; }
                tt[0, 1, i] = tt[1, 2, i] * tt[2, 0, i] - tt[2, 2, i] * tt[1, 0, i];
                tt[1, 1, i] = tt[2, 2, i] * tt[0, 0, i] - tt[0, 2, i] * tt[2, 0, i];
                tt[2, 1, i] = tt[0, 2, i] * tt[1, 0, i] - tt[1, 2, i] * tt[0, 0, i];
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++) { tt[j + 3, k + 3, i] = tt[j, k, i]; }
                }
            }


            int nsp = spinit.Count;

            ///fix0 is an index of supported & loaded points
            System.Collections.ArrayList fix0 = new System.Collections.ArrayList();
            //add supported points' indices into fix0
            for (int i = 0; i < nk; i++)
            {
                for (int j = 0; j < nsp; j++)
                {
                    double t = new double();
                    bool scparam = spinit[j].ClosestPoint(gpinit[i], out t);
                    Point3d scp = spinit[j].PointAt(t);
                    if ((gpinit[i] - scp).IsTiny())
                    {
                        fix0.Add(i);
                    }
                }
            }

            fix0.Sort();//maybe unnecessary??
            //convert to integer array "fix1"
            int[] fix1 = (int[])fix0.ToArray(typeof(int));
            //convert to HashSet(HashSet can remove overlapping values and sort automatically and quickly)
            HashSet<int> fix2 = new HashSet<int>(fix1);
            //convert to integer array "fix"
            int[] fix = new int[fix2.Count];
            //int[] fix = new int[fix2.Count];
            fix2.CopyTo(fix, 0);

            int nfix = fix.Length;
            int ifix = 0;
            int[] free = new int[nk - nfix];
            int nfree = 0;
            for (int i = 0; i < nk; i++)
            {
                if (ifix != nfix)
                {
                    if (i == fix[ifix])
                    {
                        ifix++;
                    }
                    else
                    {
                        free[nfree] = i;
                        nfree++;
                    }
                }

                else
                {
                    free[nfree] = i;
                    nfree++;
                }
            }

            int ifree = 0;
            int[,] ird = new int[3, nk];

            for (int i = 0; i < nk; i++)
            {
                ird[2, i] = 1;
            }

            for (int i = 0; i < nfix; i++)
            {
                ird[0, fix[i]] = -1;
                ird[1, fix[i]] = -1;
            }

            for (int i = 0; i < nk; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if (ird[j, i] == 0)
                    {
                        ird[j, i] = ifree;
                        ifree++;
                    }
                    else
                    {
                        ird[j, i] = nfree * 2;
                    }
                }
            }

            int[,] ir = new int[6, nm];
            for (int i = 0; i < nm; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    ir[j, i] = ird[j, istart[i]];
                    ir[j + 3, i] = ird[j, iend[i]];
                }
            }

            double[,] tp = new double[nm, nfree * 2];
            for (int j = 0; j < 6; j++)
            {
                for (int i = 0; i < nm; i++)
                {
                    if (ir[j, i] < nfree * 2)
                    {
                        tp[i, ir[j, i]] -= tt[j, 0, i];
                        tp[i, ir[j, i]] += tt[j, 3, i];
                    }
                }
            }
            for (int j = 0; j < nfree * 2; j++)
            {
                for (int i = 0; i < nm; i++)
                {
                    { tp[i, j] /= dll[i]; }
                }
            }

            double[,] dk = new double[6, 6];
            dk[0, 0] = .1e1;
            dk[0, 3] = -.1e1;
            dk[3, 0] = -.1e1;
            dk[3, 3] = .1e1;

            double[,,] ddk = new double[6, 6, nm];
            for (int ii = 0; ii < nm; ii++)
            {
                double[,] ff1 = new double[6, 6];
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        for (int k = 0; k < 6; k++)
                        {
                            ff1[i, k] += tt[i, j, ii] * dk[j, k];
                        }
                    }
                }

                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        for (int k = 0; k < 6; k++)
                        {
                            ddk[i, k, ii] += (ff1[i, j] * tt[k, j, ii]);
                        }
                    }
                }
            }

            ///ddk=ddk*E*A/L, E=A=1.0
            for (int i = 0; i < nm; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    for (int k = 0; k < 6; k++)
                    {
                        ddk[j, k, i] /= dll[i];
                    }
                }
            }

            double[,] ak = new double[nfree * 2, nfree * 2];
            for (int j = 0; j < 6; j++)
            {
                for (int k = 0; k < 6; k++)
                {
                    for (int i = 0; i < nm; i++)
                    {
                        if (ir[k, i] < nfree * 2)
                        {
                            if (ir[j, i] < nfree * 2)
                            { ak[ir[j, i], ir[k, i]] += ddk[j, k, i]; }
                        }
                    }
                }
            }

            double[] pp = new double[nfree * 2];
            for (int i = 0; i < lpinit.Count; i++)
            {
                pp[ird[0, lpinit[i]]] = lvinit[i].X;
                pp[ird[1, lpinit[i]]] = lvinit[i].Y;
            }

            var ak2 = ln.Matrix<double>.Build.DenseOfArray(ak);
            var pp2 = ln.Vector<double>.Build.DenseOfArray(pp);
            var ws = ln.Vector<double>.Build.Dense(nfree * 2);

            ws = ak2.Inverse() * pp2;

            double[] st = new double[nm];
            for (int i = 0; i < nm; i++)
            {
                for (int j = 0; j < nfree * 2; j++)
                {
                    st[i] += tp[i, j] * ws[j];
                }
            }

            double[] af = new double[nm];
            ///af=st*E*A, E=1.0, A=1.0
            for (int i = 0; i < nm; i++)
            {
                af[i] = st[i] * 1.0 * 1.0;
                q[i] = af[i] / dll[i];
            }
        }
    }
}