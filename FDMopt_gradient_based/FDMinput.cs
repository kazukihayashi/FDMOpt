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
    public class FDMinput
    {
        public const double ee = 1;

        public static void Input(ref List<Point3d> gp, ref List<Curve> tp, ref List<Curve> sp, ref List<int> lp, ref List<Vector3d> lv, ref int[] fix, ref int[] free, ref ln.Vector<double> Pfix, ref ln.Matrix<double> C, ref int[] istart, ref int[] iend, ref double[] rfload)
        {
            int nk = gp.Count;//number of nodes
            int nm = tp.Count;//number of members

            C = ln.Matrix<double>.Build.Dense(nm, nk);
            istart = new int[nm];
            iend = new int[nm];

            for (int i = 0; i < nm; i++)
            {
                istart[i] = Rhino.Collections.Point3dList.ClosestIndexInList(gp, tp[i].PointAtStart);
                iend[i] = Rhino.Collections.Point3dList.ClosestIndexInList(gp, tp[i].PointAtEnd);
                C[i, istart[i]] = -1;
                C[i, iend[i]] = 1;
            }

            ///definition of support input
            ///supports are defined by curve objects
            int nsp = sp.Count;

            ///fix0 is an index of supported & loaded points
            System.Collections.ArrayList fix0 = new System.Collections.ArrayList();
            //add supported points' indices into fix0
            for (int i = 0; i < nk; i++)
            {
                for (int j = 0; j < nsp; j++)
                {
                    double t = new double();
                    bool scparam = sp[j].ClosestPoint(gp[i], out t);
                    Point3d scp = sp[j].PointAt(t);
                    if ((gp[i] - scp).IsTiny())
                    {
                        fix0.Add(i);
                    }
                }
            }

            //add loaded points' indices into fix0
            int nl = lp.Count;//number of nodal loads
            for (int i = 0; i < nl; i++)
            {
                fix0.Add(lp[i]);
            }

            fix0.Sort();//maybe unnecessary??
            //convert to integer array "fix1"
            int[] fix1 = (int[])fix0.ToArray(typeof(int));
            //convert to HashSet(HashSet can remove overlapping values and sort automatically and quickly)
            HashSet<int> fix2 = new HashSet<int>(fix1);
            //convert to integer array "fix"
            fix = new int[fix2.Count];
            //int[] fix = new int[fix2.Count];
            fix2.CopyTo(fix, 0);

            int nfix = fix.Length;
            int ifix = 0;
            free = new int[nk - nfix];
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

            Pfix = ln.Vector<double>.Build.Dense(nfix * 3);
            int iload = 0;

            for (int i = 0; i < nfix; i++)
            {
                if (iload == nl)
                {
                    break;
                }
                if (fix[i] == lp[iload])
                {
                    Pfix[i] = lv[iload].X;
                    Pfix[nfix + i] = lv[iload].Y;
                    Pfix[nfix * 2 + i] = lv[iload].Z;
                    iload++;
                }
            }

            rfload = new double[nl * 3];
        }
    }
}