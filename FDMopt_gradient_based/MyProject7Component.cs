using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Special;
using Rhino.Geometry;
using NLoptNet;
using LibOptimization;
using ln=MathNet.Numerics.LinearAlgebra;

namespace ForceDensityMethod
{
    public class FDMmain : GH_Component
    {
        /// By adding these 2 lines, we can open the console window
        /// ex) line 1: AllocConsole();  
        ///     line 2: WriteLine("put your message here")
        [System.Runtime.InteropServices.DllImport("kernel32.dll")]
        private static extern bool AllocConsole();                 

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>

        public int iterating = 0;

        public List<Point3d> gp = new List<Point3d>();
        public List<Curve> tp = new List<Curve>();
        public List<Curve> sp = new List<Curve>();
        public List<int> lp = new List<int>();
        public List<Vector3d> lv = new List<Vector3d>();
        public int seed = new int();

        private List<Point3d> gpinit = new List<Point3d>();
        private List<Curve> tpinit = new List<Curve>();
        private List<Curve> spinit = new List<Curve>();
        private List<int> lpinit = new List<int>();
        private List<Vector3d> lvinit = new List<Vector3d>();

        public double compliance;
        public double[] rfload;
        public double[] cs;
        private double[] q;
        public int[] fix;
        public int[] free;
        public int nfix;
        public int nfree;
        public ln.Vector<double> Pfix;
        public ln.Matrix<double> C;
        public int[] istart;
        public int[] iend;
        public double[] scompliance;
        public List<double[]> srfload;
        public int iConFun;

        public FDMmain()
            : base("ForceDensityMethod", "FDMopt",
                "Opt using FDM",
                "Math", "Script")
        {}

        public override void CreateAttributes()
        {
            base.m_attributes = new CustomComponentAttributes(this);
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("GeometryPoint", "P", "Nodes to determine the geometry of model", GH_ParamAccess.list);
            pManager.AddCurveParameter("Topology", "C", "Linear curves of model", GH_ParamAccess.list);
            pManager.AddCurveParameter("Support", "S", "Index of pin-supported nodes of model", GH_ParamAccess.list);
            pManager.AddIntegerParameter("LoadedPoint", "LP", "Index of nodal load points", GH_ParamAccess.list);
            pManager.AddVectorParameter("LoadVector", "LV", "Vector of nodal loads", GH_ParamAccess.list);
            pManager.AddIntegerParameter("RandomSeed", "I", "Randomize initial values of force density", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("GeometryPoint", "P", "Nodes to determine the geometry of model", GH_ParamAccess.list);
            pManager.AddCurveParameter("Topology", "C", "Linear curves of model", GH_ParamAccess.list);
            pManager.AddNumberParameter("ForceDensity", "FD", "List of force densities of all members", GH_ParamAccess.list);
            pManager.AddNumberParameter("Compliance", "C", "External work that nodal loads did", GH_ParamAccess.item);
            pManager.AddNumberParameter("CrossSection", "A", "List of cross-sectional areas of all members", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            if (iterating==0)
            {
                if (!DA.GetDataList(0, gpinit)) { return; }
                if (!DA.GetDataList(1, tpinit)) { return; }
                if (!DA.GetDataList(2, spinit)) { return; }
                if (!DA.GetDataList(3, lpinit)) { return; }
                if (!DA.GetDataList(4, lvinit)) { return; }

                ///make HARD copy where original and copy are independent
                //gp = new List<Point3d>();
                foreach (Point3d p in gpinit)
                {
                    gp.Add(new Point3d(p));
                }

                //tp = new List<Curve>();
                foreach (Curve c in tpinit)
                {
                    tp.Add(c);
                }

                //sp = new List<Curve>();
                foreach (Curve c in spinit)
                {
                    sp.Add(c);
                }

                //lp = new List<int>();
                foreach (int i in lpinit)
                {
                    lp.Add(i);
                }

                //lv = new List<Vector3d>();
                foreach (Vector3d v in lvinit)
                {
                    lv.Add(new Vector3d(v));
                }

                int nl = lp.Count;//number of nodal loads
                int nl2 = lv.Count;//number of loaded vectors
                if (nl != nl2)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The number of LP should be equal to that of LV");
                    return;
                }

                double[] q = new double[tpinit.Count];
                FDMinit.Init(ref q, ref gpinit, ref tpinit, ref spinit, ref lpinit, ref lvinit);
                for (int i = 0; i < gpinit.Count; i++)
                {
                    if (gpinit[i].Z != 0) { break; }
                    if (i == gpinit.Count - 1) { FDMinit.Init2D(ref q, ref gpinit, ref tpinit, ref spinit, ref lpinit, ref lvinit); }
                }
                FDMinput.Input(ref gp, ref tp, ref sp, ref lp, ref lv, ref fix, ref free, ref Pfix, ref C, ref istart, ref iend, ref rfload);

                DA.SetDataList(0, gpinit);
                DA.SetDataList(1, tpinit);
                DA.SetDataList(2, q);
                DA.SetData(3, compliance);
                DA.SetDataList(4, cs);
            }

            if (iterating==1)
            {
                if (!DA.GetData(5, ref seed)) { return; }

                double[] q = new double[tpinit.Count];
                FDMinit.Init(ref q, ref gpinit, ref tpinit, ref spinit, ref lpinit, ref lvinit);
                for (int i = 0; i < gpinit.Count; i++)
                {
                    if (gpinit[i].Z != 0) { break; }
                    if (i == gpinit.Count - 1) { FDMinit.Init2D(ref q, ref gpinit, ref tpinit, ref spinit, ref lpinit, ref lvinit); }
                }

                Opt(ref q, ref gp, ref tp, ref sp, ref lp, ref lv, ref cs, ref compliance, ref rfload, ref seed);

                DA.SetDataList(0, gp);
                DA.SetDataList(1, tp);
                DA.SetDataList(2, q);
                DA.SetData(3, compliance);
                DA.SetDataList(4, cs);
            }
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                return MyProject7.Properties.Resources.FDMicon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{540c717f-1d3d-4dae-967b-ad5b36ba1e38}"); }
        }


        /// <summary>
        /// NLopt script MUST be written here.
        /// If you try to call this script from another cs.file,
        /// maybe you should fail to load it.
        /// </summary>

        internal double SObjFun(double[] q, double[] scompliance)
        {
            FDMfunc.Sens(ref q, ref gp, ref tp, ref lp, ref cs, ref compliance, ref rfload, ref fix, ref free, ref Pfix, ref C, ref istart, ref iend, ref scompliance, ref srfload);
            AllocConsole();
            Console.WriteLine("----------------------------------------------");
            Console.WriteLine("ObjFun={0}", compliance);
            return compliance;
        }

        //List<Func<double[], double[], double>> SConFun = new List<Func<double[], double[], double>>();
        internal double sConFun0(double[] q, double[] srfload_elem)
        {
            FDMfunc.Sens(ref q, ref gp, ref tp, ref lp, ref cs, ref compliance, ref rfload, ref fix, ref free, ref Pfix, ref C, ref istart, ref iend, ref scompliance, ref srfload);
            srfload_elem = srfload[0];
            Console.WriteLine("ConFun0={0}", rfload[0]);
            Console.WriteLine("ConFun1={0}", rfload[1]);
            Console.WriteLine("ConFun2={0}", rfload[2]);
            return rfload[0];
        }

        internal double sConFun1(double[] q, double[] srfload_elem2)
        {
            FDMfunc.Sens(ref q, ref gp, ref tp, ref lp, ref cs, ref compliance, ref rfload, ref fix, ref free, ref Pfix, ref C, ref istart, ref iend, ref scompliance, ref srfload);
            srfload_elem2 = srfload[1];
            return rfload[1];
        }


        public void Opt(ref double[] q, ref List<Point3d> gp, ref List<Curve> tp, ref List<Curve> sp, ref List<int> lp, ref List<Vector3d> lv, ref double[] cs, ref double compliance, ref double[] rfload, ref int seed)
        {
            FDMinput.Input(ref gp, ref tp, ref sp, ref lp, ref lv, ref fix, ref free, ref Pfix, ref C, ref istart, ref iend, ref rfload);
            uint nm = (uint)tpinit.Count;
            //FDMfunc.Sens(ref q, ref gp, ref tp, ref lp, ref cs, ref compliance, ref rfload, ref fix, ref free, ref Pfix, ref C, ref istart, ref iend, ref scompliance, ref srfload);

            using (var solver = new NLoptSolver(NLoptAlgorithm.LD_SLSQP, nm, .1e-5, 100))
            {
                var qmin = new double[nm];
                for (int i = 0; i < nm; i++)
                {
                    qmin[i] = q[i]-.1e3;
                }
                var qmax = new double[nm];
                for (int i = 0; i < nm; i++)
                {
                    qmax[i] = q[i]+.1e3;
                }
                solver.SetLowerBounds(qmin);
                solver.SetUpperBounds(qmax);


                solver.SetMinObjective(SObjFun);

                solver.AddEqualZeroConstraint(sConFun0, .1e-5);
                solver.AddEqualZeroConstraint(sConFun1, .1e-5);

                Random Random = new Random(seed);
                for (int i = 0; i < nm; i++)
                { q[i] += (Random.NextDouble()-0.5)*10; }

                double? finalScore = 0;
                var result = solver.Optimize(q, out finalScore);

                FDMfunc.Sens(ref q, ref gp, ref tp, ref lp, ref cs, ref compliance, ref rfload, ref fix, ref free, ref Pfix, ref C, ref istart, ref iend, ref scompliance, ref srfload);
                Console.WriteLine("============================================");
                Console.WriteLine(result);
                for (int i = 0; i < nm; i++)
                {
                    Console.WriteLine("q[{0}]={1}", i, q[i]);
                }
                Console.WriteLine("compliance={0}", compliance);
                Console.WriteLine("rfload[0]={0}", rfload[0]);
                Console.WriteLine("rfload[1]={0}", rfload[1]);
            }
        }
    }
}

