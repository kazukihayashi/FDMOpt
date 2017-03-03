using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace ForceDensityMethod
{
    public class FDMobject
    {
        public FDMobject()
        {
            this.gp = new List<Point3d>();
            this.tp = new List<Curve>();
        }
        public double[] q;
        public List<Point3d> gp;
        public List<Curve> tp;
        public double[] cs;
        public double ObjFun;
        public double[] ConFun;
    }
}
