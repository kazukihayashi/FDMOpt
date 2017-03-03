using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Forms;
using System.Threading.Tasks;
using System.ComponentModel;
using System.Diagnostics;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Special;
using Rhino.Geometry;
using NLoptNet;
using ln = MathNet.Numerics.LinearAlgebra;

namespace ForceDensityMethod
{
    public class CustomComponentAttributes : Grasshopper.Kernel.Attributes.GH_ComponentAttributes
    {
        // custom attribute to override double click mouse event on component and open a WPF window

        public CustomComponentAttributes(IGH_Component component)
            : base(component)
        {
            GhComponent = (FDMmain)component;
        }

        FDMmain GhComponent;

        [STAThread]
        public override Grasshopper.GUI.Canvas.GH_ObjectResponse RespondToMouseDoubleClick(Grasshopper.GUI.Canvas.GH_Canvas sender, Grasshopper.GUI.GH_CanvasMouseEvent e)
        {
            GhComponent.iterating = 1;
            Grasshopper.Instances.ActiveCanvas.Document.NewSolution(true);
            GhComponent.iterating = 2;
            return base.RespondToMouseDoubleClick(sender, e);
        }
    }
}
