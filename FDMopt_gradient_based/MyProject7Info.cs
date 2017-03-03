using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace MyProject7
{
    public class MyProject7Info : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "FDM_gradient_based";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                //Return a 24x24 pixel bitmap to represent this GHA library.
                return null;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("0549a8eb-7d8f-48a7-97a6-f5cac3baa79c");
            }
        }

        public override string AuthorName
        {
            get
            {
                //Return a string identifying you or your company.
                return "";
            }
        }
        public override string AuthorContact
        {
            get
            {
                //Return a string representing your preferred contact details.
                return "";
            }
        }
    }
}
