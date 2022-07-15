using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace uQlustCore
{
    public class RetrivalInput:BaseCInput
    {
        [Description("Directory where structures for retrival are stored")]
        public string baseDir;
        [Description("Directory where structures to which similar sturctures must be found")]
        public string retrivalDir;
        [Description("Profile name")]
        public string profileName;
        [Description("Number of retrived structures")]
        public int numToRetrive = 10;
        [Description("Distance")]
        public DistanceMeasures measure=DistanceMeasures.RMSD;
        public PDB.PDBMODE atoms = PDB.PDBMODE.ONLY_CA;
    }
}
