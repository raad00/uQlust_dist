using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace uQlustCore
{
    [Serializable]
    public class HNNCInput : BaseCInput
    {
        [Description("Clustering method (Rpart or Hash")]
        public bool Rpart = false;
        [Description("Labels for items in train file")]
        public string labelsFile = "";
        [Description("Test file")]
        public string testFile= "";
        [Description("Binary file name")]
        public string binaryfile = "";
        [Description("Distance calculation method to be used (if true mSim else mLSH-c5")]
        public bool mSim = true;

        public int retrivalSize = 20;
    }
}
