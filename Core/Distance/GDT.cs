using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using uQlustCore.PDB;

namespace uQlustCore.Distance
{
    class GDT:MaxSub
    {
        List<int> segments = new List<int>() { 3, 5, 7 };
        public float Threshold = 3.5f;

         public GDT(DCDFile dcd, string alignFile, bool flag, string refJuryProfile = null)
            : base(dcd, alignFile, flag, refJuryProfile)
        {
            order = true;
            maxSimilarity = 100.0;
        }
        public GDT(string dirName, string alignFile, bool flag, string refJuryProfile = null)
            : base(dirName, alignFile, flag,refJuryProfile)
        {
            this.dirName = dirName;
            this.alignFile = alignFile;
            this.flag = flag;
            this.refJuryProfile = refJuryProfile;

            order = true;
            maxSimilarity = 100.0;
        }
        public GDT(List<string> fileNames, string alignFile, bool flag, string refJuryProfile = null)
            : base(fileNames, alignFile, flag, refJuryProfile)
        {
            
            order = true;
            maxSimilarity = 100.0;
        }
        public override string ToString()
        {
            return "GDT";
        }
        public override int GetDistance(string refStructure, string modelStructure)
        {
            KeyValuePair<List<int>,float[,]> bestpair=default(KeyValuePair<List<int>,float[,]>);
            if (!pdbs.molDic.ContainsKey(refStructure) || !pdbs.molDic.ContainsKey(modelStructure))
                return errorValue;          

            posMOL locPosMol = Optimization.PrepareData(pdbs.molDic[refStructure], pdbs.molDic[modelStructure]);
            
            foreach(var item in segments)
            {
                if (locPosMol.posmol1.GetLength(0) < item)
                    continue;
                    KeyValuePair<List<int>, float[,]> seg = FindLongestSegment(item, Threshold, locPosMol.posmol1, locPosMol.posmol2);
                    if (seg.Key != null && (bestpair.Equals(default(KeyValuePair<List<int>, float[,]>)) || bestpair.Key.Count < seg.Key.Count))
                        bestpair = seg;
            }
          
                if (!bestpair.Equals(default(KeyValuePair<List<int>, float[,]>)))
                    if (bestpair.Key == null)
                        return 0;
                    else
                        return bestpair.Key.Count * 100 / pdbs.molDic[refStructure].molLength;
                    
            
            return 0;
          
        }
    }
}
