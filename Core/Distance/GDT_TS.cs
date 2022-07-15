using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using uQlustCore.PDB;

namespace uQlustCore.Distance
{
    public class GDT_TS : Rmsd
    {
        List<float> distances = new List<float>() { 1.0f, 2.0f, 4.0f, 8.0f };
        GDT gdt = null;
        public GDT_TS(DCDFile dcd, string alignFile, bool flag, string refJuryProfile = null)
            : base(dcd, alignFile, flag, PDBMODE.ONLY_CA, refJuryProfile)
        {
            order = true;
            maxSimilarity = 100.0;
            gdt = new GDT(dcd, alignFile, flag, refJuryProfile);
            
        }
        public GDT_TS(string dirName, string alignFile, bool flag, string refJuryProfile = null)
            : base(dirName, alignFile, flag, PDBMODE.ONLY_CA, refJuryProfile)
        {
            this.dirName = dirName;
            this.alignFile = alignFile;
            this.flag = flag;
            this.refJuryProfile = refJuryProfile;

            order = true;
            maxSimilarity = 100.0;
            gdt = new GDT(dirName, alignFile, flag, refJuryProfile);
        }
        public GDT_TS(List<string> fileNames, string alignFile, bool flag, string refJuryProfile = null)
            : base(fileNames, alignFile, flag, PDBMODE.ONLY_CA, refJuryProfile)
        {
            
            order = true;
            maxSimilarity = 100.0;
            gdt = new GDT(fileNames, alignFile, flag, refJuryProfile);
        }   
        public override void InitMeasure()
        {
            gdt.InitMeasure();
            structNames = gdt.structNames;
        }
        public override int GetDistance(string refStructure, string modelStructure)
        {
            List<int> res=new List<int>();
            int final = 0;
            foreach (var item in distances)
            {
                gdt.Threshold = item;
               
                int ww = gdt.GetDistance(refStructure, modelStructure);
                if(ww!=errorValue)              
                    res.Add(ww);
            }
            if (res.Count > 0)
            {
                foreach (var item in res)
                    final += item;

                final = 100-(int)(((float)final) / res.Count);
            }
            else final = errorValue;

            return final;
        }


    }
}
