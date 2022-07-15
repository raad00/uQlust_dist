using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace uQlustCore.Distance
{
    class Tanimoto : HammingBase
    {
        public Tanimoto(string dirName, string alignFile, bool flag, string profileName) : base(dirName, alignFile, flag, profileName)
        {
        }
        public Tanimoto(DCDFile dcd, string alignFile, bool flag, string profileName, string refJuryProfile = null)        
            :base(dcd, alignFile, flag, profileName, refJuryProfile)        
        {
        }
        public Tanimoto(List<string> fileNames, string alignFile, bool flag, string profileName, string refJuryProfile = null) 
            :base(fileNames, alignFile, flag,  profileName,refJuryProfile)
        {

        }
        public Tanimoto(string dirName, string alignFile, bool flag, string profileName, string refJuryProfile = null) 
            :base(dirName,alignFile,flag,profileName,refJuryProfile) 
        {

        }
        public Tanimoto(string profilesFile, bool flag, string profileName, string refJuryProfile)
            :base(profilesFile,flag,profileName,refJuryProfile)
        {
        }
        public Tanimoto(Alignment al, bool reference) : base(al, reference)
        { }
        public override int GetDistance(string refStructure, string modelStructure)
        {
            int dist = 0;
            if (!stateAlign.ContainsKey(refStructure))
                throw new Exception("Structure: " + refStructure + " does not exists in the available list of structures");

            if (!stateAlign.ContainsKey(modelStructure))
                throw new Exception("Structure: " + modelStructure + " does not exists in the available list of structures");

            List<byte> mod1 = stateAlign[refStructure];
            List<byte> mod2 = stateAlign[modelStructure];
            int common = 0;
            int all = 0;
            for (int j = 0; j < stateAlign[refStructure].Count; j++)
            {
                if (mod1[j] == 1)
                {
                    if (mod1[j] == mod2[j])
                        common++;
                    all++;
                }
                if (mod2[j] == 1)
                    all++;

            }
            dist=1000-(int)Math.Floor((double)(common) / (all - common) * 1000);
            return dist;
        }
    }
}
