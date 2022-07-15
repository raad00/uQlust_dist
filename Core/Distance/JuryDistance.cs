using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using uQlustCore;

namespace uQlustCore.Distance

{
    public class JuryDistance: HammingBase//,IDistance
    {

        public JuryDistance(DCDFile dcd, string alignFile, bool flag, string profileName, string refJuryProfile = null)        
            :base(dcd, alignFile, flag, profileName, refJuryProfile)        
        {
        }
        public JuryDistance(List<string> fileNames, string alignFile, bool flag, string profileName, string refJuryProfile = null) 
            :base(fileNames, alignFile, flag,  profileName,refJuryProfile)
        {

        }
        public JuryDistance(string dirName,string alignFile,bool flag,string profileName,string refJuryProfile=null) 
            :base(dirName,alignFile,flag,profileName,refJuryProfile) 
        {
          
        }
        public JuryDistance(string profilesFile, bool flag , string profileName, string refJuryProfile)
            :base(profilesFile,flag,profileName,refJuryProfile)
        {
        }
        public JuryDistance(Alignment al, bool reference) : base(al, reference)
        { }
        public void WriteStates()
        {
            StreamWriter wr = new System.IO.StreamWriter("symbols.txt");

            foreach(var item in lStates)
            {
                foreach (var iStates in states.Keys)
                    if (item.ContainsKey(iStates))
                        wr.Write(item[iStates] + " ");
                    else
                        wr.Write("0 ");
                wr.WriteLine();
            }

        }
      

        public override string ToString()
        {
            return "HAMMING";
        }
        public override bool SimilarityThreshold(float threshold, float dist)
        {
            if (dist < threshold)
                return true;
            return false;
        }


        public override int GetDistance(string refStructure, string modelStructure)
        {
            double lDist = 0;
            Dictionary <string,int> dicStates=new Dictionary<string,int>();

            if (!stateAlign.ContainsKey(refStructure) || !stateAlign.ContainsKey(modelStructure))
                return errorValue;

            List<byte> refL = stateAlign[refStructure];
            List<byte> modL = stateAlign[modelStructure];

         
            for (int j = 0; j < refL.Count; j++)
            {
                if (refL[j] == 0 || modL[j] == 0)
                {
                    lDist += 1;
                    continue;
                }
                if (weights != null && refL[j] != modL[j])
                {
                    
                    if (weights.ContainsKey(refL[j]))
                    {
                        if(weights[refL[j]].ContainsKey(modL[j]))
                            lDist += weights[refL[j]][modL[j]];
                    }

                }


            }
            return (int)(lDist * 100);

        }
         
    }
}
