using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using uQlustCore;

namespace uQlustCore.Distance
{
    public class HammingDistance:HammingBase
    {
        public HammingDistance(string dirName, string alignFile, bool flag, string profileName):
            base(dirName, alignFile, flag, profileName)
        {
            this.dirName = dirName;
            this.alignFile = alignFile;
            this.flag = flag;
            this.profileName = profileName;
        }
        public HammingDistance(Alignment al, bool reference) : base(al, reference)
        { }
/*        public override void InitMeasure()
        {
            if (fileNames == null)
                InitMeasure(dirName, alignFile, flag, profileName);
            else
                InitMeasure(fileNames, alignFile, flag, profileName);

        }*/

      /*  public void InitMeasure(string dirName, string alignFile, bool flag,string profileName)            
        {
            base.InitMeasure(dirName, alignFile, flag, profileName);
        }
        public void InitMeasure(List<string> fileNames, string alignFile, bool flag,string profileName)
            
        {
            base.InitMeasure(fileNames, alignFile, flag, profileName);
        }*/
        public override int GetDistance(string refStructure, string modelStructure)
        {
            int dist = 0;
            if(!stateAlign.ContainsKey(refStructure))
                    throw new Exception("Structure: "+refStructure+" does not exists in the available list of structures");

            if(!stateAlign.ContainsKey(modelStructure))
                    throw new Exception("Structure: "+modelStructure+" does not exists in the available list of structures");

            List<byte> mod1 = stateAlign[refStructure];
            List<byte> mod2 = stateAlign[modelStructure];
            for (int j = 0; j < stateAlign[refStructure].Count; j++)
            {
                if(mod1[j]!=mod2[j])// || mod1[j]==0)
                    dist++;
            }

            return dist;
        }
        public override string ToString()
        {
            return "HAMMING";
        }
/*        public override void CalcDistMatrix(List<string> structNames)
        {
            double[,] distTab = new double[structNames.Count, structNames.Count];


            Dictionary<string, List<int>>[] columns = MakeColumnsLists(structNames);

            int s = 0;
            hashIndex = new Dictionary<string, int>(structNames.Count);
            foreach (var item in structNames)
            {
                if (!stateAlign.ContainsKey(item))
                    throw new Exception("There is no alignment for the structure: " + item);
                if (!hashIndex.ContainsKey(item))
                    hashIndex.Add(item, s++);
            }

            for (int i = 0; i < structNames.Count; i++)
            {
                for (int j = 0; j < columns.Length; j++)
                {
                    string state = stateAlign[structNames[i]][j];
                    if (!weights.Keys.Contains(state))
                        continue;
                    foreach (var item in weights.Keys)
                    {
                        if (item == state)
                            continue;

                        if (columns[j].ContainsKey(item) && weights[state].ContainsKey(item))
                        {
                            List<int> commonStruct = columns[j][item];

                            for (int n = 0; n < commonStruct.Count; n++)
                                if (i <= commonStruct[n])
                                    distTab[i, commonStruct[n]] ++;
                                else
                                    break;
                        }
                    }
                }
            }


            distanceMatrix = new int[structNames.Count * (structNames.Count + 1) / 2];
            for (int i = 0; i < structNames.Count; i++)
                for (int j = 0; j < structNames.Count; j++)
                    if (hashIndex[structNames[i]] <= hashIndex[structNames[j]])
                    {

                        int dd = FindIndex(hashIndex[structNames[i]]);
                        distanceMatrix[dd + hashIndex[structNames[j]] - hashIndex[structNames[i]]] = (int)(distTab[hashIndex[structNames[i]], hashIndex[structNames[j]]] * 100);
                    }


        }*/

    }

}
