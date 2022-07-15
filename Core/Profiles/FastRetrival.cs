using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace uQlustCore.Profiles
{
    class FastRetrival:HashCluster
    {
        public FastRetrival()
        {

        }
        public FastRetrival(string dirName, string alignFile, HashCInput input):base(dirName,alignFile,input)
        {
            
        }
        public void InitFastRetrival()
        {
           InitHashCluster();
        }
        void PrepareFastRetrival(List<string> allStr)
        {
            PrepareClustering(allStr);
        }
        public List<string> RunFastRetrival(List<string> dataBase,List<string> retrivalList,int nearest)
        {
            List<string> bestRes = new List<string>();
            int[] distArray; 
            int[] index; 
            StreamWriter wr = new StreamWriter("retrival.txt");
            List<string> toRemove = new List<string>();
            for (int i = 0; i < dataBase.Count; i++)
                if (!File.Exists("C:\\Projects\\S35\\" + dataBase[i]))
                    toRemove.Add(dataBase[i]);

            foreach (var item in toRemove)
            {
                dataBase.Remove(item);
                retrivalList.Remove(item);
            }
            distArray = new int[dataBase.Count];
            index = new int[dataBase.Count];
            foreach (var item in retrivalList)
            {
                if (!structToKey.ContainsKey(item))
                    continue;

                string toRetrivProfile = structToKey[item];
                for (int i = 0; i < index.Length; i++)
                    index[i] = i;

                int count = 0;   
                foreach (var vecBase in dataBase)
                {
                    string baseProfile = structToKey[vecBase];

                    int dist = 0;
                    for (int i = 0; i < baseProfile.Length; i++)
                        dist += Math.Abs(baseProfile[i] - toRetrivProfile[i]);

                    distArray[count++]=dist;
                }
              //  wr.WriteLine("Next");
                Array.Sort(distArray, index);
                int counter = 0;
                int c = 0;
                //while(distArray[c]==0)
                for(int i=0;counter<nearest;i++)
                {
                    if (dataBase[index[c]] == item)
                    {
                        c++;
                        continue;
                    }
                    bestRes.Add(dataBase[index[c]]);
                    wr.WriteLine(distArray[c]+" "+item + " " + dataBase[index[c]]);
                    counter++;
                    c++;
                }
                //wr.WriteLine(c);

            }
            wr.Close();
            return bestRes;
        }
    }
}
