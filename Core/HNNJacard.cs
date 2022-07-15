using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace uQlustCore
{
    [Serializable]
    class HNNJacard:MinStateHash//RootHash
    {
        //Dictionary<int, List<int>>[] locBase = null;
        //protected int[] minStateDataBase = null;

        //protected int[,][] oTab = null;
        public HNNJacard(int binSizeG, HashCluster hk, ClusterOutput outp, HNNCInput opt):base(binSizeG,hk, outp, opt)
        {

        }
         public HNNJacard(int binSizeG, HNNCInput opt):base(binSizeG,opt)
        {

        }
        int JacardIndex(int threadNum, string [] keys, int[] index, int num)
        {
            for (int j = 0; j < dataBaseKeys.Length; j++)
            {
                double common = 0;
                double all = 0;
                string locString = caseBase[dataBaseKeys[j]];
                for (int i = 0; i < keys.Length; i++)
                {
                    int v = Convert.ToChar(keys[i]);
                    if (v == '2' && v == locString[i])
                        common++;
                    if (v == '2' || locString[i] == '2')
                        all++;
                }
                dist[threadNum][j] =100- (int)(common / all*100);
            }

            return dataBaseKeys.Length;

        }
      
        int JacardIndex(int threadNum, int[] keys, int[] index, int num)
        {
            int counter = 0;
            
            int[] locDist = dist[threadNum];
            
            for (int j = 0; j < keys.Length; j++)
            {
              
                int[] aux = oTab[j, keys[j]];
                if (aux!=null)
                {
                    counter++;
                    for (int i = 0; i <aux.Length; i++)
                        locDist[aux[i]]++;
                }                

            }
            for (int i = 0; i < locDist.Length; i++)
            {
                int x = minStateDataBase[i] + counter - locDist[i];               
                locDist[i] = 100-(int)((100.0 * locDist[i]) /x);
            }



            return dataBaseKeys.Length;

        }

        public override string CalcDist(int threadNum, int[] keys, int[] index, int num)
        {
            JacardIndex(threadNum, keys, index, num);

            //FastSort(dist[threadNum], index);
            //Array.Sort(dist[threadNum], index);
            //Array.Reverse(dist[threadNum],0,resSize);
//            Array.Reverse(index,0,resSize);
            string w = "";
            for (int i = 0; i < resSize; i++)
                //w += dataBaseKeys[index[i]]+"-"+dist[threadNum][i] + ":";
                w += dataBaseKeys[index[i]] + ":";
            w += dataBaseKeys[index[resSize]];

            return w;

        }
        public override string CalcDist(int threadNum, string [] keys, int[] index, int num)
        {

            JacardIndex(threadNum, keys, index, num);

            Array.Sort(dist[threadNum], index);
            //Array.Reverse(dist[threadNum]);
            //Array.Reverse(index);
            string w = "";
            for (int i = 0; i < resSize; i++)
                //w += dataBaseKeys[index[i]]+"-"+dist[threadNum][i] + ":";
                w += dataBaseKeys[index[i]] + ":";
            w += dataBaseKeys[index[resSize]];

            return w;

        }
    }
}
