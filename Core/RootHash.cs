using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.IO;

namespace uQlustCore
{

    [Serializable]
    class RootHash :HNN
    {
        protected int binSizeG = 0;
        [NonSerialized]
        protected ManualResetEvent[] resetEvents;
        [NonSerialized]
        protected Dictionary<string, List<int>>[] hashDataBase = null;
        [NonSerialized]
        protected Dictionary<string, int[]>[] hashDataBaseTab = null;
        protected string[] dataBaseKeys = null;
        [NonSerialized]
        protected Dictionary<string, string> results = new Dictionary<string, string>();
        [NonSerialized]
        protected int[][] dist;
        public int resSize = 20;
        string baseFileName = "";
       
        public RootHash(int binSizeG, HashCluster hk, ClusterOutput outp, HNNCInput opt):base(hk, outp, opt)
        {
            this.binSizeG = binSizeG;
            baseFileName = opt.binaryfile;
 
            // number of bits to store the universe
        }
        public RootHash(int binSizeG,  HNNCInput opt) : base(opt)
        {
            this.binSizeG = binSizeG;
            baseFileName = opt.binaryfile;
  
            // number of bits to store the universe
        }
        protected List<string> GetKeyHashes(string item)
        {
            List<string> aux = new List<string>();
            for (int i = 0, counter = 0; i < item.Length; i += binSizeG, counter++)
            {
                string keyHash = "";
                if (i + binSizeG > item.Length)
                {
                    keyHash = item.Substring(i, item.Length - i);
                }
                else
                    keyHash = item.Substring(i, binSizeG);
                aux.Add(keyHash);
            }

            return aux;
        }
        public virtual int[] GetKeyHashesTabInt(string item)
        {
            int[] aux = new int[(int)Math.Ceiling(((double)item.Length) / binSizeG)];

            for (int i = 0, counter = 0; i < item.Length; i += binSizeG, counter++)
            {
                aux[counter] = Convert.ToInt32(item[i] - '0');
            }

            return aux;
        }
        protected virtual string[] GetKeyHashesTab(string item)
        {
            string[] aux = new string [(int)Math.Ceiling(((double)item.Length)/binSizeG)];
            for (int i = 0, counter = 0; i < item.Length; i += binSizeG, counter++)
            {
                string keyHash = "";
                if (i + binSizeG > item.Length)
                {
                    keyHash = item.Substring(i, item.Length - i);
                }
                else
                    keyHash = item.Substring(i, binSizeG);
                aux[counter]=keyHash;
            }

            return aux;
        }



        public virtual void CreateBase(Dictionary<string, string> dataBase)
        {

            dataBaseKeys = dataBase.Keys.ToArray();

            double x = dataBase[dataBaseKeys[0]].Length;
            x /= binSizeG;
            int tabSize = (int)x;

            if (x != Math.Round(x))
                tabSize++;

            hashDataBase = new Dictionary<string, List<int>>[tabSize];
            for (int i = 0; i < hashDataBase.Length; i++)
                hashDataBase[i] = new Dictionary<string, List<int>>();


            for (int j = 0; j < dataBaseKeys.Length; j++)
            {                
                List<string> aux = GetKeyHashes(dataBase[dataBaseKeys[j]]);
                for (int i = 0; i < aux.Count; i++)
                {
                    Dictionary<string, List<int>> locHash = hashDataBase[i];
                    if (locHash.ContainsKey(aux[i]))
                        locHash[aux[i]].Add(j);
                    else
                    {
                        List<int> xaux = new List<int>();
                        xaux.Add(j);
                        locHash.Add(aux[i], xaux);
                    }
                }
            }

            hashDataBaseTab = new Dictionary<string, int[]>[hashDataBase.Length];
            for(int i=0;i<hashDataBase.Length;i++)
            {
                hashDataBaseTab[i] = new Dictionary<string, int[]>();
                foreach(var item in hashDataBase[i])
                {
                    hashDataBaseTab[i].Add(item.Key, item.Value.ToArray());
                }
            }
        }
        public override void Preprocessing(bool load)
        {
            if (!load)
            {
                caseBase = new Dictionary<string, string>();
                foreach (var item in hk.dicFinal)
                    for (int i = 0; i < item.Value.Count; i++)
                        caseBase.Add(hk.structNames[item.Value[i]], item.Key);
                //     caseBase.Add(hk.structNames[item.Value[0]], item.Key);




                dist = new int[set.numberOfCores][];

                CreateBase(caseBase);
                //testData = caseBase;

                for (int i = 0; i < dist.Length; i++)
                    dist[i] = new int[caseBase.Count];

              

            }
        }
        public virtual string CalcDist(int threadNum, int[] keys, int[] index, int num)
        {
            return "";
        }
        public virtual string CalcDist(int threadNum, string [] keys, int[] index, int num)
        {
            int[] locDist = dist[threadNum];
            
            for (int i = 0; i < keys.Length; i++)
            {
                int[] baseHash = null;
                if (hashDataBaseTab[i].TryGetValue(keys[i], out baseHash))
                //if (auxBase[i].TryGetValue(keys[i], out baseHash))
                {
                    for(int j=0;j<baseHash.Length;j++)
                    //foreach (var inx in baseHash)
                        locDist[baseHash[j]]++;
                }

            }
            //Array.Sort(locDist, index);
            //Array.Reverse(locDist);
            //Array.Reverse(index);
            string w = "";
            for (int i = 0; i < resSize; i++)
                w += dataBaseKeys[index[i]] + ":";
            w += dataBaseKeys[index[resSize]];

            return w;

        }
        public virtual void CalcHashDist(object o)
        {
            
            ThreadParam pp = (ThreadParam)o;
            int threadNum = pp.num;
            int start = pp.start;
            int stop = pp.stop;

          
            int[] index = new int[dist[threadNum].Length];
            int[] locDist = dist[threadNum];
            for (int n = start; n < stop; n++)
            {
                string testItem = testList[n];

                if (!testData.ContainsKey(testItem))
                    continue;

                for (int i = 0; i < locDist.Length; i++)
                {
                    locDist[i] = 0;
                    index[i] = i;
                }
                int num = 0;

                string w;
                if (opt.mSim)
                {
                    int[] keys = GetKeyHashesTabInt(testData[testItem]);
                    w = CalcDist(threadNum, keys, index, num);
                }
                else
                {
                    string []keys = GetKeyHashesTab(testData[testItem]);
                    w = CalcDist(threadNum, keys, index, num);
                }


                lock (results)
                {
                    results.Add(testItem, w);
                }


            }
            resetEvents[threadNum].Set();

        }
       /* protected void FastSortNew(int[] dist, int[] index)
        {
            //int[] aux = new int[dist.Length];
           
            int[] pos = new int[101];
           
            int i = 0;
            for (i = 0; i < dist.Length; i++)
            {
                int v = dist[i];
                auxMem[v][pos[v]] = i;
                pos[v]++;
            }
            int counter = 0;
            for(i=0;i<auxMem.Length;i++)
                if(pos[i]>0)
                {
                    int[] tmp = auxMem[i];
                    for(int n=0;n<pos[i];n++)
                        index[counter++] = tmp[n];
                    if (counter > resSize)
                        break;
                }

        }*/
        protected void FastSort(int[] dist, int[] index)
        {
            //int[] aux = new int[dist.Length];
            int[] aux = new int[101];
            int i = 0;
            for (i = 0; i < dist.Length; i++)
                aux[dist[i]]++;

            i = 0;
            int num = 0, border = 0;
            while (border < resSize && num < aux.Length)
                border += aux[num++];

            int[] sorted = new int[border];
            int counter = 0;
            for (i = 0; i < dist.Length; i++)
            {
                if (dist[i] < num)
                {
                    sorted[counter] = dist[i];
                    index[counter++] = i;
                }               
            }

            Array.Sort(sorted, index);


        }

        void Retrival(List<string> testList)
        {
            this.testList = testList;

            resetEvents = new ManualResetEvent[dist.Length];
            for (int n = 0; n < resetEvents.Length; n++)
            {
                ThreadParam pp = new ThreadParam();
                pp.num = n;
                pp.start = (int)(n * testList.Count / resetEvents.Length);
                pp.stop = (int)((n + 1) * testList.Count / resetEvents.Length);
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(CalcHashDist), (object)pp);
            }

            for (int n = 0; n < resetEvents.Length; n++)
                resetEvents[n].WaitOne();
        }
        public override Dictionary<string, string> HNNTest(List<string> testList)
        {           
            Retrival(testList);
            return results;

        }

    }
}
