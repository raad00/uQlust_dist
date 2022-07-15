using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace uQlustCore
{
    [Serializable]
    class ColMinStateHash:RootHash
    {
        protected int[] minStateDataBase = null;
        //[NonSerialized]

        protected KeyValuePair<int,int[]>[,] oTab = null;
        Dictionary<int, KeyValuePair<int, int[]>>[] colBaseInt;


        public ColMinStateHash(int binSizeG, HashCluster hk, ClusterOutput outp, HNNCInput opt) : base(binSizeG, hk, outp, opt)
        {

        }
        public ColMinStateHash(int binSizeG, HNNCInput opt) : base(binSizeG, opt)
        {
            set = new Settings();
            set.Load();
        }
        public override void Preprocessing(bool load)
        {
            base.Preprocessing(load);
            if (load)
            {
                ColMinStateHash ob;
                try
                {
                    ob = (ColMinStateHash)GeneralFunctionality.LoadBinary(opt.binaryfile);
                    minStateDataBase = ob.minStateDataBase;
                    colBaseInt = ob.colBaseInt;
                    oTab = ob.oTab;
                    dataBaseKeys = ob.dataBaseKeys;
                    dist = new int[set.numberOfCores][];
                    for (int i = 0; i < dist.Length; i++)
                        dist[i] = new int[dataBaseKeys.Length];
                   
                }
                catch (Exception ex)
                {
                    throw new Exception("Error while loading databasez: "+opt.binaryfile);
                }
            }
        }
        int ConvertToInt(string item)
        {
            int v = 0;

            for (int k = 0; k < item.Length; k++)
            {
                int w = item[k] - '1';
                v += w*(int)Math.Pow(2, k);
            }
            return v;
        }
        int  GetKeyInt(string item)
        {                                              
            return ConvertToInt(item);
        }
        public override void CreateBase(Dictionary<string, string> dataBase)
        {
            Dictionary<int, KeyValuePair<int, List<int>>>[] colBase = null;
            int remBinSize = binSizeG;
            binSizeG = 1;
            base.CreateBase(dataBase);
            binSizeG = remBinSize;

            minStateDataBase = new int[dataBaseKeys.Length];
            for (int j = 0; j < minStateDataBase.Length; j++)
                minStateDataBase[j] = 0;

            for (int j = 0; j < hashDataBase.Length; j++)
            {
                string remKey = "";
                int min = int.MaxValue;
                Dictionary<string, List<int>> locHash = hashDataBase[j];
                if (hashDataBase[j].Keys.Count > 1)
                {
                    foreach (var item in hashDataBase[j].Keys)
                    {
                        if (locHash[item].Count < min)
                        {
                            min = locHash[item].Count;
                            remKey = item;
                        }
                    }

                    List<int> aux = hashDataBase[j][remKey];
                    locHash.Clear();
                    locHash.Add(remKey, aux);
                }
                else
                    locHash.Clear();
            }





            colBase = new Dictionary<int, KeyValuePair<int, List<int>>>[hashDataBase.Length / binSizeG + 1];


            for (int i = 0; i < colBase.Length; i++)
                colBase[i] = new Dictionary<int, KeyValuePair<int, List<int>>>();

            //dataBaseKeys = new List<string>(dataBase.Keys);

            for (int j = 0; j < dataBaseKeys.Length; j++)
            {
                List<string> aux = GetKeyHashes(dataBase[dataBaseKeys[j]]);
                
                for (int i = 0; i < aux.Count; i++)
                {
                    int auxInt = GetKeyInt(aux[i]);
                    int num = 0;
                    for (int n = 0; n < aux[i].Length; n++)
                        if (hashDataBase[i * aux[i].Length + n].ContainsKey(aux[i][n].ToString()))
                            num++;

                    if (num == 0)
                        continue;


                    minStateDataBase[j] += num;

                    if (colBase[i].ContainsKey(auxInt))
                        colBase[i][auxInt].Value.Add(j);
                    else
                    {
                        List<int> xaux = new List<int>();
                        xaux.Add(j);
                        colBase[i].Add(auxInt, new KeyValuePair<int, List<int>>(num, xaux));
                    }
                }

            }
            colBaseInt = new Dictionary<int, KeyValuePair<int, int[]>>[colBase.Length];
            oTab = new KeyValuePair<int,int[]>[colBase.Length,(int)Math.Pow(2,binSizeG)];
            for (int i = 0; i < colBase.Length; i++)
            {
                for(int s=0;s<oTab.GetLength(1);s++)
                {
                    oTab[i, s] = default(KeyValuePair<int,int[]>);
                }
                colBaseInt[i] = new Dictionary<int, KeyValuePair<int, int[]>>();
                foreach (var item in colBase[i])
                {
                    KeyValuePair<int, int[]> aux = new KeyValuePair<int, int[]>(item.Value.Key, item.Value.Value.ToArray());
                    colBaseInt[i].Add(item.Key, aux);
                    oTab[i, item.Key] = colBaseInt[i][item.Key];
                }
            }
            GeneralFunctionality.SaveBinary(opt.binaryfile, this);
        }

        /*        public override void CreateBase(Dictionary<string, string> dataBase)
                {
                    Dictionary<string, KeyValuePair<int, List<int>>>[] colBase = null;
                    int remBinSize = binSizeG;
                    binSizeG = 1;
                    base.CreateBase(dataBase);
                    binSizeG = remBinSize;

                    minStateDataBase = new int[dataBaseKeys.Length];
                    for (int j = 0; j < minStateDataBase.Length; j++)
                        minStateDataBase[j] = 0;

                    for (int j = 0; j < hashDataBase.Length; j++)
                    {
                        string remKey = "";
                        int min = int.MaxValue;
                        Dictionary<string, List<int>> locHash = hashDataBase[j];
                        if (hashDataBase[j].Keys.Count > 1)
                        {
                            foreach (var item in hashDataBase[j].Keys)
                            {
                                if (locHash[item].Count < min)
                                {
                                    min = locHash[item].Count;
                                    remKey = item;
                                }
                            }                 

                            List<int> aux = hashDataBase[j][remKey];
                            locHash.Clear();
                            locHash.Add(remKey, aux);
                        }
                        else
                            locHash.Clear();
                    }





                    colBase = new Dictionary<string, KeyValuePair<int, List<int>>>[hashDataBase.Length / binSizeG + 1];


                    for (int i = 0; i < colBase.Length; i++)
                        colBase[i] = new Dictionary<string, KeyValuePair<int, List<int>>>();

                    //dataBaseKeys = new List<string>(dataBase.Keys);

                    for (int j = 0; j < dataBaseKeys.Length; j++)
                    {
                        List<string> aux = GetKeyHashes(dataBase[dataBaseKeys[j]]);

                        for (int i = 0; i < aux.Count; i++)
                        {
                            int num = 0;
                            for (int n = 0; n < aux[i].Length; n++)
                                if (hashDataBase[i * aux[i].Length + n].ContainsKey(aux[i][n].ToString()))
                                    num++;

                            if (num == 0)
                                continue;


                           minStateDataBase[j]+=num;

                            if (colBase[i].ContainsKey(aux[i]))
                                colBase[i][aux[i]].Value.Add(j);
                            else
                            {
                                List<int> xaux = new List<int>();
                                xaux.Add(j);
                                colBase[i].Add(aux[i], new KeyValuePair<int, List<int>>(num, xaux));
                            }
                        }

                    }
                    colBaseInt = new Dictionary<string, KeyValuePair<int, int[]>>[colBase.Length];
                    for(int i=0;i<colBase.Length;i++)
                    {
                        colBaseInt[i] = new Dictionary<string, KeyValuePair<int, int[]>>();
                        foreach(var item in colBase[i])
                        {
                            KeyValuePair<int, int[]> aux = new KeyValuePair<int, int[]>(item.Value.Key, item.Value.Value.ToArray());
                            colBaseInt[i].Add(item.Key, aux);
                        }
                    }
                    GeneralFunctionality.SaveBinary(opt.binaryfile, this);
                }
        */

        public override string CalcDist(int threadNum, string [] keys, int[] index, int num)
        {
            int[] locDist = dist[threadNum];
            int numCount = 0;
            for (int i = 0; i < keys.Length; i++)
            {
                int v = ConvertToInt(keys[i]);                
                KeyValuePair<int,int[]> x = oTab[i, v];
                 int val = x.Key;
                if (val!=0)
                {
                    int[] baseHash = x.Value;                   
                    for(int j=0;j<baseHash.Length;j++)
                        locDist[baseHash[j]] += val;
                    numCount += val;
                }
            }

            for (int i = 0; i < locDist.Length; i++)
            {
                int v = locDist[i];
                int x = minStateDataBase[i] + numCount - v;
                locDist[i] = 100-(int)((100.0 * v) / x);
            }

            //Array.Sort(locDist, index);
            FastSort(locDist, index);

            string w = "";
            for (int i = 0; i < resSize; i++)
                w += dataBaseKeys[index[i]] + "-"+locDist[index[i]]+":";
            w += dataBaseKeys[index[resSize]]+"-"+locDist[index[resSize]];

            return w;
        }
    }
}
