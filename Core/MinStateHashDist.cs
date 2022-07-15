using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using System.Text;
using System.Threading.Tasks;

namespace uQlustCore
{
    [Serializable]
    class MinStateHash:RootHash
    {
        
        protected int[] minStateDataBase = null;   
        [NonSerialized]
        protected Dictionary<string, KeyValuePair<int, List<int>>>[] locBase = null;
        [NonSerialized]
        protected Dictionary<string, KeyValuePair<int, int[]>>[] locBaseTab = null;       
        protected int [,][] oTab = null;
        public MinStateHash(int binSizeG, HashCluster hk, ClusterOutput outp, HNNCInput opt):base(binSizeG,hk, outp, opt)
        {
            
        }
        public MinStateHash(int binSizeG,HNNCInput opt):base(binSizeG,null,null,opt)
        {
            
        }
       
        public override void Preprocessing(bool load)
        {
            base.Preprocessing(load);
            if(load)
            {
                MinStateHash ob=(MinStateHash)GeneralFunctionality.LoadBinary(opt.binaryfile);
                minStateDataBase = ob.minStateDataBase;
                oTab = ob.oTab;
                dataBaseKeys = ob.dataBaseKeys;
                dist = new int[set.numberOfCores][];
                for (int i = 0; i < dist.Length; i++)
                    dist[i] = new int[dataBaseKeys.Length];
            }
        }

        public override void CreateBase(Dictionary<string, string> dataBase)
        {
            int remBin = binSizeG;
            binSizeG = 1;

            base.CreateBase(dataBase);
            binSizeG = remBin;
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

                    foreach (var item in locHash[remKey])
                        minStateDataBase[item]++;

                    List<int> aux = hashDataBase[j][remKey];
                    locHash.Clear();
                    locHash.Add(remKey, aux);
                }
                else
                    locHash.Clear();
            }
           // binSizeG = 1;


            locBase = new Dictionary<string, KeyValuePair<int, List<int>>>[hashDataBase.Length / binSizeG + 1];

            for (int i = 0; i < locBase.Length; i++)
                locBase[i] = new Dictionary<string, KeyValuePair<int, List<int>>>();

            dataBaseKeys = dataBase.Keys.ToArray();
            //List<string> triples = GenerateTriplets();
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

                    if (locBase[i].ContainsKey(aux[i]))
                        locBase[i][aux[i]].Value.Add(j);
                    else
                    {
                        List<int> xaux = new List<int>();
                        xaux.Add(j);
                        locBase[i].Add(aux[i], new KeyValuePair<int, List<int>>(num, xaux));
                    }
                }

            }

            locBaseTab =new Dictionary<string, KeyValuePair<int, int[]>>[locBase.Length];
            oTab = new int [locBase.Length, 3][];
            for(int i=0;i<locBase.Length;i++)
            {
                oTab[i, 0] = null;
                oTab[i, 1] = null;
                oTab[i, 2] = null;

                locBaseTab[i] = new Dictionary<string, KeyValuePair<int, int[]>>();
                foreach(var item in locBase[i])
                {
                    KeyValuePair<int, List<int>> aux = item.Value;
                    KeyValuePair<int, int[]> newAux = new KeyValuePair<int, int[]>(aux.Key,aux.Value.ToArray());
                    locBaseTab[i].Add(item.Key, newAux);
                    oTab[i, Convert.ToInt32(item.Key[0]-'0')] = locBaseTab[i][item.Key].Value;
                }
            }
            GeneralFunctionality.SaveBinary(opt.binaryfile, this);

        }
        public override int[] GetKeyHashesTabInt(string item)
        {
            int[] aux = new int[(int)Math.Ceiling(((double)item.Length) / binSizeG)];
                       
            for (int i = 0, counter = 0; i < item.Length; i += binSizeG, counter++)
            {
                aux[counter] =Convert.ToInt32(item[i]-'0');
            }

            return aux;
        }
        public override string CalcDist(int threadNum, int[] keys, int[] index, int num)
        {
            int[] aux;
            int[] locDist = dist[threadNum];
            int counter = 0;
            int i = 0;
            for (i = 0; i < keys.Length; i++)
            {
                aux = oTab[i++,keys[i]];
                if (aux!=null)
                {
                    for (int j = 0; j < aux.Length; j++) 
                   //     foreach(var item in aux)
                            //locDist[item]++;
                          locDist[aux[j]] ++;
                    counter++;
                }

            }
            for (i = 0; i < locDist.Length; i++)
                locDist[i] = counter+minStateDataBase[index[i]] - 2 * locDist[i];

            
            Array.Sort(dist[threadNum], index);
            //FastSort(locDist, index);

            string w = "";
            for (i = 0; i < resSize; i++)
                w += dataBaseKeys[index[i]] + ":";
            w += dataBaseKeys[index[resSize]];

            return w;
        }

      

        public override string CalcDist(int threadNum, string [] keys, int[] index, int num)
        {
            KeyValuePair<int, int[]> aux;
            int[] locDist = dist[threadNum];
            Dictionary<string, KeyValuePair<int, int[]>> b = null;
            string k = "";
            int i = 0;
            //foreach(var item in keys)
            for (i = 0; i < keys.Length; i++)
            {
                k = keys[i] ;
                b = locBaseTab[i];
                if (b.ContainsKey(k))
                {
                    aux = b[k];
                    int x = aux.Key;
                    int[] common = aux.Value;
                    for (int j = 0; j < common.Length; j++)
        //            foreach(var item in common)
                        locDist[common[j]] += x;
                }

            }
            for (i = 0; i < locDist.Length; i++)
                locDist[i] = minStateDataBase[index[i]] - 2 * locDist[i];

            //Array.Sort(dist[threadNum], index);
            string w = "";
            for (i = 0; i < resSize; i++)
                w += dataBaseKeys[index[i]] + ":";
            w += dataBaseKeys[index[resSize]];

            return w;
        }


    }
}
