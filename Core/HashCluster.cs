using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using uQlustCore;
using System.Threading;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Data;
using uQlustCore.Interface;
using uQlustCore.Profiles;
//using Accord.Statistics.Testing;

namespace uQlustCore
{
    public class HashCluster:IProgressBar
    {
        public Dictionary<string, List<byte>> stateAlign;

        struct threadParam
        {
            public int threadNum;
            public double threshold;
            public Dictionary<int, int>[] referenceDist;
            public List<int> candidates;
            public int referenceIndex;
            public int aux1;
            public int aux2;
        };
        public List<string> stateAlignKeys { get { return new List<string>(stateAlign.Keys); } }


        public Dictionary<string, List<int>> dicC = null;
        public Dictionary<string, List<int>> dicFinal = null;
        public bool consensusProjection = true;
        jury1D jury=null;
        List<string> allItems = new List<string>();

        //List<RangeWeight> rangeWeights = null;
        public List<int> selectedColumnsHash = new List<int>();
        public Dictionary<string, string> structToKey = null;
        DCDFile dcd = null;
        public List<string> structNames = null;
        protected int maxV = 100;
        protected int currentV = 0;

        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }


        Dictionary<string, List<byte>> stateAlignReg;
        Settings dirSettings = new Settings();
        public Alignment al;
        protected string alignFile=null;
        string dirName;
       protected Dictionary<byte, int>[] consensusStates;
        Dictionary<byte, int>[] consensusStatesReg;
        protected Dictionary<byte, int>[] columns = null;
        Dictionary<byte, int>[] columnsReg = null;

        Dictionary<string, List<int>> threadingKeys;
        public List<int> validIndexes = new List<int>();

        //Neaded for thread Tasks
        ManualResetEvent[] resetEvents = null;
        List<int>[] threadStructures;
        protected List<string> allStructures;
        Dictionary<byte, int>[] col = null;
        char[][] keyTab;
        Dictionary<string, List<byte>> auxDic;
        int[] distTab;       
        
//############################


        protected HashCInput input = null;
        int threadNumbers;
        public int refPoints;
        public HashCluster(string dirName, string alignFile,HashCInput input)
        {
            this.dirName = dirName;
            this.refPoints = input.refPoints;
            this.alignFile = alignFile;
            this.input = input;
            dirSettings.Load();
            consensusProjection = input.useConsensusStates;
            threadNumbers = dirSettings.numberOfCores;
            maxV += refPoints * 20;
        }
        public HashCluster(DCDFile dcd, HashCInput input)
        {
            this.dcd = dcd;
            this.input = input;
            currentV = 0;
            maxV = 4;
            dirSettings.Load();
            threadNumbers = dirSettings.numberOfCores;
        }
        public HashCluster()
        {

        }
        public HashCluster(Alignment _al, HashCInput input)
        {
            al = _al;
            this.input = input;
            dirSettings.Load();
            consensusProjection = input.useConsensusStates;
            threadNumbers = dirSettings.numberOfCores;
            stateAlign = al.GetStateAlign();
            maxV += refPoints * 20;
        }
        public HashCluster(Dictionary<string, List<byte>> profiles, Alignment _al, HashCInput input)
        {
            this.input = input;
            this.refPoints = input.refPoints;
            stateAlign = profiles;
            al = _al;            
            dirSettings.Load();
            consensusProjection = input.useConsensusStates;
            threadNumbers = dirSettings.numberOfCores;
            maxV += refPoints * 20;
        }
        public override string ToString()
        {
            if(input!=null && !input.rpart)
                return "phiClust:Hash";

            return "phiClust:Rpart";
        }
        public virtual void InitHashCluster()
        {
            if (!input.floatingPointsProfiles)
            {
                if (dcd != null)
                    PrepareForDCD(dcd, input);
                else
              if (alignFile != null && alignFile.Length > 0)
                {

                    al = new Alignment();
                    al.Prepare(alignFile, input.profileName);
                    // al.Prepare(profiles, profName, input.profileName);

                    al.MyAlign(alignFile);
                    //al.CombineAll();
                    stateAlign = al.GetStateAlign();

                    List<string> keys = al.r.masterNode.Keys.ToList();

                    //  if (al.r.masterNode[keys[0]].rangeWeights != null)
                    //    rangeWeights = al.r.masterNode[keys[0]].rangeWeights;
                    //AddErrors(al.errors);
                    //           stateAlign = al.GetStateAlign();
                    if (input.regular)
                    {
                        if (input.profileName == input.profileNameReg)
                            stateAlignReg = stateAlign;
                        else
                        {
                            Alignment alReg = new Alignment();
                            if (alignFile != null)
                                alReg.Prepare(alignFile, input.profileNameReg);
                            else
                                alReg.Prepare(dirName, dirSettings, input.profileNameReg);
                            alReg.MyAlign(alignFile);
                            stateAlignReg = alReg.GetStateAlign();
                        }
                    }
                }

                if (stateAlign == null)
                    PrepareForPDB(dirName, alignFile, input);
            }

        }
        public double ProgressUpdate()
        {
            double res = 0;
            if (al != null)
                res = al.ProgressUpdate();
            return StartProgress + (EndProgress - StartProgress) * (res * 0.45 + 0.55 * (double)currentV / maxV);
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;            
        }
        void PrepareForDCD(DCDFile dcd, HashCInput input)
        {
            dirSettings.Load();
            this.dcd = dcd;
            al = new Alignment();
            al.Prepare(dcd, dirSettings, input.profileName);
            al.MyAlign(alignFile);
            stateAlign = al.GetStateAlign();
            if (input.regular)
            {
                if (input.profileName == input.profileNameReg)
                    stateAlignReg = stateAlign;
                else
                {

                    Alignment alReg = new Alignment();
                    alReg.Prepare(dcd, dirSettings, input.profileNameReg);
                    alReg.MyAlign(alignFile);
                    stateAlignReg = alReg.GetStateAlign();
                }
            }
            this.input = input;

        }
        void PrepareForPDB(string dirName, string alignFile,HashCInput input)
        {
            dirSettings.Load();
            this.dirName = dirName;
            al = new Alignment();
            if (alignFile != null)
                 al.Prepare(alignFile, input.profileName);
            else
            {
                al.Prepare(dirName, dirSettings, input.profileName);
            
            }
            al.MyAlign(alignFile);
            stateAlign = al.GetStateAlign();
            if (input.regular)
            {
                if (input.profileName == input.profileNameReg)
                    stateAlignReg = stateAlign;
                else
                {
                    Alignment alReg = new Alignment();

                    if(alignFile!=null)
                        alReg.Prepare(alignFile, input.profileNameReg);
                    else
                        alReg.Prepare(dirName, dirSettings, input.profileNameReg);
                    alReg.MyAlign(alignFile);
                    stateAlignReg = alReg.GetStateAlign();
                }
            }
            this.input = input;

        }
        public double[,] KeyTable(double[,] table)
        {
            double [,] resTable=new double [table.GetLength(0),table.GetLength(1)];
            int[] aux = new int[table.GetLength(0)];
            for (int i = 0; i < table.GetLength(1); i++)
            {
                for (int j = 0; j < aux.Length; j++)
                    aux[j] = 0;
                for (int j = 0; j < table.GetLength(0); j++)                
                    aux[(int)table[j,i]]++;
                int sum=0;
                int remIndex=0;
                for (int j = 0; j < aux.Length; j++)
                {
                    sum += aux[j];
                    if (sum > table.GetLength(0) / 2)
                    {
                        remIndex = j;
                        break;
                    }
                }
                for (int j = 0; j < table.GetLength(0); j++)
                {
                    if (table[j,i] <= remIndex)
                        resTable[j,i] = 0;
                    else
                        resTable[j,i] = 1;
                }

            }


            return resTable;
        }
        private Dictionary<string, List<int>> MetaColumns(List<string> lStruct, int k, double percent)
        {
            if (columns == null)
                return null;

            Dictionary<string, List<int>> metaColumns = new Dictionary<string, List<int>>();
            Dictionary<string, int>[] states;
            Dictionary<string, int> keysNum = new Dictionary<string, int>();
            int[] indexTable = new int[columns.Length/3];

            for (int n = 0; n < lStruct.Count; n++)
                keysNum.Add(lStruct[n], n);

            for (int i = 0; i < indexTable.Length; i++)            
                indexTable[i] = i;
            
            int counter = 0;
            states = new Dictionary<string, int>[columns.Length / 3];
            for (int i = 0; i < states.Length; i++)
                states[i] = new Dictionary<string, int>();

                foreach (var item in lStruct)
                {
                    
                    counter = 0;
                    for (int i = 1; i < columns.Length - 1; i += 3)
                    {
                        byte state = stateAlign[item][i - 1];

                        if (stateAlign[item][i - 1] == stateAlign[item][i] && stateAlign[item][i] == stateAlign[item][i + 1])
                            states[counter].Add(item,0);
                        else
                            states[counter].Add(item,1);
                        counter++;

                    }
                }

            Dictionary<string, List<int>> bestSolution = new Dictionary<string, List<int>>();
            Dictionary<int, int> solutionsList = new Dictionary<int, int>();
            int positionLeft = 0;
            int positionRight = counter;
            int positionMiddle;
            do
            {
                metaColumns.Clear();
                positionMiddle = positionLeft + (positionRight - positionLeft) / 2;
                StringBuilder keyB = new StringBuilder(positionMiddle);
                foreach(var item in lStruct)
                //foreach (var item in states.Keys)
                {
                    keyB.Clear();                    
                    for (int i = 0; i < positionMiddle; i++)
                        keyB.Append(states[indexTable[i]][item]);

                    string key = keyB.ToString();
                    if (!metaColumns.ContainsKey(key))
                        metaColumns.Add(key, new List<int>());

                    metaColumns[key].Add(keysNum[item]);

                }
                if (metaColumns.Keys.Count < k)
                {
                    positionLeft = positionMiddle;              
                }
                else
                {
                    List<KeyValuePair<string, List<int>>> toSort = metaColumns.ToList();

                    toSort.Sort((nextPair,firstPair) => { return firstPair.Value.Count.CompareTo(nextPair.Value.Count); });

                    int sum=0;
                    for (int i = 0; i < k; i++)
                        sum += toSort[i].Value.Count;

                    double res;

                    res = ((double)sum )/ stateAlign.Keys.Count;

                    if (res >= percent)
                        positionLeft = positionMiddle;
                    else
                        positionRight = positionMiddle;
                }
            }
            while (positionRight-positionLeft>1);

            
            return metaColumns;
        }
       //Return 1 increase number of columns , 0 decrease
        //Return 0 if selction should be stoped, otherwise it returns number of columns that should be used
        double CalcMatching(List<List<byte>> cluster)
        {
            double res = 0;
            Dictionary<byte, int>[] freq = null;
            freq = CalcFreq(cluster);

            for (int i = 0; i < freq.Length; i++)
            {
                if (freq[i].Keys.Count == 1)
                    res++;
            }
            res /= freq.Length;

            return res;
        }

        double CalcMatching(List<string> cluster)
        {
            double res = 0;
            Dictionary<char, int>[] freq = null;
            freq=CalcFreq(cluster);

            for (int i = 0; i < freq.Length; i++)
            {
                if (freq[i].Keys.Count == 1)
                    res++;
            }
            res /= freq.Length;

            return res;      
        }
        private int SelectionNMatchHeuristic(int k, double percent, Dictionary<string, List<int>> dic)
        {
            double res = 0;
            int counter = 0;
            foreach (var item in dic)
            {
                if (item.Value.Count > 1)
                {
                    List<List<byte>> clustProfiles = new List<List<byte>>(item.Value.Count);
                    foreach (var p in item.Value)
                        clustProfiles.Add(auxDic[allStructures[p]]);

                    res += CalcMatching(clustProfiles);
                    counter++;
                }

            }
            if(counter>0)
                res /= counter;

            if (res < percent)
                return 0;
            return 1;


        }

        private int SelectionHeuristic(int k, Dictionary<string, List<int>> dic)
        {
            if (dic.Keys.Count < k)
                return 1;
            else
                return 0;
        }
        private int SelectionHeuristic(int k, double percent, Dictionary<string, List<int>> dic)
        {
            if (dic.Keys.Count < k)
                return 0;
            else
            {
                List<KeyValuePair<string, List<int>>> toSort = dic.ToList();

                toSort.Sort((nextPair, firstPair) => { return firstPair.Value.Count.CompareTo(nextPair.Value.Count); });

                int sum = 0;
                for (int i = 0; i < k; i++)
                    sum += toSort[i].Value.Count;

                double res;

                res = ((double)sum) / stateAlign.Keys.Count;

                if (res >= percent)
                    return 0;
                else
                    return 1;
            }

        }
        private Dictionary<byte, int>[] CalcFreq(List<List<byte>> profiles)
        {
            Dictionary<byte, int>[] freq = null;

            freq = new Dictionary<byte, int>[profiles[0].Count];

            for (int i = 0; i < freq.Length; i++)
                freq[i] = new Dictionary<byte, int>();

            foreach (var item in profiles)
            {
                for (int i = 0; i < item.Count; i++)
                {
                    if (freq[i].ContainsKey(item[i]))
                        freq[i][item[i]]++;
                    else
                        freq[i].Add(item[i], 1);
                }
            }

            return freq;
        }

        private Dictionary<char,int>[] CalcFreq(List<string> profiles)
        {
            Dictionary<char, int>[] freq = null;

            freq = new Dictionary<char, int>[profiles[0].Length];

            for (int i = 0; i < freq.Length; i++)
                freq[i] = new Dictionary<char, int>();

            foreach (var item in profiles)
            {
                for (int i = 0; i < item.Length; i++)
                {
                    if (freq[i].ContainsKey(item[i]))
                        freq[i][item[i]]++;
                    else
                        freq[i].Add(item[i], 1);
                }
            }

            return freq;
        }
        private List<KeyValuePair<int, double>> HashKeyEntropy(Dictionary<string, List<int>> dic)
        {
            List<KeyValuePair<int,double>> entropy=null;
            //double[] freq = null;
            Dictionary<char, int>[] freq = null; 
            foreach(var item in dic.Keys)
            {
                entropy = new List<KeyValuePair<int,double>>(item.Length);
                break;
            }
            List<string> profiles = dic.Keys.ToList();
            freq = CalcFreq(profiles);
            for (int i = 0; i < freq.Length;i++)
            {
                double sum=0;
                foreach(var item in freq[i])
                {
                    double prob;
                    if(item.Key=='0')
                        prob= 1.5*((double)item.Value)/ dic.Keys.Count;
                    else
                        prob=0.5* ((double)item.Value)/ dic.Keys.Count;
                    sum += prob * Math.Log(prob);
                }
                entropy.Add(new KeyValuePair<int, double>(i,-sum));
                
            }

            entropy.Sort((firstPair, nextPair) =>
            {
                return nextPair.Value.CompareTo(firstPair.Value);
                //return firstPair.Value.CompareTo(nextPair.Value);
            });


            return entropy;
        }
        protected virtual double[] CalcEntropy(Dictionary<string,List<int>> dic)
        {
            List<string> keys = new List<string>(dic.Keys);


            double[] entropy = new double[keys[0].Length];
            for (int i = 0; i < keys[0].Length; i++)
            {
                int counter = 0;
                for (int j = 0; j < keys.Count; j++)
                    if (keys[j][i] == '0')
                        counter++;

                entropy[i] = 0;
                double prob = 0;
 //               entropy[i] = counter;
                prob=((double)counter) / keys.Count;
                if(prob>0)
                    entropy[i]= -prob * Math.Log10(prob);
                if((1-prob)>0)
                    entropy[i] += -(1-prob) * Math.Log10(1-prob);

            }
            return entropy;

        }

        protected virtual double [] CalcEntropy(Dictionary<byte,int> []locColumns)
        {
            if (locColumns == null)
                return null;

            double[] entropy = new double[locColumns.Length];
            for (int i = 0; i < locColumns.Length; i++)
            {
                int counter = 0;
                foreach (var item in locColumns[i])
                    counter += item.Value;

                entropy[i] = 0;
                double prob = 0;
                foreach (var item in locColumns[i])
                {
                    prob = ((double)item.Value) / counter;
                    if (prob > 0)
                        entropy[i] += -prob * Math.Log10(prob);
                }

            }
            return entropy;

        }
        private void AddChildrens(HClusterNode parent, Dictionary<string, string> dic,int []indexes,int k)
        {
            Dictionary<string, List<string>> hashClusters = new Dictionary<string, List<string>>();
            foreach(var item in parent.setStruct)
            {
                string key = "";
                for (int n = 0; n < k; n++)                
                    key += dic[item][indexes[n]];
                if (!hashClusters.ContainsKey(key))
                    hashClusters.Add(key, new List<string>());
                
                hashClusters[key].Add(item);

            }
            if (hashClusters.Keys.Count > 1)
            {
                parent.joined = new List<HClusterNode>();
                foreach (var item in hashClusters)
                {
                    HClusterNode aux = new HClusterNode();
                    aux.parent = parent;
                    aux.setStruct = item.Value;
                    aux.levelDist = indexes.Length - k;
                    aux.realDist = aux.levelDist;
                    parent.joined.Add(aux);
                }
            }
            else
            {
                parent.joined = null;
                if (k == indexes.Length-1)
                    parent.levelDist = 0;
            }
        }
        public ClusterOutput DendrogHashEntropy(Dictionary<string, List<int>> dic, List<string> structures)
        {
            int[] indexes = new int[columns.Length];
            HClusterNode root = new HClusterNode();
            Queue<KeyValuePair<HClusterNode, int>> queue = new Queue<KeyValuePair<HClusterNode, int>>();
            for (int j = 0; j < indexes.Length; j++)
                indexes[j] = j;

            double[] entropy = CalcEntropy(columns);
            Array.Sort(entropy, indexes);
            Dictionary<string, string> dicStructKey = new Dictionary<string,string>();

            foreach(var item in dic.Keys)
                foreach(var str in dic[item])
                    dicStructKey.Add(structures[str],item);

            root.setStruct = new List<string>(structures);
            root.parent = null;
            root.levelDist = indexes.Length;
            root.realDist = root.levelDist;
            queue.Enqueue(new KeyValuePair<HClusterNode,int>(root,1));
            while(queue.Count!=0)
            {
                KeyValuePair<HClusterNode, int> aux = queue.Dequeue();
                AddChildrens(aux.Key,dicStructKey,indexes,aux.Value);
                if (aux.Value + 1 < indexes.Length)
                {
                    if (aux.Key.joined != null)
                    {
                        foreach (var item in aux.Key.joined)
                            queue.Enqueue(new KeyValuePair<HClusterNode, int>(item, aux.Value + 1));
                    }
                    else
                        queue.Enqueue(new KeyValuePair<HClusterNode, int>(aux.Key, aux.Value + 1));
                }
                
            }
            ClusterOutput cOut = new ClusterOutput();
            cOut.hNode = root;

            return cOut;
        }
        private string CuttKey(string keyToCutt, bool []columnAvoid)
        {
            StringBuilder keyB = new StringBuilder();
            string key = "";
            for (int i = 0; i < keyToCutt.Length; i++)
            {
                if (columnAvoid[i])
                    continue;

                keyB.Append(keyToCutt[i]);
            }
            key = keyB.ToString();

            return key;
        }

        void ThreadingAddToClusters(object o)
        {
            List<string> outList = new List<string>();
            object[] array = o as object[];
            int threadNum = (int)array[0];
            int start = (int)array[1];
            int stop = (int)array[2];
            int maxPosition = (int)array[3];
            List<string> keyQueryLocal = (List<string>)array[4];
            List<string> caseKeyLocal = (List<string>)array[5];            
            Dictionary<string, List<string>> finalList = (Dictionary<string, List<string>>)array[6];
            List<string> caseKey = (List<string>)array[7];
            List<string> keyQuery = (List<string>)array[8];
            for (int i = start; i < stop; i++)
            {
                int position = 0;
                int positionLeft = 0;
                int positionRight = maxPosition;

                do
                {
                    position = (positionRight + positionLeft) / 2;
                    outList.Clear();
                    string key = "";
                    string refKey = keyQueryLocal[i].Substring(0,position);
                   // foreach (var it in caseKeyLocal)
                        for (int j = 0; j < caseKeyLocal.Count;j++ )
                        {
                            key = caseKeyLocal[j].Substring(0, position);
                            if (key == refKey)
                                outList.Add(caseKey[j]);
                        }

                    if (outList.Count == 1)
                        positionLeft = positionRight;
                    else
                        if (outList.Count > 1)
                            positionLeft = position;
                        else
                            positionRight = position;

                }
                while (positionRight - positionLeft > 1);
                lock (finalList)
                {
                    if (!finalList.ContainsKey(keyQuery[i]))
                        finalList.Add(keyQuery[i], outList);
                }
                outList = new List<string>();

            }
            resetEvents[threadNum].Set();
        }

        public Dictionary<string, List<string>> AddToClusters(List<string> caseKeys, List<string> keyQuery)
        {
            List<string> outList = new List<string>();
            Dictionary<string, List<string>> finalList = new Dictionary<string, List<string>>();
            if (columns == null)
                return null;
            DebugClass.WriteMessage("Select");
            bool[] columnAvoid = new bool[columns.Length];

            //int[] indexes = GetIndexOrder(columns);

            List<string> caseKeyLocal = new List<string>(caseKeys.Count);
            StringBuilder keyB = new StringBuilder();
            foreach (var item in caseKeys)
            {
                keyB.Clear();
                for (int i = 0; i < validIndexes.Count; i++)
                    keyB.Append(item[i]);// item[validIndexes[i]]);

                caseKeyLocal.Add(keyB.ToString());

            }

            List<string> keyQueryLocal = new List<string>(keyQuery.Count);
            foreach (var item in keyQuery)
            {
                keyB.Clear();
                for (int i = 0; i < validIndexes.Count; i++)
                    keyB.Append(item[i]);//[validIndexes[i]]);

                keyQueryLocal.Add(keyB.ToString());

            }


            Dictionary<string, List<int>> hashClusters = new Dictionary<string, List<int>>();

            threadNumbers = 3;
            resetEvents = new ManualResetEvent[threadNumbers];
            for (int n = 0; n < threadNumbers; n++)
            {
                int p = n;
                int start = n * keyQueryLocal.Count / threadNumbers;
                int stop = (n + 1) * keyQueryLocal.Count / threadNumbers;
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(ThreadingAddToClusters), new object[] { p, start, stop, validIndexes.Count, keyQueryLocal, caseKeyLocal, finalList, caseKeys, keyQuery });


            }
            for (int i = 0; i < threadNumbers; i++)
                resetEvents[i].WaitOne();


            /*            foreach (var item in keyQueryLocal)
                        {
                            int position = 0;
                            int positionLeft = 0;
                            int positionRight = indexes.Length;

                            do
                            {
                                position = (positionRight + positionLeft) / 2;
                                outList.Clear();
                                string key = "";
                                string refKey = item.Substring(0,position);
                                foreach (var it in caseKeyLocal)
                                {
                                    key = it.Substring(0, position);
                                    if (key == refKey)
                                        outList.Add(it);
                                }

                                if (outList.Count == 1)
                                    positionLeft = positionRight;
                                else
                                    if (outList.Count > 1)
                                        positionLeft = position;
                                    else
                                        positionRight = position;

                            }
                            while (positionRight - positionLeft > 1);

                            finalList.Add(item, outList);
                            outList = new List<string>();

                        }*/
            return finalList;
        }
        private Dictionary<string, List<int>> SelectColumnsByEntropyTest(Dictionary<string, List<int>> dic, int k, double prec)
        {
            if (columns == null)
                return null;
            bool[] columnAvoid = new bool[columns.Length];
            for (int i = 0; i < columns.Length; i++)
                columnAvoid[i] = false;

            int[] indexes = new int[columnAvoid.Length];

            for (int j = 0; j < indexes.Length; j++)
                indexes[j] = j;

            double[] entropy = CalcEntropy(columns);

            Array.Sort(entropy, indexes);
            // Array.Reverse(entropy);
            //Array.Reverse(indexes);

            Dictionary<string, List<int>> hashClusters = new Dictionary<string, List<int>>();

            int n = 0;
            while (entropy[n] < 0.05)
                columnAvoid[indexes[n++]] = true;

            int position = 0;
            int positionLeft = 0;
            int positionRight = indexes.Length;
            Dictionary<int, double> xx = new Dictionary<int, double>();
            do
            {
                position = (positionRight + positionLeft) / 2;

                for (int j = 0; j < position; j++)
                    columnAvoid[indexes[j]] = true;

                hashClusters.Clear();
                xx.Clear();
                for (int i = 0; i < columnAvoid.Length; i++)
                {
                    if (columnAvoid[i])
                        continue;
                    xx.Add(i, entropy[i]);
                }
                foreach (var item in dic.Keys)
                {
                    StringBuilder keyB = new StringBuilder();
                    string key = "";
                    for (int i = 0; i < item.Length; i++)
                    {
                        if (columnAvoid[i])
                            continue;

                        keyB.Append(item[i]);
                    }
                    key = keyB.ToString();
                    if (!hashClusters.ContainsKey(key))
                        hashClusters.Add(key, new List<int>());

                    hashClusters[key].AddRange(dic[item]);
                }
                for (int j = 0; j < position; j++)
                    columnAvoid[indexes[j]] = false;

                /*if (SelectionNMatchHeuristic(k, prec, hashClusters) == 1)
                    positionLeft = position;
                else
                    positionRight = position;*/

                if (SelectionHeuristic(k, prec, hashClusters) == 1)
                    positionLeft = position;
                else
                    positionRight = position;



            }
            while (positionRight - positionLeft > 1);

            for (int i = 0; i < position; i++)
                validIndexes.Add(indexes[i]);


            return hashClusters;

        }
        Dictionary<string,List<int>> CheckSeparation(Dictionary<string, List<int>> dic)
        {
            List<string> k = new List<string>(dic.Keys);
            Dictionary<string, List<string>> locCluster = new Dictionary<string,List<string>>();            
            List<int> featuresOrder = new List<int>();
            List<double> distInside = new List<double>();
            List<double> distOutside = new List<double>();


            for (int q = 0; q < k[0].Length; q++)
            {
                distInside.Clear();
                distOutside.Clear();
                locCluster.Clear();
                Random r = new Random();
                double remColSm = double.MaxValue;
                int remIndex = -1;
                for (int j = 0; j < k[0].Length; j++)
                {
                    if (featuresOrder.Contains(j))
                        continue;
                    locCluster.Clear();
                    foreach (var cluster in dic)
                    {
                        StringBuilder newKey = new StringBuilder();                        
                        if (featuresOrder.Count > 0)
                        {
                            foreach (var feat in featuresOrder)
                                newKey.Append(cluster.Key[feat]);                                
                        }
                        newKey.Append(cluster.Key[j]);
                        string ckey = newKey.ToString();
                        if (locCluster.ContainsKey(ckey))
                            locCluster[ckey].Add(cluster.Key);
                        else
                        {
                            List<string> aux = new List<string>();
                            aux.Add(cluster.Key);
                            locCluster.Add(ckey, aux);
                        }
                    }
                    
                    List<string> locKey = new List<string>(locCluster.Keys);
                    distInside.Clear();
                    int all = 0;
                    foreach (var cluster in locCluster)
                    {
                        all += cluster.Value.Count;
                    }
                    foreach (var cluster in locCluster)
                    {

                        for (int i = 1; i < cluster.Value.Count; i++)
                        {
                            bool test = true;
                            double f = 60000.0 / all;
                            if (f < 1)
                                if (r.NextDouble() > f)
                                    test = false;

                            if (test)
                            {
                                double sum = 0;
                                for (int n = 0; n < cluster.Value[0].Length; n++)
                                {
                                    int x = cluster.Value[0][n] - cluster.Value[i][n];
                                    sum += x * x;
                                }
                                distInside.Add(sum);
                            }
                        }
                    }
                    distOutside.Clear();
                   /* double clSize = locKey.Count * (locKey.Count + 1) / 2.0;
                    
                    for (int n = 0; n < locKey.Count; n++)
                        for (int i = n + 1; i < locKey.Count; i++)
                        {
                            bool test = true;
                            if (clSize >  distInside.Count)
                            {
                                
                                double f = (distInside.Count) / clSize;
                                if (r.NextDouble() > f)
                                    test = false;
                            }
                            if (test)
                            {
                                double sum = 0;
                                for (int s = 0; s < locCluster[locKey[n]][0].Length; s++)
                                {
                                    int x = locCluster[locKey[n]][0][s] - locCluster[locKey[i]][0][s];
                                    sum += x * x;
                                }
                                distOutside.Add(sum);
                            }
                        }*/

                    double[] samplesIn = distInside.ToArray();
                    //double[] samplesOut = distOutside.ToArray();

/*                    if (samplesOut.Length > 10 && samplesIn.Length > 10)
                    {
                        var test = new TwoSampleKolmogorovSmirnovTest(samplesIn, samplesOut);

                        if (test.Statistic > remColSm)
                        {
                            remColSm = test.Statistic;
                            remIndex = j;
                        }
                    }
                    else*/
                    {
                        double avrIn = 0;
                        foreach (var xx in samplesIn)
                            avrIn += xx;
                        avrIn /= samplesIn.Length;
                        double avrOut = 0;
                       // foreach (var xx in samplesOut)
                            //avrOut += xx;
                        //avrOut /= samplesOut.Length;
                        if(Math.Abs(avrOut-avrIn)<remColSm)
                        {
                            remColSm = Math.Abs(avrOut - avrIn);
                            remIndex = j;
                        }
                    }
                }
                featuresOrder.Add(remIndex);
                if (locCluster.Count > input.relClusters)
                    break;
            }
            Dictionary<string, List<int>> outDic = new Dictionary<string, List<int>>();

            foreach (var cluster in dic)
            {
                StringBuilder newKey = new StringBuilder();
                if (featuresOrder.Count > 0)
                {
                    foreach (var feat in featuresOrder)
                        newKey.Append(cluster.Key[feat]);
                }
                string ckey = newKey.ToString();
                if (outDic.ContainsKey(ckey))                    
                    outDic[ckey].AddRange(cluster.Value);
                else
                {
                    List<int> aux = new List<int>(cluster.Value);                    
                    outDic.Add(ckey,aux);
                }
            }

            return outDic;
         //   var test = new TwoSampleKolmogorovSmirnovTest(sample1, sample2);
        }
        double CheckAcc(Dictionary<string, List<int>> dic, List<string> structures,Dictionary<string,string>labels)
        {
            double good = 0;
            double all = 0;
            Dictionary<string, int> resDic = new Dictionary<string, int>();
            foreach(var item in dic)
            {
                resDic.Clear();
                foreach(var clItem in item.Value)
                {
                    try
                    {
                        if (labels.ContainsKey(structures[clItem]))
                        {
                            if (resDic.ContainsKey(labels[structures[clItem]]))
                                resDic[labels[structures[clItem]]]++;
                            else
                                resDic.Add(labels[structures[clItem]], 1);
                        }
                    }
                    catch(Exception ex)
                    {
                        Console.WriteLine("ups");
                    }

                }
                if (resDic.Count > 0)
                {
                    List<string> keys = new List<string>(resDic.Keys);
                    try
                    {
                        if (resDic.Count > 1)
                            keys.Sort(delegate (string x, string y)
                            {
                                return resDic[x].CompareTo(resDic[y]);
                            });
                        keys.Reverse();
                        good += resDic[keys[0]];
                        all += item.Value.Count;
                    }
                    catch (Exception ex)
                    {
                        Console.Write("ips2");
                    }
                }
            }
            return good/all;
        }
        void AddColumnsOneByOne(Dictionary<string, List<int>> dic, List<string> structures)
        {
            List<string> keys = new List<string>(dic.Keys);
            //StreamReader labFile = new StreamReader("C:\\Projects\\labels.txt");
            StreamReader labFile = new StreamReader("C:\\Projects\\labels_rna");
            //StreamReader labFile = new StreamReader("C:\\Projects\\labels_rfam");
            Dictionary<string, string> labDic = new Dictionary<string, string>();
            string line = labFile.ReadLine();
            while(line!=null)
            {
                line = line.TrimEnd();
                string[] aux = line.Split(' ');
                if (!labDic.ContainsKey(aux[0]))
                    labDic.Add(aux[0], aux[1]);
                else
                    DebugClass.WriteMessage("key already exists " + aux[0] + " value: "+labDic[aux[0]]+" new value: "+aux[1]);
                line = labFile.ReadLine();
            }
            labFile.Close();

            Dictionary<string, List<int>> hashCluster = new Dictionary<string, List<int>>();
            int prev = 0;
            double res = 0;
            double max = 0;
            int locIndex = -1;
            List<int> indexList = new List<int>();            
            for (int i = 1; i < keys[0].Length; i++)
            {
                hashCluster.Clear();
                string key = "";
                StringBuilder tmpKey = new StringBuilder();
                max = 0;
                for (int j = 0; j < keys[0].Length; j++)
                {
                    if (indexList.Count > 0 && indexList.Contains(j))
                        continue;
                    foreach (var item in dic.Keys)
                    {
                        try
                        {                            
                            if (indexList.Count > 0)
                            {
                                foreach (var idx in indexList)
                                        tmpKey.Append(item[idx]);
                                        //key += item.Substring(idx, 1);
                            }
                            tmpKey.Append(item[j]);
                            key = tmpKey.ToString();
                            tmpKey.Clear();
                            if (!hashCluster.ContainsKey(key))
                                hashCluster.Add(key, new List<int>());

                            hashCluster[key].AddRange(dic[item]);
                        }
                        catch (Exception sx)
                        {
                            Console.WriteLine("hhhf");
                        }
                    }
                    res = CheckAcc(hashCluster, structures, labDic);
                    if (res>max)
                    {
                        max = res;
                        locIndex = j;
                    }
                    hashCluster.Clear();
                }
                indexList.Add(locIndex);
                //indexList.Remove(locIndex);
            }               
        }
    Dictionary<string, List<int>> DropColumnsOneByOne(Dictionary<string, List<int>> dic,List<string> structures)
        {
            List<string> keys = new List<string>(dic.Keys);
            Dictionary<string, List<int>> hashCluster=null;
            int prev = 0;
            for (int i=1;i<keys[0].Length;i++)
            {
                hashCluster = ClusterByDropColumns(dic, i);
                if(hashCluster.Count!=prev)
                {
                    ClusterOutput output=DictionaryToClusterOutput(hashCluster, structures);
                    string fileName = "Clust_" + i;
                    output.SaveTxt(fileName);
                    prev = hashCluster.Count;
                }
            }

            return hashCluster;
        }
        Dictionary<string, List<int>> ClusterByDropColumns(Dictionary<string, List<int>> dic,int refPosition)
        {
            Dictionary<string, List<int>> hashClusters = new Dictionary<string, List<int>>();                        

            hashClusters.Clear();
            string key = "";
            foreach (var item in dic.Keys)
            {
                key = item.Substring(0, refPosition);
                if (!hashClusters.ContainsKey(key))
                    hashClusters.Add(key, new List<int>());

                hashClusters[key].AddRange(dic[item]);
            }
            return hashClusters;
        }
        Dictionary<string, List<int>> SelectColumnsByEntropy(Dictionary<string, List<int>> dic, int k, double prec)
        {
            if (columns == null)
                return null;

            /*            bool [] columnAvoid=new bool [columns.Length];
                        for (int i = 0; i < columns.Length; i++)
                            columnAvoid[i] = false;


                       // if (dic.Keys.Count < k)
                         //   k = dic.Keys.Count -1;

                        int[] indexes = new int[columnAvoid.Length];

                        for (int j = 0; j < indexes.Length; j++)
                            indexes[j] = j;

                        double[] entropy = CalcEntropy(columns);*/

            double[] entropy = CalcEntropy(dic);


            int[] indexes = new int[entropy.Length];

            for (int j = 0; j < indexes.Length; j++)
                indexes[j] = j;                        

           // Array.Sort(entropy, indexes);
           //Array.Reverse(entropy);
           //Array.Reverse(indexes);
            Dictionary<string, List<int>> hashClusters = new Dictionary<string, List<int>>();

            int position=0;

            int positionLeft = 0;
            int positionRight = indexes.Length;
            do
            {
                position =(positionRight + positionLeft) / 2;

                hashClusters.Clear();

/*                selectedColumnsHash.Clear();
                for (int i = 0; i < position; i++)
                    selectedColumnsHash.Add(indexes[i]);*/
                StringBuilder keyB = new StringBuilder();
                foreach (var item in dic.Keys)
                {
                    keyB.Clear();
                    string key = "";
                    for (int i = 0; i < position; i++)
                        keyB.Append(item[indexes[i]]);
                       
                    key = keyB.ToString();
                    if (!hashClusters.ContainsKey(key))
                        hashClusters.Add(key, new List<int>());

                    hashClusters[key].AddRange(dic[item]);
                }

                if (dic.Count < k)
                {
                    if (SelectionHeuristic(dic.Count, hashClusters) != 1)
                        positionRight = position;
                    else
                        positionLeft = position;
                }
                else
                    if (hashClusters.Count == dic.Count)
                        positionRight = position;
                    else
                        if (SelectionHeuristic(k, prec, hashClusters) != 1)
                            positionLeft = position;
                        else
                            positionRight = position;


                
                
            }
            while (positionRight-positionLeft>1);
            DebugClass.WriteMessage("Mean selected:");
            double sum = 0;
            for (int i = 0; i < position; i++)
                sum += entropy[i];
            sum /= position;
            DebugClass.WriteMessage("mean: "+sum);
            sum = 0;

            DebugClass.WriteMessage("Mediana: "+entropy[entropy.Length/2]);

            DebugClass.WriteMessage("");


            return hashClusters;

        }
        private Dictionary<string, List<int>> PrepareNewKeys(Dictionary<string, List<int>> dic)
        {
            if (columns == null)
                return null;
            bool[] columnAvoid = new bool[columns.Length];
            for (int i = 0; i < columns.Length; i++)
                columnAvoid[i] = false;

            int[] indexes = new int[columnAvoid.Length];

            for (int j = 0; j < indexes.Length; j++)
                indexes[j] = j;
            Random ra=new Random();
            int referencePoints = 10;
            Dictionary<string, List<int>>[] dic2 = new Dictionary<string, List<int>>[referencePoints];
            char [,][]profileKeys = new char [referencePoints + 1,allStructures.Count][];
            



            foreach (var item in dic)
                foreach (var index in item.Value)
                {
                    profileKeys[0, index] = new char[columns.Length];
                    for (int i = 0; i < columns.Length; i++)
                        profileKeys[0, index][i] = item.Key[i];
                }
        
            for (int n = 0; n < referencePoints; n++)
            {                
                List<byte> auxXX = auxDic[allStructures[ra.Next(allStructures.Count())]];
                
                for (int i = 0; i < auxXX.Count; i++)
                {
                    consensusStates[i].Clear();
                    consensusStates[i].Add(auxXX[i], 0);
                }
                dic2[n] = PrepareKeys(allStructures, true, false);

                foreach (var item in dic2[n])
                    foreach (var index in item.Value)
                    {
                        profileKeys[n + 1, index] = new char[columns.Length];
                        for (int i = 0; i < columns.Length; i++)
                            profileKeys[n + 1, index][i] = item.Key[i];
                    }
            }
            Dictionary<string, List<int>> newDic = new Dictionary<string, List<int>>();
            for (int t = 0; t < columns.Length;t++)
            {
                int sum = 0;
                foreach(var item in dic)
                {
                    if (item.Key[t] == '0')
                        sum++;
                }
                int remIndex = 0 ;
                int maxValue=sum;
                for (int n = 0; n < referencePoints; n++)
                {
                    int v = 0;
                    foreach (var it in dic2[n])
                        if (it.Key[t] == '0')
                            v++;
                    if(v>maxValue)
                    {
                        remIndex = n;
                        maxValue = v;
                    }
                }
                if (maxValue > sum)
                {
                    for (int i = 0; i < allStructures.Count;i++ )                  
                        profileKeys[0, i][t] = profileKeys[remIndex, i][t];                    
                }
            }
            
            for (int i = 0; i < profileKeys.GetLength(1);i++)
            {
                string key = new string(profileKeys[0, i]);
                if (newDic.ContainsKey(key))
                    newDic[key].Add(i);
                else
                {
                    List<int> aux=new List<int>();
                    aux.Add(i);
                    newDic.Add(key,aux);

                }
            }
            return newDic;
        }

        public ClusterOutput RunHashCluster()
        {
            ClusterOutput output = null;
            if (!input.floatingPointsProfiles)
            {
                List<string> lista = new List<string>(stateAlign.Keys);

                //ClusterOutput output=AutomaticCluster(lista, false);
                output = RunHashCluster(lista);
            }
            else
            {
                UserDefinedProfile<double> profiles = new UserDefinedProfile<double>();
                Dictionary<string,protInfo<double>> prof=profiles.GetProfile(null, this.alignFile, null);
                Dictionary<string, List<double>> data = new Dictionary<string, List<double>>();
                foreach (var item in prof)
                    data.Add(item.Key, item.Value.profile);

                List<string> keys = new List<string>(data.Keys);
                Dictionary<string,List<int>> cl=Rpart(null, data, keys, false);

                output= DictionaryToClusterOutput(cl,keys,false);

            }
            return output;
        }
        public ClusterOutput RunHashDendrog()
        {
            List<string> lista = new List<string>(stateAlign.Keys);

            
            
            ClusterOutput output = Cluster(lista,true);
            return output;
        }
        public Dictionary<string, List<int>> HashEntropyCombine(Dictionary<string, List<int>> dic, List<string> structures, int nodesNumber)
        {
            //dic = PrepareNewKeys(dic);
            List<KeyValuePair<int, double>> entropy = HashKeyEntropy(dic);

          

            int thresholdA = 0;
            int thresholdB = entropy.Count;
            int mPoint = (thresholdB + thresholdA) / 2;
            int iter = (int)Math.Ceiling(Math.Log(thresholdB, 2) / Math.Log(2, 2));
            int remCurrent = currentV;
            double step = 40 / iter;
            Dictionary<string, List<int>> localDic = new Dictionary<string, List<int>>(dic.Keys.Count);
            int iterNumber = 0;
            double prevStep = 0;
            do
            {
                //Create new keys
                localDic.Clear();
                foreach (var item in dic)
                {
                   
                    char[] tab = new char[mPoint];
                    for (int i = 0; i < mPoint; i++)
                        tab[i] = item.Key[entropy[i].Key];

                    string newKey = new string(tab);

                    if(!localDic.ContainsKey(newKey))                    
                        localDic.Add(newKey,new List<int>());

                    //foreach (var item in dic.Keys)
                     for (int i = 0; i < item.Value.Count; i++)
                            localDic[newKey].Add(item.Value[i]);                   


                }

               /* if(SelectionNMatchHeuristic(nodesNumber,0.3,localDic)==1)
                    thresholdA = mPoint;                                
                else
                    thresholdB = mPoint;*/
                if(localDic.Count<nodesNumber)               
                    thresholdA = mPoint;                                   
                else
                    if(localDic.Count>=nodesNumber)                    
                        thresholdB = mPoint;                                            
                   
                    
                mPoint = (thresholdB + thresholdA) / 2;
                iterNumber++;
                currentV += (int)(iterNumber * step - prevStep);
                prevStep = iterNumber * step;              
            }
            while((thresholdB-thresholdA)>2);

            currentV += 40 - (currentV - remCurrent);


            //ClusterOutput clOut = new ClusterOutput();

            return localDic;


        }

        public ClusterOutput RunHashCluster(List<string> structNames)
        {
            ClusterOutput outC;
            List<string> sData = new List<string>(structNames);
            // MakeBackgroundCluster();
            //outC = AutomaticCluster(structNames, false);
    
            outC = Cluster(sData);
            outC.runParameters = input.GetVitalParameters();
            return outC;


        }

        int FindNextRefPoint(List<KeyValuePair<int,int>> sortedDist,Dictionary <int,int>[] refPoints,List<int>refIndex)
        {
            FileStream dd = new FileStream("referenceDist.txt",FileMode.Append);
            StreamWriter wr = new StreamWriter(dd);
            if (refIndex.Count - 1 == 0)
            {
                wr.Close();
                dd.Close();
                return sortedDist[sortedDist.Count - 1].Value;
            }

            int[] newDist = new int[sortedDist.Count];

            for (int i = 0; i < sortedDist.Count; i++)
                newDist[sortedDist[i].Value] = sortedDist[i].Key;
            for (int i = 0; i < refIndex.Count - 1; i++)
            {                
                for (int j = 0; j < sortedDist.Count; j++)
                   newDist[j] += refPoints[i][j];
            }

            var sorted = newDist.Select((x, i) => new KeyValuePair<int, int>(x, i)).OrderBy(x => x.Key).ToList();
            wr.WriteLine("Next");
            foreach(var item in sorted)
            {
                wr.WriteLine(item.Key);
            }
            wr.Close();
            dd.Close();
            for(int i=0;i<refIndex.Count-1;i++)
              if(!refIndex.Contains(sorted[sorted.Count - i-1].Value))
                return sorted[sorted.Count - i-1].Value; ;

            return sorted[sorted.Count-1].Value;
            //return refPoints[refPoints.Length - 1][refPoints[0].Count - 1];
        }
        KeyValuePair<List<KeyValuePair<int, int>>, Dictionary<int, int>[]> GetRefDistances2(Dictionary<string, List<int>> dic, List<string> structures)
        {
            return new KeyValuePair<List<KeyValuePair<int, int>>, Dictionary<int, int>[]>();
        }
        KeyValuePair<List<KeyValuePair<int,int>>, Dictionary<int, int>[]> GetRefDistances(Dictionary<string, List<int>> dic, List<string> structures)
        {
            
               string[] hashKeys = new string[structures.Count];
            foreach (var item in dic)
            {
                for (int i = 0; i < item.Value.Count; i++)
                    hashKeys[item.Value[i]] = item.Key;
            }


            jury = new jury1D();
            List<KeyValuePair<string, double>> juryOut = jury.JuryOptWeights(structures).juryLike;
            int[] tabDist = new int[hashKeys.Length];

            for (int i = 0; i < hashKeys.Length; i++)
            {
                tabDist[i] = 0;
                for (int n = 0; n < hashKeys[i].Length; n++)
                    if ((int)hashKeys[i][n] != '0')
                        tabDist[i]++;
            }
            var sorted = tabDist.Select((x, i) => new KeyValuePair<int, int>(x, i)).OrderBy(x => x.Key).ToList();


            Dictionary<int, int>[] referenceDist = null;
            consensusStates = new Dictionary<byte, int>[columns.Length];
            for (int i = 0; i < columns.Length; i++)
                consensusStates[i] = new Dictionary<byte, int>();

            int stepRef = 0;

            if (refPoints > 0)
            {
                stepRef = sorted.Count / refPoints;
                referenceDist = new Dictionary<int, int>[refPoints];
            }
            for (int n = 0; n < refPoints; n++)
            {

                //string profileKey = hashKeys[sorted[stepRef * (n + 1) - 1].Value];
                string profileKey = juryOut[stepRef * (n + 1) - 1].Key;
                //   string profileKey= hashKeys[FindNextRefPoint(n,sorted,referenceDist)];
                List<byte> auxXX = auxDic[allStructures[dic[profileKey][0]]];
                for (int i = 0; i < auxXX.Count; i++)
                {
                    consensusStates[i].Clear();
                    consensusStates[i].Add(auxXX[i], 0);
                }


                Dictionary<string, List<int>> dic2 = PrepareKeys(allStructures, true, false);
                referenceDist[n] = new Dictionary<int, int>(hashKeys.Length);

                foreach (var item in dic2)
                {
                    int sumX = 0;
                    for (int k = 0; k < item.Key.Length; k++)
                        if ((int)item.Key[k] != '0')
                            sumX++;

                    foreach (var num in item.Value)
                        referenceDist[n].Add(num, sumX);


                }

            }
            return new KeyValuePair<List<KeyValuePair<int, int>>, Dictionary<int, int>[]>(sorted,referenceDist);
        }
        protected Dictionary<string, List<int>> FastCombineKeysNew(Dictionary<string, List<int>> dic, List<string> structures, bool flagDecision)
        {
            Dictionary<string, List<int>> outDic = null;

            if (dic.Count < input.relClusters)
                return dic;
            input.refPoints = 10;
            KeyValuePair<List<KeyValuePair<int, int>>, Dictionary<int, int>[]> dist = GetRefDistances(dic, structures);
            var sorted = dist.Key;
            var referenceDist = dist.Value;
            int thresholdA = 0;
            int thresholdB = sorted[sorted.Count - 1].Key;

            int iter = (int)Math.Ceiling(Math.Log(thresholdB, 2) / Math.Log(2, 2));
            int remCurrent = currentV;
            double step = 40 / iter;
            //int mPoint = (thresholdB+thresholdA)/2;
            int mPoint = sorted[sorted.Count / 2].Key;
            DebugClass.WriteMessage("rel Clusters" + input.relClusters);
            int xx = sorted.Count / input.relClusters;
            DebugClass.WriteMessage("xx=" + xx);
            /* for(int i=xx;i<sorted.Count;i++)
                 if(sorted[i].Key!=sorted[xx].Key)
                 {
                     mPoint=sorted[i].Key;
                     break;
                 }*/
            //mPoint=input.relClusters
            bool decision = false;
            int iterNumber = 0;
            double prevStep = 0;
            do
            {
                //Check clusters size
                DebugClass.WriteMessage("mPoint=" + mPoint);
                outDic = null;
                GC.Collect();
                outDic = ClustDistNew(dic, referenceDist, sorted, mPoint, 1);
                //outDic = ClustDist(dic,sorted, mPoint, 1);
                DebugClass.WriteMessage("dicSize=" + outDic.Count);
                //List<KeyValuePair<string, List<int>>> toSort = outDic.AsParallel().WithDegreeOfParallelism(s.numberOfCores).ToList();
                List<KeyValuePair<string, List<int>>> toSort = outDic.ToList();
                if (input.relClusters > toSort.Count)
                {
                    thresholdB = mPoint;
                    mPoint = (thresholdB + thresholdA) / 2;
                    iterNumber++;
                    currentV += (int)(iterNumber * step - prevStep);
                    prevStep = iterNumber * step;
                    continue;
                }
                DebugClass.WriteMessage("OK");
                if (flagDecision)
                {
                    toSort.Sort((nextPair, firstPair) => { return firstPair.Value.Count.CompareTo(nextPair.Value.Count); });
                    double sum = 0;
                    for (int i = 0; i < input.relClusters; i++)
                        sum += toSort[i].Value.Count;

                    sum /= structures.Count;

                    decision = sum >= input.perData / 100.0;
                    if (decision)
                        thresholdB = mPoint;
                    else
                        thresholdA = mPoint;
                }
                else
                {
                    decision = outDic.Keys.Count <= input.relClusters;
                    if (decision)
                        thresholdB = mPoint;
                    else
                        thresholdA = mPoint;

                }
                mPoint = (thresholdB + thresholdA) / 2;
                iterNumber++;
                currentV += (int)(iterNumber * step - prevStep);
                prevStep = iterNumber * step;
            }
            while ((thresholdB - thresholdA) > 1);

            currentV += 40 - (currentV - remCurrent);
            DebugClass.WriteMessage("NEXT");
            return outDic;
        }
        protected Dictionary<string, List<int>> FastCombineKeysOLD(Dictionary<string, List<int>> dic, List<string> structures, bool flagDecision)
        {            
            Dictionary<string, List<int>> outDic = null;
            List<string> hashKeys = new List<string>(dic.Keys);
            bool [] avoid = new bool  [hashKeys.Count];
              double sum=0;
            bool end = false;
            List<List<string>> clusters = new List<List<string>>();

            int[] tabDist = new int[hashKeys.Count];

            for (int i = 0; i < hashKeys.Count; i++)
            {
                tabDist[i] = 0;
                for (int n = 0; n < hashKeys[i].Length; n++)
                    if ((int)hashKeys[i][n] != '0')
                        tabDist[i]++;
            }
            var sorted = tabDist.Select((x, i) => new KeyValuePair<int, int>(x, i)).OrderBy(x => x.Key).ToList();

            Dictionary<int, int>[] referenceDist = null;
            consensusStates = new Dictionary<byte, int>[columns.Length];
            for (int i = 0; i < columns.Length; i++)                
                    consensusStates[i] = new Dictionary<byte, int>();

            int stepRef = 0;

            if (refPoints > 0)
            {
                stepRef = sorted.Count / refPoints;
                referenceDist = new Dictionary<int, int>[refPoints];
            }
            for (int n = 0; n < refPoints; n++)
            {                
                string profileKey = hashKeys[sorted[stepRef*(n+1) - 1].Value];
                List<byte> auxXX = auxDic[allStructures[dic[profileKey][0]]];
                for (int i = 0; i < auxXX.Count; i++)
                {
                    consensusStates[i].Clear();
                    consensusStates[i].Add(auxXX[i], 0);
                }


                Dictionary<string, List<int>> dic2 = PrepareKeys(allStructures, true, false);
                referenceDist[n] = new Dictionary<int, int>(hashKeys.Count);

                foreach (var item in dic2)
                {
                    int sumX = 0;
                    for (int k = 0; k < item.Key.Length; k++)
                        if ((int)item.Key[k] != '0')
                            sumX++;

                    foreach (var num in item.Value)
                        referenceDist[n].Add(num, sumX);


                }

            }
            int thresholdA = 0;
            int thresholdB = sorted[sorted.Count - 1].Key;

            int iter = (int)Math.Ceiling(Math.Log(thresholdB, 2) / Math.Log(2, 2));
            int remCurrent = currentV;
            double step = 40 / iter;
            //int mPoint = (thresholdB+thresholdA)/2;
            int mPoint = sorted[sorted.Count / 2].Key;
            DebugClass.WriteMessage("rel Clusters" + input.relClusters);
            int xx = tabDist.Length / input.relClusters;
            DebugClass.WriteMessage("xx=" + xx);
            /* for(int i=xx;i<sorted.Count;i++)
                 if(sorted[i].Key!=sorted[xx].Key)
                 {
                     mPoint=sorted[i].Key;
                     break;
                 }*/
            //mPoint=input.relClusters
            bool decision = false;
            int iterNumber = 0;
            double prevStep = 0;
            do
            {
                //Check clusters size
                DebugClass.WriteMessage("mPoint=" + mPoint);
                outDic = null;
                GC.Collect();
               // outDic = ClustDistNew(dic,referenceDist, sorted, mPoint, 1);
                outDic = ClustDist(dic,sorted, mPoint, 1);
                DebugClass.WriteMessage("dicSize=" + outDic.Count);
                //List<KeyValuePair<string, List<int>>> toSort = outDic.AsParallel().WithDegreeOfParallelism(s.numberOfCores).ToList();
                List<KeyValuePair<string, List<int>>> toSort = outDic.ToList();
                if (input.relClusters > toSort.Count)
                {
                    thresholdB = mPoint;
                    mPoint = (thresholdB + thresholdA) / 2;
                    iterNumber++;
                    currentV += (int)(iterNumber * step - prevStep);
                    prevStep = iterNumber * step;
                    continue;
                }

                if (flagDecision)
                {
                    toSort.Sort((nextPair, firstPair) => { return firstPair.Value.Count.CompareTo(nextPair.Value.Count); });
                    sum = 0;
                    for (int i = 0; i < input.relClusters; i++)
                        sum += toSort[i].Value.Count;

                    sum /= structures.Count;

                    decision = sum >= input.perData / 100.0;
                    if (decision)
                        thresholdB = mPoint;
                    else
                        thresholdA = mPoint;
                }
                else
                {
                    decision = outDic.Keys.Count <= input.relClusters;
                    if (decision)
                        thresholdB = mPoint;
                    else
                        thresholdA = mPoint;

                }
                mPoint = (thresholdB + thresholdA) / 2;
                iterNumber++;
                currentV += (int)(iterNumber * step - prevStep);
                prevStep = iterNumber * step;
            }
            while ((thresholdB - thresholdA) > 1);

            currentV += 40 - (currentV - remCurrent);

            return outDic;                        
        }
        private void CalcHammingDist(object o)
        {
            threadParam p = (threadParam)o;
            int num=p.threadNum;
            int begin=(int)((float)(allStructures.Count-1)/resetEvents.Length)*num;
            int end=(int)((float)(allStructures.Count-1)/resetEvents.Length)*(num+1);
            string x = allStructures[allStructures.Count - 1];
            for (int i = begin; i <end; i++)
            {
                int dist = 0;
                string y = allStructures[i];
                for (int n = 0; n < y.Length; n++)
                {
                    dist += x[n] ^ y[n];

                    if (dist > p.threshold)
                        break;
                }
                distTab[i] = dist;

            }
            resetEvents[num].Set();
        }
        int MoveToCorrectPosition(int reference,int start, int thresholdH,List<KeyValuePair<int, int>> sorted)
        {
            int val = sorted[reference].Key;
            int lastPosition = sorted.Count-1;
            int current;
            if (Math.Abs(val - sorted[start].Key) > thresholdH)
                return start;

            current = (start + lastPosition) / 2;

            while(lastPosition-start>1)
            {
                if (Math.Abs(val - sorted[current].Key) <= thresholdH)
                    start = current;
                else
                    lastPosition = current;
                current = (start + lastPosition) / 2;
            }
            if (Math.Abs(val - sorted[start].Key) > thresholdH)
                return start;
            
            return lastPosition;
        }
        private Dictionary<string, List<int>> ClustDistNew(Dictionary<string, List<int>> dic, Dictionary<int, int>[] referenceDist, List<KeyValuePair<int, int>> sorted, int thresholdH, int vT)
        {
            Dictionary<string, List<int>> outDic = new Dictionary<string, List<int>>(dic.Keys.Count);
            string[] hashKeys = new string[sorted.Count];

            foreach (var item in dic)
            {
                for (int i = 0; i < item.Value.Count; i++)
                    hashKeys[item.Value[i]] = item.Key;
            }



            bool[] avoid = new bool[sorted.Count];
            distTab = new int[sorted.Count];

            for (int i = 0; i < avoid.Length; i++)
                avoid[i] = false;


            List<int> candidates = new List<int>();

            for (int i = 0, k = 0; i < sorted.Count; i++)
            {
                if (avoid[sorted[i].Value])
                    continue;
                k = i + 1;
                if (k >= sorted.Count)
                    continue;

                while (k < sorted.Count && avoid[sorted[k].Value])
                    k++;
                if (k == sorted.Count)
                    k--;
                string keyDic = hashKeys[sorted[i].Value];
                int val = sorted[i].Key;
                candidates.Clear();
                //candidates.Add(sorted[i].Value);


                int newK = MoveToCorrectPosition(i, k, thresholdH, sorted);

                for (int j = k; j < newK; j++)
                    //if (!avoid[sorted[j].Value])
                    if (!avoid[sorted[j].Value] &&Math.Abs(referenceDist[referenceDist.Length - 1][sorted[j].Value] - referenceDist[referenceDist.Length-1][sorted[i].Value])<= thresholdH)
                            candidates.Add(sorted[j].Value);


                /* while (k < sorted.Count && Math.Abs(val - sorted[k].Key) <= thresholdH)
                 {
                     if (!avoid[sorted[k].Value])
                         candidates.Add(sorted[k].Value);
                     k++;
                 }*/
                if (candidates.Count > 0)
                {
                    int remThreadNumbers = threadNumbers;
                    if (referenceDist != null)
                    {
                        if (candidates.Count < 2 * remThreadNumbers)
                            threadNumbers = remThreadNumbers / 2;
                        if (threadNumbers < 1)
                            threadNumbers = 1;
                        for (int n = 0; n < threadNumbers; n++)
                        {
                            threadParam w=new threadParam();
                            w.threshold = thresholdH;
                            w.threadNum = n;
                            w.candidates = candidates;
                            w.referenceDist = referenceDist;
                            w.referenceIndex = sorted[i].Value;
                            resetEvents[n] = new ManualResetEvent(false);
                            ThreadPool.QueueUserWorkItem(new WaitCallback(CheckCandidates), (object)w);
                        }

                        for (int n = 0; n < threadNumbers; n++)
                            resetEvents[n].WaitOne();

                        candidates.Add(sorted[i].Value);

                        foreach (var cand in candidates)
                        {
                            if (cand >= 0)
                            {
                                List<int> inx = new List<int>();// (dic[hashKeys[sorted[k].Value]].Count);
                                if (avoid[cand])
                                    continue;
                                //Console.Write("Ups");
                                inx.AddRange(dic[hashKeys[cand]]);
                                foreach (var item in dic[hashKeys[cand]])
                                    avoid[item] = true;
                                if (!outDic.ContainsKey(keyDic))
                                {
                                    List<int> auxList = new List<int>();
                                    //auxList.AddRange(dic[keyProfile]);
                                    auxList.AddRange(inx);
                                    outDic.Add(keyDic, auxList);
                                }
                                else
                                    outDic[keyDic].AddRange(inx);
                                avoid[cand] = true;
                            }
                        }
                    }
                    threadNumbers = remThreadNumbers;
                }
               

            }
            for (int i = 0; i < sorted.Count; i++)
                if (avoid[sorted[i].Value] == false)
                    if (!outDic.ContainsKey(hashKeys[sorted[i].Value]))
                    {
                        outDic.Add(hashKeys[sorted[i].Value], dic[hashKeys[sorted[i].Value]]);
                    }


            return outDic;
        }
        void CheckCandidates(object o)//Dictionary<int, int>[] referenceDist, List<int> candidates, int referanceIndex,int threshold,int threadNum)
        {
            threadParam p = (threadParam)o;
            Dictionary<int, int>[] referenceDist = p.referenceDist;
            List<int> candidates = p.candidates;
            int referanceIndex = p.referenceIndex;
            int threadNum = p.threadNum;
            int threshold = (int)p.threshold;

            double step = ((double)candidates.Count) / threadNumbers;
            int start = (int)(step * threadNum);
            int end = (int)(step * (threadNum + 1));
            if (end > candidates.Count)
                end = candidates.Count;
            for (int x = start; x < end; x++)
            {
                bool test = true;
                for (int n = referenceDist.Length-1; n >=0 ; n--)
                    if (Math.Abs(referenceDist[n][candidates[x]] - referenceDist[n][referanceIndex]) > threshold)
                    {
                        test = false;
                        break;
                    }
                if (!test)
                    lock (candidates)
                        candidates[x] = -1;
            }
            resetEvents[threadNum].Set();
        }

            private Dictionary<string, List<int>> ClustDist(Dictionary<string, List<int>> dic, List<KeyValuePair<int, int>> sorted, int thresholdH, int vT)
        {
            Dictionary<string, List<int>> outDic = new Dictionary<string, List<int>>(dic.Keys.Count);
            List<string> hashKeys = new List<string>(dic.Keys);
            bool[] avoid = new bool[hashKeys.Count];
            List<int> index = new List<int>(hashKeys.Count);
            distTab = new int[hashKeys.Count];
          
            for (int i = 0; i < avoid.Length; i++)
                avoid[i] = false;

            resetEvents = new ManualResetEvent[threadNumbers];
            for (int i = 0; i < resetEvents.Length; i++)
                resetEvents[i] = new ManualResetEvent(false);

            List<int> aux = new List<int>(hashKeys.Count);

            //float threshold = thresholdH;
           // for (float threshold = (float)thresholdH / vT; threshold <= thresholdH;)
            {
                
                for (int i = 0, k = 0; i < hashKeys.Count; i++)
                {
                    if (avoid[sorted[i].Value])
                        continue;


                    aux.Clear();
                    if (outDic.ContainsKey(hashKeys[sorted[i].Value]))
                        //foreach (var item in outDic[hashKeys[sorted[i].Value]])
                        aux.AddRange(outDic[hashKeys[sorted[i].Value]]);
                    else
                        aux.AddRange(dic[hashKeys[sorted[i].Value]]);

                    k = i + 1;
                    allStructures.Clear();
                    index.Clear();
                    int val = sorted[i].Key;
                    while (k < sorted.Count && Math.Abs(val - sorted[k].Key) <= thresholdH / 4 + 1)
                    {
                        if (avoid[sorted[k].Value])
                        {
                            k++;
                            continue;
                        }
                        allStructures.Add(hashKeys[sorted[k].Value]);
                        index.Add(sorted[k].Value);
                        k++;
                    }
                    allStructures.Add(hashKeys[sorted[i].Value]);
                    for (int n = 0; n <threadNumbers; n++)
                    {
                        threadParam w=new threadParam();
                        w.threshold = thresholdH;
                        w.threadNum = n;
                        resetEvents[n].Reset();
                        //resetEvents[n] = new ManualResetEvent(false);
                        ThreadPool.QueueUserWorkItem(new WaitCallback(CalcHammingDist), (object)w);
                    }

                    for (int n = 0; n < resetEvents.Length; n++)
                        resetEvents[n].WaitOne();

                   // Parallel.For(0, allStructures.Count - 2, n =>
                    for (int n = 0; n < allStructures.Count - 1; n++)
                    {
                        if (distTab[n] <= thresholdH)
                        {
                            string structName = allStructures[n];
                            if (outDic.ContainsKey(structName))
                            {
                                aux.AddRange(outDic[structName]);
                                outDic.Remove(structName);
                                //foreach (var item in outDic[allStructures[n]])
                                //  aux.Add(item);
                            }
                            else
                                aux.AddRange(dic[structName]);
                            //foreach (var item in dic[allStructures[n]])
                            //  aux.Add(item);

                            avoid[index[n]] = true;
                        }

                    }//);

                    if (aux.Count > 1)
                    {
                        if (!outDic.ContainsKey(hashKeys[sorted[i].Value]))
                        {
                            outDic.Add(hashKeys[sorted[i].Value],new List<int>(aux));
                        }
                    }
                }
               // threshold += (float)thresholdH / vT;
            }
            for (int i = 0; i < sorted.Count;i++)
                if (avoid[sorted[i].Value]==false )
                    if (!outDic.ContainsKey(hashKeys[sorted[i].Value]))
                    {
                        outDic.Add(hashKeys[sorted[i].Value], dic[hashKeys[sorted[i].Value]]);
                    }


            return outDic;
        }
        Dictionary<string,List<int>> PrepareRefDictionaries(List<string> structures)
        {
            Dictionary<string, List<int>> listDic = new Dictionary<string, List<int>>();
            Dictionary<int, string> aux = new Dictionary<int, string>();

            Dictionary<string,List<int>> dic=PrepareKeys(structures, true, true);
            string[] hashKeys = new string[structures.Count];

            foreach(var item in dic)
            {
                foreach(var intV in item.Value)
                {
                    aux.Add(intV, item.Key);
                }
            }


            jury = new jury1D();
            jury.PrepareJury(al, stateAlign);
            List<KeyValuePair<string, double>> juryOut = jury.JuryOptWeights(structures).juryLike;



            foreach (var item in dic)
            {
                for (int i = 0; i < item.Value.Count; i++)
                    hashKeys[item.Value[i]] = item.Key;
            }

            List<List<string>> clusters = new List<List<string>>();


            consensusStates = new Dictionary<byte, int>[columns.Length];
            for (int i = 0; i < columns.Length; i++)
                consensusStates[i] = new Dictionary<byte, int>();

            for (int n = 0; n < refPoints; n++)
            {
                string profileKey = juryOut[(int)((n + 1) * 0.75*juryOut.Count /refPoints)].Key;
                List<byte> auxXX = auxDic[profileKey];
                for (int i = 0; i < auxXX.Count; i++)
                {
                    consensusStates[i].Clear();
                    consensusStates[i].Add(auxXX[i], 0);
                }


                Dictionary<string, List<int>> dic2 = PrepareKeys(allStructures, false, false);
                foreach (var item in dic2)
                    foreach (var intV in item.Value)
                        aux[intV] += item.Key;
            }
            foreach (var item in aux)
                if (listDic.ContainsKey(item.Value))
                    listDic[item.Value].Add(item.Key);
                else
                {
                    List<int> l = new List<int>();
                    l.Add(item.Key);
                    listDic.Add(item.Value, l);
                }

            return listDic;

        }
        KeyValuePair<List<KeyValuePair<int, int>>, Dictionary<int, int>[]> ReferenceDistances(Dictionary<string, List<double>> dic, List<string> structures)
        {
            List<string> keys = new List<string>(dic.Keys);
            int[] tabDist = new int[keys.Count];

            Random r = new Random();

            int index = r.Next(keys.Count);

            for (int i = 0; i < keys.Count; i++)
            {
                double sum = 0;
                
                for (int n = 0; n < dic[keys[i]].Count; n++)
                    sum += Math.Pow(dic[keys[index]][n] - dic[keys[i]][n], 2);
                tabDist[i] = (int)Math.Sqrt(sum);
            }
            var sorted = tabDist.Select((x, i) => new KeyValuePair<int, int>(x, i)).OrderBy(x => x.Key).ToList();

            Dictionary<int, int>[] referenceDist = null;
            int stepRef = 0;

            if (refPoints > 0)
            {
                stepRef = (int)(sorted.Count * 3 / 4.0 / refPoints);
                referenceDist = new Dictionary<int, int>[refPoints];
            }
            List<int> refIndex = new List<int>();
            refIndex.Add(sorted[0].Value);
            for (int n = 0; n < refPoints; n++)
            {

                //   string profileKey = hashKeys[sorted[stepRef * (n + 1) - 1].Value];
                refIndex.Add(FindNextRefPoint(sorted, referenceDist, refIndex));
                referenceDist[n] = new Dictionary<int, int>(sorted.Count);

                for(int s=0;s< keys.Count;s++)
                {
                    double sumX = 0;
                    List<double> aux = dic[keys[refIndex[refIndex.Count-1]]];
                    for (int k = 0; k < dic[keys[s]].Count; k++)
                    {
                        double v = aux[k] - dic[keys[s]][k];
                        sumX += v * v;
                    }
                    
                    referenceDist[n].Add(s, (int)Math.Sqrt(sumX));


                }

            }
            return new KeyValuePair<List<KeyValuePair<int, int>>, Dictionary<int, int>[]>(sorted, referenceDist);
        }
        KeyValuePair<List<KeyValuePair<int,int>>, Dictionary<int,int>[]> ReferenceDistances(Dictionary<string, List<int>> dic, List<string> structures)
        {
            string[] hashKeys = new string[structures.Count];
            

            jury = new jury1D();
            jury.PrepareJury(al, stateAlign);

            /*  List<string> dicStruct = new List<string>();
              foreach (var item in dic)
                  dicStruct.Add(structures[item.Value[0]]);
              List<KeyValuePair<string, double>> juryOut = jury.JuryOptWeights(dicStruct).juryLike;*/
            List<KeyValuePair<string, double>> juryOut = jury.JuryOptWeights(structures).juryLike;

            foreach (var item in dic)
            {
                for (int i = 0; i < item.Value.Count; i++)
                    hashKeys[item.Value[i]] = item.Key;
            }
            
            List<List<string>> clusters = new List<List<string>>();

            int[] tabDist = new int[hashKeys.Length];

            for (int i = 0; i < hashKeys.Length; i++)
            {
                tabDist[i] = 0;
                for (int n = 0; n < hashKeys[i].Length; n++)
                    if ((int)hashKeys[i][n] != '0')
                        tabDist[i]++;
            }
            var sorted = tabDist.Select((x, i) => new KeyValuePair<int, int>(x, i)).OrderBy(x => x.Key).ToList();

            Dictionary<int, int>[] referenceDist = null;
            consensusStates = new Dictionary<byte, int>[columns.Length];
            for (int i = 0; i < columns.Length; i++)
                consensusStates[i] = new Dictionary<byte, int>();

            int stepRef = 0;

            if (refPoints > 0)
            {
                stepRef = (int)(sorted.Count * 3 / 4.0 / refPoints);
                referenceDist = new Dictionary<int, int>[refPoints];
            }
            List<int> refIndex = new List<int>();
            refIndex.Add(sorted[0].Value);
            for (int n = 0; n < refPoints; n++)
            {

                //   string profileKey = hashKeys[sorted[stepRef * (n + 1) - 1].Value];
                string profileKey = juryOut[stepRef * (n + 1) - 1].Key;
                //refIndex.Add(FindNextRefPoint(sorted, referenceDist,refIndex));
                //string profileKey = hashKeys[refIndex[refIndex.Count-1]];
                //List<byte> auxXX = auxDic[allStructures[dic[profileKey][0]]];
                List<byte> auxXX = auxDic[profileKey];
                for (int i = 0; i < auxXX.Count; i++)
                {
                    consensusStates[i].Clear();
                    consensusStates[i].Add(auxXX[i], 0);
                }


                Dictionary<string, List<int>> dic2 = PrepareKeys(allStructures, false, false);
                referenceDist[n] = new Dictionary<int, int>(hashKeys.Length);

                foreach (var item in dic2)
                {
                    int sumX = 0;
                    for (int k = 0; k < item.Key.Length; k++)
                        if ((int)item.Key[k] != '0')
                            sumX++;

                    foreach (var num in item.Value)
                        referenceDist[n].Add(num, sumX);


                }
                dic2.Clear();
                dic2 = null;

            }

            return new KeyValuePair<List<KeyValuePair<int, int>>, Dictionary<int, int>[]> (sorted, referenceDist);
        }
        /*void MinStates(Dictionary<string, List<int>> dic, List<string> structures)
        {
            Dictionary<string, List<int>>[] hashDataBase = new Dictionary<byte, List<int>>[stateAlign[structures[0]].Count];

            for (int i = 0; i < hashDataBase.Length; i++)
                hashDataBase[i] = new Dictionary<byte, List<int>>();


            for(int j=0;j<structures.Count;j++)
            {
                for(int i=0; i<stateAlign[structures[j]].Count;i++)
                {
                    if (hashDataBase[i].ContainsKey(stateAlign[structures[j]][i]))
                        hashDataBase[i][stateAlign[structures[j]][i]].Add(j);
                    else
                    {
                        List<int> k = new List<int>();
                        k.Add(j);
                        hashDataBase[i].Add(stateAlign[structures[j]][i],k);
                    }
                }
            }
            for (int j = 0; j < hashDataBase.Length; j++)
            {
                byte remKey =0;
                int min = int.MaxValue;
                Dictionary<byte, List<int>> locHash = hashDataBase[j];
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
            // binSizeG = 1;           
            
            //List<string> triples = GenerateTriplets();

        }*/
        protected Dictionary<string, List<int>> Rpart(Dictionary<string, List<int>> dic, Dictionary<string, List<double>> dicDouble,List<string> structures, bool flagDecision)
        {
            Dictionary<string, List<int>> outDic = null;

            if (input.relClusters > dic.Count)
                return dic;

            KeyValuePair<List<KeyValuePair<int, int>>, Dictionary<int, int>[]> refValue;
            if (dic != null)
                refValue = ReferenceDistances(dic, structures);
            else
            {
                refValue = ReferenceDistances(dicDouble, structures);
                dic = new Dictionary<string, List<int>>();
                for (int i = 0; i < structures.Count; i++)
                {
                    List<int> x = new List<int>();
                    x.Add(i);
                    dic.Add(structures[i], x);
                }
                if (resetEvents == null)
                    resetEvents = new ManualResetEvent[threadNumbers];
            }

            Debug.Flush();
            List<KeyValuePair<int, int>> sorted = refValue.Key;
            Dictionary<int, int>[] referenceDist=refValue.Value;
            int thresholdA = 0;
            int thresholdB = sorted[sorted.Count - 1].Key;
            DebugClass.WriteMessage("mPointOK");

            int iter = (int)Math.Ceiling(Math.Log(thresholdB, 2) / Math.Log(2, 2));
            int remCurrent = currentV;
            double step = 40 / iter;
            //int mPoint = (thresholdB+thresholdA)/2;
            int mPoint = sorted[sorted.Count / 2].Key;
            DebugClass.WriteMessage("rel Clusters" + input.relClusters);
            /* for(int i=xx;i<sorted.Count;i++)
                 if(sorted[i].Key!=sorted[xx].Key)
                 {
                     mPoint=sorted[i].Key;
                     break;
                 }*/
            //mPoint=input.relClusters
            bool decision = false;
            int iterNumber = 0;
            double prevStep = 0;
            Dictionary<string, List<int>> bestRes = null;
            StreamWriter tt = new StreamWriter("zeroFile.txt");
            do
            {
                //Check clusters size
                DebugClass.WriteMessage("mPoint=" + mPoint);
                outDic = null;
                GC.Collect();
                outDic = ClustDistNew(dic, referenceDist, sorted, mPoint, 1);
                //Checking single clusters
                int zeroCount = 0;
                foreach (var xx in outDic)
                    if (xx.Value.Count == 1)
                        zeroCount++;
                tt.WriteLine("size=" + outDic.Count + " single=" + zeroCount);


                //=================================
                if (bestRes == null)
                    bestRes = outDic;
                else
                {
                    DebugClass.WriteMessage("BEST");
                    if (Math.Abs(outDic.Count - input.relClusters) < Math.Abs(outDic.Count - bestRes.Count))
                        bestRes = outDic;
                }
                //outDic = ClustDist(dic,sorted, mPoint, 1);
                DebugClass.WriteMessage("dicSize=" + outDic.Count);
                //List<KeyValuePair<string, List<int>>> toSort = outDic.AsParallel().WithDegreeOfParallelism(s.numberOfCores).ToList();
                List<KeyValuePair<string, List<int>>> toSort = outDic.ToList();
                if (input.relClusters > toSort.Count)
                {
                    thresholdB = mPoint;
                    mPoint = (thresholdB + thresholdA) / 2;
                    iterNumber++;
                    currentV += (int)(iterNumber * step - prevStep);
                    prevStep = iterNumber * step;
                    continue;
                }
                DebugClass.WriteMessage("OK1");
                if (flagDecision)
                {
                    toSort.Sort((nextPair, firstPair) => { return firstPair.Value.Count.CompareTo(nextPair.Value.Count); });
                    double sum = 0;
                    for (int i = 0; i < input.relClusters; i++)
                        sum += toSort[i].Value.Count;

                    sum /= structures.Count;

                    decision = sum >= input.perData / 100.0;
                    if (decision)
                        thresholdB = mPoint;
                    else
                        thresholdA = mPoint;
                }
                else
                {
                    decision = outDic.Keys.Count <= input.relClusters;
                    if (decision)
                        thresholdB = mPoint;
                    else
                        thresholdA = mPoint;

                }
                DebugClass.WriteMessage("OK2");
                mPoint = (thresholdB + thresholdA) / 2;
                iterNumber++;
                currentV += (int)(iterNumber * step - prevStep);
                prevStep = iterNumber * step;
              
            }
            while ((thresholdB - thresholdA) > 1);

            currentV += 40 - (currentV - remCurrent);
            DebugClass.WriteMessage("NEXT");
            tt.Close();

            //return outDic;
            return bestRes;
        }

        private List<List<string>> CombineKeys(Dictionary<string, List<int>> dic, List<string> structures)
        {
            List<List<string>> lKeys = new List<List<string>>();
            Dictionary <string,int> avoid=new Dictionary<string,int>();
            List<string> hashKeys = new List<string>();
            List<string> dd = new List<string>();

            var oList = dic.OrderBy(x => x.Value.Count).ToList();
            for (int i = oList.Count - 1; i >= 0; i--)
                hashKeys.Add(oList[i].Key);

            int [] iTab = new int[hashKeys[0].Length];
            for (int i = 0; i < hashKeys.Count; i++)
            {
                List<string> aux = new List<string>();
                if (avoid.ContainsKey(hashKeys[i]))
                    continue;

                aux.Add(hashKeys[i]);
                for (int n = 0; n < hashKeys[i].Length; n++)
                    iTab[n] = (int)hashKeys[i][n];
                for (int j = i+1; j < hashKeys.Count; j++)
                {
                    if (avoid.ContainsKey(hashKeys[j]))
                        continue;
                    int sum = 0;
                    for (int n = 0; n < hashKeys[j].Length; n++)
                        if (iTab[n] != (int)hashKeys[j][n])
                            sum++;

                    if (sum <= 1)
                    {
                        aux.Add(hashKeys[j]);
                        //if(!avoid.ContainsKey(hashKeys[j]))
                            avoid.Add(hashKeys[j],0);
                    }
                }
                lKeys.Add(aux);
            }

            List<List<string>> clusters = new List<List<string>>();
            for (int i = 0; i < lKeys.Count; i++)
            {
                clusters.Add(new List<string>());
                for (int j = 0; j < lKeys[i].Count; j++)
                {
                    foreach (var item in dic[lKeys[i][j]])
                        clusters[clusters.Count - 1].Add(structures[item]);
                }
            }

            return clusters;
        }
        void CheckRPart(Dictionary<string, List<int>> dic, List<string> structures)
        {
            Dictionary<string, List<int>> res;
            int prev = 0;
            for (int i = 3000; i < 20000; i += 1000)
            {
                input.relClusters = i;
                res = Rpart(dicC, null,structNames, true);
                if (res.Count != prev)
                {
                    DebugClass.WriteMessage("output");
                    ClusterOutput oo = DictionaryToClusterOutput(res, structures);
                    DebugClass.WriteMessage("Prepared");
                    string fileName = "Rpart_scope_ref200_" + i;
                    oo.SaveTxt(fileName);
                    DebugClass.WriteMessage("done");
                    prev = res.Count;
                }
            }
        }
        private ClusterOutput PrepareClusters(Dictionary<string, List<int>> dic, List<string> structures)
        {
            List<List<string>> clusters = new List<List<string>>(dic.Count);
            validIndexes.Clear();
            if (!input.rpart)
            {
                switch (input.selectionMethod)
                {
                    case COL_SELECTION.ENTROPY:
                        refPoints =10;
                        // dicC=PrepareRefDictionaries(structures);                        
                        //AddColumnsOneByOne(dicC, structures);
                       // dicC = DropColumnsOneByOne(dicC,structures);
                        //dicC = SelectColumnsByEntropy(dicC, input.relClusters, input.perData / 100.0);
                        break;
                    case COL_SELECTION.META_COL:
                        dicC = MetaColumns(structNames, input.relClusters, input.perData / 100.0);
                        break;
                    case COL_SELECTION.KOLMSMIR:
                        dicC = SelectColumnsByEntropy(dicC, input.relClusters*3, input.perData / 100.0);                        
                        //dicC = CheckSeparation(dicC);
                        dicC = Rpart(dicC,null, structNames, true);
                        dic = dicFinal;
                        break;
                }
                currentV+=40;
                structToKey.Clear();
                foreach (var item in dicC)
                {
                    foreach (var d in item.Value)
                        structToKey.Add(structNames[d], item.Key);
                }

                dic = dicC;
                dicFinal = dic;
                DebugClass.WriteMessage("Sorted");
            }
            else
            {
                /*FastRetrival fast = new FastRetrival();
                structToKey.Clear();*/
                //fast.dicC = dic;
                //fast.structToKey = structToKey;
                /*fast.RunFastRetrival(allStructures, allStructures, 50);*/

                DebugClass.WriteMessage("before FastCombine: " + dicC.Keys.Count);
                    // dicC = HashEntropyCombine(dicC, structures, dicC.Count);
                    /*StreamWriter wr = new StreamWriter("keys.txt");
                    foreach(var item in dicC)
                    {
                        foreach(var li in item.Value)
                        {
                            wr.WriteLine(structures[li] + " " + item.Key);
                        }
                    }

                    wr.Close();*/
                   dicFinal = Rpart(dicC,null, structNames, true);
                    //CheckRPart(dicC, structNames);
                    // dicFinal = FastCombineKeysNew(dicC, structNames, true);
                    dic = dicFinal;
                
                currentV++;
                DebugClass.WriteMessage("After FastCombine: " + dicC.Keys.Count);

            }          

            return DictionaryToClusterOutput(dicFinal, structures);
        }

        private ClusterOutput DictionaryToClusterOutput(Dictionary<string,List<int>> dic, List<string> structures,bool findRepr=true)
        {
            List<List<string>> clusters = new List<List<string>>(dic.Count);

            if (findRepr)
            {
                jury = new jury1D();
                jury.PrepareJury(al);
            }
            foreach (var item in dic.Keys)
            //foreach (var item in keyList)
            {
                List<string> aux = new List<string>(dic[item].Count);

                foreach (var st in dic[item])
                    aux.Add(structures[st]);

                //Added -------------------------------------------------------------
                if (findRepr)
                {
                    ClusterOutput outJury = jury.JuryOptWeights(aux);

                    int index = aux.FindIndex(x => x == outJury.juryLike[0].Key);

                    if (index != 0)
                    {
                        string tmp = aux[0];
                        aux[0] = aux[index];
                        aux[index] = tmp;
                    }
                    //----------------------------------------------------------------------
                    
                }
                clusters.Add(aux);
            }
            ClusterOutput res = new ClusterOutput();
            res.clusters = new clusterRes();
            res.clusters.list = clusters;

            return res;
        }
        protected string Regularization(string key, string structName)
        {
            if (!stateAlignReg.ContainsKey(structName))
                return null;
            List<byte> states = stateAlignReg[structName];
            string newKey = "";
            for (int i = 0; i < key.Length; i++)
            {
                int start = i - input.wSize / 2;
                if (i < input.wSize / 2)
                    start = i;
                if (start + input.wSize >= key.Length)
                    start = key.Length - input.wSize;
                int dist = 0;
                for (int n = start; n < start + input.wSize; n++)
                    if (consensusStatesReg[n].ContainsKey(states[n]))
                    {
                        if (consensusStatesReg[n][states[n]] != 0)
                            dist++;
                    }
                    else
                        dist++;
                if (dist < input.regThreshold)
                    newKey += "0";
                else
                    newKey += key[i];
            }
            return newKey;

        }
        public string TransformToConsensusStates(List<byte> data,int k)
        {
            string key;
            int i = 0;

            if (consensusProjection && consensusStates.Length != data.Count)
                return null;

            foreach (var item in data)
            {
                if (consensusProjection)
                    if (consensusStates[i].ContainsKey(item))
                        if (consensusStates[i][item] == 0) //if agree with consensus then 0
                            keyTab[k][i++] = '0';
                        else
                            keyTab[k][i++] = '1';

                    else
                        keyTab[k][i++] = '1';
                else
                    keyTab[k][i++] = (char)(item+'0');
            }
            key = new string(keyTab[k]);
            return key;

        }
        protected string TransformToConsensusStates(string structName, int k)
        {            
            if (!stateAlign.ContainsKey(structName))
                return null;
            return TransformToConsensusStates( stateAlign[structName],k);
        }


        private void PrepareKeysThreading(object k)
        {
            string hkey;
            foreach (var item in threadStructures[(int)k])
            {
                 hkey=TransformToConsensusStates(allStructures[item],(int)k);
                if (hkey==null)
                    continue;
                if (input.regular)
                    hkey = Regularization(hkey, allStructures[item]);
                if (hkey==null)
                    continue;
                lock (threadingKeys)
                {
                    if (!threadingKeys.ContainsKey(hkey))
                        threadingKeys.Add(hkey, new List<int>());
                    threadingKeys[hkey].Add(item);
                }
            }
            resetEvents[(int)k].Set();
        }
        void FindConsensusState(List<string> structures)
        {
            columns = MakeColumnsLists(structures, stateAlign);

            if (input.regular)
            {
                columnsReg = MakeColumnsLists(structures, stateAlignReg);
                consensusStatesReg = HammingConsensusStates(columnsReg);
            }
            if (input.jury)
                consensusStates = JuryConsensusStates(al, structures);
            else
                consensusStates = HammingConsensusStates(columns);

            //columns = null;
            //  al.r.profiles.Clear();
            GC.Collect();

        }
        public Dictionary<string, List<int>> PrepareKeys(List<string> structures, bool cleanMemory = true,bool consensus=true)
        {
            if(consensus)
                FindConsensusState(structures);
            else
                columns = MakeColumnsLists(structures, stateAlign);
            threadingKeys = new Dictionary<string, List<int>>(structures.Count);
            allStructures = new List<string>(structures);
            int localThreadN = threadNumbers;

            threadStructures = new List<int>[localThreadN];
            int part = structures.Count / localThreadN + 1;
            for (int i = 0; i < localThreadN; i++)
                threadStructures[i] = new List<int>(part);

            int count=0;
            for (int i = 0; i < structures.Count; i++)
            {
                threadStructures[count].Add(i);
                if (part <= threadStructures[count].Count)
                    count++;
            }
            resetEvents = new ManualResetEvent[localThreadN];

            currentV += 10;

            List<string> aux = new List<string>(stateAlign.Keys);
            //threadingKeys = new Dictionary<string, List<int>>();
            keyTab = new char[localThreadN][];
            for (int i = 0; i < localThreadN; i++)
                keyTab[i]=new char [stateAlign[aux[0]].Count];
           
            for (int n = 0; n <localThreadN ; n++)
            {
                int k=n;
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(PrepareKeysThreading),(object) k);               
            }
            for (int i = 0; i < localThreadN; i++)
                resetEvents[i].WaitOne();


            currentV += 10;
                if (cleanMemory)
                {
                    //stateAlign.Clear();

                    //al.Clean();
                    GC.Collect();
                }
            

            return threadingKeys;
        }
        public double CalculateDaviesBouldinIndex(List<List<string>> clusters)
        {
            //ClusterOutput clustOut;
            float[] avr = new float[clusters.Count];
            double measure = 0;
            List<string> refStructures = new List<string>();
            ClusterOutput juryOut = null;

            for (int n = 0; n < clusters.Count; n++)
            {
                juryOut = jury.JuryOptWeights(clusters[n]);
                refStructures.Add(structToKey[juryOut.juryLike[0].Key]);

            }

            for (int i = 0; i < clusters.Count; i++)
            {
                float sum = 0;
                if (clusters[i].Count == 0)
                    continue;
                
                foreach(var item in clusters[i])
                {
                    int d = 0;
                    string w = structToKey[item];
                    for (int n = 0; n < refStructures[i].Length; n++)
                        if (refStructures[i][n] != w[n])
                            d++;
                    sum += d;
                }

              

                avr[i] = sum / clusters[i].Count;

            }
            for (int i = 0; i < refStructures.Count; i++)
            {
                float max = 0;
                for (int j = 0; j < refStructures.Count; j++)
                {
                    int cDist=0;
                    float v;
                    if (i == j)
                        continue;

                    for (int n = 0; n < refStructures[i].Length; n++)
                        if (refStructures[i][n] != refStructures[j][n])
                            cDist++;

                    v = ((float)(avr[i] + avr[j])) / cDist;
                    if (v > max)
                        max = v;
                }
                measure += max;
            }

            return measure / refStructures.Count;


        }
        private double CalcDisp(List<List<string>> clusters,Dictionary<string,string> dic)
        {
            int size = 0;
            List<double> clustDisp = new List<double>();
            double finalRes = 0;
            foreach (var item in clusters)
            {
                size = 0;
                int dist = 0;

                for (int i = 0; i < item.Count;i++ )
                {
                    if (!dic.ContainsKey(item[i]))
                        continue;

                    string key1 = dic[item[i]];
                    for (int j = i + 1; j < item.Count;j++)
                    {
                        if (!dic.ContainsKey(item[j]))
                            continue;
                        string key2 = dic[item[j]];

                        for (int n = 0; n < key2.Length; n++)
                            dist += key1[n] ^ key2[n];
                        size += key2.Length;
                    }
                }
                if(size>0)
                    finalRes+=((double)dist) / size;
            }
            finalRes /= clusters.Count;

            return finalRes;

        }
        private double EvalClustering(List<List<string>> clusters)
        {
            ClusterOutput juryOut = null ;

            clusters.Sort(
            delegate(List<string> p1, List<string> p2)
            {
                return p2.Count.CompareTo(p1.Count);
            }
            );
            int avrDist = 0;

            List<string> refStructures=new List<string>();
            for (int n = 0; n < clusters.Count;n++ )
            {
                juryOut = jury.JuryOptWeights(clusters[n]);
                refStructures.Add(structToKey[juryOut.juryLike[0].Key]);

            }
            int counter=0;
            for (int n = 0; n < input.relClusters;n++ )
            {

                for (int m = n+1; m < input.relClusters;m++ )
                {
                    int dist = 0;
                    for (int i = 0; i < refStructures[n].Length; i++)
                        if (refStructures[m][i] != refStructures[n][i])
                            dist++;
                    avrDist += dist;
                    counter++;
                }
            }
            avrDist /= counter;
            int countRadStr = 0;
            for (int n = 0; n < clusters.Count;n++ )
            {            
                string refStr = refStructures[n];                
                foreach (var item in clusters[n])
                {
                    int dist = 0;
                    string str = structToKey[item];
                    for (int i = 0; i < str.Length; i++)
                        if (str[i] != refStr[i])
                            dist++;
                    if (dist < avrDist * 0.25)
                        countRadStr++;
                }

            }
            counter = 0;
            for (int n = 0; n < input.relClusters; n++)
                counter += clusters[n].Count;
            double goodness = (double)countRadStr / counter;

            goodness -= (1-input.perData / 100.0) * input.relClusters;
            
            
            return goodness;
        }
        private Dictionary<string,List<byte>> RandomCoverProfiles(int num,int keyLength,Dictionary<byte,double>stateStats)// Done only for contactMap must be generalized to any profile
        {
            Random r = new Random();
            Dictionary<string, List<byte>> profiles = new Dictionary<string, List<byte>>();
            Dictionary<string,double> localStatesStat=new Dictionary<string,double>(stateStats.Count);
            List<byte> profile;
            byte[] statesProb = new byte[100];

            int s=0;
            foreach(var item in stateStats)
            {
                int remS = s;
                for (; s < remS+item.Value * 100 && s < statesProb.Length; s++)
                    statesProb[s] = item.Key;
            }
            for(int i=0;i<num;i++)
            {
                profile = new List<byte>(keyLength);
                for (int j = 0; j < keyLength;j++ )
                {
                    int index = r.Next(100);
                    profile.Add(statesProb[index]);
                }
                string name="Ex_"+i;                

                profiles.Add(name,profile);
                    
            }

            return profiles;
        }
        KeyValuePair<Dictionary<string,List<int>>,List<string>> MakeBackgroundData(List<string> data)
        {
            GStatisticsData gData = new GStatisticsData(data);
            Dictionary<string,string> dicG=gData.GenerateData();
            List<string> bStructures=new List<string>(dicG.Keys);
            Dictionary<string,int> nameToIndex=new Dictionary<string,int>();
            for(int i=0;i<bStructures.Count;i++)
                nameToIndex.Add(bStructures[i],i);


            Dictionary<string,List<int>> localDic=new Dictionary<string,List<int>>();
            foreach(var item in dicG)
            {
                if(localDic.ContainsKey(item.Value))
                    localDic[item.Value].Add(nameToIndex[item.Key]);
                else
                {
                    List <int> index=new List<int>();
                    index.Add(nameToIndex[item.Key]);
                    localDic.Add(item.Value, index);
                }
            }
            return new KeyValuePair<Dictionary<string, List<int>>, List<string>>(localDic, bStructures);
        }

        protected void PrepareClustering(List<string> _structNames)
        {
            structNames = new List<string>(_structNames.Count);
            // ProfileStat.CalcHammingDist(stateAlign, "C:\\Projects\\tmp\\tmp\\dist_rna40_IN");
            //ProfileStat.SelectFeatures(stateAlign, "C:\\Projects\\tmp\\tmp\\dist_rna40_half_IN.dat", "C:\\Projects\\tmp\\tmp\\dist_rna40_half_OUT.dat");
            //ProfileStat.CalcHammingDist(stateAlign, "C:\\Projects\\tmp\\tmp\\dist_rna40_OUT", "rna40_GaussFragBagSelected_OUT");
            //ProfileStat.CalcHammingDist(stateAlign, "C:\\Projects\\S35_IN", "dist_S35_IN_GaussFragBagNew");
            //ProfileStat.SelectFeatures(stateAlign, "C:\\Projects\\tmp\\tmp\\dist_CATH_clusters_IN", "C:\\Projects\\tmp\\tmp\\dist_CATH_clusters_OUT");
            foreach (var item in _structNames)
            {
                if (stateAlign.ContainsKey(item))
                    structNames.Add(item);
            }
            if (structNames.Count == 0)
                throw new Exception("No structures to cluster!");
            dicC = PrepareKeys(structNames,true,consensusProjection);
            structToKey = new Dictionary<string, string>();
            foreach(var item in dicC)
            {
                foreach(var d in item.Value)
                    structToKey.Add(structNames[d], item.Key);
            }
            

        }
        public double CalcClusterConsistency(List<string> cluster)
        {
            double res=0;
            if (cluster.Count == 1)
                return 1.0;

           
            List<byte> refStr = al.r.protCombineStates[cluster[0]];

            for (int j = 0; j < refStr.Count; j++)
                for (int i = 1; i < cluster.Count; i++)
                {
                    if (refStr[j] == al.r.protCombineStates[cluster[i]][j])
                        res++;
                }
            res/=refStr.Count*(cluster.Count-1);
            return res;

        }
        protected List<double> CalcClustersConsistency(List<List<string>> clusters)
        {
            List<double> consistency=new List<double>();
            foreach(var item in clusters)
            {

                        consistency.Add(CalcClusterConsistency(item));                    

            }

            return consistency;
        }
        public static Dictionary<string, KeyValuePair<List<string>, Dictionary<byte, int>[]>> CuttFreq(List<string> allNames, List<string> cuttedNames,Dictionary<string, KeyValuePair<List<string>, Dictionary<byte, int>[]>> dicToCutt)
        {
            Dictionary<string, KeyValuePair<List<string>, Dictionary<byte, int>[]>> newDic = new Dictionary<string, KeyValuePair<List<string>, Dictionary<byte, int>[]>>();
            bool[] index = new bool[allNames.Count];

            for (int i = 0; i < allNames.Count; i++)
                if (cuttedNames.Contains(allNames[i]))
                    index[i] = true;
                else
                    index[i] = false;


            foreach (var item in dicToCutt)
            {
                Dictionary<byte, int>[] freq = new Dictionary<byte, int>[cuttedNames.Count];
                for (int i = 0, j = 0; i < item.Value.Value.Length; i++)
                    if (index[i])
                        freq[j++] = item.Value.Value[i];

                newDic.Add(item.Key, new KeyValuePair<List<string>, Dictionary<byte, int>[]>(item.Value.Key, freq));
            }
            return newDic;

        }
        public void CuttAlignment(List<string> columnNames, List<string> clusterNames)
        {
            bool[] index = new bool[columnNames.Count];
            
            Dictionary<string, int> dicStruct = new Dictionary<string, int>();
            foreach (var item in clusterNames)
                dicStruct.Add(item, 0);


            for (int i = 0; i < columnNames.Count; i++)
                if (dicStruct.ContainsKey(columnNames[i]))
                    index[i] = true;
                else
                    index[i] = false;



            Dictionary<string, List<byte>> d = al.r.protCombineStates;
            
            List<string> aux = new List<string>(d.Keys);
            foreach (var item in aux)
            {                
                List<byte> states = new List<byte>();
                for (int i = 0; i < d[item].Count; i++)
                    if (index[i])
                        states.Add(d[item][i]);

                d.Remove(item);
                d.Add(item, states);
            }                    

        }

        public virtual ClusterOutput Cluster(List<string> _structNames,bool dendrog=false)
        {
            
            structNames = new List<string>(_structNames.Count);
            ClusterOutput output;

            PrepareClustering(_structNames);
            //output = PrepareClustersJuryLike(dicC, structNames);
            
            output = PrepareClusters(dicC, structNames);
            output.clusters.consistency = CalcClustersConsistency(output.clusters.list);
            if (validIndexes.Count > 0)
                output.auxInt = validIndexes;
            else
                output.auxInt = null;

            currentV = maxV;
            return output;
        }
        protected Dictionary <byte,int>[] JuryConsensusStates(Alignment alUse,List<string> structNames)
        {
            jury1D juryLocal = new jury1D();
            juryLocal.PrepareJury(alUse,stateAlign);

            Dictionary <byte ,int>[] cons=new Dictionary<byte,int>[columns.Length];


            ClusterOutput output = juryLocal.JuryOptWeights(structNames);
            if (output == null)
                return null;
            int i = 0;
            string result;
            if (structNames.Contains("consensus"))
                result = "consensus";
            else
                result = output.juryLike[0].Key;

            foreach (var item in stateAlign[result])
            {
                cons[i] = new Dictionary<byte, int>();
                cons[i++].Add(item, 0);
            }
            return cons;
        }
        private Dictionary<byte, int>[] HammingConsensusStates(Dictionary<byte, int>[] col)
        {
            Dictionary<byte, int>[] cons = new Dictionary<byte, int> [col.Length];

            for (int i = 0; i < col.Length; i++)
            {
                var orderedStates = col[i].OrderByDescending(j => j.Value);
                cons[i] = new Dictionary<byte, int>(col[i].Keys.Count);
                int n = 0;
                foreach (var item in orderedStates)                
                    cons[i].Add(item.Key, n++);                                   
                
            }

            return cons;
        }
        private void ThreadingMakeColumns(object o)
        {
            object[] array = o as object[];
            int threadNum = (int)array[0];
            int start = (int)array[1];
            int stop = (int)array[2];
            byte locState = 0;
            for (int i = start; i < stop; i++)
                for (int j = 0; j < allStructures.Count; j++)                
                {
                    if (auxDic.ContainsKey(allStructures[j]) && auxDic[allStructures[j]].Count > 0 && auxDic[allStructures[j]].Count>i)
                        locState = auxDic[allStructures[j]][i];
                    else
                        continue;
//                    if (locState == 0)
//                        continue;
                    
                    if (!col[i].ContainsKey(locState))
                        col[i].Add(locState, 1);
                    else
                        col[i][locState]++;

                }
                resetEvents[threadNum].Set();

        }
        protected Dictionary <byte,int>[] MakeColumnsLists(List<string> structNames,Dictionary<string, List<byte>> alUse)
        {
            if (structNames.Count == 0)
                return null;            

            col = new Dictionary<byte, int>[alUse[structNames[0]].Count];
            auxDic = alUse;
            allStructures = structNames;
            for (int i = 0; i < col.Length; i++)            
                col[i] = new Dictionary<byte, int>();

            resetEvents = new ManualResetEvent[threadNumbers];
            for (int n = 0; n < threadNumbers; n++)
            {
                int p = n;
                int start = n * col.Length/threadNumbers;
                int stop = (n + 1) * col.Length / threadNumbers;
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(ThreadingMakeColumns), new object[] { p, start,stop});
                

            }
            for (int i = 0; i < threadNumbers; i++)
                resetEvents[i].WaitOne();

            
            return col;
        }


    }
}
