using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.IO;
using uQlustCore;
using uQlustCore.PDB;
using System.Globalization;
using System.Threading;
using uQlustCore.Interface;

namespace uQlustCore.Distance
{
    
    public abstract class DistanceMeasure:IProgressBar
    {
        protected int[] distanceMatrix = null;
        public Dictionary<string, int> hashIndex;// = new Dictionary<string, int>();
        protected ManualResetEvent[] resetEvents = null;
        protected Dictionary<string, Dictionary<string, float [,]>> rotStruct = new Dictionary<string, Dictionary<string, float [,]>>();
        public Dictionary<string,int> structNames;
        protected List<string> structures;
//		protected  Optimization opt = null;
        public PDBFiles pdbs = null;
        public bool order; //false if best are biggest values otherwise like in regular distance true.
        public double maxSimilarity;
        private const int ErrorValue = int.MaxValue;
        public int errorValue { get { return ErrorValue; } }
        public string dirName;
        public string alignFile;
        protected List<KeyValuePair<int, int>>[] indexList;
        public jury1D jury = null;
        protected string juryProfile = null;

        protected int maxV=1;
        protected int currentV;
        public DistanceMeasure()
        {

        }    
      
        public virtual void InitMeasure(DCDFile dcd, string alignFile, bool referenceFlag, string juryProfile = null)
        {
            this.alignFile = alignFile;
            if (referenceFlag)
            {
                if (juryProfile.Length == 0 || !File.Exists(juryProfile))
                    throw new Exception("Profile for reference structure searching has been not defined or does not exist");
                jury = new jury1D();
                jury.PrepareJury(dcd, alignFile, juryProfile);
            }
        }
        public virtual void InitMeasure(string dirName,string alignFile,bool referenceFlag,string juryProfile=null)
        {
            this.alignFile = alignFile;
            this.dirName = dirName;
            if (referenceFlag)
            {
                jury = new jury1D();
                jury.PrepareJury(dirName, alignFile, juryProfile);
                //AddErrors(jury.errors);
            }
        }
        public virtual void InitMeasure(List<string> fileNames, string alignFile, bool referenceFlag, string juryProfile = null)
        {
            this.alignFile = alignFile;
            if (referenceFlag)
            {
                jury = new jury1D();
                 jury.PrepareJury(fileNames, alignFile, juryProfile);
                //AddErrors(jury.errors);
            }
        }
        public virtual void InitMeasure(string profilesFile, bool referenceFlag, string juryProfile = null, string refJuryProfile = null)
        {
            this.alignFile = profilesFile;

            if (referenceFlag)
            {
                jury = new jury1D();
                jury.PrepareJury(profilesFile, refJuryProfile);
                //AddErrors(jury.errors);
            }
        }
        public virtual void InitMeasure(Alignment al, bool flag)
        {
        }
        public abstract void InitMeasure();
        public virtual int GetDistance(string m, string n){ return 0;}
        public virtual double GetRealValue(double v){return v/100; }
        public virtual bool SimilarityThreshold(float threshold, float dist){return true;}


        public virtual double ProgressUpdate()
        {
            return (double)currentV / maxV;
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }

        protected int FindIndex(int matIndex)
        {
            return (2*hashIndex.Count - matIndex + 1) * matIndex / 2;
        }
        public virtual int GetDistance(float[,] refStructure, string modelStructure, out int[] index) {index = null; return 0; }
        //public virtual boo
        public virtual string GetReferenceStructure(List<string> structures)
        {
            ClusterOutput clustOut;

            clustOut = jury.JuryOptWeights(structures);
            if (clustOut == null)
                return null;
            //clustOut = jury.JuryOpt(structures);

            return clustOut.juryLike[0].Key;

        }

        public string GetReferenceStructure(List<string> structures, List<string> refStruct)
        {
            List<KeyValuePair<string, double>> refList = null;
            refList = GetReferenceList(structures);

            if (refStruct != null)
            {
                List<KeyValuePair<string, int>> listK = new List<KeyValuePair<string, int>>();
                foreach (var item in refList)
                {

                   // int remIndex = 0;
                    for (int j = 0; j < refStruct.Count; j++)
                        if (refStruct[j] == item.Key)
                        {
                            return item.Key;
                            //remIndex = j;
                            //break;
                        }
                    /*KeyValuePair<string,int> aux = new KeyValuePair<string, int>(item, remIndex);
                    listK.Add(aux);
                    listK.Sort((firstPair, nextPair) =>
                    {
                        return firstPair.Value.CompareTo(nextPair.Value);
                    });
                    return listK[0].Key;*/
                }

            }
            //return refList[0].Key;
            return refList[0].Key;
        }

        public virtual List<string> GetReferenceStructureList(List<string> structures)
        {
            List<string> refStruct = new List<string>();
            ClusterOutput op;
            if (jury != null)
            {
                //op = jury.JuryOpt(structures);
                op = jury.JuryOptWeights(structures);
                if (op == null)
                    return null;
                if (op.juryLike[0].Value - op.juryLike[op.juryLike.Count - 1].Value > 0.1)
                {
                    refStruct.Add(op.juryLike[0].Key);
                    refStruct.Add(op.juryLike[op.juryLike.Count - 1].Key);
                    return refStruct;
                }
            }
            return null;
        }
        protected virtual void CalcMatrix(object o)
        {

            int num = (int)o;

            for (int i = 0; i < distanceMatrix.Length; i++)
                distanceMatrix[i] = 100;

            string line="";
            StreamReader distRead; /*= new StreamReader("C:\\Projects\\dist.dat");

            line = distRead.ReadLine();
            while (line != null)
            {
                string[] aux = line.Split(' ');

                if (aux[0].Length > 4)
                    aux[0] = "pdb" + aux[0].Substring(0, 4) + ".ent_" + Char.ToUpper(aux[0][4]) + ".pdb";
                else
                    aux[0] = "pdb" + aux[0].Substring(0, 4) + ".ent.pdb";
                if (aux[1].Length > 4)
                    aux[1] = "pdb" + aux[1].Substring(0, 4) + ".ent_" + Char.ToUpper(aux[1][4]) + ".pdb";
                else
                    aux[1] = "pdb" + aux[1].Substring(0, 4) + ".ent.pdb";
                if (hashIndex.ContainsKey(aux[0]) && hashIndex.ContainsKey(aux[1]))
                {
                    int x, y, tmp;
                    x = hashIndex[aux[0]];
                    y = hashIndex[aux[1]];
                    if (x > y)
                    {
                        tmp = y;
                        y = x;
                        x = tmp;
                    }
                    int dd = FindIndex(x);
                    if (aux.Length == 3 && aux[2].Length > 0)
                        distanceMatrix[dd + y - x] = 100 - (int)(float.Parse(aux[2], CultureInfo.InvariantCulture) /10);
                }
                line = distRead.ReadLine();
            }

            distRead.Close();*/


            //distRead = new StreamReader("C:\\Projects\\final");
            distRead = new StreamReader("C:\\Projects\\Dist\\tm_score.dat");

            line = distRead.ReadLine();
            while(line!=null)
            {
                string[] aux = line.Split(' ');

                /*if(aux[0].Length>4)
                    aux[0] = "pdb" + aux[0].Substring(0, 4) + ".ent_"+ Char.ToUpper(aux[0][4]) + ".pdb";
                else
                    aux[0]= "pdb" + aux[0].Substring(0, 4) + ".ent.pdb";
                if(aux[1].Length>4)
                    aux[1] = "pdb" + aux[1].Substring(0, 4) + ".ent_" + Char.ToUpper(aux[1][4]) + ".pdb";
                else
                    aux[1] = "pdb" + aux[1].Substring(0, 4) + ".ent.pdb";*/
                if (hashIndex.ContainsKey(aux[0]) && hashIndex.ContainsKey(aux[1]))
                {
                    int x, y,tmp;
                    x = hashIndex[aux[0]];
                    y = hashIndex[aux[1]];
                    if (x > y)
                    {
                        tmp = y;
                        y = x;
                        x = tmp;
                    }
                    int dd = FindIndex(x);
                    if (aux.Length > 2 && aux[2].Length > 0 && distanceMatrix[dd + y - x] == 100)
                    {
                        //int distC = GetDistance(aux[0], aux[1]);
                        int distOrg = 100-(int)(float.Parse(aux[2], CultureInfo.InvariantCulture) * 100);
                        //if (distOrg > 50)
                            //distOrg = distC;
                        distanceMatrix[dd + y - x] = distOrg;
                    }
                }
                line = distRead.ReadLine();
            }

            distRead.Close();

            StreamWriter bb = new StreamWriter("vvv.dat");
            for(int i=0;i<structures.Count;i++)
                for(int j=i+1;j<structures.Count;j++)
                {
                    int tmp;
                    int x = hashIndex[structures[i]];
                    int y = hashIndex[structures[j]];
                    if (x > y)
                    {
                        tmp = y;
                        y = x;
                        x = tmp;
                    }
                    int dd = FindIndex(x);

                    if (distanceMatrix[dd + y - x] == 100)
                        bb.WriteLine(structures[i] + " " + structures[j]);
                }

            bb.Close();
            resetEvents[num].Set();

        }
/*        protected virtual void CalcMatrix(object o)
        {
            int num = (int)o;
            foreach (var item in indexList[num])
            {
                int dd = FindIndex(item.Key);

                int distC = GetDistance(structures[item.Key], structures[item.Value]);
                if (distC < 0 || distC >= int.MaxValue)
                    distC = -1;

                distanceMatrix[dd + item.Value - item.Key] = distC;
                Interlocked.Increment(ref currentV);
            }
            resetEvents[num].Set();
        }*/

        public virtual void CalcDistMatrix(List<string> structures) 
        {
            currentV = 0;

            hashIndex = new Dictionary<string, int>(structures.Count);
            for (int i = 0; i < structures.Count; i++)               
                hashIndex.Add(structures[i], i);
            
            Settings set = new Settings();
            set.Load();
            int threadNumbers = set.numberOfCores;
            this.structures = new List<string>(structures);

            int part = structures.Count * (structures.Count + 1) / (2 * threadNumbers) + 1;          

            distanceMatrix = new int[structures.Count * (structures.Count + 1) / 2];
            indexList = new List<KeyValuePair<int, int>>[threadNumbers];
          
            for (int n = 0; n < threadNumbers; n++)            
                indexList[n] = new List<KeyValuePair<int, int>>(part);

            int count = 0;
            foreach (var dim1 in structures)
            {
                foreach (var dim2 in structures)
                {
                    if (hashIndex[dim1] <= hashIndex[dim2])
                    {
                        indexList[count].Add(new KeyValuePair<int, int>(hashIndex[dim1], hashIndex[dim2]));
                        if (part <= indexList[count].Count)
                            count++;
                    }
                }
            }
            threadNumbers=1;
            for (int n = 0; n < threadNumbers; n++)
                maxV += indexList[n].Count;

            resetEvents = new ManualResetEvent[threadNumbers];
            for (int n = 0; n < threadNumbers; n++)
            {
                int p;
                p = n;
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(CalcMatrix), (object)p);
                
            }
            for (int j = 0; j < threadNumbers; j++)
                resetEvents[j].WaitOne();
            DebugClass.WriteMessage("Distance finished");
            currentV = maxV;
        }
        public virtual List<KeyValuePair<string, double>> GetReferenceList(List<string> structures)
        {

            ClusterOutput clustOut;
            clustOut = jury.JuryOptWeights(structures);
            if (clustOut == null)
                return null;
            return clustOut.juryLike;

        }

        public virtual float[,] GetReferenceStructure(string refStruct,List <string> structures)
        {
            float[,] refPosition = new float[rotStruct[refStruct][structures[0]].GetLength(0), 3];
            int[] count = new int[rotStruct[refStruct][structures[0]].GetLength(0)];
            for (int i = 0; i < structures.Count;i++ )
            {
                float[,] aux = rotStruct[refStruct][structures[i]];
                for (int n = 0; n <aux.GetLength(0); n++)
                {
                    if (aux[n, 0] != float.MaxValue)
                    {                        
                        for (int j = 0; j < 3; j++)                        
                            refPosition[n, j] += aux[n, j];
            
                        count[n]++;
                    }
                }
            }
            for (int i = 0; i < refPosition.GetLength(0); i++)
                for (int j = 0; j < refPosition.GetLength(1); j++)
                    if (count[i] > 0 && refPosition[i,j]!=float.MaxValue)
                        refPosition[i, j] /= count[i];

            return refPosition;
        }
      public virtual int[][] GetDistance(List<string> refStructure, List<string> structures)
      {
          int[][] dist = new int[refStructure.Count][];
          rotStruct.Clear();
          for (int j = 0; j < refStructure.Count; j++)
          {
              dist[j] = new int[structures.Count];

              //foreach (string item in structures)
                  dist[j]= GetDistance(refStructure[j], structures);

          }

          return dist;

      }
       /* public virtual int[][] GetDistance(List<string> refStructure, List<string> structures)
        {
            int[][] dist = new int[refStructure.Count][];
            int i = 0;
			if(opt!=null)
				rotStruct.Clear();
            for (int j = 0; j < refStructure.Count; j++)
            {
                dist[j] = new int[structures.Count];
                i = 0;
				
                foreach (string item in structures)				
                    dist[j][i++] = GetDistance(refStructure[j], item);
									
            }

            return dist;

        }*/
        //private void int[] 
        /*public virtual int[] GetDistance(string refStructure, List<string> structures)
        {
            int[] dist = new int[structures.Count];
            int i = 0;
            foreach (string item in this.structures)
                dist[i++] = GetDistance(refStructure, item);

            return dist;

        }*/

        private void CalcMatrixList(object o)
        {
            int num = (int)o;

            for (int i = 0; i < indexList[num].Count; i++)
                distanceMatrix[indexList[num][i].Key] = GetDistance(structures[structures.Count - 1], structures[indexList[num][i].Key]);

            resetEvents[num].Set();
        }
        public virtual int[] GetDistance(string refStructure, List<string> structures)
        {
            int[] dist = new int[structures.Count];
            distanceMatrix = new int[structures.Count];
            Settings set = new Settings();
            set.Load();
            int threadNumbers = set.numberOfCores;
            this.structures = new List<string>(structures);
            double part = (double)this.structNames.Count / threadNumbers;
            indexList = new List<KeyValuePair<int, int>>[threadNumbers];

            for (int j = 0; j < indexList.Length; j++)
                indexList[j] = new List<KeyValuePair<int, int>>();

            int count = 0;
            for (int j = 0; j < structures.Count; j++)
            {
                if (indexList[count].Count < part*(j+1)-part*j)
                {
                    KeyValuePair<int, int> aux = new KeyValuePair<int, int>(j, 0);
                    indexList[count].Add(aux);
                }
                else
                    count++;
            }
            this.structures.Add(refStructure);
            resetEvents = new ManualResetEvent[threadNumbers];
            for (int n = 0; n < threadNumbers; n++)
            {
                int k = n;
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(CalcMatrixList), (object)k);
                
            }
            for (int j = 0; j < threadNumbers; j++)
                resetEvents[j].WaitOne();


            for (int j = 0; j < this.structures.Count - 1; j++)
                dist[j] = distanceMatrix[j];

            return dist;

        }
		public virtual int[] GetDistance(float [,] refStructure, List<string> structures)
        {
            int[] dist = new int[structures.Count];
            int i = 0;
            foreach (string item in structures)
            {
                int[] index;
                dist[i++] = GetDistance(refStructure, item,out index);
            }

            return dist;

        }
        public int GetIndex(string i,string j)
        {
            return GetIndex(hashIndex[i], hashIndex[j]);
        }
		public int GetIndex(int i,int j)
        {
            int m;

            if (i < j)
            {
                m = FindIndex(i);
                //Console.WriteLine("index=" + m+" "+(j-i));
                return FindIndex(i) + j - i;
            }
            else
            {
                m = FindIndex(j);
                //Console.WriteLine("index=" + m + " " + (j - i));
                return FindIndex(j) + i - j;
            }

        }
		
        public int GetDistance(int i, int j)
        {
            if (distanceMatrix == null)
                return -1;

            return distanceMatrix[GetIndex(i, j)];
        }

        public int GetMaxDistance()
        {
            if (distanceMatrix == null)
                return -1;

            int max = distanceMatrix[0];

            for (int i = 1; i < distanceMatrix.GetLength(0); i++)
                    if (max < distanceMatrix[i])
                        max = distanceMatrix[i];

            return max;
        }
        public KeyValuePair<int,int> FindMinimalDistance(HClusterNode clust1, HClusterNode clust2, AglomerativeType linkageType)
        {
            int dist=0;
            int index=0;
            int min = 0;
            switch (linkageType)
            {
                case AglomerativeType.SINGLE:
                    min = Int32.MaxValue;
                    foreach (var item1 in clust1.setStruct)
                    {
                        int v = hashIndex[item1];
                        foreach (var item2 in clust2.setStruct)
                        {
                            dist = GetDistance(v, hashIndex[item2]);
                            if (dist < min)
                            {
                                min = dist;
                                index = GetIndex(v, hashIndex[item2]);
                            }
                        }
                    }
                    break;
                case AglomerativeType.AVERAGE:
                    dist = GetDistance(hashIndex[clust1.refStructure], hashIndex[clust2.refStructure]);

                    index = GetIndex(hashIndex[clust1.refStructure], hashIndex[clust2.refStructure]);
                    break;
                case AglomerativeType.COMPLETE:
                    min = Int32.MinValue;
                    foreach (var item1 in clust1.setStruct)
                    {
                        int v = hashIndex[item1];
                        foreach (var item2 in clust2.setStruct)
                        {
                            dist = GetDistance(v, hashIndex[item2]);
                            if (dist > min)
                            {
                                min = dist;
                                index = GetIndex(v, hashIndex[item2]);
                            }
                        }
                    }
                    break;

            }
            return new KeyValuePair<int, int>(dist, index) ;
        }

    }

}
