using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using uQlustCore;
using uQlustCore.PDB;
using uQlustCore.dcd;
using System.Threading;

namespace uQlustCore.Distance
{
    public class Rmsd : DistanceMeasure
    {
        struct paramsThread
        {
            public int k;
           
        };
        class returnV
        {
            public float[,] x;
            public List<string>[] y;
        }
        Settings dirSettings = new Settings();
        //Dictionary<string, float[,]> aRotMol = new Dictionary<string, float[,]>();
        Dictionary<string, int> atomDic;

        protected DCDFile dcd = null;
        protected PDBMODE allAtoms;
        protected string refJuryProfile = null;
        protected bool flag;

        float[][,] pdbPos = null;
        protected List<string> fileNames = null;
        public Rmsd(DCDFile dcd, string alignFile, bool flag, PDBMODE allAtoms, string refJuryProfile = null)
        {
            order = true;
            dirSettings.Load();
            this.dcd = dcd;
            this.alignFile=alignFile;
            this.flag=flag;
            this.allAtoms=allAtoms;
            this.refJuryProfile = refJuryProfile;

        }
        public Rmsd(List <string> fileNames, string alignFile, bool flag, PDBMODE allAtoms,string refJuryProfile=null)
        {
            this.fileNames = fileNames;
            this.alignFile = alignFile;
            this.flag = flag;
            this.allAtoms = allAtoms;
            this.refJuryProfile = refJuryProfile;
        }
        public Rmsd(string dirName,string alignFile,bool flag,PDBMODE allAtoms,string refJuryProfile=null)
        {
            this.dirName = dirName;
            this.alignFile = alignFile;
            this.flag = flag;
            this.allAtoms = allAtoms;
            this.refJuryProfile = refJuryProfile;
        }
        public override void InitMeasure()
        {
            if (dcd != null)
                InitMeasure(dcd, alignFile, flag, allAtoms, refJuryProfile);
            else
                if (fileNames != null)
                    InitMeasure(fileNames, alignFile, flag, allAtoms, refJuryProfile);
                else
                    InitMeasure(dirName, alignFile, flag, allAtoms, refJuryProfile);
        }
        protected void InitMeasure(DCDFile dcd, string alignFile, bool flag, PDBMODE allAtoms, string refJuryProfile = null)
        {
            base.InitMeasure(dcd, alignFile, flag, refJuryProfile);
            DCDReader readDCD = new DCDReader(dcd);
              
             pdbs = new PDBFiles();
            readDCD.DCDPrepareReading(dcd.dcdFile, dcd.pdbFile);
            int counter = 0;            
            bool cont = true;
            do
            {
                MemoryStream mStream = new MemoryStream();
                StreamWriter wStream = new StreamWriter(mStream);
                cont=readDCD.ReadAndSavePDB(wStream);
                string name = "Model" + counter++;
                mStream.Position = 0;
                pdbs.AddPDB(mStream, allAtoms, name);
            }
            while(cont);

            readDCD.FinishDCDReading();

            pdbs.FindReferenceSeq();
            pdbs.MakeAlignment(alignFile);
            structNames = CheckAvailableStructures();

            
            atomDic=pdbs.MakeAllAtomsDic();
        }

        

        private Dictionary<string,int> CheckAvailableStructures()
        {

            Dictionary<string,int> strNames = new Dictionary<string,int>();

            if (jury == null)
            {
                foreach (var item in pdbs.molDic.Keys) 
                    strNames.Add(item,1);
                return strNames;
            }
            else
            {
                Dictionary<string, int> dic = new Dictionary<string, int>();

                foreach (var item in jury.alignKeys)
                {
                            dic.Add(item, 1);
                }
                foreach (var str in pdbs.molDic.Keys)
                    if (dic.ContainsKey(str))
                        dic[str]++;
                    else
                        dic.Add(str, 1);

                int max = -1;
                foreach (var item in dic)
                    if (max < item.Value)
                        max = item.Value;

                foreach (var item in dic)
                    if (item.Value == max)
                        strNames.Add(item.Key,1);
            }


            return strNames;
        }
        
        protected void InitMeasure(string dirName,string alignFile,bool flag,PDBMODE allAtoms,string refJuryProfile=null)
        {
            base.InitMeasure(dirName, alignFile, flag, refJuryProfile);            
            pdbs = new PDBFiles();
            string[] files;
            if(dirSettings.extension==null)
                files = Directory.GetFiles(dirName);
            else
                files = Directory.GetFiles(dirName, dirSettings.extension);

            if (files.Length == 0)
                throw new Exception("In selected directory " + dirName + " there are no files with extesion " + dirSettings.extension);

            maxV = files.Length;

            string refSeqFile = dirName + ".ref";
            pdbs.ReadRefSeq(refSeqFile);

            currentV = 0;
            pdbs.AddPDB(files, allAtoms,ref currentV);

            if (pdbs.molDic.Keys.Count == 0)
                throw new Exception("Non pdb files correctly read");
            pdbs.FindReferenceSeq();
            pdbs.MakeAlignment(alignFile);
            structNames = CheckAvailableStructures();
            order = true;
            atomDic = pdbs.MakeAllAtomsDic();
         }
        protected void InitMeasure(List <string> fileNames, string alignFile, bool flag, PDBMODE allAtoms,string refJuryProfile=null)
            
        {
            base.InitMeasure(fileNames, alignFile, flag, refJuryProfile);
            dirSettings.Load();
            pdbs = new PDBFiles();
            maxV = fileNames.Count;
            currentV = 0;
            pdbs.AddPDB(fileNames, allAtoms,ref currentV);
            if (pdbs.molDic.Keys.Count == 0)
                throw new Exception("Unable to read pdb data!!");
            pdbs.FindReferenceSeq();
            if(pdbs.molDic.Count==0)
            {
                throw new Exception("All pdb structures have been removed. No Data!");
            }
            pdbs.MakeAlignment(null);
            structNames = CheckAvailableStructures();
            order = true;
            atomDic = pdbs.MakeAllAtomsDic();
        }
               
        public override string ToString()
        {
            return "Rmsd";
        }
        public override double ProgressUpdate()
        {
            return base.ProgressUpdate();
        }

        public override bool SimilarityThreshold(float threshold,float dist)
        {
            if (dist < threshold)
                return true;

            return false;
        }
        public override List<KeyValuePair<string, double>> GetReferenceList(List<string> structures)
        {
            float[,,] refPosition = null;

            if (jury != null)
                //return jury.JuryOpt(structures).juryLike[0].Key;
                return jury.JuryOptWeights(structures).juryLike;
                //return jury.ConsensusJury(structures).juryLike;

            Dictionary<string, float[,,]> aRotMol = new Dictionary<string, float[,,]>(structures.Count);
            refPosition = CalcRefPosition(structures,aRotMol);
            List<KeyValuePair<string, double>> refList = new List<KeyValuePair<string, double>>();
            foreach (var item in aRotMol.Keys)
            {
                double sum = 0;
                int count = 0;
                for (int i = 0; i < pdbs.molDic[item].indexMol.Length; i++)
                    for (int j = 0; j < aRotMol[item].GetLength(1); j++)
                        for (int n = 0; n < aRotMol[item].GetLength(2); n++)
                            if (pdbs.molDic[item].indexMol[i] != -1)
                            {
                                sum += (refPosition[i,j, n] - aRotMol[item][i,j, n]) *
                                    (refPosition[i,j, n] - aRotMol[item][i,j, n]);
                                count++;
                            }

                if (count > 0)
                {
                    KeyValuePair<string, double> aux = new KeyValuePair<string, double>(item, Math.Sqrt(sum) / count);
                    refList.Add(aux);
                }
            }
            refList.Sort((nextPair,firstPair) =>
            {
                return nextPair.Value.CompareTo(firstPair.Value);
            });

            return refList;


        }
        private float[,,] CalcRefPosition(List<string> structures,Dictionary<string, float[,,]> aRotMol)
        {
            //float[,] refPosition = null;
            int[] counts;
            aRotMol.Clear();

            if (pdbs.molDic.Keys.Count == 1)
                pdbs.EmptyAlignment(structures[0]);
            
            float[, ,] refPosition = new float[pdbs.molDic[structures[0]].indexMol.Length, atomDic.Keys.Count, 3];
            //refPosition = new float[maxV, 3];
            counts = new int[pdbs.molDic[structures[0]].indexMol.Length];


            //Find longest structure (based on Atoms)
            string refLength="";
            int maxLength=0;
            foreach (var item in structures)
            {
                int localMax=0;
                foreach (var it in pdbs.molDic[item].mol.Chains[0].Residues)
                    localMax += it.Atoms.Count;
                if (localMax > maxLength)
                {
                    maxLength = localMax;
                    refLength = item;
                }
            }

            foreach (var item in structures)
            {
                    returnV aux=GetStructAfterRotation(refLength, item);
                    if(aux==null)
                        continue;
                    float[, ,] auxK = new float[pdbs.molDic[refLength].indexMol.Length, atomDic.Keys.Count, 3];
                    
                    
                    for (int i = 0,m=0; i < aux.y.Length; i++)
                    {
                        for (int j = 0; j < aux.y[i].Count; j++, m++)
                        {
                            for (int n = 0; n < 3; n++)
                            {
                                refPosition[i, atomDic[aux.y[i][j]], n] += aux.x[m, n];
                                auxK[i, atomDic[aux.y[i][j]], n] = aux.x[m, n];
                            }

                            counts[i]++;
                        }
                         

                    }
                    aRotMol.Add(item, auxK);
            }
            for (int i = 0; i < refPosition.GetLength(0); i++)
                if (counts[i] > 0)
                    for (int j = 0; j < refPosition.GetLength(1); j++)
                        for (int n= 0; n < refPosition.GetLength(2); n++)
                        refPosition[i, j,n] /= counts[i];


                    return refPosition;
        }

        public override string GetReferenceStructure(List<string> structures)
        {
            List<KeyValuePair<string, double>> refList=null;
                refList = GetReferenceList(structures);

            return refList[0].Key;
        }
        private returnV GetStructAfterRotation(string refStructure, string modelStructure)
        {

            float[,] transMatrix=null;
            float []cent1=new float [3];
            float []cent2=new float [3];

            if (!pdbs.molDic.ContainsKey(refStructure) || !pdbs.molDic.ContainsKey(modelStructure))
                return null;
            posMOL locPosMol=Optimization.PrepareData(pdbs.molDic[refStructure], pdbs.molDic[modelStructure],true);

            //opt.Rmsd(locPosMol.posmol1, locPosMol.posmol2, cent1, cent2, true);
           /* float [,]mm=new float [3,3];

            mm[0, 0] = 0; mm[0, 1] = 0; mm[0, 2] = 1;
            mm[1, 0] = 1; mm[1, 1] = 0; mm[1, 2] = 0;
            mm[2, 0] = 0; mm[2, 1] = 1; mm[2, 2] = 0;

            float [,]test=Optimization.MultMatrixTrans(locPosMol.posmol1, mm);*/
            
           transMatrix = Optimization.TransMatrix(locPosMol.posmol2, locPosMol.posmol1);
            //transMatrix = opt.TransMatrix(locPosMol.posmol1, test);
            if (transMatrix == null)
                return null;

//            float[,] akk = Optimization.MultMatrixTrans(test, transMatrix);
            transMatrix = Optimization.MultMatrixTrans(locPosMol.posmol2, transMatrix);
            returnV r=new returnV();
            r.x = transMatrix;
            r.y = locPosMol.atoms;
            return r;
            //return Optimization.MultMatrixTrans(test, transMatrix);
            //return akk;
            //return opt.posMolRot;
            //return opt.TransMatrix(opt.posMol1, opt.posMol2);

        }


        public override int GetDistance(string refStructure, string modelStructure)
        {
            if (!pdbs.molDic.ContainsKey(refStructure) || !pdbs.molDic.ContainsKey(modelStructure))
                return errorValue;
            if (pdbPos == null)
            {
                posMOL locPosMol = Optimization.PrepareData(pdbs.molDic[refStructure], pdbs.molDic[modelStructure]);
                return (int)(Optimization.Rmsd(locPosMol.posmol1, locPosMol.posmol2, false) * 100);
            }
            else
            {
                return (int)(Optimization.Rmsd(pdbPos[hashIndex[refStructure]], pdbPos[hashIndex[modelStructure]],false) * 100);
            }
        }

        private void CalcMatrixList(object o)
        {
            int num = (int)o;
          
            for (int i = 0; i < indexList[num].Count; i++)
                distanceMatrix[indexList[num][i].Key] = GetDistance(structures[structures.Count - 1], structures[indexList[num][i].Key]);

            resetEvents[num].Set();
        }
        protected override void CalcMatrix(object o)
        {
            //StreamWriter oo = new StreamWriter("dist.dat");
            paramsThread p = (paramsThread)o;

            int maxV = -1;
            foreach (var item in indexList[p.k])
            {
                int dd = FindIndex(item.Key);

                int distC = GetDistance(structures[item.Key], structures[item.Value]);
                if (distC < 0 || distC==errorValue)
                    distC = -1;

                distanceMatrix[dd + item.Value - item.Key] = distC;
                if (distC > maxV)
                    maxV = distC;
              //  oo.WriteLine(p.structures[item.Key] + " " + p.structures[item.Value] + " " + distC);
                Interlocked.Increment(ref currentV);
            }
            resetEvents[p.k].Set();
          //  oo.Close();
        }
        public override int[] GetDistance(string refStructure, List<string> structures)
        {
            int[] dist = new int[structures.Count];
            distanceMatrix = new int[structures.Count];
            Settings set = new Settings();
            set.Load();
            int threadNumbers = set.numberOfCores;
            this.structures = new List<string>(structures);
            int part = this.structNames.Count / threadNumbers;
            indexList = new List<KeyValuePair<int, int>>[threadNumbers];
            for (int j = 0; j < indexList.Length; j++)
            {
                indexList[j] = new List<KeyValuePair<int, int>>();
            }
            int count = 0;
            for (int j = 0; j < this.structures.Count; j++)
            {
                if (indexList[count].Count < part)
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
        private void PreparePDBPos()
        {
            //Check if all alignm have the same length
            foreach(var item in pdbs.molDic)
            {
                foreach (var itemIndex in item.Value.indexMol)
                    if (itemIndex == -1)
                        return;
            }
            pdbPos = new float[structures.Count][,];
            for (int i = 0; i < structures.Count;i++)
                pdbPos[i] = new float[pdbs.molDic[structures[0]].indexMol.Length, 3];
            float []center=new float [3]; 
            for (int i = 0; i < structures.Count;i++ )
            {
                int count=0;

                for (int j = 0; j < pdbs.molDic[structures[i]].mol.Residues.Count; j++)
                    for (int n = 0; n < pdbs.molDic[structures[i]].mol.Residues[j].Atoms.Count; n++)
                    {
                        pdbPos[i][count, 0] = pdbs.molDic[structures[i]].mol.Residues[j].Atoms[n].Position.X;
                        pdbPos[i][count, 1] = pdbs.molDic[structures[i]].mol.Residues[j].Atoms[n].Position.Y;
                        pdbPos[i][count++,2] = pdbs.molDic[structures[i]].mol.Residues[j].Atoms[n].Position.Z;
                    }

                Optimization.CenterMol(pdbPos[i],center);
            }

        }
        public override void CalcDistMatrix(List <string> structures)
        {
            currentV = 0;
            hashIndex = new Dictionary<string, int>(structures.Count);
            for(int i=0;i<structures.Count;i++)
            {
                if (!pdbs.molDic.ContainsKey(structures[i]))
                    throw new Exception("There is no pdb file or file is incorrect for the structure: " + structures[i]);
                hashIndex.Add(structures[i], i);
            }
            Settings set = new Settings();
            set.Load();
            int threadNumbers =  set.numberOfCores;
            this.structures = new List<string>(structures);

            int part = structures.Count * (structures.Count+1) / (2 * threadNumbers)+1;
         
            distanceMatrix = new int[structures.Count * (structures.Count + 1) / 2];
            indexList = new List<KeyValuePair<int, int>>[threadNumbers];
            for (int n= 0; n < threadNumbers;n++)           
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
                for (int n = 0; n < threadNumbers; n++)
                    maxV += indexList[n].Count;

              //  PreparePDBPos();
            resetEvents = new ManualResetEvent[threadNumbers];
            for(int n=0;n<threadNumbers;n++)
            {
                paramsThread p=new paramsThread();
                p.k=n;
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(CalcMatrix),(object) p);
                

            }
            for (int j = 0; j < threadNumbers; j++)
                resetEvents[j].WaitOne();

            currentV = maxV;
        }
    
    }
}
