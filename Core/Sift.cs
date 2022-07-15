using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using uQlustCore.PDB;
using System.IO;
using System.Data;
using uQlustCore.Interface;

namespace uQlustCore
{
    public class Sift:IProgressBar
    {
        PDBFiles pdbs = null;
        Settings dirSettings = new Settings();
        public Dictionary<string, List<double>> bins = new Dictionary<string, List<double>>();
        Random rand = new Random();
        public List<double> binLen = new List<double>();
        List<string> listFiles = null;
        List<KeyValuePair<string, double>> field = new List<KeyValuePair<string, double>>();

        int maxV = 1;
        int currentV = 0;

        int binNumber=40;
        int start=2;
        int stop=11;

        double maxStart=3.5;
        double maxStop=6.0;
        double minStart=6.5;
        double minStop=8.5;
        double step;

        public Sift()
        {

        }
        public ClusterOutput RunSift(List<string> fileList)
        {
            this.listFiles = fileList;

            return Shape();
        }

        public  ClusterOutput RunSift(string dirName)
        {
            Options opt=new Options();
            opt.ReadDefaultFile();

            dirSettings.Load();
            pdbs = new PDBFiles();
            string[] files;

            if(dirSettings.extension==null)
                files = Directory.GetFiles(dirName);
            else
                files = Directory.GetFiles(dirName, dirSettings.extension);

            if (files.Length == 0)
                throw new Exception("In selected directory " + dirName + " there are no files with extesion " + dirSettings.extension);

            listFiles = new List<string>(files.Length);
            foreach(var item in files)
                listFiles.Add(item);

            return Shape();

        }
        public double ProgressUpdate()
        {
            return (double)currentV/maxV;
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }
        void PrepareSift()
        {
            pdbs = new PDBFiles();
            maxV = listFiles.Count * 2+1;
            foreach (var item in listFiles)
            {
                string[] aux = item.Split('\\');
                //  if (!fileList.ContainsKey(aux[aux.Length-1]))
                //       continue;

                string strName;
                strName = pdbs.AddPDB(item, PDBMODE.ALL_ATOMS);
                if (strName != null)
                {
                    pdbs.molDic[strName].CreateSideChainContactMap(11.0f, true);
                    pdbs.molDic[strName].CleanMolData();

                }
                currentV++;
            }
            step = (double)(stop - start) / binNumber;

            for (double i = start; i < stop; i += step)
                binLen.Add(i + step);

        }
        ClusterOutput Shape()
        {
            PrepareSift();
            MakeHistogram();
            FindMaxMin();
            field.Sort((firstPair, nextPair) =>
            {
                return nextPair.Value.CompareTo(firstPair.Value);
            });

            ClusterOutput siftOut = new ClusterOutput();

            siftOut.juryLike = field;
            currentV = maxV;
            return siftOut;

        }
        double GetNormalRandom(double mean, double stdDev)
        {
            
            double u1 = rand.NextDouble();
            double u2 = rand.NextDouble();

            double r = Math.Sqrt(-2.0 * Math.Log(u1));
            double theta = 2.0 * Math.PI * u2;
            double st= r * Math.Sin(theta);

            return mean + stdDev * st;
        }
        void MakeHistogram()
        {
            
            List <double> binCount=new List<double>();
            foreach (var item in pdbs.molDic.Keys)
            {
                List<double> lContacts=new List<double>();

                foreach( var it in pdbs.molDic[item].contactMap.Keys)
                {
                    for (int j = 0; j < pdbs.molDic[item].contactMap[it].Count;j++)
                    {
                        lContacts.Add(pdbs.molDic[item].contactMap[it][j]);
                        for (int n = 0; n < 50; n++)
                            lContacts.Add(pdbs.molDic[item].contactMap[it][j] + GetNormalRandom(0, 0.25));
                    }

                }
                lContacts.Sort();//Sprawdzic czy sortuje od najmniejszego do najwiekszego
                binCount.Clear();
                int k = 0;
                for (double i = start; i < stop; i += step)
                {
                 
                    int counter = 0;
                    while (k < lContacts.Count && lContacts[k] < i + step)
                    {
                        counter++;
                        k++;
                    }
                    binCount.Add(counter);
                }

                for (int i = 0; i < binCount.Count; i++)
                {
                    double tmp = (binLen[i] - step);
                    binCount[i] /= binLen[i] * binLen[i] * binLen[i] - tmp*tmp*tmp;
                    binCount[i] /= lContacts.Count;
                    binCount[i] *= 1331;// 11 * 11 * 11;
                }
                //Average histogram
                List<double> binCopy = new List<double>(binCount);
                for (int i = 1; i < binCount.Count - 1; i++)
                {
                    binCopy[i] = (binCount[i - 1] + binCount[i] + binCount[i + 1]) / 3;
                }                                
                bins.Add(item, binCopy);
                currentV++;
            }
        }
        void FindMaxMin()
        {

            foreach(var item in bins.Keys)
            {
                double maxV = 0;
                double minV = 0;
                int maxCount=0;
                int minCount=0;
                for (int i =0;i<binLen.Count;i++)
                {

                    if (binLen[i] >= start && binLen[i] <= stop)
                    {
                        if (binLen[i] >= maxStart && binLen[i] <= maxStop)
                        {
                            maxV+=bins[item][i];
                            maxCount++;
                        }
                        if (binLen[i] >= minStart && binLen[i] <= minStop)
                        {
                            minV += bins[item][i];
                            minCount++;
                        }

                    }   
                }
                if(maxCount>0)
                    maxV /= maxCount;
                if(minCount>0)
                    minV /= minCount;

                double tres;
                if (minV > maxV)
                    tres = maxV;
                else
                    tres = minV;

                Field(tres + Math.Abs(maxV - minV) / 2, item);

            }

        }

        void Field(double level,string item)
        {
            double sumUp, sumDown;

            sumUp = sumDown = 0;
            for (int i = 0; i < binLen.Count; i++)
            {
                if (binLen[i] >= maxStart && binLen[i] <= maxStop)
                {
                    double tmp = bins[item][i] - level;
                    if(tmp>0)
                        sumUp += tmp;
                }
                if (binLen[i] >= minStart && binLen[i] <= minStop)
                {
                    double tmp = level-bins[item][i];
                    if (tmp > 0)
                        sumDown += tmp;

                }

            }
            double field1 = level * (maxStop - maxStart);
            double field2 = level * (minStop - maxStart) - sumDown;

            sumUp += sumDown;
            KeyValuePair<string, double> v = new KeyValuePair<string, double>(item, sumUp * sumUp);
            field.Add(v);

        }
    }
}
