using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using uQlustCore;
using System.IO;
using uQlustCore.Interface;
using System.Data;
using uQlustCore.Distance;

namespace uQlustCore
{
   public struct fractionParams
   {
       public DistanceMeasures distance;
       public string profileName;
       public double distThreshold;
       public int clustersNum;
   };
   public class Fraction:IProgressBar
    {
        int maxV=1;
        int currentV=0;
        int errorValue = int.MaxValue;
        DataTable results;
        Exception exc = null;
        struct BestModel
        {
            public string modelName;
            public int value;
        };
      
        List<ClusterOutput> output;
        public Fraction(List<ClusterOutput> output)
        {
            if (output == null)
                throw new Exception("No data, nothing to be done!");
            this.output = output;
            
        }
        public double ProgressUpdate()
        {
            return (double)currentV / maxV;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            List<KeyValuePair<string, DataTable>> t = new List<KeyValuePair<string, DataTable>>();
            t.Add(new KeyValuePair<string,DataTable>("Fraction of good",results));

            return t;
        }
        public Exception GetException()
        {
            return exc;
        }
        public void GetFraction(object param)
        {

            float fraction = 0;
            float avrDiff=0;
            int counter = 0;
            int end;
            DistanceMeasures distance=((fractionParams)param).distance;
            string profileName = ((fractionParams)param).profileName;
            double distThreshold = ((fractionParams)param).distThreshold;
            int clustersNum = ((fractionParams)param).clustersNum;

            results = new DataTable();
            try
            {
                DataColumn col;
                col=results.Columns.Add("Output", typeof(string));
                col.AllowDBNull = true;
                col=results.Columns.Add("Best Model", typeof(string));
                col.AllowDBNull = true;
                col=results.Columns.Add("Distance To Native", typeof(double));
                col.AllowDBNull = true;
                col=results.Columns.Add("Best found", typeof(string));
                col.AllowDBNull = true;
                col=results.Columns.Add("Distance to best model", typeof(double));
                col.AllowDBNull = true;

                maxV = output.Count;
                foreach (var item in output)
                {
                    string[] kkk = item.dirName.Split(Path.DirectorySeparatorChar);
                    List<List<string>> clusters = item.GetClusters(10);

                    List<string> list = new List<string>();

                    foreach (var cl in clusters)
                        foreach (var target in cl)
                            list.Add(target);
                  
                    string native = item.dirName + ".pdb";
                    if (item.dirName.Contains("."))
                    {
                        string[] aux = item.dirName.Split('.');
                        native = aux[0] + ".pdb";
                    }



                    string fileBestModel = item.dirName + "_" + distance.ToString() + "_bestmodel.dat";
                    BestModel bestToNative;
                    if (!File.Exists(fileBestModel))
                    {
                        bestToNative = FindBestToModel(list, distance, item.dirName, "", "", native);
                        StreamWriter st = new StreamWriter(fileBestModel);
                        st.WriteLine(bestToNative.modelName + " " + bestToNative.value);
                        st.Close();
                    }
                    else
                    {
                        StreamReader st = new StreamReader(fileBestModel);
                        string line = st.ReadLine();
                        string[] aux = line.Split(' ');
                        bestToNative.modelName = aux[0];
                        bestToNative.value = Convert.ToInt32(aux[1], System.Globalization.CultureInfo.InvariantCulture);
                        st.Close();
                    }

                    if (bestToNative.value == int.MaxValue)
                        continue;

                   // list.Clear();

                    if (clusters.Count > 0)
                        clusters.Sort(delegate(List<string> first, List<string> second)
                        { return second.Count.CompareTo(first.Count); });

                    string represent = "";
                    //And now best five
                    
                    BestModel vOut = new BestModel();
                    InitBestModel(distance, ref vOut);
                    end = clustersNum;
                    if (clusters.Count < clustersNum)
                        end = clusters.Count;
                    for (int i = 0; i < end; i++)
                    //for (int i = 0; i < clusters.Count; i++)
                    {
                        if (profileName != null)
                            represent = CLusterRepresentJury(item.dirName, clusters[i], profileName)[0].Key;
                        else
                        {
                            if (item.juryLike != null)
                                represent = item.juryLike[i].Key;
                            else
                            {
                              /*  if (clusters.Count > 1)
                                {
                                    represent = CLusterRepresent(distance,item.dirName, clusters[i],"",profileName);                             
                                }
                                else*/
                                    represent = clusters[i][0];
                            }
                        }
                        int vDist = GetDist(distance, item.dirName, bestToNative.modelName, represent);
                        //int vDist = GetDist(distance, item.dirName, native, represent);
                        CheckBest(distance, vDist, represent, ref vOut);

                    }

                    if (bestToNative.value == errorValue || vOut.value == errorValue)
                    {
                        if (bestToNative.value == errorValue && vOut.value != errorValue)
                            results.Rows.Add(item.name, bestToNative.modelName, double.NaN, vOut.modelName, vOut.value / 100.0);
                        else
                            if(bestToNative.value != errorValue && vOut.value == errorValue)
                                results.Rows.Add(item.name, bestToNative.modelName, bestToNative.value / 100.0, vOut.modelName, double.NaN);
                            else
                                results.Rows.Add(item.name, bestToNative.modelName, double.NaN, vOut.modelName, double.NaN);
                    }
                    else
                        results.Rows.Add(item.name, bestToNative.modelName, bestToNative.value / 100.0, vOut.modelName, vOut.value / 100.0);
                    if (vOut.value!=errorValue && vOut.value-bestToNative.value < distThreshold * 100)
                    //if (vOut.value != errorValue && vOut.value < distThreshold * 100)
                    {
                        fraction++;
                    }
                    if (vOut.value!=errorValue)
                    {
                        counter++;
                        avrDiff += vOut.value;
                    }
                    //dist = null;
                    GC.Collect();

                    currentV++;

                }
                avrDiff /= counter;
                fraction /= counter;
                results.Rows.Add(null, null, null, null, null);
                results.Rows.Add("AVR:"," "+Math.Round(avrDiff / 100,2),null,"Fraction:",Math.Round(fraction,2));
            }
            catch(Exception ex)
            {
                exc = ex;
            }
           //Console.WriteLine("avrDiff=" + avrDiff + " fraction=" + fraction+" pearson="+Pearson());
            currentV = maxV;
        }

        double Pearson()
        {
            double pearson = 0;
            int m = 0;
            StreamWriter st = new StreamWriter("test.txt");
            foreach (var item in output)
            {
                double x = 0, y = 0, xy = 0 ,x2=0,y2=0;
                
                List<string> targets = new List<string>();
                foreach (var t in item.juryLike)
                    targets.Add(t.Key);
               
                string native = item.dirName+ ".pdb";

                targets.Add(Path.GetFileName(native));

                List<string> targetsFull = new List<string>();
                foreach (var t in item.juryLike)
                    targetsFull.Add(item.dirName+Path.DirectorySeparatorChar+t.Key);

                targetsFull.Add(native);

                //DistanceMeasure dist= new Rmsd(targetsFull, "", false, PDB.PDBMODE.ONLY_CA);
               // DistanceMeasure dist = new MaxSub(targetsFull, "", false);
                Dictionary<string, double> dic = ReadDistanceFile("H:\\tasser_dist\\"+Path.GetFileName(item.dirName)+".outDist");
                int n = 0;
                foreach (var it in item.juryLike)
                {
                    //double dd= dist.GetDistance(Path.GetFileName(native), it.Key);
                    if (!dic.ContainsKey(it.Key))
                        continue;

                    double dd = dic[it.Key];
                    if (dd < int.MaxValue)
                    {
                       // dd = 1 - dd;
                        x += it.Value;
                        y += dd;
                        xy += it.Value * dd;
                        x2 += it.Value * it.Value;
                        y2 += dd * dd;
                        //st.WriteLine(dd+" "+it.Value);
                        n++;
                    }
                }
                if (n > 0 && (n * x2 - x * x) * (n * y2 - y * y)>0)
                {
                    double pear = (n * (xy) - x * y) / Math.Sqrt((n * x2 - x * x) * (n * y2 - y * y));
                    st.WriteLine("pearson=" + pear);
                    pearson += pear;
                    m++;
                }

            }
            st.Close();
            return pearson/m;
        }


        int GetDist(DistanceMeasure distance, string item1, string item2)
        {
            string[] aux1 = item1.Split(Path.DirectorySeparatorChar);
            string[] aux2 = item2.Split(Path.DirectorySeparatorChar);
            int v = distance.GetDistance(aux1[aux1.Length - 1], aux2[aux2.Length - 1]);
            if (v == distance.errorValue)
                v = errorValue;
            return v;
        }

       int GetDist(DistanceMeasures distance,string dirName,string item1,string item2)
       {
           DistanceMeasure dist = null;
           List<string> targets = new List<string>();


           if(!item1.Contains(dirName))
            targets.Add(dirName + Path.DirectorySeparatorChar + item1);
           else
            targets.Add(item1);
           targets.Add(dirName + Path.DirectorySeparatorChar + item2);

           switch (distance)
           {
               case DistanceMeasures.HAMMING:
                   dist = new JuryDistance(targets, "", false, "");
                   break;
               case DistanceMeasures.MAXSUB:
                   dist = new MaxSub(targets, "", false);
                   break;
               case DistanceMeasures.RMSD:
                   dist = new Rmsd(targets, "", false, PDB.PDBMODE.ONLY_CA);
                   break;
               case DistanceMeasures.GDT_TS:
                   dist = new GDT_TS(targets, "", false);
                   break;
           }
           dist.InitMeasure();          
           return GetDist(dist, item1, item2);
       }
       DistanceMeasure PrepareDistance(DistanceMeasures distance,List <string> targets,string alignFile,string profileName)
        {
            DistanceMeasure dist=null;
          
            switch (distance)
            {
                case DistanceMeasures.HAMMING:
                    dist = new JuryDistance(targets, alignFile, false, profileName);
                    break;
                case DistanceMeasures.MAXSUB:
                    dist = new MaxSub(targets, alignFile, false);
                    break;
                case DistanceMeasures.RMSD:
                    dist = new Rmsd(targets, alignFile, false, PDB.PDBMODE.ONLY_CA);
                    break;
                case DistanceMeasures.GDT_TS:
                    dist = new GDT_TS(targets, alignFile, false);
                    break;
            }

            dist.InitMeasure();
            return dist;
        }
        void InitBestModel(DistanceMeasures distance,ref BestModel model)
       {
           switch (distance)
           {
               case DistanceMeasures.HAMMING:
              
                   model.value = int.MinValue;
                   break;
               case DistanceMeasures.MAXSUB:
               case DistanceMeasures.RMSD:
               case DistanceMeasures.GDT_TS:
                   model.value = int.MaxValue;
                   break;
           }


       }
       string CLusterRepresent(DistanceMeasure distance, List<string> targets)
       {           
           return distance.GetReferenceStructure(targets);
       }

       string CLusterRepresent(DistanceMeasures distance,string dirName, List<string> targets,string alignFile, string profileName)
        {
            DistanceMeasure dist;
            List<string> fileNames=new List<string>();

            foreach (var item in targets)
                fileNames.Add(dirName + Path.DirectorySeparatorChar + item);
            dist = PrepareDistance(distance, fileNames, alignFile, profileName);
            return dist.GetReferenceStructure(targets);
        }
       List<KeyValuePair<string,double>> CLusterRepresentJury(string dirName, List<string> targets,string profileName)
       {
                   
           List<string> fileNames = new List<string>(targets.Count);

           foreach (var item in targets)
               fileNames.Add(dirName + Path.DirectorySeparatorChar + item);

           jury1D jury = new jury1D();
           jury.PrepareJury(fileNames, null, profileName); 
           ClusterOutput opt = jury.JuryOptWeights(targets);

            //StreamWriter kk;

            //if (File.Exists("cccc.txt"))
            //    kk = File.AppendText("cccc.txt");
            //else
            //    kk = File.CreateText("cccc.txt");

            //kk.WriteLine("dirNamae=" + dirName + " size=" + opt.juryLike.Count);

            //for (int i = 0; i < opt.juryLike.Count; i++)
            //    kk.WriteLine(opt.juryLike[i].Key + " " + opt.juryLike[i].Value);

            //kk.Close();
           return opt.juryLike;
       }
       void CheckBest(DistanceMeasures distance,int value,string item,ref BestModel best)
       {
           switch (distance)
           {
               case DistanceMeasures.HAMMING:        
                   if (best.value <= value)
                   {
                       best.value = value;
                       best.modelName = item;
                   }
                   break;
               case DistanceMeasures.MAXSUB:
               case DistanceMeasures.GDT_TS:
               case DistanceMeasures.RMSD:
                   if (best.value >= value)
                   {
                       best.value=value;
                       best.modelName = item;
                   }
                   break;
           }

       }
        BestModel FindBestToModel(List<string> targets,DistanceMeasures distance,string dirName,string alignFile,string profileName,string model)
        {            
            BestModel vOut=new BestModel();
            List<string> localNames = new List<string>(2);
            List<string> fileNames = new List<string>(targets.Count);

            foreach (var item in targets)
                fileNames.Add(dirName + Path.DirectorySeparatorChar + item);

            fileNames.Add(model);

            //dist=PrepareDistance(distance,fileNames,alignFile,profileName);

            string[] aux = model.Split(Path.DirectorySeparatorChar);

            InitBestModel(distance, ref vOut);
            foreach (var item in targets)
            {
               /* localNames.Clear();
                localNames.Add(model);
                localNames.Add(item);
                dist = PrepareDistance(distance, fileNames, alignFile, profileName);*/
                int value=0;
                //value = dist.GetDistance(aux[aux.Length - 1], item);
                value = GetDist(distance, dirName, model, item);
                CheckBest(distance, value, item, ref vOut);
            }

            return vOut;
        }
        Dictionary <string,double> ReadDistanceFile(string fileName)
        {
            Dictionary<string, double> dic = new Dictionary<string, double>();
            StreamReader st = new StreamReader(fileName);
            string line = st.ReadLine();

            while (line != null)
            {
                string[] tmp = line.Split(' ');
                dic.Add(Path.GetFileName(tmp[0]), Convert.ToDouble(tmp[1]));
                line = st.ReadLine();
            }

            st.Close();
            return dic;
        }    
        
       
    }
    
}
