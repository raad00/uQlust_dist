using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using uQlustCore.Interface;
using uQlustCore.Distance;
using System.IO;
using System.Data;

namespace uQlustCore
{
    class Retrival : IProgressBar
    {
        RetrivalInput input;
        DistanceMeasure dist;
        Alignment alBase = null;
        List <string> data;
        List<string> database;

        Settings set;
        int currentV = 0, hcurrentV = 0;
        int maxV = 1, hmaxV = 1;
        public Retrival(RetrivalInput input)
        {
            this.input = input;
         
            set = new Settings();
            set.Load();
        }
        public override string ToString()
        {
            return "Retrival";
        }
        public double ProgressUpdate()
        {
            double res = 0.5 * (double)hcurrentV / hmaxV;
            if (dist!=null)
                res += 0.5*dist.ProgressUpdate();


            return res ;
        }
        public void PrepareRetrival()
        {
            alBase = new Alignment();
            Alignment alRetrival = new Alignment();

            alBase.Prepare(input.baseDir, set,input.profileName);
            alRetrival.Prepare(input.retrivalDir, set, input.profileName);

            database = alBase.GetStructuresNames();
            data = alRetrival.GetStructuresNames();


            alBase.AddProfiles(alRetrival);

            dist=JobManager.CreateMeasure(input.measure,alBase, false);
            dist.InitMeasure();
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }
        public ClusterOutput RunRetrival()
        {
            
            hmaxV = data.Count;
            jury1D jury = new jury1D();
            ClusterOutput res = new ClusterOutput();
            res.retrival = new List<List<KeyValuePair<string, double>>>();
            jury.PrepareJury(alBase);
            ClusterOutput clRes=jury.JuryOptWeights(alBase.GetStructuresNames());
            List<string> juryOrder = new List<string>();
            foreach (var item in clRes.juryLike)
                juryOrder.Add(item.Key);

            Dictionary<string, int> dataBaseDic = new Dictionary<string, int>();
            foreach (var item in database)
                dataBaseDic.Add(item, 0);
            int[] index;// = new int[database.Count];
            List<string> selected = new List<string>();
            int[] distTab;
            StreamWriter wr = new StreamWriter("retrival_tmp.dat");
            for(int i=0;i<data.Count;i++)
            {
                selected.Clear();
                /*int indexF = juryOrder.FindIndex(x=>x==data[i]);
                for (int j = Math.Max(0, indexF - database.Count / 6); j < Math.Min(database.Count, indexF + database.Count / 6); j++)
                    if (dataBaseDic.ContainsKey(juryOrder[j]))
                        selected.Add(juryOrder[j]);
                distTab = dist.GetDistance(data[i], selected);*/
                distTab = dist.GetDistance(data[i], database);
                index = new int[distTab.Length];
                for (int j = 0; j < index.Length; j++)
                    index[j] = j;

                Array.Sort(distTab, index);
                wr.WriteLine("Next i=" + i);
                List<KeyValuePair<string, double>> aux = new List<KeyValuePair<string, double>>();
                for (int j = 0; j < input.numToRetrive; j++)
                {
                    aux.Add(new KeyValuePair<string, double>(database[index[j]], distTab[j]));
                    //aux.Add(new KeyValuePair<string, double>(selected[index[j]], distTab[j]));
                    //wr.WriteLine(selected[index[j]]+" " +distTab[j]);
                    wr.WriteLine(database[index[j]] + " " + distTab[j]);
                }
                wr.WriteLine();
                res.retrival.Add(aux);
                currentV++;
            }
            wr.Close();
            return res;
        }
    }
}
