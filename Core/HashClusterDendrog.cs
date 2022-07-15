using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using uQlustCore;
using System.Collections.Concurrent;
using System.Threading;
using System.Threading.Tasks;
using System.Diagnostics;
using uQlustCore.Interface;
using System.Data;
using uQlustCore.Distance;

namespace uQlustCore
{
    class HashClusterDendrog:HashCluster,IProgressBar
    {
         DistanceMeasures dMeasure;
         DistanceMeasure dist=null;
         AglomerativeType linkageType;
         PDB.PDBMODE atoms;
         bool jury1d;
         string profileName;
         string refJuryProfile;
         string dirName;
         HierarchicalCInput hier;
         hierarchicalCluster hk = null;
         public HashClusterDendrog(DCDFile dcd, HashCInput input,HierarchicalCInput dendrogOpt):base(dcd,input)
        {
            this.dMeasure=dendrogOpt.distance;
            this.linkageType=dendrogOpt.linkageType;
            this.atoms = dendrogOpt.atoms;
            this.jury1d = dendrogOpt.reference1DjuryH;
            this.profileName = dendrogOpt.hammingProfile;
            this.refJuryProfile = dendrogOpt.jury1DProfileH;
            hier = dendrogOpt;
        }
         public HashClusterDendrog(string dirName, string alignFile, HashCInput input, HierarchicalCInput dendrogOpt)
             : base(dirName, alignFile, input)
        {
            this.dMeasure=dendrogOpt.distance;
            this.linkageType=dendrogOpt.linkageType;
            this.atoms = dendrogOpt.atoms;
            this.jury1d = dendrogOpt.reference1DjuryH;
            this.profileName = dendrogOpt.hammingProfile;
            this.refJuryProfile = dendrogOpt.jury1DProfileH;
            this.dirName = dirName;
            hier = dendrogOpt;
        }
        public void InitHashClusterDendrog()
         {
             base.InitHashCluster();
         }
        public new double ProgressUpdate()
        {
            double progress = 0;

            progress = 0.5*base.ProgressUpdate();

            if (hk != null)
                progress += 0.5 * hk.ProgressUpdate() ;

            return progress;
        }
        public new Exception GetException()
        {
            return null;
        }
        public new List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }

        public string UsedMeasure()
         {
             if (dist != null)
                 return dist.ToString();

             return "NONE";
         }
         public override string ToString()
         {        
             return "uQlust:Tree";
         }
         public ClusterOutput RunHashDendrogCombine()
         {
             ClusterOutput output = DendrogUsingMeasures(stateAlignKeys);
             return output;
         }
        public Dictionary<string,int> ReadLeafs()
        {
            Dictionary<string, int> leafs = new Dictionary<string, int>();
            StreamReader r = new StreamReader("C:\\Projects\\results_leavs.dat");

            string line = r.ReadLine();

            while(line!=null)
            {
                if(line.Contains("cluster"))
                {
                    string[] aux = line.Split(' ');
                    if (!leafs.ContainsKey(aux[0]))
                        leafs.Add(aux[0], 1);
                }
                line = r.ReadLine();
            }


            r.Close();
            return leafs;
        }
        Dictionary<string,List<int>> SelectClusters(Dictionary<string,int> rep, Dictionary<string,List<int>> clusters)
        {
            Dictionary<string, List<int>> res = new Dictionary<string, List<int>>();
            Dictionary<int,string> toString = new Dictionary<int,string>();

            for (int i = 0; i < allStructures.Count; i++)
            {
                string name = allStructures[i];
                name = name.Replace(".pdb", "");
                name = name.Replace("pdb", "");
                string[] aux = name.Split('.');
                string[] xx = name.Split('_');
                string ww = aux[0] + xx[xx.Length - 1].ToLower();
                toString.Add(i, ww);
            }
            StreamWriter cc = new StreamWriter("represent.dat");
            foreach(var item in clusters)
            {
                foreach(var cl in item.Value)
                    if(rep.ContainsKey(toString[cl]))
                    {
                        cc.WriteLine(toString[cl]);
                        List<int> list = new List<int>();
                        list.Add(cl);
                        foreach(var k in item.Value)
                        {
                            if (k != cl)
                                list.Add(k);
                        }
                        res.Add(item.Key, list);
                        break;
                    }
            }
            cc.Close();

            return res;
        }
         public ClusterOutput DendrogUsingMeasures(List<string> structures)
         {
             jury1D juryLocal = new jury1D();
             juryLocal.PrepareJury(al);
             
             ClusterOutput outC = null;
             Dictionary<string, List<int>> dic;
             //Console.WriteLine("Start after jury " + Process.GetCurrentProcess().PeakWorkingSet64);
             maxV = refPoints * 20*4;
             currentV = 0;
             dic = PrepareKeys(structures,false);
                          
             //DebugClass.DebugOn();
           //  input.relClusters = input.reqClusters;
           //  input.perData = 90;
             if (dic.Count > input.relClusters)
             {
                 if (!input.rpart)
                     dic = HashEntropyCombine(dic, structures, input.relClusters);
                 else
                    dic = Rpart(dic,null, structures, false);
                    //dic = FastCombineKeysNew(dic, structures, false);
            }
            //Dictionary<string, int>  xx=ReadLeafs();
            //dic = SelectClusters(xx, dic);
           maxV = 3;
            currentV = 1;
             //Console.WriteLine("Entropy ready after jury " + Process.GetCurrentProcess().PeakWorkingSet64);
             DebugClass.WriteMessage("Entropy ready");
             //Alternative way to start of UQclust Tree must be finished
             //input.relClusters = 10000;

             
             //dic = FastCombineKeys(dic, structures, true);
             DebugClass.WriteMessage("dic size" + dic.Count);
             currentV++;
             
             //Console.WriteLine("Combine ready after jury " + Process.GetCurrentProcess().PeakWorkingSet64);
             DebugClass.WriteMessage("Combine Keys ready");
             Dictionary<string, string> translateToCluster = new Dictionary<string, string>(dic.Count);
             List<string> structuresToDendrogram = new List<string>(dic.Count);
             List<string> structuresFullPath = new List<string>(dic.Count);
             DebugClass.WriteMessage("Number of clusters: "+dic.Count);
             int cc = 0;
             List<string> order = new List<string>(dic.Keys);
             order.Sort(delegate(string a, string b)
             {

                 if (dic[b].Count == dic[a].Count)
                     for (int i = 0; i < a.Length; i++)
                         if (a[i] != b[i])
                             if (a[i] == '0')
                                 return -1;
                             else
                                 return 1;


                 return dic[b].Count.CompareTo(dic[a].Count);
             });
             foreach (var item in order)
             {
                 if (dic[item].Count > 2)
                 {
                     List<string> cluster = new List<string>(dic[item].Count);
                     foreach (var str in dic[item])
                         cluster.Add(structures[str]);


                     ClusterOutput output = juryLocal.JuryOptWeights(cluster);

                     structuresToDendrogram.Add(output.juryLike[0].Key);
                     if(alignFile==null)
                        structuresFullPath.Add(dirName + Path.DirectorySeparatorChar + output.juryLike[0].Key);
                     else
                         structuresFullPath.Add(output.juryLike[0].Key);
                     translateToCluster.Add(output.juryLike[0].Key, item);
                 }
                 else
                 {
                     structuresToDendrogram.Add(structures[dic[item][0]]);
                     if(alignFile==null)
                        structuresFullPath.Add(dirName + Path.DirectorySeparatorChar + structures[dic[item][0]]);
                     else
                         structuresFullPath.Add(structures[dic[item][0]]);
                     translateToCluster.Add(structures[dic[item][0]], item);
                 }
                 cc++;
             }
             currentV++;
             DebugClass.WriteMessage("Jury finished");
             switch (dMeasure)
             {
                 case DistanceMeasures.HAMMING:
                     if (refJuryProfile == null || !jury1d)
                         throw new Exception("Sorry but for jury measure you have to define 1djury profile to find reference structure");
                     else
                         dist = new JuryDistance(structuresFullPath, alignFile, true, profileName, refJuryProfile);
                     break;


                 case DistanceMeasures.COSINE:
                     dist = new CosineDistance(structuresFullPath, alignFile, jury1d, profileName, refJuryProfile);
                     break;

                 case DistanceMeasures.RMSD:
                     dist = new Rmsd(structuresFullPath, "", jury1d, atoms, refJuryProfile);
                     break;
                 case DistanceMeasures.MAXSUB:
                     dist = new MaxSub(structuresFullPath, "", jury1d, refJuryProfile);
                     break;
             }

            // return new ClusterOutput();
             DebugClass.WriteMessage("Start hierarchical");
             //Console.WriteLine("Start hierarchical " + Process.GetCurrentProcess().PeakWorkingSet64);
             currentV = maxV;
             hk = new hierarchicalCluster(dist, hier, dirName);
             dist.InitMeasure();
            
             //Now just add strctures to the leaves             
             outC = hk.HierarchicalClustering(structuresToDendrogram);
             DebugClass.WriteMessage("Stop hierarchical");
             List<HClusterNode> hLeaves = outC.hNode.GetLeaves();

             foreach(var item in hLeaves)
             {
                 if (translateToCluster.ContainsKey(item.setStruct[0]))
                 {
                     foreach (var str in dic[translateToCluster[item.setStruct[0]]])
                         if (item.setStruct[0] != structures[str])
                            item.setStruct.Add(structures[str]);

                    item.consistency = CalcClusterConsistency(item.setStruct);
                }
                 else
                     throw new Exception("Cannot add structure. Something is wrong");
             }
             outC.hNode.RedoSetStructures();
             outC.runParameters = hier.GetVitalParameters();
             outC.runParameters += input.GetVitalParameters();
             return outC;
         }
         

    }
}
