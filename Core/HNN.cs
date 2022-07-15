using System;
using System.Text;
using System.Collections.Generic;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using uQlustCore.Profiles;
using System.Linq;
using System.IO;
using uQlustCore.Interface;
using System.Text.RegularExpressions;

namespace uQlustCore
{
    [Serializable]
    public class HNN:ISerialize
    {
        [NonSerialized]
        Dictionary<string, string> classLabels;
        [NonSerialized]
        protected Dictionary<string, string> caseBase = new Dictionary<string,string>();
        [NonSerialized]
        protected Dictionary<string, string> testData = new Dictionary<string, string>();
        [NonSerialized]
        Dictionary<string, string> labelToBaseKey = new Dictionary<string, string>();
        [NonSerialized]
        public List<string> validateList = new List<string>();
        [NonSerialized]
        public List<string> testList = new List<string>();
        [NonSerialized]
        public List<string> clusterLabels = null;
        [NonSerialized]
        protected HashCluster hk=null;
        [NonSerialized]
        ClusterOutput outP=null;
        
        public ClusterOutput outCl { get { return outP; } set { outP = value; } }
        [NonSerialized]
        protected Settings set;
        [NonSerialized]
        protected HNNCInput opt;

        public HNN(HashCluster hk,ClusterOutput outp,HNNCInput opt)        
        {
            this.hk=hk;
            this.opt = opt;
            set = new Settings();
            set.Load();
            PrepareCaseBaseLabels(outp);
        }
        public HNN(Dictionary<byte, int>[] columns, HNNCInput opt)
        {
            this.opt = opt;
        }
        public HNN(HNNCInput opt)
        {
            this.opt = opt;
        }
        public override string ToString()
        {
            return "HNN";
        }
        public void ISaveBinary(string fileName)
        {
            GeneralFunctionality.SaveBinary(fileName, this);
        }
        public ISerialize ILoadBinary(string fileName)
        {
            return GeneralFunctionality.LoadBinary(fileName);
        }
        static Dictionary<string,string> ReadClassLabels(string fileName)
        {
            Dictionary<string, string> labels = new Dictionary<string, string>();
            StreamReader wr = new StreamReader(fileName);
            string line = wr.ReadLine();
            while(line!=null)
            {
                string[] aux = line.Split(' ');
                if (aux.Length == 2)
                    labels.Add(aux[0], aux[1]);

                line = wr.ReadLine();
            }
            wr.Close();

            return labels;
        }
        void PrepareCaseBaseLabels(ClusterOutput outp)
        {
            if (outp == null)
                return;

            classLabels = null;
            if (opt.labelsFile.Length > 0)
            {
                classLabels = ReadClassLabels(opt.labelsFile);
               
            }
            clusterLabels = new List<string>();
            if (classLabels == null || classLabels.Count == 0)
            {
                int num = 1;
                this.classLabels = new Dictionary<string, string>();
                foreach (var cluster in outp.clusters.list)
                {
                    string label = "cluster_" + num++;
                    foreach(var item in cluster)                    
                        this.classLabels.Add(item, label);                    
                }
            }
            Dictionary<string, int> classDic = new Dictionary<string, int>();
            foreach (var item in hk.dicFinal)
            {
                if (item.Value.Count > 0)
                {
                    classDic.Clear();
                    foreach (var index in item.Value)
                    {

                        string structV = hk.structNames[index];
                        if (this.classLabels.ContainsKey(structV))
                        {
                            if (classDic.ContainsKey(this.classLabels[structV]))
                                classDic[this.classLabels[structV]]++;
                            else
                                classDic.Add(this.classLabels[structV], 1);
                        }
                    }
                    foreach (var index in item.Value)
                        if (labelToBaseKey.ContainsKey(hk.structNames[index]))
                            //throw new Exception("Key " + hk.structNames[index] + " already exists");
                            Console.Write("Key " + hk.structNames[index] + " already exists");
                        else
                            labelToBaseKey.Add(hk.structNames[index], item.Key);
                    if (classDic.Count == 0)
                        continue;
                    List<string> classLab = new List<string>(classDic.Keys);

                    classLab.Sort((x, y) => classDic[x].CompareTo(classDic[y]));
                    caseBase.Add(item.Key, classLab[0]);
                    clusterLabels.Add(classLab[0]);

                }
                else
                    if (labelToBaseKey.ContainsKey(hk.structNames[item.Value[0]]))
                    throw new Exception("Key " + hk.structNames[item.Value[0]] + " already exists");
                else
                    labelToBaseKey.Add(hk.structNames[item.Value[0]], item.Key);
            }


        }
        public double HNNValidate(List<string> validList)
        {
            int good = 0;
            int all = 0;
            double acc=0;
            Dictionary<string, int> classDic = new Dictionary<string, int>();
            foreach (var vItem in validList)
            {
                if(labelToBaseKey.ContainsKey(vItem) && classLabels.ContainsKey(vItem))
                {
                    if (classLabels[vItem] == caseBase[labelToBaseKey[vItem]])
                        good++;
                    all++;
                }
                
            }
            acc = (double)good / all;
            return acc;
        }
        public List<string> ITest(string fileName)
        {
            Dictionary<string, protInfo<byte>> test = new Dictionary<string, protInfo<byte>>();
            if (set.mode == INPUTMODE.USER_DEFINED)
            {
                UserDefinedProfile<byte> p = new UserDefinedProfile<byte>();
                List<string> nodes = new List<string>(hk.al.r.masterNode.Keys);
                test = p.GetProfile(hk.al.r.masterNode[nodes[0]], fileName,null);
            }
            
            foreach (var item in test)
            {
                List<string> aux = new List<string>();
                foreach (var state in item.Value.profile)
                    aux.Add(state.ToString());
                
                hk.al.r.AddItemsCombineStates(item.Key, aux);
            }
            
            return new List<string>(test.Keys);
            //Dictionary<string, string> res = HNNTest(new List<string>(test.Keys));

            //return res;
        }
        protected Dictionary<string,string> BasesForTest(List<string> testList)
        {
            Dictionary<string, string> aux = new Dictionary<string, string>();
            

            Dictionary<string, List<int>> kk = hk.PrepareKeys(testList, true, false);
            List<string> keyList = new List<string>(kk.Keys);
            List<int> validIndexes = null;
            if (hk.validIndexes.Count > 0)
            {
                validIndexes = hk.validIndexes;

                for (int j = 0; j < keyList.Count; j++)
                {
                    StringBuilder newOrder = new StringBuilder();
                    for (int i = 0; i < validIndexes.Count; i++)
                        newOrder.Append(keyList[j][validIndexes[i]]);

                    string newKey = newOrder.ToString();
                    if (kk.ContainsKey(newKey))
                        kk[newKey].AddRange(kk[keyList[j]]);
                    else
                    {
                        List<int> indx = new List<int>();
                        indx.AddRange(kk[keyList[j]]);
                        kk.Add(newKey, indx);
                    }
                    kk.Remove(keyList[j]);
                }
            }
            List<string> caseKeys = new List<string>(kk.Keys);
            foreach (var item in caseKeys)
            {
                /*  if (caseBase.ContainsKey(item))
                  {
                      foreach (var v in kk[item])
                          if (!res.ContainsKey(testList[v]))
                              res.Add(testList[v], caseBase[item]);
                  }
                  else*/
                {
                    foreach (var v in kk[item])
                        if (!aux.ContainsKey(item))
                            aux.Add(item, testList[v]);
                }
            }
            return aux;
        }
        public Dictionary<string, string> GetTestProfileFromPDB(Options opt)
        {
            Alignment al;
            al = new Alignment();
            al.Prepare(opt.hnn.testFile, set, opt.hash.profileName);
            al.MyAlign(null);
            Dictionary<string, List<byte>> aux = al.GetStateAlign();
            foreach (var item in aux)
            {
                StringBuilder st = new StringBuilder();
                foreach (var s in item.Value)
                    st.Append((char)(s+'0'));

                testData.Add(item.Key,st.ToString());
            }
            return testData;
        }

        public virtual void PairwiseDistance(string listFileName, string distOut)
        {
            StreamReader str = new StreamReader(listFileName);
            string line = str.ReadLine();
            List<string> lista = new List<string>();
            while(line!=null)
            {
                lista.Add(line);
                line = str.ReadLine();
            }
            str.Close();
            StreamWriter wr = new StreamWriter(distOut);
            foreach (var item in lista)
            {
                double common = 0;
                double all = 0;
                int hamming = 0;
                string[] aux = item.Split(' ');
                if (!caseBase.ContainsKey(aux[0]))
                    continue;
                string str1 = caseBase[aux[0]];
                if (!caseBase.ContainsKey(aux[1]))
                    continue;

                string str2 = caseBase[aux[1]];
                double iJacard = 0;
                for (int i = 0; i < str1.Length; i++)
                {
                    if (str1[i] == '2' && str1[i] == str2[i])
                        common++;
                    if (str1[i] == '2' || str2[i] == '2')
                        all++;
                    if (str1[i] != str2[i])
                        hamming++;
                }
                iJacard = common / all;
                wr.WriteLine(item + " " + iJacard+" "+hamming);
            }
            wr.Close();
        }
        public virtual Dictionary<string, string> HNNTest(List<string> testList)
        {
            Dictionary<string, string> res = new Dictionary<string, string>();
            Dictionary<string,string> aux= BasesForTest(testList);
            if (aux.Count > 0)
            {
                StringBuilder final = new StringBuilder();
                Dictionary<string, List<string>> keys = hk.AddToClusters(new List<string>(caseBase.Keys), new List<string>(aux.Keys));
                foreach (var item in keys.Keys)
                {
                    final.Clear();
                    if (keys[item].Count > 0)
                    {
                        Dictionary<string, int> c = new Dictionary<string, int>();
                        for (int i = 0; i < keys[item].Count; i++)
                            if (!c.ContainsKey(caseBase[keys[item][i]]))
                                c.Add(caseBase[keys[item][i]], 1);
                        List<string> xx = new List<string>(c.Keys);
                        if (xx.Count > 100)
                        {
                            jury1D j = new jury1D(null, hk.al);
                            j.PrepareJury(hk.al);
                            ClusterOutput r = j.JuryOptWeights(xx, null);
                            int remN = -1;
                            for (int n = 0; n < r.juryLike.Count; n++)
                            {
                                if (r.juryLike[n].Key == caseBase[item])
                                {
                                    remN = n;
                                    break;
                                }
                            }
                            List<double> tmp = new List<double>();
                            for (int n = 0; n < r.juryLike.Count; n++)
                            {
                                double s = Math.Abs(r.juryLike[remN].Value - r.juryLike[n].Value);
                                tmp.Add(s);
                            }

                            var sorted = tmp.Select((x, i) => new { Value = x, OriginalIndex = i }).OrderBy(x => x.Value).ToList();

                            for (int i = 0; i < r.juryLike.Count - 1; i++)
                                final.Append(r.juryLike[sorted[i].OriginalIndex].Key + ":");
                            final.Append(r.juryLike[sorted[xx.Count - 1].OriginalIndex].Key);


                        }
                        else
                        {
                            for (int i = 0; i < xx.Count - 1; i++)
                                final.Append(xx[i] + ":");
                            final.Append(xx[xx.Count - 1]);
                        }
                    }
                    if (final.Length == 0)
                        res.Add(aux[item], "NOT CLASSIFIED");
                    else
                        res.Add(aux[item], final.ToString());
                }
            }

            return res;
        }
        public virtual void Preprocessing(bool Load)
        {

        }
       
       
        public void SelectFeatures()
        {
            StreamReader wr = new StreamReader("C:\\Projects\\cath4.2-domain-list.txt");
            Dictionary<string, List<string>> cathClusters = new Dictionary<string, List<string>>();
            string line = wr.ReadLine();
            while(line!=null)
            {
                line = line.Replace('\t', ' ');
                line=Regex.Replace(line, @"\s+", " ");
                string[] aux = line.Split(' ');
                if(cathClusters.ContainsKey(aux[1] + "." + aux[2] + "." + aux[3] + "." + aux[4]))                
                    cathClusters[aux[1] + "." + aux[2] + "." + aux[3] + "." + aux[4]].Add(aux[0]);
                else
                {
                    List<string> x = new List<string>();
                    x.Add(aux[0]);
                    cathClusters.Add(aux[1] + "." + aux[2] + "." + aux[3] + "." + aux[4], x);
                }
                
                line = wr.ReadLine();
            }
            wr.Close();
            List<int> keep = new List<int>();
            StreamWriter re = new StreamWriter("features");
            int[] select = new int[hk.validIndexes.Count];
            foreach (var item in cathClusters)
            {
                List<string> clusterKeys = new List<string>();
                List<string> nameKey = new List<string>();
                foreach (var lab in item.Value)
                    if (labelToBaseKey.ContainsKey(lab))
                    {

                        clusterKeys.Add(labelToBaseKey[lab]);
                        nameKey.Add(lab);
                    }

                if (clusterKeys.Count == 0)
                    continue;

                re.WriteLine("Cluster size: " + clusterKeys.Count);
                int counter = 0;
                
                keep.Clear();
                for (int i = 0; i < clusterKeys[0].Length; i++)
                {
                    bool test = false;
                    
                    for (int j=1;j<clusterKeys.Count;j++)
                    {
                            if (clusterKeys[0][i] != clusterKeys[j][i])
                            {
                                test = true;
                                break;
                            }

                    }
                    if(test)
                    {
                        select[i]++;
                        counter++;                        
                    }
                    else
                        keep.Add(i);
                }
                string refKey = "";
                for (int i = 0; i < keep.Count; i++)
                    refKey += labelToBaseKey[nameKey[0]][keep[i]];

                StringBuilder key = new StringBuilder();
                Dictionary<string, int> locCl = new Dictionary<string, int>();
                foreach (var it in labelToBaseKey.Keys)
                {
                    key.Clear();
                    for(int i=0;i<keep.Count;i++)
                    {
                        key.Append(labelToBaseKey[it][keep[i]]);
                    }
                    string xx = key.ToString();
                    if (locCl.ContainsKey(xx))
                        locCl[xx]++;
                    else
                        locCl.Add(xx,1);
                }

                if (clusterKeys.Count > locCl[refKey])
                    Console.Write("UPPP");
                re.WriteLine("liczba wyrzuconych: " + counter);
                re.WriteLine("Rozmiar clustra dla wszystkich: " + locCl[refKey]);
                re.WriteLine();
            }
            re.Close();
        }
    }
}
