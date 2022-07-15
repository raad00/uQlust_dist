using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using uQlustCore.PDB;
using uQlustCore;
using uQlustCore.Distance;
using System.Threading;
using System.Reflection;
using System.Text.RegularExpressions;
namespace uQlustCore.Profiles
{
    public class FragBagProfile:ContactProfile
    {
        protected string fragBagProfile = "FragBag profile ";
        List<float [,]> fragBagLibrary=new List<float[,]>();
        Dictionary<char, string> chainId = new Dictionary<char, string>();
        float[] distData;
        int[] index;
        int[] profile;
        int[] zeroCount;
        float[][,] chunk = null;
        float[][] center;
        List<Residue> res;
        int count = 0;

        ManualResetEvent[] resetEvents = null;
        public FragBagProfile()
        {
            dirSettings.Load();
           
            destination = new List<INPUTMODE>();
            if(dirSettings.mode==INPUTMODE.PROTEIN)
                destination.Add(INPUTMODE.PROTEIN);
            else
                destination.Add(INPUTMODE.RNA);
            profileName = "FragBag";
            contactProfile = "FragBag profile ";
            AddInternalProfiles();            
            maxV = 1;

        }
        public override void AddInternalProfiles()
        {
            profileNode node = new profileNode();

            node.profName = profileName;
            node.internalName = profileName;
            for (int i = 0; i < 255; i++)
                node.AddStateItem(i.ToString(), i.ToString(),true);

             InternalProfilesManager.AddNodeToList(node, typeof(FragBagProfile).FullName);

        }
        public override bool CheckIfAvailable()
        {
            if (!Directory.Exists(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) +
                    Path.DirectorySeparatorChar+"fragLib"))
                if (!Directory.Exists("C:\\Projects\\UQlast\\fragLib"))
                    throw new Exception("Directory fragLib not exists. Profile FragBag cannot be used!");

            ReadLibrary("fragBagProtein.txt");
            return true;
        }
        private void ReadCathFile(string fileName,Dictionary<string,Dictionary<string,string> > aux)
        {
            
            StreamReader readFile = new StreamReader(fileName);
            string[] tmp;
            string line = readFile.ReadLine();
            string name="",seq;
            string pdbName, chainName;
            while(line!=null)
            {
                if(line.Contains(">"))
                {
                    tmp = line.Split('|');
                    name =tmp[tmp.Length - 1];
                }
                else
                {
                    seq = line;
                    pdbName = "pdb"+name.Substring(0, 4)+".ent";
                    chainName = name.Substring(4, name.Length - 4);

                    if (!aux.ContainsKey(pdbName))
                    {
                        Dictionary<string, string> dic = new Dictionary<string, string>();
                        dic.Add(chainName, seq);

                        aux.Add(pdbName, dic);
                    }
                    else
                        aux[pdbName].Add(chainName,seq);
                }
                line = readFile.ReadLine();
            }

            readFile.Close();

            
        }

     /*  public override int Run(object processParams)
       {
                    Dictionary<int, List<int>> contacts = new Dictionary<int, List<int>>();
                    string fileName = ((ThreadFiles)(processParams)).fileName;

                    List<string> files = threadingList[((ThreadFiles)(processParams)).threadNumber];// CheckFile(listFile);
                    if (files.Count == 0)
                        return 0;
                    StreamWriter wr;

                    if (File.Exists(fileName))
                        wr = File.AppendText(fileName);
                    else
                        wr = new StreamWriter(fileName);

                    if (wr == null)
                        throw new Exception("Cannot open file: " + fileName);

                    PDBFiles pdbs = new PDBFiles();
                    PDBFiles pdbsNew = new PDBFiles();
                    //Needed only for CATH
                    string []cathFiles=Directory.GetFiles("F:\\CATH");
                    Dictionary<string, Dictionary<string, string>> aux = new Dictionary<string, Dictionary<string, string>>();
                    foreach(var item in cathFiles)
                    {
                        if(item.Contains("COMBS"))
                            ReadCathFile(item, aux);
                    }



                    //-----------------------


                    maxV = files.Count;
                    bool test=false;
                    try
                    {
                        foreach (var item in aux.Keys)
                        {
                            string strName;

//if (!test && item.Contains("1vij"))
                                test = true;

                            if (!test)
                                continue;
                            //strName = pdbs.AddPDB(item, PDBMODE.CA_CB);
                            strName = pdbs.AddPDB("F:\\cath_pdb\\"+item, PDBMODE.ONLY_CA,CHAIN_MODE.ALL);
                  
                         //   if(aux.ContainsKey(pdbs.molDic))

                            if (strName != null)
                            {
                                if (!aux.ContainsKey(strName))
                                    Console.WriteLine("UPS");
                                List<Chain> toRemove = new List<Chain>();
                                chainId.Clear();
                                pdbsNew.molDic.Clear();
                                MolData mol=new MolData();
                                mol.mol = new Molecule();
                                pdbsNew.molDic=new Dictionary<string,MolData>();
                                pdbsNew.molDic.Add(strName,mol);
                                char chStart = '0';
                                foreach (var chains in aux[strName].Keys)
                                {
                                    if (!chains.Contains("_"))
                                        continue;
                                    Console.WriteLine(chains);
                                    Chain cc=pdbs.molDic[strName].mol.CuttMoleculeToSEQ(aux[strName][chains], chains.Substring(0, 1)[0]);                                    
                                    if(cc!=null)
                                    {
                                        cc.ChainIdentifier = chStart;
                                        chainId.Add(chStart, chains);
                                        chStart++;
                                        pdbsNew.molDic[strName].mol.Chains.Add(cc);
                                    }
                                }                
                                if(pdbsNew.molDic[strName].mol.Chains.Count>0)
                                    MakeProfiles(strName, pdbsNew.molDic[strName], wr);
                            }
                            pdbs.molDic.Clear();
                            currentProgress++;
                        }
                    }
                    catch (Exception ex)
                    {
                        // wr.Close();
                        throw new Exception(ex.Message);
                    }

                    wr.Close();
                    currentProgress = maxV;
                    return 0;
                }*/
        public override int Run(object processParams)
        {
            Dictionary<int, List<int>> contacts = new Dictionary<int, List<int>>();
            string fileName = ((ThreadFiles)(processParams)).fileName;

            List<string> files = threadingList[((ThreadFiles)(processParams)).threadNumber];// CheckFile(listFile);
            if (files.Count == 0)
                return 0;
            StreamWriter wr;

            if (File.Exists(fileName))
                wr = File.AppendText(fileName);
            else
                wr = new StreamWriter(fileName);

            if (wr == null)
                throw new Exception("Cannot open file: " + fileName);

            PDBFiles pdbs = new PDBFiles();
            chunk = new float[dirSettings.numberOfCores][,];
            center = new float[dirSettings.numberOfCores][];
            for (int i = 0; i < dirSettings.numberOfCores; i++)
            {
                chunk[i] = new float[fragBagLibrary[0].GetLength(0), fragBagLibrary[0].GetLength(1)];
                center[i] = new float[3];
            }
            maxV = files.Count;
            try
            {
                resetEvents = new ManualResetEvent[dirSettings.numberOfCores];
                for (int i = 0; i < dirSettings.numberOfCores; i++)
                    resetEvents[i] = new ManualResetEvent(false);
                foreach (var item in files)
                {
                    string strName;
                    //strName = pdbs.AddPDB(item, PDBMODE.CA_CB);
                    if(dirSettings.mode==INPUTMODE.PROTEIN)
                        strName = pdbs.AddPDB(item, PDBMODE.ONLY_CA,CHAIN_MODE.ALL);
                    else
                        strName = pdbs.AddPDB(item, PDBMODE.ONLY_P, CHAIN_MODE.ALL);
                  
                 //   if(aux.ContainsKey(pdbs.molDic))

                    if (strName != null)
                    {
                        MakeProfiles(strName, pdbs.molDic[strName], wr);
                    }
                    pdbs.molDic.Clear();
                    ErrorBase.ClearErrors();
                    Interlocked.Increment(ref currentProgress);
                }
            }
            catch (Exception ex)
            {
                // wr.Close();
                throw new Exception(ex.Message);
            }

            wr.Close();
           // currentProgress = maxV;
            return 0;
        }

        private void ThreadLibrary(object o)
        {
            int[] inx = (int [])o;

            
            

            if (fragBagLibrary[0].GetLength(0) + inx[2] > res.Count)
            {
                resetEvents[inx[3]].Set();
                return;
            }
            
            for (int n = 0; n < fragBagLibrary[0].GetLength(0); n++)
            {
                chunk[inx[3]][n, 0] = res[inx[2] + n].Atoms[0].Position.X;
                chunk[inx[3]][n, 1] = res[inx[2] + n].Atoms[0].Position.Y;
                chunk[inx[3]][n, 2] = res[inx[2] + n].Atoms[0].Position.Z;
            }
            Optimization.CenterMol(chunk[inx[3]], center[inx[3]]);

            for (int j = inx[0]; j < inx[1]; j++)
                distData[j] = Optimization.Rmsd(fragBagLibrary[j], chunk[inx[3]], false);

            resetEvents[inx[3]].Set();
        }
        protected override void MakeProfiles(string strName, MolData molDic, StreamWriter wr)
        {
            int tNumb = dirSettings.numberOfCores;
         
            if (molDic != null)
            {
            
                foreach (var chain in molDic.mol.Chains)
                {
                    for (int i = 0; i < profile.Length; i++)
                        profile[i] = 0;
                    res = chain.Residues;
                    if(dirSettings.mode==INPUTMODE.PROTEIN)
                        if (res.Count <= 10)
                            continue;
                    if(dirSettings.mode == INPUTMODE.RNA)
                        if (res.Count <= 5)
                            continue;
                    for (int i = 0; i < res.Count; i++)
                    {
                        for (int j = 0; j < fragBagLibrary.Count; j++)
                        {
                            distData[j] = float.MaxValue;
                            index[j] = j;
                        }

                        for (int n = 0; n < tNumb; n++)
                        {
                            int[] pp = new int[4];
                            pp[0] = (int)(n * fragBagLibrary.Count / Convert.ToDouble(tNumb, System.Globalization.CultureInfo.InvariantCulture));
                            pp[1] = (int)((n + 1) * fragBagLibrary.Count / Convert.ToDouble(tNumb, System.Globalization.CultureInfo.InvariantCulture));
                            pp[2] = i;
                            pp[3] = n;
                            resetEvents[n].Reset();
                            
                            ThreadPool.QueueUserWorkItem(new WaitCallback(ThreadLibrary), (object)pp);

                        }
                        for (int n = 0; n < tNumb; n++)
                            resetEvents[n].WaitOne();


                        Array.Sort(distData, index);
                        //Console.WriteLine(distData[0] + " " + distData[1]);
                        if (distData[0] != float.MaxValue)// && distData[0]<1)
                        {
                            if (profile[index[0]] > Byte.MaxValue)
                                profile[index[0]] = Byte.MaxValue;
                            else
                                profile[index[0]]++;                             
                        }
                        
                    }
                    if (profile.Length > 0)
                    {
                        if(chainId.Count>0)
                            wr.WriteLine(">" + strName + "|" + chainId[chain.ChainIdentifier]);
                        else
                            wr.WriteLine(">" + strName);
                        wr.Write(fragBagProfile);
                        foreach (var item in profile)
                        {
                            byte value = 0;
                            if (item < Byte.MaxValue)
                                value = Convert.ToByte(item);
                            else
                                value = Byte.MaxValue;
                                wr.Write(value + " ");
                        }
                        wr.WriteLine();
                        count++;
                        for (int i = 0; i < profile.Length; i++)
                            if (profile[i] == 0)
                                zeroCount[i]++;
                    }
                }
            }
        }
        public override void JoinFiles(string fileName)
        {
            StreamWriter wr;
            Dictionary<int, int>[] columns=null;
                if (File.Exists(GetProfileFileName(fileName)))
                    wr = File.AppendText(GetProfileFileName(fileName));
                else
                    wr = new StreamWriter(GetProfileFileName(fileName));

                if (wr == null)
                    throw new Exception("Cannot open file: " + GetProfileFileName(fileName));

                for (int i = 0; i < threadNumbers; i++)
                {
                    string fileN = GetProfileFileName(fileName) + "_" + i;
                    using (StreamReader rr = new StreamReader(fileN))
                    {
                        if (rr == null)
                            throw new Exception("Cannot open file: " + fileN);

                        string line = rr.ReadLine();
                    if (line == null)
                    {
                        wr.Close();
                        rr.Close();
                        return;
                    }
                        while (line != null)
                        {
                            if(!line.Contains(">"))
                            {
                                line=line.Trim();
                                line = Regex.Replace(line, @"\s+", " ");
                                string []aux=line.Split(' ');                             
                                if (columns == null)
                                {
                                    columns = new Dictionary<int, int>[aux.Length - 2];
                                    for (int n = 0; n < columns.Length; n++)
                                        columns[n] = new Dictionary<int, int>();
                                }

                                for (int n = 0; n < aux.Length-2; n++)                                    
                                    if (!columns[n].ContainsKey(Convert.ToInt32(aux[n + 2], System.Globalization.CultureInfo.InvariantCulture)))
                                        columns[n].Add(Convert.ToInt32(aux[n + 2], System.Globalization.CultureInfo.InvariantCulture), 0);


                            }
                            line = rr.ReadLine();
                        }
                        rr.Close();
                    }
                }
            if (columns == null)
            {
                wr.Close();
                return;
            }
                double[] devStd = new double[columns.Length];
                double[,] mean = new double[columns.Length,3];
                for (int n = 0; n < columns.Length;n++ )
                    if(columns[n].Keys.Count>10)
                    {
                        double sum = 0, sum2 = 0 ;
                        mean[n, 2] = double.MinValue;
                        mean[n, 1] = double.MaxValue;
                        foreach (var item in columns[n])
                        {
                            sum += item.Key;
                            sum2 += item.Key * item.Key;
                            if (item.Key > mean[n, 2])
                                mean[n, 2] = item.Key;
                            if (item.Key < mean[n, 1])
                                mean[n, 1] = item.Key;
                        }
                        mean[n,0]=sum/columns[n].Count;                        
                        devStd[n] = Math.Sqrt(sum2 / columns[n].Count - mean[n,0] * mean[n,0]);
                    }

                    for (int i = 0; i < threadNumbers; i++)
                    {
                        string fileN = GetProfileFileName(fileName) + "_" + i;
                        using (StreamReader rr = new StreamReader(fileN))
                        {
                            if (rr == null)
                                throw new Exception("Cannot open file: " + fileN);

                            string line = rr.ReadLine();
                            while (line != null)
                            {
                                string naa = line;
                                if (!line.Contains(">"))
                                {
                                    
                                    line = line.Trim();
                                    
                                    line = Regex.Replace(line, @"\s+", " ");

                                    string[] aux = line.Split(' ');
                                  /*  if(aux.Length==zeroCount.Length)
                                        for (int n = 0; n < zeroCount.Length; n++)
                                        {
                                            if (devStd[n] > 0)
                                            {
                                                bool flag = false;
                                                double step = (mean[n, 2] - mean[n, 1]) / 10;
                                                for (int m = 1; m <= 10; m++)
                                                    if (Convert.ToInt32(aux[n + 2]) < (mean[n, 1] + step * m))
                                                    {
                                                        aux[n + 2] = m.ToString();
                                                        flag = true;
                                                        break;
                                                    }
                                                if (!flag)
                                                    aux[n + 2] = "0";
                                            }
                                            if (zeroCount[n] == count)
                                                aux[n + 2] = " ";
                                        }*/

                                    line = String.Join(" ", aux);
                                    line = line.Trim();
                                    line = Regex.Replace(line, @"\s+", " ");

                                }
                                wr.WriteLine(line);
                                line = rr.ReadLine();
                            }
                            rr.Close();
                            File.Delete(fileN);
                        }
                    }
                wr.Close();

        }
        public override List<byte> CreateNewProfile(profileNode node, string[] profile)
        {
            List<byte> newProfile = new List<byte>(profile.Length);

            for (int i = 0; i < profile.Length; i++)
            {
                string state = profile[i];
                    if (node.ContainsState(state))
                    if (Convert.ToInt32(state, System.Globalization.CultureInfo.InvariantCulture) > 255)
                        state = "255";
                if (node.ContainsState(state))
                    newProfile.Add(node.GetCodedState(node.states[state],true));
                else
                    ErrorBase.AddErrors("Unknow state " + state + " in " + node.profName + " profile!");
            }
            return newProfile;
        }
        public static Dictionary<string, protInfo<byte>> RearangeColumnOrder(Dictionary<string, protInfo<byte>> dic,string fileName)
        {
            Dictionary<string, protInfo<byte>> res = new Dictionary<string, protInfo<byte>>();
            List<string> keys = new List<string>(dic.Keys);
            StreamReader f = new StreamReader(fileName);

            string line = f.ReadLine();
            f.Close();

            string[] aux = line.Split();
            int[] index = new int[aux.Length];
            for(int i=0;i<aux.Length;i++)
            {
                index[i] = Convert.ToInt32(aux[i]);
            }

            for (int j = 0; j < keys.Count; j++)
            {
                protInfo<byte> xx = dic[keys[j]];
                List<byte> newProfile = new List<byte>();
                for (int i = 0; i < xx.profile.Count; i++)
                    newProfile.Add(xx.profile[index[i]]);

                xx.profile = newProfile;
                res.Add(keys[j], xx);
            }
            return res;
        }

        public static Dictionary<string, protInfo<byte>> RearangeColumnOrder(Dictionary<string, protInfo<byte>> dic)
        {
            Dictionary<string, protInfo<byte>> res = new Dictionary<string, protInfo<byte>>();

            List<double> sumColumn = new List<double>();

            List<string> keys = new List<string>(dic.Keys);

            for (int i = 0; i < dic[keys[0]].profile.Count; i++)
            {
                double sum = 0;
                foreach (var item in keys)
                    if (dic[item].profile[i] > 0)
                        sum += dic[item].profile[i];
                        //sum++;

                sumColumn.Add(sum);
            }
            int[] index = Enumerable.Range(0, sumColumn.Count).ToArray<int>();

            Array.Sort<int>(index, (b, a) => sumColumn[a].CompareTo(sumColumn[b]));

            for (int j = 0; j < keys.Count; j++)
            {
                protInfo<byte> xx = dic[keys[j]];
                List<byte> newProfile = new List<byte>();
                for (int i = 0; i < xx.profile.Count; i++)
                    newProfile.Add(xx.profile[index[i]]);

                xx.profile = newProfile;
                res.Add(keys[j], xx);
            }
            return res;
        }
        public override Dictionary<string, protInfo<byte>> GetProfile(profileNode node, string listFile, DCDFile dcd = null)
        {
            Dictionary<string, protInfo<byte>> res = ReadProfile(node, listFile, dcd);
            /*            StreamReader r = new StreamReader("sekw");
                        string line = r.ReadLine();
                        int []index = new int[399];
                        int counter = 0;
                        while (line != null)
                        {
                            index[counter] = Convert.ToInt32(line);
                            line = r.ReadLine();
                            counter++;
                        }*/
            //2states
            //int[] index =new int[20] { 0, 74, 57, 375, 343, 174, 33, 371, 227, 330, 9, 327, 141, 160, 256, 216, 138, 385, 51, 268 };// for proteins

            //4states
            //int[] index = new int[30] { 0, 188, 126, 4, 33, 322, 389, 84, 111, 27, 338, 45, 35, 375, 127, 235, 223,25 , 360, 232,367,28,329,121,61,258,91,225,330,301 };// for proteins

            //int[] index = new int[30] { 256, 38, 374, 54, 162, 301, 178, 363, 61, 275, 360, 302, 117, 268, 174, 40, 141, 20, 390, 51,167,338,156,28,318,356,298,10,97,158};//newr version

            /*List<string> keys = new List<string>(res.Keys);

            foreach(var item in keys)
            {
                List<byte> newProfile = new List<byte>();
                for (int i = 0; i < index.Length; i++)
                    newProfile.Add(res[item].profile[index[i]]);
                protInfo<byte> xx = res[item];
                xx.profile = newProfile;
                res[item] = xx;
            }*/
            //res = RearangeColumnOrder(res);
            // res = RearangeColumnOrder(res, "C:\\Projects\\listIndex");
            //res = ProfileStat.RearangeStates(res, 0.51);
            /*StreamWriter yy = new StreamWriter("S35_full_FragBag_org");
            foreach(var item in res)
            {
                yy.WriteLine(">" + item.Key);
                yy.Write("Selected_frag profile ");
                for(int i = 0; i < item.Value.profile.Count-1; i++)
                    yy.Write(item.Value.profile[i] + " ");

                yy.WriteLine(item.Value.profile[item.Value.profile.Count - 1]);
            }
            yy.Close();*/

            res = ProfileStat.TwoStates(res);
            //res = ProfileStat.RearangeStates(res, 0.2);
            //res = ProfileStat.SelectFeatures(res, 30);
            //ProfileStat.CalcHammingDist(res, "C:\\Projects\\tmp\\tmp\\dist_in_fragBag_CATH");
            return res;
        }
        public void ReadLibrary(string libraryName)
        {
            string[] files;
            string location = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) +
                Path.DirectorySeparatorChar + "fragLib";

            if (Directory.Exists(location))
                files = Directory.GetFiles(location);
            else
                throw new Exception("Directory fragLib not exists. Profile fragBag cannot be used!");

            fragBagLibrary.Clear();

            libraryName = location + Path.DirectorySeparatorChar + libraryName;
            if (!File.Exists(libraryName))
                throw new Exception("Cannot find: " + libraryName + " Profile fragBag cannot be used!");

            using (StreamReader s = new StreamReader(libraryName))
            {
                string line = s.ReadLine();
                float[] cent = new float[3];
                while (line != null)
                {
                    if (line.Contains("----"))
                    {
                        line = s.ReadLine();
                        List<float[]> str = new List<float[]>();
                        while (line != null && !line.Contains("***"))
                        {
                            try
                            {
                                line = line.Trim();
                                line = Regex.Replace(line, @"\s+", " ");
                                string[] aux = line.Split(' ');
                                float[] p = new float[3];
                                for (int i = 0; i < 3; i++)
                                    p[i] = (float)Convert.ToDouble(aux[i], System.Globalization.CultureInfo.InvariantCulture);

                                str.Add(p);
                            }
                            catch(Exception ex)
                            {
                                throw new Exception("Error in frag bag library: " + line);
                            }
                            line = s.ReadLine();
                        }
//                        if (str.Count > 5)
//                            Console.WriteLine("UPS");
                        float[,] auxTab = new float[str.Count, 3];
                        for (int i = 0; i < str.Count; i++)
                        {
                            for (int j = 0; j < str[i].Length; j++)
                                auxTab[i, j] = str[i][j];
                        }

                        Optimization.CenterMol(auxTab, cent);
                        fragBagLibrary.Add(auxTab);
                    }
                    line = s.ReadLine();
                }
            }


            if (fragBagLibrary.Count == 0)
                throw new Exception("Frag Bag Library is empty!");
            distData = new float[fragBagLibrary.Count];
            index = new int[fragBagLibrary.Count];
            zeroCount = new int[fragBagLibrary.Count];
            profile = new int[fragBagLibrary.Count];
        }
    }
}
