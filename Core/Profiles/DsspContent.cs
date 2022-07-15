using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.IO;
using System.Reflection;
using System.Runtime.InteropServices;

namespace uQlustCore.Profiles
{
    class DsspContent: DsspInternalProfile
    {
        
        protected static List<char> states = new List<char>() { 'H', 'E','C'};
        public DsspContent()
        {
            SSprofile = "DSSPContent ";
            SAprofile = "DSSPContent ";

            profileName = "dssContent.native";
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.PROTEIN);
            AddInternalProfiles();
            //Check if dll exists


        }

        public override int Run(object processParams)
        {
            string[] aux;

            System.OperatingSystem osinfo = System.Environment.OSVersion;

            if (osinfo.VersionString.Contains("Win"))
            {
                if (!Environment.Is64BitOperatingSystem)
                {
                    this.ex = new Exception("Dssp profile is not working on 32 bit windows system!");
                    return 0;
                }
            }

            StreamWriter wr;
            wr = new StreamWriter(((ThreadFiles)processParams).fileName);

            if (wr == null)
            {
                this.ex = new Exception("Cannot open file: " + ((ThreadFiles)processParams).fileName);
                return 0;
            }

            List<string> auxFiles = threadingList[((ThreadFiles)processParams).threadNumber];

            try
            {

                foreach (var item in auxFiles)
                {
                    //  wrapper.timeSp = 0;

                    if (!File.Exists(item))
                    {
                        ErrorBase.AddErrors("File " + item + "does not exist");
                        continue;
                    }

                    //wrapper.Run(item,item.Length);
                    //timeSp += wrapper.timeSp;
                    IntPtr dsspExt = IntPtr.Zero;
                    IntPtr pSS = IntPtr.Zero;
                    string SS = "";
                    int error = 0;


                    dsspExt = PrepProtein();
                    DebugClass.WriteMessage("PDB:" + item);
                    error = ReadProt(item, dsspExt);
                    if (error == 0)
                    {
                        pSS = GetSS(dsspExt);
                        SS = Marshal.PtrToStringAnsi(pSS);
                        SS=SS.Replace('B', 'E');
                        SS = SS.Replace('G', 'H');
                        SS = SS.Replace('I', 'C');
                        SS = SS.Replace('S', 'C');
                        SS = SS.Replace('T', 'C');

                        Dictionary<char, int> counter = new Dictionary<char, int>();
                        for (int i = 0; i < SS.Length; i++)
                        {
                            if (counter.ContainsKey(SS[i]))
                                counter[SS[i]]++;
                            else
                                counter.Add(SS[i], 1);
                        }
                        List<char> keyList = new List<char>(counter.Keys);
                        foreach (var it in keyList)
                        {
                            double res = ((double)counter[it]) / SS.Length;
                            res *= 10;
                            counter[it] = (int)Math.Floor(res);
                        }
                       
                        aux = item.Split(Path.DirectorySeparatorChar);
                        if (SS != null && SS.Length > 0)
                        {
                            wr.WriteLine(">" + aux[aux.Length - 1]);
                            string txt = "";
                            for (int i = 0; i < states.Count - 1; i++)
                            {
                                if (counter.ContainsKey(states[i]))
                                    txt += counter[states[i]] + " ";
                                else
                                    txt += "0 ";
                            }
                            if (counter.ContainsKey(states[states.Count - 1]))
                                txt += counter[states[states.Count - 1]];
                            else
                                txt += "0";
                            wr.WriteLine(SSprofile + txt);
                        }
                    }
                    else
                    {
                        ErrorBase.AddErrors(item + ": " + DSSPErrors(error));
                    }
                    if (dsspExt != IntPtr.Zero)
                        DisposeMProtein(dsspExt);
                    if (pSS != IntPtr.Zero)
                        DisposeBuffor(pSS);

                    Interlocked.Increment(ref currentProgress);
                }

                wr.Close();
            }
            catch (Exception ex)
            {
                this.ex = ex;
            }
            //currentV = maxV;
            return 0;
        }

        public override void AddInternalProfiles()
        {
            profileNode node = new profileNode();

            node.profName = "dsspContent";
            node.internalName = "dsspContent";
            for (int i = 0; i <=10; i++)
                node.AddStateItem(i.ToString(), i.ToString());

            InternalProfilesManager.AddNodeToList(node, typeof(DsspContent).FullName);

        }
        public override void RemoveInternalProfiles()
        {
            InternalProfilesManager.RemoveNodeFromList("SA");
            InternalProfilesManager.RemoveNodeFromList("SS");
        }

    }
}
