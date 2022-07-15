using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.IO;
using System.Runtime.InteropServices;
using System.Reflection;
using System.Data;

namespace uQlustCore.Profiles
{
    class RNASSProfile : InternalProfileBase<byte>
    {
        [DllImport("RNA")]
        static public extern IntPtr RNAProfiles(string fileName);
        [DllImport("RNA")]
        static public extern IntPtr GetSS(IntPtr r);
        [DllImport("RNA")]
        static public extern int GetError(IntPtr r);
        [DllImport("RNA")]
        static public extern IntPtr GetLW(IntPtr r);
        [DllImport("RNA")]
        static public extern IntPtr GetSEQ(IntPtr r);
        [DllImport("RNA")]
        static public extern void RNAFree(IntPtr r);

        static string SSprofile = "SS profile ";
        static string LWprofile = "LW profile ";
        static List<string> SSList = new List<string>()
        {
            "S","H"
        };
        static List<string> LWList = new List<string>()
        {
            ".","+",":","A","B","C","E","F","G","H","I","K","L","M","O"
        };


        public RNASSProfile()
        {
            profileName = "RNA_profiles";           
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.RNA);
            AddInternalProfiles();
            //Check if dll exists


        }
        string RNAViewErrors(int num)
        {
            switch(num)
            {
                case 1:
                    return "Unable to open file!";
                case 2:
                    return "residue name field empty";
                case 3:
                    return "xyz coordinate over %s limit. reset origin to geometrical center";
                case 4:
                    return "xyz coordinate under %s limit. reset origin to minimum xyz coordinates";
                case 5:
                    return "Error! the base is not in the library.";
                case 6:
                    return "Warning: structure with overlapped base-pairs";
                case 7:
                    return "Too many possible H-bonds between two bases";
            }
            return null;
        }
        public new double ProgressUpdate()
        {
            return (double)currentProgress / maxV;
        }

        public override void Run(DCDFile dcd)
        {
            throw new NotImplementedException();
        }
        public override bool CheckIfAvailable()
        {
            return false;
            string path;
            System.OperatingSystem osinfo = System.Environment.OSVersion;
            if (osinfo.VersionString.Contains("Win"))
            {
                if (Environment.Is64BitOperatingSystem)
                {
                    path = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) +
                        Path.DirectorySeparatorChar + "RNA.dll";
                    if (!File.Exists(path))
                        throw new Exception("Cannot find RNA.dll");
                }
                else
                    throw new Exception("For this system rna is not available");

            }
            else
            {
                path = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) +
                    Path.DirectorySeparatorChar + "libRNA.so";
                if (!File.Exists(path))
                    throw new Exception("Cannot find libRNA.so in " + path);
            }
            string dir;

            dir = Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location);

            if (Directory.Exists(dir + Path.DirectorySeparatorChar + "BASEPARS"))
                dir = dir + Path.DirectorySeparatorChar + "BASEPARS";
            else
                throw new Exception("Cannot find BASEPARS directory!");

            dir += Path.DirectorySeparatorChar;
            System.Environment.SetEnvironmentVariable("BDIR", dir);
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
            string currentDir = Directory.GetCurrentDirectory();
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
                    IntPtr pSS = IntPtr.Zero;
                    IntPtr pSEQ = IntPtr.Zero;
                    IntPtr pLW = IntPtr.Zero;
                    string SS = "";
                    string SEQ = "";
                    string LW = "";
                    int error = 0;

                    IntPtr ptr=RNAProfiles(item);
                    int errorNum = 0;
                    if (errorNum > 0)
                    {
                        ErrorBase.AddErrors(item + ": " + RNAViewErrors(error));
                        continue;
                    }
                    pSS = GetSS(ptr);
                    pLW = GetLW(ptr);
                    pSEQ = GetSEQ(ptr);
                    if (pSS == IntPtr.Zero || pLW == IntPtr.Zero || pSEQ == IntPtr.Zero)
                    {
                        ErrorBase.AddErrors("Could not create profile for: "+item+" structure");
                        RNAFree(ptr);
                        continue;
                    }
                    DebugClass.WriteMessage("PDB:" + item);
                    SS = Marshal.PtrToStringAnsi(pSS);
                    SEQ = Marshal.PtrToStringAnsi(pSEQ);
                    LW = Marshal.PtrToStringAnsi(pLW);
                    aux = item.Split(Path.DirectorySeparatorChar);
                    if (SS != null && SS.Length > 0)
                    {
                        wr.WriteLine(">" + aux[aux.Length - 1]);
                        wr.WriteLine(SSprofile + SS);
                        wr.WriteLine(LWprofile + LW);
                        wr.WriteLine(SEQprofile + SEQ);
                    }
                    RNAFree(ptr);
                    Interlocked.Increment(ref currentProgress);
                }

                wr.Close();
            }
            catch (Exception ex)
            {
                this.ex = ex;
            }
            //currentV = maxV;
            Directory.SetCurrentDirectory(currentDir);
            return 0;
        }

        public override Dictionary<string, protInfo<byte>> GetProfile(profileNode node, string fileName, DCDFile dcd = null)
        {
            Dictionary<string, protInfo<byte>> dic = new Dictionary<string, protInfo<byte>>();
            currentProgress = 0;

            if (!File.Exists(GetProfileFileName(fileName)))
                throw new Exception("Profile rnaSS does not exists!");

            string profileName;
            if (node.ContainsState("+"))
                profileName = LWprofile;
            else
                profileName = SSprofile;
            StreamReader wr = new StreamReader(GetProfileFileName(fileName));

            List<string> files = GetFileList(fileName);
            maxV = files.Count;
            Dictionary<string, int> dicFile = new Dictionary<string, int>();

            foreach (var item in files)
            {
                string[] aux = item.Split(Path.DirectorySeparatorChar);
                if (!dicFile.ContainsKey(aux[aux.Length - 1]))
                    dicFile.Add(aux[aux.Length - 1], 0);
            }

            protInfo<byte> info;
            string line = wr.ReadLine();
            string name = "";
            string seq = "";
            List<string> profile = new List<string>();
            List<byte> newProfile = new List<byte>();
            while (line != null)
            {
                if (line.Contains('>'))
                {
                    if (name.Length > 0)
                    {
                        info = new protInfo<byte>();
                        info.sequence = seq;
                        info.profile = newProfile;
                        if (dicFile.ContainsKey(name))
                            if (!dic.ContainsKey(name))
                                dic.Add(name, info);
                            else
                                ErrorBase.AddErrors("There are two the same names " + name + " in the profile file " + GetProfileFileName(fileName));
                        currentProgress++;
                    }

                    name = line.Remove(0, 1);
                }
                if (line.Contains(SEQprofile))
                    seq = line.Remove(0, SEQprofile.Length);

                if (line.Contains(profileName))
                {
                    string tmp = line.Remove(0, profileName.Length);
                    string[] aux;
                    profile.Clear();
                    if (tmp.Contains(" "))
                    {
                        tmp = tmp.TrimEnd();
                        aux = tmp.Split(' ');
                        foreach (var item in aux)
                            profile.Add(item);
                    }
                    else
                        for (int i = 0; i < tmp.Length; i++)
                            profile.Add(tmp[i].ToString());

                    newProfile = new List<byte>();
                    for (int i = 0; i < profile.Count; i++)
                        if (node.ContainsState(profile[i].ToString()))
                            newProfile.Add(node.GetCodedState(node.states[profile[i].ToString()]));
                        else
                            if (profile[i].ToString() != "-" && profile[i].ToString() != "X")
                                throw new Exception("Unknow state: " + profile[i].ToString() + " in structure " + name + " in profile " + node.profName + " profile!");
                            else
                                newProfile.Add(0);



                }

                line = wr.ReadLine();
            }

            info = new protInfo<byte>();
            info.sequence = seq;
            info.profile = newProfile;
            if (!dic.ContainsKey(name))
                dic.Add(name, info);
            else
                ErrorBase.AddErrors("There are two the same names " + name + " in the profile file " + GetProfileFileName(fileName));

            wr.Close();
            currentProgress = maxV;
            return dic;
        }



        public override void AddInternalProfiles()
        {
            profileNode node = new profileNode();

            node.profName = "RNA_SS";
            node.internalName = "RNA_SS";
            foreach (var item in SSList)
                node.AddStateItem(item, item);

            InternalProfilesManager.AddNodeToList(node, typeof(RNASSProfile).FullName);

            node = new profileNode();

            node.profName = "RNA_LW";
            node.internalName = "RNA_LW";
            foreach (var item in LWList)
                node.AddStateItem(item, item);

            InternalProfilesManager.AddNodeToList(node, typeof(RNASSProfile).FullName);
        }
        public override void RemoveInternalProfiles()
        {
            InternalProfilesManager.RemoveNodeFromList("RNA_SS");
            InternalProfilesManager.RemoveNodeFromList("RNA_LW");
        }
        public override List<INPUTMODE> GetMode()
        {
            List<INPUTMODE> oList = new List<INPUTMODE>();
            oList.Add(INPUTMODE.RNA);
            return oList;
        }
    }

}
