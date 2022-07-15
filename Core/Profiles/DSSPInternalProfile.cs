using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;
using System.Runtime.InteropServices;
using System.Reflection;
using System.Data;

namespace uQlustCore.Profiles
{
    public class DsspInternalProfile:InternalProfileBase<byte>
    {

        [DllImport("DSSP")]
        static public extern IntPtr PrepProtein();
        [DllImport("DSSP")]
            static public extern int ReadProt(string name, IntPtr obj);
        [DllImport("DSSP")]
        static public extern IntPtr DisposeMProtein(IntPtr obj);
        [DllImport("DSSP")]
        static public extern IntPtr DisposeBuffor(IntPtr obj);
        [DllImport("DSSP")]
        static public extern IntPtr GetSS(IntPtr obj);
        [DllImport("DSSP")]
        static public extern IntPtr GetSA(IntPtr obj);
        [DllImport("DSSP")]
        static public extern IntPtr GetSEQ(IntPtr obj);

        protected string SSprofile = "SS profile ";
        protected string SAprofile = "SA profile ";
        
        //static string PSIprofile = "PSI profile ";
        //static string PHIprofile = "PHI profile ";

        static Dictionary<string, int> surface = new Dictionary<string, int>()
        {
            {"A",115},{"R",225},{"D",150},{"N",160},{"C",135},
            {"E",190},{"Q",180},{"G",75},{"H",195},{"I",175},
            {"L",170},{"K",200},{"M",185},{"F",210},{"P",145},
            {"S",115},{"T",140},{"W",255},{"Y",230},{"V",155},{"X",1}
        };
        static List<string> SSList = new List<string>()
        {
            "H","E","I","S","B","T","C","G"
        };
        //InternalProfilesManager manager = new InternalProfilesManager();
    
        public DsspInternalProfile()
        {
            profileName = "dssProfile.native";
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.PROTEIN);
            AddInternalProfiles();
            //Check if dll exists


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
                    path=Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) +
                        Path.DirectorySeparatorChar + "DSSP.dll";
                    if (!File.Exists(path))
                        throw new Exception("Cannot find DSSP.dll");
                }
                else
                    throw new Exception("For this system dssp is not available");

            }
            else
            {
                path=Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) +
                    Path.DirectorySeparatorChar + "libDSSP.so";
                if (!File.Exists(path))
                        throw new Exception("Cannot find libDSSP.so in "+path);
            }

            

        }
        protected string DSSPErrors(int num)
        {
            switch (num)
            {
                case 1:
                    return "DSSP: inconsistent residue types in atom records for residue";                    
                case 2:
                    return "DSSP: inconsistent residue sequence numbers";                    
                case 3:
                    return "DSSP: first and second residue are the same";                    
                case 4:
                    return "DSSP: empty protein, or no valid complete residues";
                case 5:
                    return "DSSP: Only cysteine residues can form sulphur bridges";
                case 6:
                    return "DSSP: Some other error";

            }
            return "DSSP: Some other error";
        }

        public override  int Run(object processParams)
        {
            string[] aux;

            System.OperatingSystem osinfo=System.Environment.OSVersion;

            if (osinfo.VersionString.Contains("Win"))
            {
                if (!Environment.Is64BitOperatingSystem)
                {
                    this.ex = new Exception("Dssp profile is not working on 32 bit windows system!");
                    return 0;
                }
            }

            StreamWriter wr;
            wr = new StreamWriter((( ThreadFiles)processParams).fileName);

            if (wr == null)
            {
                this.ex= new Exception("Cannot open file: " + ((ThreadFiles)processParams).fileName);
                return 0;
            }
            
            List<string> auxFiles = threadingList[(( ThreadFiles)processParams).threadNumber];

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
                    IntPtr pSA = IntPtr.Zero;
                    IntPtr pSEQ = IntPtr.Zero;
                    string SS = "";
                    string SA = "";
                    string SEQ = "";
                    int error = 0;


                    dsspExt = PrepProtein();
                    DebugClass.WriteMessage("PDB:" + item);
                    error = ReadProt(item, dsspExt);
                    if (error == 0)
                    {
                        pSS = GetSS(dsspExt);
                        pSA = GetSA(dsspExt);
                        pSEQ = GetSEQ(dsspExt);
                        SS = Marshal.PtrToStringAnsi(pSS);
                        SA = Marshal.PtrToStringAnsi(pSA);
                        SEQ = Marshal.PtrToStringAnsi(pSEQ);
                        aux = item.Split(Path.DirectorySeparatorChar);
                        if (SS != null && SS.Length > 0)
                        {
                            wr.WriteLine(">" + aux[aux.Length - 1]);
                            wr.WriteLine(SSprofile + SS);
                            wr.WriteLine(SAprofile + ConvertSAProfile(SA, SEQ));
                            //wr.WriteLine(PSIprofile + wrapper.PSI);
                            //wr.WriteLine(PHIprofile + wrapper.PHI);
                            wr.WriteLine(SEQprofile + SEQ);

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
                    if (pSA != IntPtr.Zero)
                        DisposeBuffor(pSA);
                    if (pSEQ != IntPtr.Zero)
                        DisposeBuffor(pSEQ);

                    Interlocked.Increment(ref currentProgress);                    
                }

                wr.Close();
            }
            catch(Exception ex)
            {
                this.ex = ex;               
            }
            //currentV = maxV;
            return 0;
        }
        private string ConvertSAProfile(string saProfile,string seq)
        {
            string convert="";
            string[] aux;

            saProfile = saProfile.Remove(saProfile.Length - 1, 1);
            aux = saProfile.Split(' ');
            int i=0;
            foreach(var item in aux)
            {
                int con = Convert.ToInt32(item);

                double surf = ((double)con) / surface[seq[i].ToString()];

                if (surf > 1)
                    surf = 1;

                if (seq[i] == 'X')
                    surf = -1;

                surf = (int)((surf+0.05)*10);

                if (surf > 9)
                    surf = 9;
                if (surf >= 0)
                    convert += surf.ToString();
                else
                    convert += "X";
                i++;
                if (i < seq.Length)
                    convert += " ";

                
            }
            return convert;
        }
        public override Dictionary<string, protInfo<byte>> GetProfile(profileNode node, string fileName,DCDFile dcd=null)
        {
            Dictionary<string,protInfo<byte>> dic=new Dictionary<string,protInfo<byte>>();
            currentProgress = 0;

            if (!File.Exists(GetProfileFileName(fileName)))
                throw new Exception("Profile dssp does not exists!");

            string profileName;
            if (node.ContainsState("H"))
                profileName = SSprofile;
            else
                profileName = SAprofile;
            StreamReader wr = new StreamReader(GetProfileFileName(fileName));
            
            List<string> files = GetFileList(fileName);
            maxV = files.Count; 
            Dictionary<string, int> dicFile = new Dictionary<string, int>();

            foreach (var item in files)
            {
                string[] aux = item.Split(Path.DirectorySeparatorChar);
                if (!dicFile.ContainsKey(aux[aux.Length - 1]))
                    dicFile.Add(aux[aux.Length-1],0);
            }

            protInfo<byte> info;
            string line = wr.ReadLine();
            string name="";
            string seq="";
            List <string> profile =new List<string>();
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
                        if( dicFile.ContainsKey(name))
                            if (!dic.ContainsKey(name))
                                dic.Add(name, info);
                            else
                                ErrorBase.AddErrors("There are two the same names " + name + " in the profile file " + GetProfileFileName(fileName));
                        currentProgress++;
                    }

                    name = line.Remove(0, 1);
                }
                if(line.Contains(SEQprofile))                
                    seq = line.Remove(0, SEQprofile.Length);

                if (line.Contains(profileName))
                {
                    string tmp = line.Remove(0, profileName.Length);
                    string[] aux;
                    profile.Clear();
                    if (tmp.Contains(" "))
                    {
                        tmp=tmp.TrimEnd();
                        aux = tmp.Split(' ');
                        foreach(var item in aux)
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
                                throw new Exception("Unknow state: " + profile[i].ToString() + " in structure " +name+" in profile "+ node.profName + " profile!");
                            else
                                newProfile.Add(0);



                }

                line = wr.ReadLine();
            }

            info = new protInfo<byte>();
            info.sequence = seq;
            info.profile = newProfile;
            if(!dic.ContainsKey(name))
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

            node.profName = "SS";
            node.internalName = "SS";
            foreach(var item in SSList)
                node.AddStateItem(item,item);

             InternalProfilesManager.AddNodeToList(node, typeof(DsspInternalProfile).FullName);

            node = new profileNode();

            node.profName = "SA";
            node.internalName = "SA";
            for(int i=0;i<10;i++)
                node.AddStateItem(i.ToString(), i.ToString());

             InternalProfilesManager.AddNodeToList(node, typeof(DsspInternalProfile).FullName);
        }
        public override void RemoveInternalProfiles()
        {
            InternalProfilesManager.RemoveNodeFromList("SA");
            InternalProfilesManager.RemoveNodeFromList("SS");
        }
        public override List<INPUTMODE> GetMode()
        {
            List<INPUTMODE> oList = new List<INPUTMODE>();
            oList.Add(INPUTMODE.PROTEIN);

            return oList;
        }
    }
}
