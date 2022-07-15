using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;
using System.Threading;
using System.IO;
using uQlustCore.PDB;
using uQlustCore.dcd;

namespace uQlustCore.Profiles
{
    public class ContactProfile:InternalProfileBase<byte>
    {
        protected string contactProfile="CONTACT profile ";
        static double kPI = 4 * Math.Atan(1.0);
        //static string PhiProfile = "PHI profile ";
        //static string PsiProfile = "PSI profile ";
        protected string ssProfile = "SS profile ";
        private short[] residuePhi;
        private short[] residuePsi;
        static Dictionary<string, string> mesoStates = new Dictionary<string, string>(){
            {"180 -165","Aa"}, {"180 -135","Ab"}, {"180 -105","Ac"}, {"180 -75","Ad"}, {"180 -45","Ae"}, {"180 -15","Af"},
            {"180 15","Ag"},  {"180 45","Ah"},   {"180 75","Ai"},   {"180 105","Aj"}, {"180 135","Ak"}, {"180 165","Al"},
            {"-150 -165","Ba"},{"-150 -135","Bb"},{"-150 -105","Bc"},{"-150 -75","Bd"},{"-150 -45","Be"},{"-150 -15","Bf"},
            {"-150 15","Bg"},  {"-150 45","Bh"},  {"-150 75","Bi"},  {"-150 105","Bj"},{"-150 135","Bk"},{"-150 165","Bl"},
            {"-120 -165","Ca"},{"-120 -135","Cb"},{"-120 -105","Cc"},{"-120 -75","Cd"},{"-120 -45","Ce"},{"-120 -15","Cf"},
            {"-120 15","Cg"},  {"-120 45","Ch"},  {"-120 75","Ci"},  {"-120 105","Cj"},{"-120 135","Ck"},{"-120 165","Cl"},
            {"-90 -165","Da"}, {"-90 -135","Db"}, {"-90 -105","Dc"}, {"-90 -75","Dd"}, {"-90 -45", "De"},{"-90 -15","Df"},
            {"-90 15","Dg"},   {"-90 45","Dh"},   {"-90 75","Di"},   {"-90 105","Dj"}, {"-90 135","Dk"}, {"-90 165","Dl"},
            {"-60 -165","Ea"}, {"-60 -135","Eb"}, {"-60 -105","Ec"}, {"-60 -75","Ed"}, {"-60 -45","Ee"}, {"-60 -15","Ef"},
            {"-60 15","Eg"},   {"-60 45","Eh"},   {"-60 75","Ei"},   {"-60 105","Ej"}, {"-60 135","Ek"}, {"-60 165","El"},
            {"-30 -165","Fa"}, {"-30 -135","Fb"}, {"-30 -105","Fc"}, {"-30 -75","Fd"}, {"-30 -45","Fe"}, {"-30 -15","Ff"},
            {"-30 15","Fg"},   {"-30 45","Fh"},   {"-30 75","Fi"},   {"-30 105","Fj"}, {"-30 135","Fk"}, {"-30 165","Fl"},
            {"0 -165","Ja"},   {"0 -135","Jb"},   {"0 -105","Jc"},   {"0 -75","Jd"},   {"0 -45","Je"},   {"0 -15","Jf"},
            {"0 15","Jg"},     {"0 45","Jh"},     {"0 75","Ji"},     {"0 105","Jj"},   {"0 135","Jk"},   {"0 165","Jl"},
            {"30 -165","Ha"},  {"30 -135","Hb"},  {"30 -105","Hc"},  {"30 -75","Hd"},  {"30 -45","He"},  {"30 -15","Hf"},
            {"30 15","Hg"},    {"30 45","Hh"},    {"30 75","Hi"},    {"30 105","Hj"},  {"30 135","Hk"},  {"30 165","Hl"},
            {"60 -165","Ia"},  {"60 -135","Ib"},  {"60 -105","Ic"},  {"60 -75","Id"},  {"60 -45","Ie"},  {"60 -15","If"},
            {"60 15","Ig"},    {"60 45","Ih"},    {"60 75","Ii"},    {"60 105","Ij"},  {"60 135","Ik"},  {"60 165","Il"},
            {"90 -165","Ja"},  {"90 -135","Jb"},  {"90 -105","Jc"},  {"90 -75","Jd"},  {"90 -45","Je"},  {"90 -15","Jf"},
            {"90 15","Jg"},    {"90 45","Jh"},    {"90 75","Ji"},    {"90 105","Jj"},  {"90 135","Jk"},  {"90 165","Jl"},
            {"120 -165","Ka"}, {"120 -135","Kb"}, {"120 -105","Kc"}, {"120 -75","Kd"}, {"120 -45","Ke"}, {"120 -15","Kf"},
            {"120 15","Kg"},   {"120 45","Kh"},   {"120 75","Ki"},   {"120 105","Kj"}, {"120 135","Kk"}, {"120 165","Kl"},
            {"150 -165","La"}, {"150 -135","Lb"}, {"150 -105","Lc"}, {"150 -75","Ld"}, {"150 -45","Le"}, {"150 -15","Lf"},
            {"150 15","Lg"},   {"150 45","Lh"},   {"150 75","Li"},   {"150 105","Lj"}, {"150 135","Lk"}, {"150 165", "Ll"}};

        static Dictionary<string, int> turns = new Dictionary<string, int>(){
             {"EfDf",1},{"EeEf",1},{"EfEf",1},{"EfDg",1},{"EeDg",1},{"EeEe",1},{"EfCg",1},{"EeDf",1},{"EkJf",1},{"EkIg",1},{"EfEe",1},{"EkJg",1},
             {"EeCg",1},{"DfDf",1},{"EfCf",1},{"DgDf",1},{"DfDg",1},{"IhIg",1},{"EfDe",1},{"EkIh",1},{"DgCg",1},{"DfCg",1},{"IbDg",1},{"DfEe",1},
             {"FeEf",1},{"IbEf",1},{"DfEf",1},{"IhJf",1},{"IhJg",1},{"IgIg",1},{"EfCh",1},{"DgEe",1},{"DgEf",1},{"EeEg",1},{"IhIh",1},{"EeDe",1},
             {"IgJg",1},{"EkKf",1},{"EeCh",1},{"IbDf",1},{"DgDg",1},{"EgDf",1},{"FeDg",1},{"ElIg",1},{"IgIh",1},{"DfDe",1},{"EjIg",1},{"EeCf",1},
             {"DfCh",1},{"DgCf",1},{"DfCf",1},{"DeEe",1},{"DkIh",1},{"FeDf",1},{"EkIf",1},{"EeDh",1},{"DgCh",1},{"IgJf",1},{"EjJg",1},{"FeEe",1},
             {"DlIh",1},{"EgCg",1},{"ElIh",1},{"EjJf",1},{"FeCg",1},{"DlIg",1},{"IbCg",1},{"EfEg",1},{"EkJe",1},{"FkJf",1},{"ElJg",1},{"DgDe",1},
             {"DlJg",1},{"EgCf",1},{"IaEf",1},{"FkIg",1},{"JaEf",1},{"EjIh",1},{"EgEf",1},{"DkJg",1},{"DeEf",1},{"EeCi",1},{"JgIh",1},{"IcEf",1},
             {"EkKe",1},{"DkIg",1},{"IbEe",1},{"EgDg",1},{"EeFe",1},{"EjKf",1},{"IaDf",1},{"HhIg",1},{"HbDg",1},{"ElJf",1},{"EfDh",1},{"IcDf",1},
             {"EfBh",1},{"IcDg",1},{"IcCg",1},{"FkJg",1},{"FeCh",1},{"IgKf",1},{"FdDg",1},{"EkHh",1},{"DfDh",1},{"DgBh",1},{"DfBh",1},{"DeDf",1},
             {"DfFe",1},{"EfFe",1},{"EgEe",1},{"EgDe",1},{"DkJf",1},{"JgJg",1},{"IbEg",1},{"IbCh",1},{"EfBg",1},{"DgCe",1},{"JlEf",1},{"CgCg",1},
             {"HhJf",1},{"EeBi",1},{"DfBi",1},{"IhIf",1},{"FeEg",1},{"FdEf",1},{"EdEf",1},{"DlJf",1},{"DhCg",1},{"JgIg",1},{"IeBg",1},{"FjIg",1},
             {"FdCh",1},{"EdEe",1},{"JfIh",1},{"JaEe",1},{"HhJg",1},{"HbEf",1},{"HbCh",1},{"FkIh",1},{"FjJf",1},{"ElJe",1},{"DhDf",1},{"CgDf",1}};


        static Dictionary<string, int> PII = new Dictionary<string, int>() { { "Dk", 1 }, { "Dl", 1 }, { "Ek", 1 }, { "El", 1 } };
        static Dictionary<string, int> helix = new Dictionary<string, int>() { { "De", 1 }, { "Df", 1 }, { "Ed", 1 }, { "Ee", 1 }, { "Ef", 1 }, { "Fd", 1 }, { "Fe", 1 } };
        static Dictionary<string, int> strand = new Dictionary<string, int>() { { "Bj", 1 }, { "Bk", 1 }, { "Bl", 1 }, { "Cj", 1 }, { "Ck", 1 }, { "Cl", 1 }, { "Dj", 1 }, { "Dk", 1 }, { "Dl", 1 } };


        protected Settings dirSettings = new Settings();
        //InternalProfilesManager manager = new InternalProfilesManager();
        
        public ContactProfile()
        {
            dirSettings.Load();
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.PROTEIN);
            profileName = "Contact";
            AddInternalProfiles();
            maxV = 1;

        }
        private string BuildKey(int phi, int psi)
        {
            string key = "";

            int delta = 30;

            if (psi > 180 || phi > 180)
                return null;

            phi = (int)Math.Round(((double)phi / delta)) * delta;
            psi =-15+ (int)(Math.Round(((double)psi+15) / delta) * delta);

            key = phi.ToString() +" "+ psi.ToString();

            return key;
        }
        public override void Run(DCDFile dcd)
        {
            DCDReader dcdReader = new DCDReader(null);
            string modelName=null;
            int counter = 0;
            string fileName;
            StreamWriter wr;
            PDBFiles pdbs = new PDBFiles();
            Dictionary<string,int> modelList=null;
            fileName = GetProfileFileName(dcd);

            if (File.Exists(fileName))
            {
                modelList = CheckFile(dcd);
            
                wr = File.AppendText(fileName);

            }
            else
                wr = new StreamWriter(fileName);

            if (wr == null)
                throw new Exception("Cannot open file: " + fileName);

            dcdReader.DCDPrepareReading(dcd.dcdFile, dcd.pdbFile);
            bool cont = true;
            do
            {
                MemoryStream mStream = new MemoryStream();
                StreamWriter write = new StreamWriter(mStream);
                cont = dcdReader.ReadAndSavePDB(write);

                modelName = "Model" + counter++;
                if (modelList == null || modelList.ContainsKey(modelName))
                    continue;

                mStream.Position = 0;
                if (this.GetType() == typeof(ContactProfile))
                    pdbs.AddPDB(mStream, PDBMODE.ALL_ATOMS, modelName);
                else
                    pdbs.AddPDB(mStream, PDBMODE.ONLY_CA, modelName);
                if(pdbs.molDic.ContainsKey(modelName))
                    MakeProfiles(modelName,pdbs.molDic[modelName], wr);

            }
            while (cont);

            dcdReader.FinishDCDReading();

            wr.Close();            
        }
        private bool CheckIfAllNeededAtomsExists(MolData molDic)
        {
            List<Residue> rList = molDic.mol.Chains[0].Residues;
            for (int i = 0; i < rList.Count; i++)
            {
                int counter = 0;
                for (int j = 0; j < rList[i].Atoms.Count; j++)
                {
                    switch (rList[i].Atoms[j].AtomName)
                    {
                        case "CA":
                        case "N":
                        case "C":
                            counter++;
                            break;
                    }
                }
                if (counter < 3)
                    return false;
            }
            return true;
        }

        public bool CalcAllDihedrals(MolData molDic)
        {
            Dictionary<string, Point3D> pfAtoms = new Dictionary<string, Point3D>() { { "CA", null }, { "C", null }, { "N", null }, { "PrevC", null }, { "NextN", null } };
            short result = 360;

            
            residuePsi[0]= result;

            if (!CheckIfAllNeededAtomsExists(molDic))
                return false;

            for (int i = 0; i < molDic.mol.Residues.Count; i++)
            {

                foreach (var item in molDic.mol.Residues[i].Atoms)
                {
                    if (pfAtoms.ContainsKey(item.AtomName))
                        pfAtoms[item.AtomName] = item.Position;
                }

                if (i > 0)
                    residuePhi[i]=(short)Math.Round(CalculateDihedralAngles(pfAtoms["PrevC"], pfAtoms["N"], pfAtoms["CA"], pfAtoms["C"]));
                else
                    residuePhi[i] = 360;
                if (i < molDic.mol.Residues.Count - 1)
                {
                    foreach (var item in molDic.mol.Residues[i + 1].Atoms)
                    {
                        if (item.AtomName == "N")
                        {
                            pfAtoms["NextN"] = item.Position;
                            break;
                        }
                    }
                    residuePsi[i] = (short)Math.Round(CalculateDihedralAngles(pfAtoms["N"], pfAtoms["CA"], pfAtoms["C"], pfAtoms["NextN"]));

                }
                else
                    residuePsi[i] = 360;

                pfAtoms["PrevC"] = pfAtoms["C"];
            }
            return true;
        }
        private double CalculateDihedralAngles(Point3D p1, Point3D p2, Point3D p3, Point3D p4)
        {
            Point3D v1 = p1 - p2;
            Point3D v2 = p4 - p3;
            Point3D v3 = p2 - p3;

            Point3D p = CrossProduct(v3, v1);
            Point3D n = CrossProduct(v3, v2);
            Point3D z = CrossProduct(v3, n);

            double u = DotProduct(n, n);
            double v = DotProduct(z, z);

            double result = 360;
            if (u > 0 && v > 0)
            {
                u = DotProduct(p, n) / Math.Sqrt(u);
                v = DotProduct(p, z) / Math.Sqrt(v);
                if (u != 0 && v != 0)
                    result = Math.Atan2(v, u) * 180 / kPI;


            }

            return result;


        }

        private Point3D CrossProduct(Point3D vec1, Point3D vec2)
        {
            Point3D resPoint = new Point3D();

            resPoint.X = vec1.Y * vec2.Z - vec1.Z * vec2.Y;
            resPoint.Y = vec1.Z * vec2.X - vec1.X * vec2.Z;
            resPoint.Z = vec1.X * vec2.Y - vec1.Y * vec2.X;

            return resPoint;
        }
        private float DotProduct(Point3D vec1, Point3D vec2)
        {
            return vec1.X * vec2.X + vec1.Y * vec2.Y + vec1.Z * vec2.Z;
        }

        protected virtual void MakeProfiles(string strName,MolData molDic, StreamWriter wr)
        {
            Dictionary<int, List<int>> contacts = new Dictionary<int, List<int>>();
            if (molDic != null)
            {
                residuePsi = new short[molDic.mol.Residues.Count];
                residuePhi = new short[molDic.mol.Residues.Count];

                if (!CalcAllDihedrals(molDic))
                    return;
                molDic.CreateSideChainContactMap(8.5f, false);
                //pdbs.molDic[strName].CreateContactMap(9.5f, "CB");


                foreach (var contItem in molDic.contactMap.Keys)
                {
                    if (!contacts.ContainsKey(contItem))
                        contacts.Add(contItem, new List<int>());


                    foreach (var itemList in molDic.contactMap[contItem])
                    {
                        contacts[contItem].Add((int)itemList);
                        if (!contacts.ContainsKey((int)itemList))
                        {
                            contacts.Add((int)itemList, new List<int>());
                            contacts[(int)itemList].Add(contItem);
                        }
                        else
                            if (!contacts[(int)itemList].Contains(contItem))
                                contacts[(int)itemList].Add(contItem);

                    }
                }

                int num;
                string profile = "";
                int len = molDic.mol.Chains[0].chainSequence.Length;
                for (int i = 0; i < len; i++)
                {

                    if (contacts.ContainsKey(i))
                    {
                        num = contacts[i].Count;
                        if (num > 9)
                            num = 9;
                    }
                    else
                        num = 0;
                    profile += num;
                    if (i < len - 1)
                        profile += " ";
                }
                string psiProfile = "";
                string phiProfile = "";
                string ss = "";
                List<string> meso = new List<string>();
                List<string> mesoDF = new List<string>();
                for (int i = 0; i < molDic.mol.Residues.Count; i++)
                {
                    string key;
                    psiProfile += residuePsi[i];
                    phiProfile += residuePhi[i];
                    key = BuildKey(residuePhi[i], residuePsi[i]);
                    meso.Add(key);
                    if (i < molDic.mol.Residues.Count - 1)
                    {
                        ss += " ";
                        psiProfile += " ";
                        phiProfile += " ";
                    }

                }
                char[] ssStates = new char[meso.Count];
                for (int i = 0; i < meso.Count; i++)
                    ssStates[i] = 'C';
                for (int i = 0; i < meso.Count; i++)
                {
                    if (meso[i] == null)
                        continue;


                    if (mesoStates.ContainsKey(meso[i]))
                    {
                        //mesoDF.Add(mesoStates[meso[i]]);
                        if (PII.ContainsKey(mesoStates[meso[i]]))
                            ssStates[i] = 'P';

                        if (i < meso.Count - 1)
                        {
                            if (meso[i + 1] != null && mesoStates.ContainsKey(meso[i + 1]))
                            {
                                string combStates = mesoStates[meso[i]] + mesoStates[meso[i + 1]];
                                if (turns.ContainsKey(combStates))
                                {
                                    ssStates[i] = ssStates[i + 1] = 'T';
                                }
                            }
                        }
                    }
                }
                for (int i = 0; i < meso.Count; i++)
                {
                    int counter = 0;
                    if (meso[i] == null)
                        continue;
                    if (mesoStates.ContainsKey(meso[i]))
                    {
                        if (helix.ContainsKey(mesoStates[meso[i]]) && i < meso.Count - 5)
                        {

                            for (int k = 0; k < 5; k++)
                                if (meso[i + k] != null && mesoStates.ContainsKey(meso[i + k]) && helix.ContainsKey(mesoStates[meso[i + k]]))
                                    counter++;
                                else
                                    break;
                            if (counter == 5)
                                for (int k = 0; k < 5; k++)
                                    ssStates[i + k] = 'H';
                        }
                        counter = 0;
                        if (i < meso.Count - 3)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                if (meso[i + k] != null && mesoStates.ContainsKey(meso[i + k]) && strand.ContainsKey(mesoStates[meso[i + k]]))
                                    counter++;
                                else
                                    break;
                            }
                            if (counter == 3)
                                for (int k = 0; k < 3; k++)
                                    if (ssStates[i + k] == 'C' || ssStates[i + k] == 'P')
                                        ssStates[i + k] = 'E';
                        }
                    }

                }

                ss = new string(ssStates);

                if (profile.Length > 0)
                {
                    wr.WriteLine(">" + strName);
                    wr.WriteLine(contactProfile + profile);
                    //wr.WriteLine(PsiProfile+psiProfile);
                    //wr.WriteLine(PhiProfile+phiProfile);
                    wr.WriteLine(ssProfile + ss);
                    wr.WriteLine(SEQprofile + molDic.mol.Chains[0].chainSequence);

                }

                molDic.CleanMolData();
            }

        }
        public override int Run(object processParams)
        {
            Dictionary<int, List<int>> contacts= new Dictionary<int, List<int>>();
            string fileName = ((ThreadFiles)(processParams)).fileName;

            List<string> files = threadingList[((ThreadFiles)(processParams)).threadNumber];// CheckFile(listFile);
            if (files.Count == 0)
            {
                return 0;
            }
            else Console.WriteLine("Jest: ", files.Count);
            StreamWriter wr;
                                 
            if (File.Exists(fileName))
                wr = File.AppendText(fileName);
            else
                wr = new StreamWriter(fileName);

            if (wr == null)
                throw new Exception("Cannot open file: " + fileName);

            PDBFiles pdbs = new PDBFiles();

            string refSeqFile = ((ThreadFiles)(processParams)) + ".ref";
            pdbs.ReadRefSeq(refSeqFile);
            //maxV = files.Count;
            foreach (var item in files)
            {
                try
                {

                    string strName;
                    //strName = pdbs.AddPDB(item, PDBMODE.CA_CB);
                    if (dirSettings.mode == INPUTMODE.PROTEIN)
                        if (this.GetType() == typeof(ContactProfile))
                            strName = pdbs.AddPDB(item, PDBMODE.ALL_ATOMS);
                        else
                            strName = pdbs.AddPDB(item, PDBMODE.ONLY_CA);
                    else
                        strName = pdbs.AddPDB(item, PDBMODE.ONLY_P_AND_C4);
                    if (strName != null)
                        MakeProfiles(strName, pdbs.molDic[strName], wr);

                    Interlocked.Increment(ref currentProgress);
                }
                catch (Exception ex)
                {
                    ErrorBase.AddErrors("Incorrect profile for " + item);
                }

            }

            wr.Close();
            currentProgress = maxV;
            return 0;
        }
        public virtual List<byte> CreateNewProfile(profileNode node,string [] profile)
        {
            List <byte>newProfile = new List<byte>(profile.Length);

            for (int i = 0; i < profile.Length; i++)
            {
                string state = profile[i];
                if (node.ContainsState(state))
                {
                   /* byte st= node.GetCodedState(node.states[state]);

                    for (int j = 1; j <= 10;j++) // node.states.Count; j++)
                        if (st == j)
                            newProfile.Add(2);
                        else
                            newProfile.Add(1);*/

                    newProfile.Add(node.GetCodedState(node.states[state]));
                }
                else
                    ErrorBase.AddErrors("Unknow state " + state + " in " + node.profName + " profile!");
            }
            return newProfile;
        }
        protected Dictionary<string, protInfo<byte>> ReadProfile(profileNode node, string listFile, DCDFile dcd = null)
        {
            Dictionary<string, protInfo<byte>> dic = null;
            StreamReader wr;
            DebugClass.WriteMessage("start reading file" + listFile);
            string useFile = listFile;


            List<string> names = GetFileList(listFile);
            if (names == null)
                DebugClass.WriteMessage("After GetFileList null");
            else
                DebugClass.WriteMessage("After GetFileList " + names.Count);
            maxV = names.Count;
            Dictionary<string, int> dicNames = new Dictionary<string, int>(names.Count);

            foreach (var item in names)
            {
                string[] aux = item.Split(Path.DirectorySeparatorChar);
                if (!dicNames.ContainsKey(aux[aux.Length - 1]))
                    dicNames.Add(aux[aux.Length - 1], 0);
            }
            DebugClass.WriteMessage("After loop");
            string[] dotFlag = names[0].Split('|');
            bool flag = false;
            if (dotFlag.Length == 2)
                flag = true;
            DebugClass.WriteMessage("Before profile");
            if (dcd == null)
            {
                wr = new StreamReader(GetProfileFileName(listFile));
            }
            else
                wr = new StreamReader(GetProfileFileName(dcd));

            DebugClass.WriteMessage("start ddd reading");

            protInfo<byte> info;
            string line = wr.ReadLine();
            string name = "";
            string seq = "";
            List<string> profile; ;
            List<byte> newProfile = null;

            string profileName;
            if (node.ContainsState("H"))
                profileName = ssProfile;
            else
                profileName = contactProfile;
            //Check number of elements in file
            int lineLength = 0;


            while (line != null)
            {
                if (lineLength < line.Length)
                    lineLength = line.Length;
                line = wr.ReadLine();

            }
            wr.BaseStream.Position = 0;
            wr.DiscardBufferedData();
            profile = new List<string>(lineLength);
            dic = new Dictionary<string, protInfo<byte>>(names.Count);
            line = wr.ReadLine();
            string remName = "";
            DebugClass.WriteMessage("Starrrrr");
            while (line != null)
            {
                if (line[0] == '>')
                {
                    if (name.Length > 0)
                    {
                        if (dicNames.ContainsKey(name))
                        {
                            info = new protInfo<byte>();
                            info.sequence = seq;
                            info.profile = newProfile;
                            //string newName = line;
                            //newName = newName.Remove(0, 1);                            
                            if (!dic.ContainsKey(remName))
                                dic.Add(remName, info);
                            else
                                ErrorBase.AddErrors("!!In the generated profile file structure " + remName + " exists more then once!!!\nOnly the first profile will be considered!");
                            currentProgress++;
                        }
                    }

                    name = line.Remove(0, 1);
                    remName = name;
                    if (name.Contains("|") && !flag)
                    {
                        string[] tmp = name.Split('|');
                        name = tmp[0];

                    }
                }
                if (dicNames.ContainsKey(name))
                {
                    if (line.Contains(SEQprofile))
                        seq = line.Remove(0, SEQprofile.Length);
                    else
                        if (line.Contains(profileName))
                    {
                        //StringBuilder tmp = new StringBuilder(line);
                        line = line.Remove(0, profileName.Length);
                        //tmp = tmp.Remove(0, profileName.Length);
                        //StringBuilder tmp = new StringBuilder (line.Remove(0, profileName.Length));
                        string[] aux;
                        if (line.Contains(' '))
                            aux = line.Split(' ');
                        else
                        {
                            aux = new string[line.Length];
                            for (int i = 0; i < line.Length; i++)
                                aux[i] = line[i].ToString();
                        }
                        newProfile = CreateNewProfile(node, aux);

                    }
                }
                line = wr.ReadLine();
            }
            if (dicNames.ContainsKey(name))
            {
                info = new protInfo<byte>();
                info.sequence = seq;
                info.profile = newProfile;
                dic.Add(remName, info);
            }

            wr.Close();
            DebugClass.WriteMessage("Reading finished");
            return dic;


        }
        public override Dictionary<string, protInfo<byte>> GetProfile(profileNode node, string listFile,DCDFile dcd=null)
        {
            return ReadProfile(node, listFile, dcd);
        }
        public override void AddInternalProfiles()
        {
            profileNode node = new profileNode();

            node.profName = "Contact";
            node.internalName = "Contact";
            for (int i = 0; i < 10; i++)
                node.AddStateItem(i.ToString(), i.ToString());

            //manager.AddNodeToList(node, "P");
            node.AddStateItem("T", "T");
            node.AddStateItem("C", "C");
            node.AddStateItem("P", "P");
            node.AddStateItem("H", "H");
            node.AddStateItem("E", "E");

             InternalProfilesManager.AddNodeToList(node, typeof(ContactProfile).FullName);

        }
        public override void RemoveInternalProfiles()
        {
            InternalProfilesManager.RemoveNodeFromList("Contact");
            InternalProfilesManager.RemoveNodeFromList("SS Dihedral");

        }
        public override List<INPUTMODE> GetMode()
        {
            List<INPUTMODE> oList = new List<INPUTMODE>();
            oList.Add(INPUTMODE.PROTEIN);

            return oList;
        }
    }
}