using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using uQlustCore;
using System.Collections.Concurrent;
using System.Threading;
using System.Diagnostics;

namespace uQlustCore.PDB
{


    public enum PDBMODE { ONLY_CA, ONLY_CB, ALL_ATOMS, CA_CB,ONLY_SEQ, ONLY_P,ONLY_P_AND_C4 };
    public enum CHAIN_MODE { SINGLE, ALL };
    public class IncorrectSideChainException : Exception
    {
        public IncorrectSideChainException(String s):base(s)
        {            
        }
    }
    public class Params
    {
        public string fileName;
        public int k;
    }

    public class MolData
    {
        static double kPI = 4 * Math.Atan(1.0);
        public Molecule mol=null;
        public int molLength;
        public int[] indexMol;
        public Dictionary <int,List<float>> contactMap=null;
        public int[] fullContactMap = null;

        private static float[] center= new float[3];
        

        //Flag true only CA will be readed
        public MolData(string fileName,PDBMODE flag,INPUTMODE mode,CHAIN_MODE chainFlag=CHAIN_MODE.SINGLE)
        {
            using (StreamReader rr = new StreamReader(fileName))            {
                switch (mode)
                {
                    case INPUTMODE.PROTEIN:
                        mol = new Molecule( flag,chainFlag);
                        break;
                    case INPUTMODE.RNA:
                        mol = new MoleculeRNA(flag);
                        break;
                }
                
                bool res = mol.ReadMolecule(rr);
                if (!res)
                    return;
                if (mol.Chains.Count == 0)
                {
                    ErrorBase.AddErrors("Error in reading file: " + fileName + "\nCannot find residues, file will not be considered!");
                    rr.Close();
                    return;
                }
               molLength = mol.Chains[0].chainSequence.Length;         
            }
            //CenterMol();
        }
        public MolData()
        {
            
        }
        public bool ReadMolData(MemoryStream stream, PDBMODE flag,string modelName)
        {
            StreamReader rr = new StreamReader(stream);
            mol = new Molecule(flag);
            bool res = mol.ReadMolecule(rr);
            if (!res)
            {
                return false;
            }

            if (mol.Chains.Count == 0)
            {
                ErrorBase.AddErrors("Error in reading file: " + modelName + "\nCannot find residues, file will not be considered!");
                return false;
            }
            molLength = mol.Chains[0].chainSequence.Length;

            return true;
            //CenterMol();
        }
        public void CleanMolData()
        {
            mol = null;
        }
        private Point3D GetAtomPosition(int indexM,string atomName)
        {
            for (int n = 0; n < mol.Residues[indexM].Atoms.Count; n++)
                if (mol.Residues[indexM].Atoms[n].AtomName == atomName)
                    return mol.Residues[indexM].Atoms[n].Position;

            return null;
        }
        public void CreateFullContactMap(float threshold,byte [] tab,string atomName)
        {
            float distance;
            int nextStart = 12;
            //int nextStart = 6;
            threshold *= threshold;
            int counter = 0;
            int thr1 = 4 * 4;
            int thr2 = 6 * 6;
            int thr3 = 8 * 8;
            for (int j = 0; j < indexMol.Length; j++)
            {
                if (indexMol[j] != -1)
                {
                    Point3D atom1 = GetAtomPosition(indexMol[j],atomName) ;

                   
                    for (int i = j + nextStart; i < indexMol.Length; i++)
                    {
                        if (i == j + 3)
                            continue;
                        distance = 0;

                        if (indexMol[i] != -1)
                        {
                            Point3D atom2 = GetAtomPosition(indexMol[i],atomName);
                            if (atom1 != null && atom2 != null)
                            {
                                float tmp = atom1.X - atom2.X;
                                distance += tmp * tmp;
                                tmp = atom1.Y - atom2.Y;
                                distance += tmp * tmp;
                                tmp = atom1.Z - atom2.Z;
                                distance += tmp * tmp;

                                if (i == j + 2)
                                {
                                    if (distance < thr1)
                                        tab[counter++] = 1;
                                    else
                                        if (distance < thr2)
                                            tab[counter++] = 2;
                                        else
                                            if (distance < thr3)
                                                tab[counter++] = 3;
                                            else
                                                tab[counter++] = 4;


                                }
                            }
                            else
                                distance = float.MaxValue;

                            if (distance < threshold)
                                tab[counter++] = 1;
                            else
                                tab[counter++] = 0;
                        }
                        else
                            tab[counter++] = 0;
                    }
                }
                else
                {
                    for (int i = j + nextStart; i < indexMol.Length; i++)
                    {
                        if (i == j + 3)
                            continue;
                        tab[counter++] = 0;
                    }
                }
            }
            


        }
        Point3D AvrAtomsPosition(List<Atom> atoms)
        {
            Point3D pos = new Point3D();
            int count=0;
            float x=0, y=0, z = 0;
            for (int i = 0; i < atoms.Count; i++)
            {
                if (atoms[i].AtomName == "C" || atoms[i].AtomName == "CA" || atoms[i].AtomName == "N" || atoms[i].AtomName == "O")
                    continue;                
                x += atoms[i].Position.X;
                y += atoms[i].Position.Y;
                z += atoms[i].Position.Z;
                count++;
            }
            x /= count;
            y /= count;
            z /= count;


            pos.X=x;
            pos.Y=y;
            pos.Z=z;

            return pos;

        }
        /*private bool CheckIfAllNeededAtomsExists()
        {
            List<Residue> rList = mol.Chains[0].Residues;
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
        public bool CalcAllDihedrals()
        {
            Dictionary <string,Point3D> pfAtoms=new Dictionary<string,Point3D> () {{"CA",null},{"C",null},{"N",null},{"PrevC",null},{"NextN",null}};
            short result = 360;
           
            
            
            mol.Residues[0].Psi=result;

            if (!CheckIfAllNeededAtomsExists())
                return false;

            for(int i=0;i<mol.Residues.Count;i++)
            {
                
                foreach(var item in mol.Residues[i].Atoms)
                {
                    if(pfAtoms.ContainsKey(item.AtomName))
                        pfAtoms[item.AtomName]=item.Position;
                }
                
                if(i>0)
                    mol.Residues[i].Phi=(short)Math.Round(CalculateDihedralAngles(pfAtoms["PrevC"],pfAtoms["N"],pfAtoms["CA"],pfAtoms["C"]));
                else
                    mol.Residues[i].Phi=360;
                if (i < mol.Residues.Count - 1)
                {
                    foreach (var item in mol.Residues[i+1].Atoms)
                    {
                        if (item.AtomName == "N")
                        {
                            pfAtoms["NextN"] = item.Position;
                            break;
                        }
                    }
                    mol.Residues[i].Psi = (short)Math.Round(CalculateDihedralAngles(pfAtoms["N"], pfAtoms["CA"], pfAtoms["C"], pfAtoms["NextN"]));

                }
                else
                    mol.Residues[i].Psi = 360;

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

            resPoint.X=vec1.Y * vec2.Z - vec1.Z * vec2.Y;
            resPoint.Y=vec1.Z * vec2.X - vec1.X * vec2.Z;
            resPoint.Z=vec1.X * vec2.Y - vec1.Y * vec2.X;

            return resPoint;
        }
        private float DotProduct(Point3D vec1, Point3D vec2)
        {
            return vec1.X * vec2.X + vec1.Y * vec2.Y + vec1.Z * vec2.Z;
        }*/
        public void CreateContactMap(float threshold, bool flagDistance,string atomName)
        {
            float distance;
            Point3D [] sideResPosition;

            sideResPosition = new Point3D [mol.Residues.Count];
            float localThreshold = threshold * threshold;

            
                
            for (int i = 0; i < sideResPosition.Length; i++)
            {
                if (atomName!=null)
                    sideResPosition[i] = GetAtomPosition(i,atomName);
                else
                    sideResPosition[i] = AvrAtomsPosition(mol.Residues[i].Atoms);
            }
            contactMap = new Dictionary<int, List<float>>();//(mol.Residues.Count*(mol.Residues.Count+1)/2);

            for (int j = 0; j < mol.Residues.Count; j++)
            {
                List<float> rijList = new List<float>();
                for (int i = j + 6; i < mol.Residues.Count; i++)
                {                                        
                    distance = 0;
                    float tmp = sideResPosition[j].X-sideResPosition[i].X;
                    distance += tmp * tmp;
                    tmp =sideResPosition[j].Y-sideResPosition[i].Y;
                    distance += tmp * tmp;
                    tmp = sideResPosition[j].Z-sideResPosition[i].Z;
                    distance += tmp * tmp;
                    if (distance < localThreshold)
                        if(flagDistance)
                            rijList.Add(distance);
                        else
                            rijList.Add(i);
                }
                if (rijList.Count > 0)
                    contactMap.Add(j, rijList);
            }
        }
        public void CreateCAContactMap(float threshold,bool flagDistance)
        {
            CreateContactMap(threshold, flagDistance, "CA");
        }
        //if flag distance is true creates conatct map with distances otherwise 
        //indexes of resiudes are stored
        public void CreateSideChainContactMap(float threshold,bool flagDistance)
        {
            CreateContactMap(threshold, flagDistance, null);
        }

        private void CenterMol()
        {
            center[0] = 0; center[1] = 0; center[2] = 0;
            for (int i = 0; i < mol.Atoms.Count; i++)
            {
                center[0] += mol.Atoms[i].Position.X;
                center[1] += mol.Atoms[i].Position.Y;
                center[2] += mol.Atoms[i].Position.Z;
            }
            center[0]/=mol.Atoms.Count;
            center[1]/=mol.Atoms.Count;
            center[2]/=mol.Atoms.Count;
            for (int i = 0; i < mol.Atoms.Count; i++)
            {
                mol.Atoms[i].Position.X=mol.Atoms[i].Position.X - center[0];
                mol.Atoms[i].Position.Y=mol.Atoms[i].Position.Y - center[1];
                mol.Atoms[i].Position.Z=mol.Atoms[i].Position.Z - center[2];
            }

        }
		
    }

    public class PDBFiles
    {
        public Dictionary<string, MolData> molDic = new Dictionary<string, MolData>();
        string refStuctName;
        public string refSeq = null;
        ManualResetEvent[] resetEvents;
        List<string>[] auxFiles;
        INPUTMODE mode;
        public PDBFiles()
        {
            Settings set = new Settings();
            set.Load();
            mode = set.mode;
        }
        public Dictionary<string, int> MakeAllAtomsDic()
        {
            Dictionary<string,int> dic=new Dictionary<string,int>();
            int num = 0;


            foreach(var item in molDic)
            {
                foreach (var res in item.Value.mol.Residues)
                    foreach (var it in res.Atoms)
                        if (!dic.ContainsKey(it.AtomName))
                            dic.Add(it.AtomName,num++);

            }


            return dic;
        }
        public void AddPDB(List<string> structures, PDBMODE flag,ref int currentProgress)
        {
            molDic = new Dictionary<string, MolData>(structures.Count);
            foreach (var item in structures)
            {
                AddPDB(item, flag);
                currentProgress++;
            }            
        }

        public void AddPDB(string []structures,PDBMODE flag,ref int currentProgress)
        {
            if (mode == INPUTMODE.RNA)
                flag = PDBMODE.ONLY_P;
            
            foreach (var item in structures)
            {
                AddPDB(item, flag);
                currentProgress++;
            }
        }
        public string AddPDB(string fileName,PDBMODE flag,CHAIN_MODE flagChain=CHAIN_MODE.SINGLE)
        {
            string name = Path.GetFileName(fileName);
            if(molDic!=null && molDic.ContainsKey(name))            
                return name;
            try
            {
                MolData molD = new MolData(fileName, flag,mode,flagChain);
                if ((molD.mol == null || molD.mol.Chains == null || molD.mol.Chains.Count == 0))
                {
                    ErrorBase.AddErrors("PDB reading: file " + fileName + " is removed from consideration it looks that it has wrong format");
                    return null;
                }
                molDic.Add(name, molD);

                return name;

            }
            catch (IncorrectSideChainException ex)
            {
                ErrorBase.AddErrors("PDB reading: file " + fileName + " is removed from consideration because\n" + ex.Message);
            }
            catch (Exception ee)
            {
                ErrorBase.AddErrors("PDB reading: file " + fileName + " is removed from consideration because\n" + ee.Message);
            }
            return null;
        }
        public void ReadRefSeq(string fileName)
        {
            refSeq = Alignment.ReadRefSeq(fileName);
        }
        public string AddPDB(MemoryStream stream, PDBMODE flag,string modelName)
        {
            try
            {
                MolData molD = new MolData();
                if (!molD.ReadMolData(stream, flag,modelName))
                {
                    return null;
                }
                molDic.Add(modelName, molD);

                return modelName;

            }
            catch (IncorrectSideChainException ex)
            {
                ErrorBase.AddErrors("PDB reading: file " + modelName+ " is removed from consideration because\n" + ex.Message);
            }
            return null;
        }
       
        public void FindReferenceSeq()
        {
            int remLen=0;
            foreach (var item in molDic.Keys)
            {               
                if(remLen<molDic[item].molLength)
                {
                    remLen = molDic[item].molLength;
                    refStuctName = item;                 
                }
            }
            if(refSeq==null || refSeq.Length<molDic[refStuctName].mol.Chains[0].chainSequence.Length)
                refSeq = molDic[refStuctName].mol.Chains[0].chainSequence;
        }
        public void EmptyAlignment(string item)
        {
            int len = molDic[item].mol.Chains[0].chainSequence.Length;
            molDic[item].indexMol = new int[len];
            for (int j = 0, count = 0; j <len; j++)
                    molDic[item].indexMol[j] = count++;

        }
        private void AddAlignment(string item)
        {
            MAlignment align = new MAlignment(molDic[refStuctName].mol.Chains[0].chainSequence.Length);
            molDic[item].indexMol = new int[molDic[refStuctName].mol.Chains[0].chainSequence.Length];
            string ss = align.Align(molDic[refStuctName].mol.Chains[0].chainSequence, molDic[item].mol.Chains[0].chainSequence).seq2;
            for (int j = 0, count = 0; j < ss.Length; j++)
            {
                if (ss[j] != '-')
                    molDic[item].indexMol[j] = count++;
                else
                    molDic[item].indexMol[j] = -1;
            }

        }
        private void PDBSAlign(object o)
        {
            Params p=(Params) o;
            List<string> toRemove = new List<string>();
            MAlignment align = new MAlignment(refSeq.Length);
            DebugClass.WriteMessage("Started");



            foreach (var item in auxFiles[p.k])
            {
                molDic[item].indexMol = new int[refSeq.Length];
                MAlignment.alignSeq alignRes = align.Align(refSeq, molDic[item].mol.Chains[0].chainSequence);


                if (alignRes.seq1.Contains('-'))
                {
                    ErrorBase.AddErrors("Reference structure " + refStuctName + " cannot be used as reference to " + item + " structure because after alignment gaps to referense structure must be added what is not allowed");
                    toRemove.Add(item);
                    continue;
                }
                //string ss = align.Align(molDic[refStuctName].mol.Chains[0].chainSequence, molDic[item].mol.Chains[0].chainSequence).seq2;
                string ss = alignRes.seq2;
                for (int j = 0, count = 0; j < ss.Length; j++)
                {
                    if (ss[j] != '-')
                        molDic[item].indexMol[j] = count++;
                    else
                        molDic[item].indexMol[j] = -1;
                }
            }
                DebugClass.WriteMessage("Almost done");
                foreach (var item in toRemove)
                    molDic.Remove(item);
           
            resetEvents[p.k].Set();
        }
        
        public void MakeAlignment(string alignFile)
        {
            List<string> toRemove = new List<string>();

            if (molDic.Count == 0)
                return;

            if (molDic.Count == 1)
            {
                foreach (var item in molDic.Keys)
                {
                    molDic[item].indexMol = new int[molDic[item].mol.Chains[0].chainSequence.Length];
                    for (int i = 0; i < molDic[item].indexMol.Length; i++)
                        molDic[item].indexMol[i] = i;

                }
            }
            
            if (alignFile != null && alignFile.Length > 0)
            {
                Dictionary<string, string> al = Alignment.ReadAlignment(alignFile);
                foreach (var item in molDic.Keys)
                {
                    string ss;
                    if (al.ContainsKey(item))
                    {
                        ss = al[item];
                        molDic[item].indexMol = new int[ss.Length];

                        for (int j = 0, count = 0; j < ss.Length; j++)
                        {
                            if (ss[j] != '-')
                                molDic[item].indexMol[j] = count++;
                            else
                                molDic[item].indexMol[j] = -1;
                        }
                    }
                }
            }
            else
                if (molDic.Count == 2)
                {
                    using (MAlignment align = new MAlignment(refSeq.Length))
                    {

                        MAlignment.alignSeq alignRes;
                        List<string> molItems = new List<string>(molDic.Keys);
                        string str1, str2;
                        if (molDic[molItems[0]].mol.Chains[0].chainSequence.Length > molDic[molItems[1]].mol.Chains[0].chainSequence.Length)
                        {
                            str1 = molItems[0];
                            str2 = molItems[1];
                        }
                        else
                        {
                            str1 = molItems[1];
                            str2 = molItems[0];
                        }
                        alignRes = align.Align(molDic[str1].mol.Chains[0].chainSequence, molDic[str2].mol.Chains[0].chainSequence);
                        molDic[str1].indexMol = new int[alignRes.seq1.Length];
                        molDic[str2].indexMol = new int[alignRes.seq1.Length];
                        //string ss = align.Align(molDic[refStuctName].mol.Chains[0].chainSequence, molDic[item].mol.Chains[0].chainSequence).seq2;
                        string ss = alignRes.seq2;
                        for (int j = 0, count = 0; j < ss.Length; j++)
                        {
                            if (ss[j] != '-')
                                molDic[str2].indexMol[j] = count++;
                            else
                                molDic[str2].indexMol[j] = -1;
                        }
                        ss = alignRes.seq1;
                        for (int j = 0, count = 0; j < ss.Length; j++)
                        {
                            if (ss[j] != '-')
                                molDic[str1].indexMol[j] = count++;
                            else
                                molDic[str1].indexMol[j] = -1;
                        }
                    }
                }
                else
                {
                    Settings set = new Settings();
                    set.Load();
                    int threadNumbers = set.numberOfCores;
                //return;
                    auxFiles = new List<string>[threadNumbers];
                    resetEvents = new ManualResetEvent[threadNumbers];
                    List<string> allFiles = new List<string>(molDic.Keys);
                    for(int i=0;i<threadNumbers;i++)
                    {
                        auxFiles[i] = new List<string>((i + 1) * molDic.Count / threadNumbers - i * molDic.Count / threadNumbers);
                        for (int j = i * allFiles.Count / threadNumbers; j < (i + 1) * allFiles.Count / threadNumbers; j++)
                            auxFiles[i].Add(allFiles[j]);
                    }
                    DebugClass.WriteMessage("Alignment start =" + threadNumbers);
                    //PDBSAlign(auxFiles)

                    for (int i = 0; i < threadNumbers; i++)
                    {
                        //PDBSAlign(auxFiles[i]);
                        Params p = new Params();
                        p.k=i;
                        resetEvents[p.k] = new ManualResetEvent(false);
                        ThreadPool.QueueUserWorkItem(new WaitCallback(PDBSAlign), (object)p);
                        //  startProg = Task.Factory.StartNew(() => RunMakeProfiles(fileName,auxFiles,k));
                        // runnigTask[i] = startProg;
                    }
                    for (int i = 0; i < threadNumbers; i++)
                        resetEvents[i].WaitOne();
                    //WaitHandle.WaitAll(resetEvents);
                    DebugClass.WriteMessage("Alignment stop");
                }
        }
               

    }


}
