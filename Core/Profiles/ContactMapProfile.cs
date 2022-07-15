using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;
using uQlustCore.PDB;
using uQlustCore.dcd;
using System.Threading;

namespace uQlustCore.Profiles
{
  
    public class ContactMapProfile:ContactProfile
    {
        protected int[] contOne = null;
        
        protected byte [][] contact;
        protected char[][] contactToString;
        protected List<string>[] auxFiles;
       //protected Settings dirSettings = new Settings();
       //InternalProfilesManager manager = new InternalProfilesManager();      
       protected PDBFiles pdbs;
       protected ManualResetEvent[] resetEvents;
       public ContactMapProfile()
        {
            dirSettings.Load();
            profileName = "ContactMap";
            contactProfile = "ContactMap profile ";

            AddInternalProfiles();
            resetEvents = new ManualResetEvent[threadNumbers];
            maxV = 1;

        }
        protected void ReadPdbs(List <string> structures)
        {
           try
           {
               pdbs = new PDBFiles();
               foreach (var item in structures)
               {
                    pdbs.AddPDB(item, PDBMODE.ONLY_SEQ);
                    currentProgress++;
               }
             //  pdbs.AddPDB(structures, PDBMODE.ONLY_CA,ref currentV);
               //pdbs.AddPDB(structures, PDBMODE.ONLY_SEQ, ref currentV);
           }
           catch (Exception ex)
           {
               throw new Exception(ex.Message);
           }
           
        }
        protected void RunMakeProfiles(object o)
        {
            Params p = (Params)o;
            string fileN = GetProfileFileName(p.fileName) + "_" + p.k;

            using (StreamWriter wr = new StreamWriter(fileN))
            {
                foreach (var item in auxFiles[p.k])
                {
                    string nn = Path.GetFileName(item);
                    if(pdbs.molDic.ContainsKey(nn))
                        MakeProfiles(item, pdbs.molDic[nn], wr, p.k);

                    Interlocked.Increment(ref currentProgress);
                }

                wr.Close();
            }
            foreach (var item in auxFiles[p.k])
                pdbs.molDic.Remove(item);

            GC.Collect();
            resetEvents[p.k].Set();
        }
       public override void RunThreads(string fileName)
       {
           List<string> files = CheckFile(fileName);

           if (files.Count == 0)
               return;
           
           //Task[] runnigTask = new Task[threadNumbers];
           //Task startProg;
           maxV = files.Count*2;
           ReadPdbs(files);
           if (pdbs.molDic.Count == 0)
               return;
           string aux = Path.GetDirectoryName(files[0]).TrimEnd(Path.DirectorySeparatorChar);

           string refSeqFile = aux+".ref";
           if(File.Exists(refSeqFile))
            pdbs.ReadRefSeq(refSeqFile);
           else
               pdbs.FindReferenceSeq();
           pdbs.MakeAlignment(null);

           if (contOne == null)
           {
               string molDicKey="";
               foreach (var item in pdbs.molDic.Keys)
               {
                   molDicKey = item;
                   break;
               }
               int len = pdbs.molDic[molDicKey].indexMol.Length;
               contOne = new int[len * (len + 1) / 2];
               contact = new byte[threadNumbers][];
               contactToString=new char [threadNumbers][];

               for(int i=0;i<threadNumbers;i++)
               {
                   contact[i]=new byte [len * (len + 1) / 2];
                   contactToString[i]=new char [len * (len + 1)];
               }
               
               for (int i = 0; i < contOne.Length; i++)
                   contOne[i] =0;
               for (int i = 0; i < threadNumbers; i++)
                   for (int j = 0; j < contOne.Length; j++)
                       contact[i][j] = 0;
              
           }
           
           auxFiles =new List<string>[threadNumbers];
           List<string> allFiles = new List<string>(pdbs.molDic.Keys);

           

           for (int i = 0; i < threadNumbers; i++)
           {
               auxFiles[i] = new List<string>((i + 1) * pdbs.molDic.Count / threadNumbers - i * pdbs.molDic.Count / threadNumbers);               
               for (int j = i * allFiles.Count / threadNumbers; j < (i + 1) * allFiles.Count / threadNumbers; j++)
                       auxFiles[i].Add(files[j]);
           }
           for (int i = 0; i < threadNumbers; i++)
           {
               Params p = new Params();
               p.fileName = fileName;
               p.k = i;

               resetEvents[p.k] = new ManualResetEvent(false);
               ThreadPool.QueueUserWorkItem(new WaitCallback(RunMakeProfiles), (object)p);
           }
           for (int i = 0; i < threadNumbers; i++)
               resetEvents[i].WaitOne();

           //JoinFiles(fileName);

           CuttProfiles(fileName);
           currentProgress = maxV;
       }

       public override int Run(object processParams)
       {
           return 0;
       }

        private void CuttProfiles(string fileName)
       {
           using (StreamWriter wr = new StreamWriter(GetProfileFileName(fileName), true))
           {
               if (wr == null)
                   throw new Exception("Cannot open file: " + GetProfileFileName(fileName));
               string profileName = contactProfile;
               profileName = profileName.Trim();
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
                           if (!line.Contains(">"))
                           {
                               string[] aux = line.Split(' ');
                               wr.Write(profileName);
                               for (int j = 0; j < aux.Length; j++)
                               {
                                   if (contOne[j] > 0)
                                       wr.Write(" "+aux[j]);
                               }
                               if (contOne[aux.Length - 1] > 0)
                                   wr.WriteLine(" "+aux[aux.Length - 1]);
                               else
                                   wr.WriteLine();
                           }
                           else
                           {
                               wr.WriteLine(line);
                           }
                           line = rr.ReadLine();
                       }
                       rr.Close();
                       File.Delete(fileN);
                   }
               }
               wr.Close();
           }
       }
       protected virtual void GenerateContactMap(MolData mol,int k)
        {
            mol.CreateFullContactMap(8.5f, contact[k],"CA");               
        }
       protected void MakeProfiles(string strName,MolData molDic,StreamWriter wr,int k)
       {
           if (molDic != null)
           {
               PDBFiles pdb = new PDBFiles();
               pdb.AddPDB(strName, PDBMODE.ONLY_CA);
               DebugClass.WriteMessage("Make Started");
               //int []contact=pdbs.molDic[strName].CreateFullContactMap(8.5f);
               //molDic.CreateFullContactMap(9.5f,contact[k]);//correct
               List<string> cc=new List<string>(pdb.molDic.Keys);
               pdb.molDic[cc[0]].indexMol = molDic.indexMol;
               GenerateContactMap(pdb.molDic[cc[0]],k);
               //pdb.molDic[cc[0]].CreateFullContactMap(8.5f, contact[k]);               
               //molDic.CreateFullContactMap(8.5f,contact[k]);
               //pdbs.molDic[strName].CreateContactMap(9.5f, "CB");

            
               for (int i = 0; i < contact[k].Length; i++)
                   if (contact[k][i] == 1)
                       contOne[i]++;
               int j = 0;
               for (int i = 0; i < contact[k].Length - 1; i++)
               {
                   if (contact[k][i] == 1)
                       contactToString[k][j++] = '1';
                   else
                       contactToString[k][j++] = '0';

                   contactToString[k][j++] = ' ';
               }
               if (contact[k][contact[k].Length - 1] == 1)
                   contactToString[k][j] = '1';
               else
                   contactToString[k][j] = '0';

               string all = new string(contactToString[k], 0, j);

               wr.WriteLine(">" + cc[0]);
               //for (int i = 0; i < contact[k].Length-1; i++)               
                 //  wr.Write(contact[k][i]+" ");
               //wr.Write(contact[k][contact.Length-1]);
               //all = all.Trim();
               wr.WriteLine(all.Trim());
               //wr.WriteLine();
               //molDic.CleanMolData();
               DebugClass.WriteMessage("Make finished");
           }
          
       }

       public override void AddInternalProfiles()
       {
           profileNode node = new profileNode();

           node.profName = profileName;
           node.internalName = profileName;

           node.AddStateItem("0", "0");
           node.AddStateItem("1", "1");
           node.AddStateItem("2", "2");
           node.AddStateItem("3", "3");
            InternalProfilesManager.AddNodeToList(node, typeof(ContactMapProfile).FullName);           

       }
       public override void RemoveInternalProfiles()
       {
           InternalProfilesManager.RemoveNodeFromList("ContactMap");
       }
    }
}
