using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;
using uQlustCore.PDB;
using uQlustCore.dcd;

namespace uQlustCore.Profiles
{
    public class CAProfiles:ContactProfile
    {
        //Settings dirSettings = new Settings();
        //InternalProfilesManager manager = new InternalProfilesManager();
        protected static List<char> states =new List<char>() { 'H', 'J', 'K', 'L', 'E', 'R', 'T', 'Y', 'U', 'N' };
    public CAProfiles()
        {
            dirSettings.Load();
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.PROTEIN);
           // pdbs = new PDBFiles();
            profileName = "SS_CA_SA_CA";
            contactProfile = "SA CA profile ";
            ssProfile = "SS CA profile ";

            AddInternalProfiles();       

        }
       protected KeyValuePair<char[],string> GenerateStates(MolData molDic)
        {
            Dictionary<int, List<int>> contacts = new Dictionary<int, List<int>>();
            char[] statesTab = Enumerable.Repeat<char>('N', molDic.mol.Residues.Count).ToArray();
            double[] dist2;
            molDic.CreateCAContactMap(8.5f, false);
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
            dist2 = new double[molDic.mol.Residues.Count - 2];
            for (int i = 0; i < molDic.mol.Residues.Count - 2; i++)
            {

                Atom aux1 = molDic.mol.Residues[i].Atoms[0];
                Atom aux2 = molDic.mol.Residues[i + 2].Atoms[0];
                double sum = (aux1.Position.X - aux2.Position.X) * (aux1.Position.X - aux2.Position.X);
                sum += (aux1.Position.Y - aux2.Position.Y) * (aux1.Position.Y - aux2.Position.Y);
                sum += (aux1.Position.Z - aux2.Position.Z) * (aux1.Position.Z - aux2.Position.Z);
                sum = Math.Sqrt(sum);
                dist2[i] = sum;
            }
            for (int i = 0; i < molDic.mol.Residues.Count - 4; i++)
            {
                Atom aux1 = molDic.mol.Residues[i].Atoms[0];
                Atom aux2 = molDic.mol.Residues[i + 4].Atoms[0];
                double sum = (aux1.Position.X - aux2.Position.X) * (aux1.Position.X - aux2.Position.X);
                sum += (aux1.Position.Y - aux2.Position.Y) * (aux1.Position.Y - aux2.Position.Y);
                sum += (aux1.Position.Z - aux2.Position.Z) * (aux1.Position.Z - aux2.Position.Z);
                sum = Math.Sqrt(sum);
                if (dist2[i] > 4 && dist2[i] < 8)
                {
                    if (dist2[i] < 6)
                    {
                        if (sum > 4 && sum < 14)
                        {
                            if (sum < 7)
                                statesTab[i] = 'H';
                            else
                                if (sum < 9)
                                statesTab[i] = 'J';
                            else
                                    if (sum < 11)
                                statesTab[i] = 'K';
                            else
                                        if (sum < 13)
                                statesTab[i] = 'L';
                        }
                        else
                            statesTab[i] = 'U';//unphysical
                    }
                    else
                        if (sum > 4 && sum < 14)
                    {
                        if (sum < 7)
                            statesTab[i] = 'E';
                        else
                            if (sum < 9)
                            statesTab[i] = 'R';
                        else
                                if (sum < 11)
                            statesTab[i] = 'T';
                        else
                                    if (sum < 13)
                            statesTab[i] = 'Y';
                    }

                }
                else
                    statesTab[i] = 'U';//unphysical                                  

            }

            return new KeyValuePair<char[], string>(statesTab, profile);
        }
        protected override void MakeProfiles(string strName, MolData molDic,StreamWriter wr)
         {
              Dictionary<int, List<int>> contacts = new Dictionary<int, List<int>>();
            
             if(molDic!=null)             
             {

                KeyValuePair<char [],string> x = GenerateStates(molDic);
                string ss = new string(x.Key);

                 if (x.Value.Length > 0)
                 {
                     wr.WriteLine(">" + strName);
                     wr.WriteLine(contactProfile + x.Value);                
                     wr.WriteLine(ssProfile + ss);
                     wr.WriteLine(SEQprofile + molDic.mol.Chains[0].chainSequence);

                 }

                 molDic.CleanMolData();
             }

         }
         public override void AddInternalProfiles()
         {
             profileNode node = new profileNode();

             node.profName = "SA CA";
             node.internalName = "SA CA";
             for (int i = 0; i < 10; i++)
                 node.AddStateItem(i.ToString(), i.ToString());

              InternalProfilesManager.AddNodeToList(node, typeof(CAProfiles).FullName);

             node = new profileNode();

             node.profName = "SS CA";
             node.internalName = "SS CA";
             node.AddStateItem("H", "H");
             node.AddStateItem("J", "J");
             node.AddStateItem("K", "K");
             node.AddStateItem("L", "L");
             node.AddStateItem("E", "E");
             node.AddStateItem("R", "R");
             node.AddStateItem("T", "T");
             node.AddStateItem("Y", "Y");
             node.AddStateItem("U", "U");
             node.AddStateItem("N", "N");
              InternalProfilesManager.AddNodeToList(node, typeof(CAProfiles).FullName);

         }
         public override void RemoveInternalProfiles()
         {
             InternalProfilesManager.RemoveNodeFromList("SS CA");
             InternalProfilesManager.RemoveNodeFromList("SA CA");
         }
         public override List<INPUTMODE> GetMode()
         {
             return base.GetMode();
         }
    }
}
