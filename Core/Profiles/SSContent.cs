using System;
using System.Collections.Generic;
using uQlustCore.PDB;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace uQlustCore.Profiles
{
    class SSComposition:CAProfiles
    {
        public SSComposition()
        {
            dirSettings.Load();
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.PROTEIN);
            // pdbs = new PDBFiles();
            profileName = "SS_Comp";
            contactProfile = "SS_Comp profile ";
            ssProfile = "SS_Comp profile ";

            AddInternalProfiles();

        }
        public override void AddInternalProfiles()
        {
            profileNode node = new profileNode();

            node.profName = "SS_Comp";
            node.internalName = "SS_Comp";
         
            for (int i = 0; i < 10; i++)
                node.AddStateItem(i.ToString(), i.ToString());

            InternalProfilesManager.AddNodeToList(node, typeof(SSComposition).FullName);

        }
        protected override void MakeProfiles(string strName, MolData molDic, StreamWriter wr)
        {
            Dictionary<int, List<int>> contacts = new Dictionary<int, List<int>>();

            if (molDic != null)
            {

                KeyValuePair<char[], string> x = GenerateStates(molDic);
                string ss = new string(x.Key);
                Dictionary<char, int> counter = new Dictionary<char, int>();
                for(int i=0;i<ss.Length;i++)
                {
                    if (counter.ContainsKey(ss[i]))
                        counter[ss[i]]++;
                    else
                        counter.Add(ss[i], 1);
                }
                List<char> keyList = new List<char>(counter.Keys);
                foreach (var item in keyList)
                {
                    double res= ((double)counter[item])/ ss.Length;                    
                    res *= 10;
                    counter[item] = (int)Math.Floor(res);
                }
                if (ss.Length > 0)
                {
                    wr.WriteLine(">" + strName);
                    string txt = "";
                    for(int i=0;i<states.Count-1;i++)
                    {
                        if (counter.ContainsKey(states[i]))                        
                            txt += counter[states[i]] + " ";                        
                        else
                            txt += "0 ";
                    }
                    if (counter.ContainsKey(states[states.Count-1]))
                        txt += counter[states[states.Count - 1]];
                    else
                        txt += "0";
                    wr.WriteLine(ssProfile + txt);

                }

                molDic.CleanMolData();
            }

        }
        public override void RemoveInternalProfiles()
        {
            InternalProfilesManager.RemoveNodeFromList("SS_Comp");
        }

    }
}
