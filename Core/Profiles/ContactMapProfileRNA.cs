using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using uQlustCore.PDB;

namespace uQlustCore.Profiles
{
    public class ContactMapProfileRNA:ContactMapProfile
    {
        public ContactMapProfileRNA()
        {
            dirSettings.Load();
            profileName = "ContactMapRNA";
            contactProfile = "ContactMapRNA profile ";
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.RNA);
            AddInternalProfiles();
            resetEvents = new ManualResetEvent[threadNumbers];
            maxV = 1;

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
            InternalProfilesManager.AddNodeToList(node, typeof(ContactMapProfileRNA).FullName);

        }
        protected override void GenerateContactMap(MolData mol, int k)
        {
            mol.CreateFullContactMap(15.5f, contact[k],"P");
        }
    }
}
