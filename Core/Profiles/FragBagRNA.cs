using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Reflection;
namespace uQlustCore.Profiles
{
    class FragBagRNA:FragBagProfile
    {
        public FragBagRNA()
        {
            dirSettings.Load();
            fragBagProfile = "FragBagRNA profile ";
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.RNA);
            profileName = "FragBagRNA";
            contactProfile = "FragBagRNA profile ";
            AddInternalProfiles();
            maxV = 1;

        }
        public override bool CheckIfAvailable()
        {
            return false;
            if (!Directory.Exists(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) +
                    Path.DirectorySeparatorChar + "fragLib"))
                if (!Directory.Exists("C:\\Projects\\UQlast\\fragLib"))
                    throw new Exception("Directory fragLib not exists. Profile FragBag cannot be used!");

            //ReadLibrary("fragBagRNAv2.txt");
            ReadLibrary("fragBagRNA.txt");
        }
        public override void AddInternalProfiles()
        {
            profileNode node = new profileNode();

            node.profName = "FragBagRNA";
            node.internalName = "FragBagRNA";
            for (int i = 0; i < 255; i++)
                node.AddStateItem(i.ToString(), i.ToString(), true);

            InternalProfilesManager.AddNodeToList(node, typeof(FragBagRNA).FullName);

        }
        public override Dictionary<string, protInfo<byte>> GetProfile(profileNode node, string listFile, DCDFile dcd = null)
        {
            Dictionary<string, protInfo<byte>> res = ReadProfile(node, listFile, dcd);
            //int[] index = new int[29] { 21, 64, 29, 2, 33, 57, 77, 84, 14, 44, 81, 0, 90, 74, 30, 71, 76, 8, 48, 61,86,41,3,7,15,20,23,25,51};//for rna
            //int[] index = new int[20] { 90, 203, 95, 39, 59, 97, 12, 49, 122, 195, 205, 91, 221, 187, 230, 198, 84, 96, 124, 229};//for rfam
            // int[] index = new int[10] { 104,195,81,59,91,95,150,186,2,11};//for ribos
            //int[] index = new int[30] { 158,172,224,120,94,20,26,44,248,173,75,92,201,170,243,118,195,182,116,133,217,137,25,177,188,105,89,189,202,56}; //tylko na rna
            //int[] index = new int[30] { 75,118,198,224,219,173,218,145,100,44,124,152,137,94,182,20,80,17,202,92,165,166,238,247,97,129,172,216,133,248 };//Wersja na calosi rna40
            //int[] index = new int[30] { 21, 9, 33, 2, 79, 29, 90, 64, 0, 81, 86, 14, 44, 84, 74, 30, 76, 3, 82, 77, 57, 34, 26, 85, 61,8, 71, 20, 47, 56 };//Wersja na calosi rna40 biblioteka ver1
            //int[] index = new int[30] { 75, 118, 124, 173, 20, 100, 44, 198, 145, 219, 218, 152, 224, 202, 17, 137, 80, 182, 92, 165, 94, 166, 97, 216, 89, 247, 34, 125, 172, 238 };//half rna40
            /*List<string> keys = new List<string>(res.Keys);

            foreach(var item in keys)
            {
                List<byte> newProfile = new List<byte>();
                for (int i = 0; i < index.Length; i++)
                    newProfile.Add(res[item].profile[index[i]]);
                protInfo xx = res[item];
                xx.profile = newProfile;
                res[item] = xx;
            }*/
            //res = RearangeColumnOrder(res);
            // res = RearangeColumnOrder(res, "C:\\Projects\\listIndex");
            res = ProfileStat.RearangeStates(res, 0.51);
            //res = ProfileStat.SelectFeatures(res,30);
            return res;
        }

    }
}
