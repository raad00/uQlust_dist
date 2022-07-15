using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Text.RegularExpressions;
using System.Globalization;

namespace uQlustCore.Profiles
{
    
    public class UserDefinedProfile<T> :InternalProfileBase<T> where T: struct
    {
        //string profileName;

//        InternalProfilesManager manager = new InternalProfilesManager();
        private static string ProfileName = "User defined profile";

        public UserDefinedProfile()
        {           
          
            AddInternalProfiles();
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.PROTEIN);
            destination.Add(INPUTMODE.RNA);
            destination.Add(INPUTMODE.USER_DEFINED);
        }
        public static string GetName()
        {
            return ProfileName;
        }
        public override int Run(object processParams)
        {
            return 0;
        }
        public override void Run(DCDFile dcd)
        {
            throw new NotImplementedException();
        }
        public override void RunThreads(string fileName)
        {

        }
        public override Dictionary<string, protInfo<T>> GetProfile(profileNode node, string fileName, DCDFile dcd)
        {
            Dictionary<string, protInfo<T>> dic = new Dictionary<string, protInfo<T>>();
            StreamReader wr;
            DebugClass.WriteMessage("profile" + fileName);
            wr = new StreamReader(fileName);
            
            protInfo<T> info;
            string line = wr.ReadLine();
            string name = "";
            int prevSize = -1;
            string seq = "";
            List<string> profile = new List<string>();
            List<T> newProfile = new List<T>();
            while (line != null)
            {
                if (line.Contains(">"))
                {
                    if (name.Length > 0)
                    {
                        info = new protInfo<T>();
                        info.sequence = null;
                        info.profile = newProfile;
                        if (dic.ContainsKey(name))
                        {
                            //wr.Close();
                            //throw new Exception("The nameof profile must be unique, name: " + name + " already exists in " + fileName);
                        }
                        else
                        {
                            if (prevSize > 0)
                            {
                                if (newProfile.Count == prevSize)
                                    dic.Add(name, info);
                            }
                            else
                            {
                                prevSize = newProfile.Count;
                                dic.Add(name, info);
                            }

                        }
                    }
                    newProfile = new List<T>();
                    name = line.Replace(">", "");
                    line = wr.ReadLine();
                }
                if (line != null && line.Contains(" profile "))
                {

                    int index = line.IndexOf(" profile ");
                    if (index == -1)
                    {
                        wr.Close();
                        throw new Exception("Incorrect file format you have to have 'profile' word!");
                    }

                    string cLine = line.Remove(0, " profile ".Length + index);

                    cLine = Regex.Replace(cLine, @"\s+", " ");
                    if (line.Length == 0)
                        continue;
                    cLine = cLine.Trim();
                    string[] aux;
                    if (cLine.Contains(' '))
                        aux = cLine.Split(' ');
                    else
                    {
                        char[] charArray = cLine.ToCharArray();
                        aux = cLine.Select(x => x.ToString()).ToArray();
                    }

                    profile.Clear();

                    foreach (var item in aux)
                        profile.Add(item);


                    newProfile = new List<T>();
                    T z = default(T);
                    if (z is byte)
                        for (int i = 0; i < profile.Count; i++)
                            if (node.ContainsState(profile[i].ToString()))
                                newProfile.Add((T)Convert.ChangeType(node.GetCodedState(node.states[profile[i].ToString()]), typeof(T)));
                            else
                                if (profile[i].ToString() != "-" && profile[i].ToString() != "X")
                                throw new Exception("Unknow state " + profile[i].ToString() + " in " + node.profName + " profile!");
                            else
                                newProfile.Add((T)Convert.ChangeType(0.0, typeof(T)));
                    else
                        for (int i = 0; i < profile.Count; i++)
                            newProfile.Add((T)Convert.ChangeType(profile[i], typeof(T), CultureInfo.InvariantCulture));


                }
                line = wr.ReadLine();
            }
            if (newProfile.Count > 0 && !dic.ContainsKey(name))
            {
                info = new protInfo<T>();
                info.sequence = seq;
                info.profile = newProfile;
                dic.Add(name, info);
            }
            DebugClass.WriteMessage("number of profiles " + dic.Keys.Count);

            wr.Close();            

            return dic;
        }
        public override void AddInternalProfiles()
        {
            profileNode node = new profileNode();

            node.profName = ProfileName;
            node.internalName = ProfileName;           
            InternalProfilesManager.AddNodeToList(node, this.GetType().FullName);
        }
        public override void RemoveInternalProfiles()
        {            
            InternalProfilesManager.RemoveNodeFromList(ProfileName);
        }
        public override List<INPUTMODE> GetMode()
        {
            List<INPUTMODE> oList = new List<INPUTMODE>();
            oList.Add(INPUTMODE.PROTEIN);
            oList.Add(INPUTMODE.RNA);
            oList.Add(INPUTMODE.USER_DEFINED);
            return oList;
        }
    }
}
