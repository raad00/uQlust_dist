using System;
using System.IO;
using System.Diagnostics;
using System.Collections;
using System.Collections.Generic;
using System.Text.RegularExpressions;
using System.Text;
using System.Data;
using uQlustCore.Interface;

namespace uQlustCore
{
    public class Alignment : IProgressBar
    {
        Settings dirSettings = new Settings();
        Dictionary<string, Dictionary<string, string>> align=new Dictionary<string, Dictionary<string, string>>();
        int gcCounter = 0;
        string refSeq = null;
        int maxV;
		public ProfileTree r;
             

        public Alignment()
        {
            r = new ProfileTree();
            dirSettings.Load();
        }

        public double ProgressUpdate()
        {
            return r.ProgressUpdate();
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }


        public void Prepare(DCDFile dcd,Settings dirSettings, string profName)
        {
            StartAlignment(dirSettings, profName, dcd, null, null);
        }       
		public void Prepare(string pathName,Settings dirSettings,string profName)
		{
            StartAlignment(dirSettings, profName, null, pathName, null);
        }
        public void Prepare(List<string>fileNames, Settings dirSettings, string profName)
        {
            StartAlignment(dirSettings, profName, null, null, fileNames);
        }
        public void Prepare(List<string> names,string alignFile,string profName)
        {

        }
        public void Clean()
       {
           align.Clear();
           r.protCombineStates.Clear();
           r.profiles.Clear();          
       }
        public void Prepare(List<KeyValuePair<string,string>> profilesStr,string profName,string profFile)
        {
            r = new ProfileTree();
            r.LoadProfiles(profFile);            
            
            foreach(var item in profilesStr)
            {
                List<string> aux = new List<string>(item.Value.Length);
                for (int i = 0; i < item.Value.Length; i++)
                    aux.Add(item.Value[i].ToString());

               r.AddItemsCombineStates(item.Key,aux);
                                
            }
        }
        public void Prepare(string profilesFile, string profName)
        {
            r.LoadProfiles(profName);
            r.listFile = profilesFile;
            DebugClass.WriteMessage("profiles gen started "+profName);
            r.MakeProfiles();
            DebugClass.WriteMessage("Prfofiles end");
        }
        public void AddProfiles(Alignment al)
        {
            r.JoinProfiles(al.r.profiles);
        }
        private void StartAlignment(Settings dirSettings, string profName, DCDFile dcd, string dirName, List<string> fileNames )
        {

            DebugClass.WriteMessage("Start align");
            string refFile = null;
            this.dirSettings = dirSettings;            
           // r = new ProfileTree();
            r.LoadProfiles(profName);
            if (dcd != null)
            {
                DebugClass.WriteMessage("profiles gen started");
                r.PrepareProfiles(dcd);
            }
            else
                if (dirName != null)
                {
                    refFile = dirName + ".ref";
                    DebugClass.WriteMessage("profiles gen started");
                    r.PrepareProfiles(dirName);
                    DebugClass.WriteMessage("finished");
                    refSeq=ReadRefSeq(refFile);
                }
                else
                {
                    maxV = fileNames.Count;
                    string name = fileNames[0];
                    if(fileNames[0].Contains("|"))
                    {
                        string[] aux = fileNames[0].Split('|');
                        name = aux[0];
                    }
                    refFile = Directory.GetParent(name).ToString() + ".ref";
                    DebugClass.WriteMessage("profiles gen started");
                    r.PrepareProfiles(fileNames);
                    DebugClass.WriteMessage("finished");
                    refSeq = ReadRefSeq(refFile);
                }
            DebugClass.WriteMessage("Prfofiles end"); 


        }
		public Dictionary<string,Dictionary<string,string>> GetAlignment()
		{
			return align;
		}
        public List<string> GetStructuresNames()
        {
            List<string> names;
            List<string> prof = new List<string>(r.profiles.Keys);
            names = new List<string>(r.profiles[prof[0]].Keys);
            return names;
        }
		public Dictionary<string, List<byte>> GetStateAlign()
		{
			return r.protCombineStates;
		}
        public static string ReadRefSeq(string fileName)
        {
            string rSeq;
            if (!File.Exists(fileName))
                return null;

            StreamReader stR = new StreamReader(fileName);
            rSeq = stR.ReadLine();
            stR.Close();
            return rSeq;
        }
        public static Dictionary<string, string> ReadAlignment(string fileName)
        {
            Dictionary<string, string> alignLoc = new Dictionary<string, string>();
            StreamReader file_in = null;
            string line, remName = "";

            try
            {
                file_in = new StreamReader(fileName);
                line = file_in.ReadLine();
                while (line != null)
                {
                    if (line.Contains(">"))
                    {
                        string name = line.Substring(1, line.Length-1);
                        string profile = "";
                        line = file_in.ReadLine();
                        while (line!=null && !(line.Contains(">")))
                        {
                            profile+=line;
                            line = file_in.ReadLine();                                                        
                        }
                        if (alignLoc.Count == 0)
                            remName = name;
                        else
                        {
                            profile.Replace("\n", "");
                            profile.Replace(" ", "");
                            if (profile.Length != alignLoc[remName].Length)
                                throw new Exception("Alignment incorrect for " + remName + " and "+ name + "\nDifferent number of symbols in the alignment!");
                        }
                        alignLoc.Add(name, profile);
                    }
                    else
                        line = file_in.ReadLine();
                }
            }
            finally
            {
                if (file_in != null)
                    file_in.Close();
            }

            return alignLoc;

        }

      
		public void MyAlign(string alignFile)
		{
            bool test=false;
            //Check if there is sequence that could be aligned        
            align = new Dictionary<string, Dictionary<string, string>>();    
            foreach (var item in r.profiles)
            {
                align.Add(item.Key, new Dictionary<string, string>());

                if((alignFile==null || alignFile.Length==0))
                {
                    refSeq = "";
                    foreach (var pp in item.Value)
                    {
                        if (pp.Value.profile == null)
                            throw new Exception("For " + pp.Key + " no profile has been found!");

                        if (refSeq == null || pp.Value.sequence.Length > refSeq.Length)                        
                            refSeq = pp.Value.sequence;                                             
                    }                   
                    MAlignment al = new MAlignment(refSeq.Length);
                    foreach (var pp in item.Value)
                    {
                        if (pp.Value.sequence!=null && refSeq.Length != pp.Value.sequence.Length)
                        {
                            string seq = (al.Align(refSeq, pp.Value.sequence)).seq2;
                            if (seq.Length > 0)
                                align[item.Key][pp.Key] = seq;
                        }
                    
                    }                

                }
                if (align.Count == 0)
                    throw new Exception("There is no alignment");
            }
            try
            {
                if (dirSettings.mode != INPUTMODE.USER_DEFINED)
                    AlignProfiles();
                else
                    OmitAlignProfiles();
            }
            catch (Exception)
            {
                throw new Exception("Combining alignments went wrong!");
            }
			
		}
		
        private void AlignProfile(string protName,string profileName,List<byte>prof)
        {
            List<byte> ll = new List<byte>(prof.Count);
            int m = 0;
            if (align == null || align[profileName] == null || !align[profileName].ContainsKey(protName))
            {
                foreach (var item in prof)
                    ll.Add(item);
            }
            else
            {
                string alignProfile = align[profileName][protName];
                
                for (int i = 0; i < alignProfile.Length; i++)
                {
                    if (alignProfile[i] == '-')
                    {
                        ll.Add(0);
                        continue;
                    }
                    if (m < prof.Count)
                        ll.Add(prof[m++]);
                    else
                        ErrorBase.AddErrors("Profile " + profileName + " for " + protName + " seems to be incorect");

                }
            }
            protInfo<byte> tmp = new protInfo<byte>();
            if (r.profiles[profileName].ContainsKey(protName))
            {
                tmp = r.profiles[profileName][protName];
                tmp.alignment = ll;
                r.profiles[profileName][protName] = tmp;
            }
           

        }

        public void OmitAlignProfiles()
        {
            List<string> keys = new List<string>(r.profiles.Keys);
            r.protCombineStates = new Dictionary<string, List<byte>>();
            foreach (var item in r.profiles[keys[0]].Keys)
                r.protCombineStates.Add(item,r.profiles[keys[0]][item].profile);

            r.MakeDummyCombineCoding();
        }
		public void AlignProfiles()
		{
            List <string> profName=new List<string>(r.profiles.Keys);
            List<string> protNames = new List<string>(r.profiles[profName[0]].Keys);
            foreach (string protName in protNames)
            {
                if(align[profName[0]].ContainsKey(protName))
                    if (((string)align[profName[0]][protName]).Length < 5)
                        continue;
                bool test = true;
                foreach (var item in r.profiles)
                    if (!item.Value.ContainsKey(protName))
                        test = false;

                if (!test)
                    continue;


                foreach (var item in r.profiles)
                    AlignProfile(protName, item.Key, item.Value[protName].profile);                    
                
               
             
              r.CombineProfiles(protName, r.profiles);
             
               align.Remove(protName);
               gcCounter++;
               if (gcCounter > 5000)
               {
                   GC.Collect();
                   gcCounter = 0;
               }
            }
          
            GC.Collect();

           
		}
        public int ProfileLength()
        {
            return r.GetProfileLength();
        }
        public Dictionary<byte, double> StatesStat()
        {
            return r.StatesStat();
        }
        
		
	}
}

