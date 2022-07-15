using System;
using System.Collections;
using System.IO;
using System.Collections.Specialized;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using uQlustCore;
using System.Data;
using System.Threading;
using uQlustCore.Interface;
namespace uQlustCore
{
	public class jury1D:IProgressBar
	{
        string currentProfile="";
        //public static string defaultProfileFile = "defaultJury1d.profile";

        Dictionary<string, List<byte>> stateAlign;
                List<KeyValuePair<string, double>> res1Djury;
        public List<string> alignKeys
        {
            get
            {
                if (stateAlign != null && stateAlign.Keys != null)
                    return new List<string>(stateAlign.Keys);
                else
                    return null;
            }
        }
        Alignment al;
       
        Dictionary<byte, Dictionary<byte,double>> weights;
        public Dictionary<byte, int>[] columns = null;
        List<string> allStructures;
        ManualResetEvent[] resetEvents = null;
        long maxV, currentV;
        Settings set = new Settings();
        int threadNumbers;

        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }

        public jury1D()
        {
            maxV = 1;
            currentV = 0;
            al = new Alignment();
            set.Load();
            threadNumbers = set.numberOfCores;
        }
        public jury1D(Dictionary<byte, int>[] freq1,Alignment al)
        {
            columns = freq1;          
            this.al = al;
        }
        public double ProgressUpdate()
        {
            double progress = al.ProgressUpdate();
            return StartProgress+(EndProgress - StartProgress) * (progress * 0.75 + 0.25 * ((double)currentV / maxV));
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }

        public void PrepareJury(DCDFile dcd, string alignFile, string profile = null)
        {
            currentProfile = profile;
            al.Prepare(dcd, set, currentProfile);
            CommonConstr(alignFile, profile);
        }

        public void PrepareJury(string dirName,string alignFile,string profile=null)
        {
            if (alignFile == null)
            {
                currentProfile = profile;
                al.Prepare(dirName, set, currentProfile);


                CommonConstr(alignFile, profile);
                if (al.GetAlignment().Count == 0)                                
                    return;               
            }
            else
                JuryProfileFile(alignFile, profile);

        }
        public void PrepareJury(List<string> fileNames, string alignFile,string profile=null)
        {
            currentProfile = profile;

            if (alignFile == null)
                al.Prepare(fileNames, set, currentProfile);
            else
                al.Prepare(alignFile, profile);

            CommonConstr(alignFile, profile);

        }
        public void PrepareJury(List<KeyValuePair<string,string>> profiles,string profName, string profileFile)
        {
            al.Prepare(profiles, profName, profileFile);
            stateAlign = al.GetStateAlign();
            weights = al.r.GenerateWeights(wOpertion.MULT);
           // maxProgress = 0.75;
        }
        //all profiles in profileFile should be aligned
        public void PrepareJury(string profilesFile, string profile)
        {
            JuryProfileFile(profilesFile, profile);
            //weights = al.r.GenerateWeights(wOpertion.MULT);
            //maxProgress = 0.75;
        }
        private void JuryProfileFile(string profilesFile, string profile)
        {
            al.Prepare(profilesFile, profile);
            CommonConstr(profilesFile, profile);        
            /*stateAlign = new Dictionary<string, List<byte>>();
            foreach (var itemK in al.r.profiles)
            {
                Dictionary<string, protInfo> aux = itemK.Value;
                foreach (var item in aux)  
                    stateAlign.Add(item.Key,item.Value.profile);
                
                break;
            }*/
            //CommonConstr(profilesFile, profile);        
        }
        private void CommonConstr(string alignFile, string profile)
        {
            al.MyAlign(alignFile);
            stateAlign = al.GetStateAlign();
            //AddErrors(al.errors);
            weights = al.r.GenerateWeights(wOpertion.MULT);
            
          
            if (weights.Keys.Count == 0)
                throw new Exception("Weights in profile: " + profile + " have not been defined. 1DJury cannot work!");


            
            //maxProgress = 0.75;

        }
        public void PrepareJury(Alignment al,Dictionary<string,List<byte>> profDic)
        {
            this.al = al;
            stateAlign = profDic;
            weights = al.r.GenerateWeights(wOpertion.MULT);
        }
        public void PrepareJury(Alignment al)
        {
            this.al = al;
            stateAlign = al.GetStateAlign();

            weights = al.r.GenerateWeights(wOpertion.MULT);

        }
        public override string ToString()
        {
            return "1DJury";
        }
        public Dictionary<byte, List<int>>[] MakeColumnsLists(List<string> structNames)
        {
            byte locState=0;
            Dictionary<byte, List<int>>[] columns = null;

            if (structNames.Count == 0)
                return null;

            columns = new Dictionary<byte, List<int>>[stateAlign[structNames[0]].Count];
            for (int i = 0; i < columns.Length; i++)
            {

                
                Dictionary<byte, List<int>> dicCol = new Dictionary<byte, List<int>>();
                for (int j = 0; j < structNames.Count; j++)
                {
                    if (stateAlign.ContainsKey(structNames[j]) && stateAlign[structNames[j]].Count > 0)
                        locState = stateAlign[structNames[j]][i];
                    else
                        continue;
                    if (locState == 0)
                        continue;

                    if (!dicCol.ContainsKey(locState))
                    {
                        List<int> lista = new List<int>();
                        lista.Add(j);
                        dicCol.Add(locState, lista);
                    }
                    else
                        dicCol[locState].Add(j);
                }

                columns[i] = dicCol;
            }
            
            return columns;
        }
        private void ThreadingMakeColumns(object o)
        {
            byte locState;
            ThreadParam pp = (ThreadParam)o;
            int threadNum = pp.num;
            int start = pp.start;
            int stop = pp.stop;

            for (int i = start; i < stop; i++)
            {
                lock (columns[i])
                {
                   // int counter = 0;
                    foreach (var name in allStructures)
                    {
                        List<byte> cc = stateAlign[name];

                        locState = cc[i];
                        if (locState == 0)
                            continue;
                        if (!columns[i].ContainsKey(locState))
                        {
                            columns[i].Add(locState, 1);
                            continue;
                        }
                        
                        //columns[i].AddOrUpdate(locState,0, (key, value) => value + 1);
                        columns[i][locState]++;//, 0, (key, value) => value + 1);
                    }

                }                
                Interlocked.Increment(ref currentV);
            }
            resetEvents[threadNum].Set();
        }
        public Dictionary<byte, int>[] MakeColumns(List<string> structNames)
        {           

            if (structNames.Count == 0)
                return null;

            maxV = stateAlign[structNames[0]].Count + stateAlign.Count;

            columns = new Dictionary<byte, int>[stateAlign[structNames[0]].Count];
            allStructures = structNames;

            resetEvents = new ManualResetEvent[threadNumbers];

            for (int i = 0; i < columns.Length; i++)
            {
                columns[i] = new Dictionary<byte, int>();
           //     foreach (var item in weights.Keys)
             //       columns[i].Add(item,0);
            }
            for (int n = 0; n < threadNumbers; n++)
            {
                ThreadParam pp = new ThreadParam();
                pp.num = n;
                pp.start = (int)(n * columns.Length / Convert.ToDouble(threadNumbers));
                pp.stop=(int)( (n + 1) * columns.Length / Convert.ToDouble(threadNumbers));
                //int p = n;
                //int start = (int)(n * columns.Length /Convert.ToDouble(threadNumbers));
                //int stop =(int)( (n + 1) * columns.Length / Convert.ToDouble(threadNumbers));
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(ThreadingMakeColumns), (object) pp);
                

            }
            for (int n = 0; n < threadNumbers; n++)
                resetEvents[n].WaitOne();
            return columns;
        }
        public List<byte> GetStructureStates(string structName)
        {

            if(stateAlign.ContainsKey(structName))
                return stateAlign[structName];

            return null;
        }
        private void ThreadingScoreCalc(object o)
        {
            double score;
            byte locState;

            //object[] array = o as object[];
            //int threadNum = (int)array[0];
            //int start = (int)array[1];
            //int stop = (int)array[2];
            ThreadParam k = (ThreadParam)o;
            int threadNum = k.num;
            int start = k.start;
            int stop = k.stop;
            for (int j = start; j < stop;j++ )
            {
                score = 0;
                string name = allStructures[j];
                List<byte> listStates = stateAlign[name];
                for (int i = 0; i < listStates.Count; i++)
                {
                    locState = listStates[i];
                    if (locState == 0)
                        continue;
                    if (weights.ContainsKey(locState))
                    {
                        Dictionary<byte, double> auxDic = weights[locState];
                        foreach (var item in columns[i])
                            if (auxDic.ContainsKey(item.Key))
                                score += auxDic[item.Key] * item.Value;
                    }
                }
                score /= allStructures.Count * listStates.Count;

                KeyValuePair<string, double> v = new KeyValuePair<string, double>(name, score);
                lock (res1Djury)
                {
                    res1Djury.Add(v);
                }
                Interlocked.Increment(ref currentV);
            }
            resetEvents[threadNum].Set();
            
        }
        public ClusterOutput JuryOptWeights(List<string> structNames,Dictionary<byte,int>[]locColumns=null)
        {
            res1Djury = new List<KeyValuePair<string, double>>(structNames.Count);
            Dictionary<byte, int>[] columns = null;

            List<string> aux = new List<string>(structNames);
            foreach (var item in structNames)
            {
                if (!stateAlign.ContainsKey(item) || stateAlign[item].Count == 0)
                    aux.Remove(item);
            }
            if (locColumns == null)
                columns = MakeColumns(aux);
            else
            {
                columns = locColumns;
                this.columns = locColumns;
            }

            
            currentV++;
            if (columns == null)
                return null;
            
            allStructures = new List<string>(aux);
            
            if (weights == null || weights.Count==0)
                return null;

            

            resetEvents = new ManualResetEvent[threadNumbers];

            /*weights["0"]["0"] = 1;
            weights["1"]["1"] = columns.Length / 10;//Do usuniecia!!!!!!!!!!!!!!!
            weights["1"]["0"] = 0;
            weights["0"]["1"] = 0;*/



            for (int n = 0; n < threadNumbers; n++)
            {
                ThreadParam pp = new ThreadParam();
                //int p = n;
                //int start = (int)(n * allStructures.Count / Convert.ToDouble(threadNumbers));
                //int stop = (int)((n + 1) * allStructures.Count / Convert.ToDouble(threadNumbers));
                pp.num = n;
                pp.start = (int)(n * allStructures.Count / Convert.ToDouble(threadNumbers));
                pp.stop = (int)((n + 1) * allStructures.Count / Convert.ToDouble(threadNumbers));
                resetEvents[n] = new ManualResetEvent(false);
                //ThreadPool.QueueUserWorkItem(new WaitCallback(ThreadingScoreCalc), new object[] { p, start, stop });
                ThreadPool.QueueUserWorkItem(new WaitCallback(ThreadingScoreCalc), (object)pp);

            }
            for (int n = 0; n < threadNumbers; n++)
                resetEvents[n].WaitOne();
            
            currentV++;
            res1Djury.Sort((firstPair, nextPair) =>
            {
                return nextPair.Value.CompareTo(firstPair.Value);
            });

            ClusterOutput juryRes = new ClusterOutput();
            juryRes.runParameters = "Profile: " + this.currentProfile;
            juryRes.juryLike = res1Djury;
            currentV = maxV;
            
            return juryRes;
			

        }
        public ClusterOutput ConsensusJury(List<string> structNames)
        {
            List<KeyValuePair<string, double>> distCons = new List<KeyValuePair<string, double>>();
            byte locState;
            List<Dictionary<byte, int>> cons = new List<Dictionary<byte,int>>(structNames.Count);
            List<byte> finalCons = new List<byte>();

            foreach (string name in structNames)
            {
                for (int i = 0; i < stateAlign[name].Count; i++)
                     cons.Add(new Dictionary<byte, int>());
                break;
            }

            foreach (string name in structNames)
            {
                if (!stateAlign.ContainsKey(name))
                    continue;
                for (int i = 0; i < stateAlign[name].Count; i++)
                {
                    locState = stateAlign[name][i];
                    if (cons[i].ContainsKey(locState))
                        cons[i][locState]++;
                    else
                        cons[i].Add(locState, 1);
                }

            }
            foreach (var item in cons)
            {
                var items = from pair in item
                            orderby pair.Value descending
                            select pair;
                // Display results.
                foreach (KeyValuePair<byte, int> pair in items)
                {
                    finalCons.Add(pair.Key);
                    break;
                }

            }
            foreach (string name in structNames)
            {
                if (!stateAlign.ContainsKey(name))
                    continue;
                float dist = 0;
                for (int i = 0; i < stateAlign[name].Count; i++)
                {
                    if (stateAlign[name][i] == finalCons[i])
                        dist++;
                }
                KeyValuePair<string, double> v = new KeyValuePair<string, double>(name, dist);
                distCons.Add(v);
            }
            distCons.Sort((firstPair, nextPair) =>
            {
                return firstPair.Value.CompareTo(nextPair.Value);
            });

            ClusterOutput juryRes = new ClusterOutput();

            juryRes.juryLike = distCons;
            return juryRes;

        }
/*        [Obsolete("Instead of this method use JuryOptWeights, this method has fixed weights values")]
		public ClusterOutput JuryOpt(List <string> structNames)
		{
			string locState;
			List<KeyValuePair<string,double>> res1Djury=new List<KeyValuePair<string,double>>();
			Dictionary <string, float> score= new Dictionary<string, float>();
            Dictionary<string, int>[] columns = null;


            columns = MakeColumns(structNames);

            if (columns == null)
                return null;

			foreach (string name in structNames)
				score[name]=0;
			
			
			foreach (string name in structNames)
			{
                for (int i = 0; i < stateAlign[name].Count; i++)
				{
                    locState = stateAlign[name][i];
					if(locState=="-")
						continue;
					if(locState.Contains("H") || locState.Contains("E"))
						score[name]+=weightHE*columns[i][locState];		
                    else
                        score[name] += weightC * columns[i][locState];		
					
					if(locState.Contains("H") || locState.Contains("E"))
					{
						//weight=0.5f;
						locState=locState.Replace("H","C");
						locState=locState.Replace("E","C");
						if(columns[i].ContainsKey(locState))
							score[name]+=weightC*columns[i][locState];														
					}
                    if (stateAlign[name][i].Contains("C"))
				    {
                        locState = stateAlign[name][i];					
						locState=locState.Replace("C","H");
						if(columns[i].ContainsKey(locState))
							if(score.ContainsKey(name))
								score[name]+=weightC*columns[i][locState];		
							else 
								score[name]=weightC*columns[i][locState];

                        locState = stateAlign[name][i];					
						locState=locState.Replace("C","E");
						if(columns[i].ContainsKey(locState))
							if(score.ContainsKey(name))
								score[name]+=weightC*columns[i][locState];		
							else 
								score[name]=weightC*columns[i][locState];
							
					}
				}
                score[name] /= structNames.Count * stateAlign[name].Count;
				
				KeyValuePair<string,double> v=new KeyValuePair<string, double>(name,score[name]);				
				res1Djury.Add(v);
			}			
			res1Djury.Sort((firstPair,nextPair) =>
    		{
        		return nextPair.Value.CompareTo(firstPair.Value);
    		});

			//Console.WriteLine("Next");
			//foreach (KeyValuePair<string,double> s in res1Djury)
			//	Console.WriteLine("Result "+s.Value+" VALUE="+s.Key);
            ClusterOutput juryRes = new ClusterOutput();

            juryRes.juryLike = res1Djury;
			return juryRes;
			
		}	*/
							
	}
	
}

