using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace uQlustCore.Profiles
{
    public static class ProfileStat
    { 

        public static Dictionary<string, protInfo<byte>> TwoStates(Dictionary<string, protInfo<byte>> dic)
        {
            Dictionary<string, protInfo<byte>> res = new Dictionary<string, protInfo<byte>>();
            List<string> keys = new List<string>(dic.Keys);
            for (int j = 0; j < keys.Count; j++)
            {
                protInfo<byte> xx = dic[keys[j]];
                List<byte> newProfile = new List<byte>();
                for (int i = 0; i < xx.profile.Count; i++)
                {
                    if (xx.profile[i] > 0)
                        newProfile.Add(2);
                    else
                        newProfile.Add(1);
                }
                xx.profile = newProfile;
                res.Add(keys[j], xx);
            }

            return res;
        }
        public static Dictionary<string, protInfo<byte>> RearangeStates(Dictionary<string, protInfo<byte>> dic,double percent)
        {
            Dictionary<string, protInfo<byte>> res = new Dictionary<string, protInfo<byte>>();
            int[] freq = new int[byte.MaxValue];
            byte[] table = new byte[dic.Keys.Count];
            List<string> keys = new List<string>(dic.Keys);
            int[,,] interval = new int[dic[keys[0]].profile.Count,(int)Math.Floor(1.0 / percent)+1, 2];

            for (int j = 0; j < dic[keys[0]].profile.Count; j++)
            {
                for (int i = 0; i < freq.Length; i++)
                    freq[i] = 0;

                for (int i = 0; i < keys.Count; i++)
                {
                    table[i] = dic[keys[i]].profile[j];
                    freq[table[i]]++;
                }
                double sum = 0;
                int remIndex = 0;
                int intervCount = 0;
                for (int i = 0; i < freq.Length; i++)
                {
                    sum += freq[i];
                    if (sum/keys.Count>percent)
                    //if (i==0)//For two states only
                    {
                        sum = 0;
                        interval[j, intervCount, 0] = remIndex;
                        interval[j, intervCount, 1] = i+1;
                        intervCount++;
                        remIndex = i+1;
                    }
                }
                if(sum!=0)
                {
                    interval[j, intervCount, 0] = remIndex;
                    interval[j, intervCount, 1] = int.MaxValue;
                    intervCount++;
                }
                interval[j, intervCount - 1, 1] = int.MaxValue;
            }
            for(int j=0;j<keys.Count;j++)
            { 
                protInfo<byte> xx = dic[keys[j]];
                List<byte> newProfile = new List<byte>();
                for (int i = 0; i < xx.profile.Count; i++)
                    for (byte k = 0; k < interval.GetLength(1); k++)
                        if (xx.profile[i] >= interval[i,k, 0] && xx.profile[i] < interval[i,k, 1])
                        {
                            /*for (int s = 0; s < interval.GetLength(1); s++)
                                if (s == k)
                                    newProfile.Add(2);
                                else
                                    newProfile.Add(1);*/
                            
                                newProfile.Add((byte)(k + 1));
                            break;
                        }
                        else
                        {
                            Console.WriteLine("ups");
                        }
                xx.profile = newProfile;
                res.Add(keys[j],xx);
            }


            return res;
        }
        public static Dictionary<KeyValuePair<string, string>, string> ReadFile(string fileName)
        {
            Dictionary<KeyValuePair<string, string>, string> dicOrg = new Dictionary<KeyValuePair<string, string>, string>();
            StreamReader str = new StreamReader(fileName);
            Random r = new Random();
            bool test = true;
            string line = str.ReadLine();
            while (line != null)
            {
                line = line.Trim();
                line = line.Replace("  ", " ");
                string[] aux = line.Split(' ');
                if (aux.Length == 3)
                {
                    KeyValuePair<string, string> a = new KeyValuePair<string, string>(aux[1], aux[2]);
                    if (fileName.Contains("OUT"))
                        if (r.Next(100) <= 50)
                            test = true;
                        //else
                          //  test = false;
                            
                    if (!dicOrg.ContainsKey(a) && test)
                        dicOrg.Add(a, aux[0]);
                    //else
                        //Console.Write("Ups+" + aux[1] + " " + aux[2]);

                }
                line = str.ReadLine();
            }
            str.Close();

            return dicOrg;
        }
        public static double CalcDicDist(Dictionary<string, List<byte>> dic, Dictionary<KeyValuePair<string, string>, string> cl,List<int> bestIndex,int index)
        {
            double avr = 0;
            int count = 0;
            foreach (var item in cl)
            {
                double dist = 0;
                if (dic.ContainsKey(item.Key.Key) && dic.ContainsKey(item.Key.Value))
                {
                    List<byte> m1 = dic[item.Key.Key];
                    List<byte> m2 = dic[item.Key.Value];
                    for (int i = 0; i < bestIndex.Count; i++)
                        if (m1[bestIndex[i]] != m2[bestIndex[i]])
                            dist++;
                    if (m1[index] != m2[index])
                        dist++;

                    avr += dist;
                    count++;
                }
            }
            return avr/count;
        }
        public static Dictionary<string, protInfo<byte>> SelectFeatures(Dictionary<string, protInfo<byte>> dic,int featsNumber)
        {
            Dictionary<string, protInfo<byte>> resDic =null;
            List<string> dicKey = new List<string>(dic.Keys);
            int[] freq = new int[dic[dicKey[0]].profile.Count];
            int[] idx = new int[dic[dicKey[0]].profile.Count];
            for (int i = 0; i < idx.Length; i++)
                idx[i] = i;

            for (int j = 0; j < dic[dicKey[0]].profile.Count; j++)
            {
                foreach(var item in dicKey)
                {
                    freq[j] += dic[item].profile[j];
                }
            }

            Array.Sort<int>(idx, (a, b) => freq[b].CompareTo(freq[a]));

            resDic = new Dictionary<string, protInfo<byte>>();
            foreach(var item in dicKey)
            {
                List<byte> qq = new List<byte>();
                for (int i = 0; i < featsNumber; i++)
                    qq.Add(dic[item].profile[idx[i]]);

                protInfo<byte> aux = new protInfo<byte>();
                aux.profile = qq;
                aux.sequence = dic[item].sequence;
                aux.alignment = dic[item].alignment;

                resDic.Add(item, aux);
            }

            return resDic;
        }
        public static void  SelectFeatures(Dictionary<string, List<byte>> dic, string fileNameIN, string fileNameOUT)
        {
            Dictionary<KeyValuePair<string, string>, string> dicIN = ReadFile(fileNameIN);
            Dictionary<KeyValuePair<string, string>, string> dicOUT = ReadFile(fileNameOUT);

           List<int> bestIndex = new List<int>();
            int remIndex = -1;
            List<string> dicKey = new List<string>(dic.Keys);
            for (int k = 0; k < dic[dicKey[0]].Count; k+=4)
            {
                double remRes = double.MinValue;
                for (int j = k; j < dic[dicKey[0]].Count; j+=4)
                {
                    if (bestIndex.Contains(j))
                        continue;
                    double avrIn, avrOut;
                    avrIn = CalcDicDist(dic, dicIN, bestIndex, j);
                    avrOut = CalcDicDist(dic, dicOUT, bestIndex, j);
                    double res = avrOut - avrIn;
                    if (res > remRes)
                    {
                        remRes = res;
                        remIndex = j;
                    }
                }
                for(int j=0;j<4;j++)
                    bestIndex.Add(remIndex+j);
            }

        }
        public static void CalcHammingDist(Dictionary<string, List<byte>> dic, string fileName,string fileNameOut)
        {
            StreamReader str = new StreamReader(fileName);
            StreamWriter outStr = new StreamWriter(fileNameOut);
            List<KeyValuePair<string,string>> xx = new List<KeyValuePair<string, string>>();
            Dictionary<KeyValuePair<string, string>, string> dicOrg = new Dictionary<KeyValuePair<string, string>, string>();
            string line = str.ReadLine();
            while (line != null)
            {
                line = line.Trim();
                line = line.Replace("  ", " ");
                string[] aux = line.Split(' ');
                if (aux.Length >= 2)
                {
                    if (dic.ContainsKey(aux[aux.Length-1]) && dic.ContainsKey(aux[aux.Length-2]))
                    {
                        int count = 0;
                        int val = 0, key = 0;
                        for (int i = 0; i < dic[aux[aux.Length-2]].Count; i++)
                        {
                            key += dic[aux[aux.Length-2]][i];
                            val += dic[aux[aux.Length-1]][i];
                            if (dic[aux[aux.Length-2]][i] != dic[aux[aux.Length-1]][i])
                                count++;
                        }
                        outStr.WriteLine(count + " " + aux[aux.Length-2] + " " + aux[aux.Length-1]);
                    }
                }
                line = str.ReadLine();
            }
            str.Close();            
            outStr.Close();
        }
    }
}
