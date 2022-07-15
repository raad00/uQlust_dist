using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Reflection;
using System.Linq.Expressions;
using System.Timers;
using uQlustCore;
using uQlustCore.dcd;
using uQlustCore.Distance;
using uQlustCore.Profiles;


namespace uQlustTerminal
{
    class Program
    {
        static JobManager manager = new JobManager();
       // static Timer t = new Timer();
        private static void UpdateProgress(object sender, EventArgs e)
        {
           
            Dictionary<string, double> res = manager.ProgressUpdate();

            if (res == null)
            {
//                TimeInterval.Stop();
                Console.Write("\r                                                             ");
                return;
            }
            Console.Write("\r                                                                              ");
            foreach (var item in res)
                Console.Write("\rProgress " + item.Key + " " + (item.Value*100).ToString("0.00")+"%");

        }

        public static void ErrorMessage(string message)
        {
            Console.WriteLine(message);
        }
        public static Options SetOptions(string minState,string testDir,string resFile,int retrivalSize,bool creatBase,bool mSim)
        {
            Options opt = new Options();
            opt.dataDir.Clear();
            opt.clusterAlgorithm.Clear();
            opt.profileFiles.Clear();

            opt.hash.rpart = false;
            opt.hash.selectionMethod = COL_SELECTION.ENTROPY;

            //opt.clusterAlgorithm.Add(ClusterAlgorithm.HNNBinary);
            if (creatBase)
            {
                opt.dataDir.Add(testDir);
                opt.clusterAlgorithm.Add(ClusterAlgorithm.HNNDataBase);
            }
            else
            {
                opt.profileFiles.Add("gauss_frag.profiles");
                opt.clusterAlgorithm.Add(ClusterAlgorithm.HNNBinary);
            }
            
            opt.hnn.testFile = testDir;
//            opt.hnn.labelsFile = textBox3.Text;
  
            
            opt.hash.relClusters = 100000000;
            opt.hash.perData = 100;
            opt.hash.useConsensusStates = false;
            opt.outputFile = resFile;
            opt.hnn.binaryfile = minState;
            opt.hnn.mSim = mSim;
            opt.hash.profileName = "profiles/gauss_frag.profiles";
            opt.hnn.retrivalSize = retrivalSize;
            if(!File.Exists(opt.hash.profileName))
            {
                Console.WriteLine("File gauss_frag.profiles cannot be found in the directory:profiles");
                System.Environment.Exit(-1);

            }
            return opt;
        }
        static void Main(string[] args)
        {
            bool errors = false;
            bool times = false;
            bool progress = false;
            bool createBase = false;
            bool mSim = true;
            string minStateFileName = "";
            string testDirectory = "";
            string resultFile = "";
            int retrivalSize = 20;

            Options opt = new Options();
            try
            {
                InternalProfilesManager.InitProfiles();
            }
            catch(Exception ex)
            {
                Console.WriteLine("Some of the profiles are not available: ", ex.Message);
            }
            //Console.WriteLine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location));
            /*foreach(var item in InternalProfilesManager.internalList)
            {
                Console.WriteLine("profile=" + item);
            }*/
            if (args.Length == 0)
            {
                Console.WriteLine("To create database use following options:");
                Console.WriteLine("-c file name - save database to the file");
                Console.WriteLine("-d directory or file name with list of pdb files from which database will be created");
                Console.WriteLine("-mSim  or -mLSH-c5 use method for database preparing default mSim");
                Console.WriteLine("To find the closest structures use following arguments");
                Console.WriteLine("-mSim or -mLSH-c5 use method for distance calcultion default mSim");
                Console.WriteLine("-f database file name (Save name with -c option, Load otherwise");
                Console.WriteLine("-d directory with structures for testing or file name with list of files");
                Console.WriteLine("-s result file name");
                Console.WriteLine("Following parameters are optional:");
                Console.WriteLine("-e \n\t show all errors");
                Console.WriteLine("-n \n\t set number of cores to be used");
                Console.WriteLine("-t \n\tshow time information");
                Console.WriteLine("-p \n\tShow progres bar");
                Console.WriteLine("-r   \n\tretrival size [default 20]");
                return;
            }
            Settings set = new Settings();
            set.Load();
            if (set.profilesDir == null || set.profilesDir.Length == 0)
            {
                set.profilesDir = "generatedProfiles";
            //    set.Save();
            }
            set.mode = INPUTMODE.PROTEIN;
            for (int i = 0; i < args.Length; i++)
            {
                switch (args[i])
                {
                    case "-mSim":
                        mSim = true;
                        break;
                    case "-mLSH-c5":
                        mSim = false;
                        break;
                    case "-f":
                        if(i+1>=args.Length)
                        {
                            Console.WriteLine("After -f option you have to provide configuration file");
                            return;
                        }
                        if (!File.Exists(args[i + 1]))
                        {
                            Console.WriteLine("File " + args[i + 1] + " does not exist");
                            return;
                        }
                        minStateFileName = args[i + 1];
                        i++;
                        break;
                    case "-c":
                        if (i + 1 >= args.Length)
                        {
                            Console.WriteLine("After -c file name for database is need");
                            return;
                        }
                        if (File.Exists(args[i + 1]))
                        {
                            Console.WriteLine("File "+args[i+1]+" already exists! Provide non existing file name");
                            return;
                        }
                        createBase = true;
                        minStateFileName = args[i + 1];
                        break;

                    case "-d":
                        if (i + 1 >= args.Length)
                        {
                            Console.WriteLine("After -f option you have to provide directory name or file nam with list of files");
                            return;
                        }
                        if (!Directory.Exists(args[i + 1]) && !File.Exists(args[i+1]))
                        {
                            Console.WriteLine("Directory " + args[i + 1] + " does not exist");
                            return;
                        }
                        testDirectory = args[i + 1];
                        break;
                    case "-s":
                        if (i + 1 >= args.Length)
                        {
                            Console.WriteLine("After -s option you have to provide results file name");
                            return;
                        }                       
                        resultFile = args[i + 1];
                        i++;
                        break;
                    case "-r":
                        try
                        {
                            retrivalSize = Convert.ToInt32(args[++i]);
                        }
                        catch (Exception ex)
                        {
                            Console.WriteLine("Wrong value of retrival size: " + ex.Message);
                            return;
                        }
                        break;
                    case "-n":
                        if (args.Length > i)
                        {
                                        
                            int num;
                            try
                            {
                                num = Convert.ToInt32(args[++i]);                                
                                set.numberOfCores = num;
                            }
                            catch(Exception ex)
                            {
                                Console.WriteLine("Wrong definition of number of cores: " + ex.Message);
                                return;
                            }
                        }
                        else
                            Console.WriteLine("Number of cores has been not provided");
                        break;
                    case "-e":
                        errors = true;
                        break;
                    case "-t":
                        times = true;
                        break;
                    case "-p":
                        progress = true;
                        break;
                    default:
                        if(args[i].Contains("-"))
                            Console.WriteLine("Unknown option " + args[i]);
                        break;

                }
            }
            set.Save();
            if (!createBase)
            {
                if (minStateFileName.Length == 0 || testDirectory.Length == 0 || resultFile.Length == 0)
                {
                    Console.WriteLine("One of the required file was not provided!");
                    return;
                }

            }            
            else
                if (minStateFileName.Length == 0 || testDirectory.Length == 0)
            {
                Console.WriteLine("One of the required file was not provided!");
                return;
            }
            string[] aux = null;
            try
            {
                Console.WriteLine("minState file " + minStateFileName);
                Console.WriteLine("test directory " + testDirectory);
                Console.WriteLine("resultFile " + resultFile);
                Console.WriteLine("retrivalSize " + retrivalSize);
                opt=SetOptions(minStateFileName,testDirectory,resultFile,retrivalSize,createBase,mSim);
                
                
                aux = args[0].Split('.');
                manager.opt = opt;
                manager.message = ErrorMessage;
                if (progress)
                {
                    TimeIntervalTerminal.InitTimer(UpdateProgress);
                    TimeIntervalTerminal.Start();
                }
                manager.RunJob("");
                manager.WaitAllNotFinished();
                
                UpdateProgress(null,null);
                
                if(progress)
                    TimeIntervalTerminal.Stop();
                Console.Write("\r                                                                     ");
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception : " + ex.Message);
            }
            if (!createBase && manager.clOutput.Count > 0)
            {
                foreach (var item in manager.clOutput.Keys)
                {
                    Dictionary<string,string> res = manager.clOutput[item].hNNRes;
                    ClusterOutput.SaveHnn(res, opt.outputFile);
                    //clusterOut.SCluster(clustName+"_"+opt.outputFile);
                }
            }
            if (times)
            {
                foreach (var item in manager.clOutput)
                    Console.WriteLine(item.Value.dirName + " " + item.Value.measure + " " + item.Value.time);
            }
            if (errors)
            {
                foreach (var item in ErrorBase.GetErrors())
                    Console.WriteLine(item);
            }
            Console.WriteLine();
        }
    }
}
