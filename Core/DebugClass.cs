using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace uQlustCore
{
    public static class DebugClass
    {
        static bool DEBUG = false;
        static string fileName = "logFile.dat";
        static StreamWriter logFile;

        public static void WriteMessage(string message)
        {
            if (DEBUG)
            {
                logFile.WriteLine(message);
                logFile.Flush();
            }
        }
        public static void WriteMessageInLine(string message)
        {
            if (DEBUG)
                logFile.Write(message);
        }

        public static void DebugOn()
        {
            try
            {
                logFile = new StreamWriter(fileName);
            }
            catch(IOException ex)
            {
                logFile.Close();
                logFile = new StreamWriter(fileName);

            }
            DEBUG = true;
        }
        public static void DebugOff()
        {
            DEBUG = false;
            logFile.Close();
        }
    }
}
