using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Runtime.InteropServices;
using uQlustCore;

namespace uQlustCore.dcd
{

    public class Header
    {
        public int control;
        public char[] nn = new char[4];
        public int coordNumber;
        public int simNumber;
        public int freq;
        public int atomMoveNumber;
        public float timeStep;
        public int blockSize;
        public int sizeTitle;
        
    }
    public class DCDReader
    {
        float[] posX;
        float[] posY;
        float[] posZ;
        string[] pdbLines;
        Header h;
        int atoms;
        BinaryReader reader=null;
        DCDFile dcdFile = null;
        FileStream fs=null;

        public DCDReader(DCDFile dcdFile)
        {
            this.dcdFile = dcdFile;
        }
        public override string ToString()
        {
            return " DCD";
        }
        public void DCDPrepareReading(string fileName, string pdbFile)
        {
            char[] tab = new char[80];        
            fs = new FileStream(fileName, FileMode.Open);


            //int count = Marshal.SizeOf(typeof(Header));

            int count = 100;
            byte[] readBuffer = new byte[count];
            reader = new BinaryReader(fs);
            count = 100;
            readBuffer = reader.ReadBytes(count);

            byte []test=new byte [4];
            Array.Copy(readBuffer,0,test,0,4);

            if (BitConverter.IsLittleEndian)
                Array.Reverse(test);

            h = new Header();
            int m = BitConverter.ToInt32(test, 0);
            h.control = BitConverter.ToInt32(test, 0);
            Array.Copy(readBuffer, 4, test, 0, 4);
            h.nn = System.Text.Encoding.UTF8.GetString(test).ToCharArray();
            Array.Copy(readBuffer, 8, test, 0, 4);
            if (BitConverter.IsLittleEndian)
                Array.Reverse(test);
            h.coordNumber = BitConverter.ToInt32(test, 0);
            Array.Copy(readBuffer, 12, test, 0, 4);
            if (BitConverter.IsLittleEndian)
                Array.Reverse(test);
            h.simNumber = BitConverter.ToInt32(test, 0);
            Array.Copy(readBuffer, 16, test, 0, 4);
            if (BitConverter.IsLittleEndian)
                Array.Reverse(test);
            h.freq = BitConverter.ToInt32(test, 0);
            Array.Copy(readBuffer, 40, test, 0, 4);
            if (BitConverter.IsLittleEndian)
                Array.Reverse(test);
            h.atomMoveNumber = BitConverter.ToInt32(test, 0);
            Array.Copy(readBuffer, 44, test, 0, 4);
            if (BitConverter.IsLittleEndian)
                Array.Reverse(test);
            h.timeStep= BitConverter.ToSingle(test, 0);
            Array.Copy(readBuffer, 92, test, 0, 4);
            if (BitConverter.IsLittleEndian)
                Array.Reverse(test);
            h.blockSize = BitConverter.ToInt32(test, 0);
            Array.Copy(readBuffer, 96, test, 0, 4);
            if (BitConverter.IsLittleEndian)
                Array.Reverse(test);
            h.sizeTitle = BitConverter.ToInt32(test, 0);

            int controlNumber;
            for (int i = 0; i < h.sizeTitle; i++)
                readBuffer=reader.ReadBytes(80);
            //Teraz 4 inty
            readBuffer = reader.ReadBytes(4);
            Array.Copy(readBuffer, 0, test, 0, 4);
            if (BitConverter.IsLittleEndian)
                Array.Reverse(test);
             controlNumber= BitConverter.ToInt32(test, 0);
             readBuffer = reader.ReadBytes(4);
             if (BitConverter.IsLittleEndian)
                 Array.Reverse(test);
             controlNumber = BitConverter.ToInt32(test, 0);

             readBuffer = reader.ReadBytes(4);
             if (BitConverter.IsLittleEndian)
                 Array.Reverse(test);
             atoms = BitConverter.ToInt32(test, 0);

            readBuffer = reader.ReadBytes(8);

            posX = new float[atoms];
            posY = new float[atoms];
            posZ = new float[atoms];

            ReadPDBFile(pdbFile);
        }
        public bool ReadAndSavePDB(StreamWriter sw)
        {
            string rep;
            int countAtom = 0;

            ReadFrame(reader, atoms);

            foreach (var line in pdbLines)
            {
                rep = line;
                if (line.Contains("ATOM"))
                {
                    rep = line.Substring(0, 31) + String.Format("{0,7}", posX[countAtom]).Substring(0, 7) + " " + String.Format("{0,7}", posY[countAtom]).Substring(0, 7) + " " + String.Format("{0,7}", posZ[countAtom]).Substring(0, 7) +
                        line.Substring(54, line.Length - 54);
                    countAtom++;
                    if (countAtom >= posX.Length)
                        throw new Exception("Number of atoms in pdb file is bigger than in dcd file");
                }
                sw.WriteLine(rep);

            }

            sw.Flush();
            h.coordNumber--;

            if (h.coordNumber > 0)
                return true;

            return false;
        }
        public void FinishDCDReading()
        {
            reader.Close();
            fs.Close();
        }
        public void DCDReadFile(string fileName,string pdbFile,string outDir)
        {
            DCDPrepareReading(fileName,pdbFile);
            
            for (int i = 0; i < h.coordNumber; i++)
                ReadAndSavePdbToFile(outDir+"\\"+"test_"+i+".pdb");

            FinishDCDReading();
            
        }
        void ReadFrame(BinaryReader reader,int atoms)
        {
            byte[] buff = new byte[4];
            buff = reader.ReadBytes(4);
            buff = reader.ReadBytes(4);
            int dSize;// = reader.ReadInt32();
            //dSize = reader.ReadInt32();
            for (int i = 0; i < atoms; i++)
            {
                buff = reader.ReadBytes(4);
                if (BitConverter.IsLittleEndian)
                    Array.Reverse(buff);
                posX[i]=BitConverter.ToSingle(buff, 0);
                //posX[i] = reader.ReadSingle();
            }
            buff = reader.ReadBytes(4);
            if (BitConverter.IsLittleEndian)
                Array.Reverse(buff);
            dSize = BitConverter.ToInt32(buff, 0);
//            dSize = reader.ReadInt32();
            dSize++;
            buff = reader.ReadBytes(4);
            //dSize = reader.ReadInt32();
            for (int i = 0; i < atoms; i++)
            {
                buff = reader.ReadBytes(4);
                if (BitConverter.IsLittleEndian)
                    Array.Reverse(buff);
                posY[i] = BitConverter.ToSingle(buff, 0);

//                posY[i] = reader.ReadSingle();
            }
            buff = reader.ReadBytes(4);
            buff = reader.ReadBytes(4);
            //dSize = reader.ReadInt32();
            //dSize = reader.ReadInt32();
            for (int i = 0; i < atoms; i++)
            {
                buff = reader.ReadBytes(4);
                if (BitConverter.IsLittleEndian)
                    Array.Reverse(buff);
                posZ[i] = BitConverter.ToSingle(buff, 0);
            }
             //   posZ[i] = reader.ReadSingle();
            buff = reader.ReadBytes(4);
            //dSize = reader.ReadInt32();
     
        }
        void ReadPDBFile(string pdbFile)
        {
             String line;
            using (StreamReader sr = new StreamReader(pdbFile))
            {
                    line = sr.ReadToEnd();            
            }
            if (line != null && line.Length > 0)
            {
                pdbLines = line.Split('\n');
            }
        }

        void ReadAndSavePdbToFile(string outFile)
        {            
            StreamWriter fileOut = new StreamWriter(outFile);
            ReadAndSavePDB(fileOut);
            fileOut.Close();
        }
    }
    
}
