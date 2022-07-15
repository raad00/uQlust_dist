
using System;
using System.Collections.Generic;
//using System.Windows.Media;
//using System.Windows.Media.Media3D;
using System.Globalization;
using uQlustCore;

namespace uQlustCore.PDB
{
	public class Point3D
	{
		private double x,y,z;
		
        public float X
        {
            set
            {
                x = value;
            }
            get
            {
                return (float)x;
            }
        }
        public float Y
        {
            set
            {
                y =value ;
            }
            get
            {
                return  (float)y;
            }
        }
        public float Z
        {
            set
            {
                z = value ;
            }
            get
            {
                return (float) z;
            }
        }

		public Point3D(float xp,float yp, float zp)
		{
			X=xp;Y=yp;Z=zp;
		}
        public Point3D()
        {
            x = y = z = 0;
        }
        public static Point3D operator +(Point3D p1, Point3D p2)
        {
            return new Point3D(p1.X + p2.X, p1.Y + p2.Y, p1.Z + p2.Z);
        }
        public static Point3D operator -(Point3D p1, Point3D p2)
        {
            return new Point3D(p1.X - p2.X, p1.Y - p2.Y, p1.Z - p2.Z);
        }
        public static double operator *(Point3D p1, Point3D p2)
        {
            return p1.x * p2.x+ p1.y * p2.y+ p1.z * p2.z;
        }
        public static Point3D operator *(double x, Point3D p)
        {
            return new Point3D((float)(p.X * x), (float)(p.Y * x), (float)(p.Z * x));
        }
        public static Point3D operator /(Point3D p, double x)
        {
            return new Point3D((float)(p.X / x), (float)(p.Y / x), (float)(p.Z / x));
        }
        public double Norm()
        {
            return Math.Sqrt(X * X + Y * Y + Z * Z);
        }
        public Point3D Normalize()
        {
            Point3D res = new Point3D();
            res.X = X;
            res.Y = Y;
            res.Z = Z;
            double len = Math.Sqrt(X * X + Y * Y + Z * Z);

            res.X /= (float)len;
            res.Y /= (float)len;
            res.Z /= (float)len;

            return res;
        }
    }
    public class Atom
    {
        private static Dictionary<byte, string> AtomNamesToBytes = new Dictionary<byte, string>();
        private static Dictionary<string, byte> AtomNames = new Dictionary<string, byte>();

      // private byte atomIndex;

        public string AtomName;
       
        public short[] tabParam = new short[3];
        public Point3D Position;
        protected virtual bool CheckResidue(string residueName)
        {
            if (Residue.IsAminoName(residueName))
                return true;
            
            return false;
        }
        protected virtual char ResidueIdentifier(string residueName)
        {
            return Residue.GetResidueIdentifier(residueName);
        }
        protected virtual bool CheckAtomName(string atName)
        {
            if (atName.StartsWith("H"))
                return false;

            return true;
        }
        public  string ParseAtomLine(Molecule molecule, string pdbLine, PDBMODE flag)
        {
            
            try
            {
                string atomName = pdbLine.Substring(12, 4).Trim();

                if (!CheckAtomName(atomName))
                    return "Wrong Atom name: " + atomName+" atom will be removed";

                string residueName = pdbLine.Substring(17, 3).Trim();
                if (!CheckResidue(residueName))
                    return "Incorrect residue name: "+residueName;

                this.AtomName = atomName;
                //ResidueName = ResidueIdentifier(residueName);
                tabParam[0] = (short)ResidueIdentifier(residueName);

                //ResidueSequenceNumber = Convert.ToInt16(pdbLine.Substring(22, 4));
                tabParam[2] = Convert.ToInt16(pdbLine.Substring(22, 4));

                //ChainIdentifier = (pdbLine.Substring(21, 1))[0];
                tabParam[1] = (short)(pdbLine.Substring(21, 1))[0];
                //if (ResidueName == 'O') ChainIdentifier = ' ';
                //if (tabParam[0] == 'O') ChainIdentifier = ' ';
                if (tabParam[0] == 'O') tabParam[1] =(short) ' ';
                else //if (ChainIdentifier == ' ') ChainIdentifier = '1';
                    if (tabParam[1] == ' ') tabParam[1] =(short) '1';
              
                    float x = float.Parse(pdbLine.Substring(30, 7), CultureInfo.InvariantCulture);
                    float y = float.Parse(pdbLine.Substring(38, 7), CultureInfo.InvariantCulture);
                    float z = float.Parse(pdbLine.Substring(46, 7), CultureInfo.InvariantCulture);


                    if (flag != PDBMODE.ONLY_SEQ)                    
                        Position = new Point3D(x, y, z);
                    
            }
            catch(Exception ex)
            {
                return "Error in reading ATOM line: " + ex.Message;
            }
			
            return null;
        }

    }
}
