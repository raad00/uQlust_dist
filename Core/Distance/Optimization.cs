using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using uQlustCore.PDB;

namespace uQlustCore.Distance
{
    public struct posMOL
    {
        public float[,] posmol1;
        public float[,] posmol2;
        public List<string>[] atoms;
    }
    public struct MultRes
    {
        public float[,] mult;
        public double val;
    }
    public class Optimization
    {

		public Optimization()
		{
		}
        public static float[,] GetFullStructure(MolData mol)
        {
            float[,] structure = new float [mol.indexMol.Length,3];
            for (int i = 0; i < mol.indexMol.Length; i++)
                if (mol.indexMol[i] != -1 && mol.mol.Residues[mol.indexMol[i]].ResidueName != 'X')
                {
                    //Take only first Atom (should be CA)
                    Atom atom=mol.mol.Residues[mol.indexMol[i]].Atoms[0];
                    structure[i, 0] = atom.Position.X;
                    structure[i, 1] = atom.Position.Y;
                    structure[i, 2] = atom.Position.Z;

                }
                else
                {
                    structure[i, 0] = float.MaxValue;
                    structure[i, 1] = float.MaxValue;
                    structure[i, 2] = float.MaxValue;
                }

            return structure;
        }
     /*   public static posMOL PrepareData(MolData mol1, MolData mol2, bool atomList = false)
        {
            posMOL globPosMol = new posMOL();

            globPosMol.posmol1 = new float[mol1.mol.Residues.Count, 3];
            globPosMol.posmol2 = new float[mol2.mol.Residues.Count, 3];

            for (int i = 0; i < mol1.mol.Residues.Count; i++)
            {
                globPosMol.posmol1[i, 0] = mol1.mol.Residues[i].Atoms[0].Position.X;
                globPosMol.posmol1[i, 1] = mol1.mol.Residues[i].Atoms[0].Position.Y;
                globPosMol.posmol1[i, 2] = mol1.mol.Residues[i].Atoms[0].Position.Z;
            }
            for (int i = 0; i < mol2.mol.Residues.Count; i++)
            {
                globPosMol.posmol2[i, 0] = mol2.mol.Residues[i].Atoms[0].Position.X;
                globPosMol.posmol2[i, 1] = mol2.mol.Residues[i].Atoms[0].Position.Y;
                globPosMol.posmol2[i, 2] = mol2.mol.Residues[i].Atoms[0].Position.Z;
            }

            float[] center = new float[3];

            CenterMol(globPosMol.posmol1, center);
            CenterMol(globPosMol.posmol2, center);

            return globPosMol;

        }*/
        public static posMOL PrepareData(MolData mol1,MolData mol2,bool atomList=false)
        {
            int count=0;
            posMOL globPosMol = new posMOL();
            Dictionary<string, int> atomNameMol1 = null;

            for (int i = 0; i < mol1.mol.Residues.Count; i++)
                if (mol1.mol.Residues[i].Atoms.Count > count)
                    count = mol1.mol.Residues[i].Atoms.Count;
           
           atomNameMol1 = new Dictionary<string, int>(count);

            count = 0;
            for (int i = 0; i < mol1.indexMol.Length; i++)
                if (mol1.indexMol[i] != -1 && mol2.indexMol[i] != -1)// && mol1.mol.Residues[mol1.indexMol[i]].ResidueName != 'X' && mol2.mol.Residues[mol2.indexMol[i]].ResidueName != 'X')
                {
                    atomNameMol1.Clear();
                    List<Atom> atoms = mol1.mol.Residues[mol1.indexMol[i]].Atoms;
                    for (int j = 0; j < atoms.Count;j++ )
                        atomNameMol1.Add(atoms[j].AtomName, j);

                    foreach (var item in mol2.mol.Residues[mol2.indexMol[i]].Atoms)
                        if (atomNameMol1.ContainsKey(item.AtomName))
                            count++;
                }
           // index1 = mol1.indexMol;
            //index2 = mol2.indexMol;

            if (globPosMol.posmol1==null || count != globPosMol.posmol1.GetLength(0))
            {
                globPosMol.posmol1 = new float[count, 3];
                globPosMol.posmol2 = new float[count, 3];
            }
            if (atomList)
            {
                globPosMol.atoms = new List<string>[mol1.indexMol.Length];
                for(int i=0;i<mol1.indexMol.Length;i++)
                {
                    globPosMol.atoms[i] = new List<string>();
                }
            }

            count = 0;            
            for (int i = 0; i < mol1.indexMol.Length; i++)
                if (mol1.indexMol[i] != -1 && mol2.indexMol[i] != -1)// && mol1.mol.Residues[mol1.indexMol[i]].ResidueName != 'X' && mol2.mol.Residues[mol2.indexMol[i]].ResidueName != 'X')
                {
                    atomNameMol1.Clear();
                    List<Atom> atoms = mol1.mol.Residues[mol1.indexMol[i]].Atoms;
                    for (int j = 0; j < atoms.Count; j++)
                        atomNameMol1.Add(atoms[j].AtomName, j);

                    foreach (var item in mol2.mol.Residues[mol2.indexMol[i]].Atoms)
                        if (atomNameMol1.ContainsKey(item.AtomName))
                        {
                            Atom aux = mol1.mol.Residues[mol1.indexMol[i]].Atoms[atomNameMol1[item.AtomName]];
                            globPosMol.posmol1[count, 0] = aux.Position.X;
                            globPosMol.posmol1[count, 1] = aux.Position.Y;
                            globPosMol.posmol1[count, 2] = aux.Position.Z;
                            
                            globPosMol.posmol2[count, 0] = item.Position.X;
                            globPosMol.posmol2[count, 1] = item.Position.Y;
                            globPosMol.posmol2[count, 2] = item.Position.Z;
                            if(atomList)
                                globPosMol.atoms[i].Add(item.AtomName);
                            count++;
                        }
                }

            float[] center = new float[3];
            
            CenterMol(globPosMol.posmol1, center);
            CenterMol(globPosMol.posmol2, center);

            return globPosMol;
        }
        
		static public void TransVec(float [,] mat,float [] center)
		{
            for (int i = 0; i < mat.GetLength(1); i++)
            {
                for (int j = 0; j < mat.GetLength(0); j++)
                    mat[j, i] -= center[i];
            }            
			
		}
		
        static private void FindCenter(float[,] mat,float [] center)
        {
            for (int i = 0; i < mat.GetLength(1); i++)
            {
                center[i] = 0;
                for (int j = 0; j < mat.GetLength(0); j++)
                    center[i] += mat[j, i];
                center[i] /= mat.GetLength(0);
            }
        }
        public static float [,] Transp(float [,] u)
        {
            float[,] newU = new float[u.GetLength(1), u.GetLength(0)];

            for (int i = 0; i < u.GetLength(0); i++)
                for (int j = 0; j < u.GetLength(1); j++)
                    newU[j, i] = u[i, j];

            return newU;

        }
        private int Det(float[,] u, float[,] v)
        {
            float[,] detM = MultMatrixTrans(u, v);
                        
            return Math.Sign(detM[0,0]*detM[1,1]*detM[2,2]+detM[0,1]*detM[1,2]*detM[2,0]+detM[0,2]*detM[1,0]*detM[2,1]-
                             detM[0,2]*detM[1,1]*detM[2,0]-detM[0,1]*detM[1,0]*detM[2,2]-detM[0,0]*detM[1,2]*detM[2,1]);
        }

        static public  float [] Cross(float []  x,float [] y) 
        {
            float [] z=new float[x.GetLength(0)];

            z[0] = x[1]*y[2]-x[2]*y[1];
            z[1] = -(x[0]*y[2]-x[2]*y[0]);
            z[2] = x[0]*y[1]-x[1]*y[0];

            return z;
        }

        public static float[,] MultMatrix(float[,] u, float[,] v)
        {
            float[,] mult = new float[u.GetLength(0), v.GetLength(1)];

            for (int i = 0; i < u.GetLength(0); i++)
            {                
                for (int j = 0; j < v.GetLength(1); j++)
                {
                    double sum = 0;
                    for (int k = 0; k < u.GetLength(1); k++)
                        sum += u[i, k]*v[k, j];

                    mult[i, j] =(float) sum;
                }
            }
            return mult;
        }
      
		public static float[,] MultMatrixTrans(float[,] u, float[,] v)
        {
            float[,] mult = new float[u.GetLength(0), v.GetLength(0)];

            for (int i = 0; i < u.GetLength(0); i++)
            {                
                for (int j = 0; j < v.GetLength(0); j++)
                {
                    double sum = 0;
                    for (int k = 0; k < u.GetLength(1); k++)
                        sum += ((double)u[i, k])*v[j, k];

                    mult[i, j] = (float)sum;
                }
            }
            return mult;
        }
        public static float[,] MultMatrixTrans(float[,] u, double[,] v)
        {
            float[,] mult = new float[u.GetLength(0), v.GetLength(0)];

            for (int i = 0; i < u.GetLength(0); i++)
            {
                for (int j = 0; j < v.GetLength(0); j++)
                {
                    double sum = 0;
                    for (int k = 0; k < u.GetLength(1); k++)
                        sum += ((double)u[i, k]) * v[j, k];

                    mult[i, j] = (float)sum;
                }
            }
            return mult;
        }
		 public static MultRes TransMultMatrix(float[,] u, float[,] v)
        {
            MultRes res;
            double val;
            float[,] mult = new float[u.GetLength(1), v.GetLength(1)];
             
            double G1=0,G2=0;

            for (int k = 0; k < u.GetLength(0); k++)
                for (int i = 0; i < u.GetLength(1); i++)
                {
                    G1 += u[k, i] * u[k, i];
                    G2 += v[k, i] * v[k, i];
                }

            val = (G1 + G2) / 2;
            for (int i = 0; i < u.GetLength(1); i++)
            {                
                for (int j = 0; j < v.GetLength(1); j++)
                {
                    double sum = 0;
                    for (int k = 0; k < u.GetLength(0); k++)
                        sum += ((double)u[k,i])*v[k,j];

                    mult[i, j] = (float)sum;
                }
            }
            res.mult = mult;
            res.val = val;
            return res;
        }
        static public void CenterMol(float[,] mol1, float []cent1)
        {
            Optimization.FindCenter(mol1, cent1);
            Optimization.TransVec(mol1, cent1);

        }
        public static float Rmsd(float[,] mol1, float[,] mol2, bool flag)
        {
            float rmsd;
           
            MultRes res = Optimization.TransMultMatrix(mol1, mol2);
            Qcp qcp = new Qcp(res.mult, res.val, mol1.GetLength(0));

            rmsd = (float)qcp.CalcRmsd();

            if (flag)
            {
                qcp.CalcRotMatrix();
                float[,] rot = new float[3, 3];
                /*for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        rot[i, j] = (float)qcp.rot[i, j];*/
                rot[0, 0] = (float)qcp.rot00; rot[0, 1] = (float)qcp.rot01; rot[0, 2] = (float)qcp.rot02;
                rot[1, 0] = (float)qcp.rot10; rot[1, 1] = (float)qcp.rot11; rot[1, 2] = (float)qcp.rot12;
                rot[2, 0] = (float)qcp.rot20; rot[2, 1] = (float)qcp.rot21; rot[2, 2] = (float)qcp.rot22;

               float [,]posMolRot = Optimization.MultMatrixTrans(mol2, rot);
               // posMolRot = TransMultMatrix(rot,mol2);

            }

            return rmsd;
            //return MultMatrixTrans(tmp, U);

        }
        public static float [,] TransMatrix(float [,] mol1, float [,] mol2)
		{
			MultRes res = Optimization.TransMultMatrix(mol1, mol2);           
			Qcp qcp = new Qcp(res.mult, res.val, mol1.GetLength(0));
            
            if (!qcp.CalcRotMatrix())
                return null;
            
            float[,] rot = new float[3, 3];
            rot[0, 0] = (float)qcp.rot00; rot[0, 1] = (float)qcp.rot01; rot[0, 2] = (float)qcp.rot02;
            rot[1, 0] = (float)qcp.rot10; rot[1, 1] = (float)qcp.rot11; rot[1, 2] = (float)qcp.rot12;
            rot[2, 0] = (float)qcp.rot20; rot[2, 1] = (float)qcp.rot21; rot[2, 2] = (float)qcp.rot22;
            /*for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    rot[i,j] = (float)qcp.rot[i,j];*/

            return rot;			
		}
   }

}
