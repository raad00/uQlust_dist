using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using uQlustCore;
using uQlustCore.PDB;


namespace uQlustCore.Distance
{
    public class MaxSub:Rmsd
    {

        public MaxSub(DCDFile dcd, string alignFile, bool flag, string refJuryProfile = null)
            : base(dcd, alignFile, flag, PDBMODE.ONLY_CA,refJuryProfile)
        {
            order = true;
            maxSimilarity = 0.0;
        }
        public MaxSub(string dirName, string alignFile, bool flag, string refJuryProfile = null)
            : base(dirName, alignFile, flag, PDBMODE.ONLY_CA,refJuryProfile)
        {
            this.dirName = dirName;
            this.alignFile = alignFile;
            this.flag = flag;
            this.refJuryProfile = refJuryProfile;

            order = true;
            maxSimilarity = 0.0;
        }
        public MaxSub(List<string> fileNames, string alignFile, bool flag, string refJuryProfile = null)
            : base(fileNames, alignFile, flag, PDBMODE.ONLY_CA, refJuryProfile)
        {
            
            order = true;
            maxSimilarity = 0.0;
        }
        public override string ToString()
        {
            return "MaxSub";
        }
        public override bool SimilarityThreshold(float threshold, float dist)
        {
            if (dist < threshold)
                return true;
            return false;

        }
        public override double GetRealValue(double v)
        {
            return (100 - v) / 100;
        }
        protected KeyValuePair<List<int>,float[,]> FindLongestSegment(int seedL,float threshold, float [,] mol1,float [,] mol2)
        {
            List<int> bestList = null;
            List<int> pList = new List<int>();
            List<int> cList = null;
            float[,] remTransMol = null;
            float[,] copyMol2 = null;
            copyMol2 = new float[mol2.GetLength(0), mol2.GetLength(1)];
            remTransMol = new float[mol2.GetLength(0), mol2.GetLength(1)];
            for (int i = 0; i < mol1.GetLength(0) - seedL; i++)
            {
                pList.Clear();
                for (int j = 0; j < seedL; j++)
                    pList.Add(i + j);

                KeyValuePair<List<int>, float[,]> extRes = Extend(pList, bestList, mol1, mol2, copyMol2,threshold);
                cList = extRes.Key;

                if (cList == null)
                    continue;

                if (bestList==null || cList.Count > bestList.Count)
                {
                    bestList = new List<int>(cList);
                    Buffer.BlockCopy(extRes.Value, 0, remTransMol, 0, copyMol2.Length * sizeof(float));
                }
                if (cList.Count == mol1.GetLength(0))
                    break;


            }
            KeyValuePair<List<int>, float[,]> res = new KeyValuePair<List<int>, float[,]>(bestList, remTransMol);

            return res;
        }

        public override int GetDistance(string refStructure, string modelStructure)
        {
            int seedL = 4;
            //float maxDist = 3.5f;
            float maxDist = 5.0f;

            if (!pdbs.molDic.ContainsKey(refStructure) || !pdbs.molDic.ContainsKey(modelStructure))
                return errorValue;

            
            posMOL locPosMol = Optimization.PrepareData(pdbs.molDic[refStructure], pdbs.molDic[modelStructure]);

            KeyValuePair<List<int>,float[,]> seg = FindLongestSegment(seedL, maxDist, locPosMol.posmol1, locPosMol.posmol2);
            if ( seg.Key == null)
                return 0;

            float sum = 0;
            float mDist = maxDist * maxDist;
            float rmsd = 0;
                foreach(var item in seg.Key)
                {
                    float dist = 0;
                    dist = calcDist(locPosMol.posmol1, seg.Value, item);
                    rmsd += dist;
                    dist *= dist;
                    sum += 100 / (1 + dist / mDist);
                }
            rmsd /= seg.Key.Count;
            return 100 - (int)((sum / pdbs.molDic[refStructure].molLength));//locPosMol.posmol1.GetLength(0)));
        }

        //Warning MaxSub has been reversed 0 is the best one, 1 is the worst
        /*public override int GetDistance(string refStructure, string modelStructure)
        {
            int seedL = 4;
            int maxS = 0;
            float maxDist = 3.5f;
            //float maxDist = 1.5f;
            List<int> pList = new List<int>();
            List<int> cList = null;
            float[,] remMol = null;
            float[,] remTransMol = null;
            float[,] copyMol2 = null;
            List<int> bestList = null;

         //   TestRmsd(refStructure, modelStructure);
            if (!pdbs.molDic.ContainsKey(refStructure) || !pdbs.molDic.ContainsKey(modelStructure))
                return int.MaxValue;

            posMOL locPosMol=Optimization.PrepareData(pdbs.molDic[refStructure], pdbs.molDic[modelStructure]);
            copyMol2 = new float[locPosMol.posmol2.GetLength(0), locPosMol.posmol2.GetLength(1)];

            remMol = new float[locPosMol.posmol1.GetLength(0), locPosMol.posmol1.GetLength(1)];
            remTransMol = new float[locPosMol.posmol1.GetLength(0), locPosMol.posmol1.GetLength(1)];

            //opt.CenterMols(opt.posMol1, opt.posMol2);
            if (bestList != null)
                bestList.Clear();
            for (int i = 0; i < locPosMol.posmol1.GetLength(0) - seedL; i++)
            {
                pList.Clear();
                for (int j = 0; j < seedL; j++)
                    pList.Add(i + j);

                KeyValuePair<List<int>, float[,]> extRes = Extend(pList, bestList,locPosMol.posmol1, locPosMol.posmol2, copyMol2, maxDist);
                cList = extRes.Key;

                if (cList == null)
                    continue;

                if (cList.Count > maxS)
                {
                    bestList = new List<int>(cList);
                    maxS = cList.Count;
             //       Buffer.BlockCopy(copyMol1, 0, remMol, 0, copyMol2.Length * sizeof(float));
                    Buffer.BlockCopy(extRes.Value, 0, remTransMol, 0, copyMol2.Length * sizeof(float));                    
                }
                if (cList.Count == locPosMol.posmol1.GetLength(0))
                    break;


            }
            if (remTransMol == null || bestList==null || bestList.Count==0)
                return int.MaxValue;
			
			
			
            float sum = 0;
            float mDist = maxDist * maxDist;
			float rmsd=0;
            foreach (var item in bestList)
            {
                float dist = 0;
                dist = calcDist(locPosMol.posmol1, remTransMol, item);
				rmsd+=dist;
                dist *= dist;
                sum += 100/(1 + dist/ mDist);
            }
			rmsd/=bestList.Count;
            return 100 - (int)((sum / locPosMol.posmol1.GetLength(0)));
        }*/
        private float [,] Rotate(List<int> indexList, float[,] pMol1, float[,] pMol2,float [,]copyMol)
        {
            float[,] indexMol1;
            float[,] indexMol2;
            float[,] transMatrix = null;
            float[,] rotMol = null;
			
			

			if (indexList.Count <= 2)
                return null;

			
            indexMol1 = new float[indexList.Count, 3];
            indexMol2 = new float[indexList.Count, 3];
            for (int j = 0; j < indexList.Count; j++)
            {
                for (int l = 0; l < 3; l++)
                {
                    indexMol1[j, l] = pMol1[indexList[j], l];
                    indexMol2[j, l] = pMol2[indexList[j], l];
                }
            }
            float[] center1 = new float[3];
            float[] center2 = new float[3];
            //opt.CenterMols(indexMol1, indexMol2);
            Optimization.CenterMol(indexMol1, center1);
            Optimization.TransVec(indexMol1, center1);
            Optimization.CenterMol(indexMol2, center2);
            Optimization.TransVec(indexMol2, center2);

            transMatrix = Optimization.TransMatrix(indexMol1,indexMol2);
            if (transMatrix == null)
                return null;

            Buffer.BlockCopy(pMol2, 0, copyMol, 0, pMol2.Length * sizeof(float));
            Optimization.TransVec(copyMol, center2);
            Optimization.TransVec(pMol1, center1);
           
			rotMol = Optimization.MultMatrixTrans(copyMol,transMatrix);          			

            return rotMol;
        }
        private float calcDist(float[,] d1, float[,] d2, int index)
        {
            double sum = 0;
            for (int p = 0; p < 3; p++)
            {
                double tmp = d1[index, p] - d2[index, p];
                sum += tmp * tmp;
            }
            sum = (float)Math.Sqrt(sum);

            return (float) sum;
        }
        KeyValuePair< List<int>,float[,]> Extend(List<int> cList,List<int> bestList, float[,] pMol1, float[,] pMol2,float [,]copyMol, float maxDist)
        {
            int loopK = 3;
            float[,] transMol=null;
            float[,] cMol1 = new float[pMol1.GetLength(0), pMol1.GetLength(1)];
            float[,] cMol2 = new float[pMol2.GetLength(0), pMol2.GetLength(1)];
            List<int> locList = new List<int>(cList);
            for (int i = 0; i < loopK; i++)
            {
                Buffer.BlockCopy(pMol2, 0, cMol2, 0, pMol2.Length * sizeof(float));
                Buffer.BlockCopy(pMol1, 0, cMol1, 0, pMol1.Length * sizeof(float));
                transMol=Rotate(locList, cMol1, cMol2,copyMol);

                if (transMol == null)
                    return new KeyValuePair<List<int>, float[,]>(null, null);

                locList.Clear();
                float val = (i + 1) * maxDist / loopK;
                for (int j = 0; j < cMol1.GetLength(0); j++)
                {
                    if (calcDist(transMol,cMol1,j) < val)
                        locList.Add(j);
                }
                if (locList.Count <= 1)
                    break;
            }
            if (bestList != null && bestList.Count > locList.Count || locList.Count < 2)
                return new KeyValuePair<List<int>, float[,]>(locList, transMol);

            float[,] indexMol1;
            float[,] indexMol2;
            indexMol1 = new float[locList.Count, 3];
            indexMol2 = new float[locList.Count, 3];
            for (int j = 0; j < locList.Count; j++)
            {
                for (int l = 0; l < 3; l++)
                {
                    indexMol1[j, l] = pMol1[locList[j], l];
                    indexMol2[j, l] = pMol2[locList[j], l];
                }
            }
            float[] center = new float[3];
      
            //opt.CenterMols(indexMol1, indexMol2);
            Optimization.CenterMol(indexMol1, center);
            Optimization.TransVec(pMol1, center);
            Optimization.CenterMol(indexMol2, center);
            Optimization.TransVec(pMol2, center);

            transMol=Rotate(locList, pMol1, pMol2,copyMol);
            if (transMol == null)
                return new KeyValuePair<List<int>, float[,]>(null, null);
			
            for (int j = locList.Count - 1; j >= 0; j--)
            {
				float ww=calcDist(transMol,pMol1,locList[j]);
                if (ww> maxDist)
                    locList.RemoveAt(j);
            }


            return new KeyValuePair<List<int>, float[,]>(locList, transMol);


        }

    }
}
