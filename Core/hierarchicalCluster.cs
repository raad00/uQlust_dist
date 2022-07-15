using System;
using System.Collections;
using System.Collections.Specialized;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using uQlustCore;
using System.Data;
using uQlustCore.Interface;
using uQlustCore.Distance;


namespace uQlustCore
{
    [Serializable]
	public class HClusterNode
	{
		public List <string> setStruct=new List<string>();
        public List<List<byte>> setProfiles = new List<List<byte>>();
        public List <HClusterNode> joined=null;
        [NonSerialized]
        public HClusterNode parent = null;
		public string refStructure;
        public string dirName;
        public Dictionary<byte, int>[] stateFreq = null;
        public bool fNode=false;
        public bool flagSign = false;
        public int num;
		public int iterNum=0;
		public int levelNum=0;
		public int counter=0;
        public double consistency = 0;
        public Color color=Color.Green;
        public double levelDist;
        public double realDist;
        //[NonSerialized]
        //public GraphNode gNode;
        [NonSerialized]
        public bool mark=false;
        int kMin;


        public void ClearColors()
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            st.Push(this);

            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.joined != null)
                {
                    current.color = Color.Green;
                    foreach (var item in current.joined)
                        st.Push(item);
                }

            }
        }
        public void ColorNode(string structMark,int color)
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            st.Push(this);

            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.joined != null)
                {
                    if (current.setStruct.Contains(structMark))
                        current.color = Color.Red;
                    else
                        current.color = Color.Green;
                    foreach (var item in current.joined)
                        st.Push(item);
                }

            }

        }

        public Dictionary<HClusterNode,System.Drawing.Color> MarkNodes(List<string> toMark,System.Drawing.Color color)
        {
            Dictionary<HClusterNode, System.Drawing.Color> returnList = new Dictionary<HClusterNode, System.Drawing.Color>();
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            st.Push(this);

            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.joined == null || current.joined.Count == 0)
                {
                    foreach(var item in toMark)
                        if(current.setStruct.Contains(item))
                        {
                            returnList.Add(current,color);
                            break;
                        }
                }
                else
                    if (current.joined != null)
                        foreach (var item in current.joined)
                            st.Push(item);

            }
            return returnList;
        }
        public int SearchMaxDist()
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;
            int kMax;
            kMax = -(int)this.levelDist;
            kMin = (int)this.levelDist;

            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.levelDist > kMax)
                    kMax = (int)current.levelDist;

                if(current.levelDist<kMin)
                    kMin = (int)current.levelDist;

                if (current.joined != null)
                    foreach (var item in current.joined)
                        st.Push(item);
            }
            return kMax;
        }
        public Dictionary<double, List<HClusterNode>> GetClustersByLevels()
        {
            Dictionary<double, List<HClusterNode>> dic = new Dictionary<double, List<HClusterNode>>();

            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;
            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (!dic.ContainsKey(current.levelDist))
                    dic.Add(current.levelDist, new List<HClusterNode>());
                                    
                dic[current.levelDist].Add(current);

                if (current.joined != null)
                    foreach (var item in current.joined)
                        st.Push(item);
            }



            return dic;
        }
        public List<List<string>> GetClusters(int k)
        {
            List<List<string>> results = new List<List<string>>();
            Dictionary<HClusterNode,System.Drawing.Color> currentRes;
            Dictionary<HClusterNode, System.Drawing.Color> bestRes;
            Dictionary<HClusterNode, System.Drawing.Color> newRes;
            bool NotFinished = true;
            int min, max;
            max=SearchMaxDist();
            min = kMin;
            //max = kMax;

            bool direction = true;

            if (levelDist < joined[0].levelDist)
                direction = false;

            int threshold = min+(max - min) / 2;
            currentRes = CutDendrog(threshold);
            bestRes=currentRes;            
            if(currentRes.Count!=k)
                while (NotFinished)
                {
                    if (direction)
                    {
                        if (currentRes.Count < k)
                            max = threshold;
                        else
                            min = threshold;
                    }
                    else
                        if (currentRes.Count > k)
                            max = threshold;
                        else
                            min = threshold;

                    threshold = min+(max - min) / 2;
                    newRes = CutDendrog(threshold);
                    if (max-min>1 && newRes.Count!=k)             
                        NotFinished = true;
                    else
                        NotFinished = false;

                    if (Math.Abs(newRes.Count - k) < Math.Abs(bestRes.Count-k))                                                    
                            bestRes=newRes;
                        
                    currentRes = newRes;

                }

            foreach (var item in bestRes)            
                results.Add(item.Key.setStruct);

            return results;
        }
        public Dictionary<HClusterNode,System.Drawing.Color> CutDendrog(int distThreshold)
        {
            Dictionary<HClusterNode, System.Drawing.Color> returnList = new Dictionary<HClusterNode,System.Drawing.Color>();
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.levelDist <= distThreshold)
                {
                    returnList.Add(current, System.Drawing.Color.Red);
                    if (current.joined != null)
                        foreach (var item in current.joined)
                            if(item.levelDist==current.levelDist)
                                returnList.Add(item, System.Drawing.Color.Red);

                }
                else
                    if (current.joined != null)
                        foreach (var item in current.joined)
                            st.Push(item);

            }
            return returnList;
        }

        public void SaveNode(StreamWriter file,int num)
        {
            if (file != null)
            {
                file.WriteLine("================Cluster Num=" + num + "=======================");
                if (refStructure != null && refStructure.Length != 0)
                    file.WriteLine(refStructure);
                foreach (var item in setStruct)
                    if(item!=refStructure)
                        file.WriteLine(item);
            }
        }
        private void MakeStructures()
        {
            
            int counter=0;
            if (joined != null)
            {
                setStruct.Clear();
                foreach (var item in joined)
                    counter += item.setStruct.Count;
                setStruct = new List<string>(counter);

                foreach (var item in joined)
                    foreach (var str in item.setStruct)
                        setStruct.Add(str);
            }
        }
        public void RedoSetStructures()
        {
            if(joined==null)
                return;
            foreach (var item in joined)            
                   item.RedoSetStructures();                
                            
            this.MakeStructures();

        }

        public List<HClusterNode> GetLeaves()
        {
            List<HClusterNode> listH = new List<HClusterNode>();

            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;
            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();

                if (current.joined != null)
                    foreach (var item in current.joined)
                        st.Push(item);
                else
                    listH.Add(current);
            }

            return listH;

        }
        public bool IsVisible(HClusterNode node)
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            if (this == node)
                return true;


            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.joined != null)
                    foreach (var item in current.joined)
                    {
                        if (item == node)
                            return true;
                        st.Push(item);
                    }
            }
            return false;

        }
	}

	public class hierarchicalCluster:IProgressBar
	{
		DistanceMeasure dMeasure;
        AglomerativeType linkageType;
        public string mustRefStructure = null;
        string dirName;
        int min;
        int currentV = 0;
        int maxV = 1;
        int progressRead = 0;
        HierarchicalCInput hierOpt;

		public hierarchicalCluster (DistanceMeasure dMeasure,HierarchicalCInput hier,string dirName)
		{
			this.dMeasure=dMeasure;
            this.linkageType = hier.linkageType;
            this.dirName = dirName;
            progressRead = 0;
            hierOpt = hier;

		}
        public double ProgressUpdate()
        {
            double sumProgress = 0;
            double progress = dMeasure.ProgressUpdate();

            if (progressRead == 1)
                sumProgress = 0.05 + progress * 0.7;
            else
                sumProgress = 0.05 * progress;

            return sumProgress + 0.25 * ((double)currentV / maxV);

        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }


        public override string ToString()
        {
            return "Agglomerative "+linkageType.ToString();
        }
        private int LMinimalDist(List<HClusterNode> levelNodes)
        {
            
            int min = Int32.MaxValue;
            for (int i = 0; i < levelNodes.Count; i++)
            {
                HClusterNode refStruct = levelNodes[i];
                for (int j = i + 1; j < levelNodes.Count; j++)
                {
                    int dist = dMeasure.FindMinimalDistance(refStruct, levelNodes[j], linkageType).Key;
                    if (min > dist)
                    {
                        min = dist;
                    }
                }
            }
            return min;
        }
        /*private List<List<HClusterNode>> LevelMinimalDist(List<HClusterNode> levelNodes)
        {
            List<List<HClusterNode>> res = new List<List<HClusterNode>>();
            Dictionary<HClusterNode, int> test = new Dictionary<HClusterNode, int>();

            foreach (var item in levelNodes)
                if (!test.ContainsKey(item))
                    test.Add(item, 1);

            min = LMinimalDist(levelNodes);
            for (int i = 0; i < levelNodes.Count; i++)
            {                
                List<HClusterNode> lista = new List<HClusterNode>();
                for (int j = i + 1; j < levelNodes.Count; j++)
                {
                    int dist = dMeasure.FindMinimalDistance(levelNodes[i], levelNodes[j], linkageType).Key;
                    if (min == dist)
                    {
                        if (!levelNodes[j].fNode && !levelNodes[i].fNode)
                        {
                            lista.Add(levelNodes[i]);
                            lista.Add(levelNodes[j]);
                            levelNodes[i].fNode = true;
                            levelNodes[j].fNode = true;
                        }
                    }
                }
                if (lista.Count >= 2)
                    res.Add(lista);
            }

            foreach (var item in res)
                foreach (var it in item)
                    it.fNode = true;

            return res;
        }*/

       private List<List<HClusterNode>> LevelMinimalDist(List<HClusterNode> levelNodes)
        {
            List<List<HClusterNode>> res = new List<List<HClusterNode>>();
            byte[] e = new byte[levelNodes.Count];
            Dictionary<HClusterNode, int> test = new Dictionary<HClusterNode, int>();

            foreach (var item in levelNodes)
                if (!test.ContainsKey(item))
                    test.Add(item, 1);

            min = Int32.MaxValue;
           
            for (int i = 0; i < levelNodes.Count; i++)
            {
                List<HClusterNode> lista = new List<HClusterNode>();
              
                    for (int j = i + 1; j < levelNodes.Count; j++)
                    {
                        //if (levelNodes[j].fNode)
                        //    continue;
                        int dist = dMeasure.FindMinimalDistance(levelNodes[i], levelNodes[j], linkageType).Key;
                        if (min > dist)
                        {
                            min = dist;
                           foreach (var item in lista)
                                item.fNode = false;
                            lista.Clear();
                            foreach (var item in res)
                                foreach (var it in item)
                                    it.fNode = false;
                            res.Clear();

                            lista.Add(levelNodes[i]);
                            lista.Add(levelNodes[j]);
                            levelNodes[i].fNode = true;
                            levelNodes[j].fNode = true;
                          

                        }
                        else
                            if (min == dist)
                            {
                               if (!levelNodes[j].fNode && !levelNodes[i].fNode)
                                {
                                    lista.Add(levelNodes[i]);
                                    lista.Add(levelNodes[j]);
                                    levelNodes[i].fNode = true;
                                    levelNodes[j].fNode = true;
                                }


                            }
                    }
                if (lista.Count >= 2)
                    res.Add(lista);
            }

            foreach (var item in res)
                foreach (var it in item)
                    it.fNode = true;

            return res;
        }
       
        public ClusterOutput HierarchicalClustering(List <string> structures)
		{
			List<List<HClusterNode>> level =new List<List<HClusterNode>>();
			List<HClusterNode> levelNodes=new List<HClusterNode>();		
			List<HClusterNode> rowNodes=new List<HClusterNode>();
            ClusterOutput outCl = new ClusterOutput();
            int levelCount = 0;
			bool end=false;
			HClusterNode node;

            if (structures.Count <= 1)
            {
                outCl.hNode = new HClusterNode();
                outCl.hNode.setStruct = structures;
                outCl.hNode.refStructure = structures[0];
                outCl.hNode.levelDist = 0;
                outCl.hNode.joined = null;
                return outCl;
            }

            progressRead = 1;
            dMeasure.CalcDistMatrix(structures);

            
			for(int i=0;i<structures.Count;i++)
			{
				node=new HClusterNode();
				node.refStructure=structures[i];
				node.joined=null;
				node.setStruct.Add(structures[i]);
                node.levelNum = levelCount;
                
                node.levelDist = dMeasure.maxSimilarity;
                node.realDist = dMeasure.GetRealValue(node.levelDist);
				levelNodes.Add(node);
			}
            maxV = levelNodes.Count+1;
			level.Add(levelNodes);
						
			while(!end)
			{
                
				levelNodes=new List<HClusterNode>();
                List<List<HClusterNode>> rowList = LevelMinimalDist(level[level.Count - 1]);
                if (rowList.Count > 0)
                {
                    foreach (var item in rowList)
                    {
						node=new HClusterNode();
						node.joined=item;
                        node.levelDist = min;
                        node.realDist = dMeasure.GetRealValue(min);
                        node.levelNum = level.Count;
						for(int m=0;m<item.Count;m++)
						{
							node.setStruct.AddRange(item[m].setStruct);
                            item[m].fNode = true;
							
						}																					
                        node.refStructure = dMeasure.GetReferenceStructure(node.setStruct);
                        
                        List<string> refList = new List<string>();
                        foreach (var itemJoined in node.joined)
                            refList.Add(itemJoined.refStructure);

                        node.refStructure = null;
                        if (mustRefStructure != null)
                            foreach (var itemRef in refList)
                                if (itemRef == mustRefStructure)
                                    node.refStructure = mustRefStructure;

                        if (node.refStructure == null)
                        {
                            node.refStructure = dMeasure.GetReferenceStructure(node.setStruct, refList);
                        }


                        levelNodes.Add(node);
					}
					
					
				}
				if(levelNodes.Count>0)
				{
					level.Add(levelNodes);
					for(int i=0;i<level[level.Count-2].Count;i++)
					{
						if(!level[level.Count-2][i].fNode)
							level[level.Count-1].Add(level[level.Count-2][i]);
					}
                    currentV = maxV - levelNodes.Count;

				}				
				if(level[level.Count-1].Count==1)
					end=true;
				
			}
            

            outCl.hNode = level[level.Count - 1][0];
            outCl.hNode.levelNum = 0;

            


            //At the end level num must be set properly
            Queue<HClusterNode> qq = new Queue<HClusterNode>();
            HClusterNode h;
            for (int i = 0; i < level.Count; i++)
                for (int j = 0; j < level[i].Count; j++)
                    level[i][j].fNode = true;

            for (int i = 0; i < level.Count; i++)
                for (int j = 0; j < level[i].Count; j++)
                    if (level[i][j].fNode)
                    {
                        level[i][j].levelDist = Math.Abs(level[i][j].levelDist - dMeasure.maxSimilarity);
                        level[i][j].realDist = dMeasure.GetRealValue(level[i][j].levelDist);
                        level[i][j].fNode = false;
                    }



            
            qq.Enqueue(level[level.Count - 1][0]);
            while (qq.Count != 0)
            {
                h = qq.Dequeue();

                if (h.joined != null)
                    foreach (var item in h.joined)
                    {
                        item.levelNum = h.levelNum + 1;
                        qq.Enqueue(item);

                    }
            }

            outCl.hNode.dirName = dirName;
            outCl.clusters = null;
            outCl.juryLike = null;
            currentV = maxV;
            outCl.runParameters = hierOpt.GetVitalParameters();
            return outCl;
		}
		
	}

}

