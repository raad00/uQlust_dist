using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using uQlustCore;
using System.Data;
using uQlustCore.Interface;
using uQlustCore.Distance;
namespace uQlustCore
{
    class FastDendrog:IProgressBar
    {
        HClusterNode root = new HClusterNode();
        public List<string> setStruct = new List<string>();
        List<HClusterNode> dendrogList = new List<HClusterNode>();
        List<HClusterNode> leaves = new List<HClusterNode>();
        List<HClusterNode> st = new List<HClusterNode>();
        const string clusterName = "Fast Dendrogram";
        int initNodesNum = 20;
        DistanceMeasure dMeasure;
        ClusterOutput final=new ClusterOutput();
        jury1D jury;
        AglomerativeType linkage;
        HammingConsensus consensus;
        bool useKMeans = false;
        string dirName;
        int kMeansIter;
        bool hConcensus;
        HierarchicalCInput input;
        double currentProgress = 0;
        int currentV = 0;
        int maxProgress = 1;
        int maxV = 1;
        IProgressBar progressObject=null;
        public FastDendrog(DistanceMeasure dMeasure,HierarchicalCInput input,string dirName)
        {
            this.dMeasure = dMeasure;
            this.initNodesNum = input.numInitNodes;
            this.linkage = input.linkageType;
            this.useKMeans = input.ItMeans;
            this.hConcensus = input.HConsensus;
            this.kMeansIter = input.ItNum;
            this.dirName = dirName;
            this.input = input;
          
        }
        public double ProgressUpdate()
        {

            if (progressObject != null)
                return currentProgress +1.0/maxProgress* progressObject.ProgressUpdate();

            return currentProgress;
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
            return clusterName;
        }
        ClusterOutput DivideSpaceKmeans(List<string> list)
        {
             kMeans km = new kMeans(dMeasure);
             return km.kMeansLevel(2,30, list);
        }
        ClusterOutput DivideSpaceHamming(List<string> list)
        {
            ClusterOutput output=new ClusterOutput();

            Dictionary<string, int> aux = new Dictionary<string, int>();
            ClusterOutput outC = jury.JuryOptWeights(list);

            if (outC == null)
                return null;

            List<string> clust1 = new List<string>();
            List<string> clust2 = new List<string>();
           

            //consensus.ToConsensusStates(list, jury.GetStructureStates(outC.juryLike[0].Key));
            consensus.ToConsensusStates(list,outC.juryLike[0].Key);
            foreach (var item in list)
                aux.Add(item, consensus.distanceOrdered[item]);

            if (useKMeans)
            {
                clust1.Add(outC.juryLike[0].Key);
                var item = aux.OrderByDescending(j => j.Value);
                int dist = item.First().Value;
                foreach (var ll in item)
                {
                    if (dist != ll.Value)
                        break;
                    else
                        clust2.Add(ll.Key);
                }
                Dictionary<string, double> dic = new Dictionary<string, double>();
                foreach (var it in outC.juryLike)
                    dic.Add(it.Key, it.Value);
                double min = Double.MaxValue;
                string rem = "";
                foreach (var it in clust2)
                {
                    if (dic[it] < min)
                    {
                        min = dic[it];
                        rem = it;
                    }
                }
                clust1.Add(rem);

                kMeans km = new kMeans(dMeasure);
                return km.kMeansRun(kMeansIter,list, clust1);

            }
            else
            {
               


                int i = 0;
                foreach (var item in aux.OrderBy(j => j.Value))
                {

                    if (i < list.Count() / 2)
                        clust1.Add(item.Key);
                    else
                        clust2.Add(item.Key);
                    i++;
                }


                output.clusters = new clusterRes();
                    

                output.clusters.list.Add(clust1);
                output.clusters.list.Add(clust2);
            }
            return output;
            

        }
   
        ClusterOutput DivideSpace1DJury(List<string> list)
        {
            
            ClusterOutput outC,finalOut=new ClusterOutput();

            outC=jury.JuryOptWeights(list);
            double dev = CalcStandDev(outC);
            List<string> clust1 = new List<string>();
            List<string> clust2 = new List<string>();

            dev += dev;

            if (useKMeans)
            {
                clust1.Add(outC.juryLike[0].Key);
                clust1.Add(outC.juryLike[outC.juryLike.Count - 1].Key);
                kMeans km = new kMeans(dMeasure);
                return km.kMeansRun(kMeansIter,list, clust1);

            }
            else
            {
                for (int i = 0; i < outC.juryLike.Count; i++)
                    if(i<outC.juryLike.Count / 2)
                        clust1.Add(outC.juryLike[i].Key);
                    else
                        clust2.Add(outC.juryLike[i].Key);

                finalOut.clusters = new clusterRes();                
                finalOut.clusters.list.Add(clust1);
                finalOut.clusters.list.Add(clust2);
            }
            return finalOut;
        }
        private double CalcStandDev(ClusterOutput outC)
        {
            if (outC.juryLike != null)
            {
                double sum = 0,avr,dev=0;
                for (int i = 0; i < outC.juryLike.Count; i++)
                {
                    sum += outC.juryLike[i].Value;
                }
                avr = sum / outC.juryLike.Count;
                for (int i = 0; i < outC.juryLike.Count; i++)
                {
                    double aux = avr - outC.juryLike[i].Value;
                    dev += aux * aux;
                }
                dev /= outC.juryLike.Count;
                return Math.Sqrt(dev);
            }
            return -1;
        }
        void FastCluster(HClusterNode parent)
        {
            HClusterNode c;
         //   if (parent.setStruct.Count > 2)
           // {               
                ClusterOutput outClust;

                //outClust = DivideSpaceKmeans(parent.setStruct);
                //outClust = DivideSpace1DJury(parent.setStruct);
                if (hConcensus)
                    outClust = DivideSpaceHamming(parent.setStruct);
                else
                    outClust = DivideSpace1DJury(parent.setStruct);

                if (outClust==null || outClust.clusters.list.Count <= 1)
                {
                    leaves.Add(parent);
                    return;
                }
                //dist = dMeasure.GetDistance(outClust.clusters[0][0], outClust.clusters[1][0]);                
//                if (!dMeasure.SimilarityThreshold(distThreshold,dist))
                    
                parent.joined = new List<HClusterNode>();                
                for (int i = 0; i < outClust.clusters.list.Count; i++)
                {
                  //      dist =(int)(dMeasure.GetDistance(outClust.clusters[0][0], root.setStruct[0]));
                        
                        c = new HClusterNode();
                        //    c.levelDist = dist;                        
                        c.setStruct = outClust.clusters.list[i];
                        if (c.setStruct.Count > 20)
                        {
                            parent.joined.Add(c);
                            st.Add(c);
                        }
                        else
                        {
                            leaves.Add(c);
                        }
                }
                
            //}
        }
        private void MakeDendrogs(AglomerativeType linkage)
        {
            ClusterOutput outCl;
            hierarchicalCluster dendrog = new hierarchicalCluster(dMeasure,input,dirName);
            currentV = 0;
            maxV = leaves.Count+1;
            double remProgress = currentProgress;
            for(int i=0;i<leaves.Count;i++)
            {
                HClusterNode c = leaves[i];
                dendrog.mustRefStructure = c.setStruct[0];
                outCl = dendrog.HierarchicalClustering(c.setStruct);
                dendrogList.Add(c);
                c.levelDist = outCl.hNode.levelDist;
                c.realDist = dMeasure.GetRealValue(c.levelDist);
                c.refStructure = outCl.hNode.refStructure;
                if(outCl.hNode.joined!=null)
                {
                    c.joined = new List<HClusterNode>();
                    foreach (var item in outCl.hNode.joined)
                        c.joined.Add(item);
                }
                currentV++;
                currentProgress=remProgress+1.0/maxProgress* (double)currentV / maxV;
            }
            maxV = currentV;
            currentProgress = remProgress;
        }
        /*private void CheckRefDistances()
        {
            HClusterNode current;
            hierarchicalCluster dendrog = new hierarchicalCluster(dMeasure);
            float dist;
            Stack <HClusterNode> localSt=new Stack<HClusterNode>();
            st.Clear();
            st.Push(root);
            current = root;
            while (st.Count != 0)
            {                
                current=st.Pop();
                st.Push(current.joined[0]);
                st.Push(current.joined[1]);

                localSt.Push(current);                
            }
            while(localSt.Count!=0)
            {
                    float dist2;

                    current=localSt.Pop();
                    dist = dMeasure.GetDistance(current.setStruct[0], current.joined[0].setStruct[0]);
                    dist2 = dMeasure.GetDistance(current.setStruct[0], current.joined[1].setStruct[0]);
                    //current.levelDist = (dist > dist2) ? dist : dist2;
                    //current.levelDist = (dist + current.joined[0].levelDist + dist2 + current.joined[1].levelDist) / 2;
                    current.levelDist = dMeasure.GetDistance(current.joined[0].setStruct[0], current.joined[1].setStruct[0]);
                
            }
        }*/
        private HClusterNode JoinNodes(List<HClusterNode> nodes)
        {
            HClusterNode node = new HClusterNode();
            node.joined = new List<HClusterNode>();
            node.setStruct = new List<string>();
            foreach (var item in nodes)
            {
                node.joined.Add(item);
                foreach (var itemN in item.setStruct)
                    node.setStruct.Add(itemN);
            }
            List<string> refList = null;
            if (node.joined != null)
            {
                refList = new List<string>();
                foreach (var item in node.joined)
                    refList.Add(item.refStructure);
            }
            string refStr = dMeasure.GetReferenceStructure(node.setStruct,refList);
            node.refStructure = refStr;
            for (int i = 0; i < node.setStruct.Count; i++)
                if (refStr == node.setStruct[i])
                {
                    refStr = node.setStruct[0];
                    node.setStruct[0] = node.setStruct[i];
                    node.setStruct[i] = refStr;
                    break;
                }

            return node;

        }
        private HClusterNode ConnectDendrogs(AglomerativeType linkage)
        {
            List <Dictionary<int,int>> sim = new List<Dictionary<int,int>>();
            HClusterNode rootNode;
            int maxV = 1000000000;
            int minV = maxV-1;

            while (minV != maxV && dendrogList.Count>2)
            {

                int[,] distanceM = new int[dendrogList.Count, dendrogList.Count];
                minV = maxV;
                for (int i = 0; i < dendrogList.Count; i++)
                {
                    for (int j = i + 1; j < dendrogList.Count; j++)
                    {
                         distanceM[i, j] = dMeasure.GetDistance(dendrogList[i].refStructure, dendrogList[j].refStructure);                        
                        //distanceM[i, j] = dMeasure.FindMinimalDistance(dendrogList[i], dendrogList[j],linkage);                        
                        if (distanceM[i, j] < minV)
                            minV = distanceM[i, j];
                    }
                }

                if (minV != maxV)
                {
                    sim.Clear();
                    for (int i = 0; i < dendrogList.Count; i++)
                    {
                        Dictionary<int, int> aux = new Dictionary<int, int>();
                        aux.Add(i, 0);
                        for (int j = i + 1; j < dendrogList.Count; j++)
                        {
                            if (distanceM[i, j] == minV)
                                aux.Add(j, 0);
                        }
                        if (aux.Keys.Count > 1)
                            sim.Add(aux);

                    }
                    for (int i = 0; i < sim.Count; i++)
                    {                        
                        for (int j = i + 1; j < sim.Count; j++)
                        {
                            foreach (var item in sim[j].Keys)
                                if (sim[i].ContainsKey(item))
                                {
                                    foreach (var itemCopy in sim[j].Keys)
                                        if (!sim[i].ContainsKey(itemCopy))
                                            sim[i].Add(itemCopy, 0);

                                    sim.RemoveAt(j);
                                    i = -1;
                                    j=sim.Count;
                                    break;
                                }
                        }
                    }
                    List<HClusterNode> lNodes = new List<HClusterNode>();
                    List<int> removeList = new List<int>();
                    for (int n = sim.Count-1; n >=0; n--)
                    {
                        HClusterNode node = new HClusterNode();
                        node.joined = new List<HClusterNode>();
                        node.setStruct = new List<string>();
                        lNodes.Clear();
                        foreach (var item in sim[n].Keys)
                            if(!lNodes.Contains(dendrogList[item]))
                                lNodes.Add(dendrogList[item]);

                        node = JoinNodes(lNodes);
                        node.levelDist = minV;
                        node.realDist = dMeasure.GetRealValue(minV);
                        List<int> keys = new List<int>(sim[n].Keys);
                        keys.Sort();
                        dendrogList[keys[0]] = node;

                        for (int i = keys.Count - 1; i >= 1;i-- )
                        {
                            if(!removeList.Contains(keys[i]))
                                removeList.Add(keys[i]);
                               // dendrogList.RemoveAt(keys[i]);
                        }


                    }
                    removeList.Sort();
                    for (int i = removeList.Count - 1; i >= 0; i--)
                        dendrogList.RemoveAt(removeList[i]);
                }
            }

            if (dendrogList.Count > 1)
                rootNode = JoinNodes(dendrogList);
            else
                rootNode = dendrogList[0];

            return rootNode;
        }
        public List<HClusterNode> RearangeDendrogram(HClusterNode rRoot,double dist)
        {
            List <HClusterNode> rList=new List<HClusterNode>();

            Queue<HClusterNode> lQueue = new Queue<HClusterNode>();
            lQueue.Enqueue(rRoot);
            while (lQueue.Count != 0)
            {
                HClusterNode h = lQueue.Dequeue();
                if(h.joined!=null)
                for (int i = 0; i < h.joined.Count; i++)
                {

                    if (h.joined[i].levelDist < dist )
                        rList.Add(h.joined[i]);
                    else
                        if (h.joined.Count > 0)
                            lQueue.Enqueue(h.joined[i]);
                }
            }

            return rList;
        }
        public void PrepareList()
        {
            double min = 100000;
            double max = 0;
            List<HClusterNode> newList = new List<HClusterNode>();
            List<HClusterNode> finalList=new List<HClusterNode>();
            List<HClusterNode> auxList;
            for (int i = 0; i < dendrogList.Count; i++)
            {
                if (dendrogList[i].levelDist < min)
                    min = dendrogList[i].levelDist;
                if (dendrogList[i].levelDist >max)
                    max = dendrogList[i].levelDist;

            }
            min = min+(max - min) / 2;
            for (int i = 0; i < dendrogList.Count; i++)
            {
                if (dendrogList[i].levelDist > min)
                    newList.Add(dendrogList[i]);
                else
                    finalList.Add(dendrogList[i]);
            }
            
            for (int i = 0; i < newList.Count; i++)
                {
                auxList = RearangeDendrogram(newList[i],min);
                if (auxList.Count > 0)
                    foreach (var item in auxList)
                        finalList.Add(item);
            }
            dendrogList = new List<HClusterNode>(finalList);
        }
        public ClusterOutput Run(List<string> structs)
        {
            maxProgress=5;
            currentProgress=0;
            if (hConcensus)
            {
                maxProgress++;
                consensus = new HammingConsensus(dMeasure.dirName, null, false, input.consensusProfile);
                progressObject = consensus;
                consensus.InitMeasure();
                currentProgress += 1.0 / maxProgress;
            }
            jury = new jury1D();
            progressObject = jury;
            currentProgress += 1.0 / maxProgress;
            progressObject = null;
            jury.PrepareJury(dMeasure.dirName, dMeasure.alignFile, input.jury1DProfileFast);
            currentProgress += 1.0 / maxProgress;
            ClusterOutput clOut = new ClusterOutput();



            root.setStruct = structs;
//            if(hConcensus)
//                consensus.ToConsensusStates(structs);
                     
            FastCluster(root);
            maxV = initNodesNum;
            while (st.Count>0 && (leaves.Count+st.Count)<initNodesNum)
            {
                st.Sort(
                        delegate(HClusterNode p1, HClusterNode p2)
                        {                                        
                            return p2.setStruct.Count.CompareTo(p1.setStruct.Count);
                        }             
                );

                HClusterNode node = st[0];
                st.RemoveAt(0);
                FastCluster(node);
                currentV += leaves.Count + st.Count;
            }
            currentV = maxV;
            currentProgress+=1.0/maxProgress;
            while (st.Count > 0)
            {
                HClusterNode node = st[0];
                st.RemoveAt(0);
                leaves.Add(node);
            }
            MakeDendrogs(linkage);
            currentProgress+=1.0/maxProgress;
            PrepareList();
            root = ConnectDendrogs(linkage);
            root.levelDist = root.SearchMaxDist();
            root.realDist = dMeasure.GetRealValue(root.levelDist);
            //CheckRefDistances();
            //dendrogList = RearangeDendrogram(root);
            //root = ConnectDendrogs();
            clOut.hNode = root;
            currentProgress+=1.0/maxProgress;
            return clOut;
        }
    }
}
