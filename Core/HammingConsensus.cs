using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using uQlustCore.Distance;

namespace uQlustCore
{
    class HammingConsensus : HammingBase
    {
        List<byte> consensusStates = new List<byte>();
        Dictionary<string, string> consensus = new Dictionary<string, string>();
        public Dictionary<string,int> distanceOrdered = new Dictionary<string,int>();
       
        public HammingConsensus(string dirName, string alignFile, bool flag, string consensusProfile)
            :base(dirName, alignFile,flag, consensusProfile)
        {
        }
        public HammingConsensus(List<string> fileNames, string alignFile, bool flag, string consensusProfile)
            :base(fileNames, alignFile, flag,consensusProfile)
        {
        }

        public override void InitMeasure()
        {
            if (fileNames != null)
                base.InitMeasure(fileNames, alignFile, flag, refJuryProfile);
            else
                base.InitMeasure(dirName, alignFile, flag, refJuryProfile);
        }

        public override string ToString()
        {
            return "HAMMING";
        }
        private string TransformToConsensusStates(string structName)
        {
            byte locState;
            string consensusRepr = "";

            if (!stateAlign.ContainsKey(structName))
                return null;

            for (int i = 0; i < stateAlign[structName].Count; i++)
            {
                locState = stateAlign[structName][i];
                if (locState == 0)
                    continue;

                if (consensusStates[i] == locState)
                    consensusRepr += "0";
                else
                    consensusRepr += "1";

            }
            return consensusRepr;
        }
        private int DistanceToConsensus(string structName)
        {
            int dist = 0;
            for (int i = 0; i < consensus[structName].Length; i++)
                if (consensus[structName][i] != '0')
                    dist++;

            return dist;
        }
        private void CalcAllDistances(List<string> structNames)
        {
            consensus.Clear();
            distanceOrdered.Clear();
            foreach (var item in structNames)
            {
                consensus.Add(item, TransformToConsensusStates(item));
                int dist = DistanceToConsensus(item);
                distanceOrdered.Add(item, dist);
            }

        }
        public void ToConsensusStates(List<string> structNames, string newConsensusStates)
        {
            List<string> states = new List<string>();
            if (!stateAlign.ContainsKey(newConsensusStates))
            {
                states=null;
                return;
            }

            consensusStates=stateAlign[newConsensusStates];

            CalcAllDistances(structNames);
        }
        public void ToConsensusStates(List<string> structNames, List<byte> newConsensusStates)
        {
            consensusStates = newConsensusStates;
            CalcAllDistances(structNames);
        }
    }
}
