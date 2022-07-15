
using System.Collections.Generic;
using System.Text;

namespace uQlustCore.PDB
{
   
    public class Chain
    {
        public string chainSequence;
        private char chainIdentifier;
        public List<Residue> residues;

        internal Chain(char chainIdentifier)
        {
            this.chainIdentifier = chainIdentifier;
            this.residues = new List<Residue>();
        }
        public List<Residue> CuttChain(int start,int end)
        {
            List<Residue> newResidue = new List<Residue>();
          
            for (int i = start; i <= end; i++)
                newResidue.Add(residues[i]);

            return newResidue;
            

        }
        public void ClearChain()
        {
            residues.Clear();
        }
        internal char ChainIdentifier { get { return this.chainIdentifier; } set { this.chainIdentifier = value; } }

        internal List<Residue> Residues { get { return this.residues; } }
        public void CreateChainString()
        {
            StringBuilder st = new StringBuilder(residues.Count);           
            for (int i = 0; i < residues.Count; i++)
                    st.Append(residues[i].ResidueName);

            chainSequence = st.ToString();
        }
    }
}
