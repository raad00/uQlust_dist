using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using uQlustCore.PDB;

namespace uQlustCore.Profiles
{
    class GaussIntegralsRNA:GaussIntegrals
    {
        public GaussIntegralsRNA()
        {
            dirSettings.Load();
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.RNA);
            profileName = "GaussIntegralsRNA ";
            contactProfile = profileName;
            AddInternalProfiles();
            init_reference_gauss_integrals(null);
            maxV = 1;
        }
        public override void AddInternalProfiles()
        {
            profileNode node = new profileNode();

            node.profName = "GaussIntegralsRNA";
            node.internalName = "GaussIntegralsRNA";
            for (int i = 1; i <= 10; i++)
                node.AddStateItem(i.ToString(), i.ToString());

            InternalProfilesManager.AddNodeToList(node, typeof(GaussIntegralsRNA).FullName);

        }

        List<Atom> SelectAtoms(Chain chain,int countAtoms)
        {
            List<Atom> atoms = new List<Atom>();
            if (countAtoms > 2 * 876)
            {
                //Remove C4
                foreach (var res in chain.Residues)
                {
                    foreach (var atom in res.Atoms)
                        if (atom.AtomName == "P")
                            atoms.Add(atom);
                }
            }
            else
            {
                foreach (var res in chain.Residues)
                {
                    foreach (var atom in res.Atoms)
                        atoms.Add(atom);
                }
            }
            if (atoms.Count > 875)
            {
                double step = (atoms.Count -875.0)/atoms.Count;
                double aux = atoms.Count / step;
                double sum = 0;
                List<int> toRemove = new List<int>();
                for (int i = 0; i < atoms.Count; i++,sum+= step)
                {
                    if (sum >= 1)
                    {
                        toRemove.Add(i);
                        sum = sum-1;
                    }
                }
                for (int i = toRemove.Count - 1; i >= 0; i--)
                    atoms.RemoveAt(toRemove[i]);
            }
            if (atoms.Count > 876)
                Console.Write("UPS");

            return atoms;
        }
        protected override void MakeProfiles(string strName, MolData molDic, StreamWriter wr)
        {
            List<Atom> atoms = new List<Atom>();

            foreach (var chain in molDic.mol.Chains)
            {
                int atomsNum = 0;

                for (int i = 0; i < chain.Residues.Count; i++)
                    atomsNum += chain.Residues[i].Atoms.Count;

                if (atomsNum < 12)
                    continue;

                atoms = SelectAtoms(chain, atomsNum);

                double[] git_vector = generate_gauss_integrals(atoms);
                if (git_vector == null)
                    continue;
                if (molDic.mol.Chains.Count > 1)
                    wr.WriteLine(">" + strName + "|" + chain.ChainIdentifier);
                else
                    wr.WriteLine(">" + strName);
                string name = profileName.Remove(profileName.Length - 1, 1);
                wr.Write(name);
                for (int i = 0; i < git_vector.Length; i++)
                    wr.Write(" " + git_vector[i]);

                wr.WriteLine();
            }
            molDic.CleanMolData();
        }
        public override List<INPUTMODE> GetMode()
        {
            List<INPUTMODE> oList = new List<INPUTMODE>();
            oList.Add(INPUTMODE.RNA);

            return oList;
        }

    }
}
