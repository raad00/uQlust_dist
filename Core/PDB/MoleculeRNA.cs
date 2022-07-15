using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using uQlustCore;

namespace uQlustCore.PDB
{
    class MoleculeRNA : Molecule
    {
        public MoleculeRNA(PDBMODE flag):base(flag)
        {

        }
      /*  public override bool ReadMolecule()
        {
            return this.ReadAtoms(pdbStream);
        }*/
        internal override bool ReadAtoms(StreamReader pdbStream)
        {
            this.atoms = new List<Atom>();

//            using (StreamReader pdbReader = new StreamReader(pdbStream))
  //          {
                string pdbLine = pdbStream.ReadLine();

                while (pdbLine != null)
                {
                    if (pdbLine.StartsWith("ENDMDL") || pdbLine.StartsWith("TER") || pdbLine.StartsWith("END")) break;

                    if (pdbLine.StartsWith("ATOM"))
                    {
                        if (pdbLine.Contains("\t"))
                        {
                            ErrorBase.AddErrors("ATOM line containes tab what is not allowed");
                            return false;
                        }
                        AtomRNA atom = new AtomRNA();
                        string error = atom.ParseAtomLine(this, pdbLine, flag);
                        if (error==null)
                        {
                            if(flag==PDBMODE.ONLY_P || flag==PDBMODE.ONLY_P_AND_C4)
                            {
                                if(atom.AtomName=="P")
                                    this.atoms.Add(atom);
                                if(atom.AtomName=="C4*" || atom.AtomName == "C4'")
                                    this.atoms.Add(atom);
                        }
                            else
                                this.atoms.Add(atom);
                        }
                        else
                        {
                            ErrorBase.AddErrors("Error in file: " + ((FileStream)pdbStream.BaseStream).Name + " " + error);
                        }
                                
                    }

                    pdbLine = pdbStream.ReadLine();
                }
    //        }
            return true;
        }

    }
}
