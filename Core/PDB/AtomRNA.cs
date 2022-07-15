using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace uQlustCore.PDB
{
    class AtomRNA:Atom
    {
        static Dictionary<string, int> allowedATomNames = new Dictionary<string, int>(){{"P",1},{"O1P",1},{"O2P",1},{"O5*",1},{"C5*",1},{"C4*",1},{"O4*",1},{"C3*",1},{"O3*",1},{"C2*",1},{"O2*",1},{"C1*",1},{"N1",1},
                                                                                    {"C2",1},{"N2",1},{"N3",1},{"C4",1},{"C5",1},{"C6",1},{"O6",1},{"N7",1},{"C8",1},{"N9",1},{"OP1",1},{"OP2",1},{"O5'",1},{"C5'",1},{"C4'",1},{"O4'",1},{"C3'",1},{"O3'",1},{"C2'",1},{"O2'",1},{"C1'",1},
                                                                                    };

        protected override bool CheckResidue(string residueName)
        {
            //if (Residue.IsAminoName(residueName))
            //    return true;

            return true;
        }
        protected override bool CheckAtomName(string atName)
        {
         //   if (atName.StartsWith("H"))
           //     return false;

            return true;
        }
        protected override char  ResidueIdentifier(string residueName)
        {
            return ResidueRNA.GetResidueIdentifier(residueName);
        }
    }
}
