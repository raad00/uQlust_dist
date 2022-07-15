using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace uQlustCore.PDB
{
    class ResidueRNA:Residue
    {
        static Dictionary<string, int> residueSize=new Dictionary<string,int>(){{"G",34},{"C",31},{"A",34},{"U",31}};//G,C,U,A;


        static ResidueRNA()
        {

        }
        public static new char GetResidueIdentifier(string residueName)
        {
            if (residueName.Contains("G") || residueName.Contains("g")) return 'G';
            else if (residueName.Contains("C") || residueName.Contains("c")) return 'C';
            else if (residueName.Contains("A") || residueName.Contains("a")) return 'A';
            else if (residueName.Contains("U") || residueName.Contains("u")) return 'U';
            else return 'X';
        }

    }
}
