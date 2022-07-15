using System;
using System.Collections.Generic;
using System.Text;
using uQlustCore;

namespace uQlustCore.PDB
{
    public class Residue
    {
       // static Dictionary<string, int> A,R,N,D,C,E,G,Q,H,I,L,K,M,F,P,S,T,W,Y,V,X;
       
        //static readonly Dictionary<char, Dictionary<string,int>> allAtoms;


        private char residueName;
        //private char chainIdentifier;
        //private short residueSequenceNumber;
        public short[] tabParam = new short[2];
        //private string residueIdentifier;
        private List<Atom> atoms;
       
        public Residue()
        {

        }
        static Residue()
        {

         /*   A = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5}};
            R = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CD",6},{"CG",7},{"CZ",8},{"NE",9},{"NH1",10},{"NH2",11}};
            N = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CG",6},{"OD1",7},{"ND2",8}};
            D = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CG",6},{"OD1",7},{"OD2",8}};
            C = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"SG",6}};
            E = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CD",6},{"CG",7},{"OE1",8},{"OE2",9}};
            G = new Dictionary<string, int> { { "CA", 1 }, { "N", 2 }, { "C", 3 }, { "O", 4 } };
            Q = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CD",6},{"CG",7},{"OE1",8},{"NE2",9}};
            H = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CG",6},{"ND1",7},{"CD2",8},{"CE1",9},{"NE2",10}};
            I = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CG1",6},{"CG2",7},{"CD1",8}};
            L = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CD1",6},{"CD2",7},{"CG",8}};
            K = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CD",6},{"CG",7},{"CE",8}};
            M = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CG",6},{"CE",7},{"SD",8}};
            F = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CG",6},{"CD1",7},{"CD2",8},{"CE1",9},{"CE2",10},{"CZ",11}};
            P = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CG",6},{"CD",7}};
            S = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"OG",6}};
            T = new Dictionary<string, int> { { "CB", 1 }, { "CA", 2 }, { "N", 3 }, { "C", 4 }, { "O", 5 }, { "OG1", 6 }, { "CG2", 7 } };
            W = new Dictionary<string, int> { { "CG", 1 }, { "CA", 2 }, { "N", 3 }, { "C", 4 }, { "O", 5 }, { "CB", 6 }, { "CD1", 7 }, { "CD2", 8 }, { "NE1", 9 }, { "CE2", 10 }, { "CE3", 11 }, { "CZ2", 12 }, { "CZ3", 13 }, { "CH2", 14 } };
            Y = new Dictionary<string, int> {{"CB",1},{"CA",2},{"N",3},{"C",4},{"O",5},{"CG",6},{"CD1",7},{"CD2",8},{"CE1",9},{"CE2",10},{"CZ",11},{"OH",12}};
            V = new Dictionary<string, int> { { "CB", 1 }, { "CA", 2 }, { "N", 3 }, { "C", 4 }, { "O", 5 }, { "CG1", 6 }, { "CG2", 7 } };
            X = new Dictionary<string, int> { { "CG", 1 }, { "CA", 2 }, { "N", 3 }, { "C", 4 }, { "O", 5 }, { "CB", 6 }, { "CD1", 7 }, { "CD2", 8 }, { "NE1", 9 }, { "CE2", 10 }, 
                                              { "CE3", 11 }, { "CZ2", 12 }, { "CZ3", 13 }, { "CH2", 14 },{"CE",15},{"ND1",16},{"OD1",17},{"OD2",18},
                                              {"NH1",19},{"NH2",20},{"SG",21},{"OG",22},{"SD",23}};
            allAtoms = new Dictionary<char, Dictionary<string,int>>
            {   
                {'A',A},{'R',R},{'N',N},{'D',D},{'C',C},{'E',E},{'G',G},{'Q',Q},{'H',H},{'I',I},{'L',L},{'K',K},{'M',M},
                {'F',F},{'P',P},{'S',S},{'T',T},{'W',W},{'Y',Y},{'V',V},{'X',X}         
            };            */
        }
        internal Residue(Molecule molecule, Atom atom,PDBMODE flag)
        {
          //  this.molecule = molecule;
            //this.residueName = atom.ResidueName;
            this.residueName = (char)atom.tabParam[0];
            //this.chainIdentifier = atom.ChainIdentifier;
            this.tabParam[0] = atom.tabParam[1];
            //this.residueSequenceNumber = atom.ResidueSequenceNumber;         
            //this.residueSequenceNumber = atom.tabParam[2];         
            this.tabParam[1] = atom.tabParam[2];         

            atoms = new List<Atom>();
            this.atoms.Add(atom);
            //AddAtom(atom);
            //this.residueIdentifier = Residue.GetResidueIdentifier(this.residueName);
        }
        /// <summary>
        /// Label used for atom tooltips.
        /// </summary>
        internal char ResidueName { get { return this.residueName; } }
        internal char ChainIdentifier { get { return (char)this.tabParam[0]; } }
        internal int ResidueSequenceNumber { get { return this.tabParam[1]; } }
        internal List<Atom> Atoms { get { return this.atoms; } }
       
        internal static bool IsAminoName(string residueName)
        {
            if (residueName == "ALA") return true;
            else if (residueName == "ARG") return true;
            else if (residueName == "ASP") return true;
            else if (residueName == "CYS") return true;
            else if (residueName == "GLN") return true;
            else if (residueName == "GLU") return true;
            else if (residueName == "GLY") return true;
            else if (residueName == "HIS") return true;
            else if (residueName == "ILE") return true;
            else if (residueName == "LEU") return true;
            else if (residueName == "LYS") return true;
            else if (residueName == "MET") return true;
            else if (residueName == "MSE") return true;
            else if (residueName == "PHE") return true;
            else if (residueName == "PRO") return true;
            else if (residueName == "SER") return true;
            else if (residueName == "THR") return true;
            else if (residueName == "TRP") return true;
            else if (residueName == "TYR") return true;
            else if (residueName == "VAL") return true;
            else if (residueName == "ASN") return true;
            else return false;
        }

        public static char GetResidueIdentifier(string residueName)
        {
            if (residueName == "HOH") return 'O';
            else if (residueName == "ALA") return 'A';
            else if (residueName == "ARG") return 'R';
            else if (residueName == "ASP") return 'D';
            else if (residueName == "CYS") return 'C';
            else if (residueName == "GLN") return 'Q';
            else if (residueName == "GLU") return 'E';
            else if (residueName == "GLY") return 'G';
            else if (residueName == "HIS") return 'H';
            else if (residueName == "ILE") return 'I';
            else if (residueName == "LEU") return 'L';
            else if (residueName == "LYS") return 'K';
            else if (residueName == "MET") return 'M';
            else if (residueName == "MSE") return 'M';
            else if (residueName == "PHE") return 'F';
            else if (residueName == "PRO") return 'P';
            else if (residueName == "SER") return 'S';
            else if (residueName == "THR") return 'T';
            else if (residueName == "TRP") return 'W';
            else if (residueName == "TYR") return 'Y';
            else if (residueName == "VAL") return 'V';
            else if (residueName == "ASN") return 'N';
            else return 'X';
        }
/*        public bool CheckSideChain()
        {
            int counter = 0;
            foreach(var item in atoms)
            {
                if (item.AtomName == "C" || item.AtomName== "CA" || item.AtomName == "N" || item.AtomName == "O")
                    continue;
                if (allAtoms.ContainsKey(this.residueName) && allAtoms[residueName].ContainsKey(item.AtomName))
                    counter++;    
            }
            if (!allAtoms.ContainsKey(this.residueName) && counter == 0)
                return true;

            if (counter == allAtoms[residueName].Keys.Count)
                return true;

            return false;
        }*/
       
    }
}
