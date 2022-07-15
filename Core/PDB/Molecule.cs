using System;
using System.ComponentModel;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using uQlustCore;


namespace uQlustCore.PDB
{
    public class Molecule
    {
        protected List<Atom> atoms;
        protected List<Residue> residues;
        protected List<Chain> chains;
        protected PDBMODE flag;    
        CHAIN_MODE flagChain=CHAIN_MODE.SINGLE;
        internal Molecule(PDBMODE flag,CHAIN_MODE flagChain=CHAIN_MODE.SINGLE)
        {
            this.flag = flag;
            this.flagChain=flagChain;
        }
        public Molecule()
        {
            chains = new List<Chain>();
        }
        public virtual bool ReadMolecule(StreamReader pdbStream)
        {
            bool res;
            res = this.ReadAtoms(pdbStream);
            if (!res)
                return false;
            this.CreateResidues();
            this.CreateChains();

            atoms.Clear();
            if(chains.Count==0)
                return false;
            residues = chains[0].Residues;

            if (flag == PDBMODE.ONLY_SEQ)
            {
                atoms.Clear();
                residues.Clear();
                foreach (var item in chains)
                    item.ClearChain();
            }

            return true;
        }
        
        
        public List<Atom> Atoms { get { return this.atoms; } }

        public List<Residue> Residues { get { return this.residues; } }

        public  List<Chain> Chains { get { return this.chains; } }

        public void RemoveChain(char selChain)
        {
            for(int i=0;i<chains.Count;i++)
            {
                if(chains[i].ChainIdentifier==selChain)
                {
                    chains.RemoveAt(i);
                    return;
                }
            }
        }

        public Chain CuttMoleculeToSEQ(string seq,char selChain)
        {
            int start, end;
            foreach(var item in chains)
            {
                if(item.ChainIdentifier==selChain)
                {
                    StreamWriter er = new StreamWriter("query");
                    er.WriteLine(">test");
                    er.WriteLine(seq);
                    er.Close();
                    er = new StreamWriter("baza");
                    er.WriteLine(">baza");
                    er.WriteLine(item.chainSequence);
                    er.Close();
                    ProcessStartInfo startInfo = new ProcessStartInfo();
                    startInfo.FileName = "C:\\Blast\\bin\\formatdb.exe";
                    startInfo.Arguments = "-i baza -p T";
                    startInfo.ErrorDialog = false;
                    startInfo.RedirectStandardOutput = true;
                    startInfo.UseShellExecute = false;
                    startInfo.CreateNoWindow = true;
                    try
                    {
                        Process exeProcess = Process.Start(startInfo);

                        exeProcess.WaitForExit();
                        startInfo.FileName = "C:\\Blast\\bin\\blastpgp";
                        startInfo.Arguments = "-i query -j 1 -d baza -o align";
                        exeProcess = Process.Start(startInfo);
                        exeProcess.WaitForExit();
                    }
                    catch(Exception)
                    {                      
                        return null;
                    }
                    using (StreamReader rr = new StreamReader("align"))
                    {
                        string line = "";
                        while (line != null && !line.Contains("Sbjct:"))
                        {
                            line = rr.ReadLine();
                        }
                        if (line == null)
                            return null;
                        line = line.Replace("  ", " ");
                        string[] tmp = line.Split(' ');
                        start = Convert.ToInt32(tmp[1]);
                        end = Convert.ToInt32(tmp[tmp.Length-1]);
                        while (line != null)
                        {
                            line = rr.ReadLine();
                            if (line.Contains("Score") || line.Contains("Database"))
                                break;
                            if (line.Contains("Sbjct:"))
                            {
                                line = line.Replace("  ", " ");
                                tmp = line.Split(' ');
                                end = Convert.ToInt32(tmp[tmp.Length - 1]);
                            }
                        }

                    }
                    if (end - start > item.chainSequence.Length)
                        return null;

                    Chain c = new Chain(selChain);


                    c.residues=item.CuttChain(start - 1, end - 1);
                    c.ChainIdentifier = selChain;
                    c.CreateChainString();
                    //Console.WriteLine("seq="+seq + "\nend=" +item.chainSequence);
                    return c;
                   
                }
            }
            return null;
        }

        internal virtual bool ReadAtoms(StreamReader pdbStream)
        {
            List<Atom> auxList = new List<Atom>();
            //using (StreamReader pdbReader = new StreamReader(pdbStream))
            //{
                string pdbLine = pdbStream.ReadLine();
               
                while (pdbLine != null)
                {
                    if (pdbLine.StartsWith("ENDMDL") /*|| pdbLine.StartsWith("TER") */|| pdbLine.StartsWith("END")) break;
                    if(flagChain==CHAIN_MODE.SINGLE)
                        if(pdbLine.StartsWith("TER"))
                            break;
                    if (pdbLine.StartsWith("ATOM") || pdbLine.StartsWith("HETATM"))
                    {
                        if (pdbLine.Contains("\t"))
                        {
                            ErrorBase.AddErrors("Error in file: " + ((FileStream)pdbStream.BaseStream).Name + " " + "ATOM line containes tab what is not allowed");
                            return false;
                        }
                        Atom atom = new Atom();
                        string error=atom.ParseAtomLine(this, pdbLine,flag);
                        if (error==null)
                        {
                            //Check if the atom already exists!!
                            if (flag == PDBMODE.ONLY_CA || flag == PDBMODE.CA_CB || flag == PDBMODE.ONLY_CB || flag==PDBMODE.ONLY_SEQ)
                            {
                                if (flag == PDBMODE.ONLY_CA || flag == PDBMODE.CA_CB || flag == PDBMODE.ONLY_SEQ)
                                    if (atom.AtomName == "CA")                                                                            
                                        auxList.Add(atom);                                    


                                if (flag == PDBMODE.ONLY_CB || flag == PDBMODE.CA_CB )//|| (atom.ResidueName=='G' && flag==PDBMODE.ONLY_CA))
                                    if (atom.AtomName == "CB")
                                        auxList.Add(atom);
                            }
                            else
                                auxList.Add(atom);
                        }
                        else
                            ErrorBase.AddErrors("Error in file: " + ((FileStream)pdbStream.BaseStream).Name+" "+error);
      
                    }

                    pdbLine = pdbStream.ReadLine();
                }
                this.atoms = new List<Atom>(auxList);
          //  }
            return true;
        }
        
        private void CreateResidues()
        {
            List<Residue> auxRes = new List<Residue>();
            List<Atom> localAtoms = new List<Atom>();
            //this.residues = new List<Residue>();

            Residue residue = null;
           
            foreach (Atom atom in this.atoms)
            {
                //if (residue == null ||  atom.ResidueSequenceNumber != residue.ResidueSequenceNumber||
                if (residue == null || atom.tabParam[2] != residue.ResidueSequenceNumber ||
                    atom.tabParam[1] != residue.ChainIdentifier)
                    //atom.ChainIdentifier != residue.ChainIdentifier)
                {
                    residue = new Residue(this, atom,flag);                    
                    auxRes.Add(residue);
                }
                else
                {
                    
                    bool test = false;                   
                    foreach(var item in residue.Atoms )
                        if (atom.AtomName == item.AtomName)
                        {
                            ErrorBase.AddErrors("Residue " + residue.ResidueName + " has two the same atoms " + atom.AtomName);
                            test = true;
                        }
                    if (!test)
                        residue.Atoms.Add(atom);
                        //residue.AddAtom(atom);
                    
                }
                atom.tabParam = null;

            }
            this.residues = new List<Residue>(auxRes);
            auxRes = null;

        }
       

        private void CreateChains()
        {
            this.chains = new List<Chain>();

            Chain chain = null;

            foreach (Residue residue in this.residues)
            {
                    if (chain == null || residue.ChainIdentifier != chain.ChainIdentifier)
                    {
                        chain = new Chain(residue.ChainIdentifier);
                        this.chains.Add(chain);
                    }

                    chain.Residues.Add(residue);

                    residue.tabParam = null;
            }

            for(int i=0;i<chains.Count;i++)
              chains[i].CreateChainString();

        }
    }
}
