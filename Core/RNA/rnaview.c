/* the input coordinate file either is PDF format or CIF */
//#include "stdafx.h"
#include <stdio.h>
//#include <tchar.h>
//#include <wtypes.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include "nrutil.h"
#include "rna.h"

#define SYSTEM_WINDOWS


#ifdef SYSTEM_WINDOWS
#include <wtypes.h>
#endif
#ifndef SYSTEM_WINDOWS
#define LPSTR char *
#endif

struct RNAMemory
{
	char *profile;
	char *sequence;
	char *secStruct;
};
//int errorNum = 0;



struct RNAMemory * process_single_file(char *pdbFile)
/* processing a single PDB (or CIF or RNAML) file */
{
    long i, j,  base_all;
    long type_stat[20]; /* maxmum 20 different pairs */
    long **pair_stat; /* maxmum 20 different pairs */
	struct RNAMemory *r=0;

    pair_stat = lmatrix(0, 20, 0, 40);
    for (i = 0; i < 20; i++){
        type_stat[i] =0;        
    }
      
    for (i = 0; i <20; i++) 
        for (j = 0; j <40; j++)
            pair_stat[i][j]=0;
    
    

    r=rna(pdbFile, type_stat, pair_stat, &base_all);
    free_lmatrix(pair_stat,0, 20, 0, 40);
	//chdir(dirLocation);
	return r;
    
}

    
struct RNAMemory * rna(char *pdbfile, long *type_stat, long **pair_stat, long *bs_all)
/* do all sorts of calculations */
{
    char HB_ATOM[BUF512], ALT_LIST[BUF512];
    char *ChainID, *bseq, **AtomName, **ResName, **Miscs;
    long i, j, k,m,n, ie, ib, dna_rna, num, num_residue, nres, bs_atoms;
    long *ResSeq, *RY, **seidx, num_modify, *modify_idx, nprot_atom=0;
    long **chain_idx,nchain;
    double HB_UPPER[2], **xyz;
    static long base_all;    
	struct RNAMemory *r = 0;
    
     
/* read in H-bond length upper limit etc from <misc_rna.par> */    
    hb_crt_alt(HB_UPPER, HB_ATOM, ALT_LIST);

/* read in the PDB file */
    num = number_of_atoms(pdbfile);
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    
    num = read_pdb(pdbfile,AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
                   ALT_LIST);

/* get the numbering information of each residue.
   seidx[i][j]; i = 1-num_residue  j=1,2
*/
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);
      
            
        
/* Below is only for nucleic acids ie RY >= 0*/  
    bs_atoms = 0;
    for(i = 1; i <= num_residue; i++){
        ib = seidx[i][1];
        ie = seidx[i][2];
        dna_rna = residue_ident(AtomName, xyz, ib, ie);
       
        if (dna_rna >= 0){         
            for(j = seidx[i][1]; j <= seidx[i][2]; j++){
                    //if(!strchr(user_chain, toupper(ChainID[j]) ) )continue;
                    if(toupper(ChainID[seidx[1][1]])!=toupper(ChainID[j])) continue;                    
                    bs_atoms++;
                    strcpy(AtomName[bs_atoms], AtomName[j]);
                    strcpy(ResName[bs_atoms], ResName[j]);
                    ChainID[bs_atoms] = ChainID[j];
                    ResSeq[bs_atoms] = ResSeq[j];
                    for(k = 0 ; k <=NMISC; k++)
                        Miscs[bs_atoms][k] = Miscs[j][k];
                    for(k = 1 ; k <=3; k++)
                        xyz[bs_atoms][k] = xyz[j][k];
                
                
            }
        }
    }

    bseq = cvector(1, num_residue);
    nres=0;    
    seidx=residue_idx(bs_atoms, ResSeq, Miscs, ChainID, ResName, &nres);

    chain_idx = lmatrix(1,500 , 1, 2);  /* # of chains max = 200 */    
    get_chain_idx(nres, seidx, ChainID, &nchain, chain_idx);

    bs_atoms = 0;
    for (i=1; i<=nchain; i++){ /* rid of ligand */
        if((chain_idx[i][2] - chain_idx[i][1]) <= 0)continue;
        ib=chain_idx[i][1];
        ie=chain_idx[i][2];
        
        for (k=chain_idx[i][1]; k<=chain_idx[i][2]; k++){
            ib = seidx[k][1];
            ie = seidx[k][2];
            dna_rna = residue_ident(AtomName, xyz, ib, ie);
       
            if (dna_rna >= 0){          
                for(j = ib; j <= ie; j++){
                    bs_atoms++;
                    strcpy(AtomName[bs_atoms], AtomName[j]);
                    strcpy(ResName[bs_atoms], ResName[j]);
                    ChainID[bs_atoms] = ChainID[j];
                    ResSeq[bs_atoms] = ResSeq[j];
                    for(m = 0 ; m <=NMISC; m++)
                        Miscs[bs_atoms][m] = Miscs[j][m];
                    for(m = 1 ; m <=3; m++)
                        xyz[bs_atoms][m] = xyz[j][m];
                    n=bs_atoms;                    
                }
            }
            
        }
        
    }
    
    
    nres=0;    
    seidx=residue_idx(bs_atoms, ResSeq, Miscs, ChainID, ResName, &nres);

    RY = lvector(1, num_residue);
    modify_idx = lvector(1, num_residue);    
    get_seq(nres, seidx, AtomName, ResName, ChainID, ResSeq, Miscs,
            xyz, bseq, RY, &num_modify,modify_idx);  /* get the new RY */

    r=work_horse(pdbfile, nres, bs_atoms, bseq, seidx, RY, AtomName,
               ResName, ChainID, ResSeq, Miscs, xyz,num_modify, modify_idx,  
               type_stat, pair_stat);

    base_all=base_all+nres; /* acculate all the bases */
    *bs_all=base_all;
     
    
    free_cmatrix(AtomName, 1, num, 0, 4);
    free_cmatrix(ResName, 1, num, 0, 3);
    free_cvector(ChainID, 1, num);
    free_lvector(ResSeq, 1, num);
    free_dmatrix(xyz, 1, num, 1, 3);
    free_cmatrix(Miscs, 1, num, 0, NMISC);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lmatrix(chain_idx, 1,500 , 1, 2);   
    free_lvector(RY, 1, num_residue);
    free_lvector(modify_idx, 1, num_residue);

	return r;
    
}

struct RNAMemory * work_horse(char *pdbfile, long num_residue, long num,
                char *bseq, long **seidx, long *RY, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq,char **Miscs, 
                double **xyz,long num_modify, long *modify_idx, 
                long *type_stat,long **pair_stat)
/* perform all the calculations */

{    
    
    long **bs_pairs_tot, num_pair_tot=0, num_single_base=0,*single_base, ntot;
    long i, j;
    long num_bp = 0, num_helix = 1, nout = 16, nout_p1 = 17;
    long pair_istat[17], pair_jstat[17];
    long *bp_idx, *helix_marker, *matched_idx;
    long **base_pairs, **helix_idx;
    long **bp_order, **end_list;
    long num_multi, *multi_idx, **multi_pair;
    long num_bp_best=0, **pair_num_best,*sugar_syn;
    
    double BPRS[7];
    double **orien, **org, **Nxyz, **o3_p, **bp_xyz;
    char **pair_type;
	struct RNAMemory *r=0;

	r = malloc(sizeof(struct RNAMemory));
	r->profile = 0;
	r->secStruct = 0;
	r->sequence = 0;

    multi_pair = lmatrix(1, num_residue, 1, 20); /* max 20-poles */
    multi_idx = lvector(1, num_residue);  /*max multipoles = num_residue */
    pair_type = cmatrix(1, num_residue*2, 0, 3); /* max base pairs */    
    bs_pairs_tot = lmatrix(1, 2*num_residue, 1, 2);
    single_base = lvector(1, num_residue);  /* max single base */   
    sugar_syn = lvector(1, num_residue); 
 
    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    Nxyz = dmatrix(1, num_residue, 1, 3);     /* RN9/YN1 atomic coordinates */
    o3_p = dmatrix(1, num_residue, 1, 8);     /* O3'/P atomic coordinates */

/* get the  base information for locating possible pairs later
  orien --> rotation matrix for each library base to match each residue.
  org --> the fitted sxyz origin for each library base.
*/
    base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
              ResSeq, Miscs, xyz,  orien, org, Nxyz, o3_p, BPRS);    
/* find all the base-pairs */
    r->profile=all_pairs(pdbfile, num_residue, RY, Nxyz, orien, org, BPRS,
              seidx, xyz, AtomName, ResName, ChainID, ResSeq, Miscs, bseq,
              &num_pair_tot, pair_type, bs_pairs_tot, &num_single_base,
              single_base, &num_multi, multi_idx, multi_pair, sugar_syn);
	r->sequence = malloc(sizeof(char)*(num_residue + 1));
	for (i = 0; i < num_residue; i++)
		r->sequence[i] = bseq[i + 1];
	r->sequence[num_residue] = '\0';
    
    
    if (!num_pair_tot) {   
		r->secStruct = malloc(sizeof(char)*(num_residue + 1));
		for (i = 0; i < num_residue; i++)
			r->secStruct[i] = 'S';
		r->secStruct[num_residue] = '\0';

        return r; 
    }

    ntot = 2*num_pair_tot;
/*    if(ARGC <=2){ */

    matched_idx = lvector(1, num_residue);
    base_pairs = lmatrix(1, num_residue, 1, nout_p1);

    for (i = 1; i <= num_residue; i++) {
        best_pair(i, num_residue, RY, seidx, xyz, Nxyz, matched_idx, orien,
                  org, AtomName, bseq, BPRS, pair_istat);
        if (pair_istat[1]) {        /* with paired base */
            best_pair(pair_istat[1], num_residue, RY, seidx, xyz, Nxyz,
                      matched_idx, orien, org, AtomName, bseq, BPRS,
                      pair_jstat);
            if (i == pair_jstat[1]) { /* best match between i && pair_istat[1] */
                matched_idx[i] = 1;
                matched_idx[pair_istat[1]] = 1;
                base_pairs[++num_bp][1] = i;
                for (j = 1; j <= nout; j++)
                    base_pairs[num_bp][j + 1] = pair_istat[j];
            }
        }
    }
    
    bp_idx = lvector(1, num_bp);
    helix_marker = lvector(1, num_bp);
    helix_idx = lmatrix(1, num_bp, 1, 7);
   
    re_ordering(num_bp, base_pairs, bp_idx, helix_marker, helix_idx, BPRS,
                &num_helix, o3_p, bseq, seidx, ResName, ChainID, ResSeq,
                Miscs);

    bp_order = lmatrix(1, num_bp, 1, 3);
    end_list = lmatrix(1, num_bp, 1, 3);
    bp_xyz = dmatrix(1, num_bp, 1, 9); /*   bp origin + base I/II normals: 9 - 17 */
    
    for (i = 1; i <= num_bp; i++)
        for (j = 1; j <= 9; j++)
            bp_xyz[i][j] = base_pairs[i][j + 8] / MFACTOR;
    
    pair_num_best = lmatrix(1, 3, 1, num_bp);
    r->secStruct=write_best_pairs(num_residue,num_helix, helix_idx, bp_idx, helix_marker,
                     base_pairs, seidx, ResName, ChainID, ResSeq,
                     Miscs, bseq, BPRS, &num_bp_best, pair_num_best);    

	fflush(stdout);
    
    
    free_lmatrix(multi_pair, 1, num_residue, 1, 20);  
    free_lvector(multi_idx , 1, num_residue);       
    free_cmatrix(pair_type , 1, num_residue*2, 0, 3);          
    free_lmatrix(bs_pairs_tot , 1, 2*num_residue, 1, 2);     
    free_lvector(single_base , 1, num_residue);         
    free_lvector(sugar_syn , 1, num_residue); 
    free_dmatrix(bp_xyz, 1, num_bp, 1, 9);
    free_lmatrix(pair_num_best, 1, 3, 1, num_bp);

    free_lmatrix(bp_order, 1, num_bp, 1, 3);
    free_lmatrix(end_list, 1, num_bp, 1, 3);

    free_lvector(bp_idx, 1, num_bp);
    free_lvector(helix_marker, 1, num_bp);
    free_lmatrix(helix_idx, 1, num_bp, 1, 7);

    free_lvector(matched_idx, 1, num_residue);
    free_lmatrix(base_pairs, 1, num_residue, 1, nout_p1);


    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);
    free_dmatrix(Nxyz, 1, num_residue, 1, 3);
    free_dmatrix(o3_p, 1, num_residue, 1, 8);

	return r;
}
#ifdef SYSTEM_WINDOWS
extern __declspec(dllexport) int __stdcall GetError()
#else
	int GetError()
#endif
	{
		return 0;// errorNum;
	}


#ifdef SYSTEM_WINDOWS
	extern __declspec(dllexport) struct RNAMemory * __stdcall RNAProfiles(char *fileName)
#else
	struct RNAMemory * RNAProfiles(char *fileName)
#endif
	{
		return process_single_file(fileName);		
	}

#ifdef SYSTEM_WINDOWS
	extern __declspec(dllexport) LPSTR __stdcall GetSS(struct RNAMemory *r)
#else
	LPSTR GetSS(struct RNAMemory *r)
#endif
	{
		return r->secStruct;
	}
#ifdef SYSTEM_WINDOWS
	extern __declspec(dllexport) LPSTR __stdcall GetLW(struct RNAMemory *r)
#else
	LPSTR GetLW(struct RNAMemory *r)
#endif
		{
			return r->profile;
		}

#ifdef SYSTEM_WINDOWS
		extern __declspec(dllexport) LPSTR __stdcall GetSEQ(struct RNAMemory *r)
#else
		LPSTR GetSEQ(struct RNAMemory *r)
#endif
		{
			return r->sequence;
		}
#ifdef SYSTEM_WINDOWS
extern  __declspec(dllexport) void __stdcall RNAFree(struct RNAMemory *r)
#else
		void RNAFree(struct RNAMemory *r)
#endif
	{
			if (r != 0)
			{

				free(r->profile);
				free(r->sequence);
				free(r->secStruct);
				free(r);
			}
	}

