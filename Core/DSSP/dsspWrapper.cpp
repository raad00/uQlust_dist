#include "stdafx.h"


#include <iostream>
#include <fstream>
#include <time.h>
#include "mas.h"
#include "structure.h"
#include "dsspWrapper.h"


using namespace std;
using namespace System;
using namespace System::Runtime::InteropServices;

namespace dsspWrapper
{
	void Wrapper::Run(String^ fileName,int len)
	{

		char buffor[100];
		struct dirent *ent;
		string sS;
		string sA;
		string seQ;
		string psi;
		string phi;
		for(int i=0;i<len;i++)
			buffor[i]=fileName[i];

		buffor[len]='\0';
		ifstream infile(buffor, ios_base::in | ios_base::binary);
			//ifstream infile(fileName, ios_base::in | ios_base::binary);

		if (!infile.is_open())
			return;
		vector<MChain*> chain;
		try
		{
			MProtein a(infile);
		
			clock_t start,end;
			start=clock();
			a.CalculateSecondaryStructure();
			end=clock();
			chain=a.GetChains();
		
			timeSp+=(end-start)/(CLOCKS_PER_SEC/1000);
		//Only first chain is considered
		for(int i=0;i<1;i++)
		{
			vector <MResidue *> residue=chain[i]->GetResidues();
			for(int j=0;j<residue.size();j++)
			{
				char ss;
				switch (residue[j]->GetSecondaryStructure())
				{
					case alphahelix:	ss = 'H'; break;
					case betabridge:	ss = 'B'; break;
					case strand:		ss = 'E'; break;
					case helix_3:		ss = 'G'; break;
					case helix_5:		ss = 'I'; break;
					case turn:			ss = 'T'; break;
					case bend:			ss = 'S'; break;
					case loop:			ss = 'C'; break;
				}
				char aa=kResidueInfo[residue[j]->GetType()].code;
				seQ+=aa;
				sS+=ss;
				sprintf(buffor,"%d ",(int)floor(residue[j]->Accessibility() + 0.5));
				sA+=buffor;				

			}
		
			//cout<<SS;
		}
		infile.close();
		const char *cS=sS.c_str();
		SS=gcnew String(cS);
		const char *cA=sA.c_str();
		SA=gcnew String(cA);
		const char *cSEQ=seQ.c_str();
		SEQ=gcnew String(cSEQ);
		}
		catch(exception &e)
		{
			const char *er=e.what();
			Errors=gcnew String(er);


		}

	}
}
