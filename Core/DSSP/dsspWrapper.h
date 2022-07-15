#pragma once

using namespace System;
using namespace std;

namespace dsspWrapper 
{

	public ref class Wrapper
	{
		
	public:
		System::String^ SEQ;
		System::String^ SS;
		System::String^ SA;
		System::String^ Errors;
		double timeSp;
		
	public:
		
		void Run(System::String^ fileName,int len);
	};
}
