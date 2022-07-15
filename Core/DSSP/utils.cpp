// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
#include "stdafx.h"
#define BOOST_USE_WINDOWS_H

#include "mas.h"

#include <iostream>
#include <cstdio>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
//#include <boost/timer/timer.hpp>


#include "utils.h"


using namespace std;

// --------------------------------------------------------------------

/*arg_vector::operator char* const*()
{
	m_argv.clear();
	foreach (string& s, m_args)
	{
		m_argv.push_back(s.c_str());
		if (VERBOSE > 1)
			cerr << m_argv.back() << ' ';
	}
	if (VERBOSE > 1)
		cerr << endl;

	m_argv.push_back(nullptr);
	return const_cast<char*const*>(&m_argv[0]);
}

ostream& operator<<(ostream& os, const arg_vector& argv)
{
	os << "About to execute: " << endl;
	foreach (const string& a, argv.m_args)
		os << a << ' ';
	os << endl;

	return os;
}*/

// --------------------------------------------------------------------

mas_exception::mas_exception(const string& msg)
{
	snprintf(m_msg, sizeof(m_msg), "%s", msg.c_str());
}


//// --------------------------------------------------------------------
//
//string decode(const sequence& s)
//{
//	string result;
//	result.reserve(s.length());
//	
//	foreach (aa a, s)
//		result.push_back(kAA[a]);
//
//	return result;
//}
//
//namespace {
//
//bool sInited = false;
//uint8 kAA_Reverse[256];
//
//inline void init_reverse()
//{
//	if (not sInited)
//	{
//		// init global reverse mapping
//		for (uint32 a = 0; a < 256; ++a)
//			kAA_Reverse[a] = 255;
//		for (uint8 a = 0; a < sizeof(kAA); ++a)
//		{
//			kAA_Reverse[toupper(kAA[a])] = a;
//			kAA_Reverse[tolower(kAA[a])] = a;
//		}
//	}
//}
//
//}
//
//aa encode(char r)
//{
//	init_reverse();
//
//	if (r == '.' or r == '*' or r == '~' or r == '_')
//		r = '-';
//	
//	aa result = kAA_Reverse[static_cast<uint8>(r)];
//	if (result >= sizeof(kAA))
//		throw mas_exception(boost::format("invalid residue %1%") % r);
//	
//	return result;
//}
//
//sequence encode(const string& s)
//{
//	init_reverse();
//	
//	sequence result;
//	result.reserve(s.length());
//
//	foreach (char r, s)
//	{
//		if (r == '.' or r == '*' or r == '~' or r == '_')
//			r = '-';
//		
//		aa rc = kAA_Reverse[static_cast<uint8>(r)];
//		if (rc >= sizeof(kAA))
//			throw mas_exception(boost::format("invalid residue in sequence %1%") % r);
//
//		result.push_back(rc);
//	}
//	
//	return result;
//}

// --------------------------------------------------------------------

/*#ifndef NDEBUG
stats::~stats()
{
	if (VERBOSE) 
		cerr << endl << "max: " << m_max << " count: " << m_count << " average: " << (m_cumm / m_count) << endl;
}
#endif*/

// --------------------------------------------------------------------

