/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_STRING_H
#define MECHSYS_STRING_H

// STL
#include <string>
#include <cstdarg> // for va_list, va_start, va_end
#include <cstdlib> // for vsnprintf

#define _(STRING) (STRING) // for translation

class String : public std::string
{
public:
	String(){}
	String(std::string const & Other)
	{
		this->clear();
		this->append(Other);
	}
	String(const char * Msg)
	{
		this->clear();
		this->append(Msg);
	}
	int Printf(String const & Fmt, ...)
	{
		int len;
		va_list       arg_list;
		va_start     (arg_list, Fmt);
		len=_set_msg (Fmt.c_str(), arg_list);
		va_end       (arg_list);
		return len;
	}
	int PrintfV(String const & Fmt, va_list ArgList)
	{
		return _set_msg (Fmt.c_str(), ArgList);
	}
	char const * CStr() const { return this->c_str(); }
private:
	int _set_msg(char const * Fmt, va_list ArgList)
	{
		const int size = 4048; // TODO: remove this limitation by using a loop and reallocating space
		char      buffer[size];
		int       len = std::vsnprintf(buffer, size, Fmt, ArgList);
		this->clear();
		if (len<0) this->append("String::_set_msg: INTERNAL ERROR: std::vsnprintf FAILED");
		else
		{
			buffer[len]='\0';
			if (len>size) this->append("String::_set_msg: INTERNAL ERROR: std::vsnprintf MESSAGE TRUNCATED: ");
			this->append(buffer);
		}
		return len;
	}
}; // class String

#endif // MECHSYS_STRING_H
