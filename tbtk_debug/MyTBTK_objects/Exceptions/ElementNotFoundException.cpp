#include "MyTBTK/ElementNotFoundException.h"

namespace MyTBTK{

ElementNotFoundException::ElementNotFoundException(){
}

ElementNotFoundException::ElementNotFoundException(
	const std::string& function,
	const std::string& where,
	const std::string& message,
	const std::string& hint
) : Exception(function, where, message, hint){
}

ElementNotFoundException::~ElementNotFoundException(){
}

};	//End of namespace Ygg
