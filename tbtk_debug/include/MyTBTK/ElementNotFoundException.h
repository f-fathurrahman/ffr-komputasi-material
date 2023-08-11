#ifndef COM_DAFER45_MyTBTK_ELEMENT_NOT_FOUND_EXCEPTION
#define COM_DAFER45_MyTBTK_ELEMENT_NOT_FOUND_EXCEPTION

#include "MyTBTK/Exception.h"

#include <string>

namespace MyTBTK{

class ElementNotFoundException : public Exception{
public:
	/** Constructor. */
	ElementNotFoundException();

	/** Constructor. */
	ElementNotFoundException(
		const std::string& function,
		const std::string& where,
		const std::string& message,
		const std::string& hint
	);

	/** Destructor. */
	virtual ~ElementNotFoundException();
private:
};

};	//End of namespace MyTBTK

#endif
