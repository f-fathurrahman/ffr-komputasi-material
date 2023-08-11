#ifndef COM_DAFER45_MyTBTK_INDEX_EXCEPTION
#define COM_DAFER45_MyTBTK_INDEX_EXCEPTION

#include "MyTBTK/Exception.h"

#include <string>

namespace MyTBTK{

class IndexException : public Exception{
public:
	/** Constructor. */
	IndexException();

	/** Constructor. */
	IndexException(
		const std::string& function,
		const std::string& where,
		const std::string& message,
		const std::string& hint
	);

	/** Destructor. */
	virtual ~IndexException();
private:
};

};	//End of namespace MyTBTK

#endif
