#ifndef STOLEN_DATETIME_FUNCTION
#define STOLEN_DATETIME_FUNCTION
#include<ctime>
const std::string currentDateTime(){
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
	return buf;
};
#endif//STOLEN_DATETIME_FUNCTION
