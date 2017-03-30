#include "../headers/field_gain.h"

HAC_USING_LIBRARY_NAMESPACE;

std::ostream HAC_EXPORT &operator << (std::ostream& _stream, const field_gain &_field_gain)
{
	if (_field_gain.rays.empty())
	{
		return _stream;
	}

	for (auto iter = _field_gain.rays.begin(); iter != _field_gain.rays.end(); iter++)
	{
		_stream << iter->second.angle_rcv << " " << iter->first << "\n";
	}

	return _stream;
}