#ifndef HAC_FIELD_CLUTTER_H
#define HAC_FIELD_CLUTTER_H

#include "main.h"
#include "bands.h"

#include <map>

HAC_LIBRARY_NAMESPACE
{
	struct HAC_EXPORT field_clutter	//!	выходная структура распределённой помехи
	{
		float  P;					//!< приведенная мощность шумов моря на входе антенны, Па^2 в полосе 1 Гц на частоте 1 кГц  
		std::map<float, float> Q;	//!< нормированное (по всем лучам) распределение интенсивности шумов моря по верт. углу на входе антенны, ед.(первый индекс - угол скольжения луча)

		float const get_pwr(space_band_v const&)const;	//!< возвращет мощность помехи в заданном пространственном угле, Па^2 в полосе 1 Гц на частоте 1 кГц 
	};
};

#endif /*HAC_FIELD_CLUTTER_H*/
