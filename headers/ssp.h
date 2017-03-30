#ifndef HAC_SSP_H
#define HAC_SSP_H

#include "../headers/main.h"
#include "../headers/coord_2d.h"

#include "../../GISDataBase/headers/det_hydro.h"

#include <map>

HAC_LIBRARY_NAMESPACE
{
	struct sonic	//! элемент среза ВРСЗ
	{
		float r;	//!< удаленность от источника, м
		float h;	//!< глубины горизонтов, м
		float c;	//!< скорость звука, м/с
		float cz;	//!< производная скорости звука по вертикальной координате
	};

	struct ss_at_point	//! структура скорости звука и его градиента в точке пространства
	{
		float c;	//!< скорость звука, м/с
		coord_2d grad_c;	//!< градиент скорости звука
		float crr;	//!< вторая производная скорости звука по горизонтальной координате
		float crz;	//!< производная скорости звука по горизонтальной координате и глубине
		float czz;	//!< вторая производная скорости звука по глубине

//		unsigned int Layer;
	};

	struct HAC_EXPORT sspmod	//!< описание ВРСЗ
	{
		sspmod();

		enum dist_numbers : char
		{
			one,
			some
		} dist_number;

		std::map<float, sonic> depth_ssp;	// 1 - глубина
		std::map<float, sonic>::iterator Layer;
		std::map<float, std::map<float, sonic>> range_depth_ssp;	// 1 - гор. коор, 2 - верт. коорд
		std::map<float, std::map<float, sonic>>::const_iterator iter_range_depth_ssp_low;
		std::map<float, std::map<float, sonic>>::const_iterator iter_range_depth_ssp_up;
 		sonic sonic_l;
		sonic sonic_u;
		float z_up, z_low;

		error const set_ssp(std::map<float, GDB::det::SSP> const &ssp_slices);	// заполнение данных 
		ss_at_point get_ssp(const float &_r, const float &_z);	//!< получение данных в данной точке
	};
}

#endif /*HAC_SSP_H*/
