#ifndef HAC_FIELD_REVERB_H
#define HAC_FIELD_REVERB_H

#include "main.h"
#include "bands.h"

#include <map>

HAC_LIBRARY_NAMESPACE
{
	struct HAC_EXPORT field_reverb	//!	выходная структура реверберационной помехи
	{
		struct field_reverb_slice	//! реализация реверберации во времени
		{
		private:
			std::multimap<float, float> Q_time;	//!< массив (по времени, c) интенсивностей реверберационной помехи в данном направлении, в единицах от мощности излучения в данном направлении.
						// время  относительная мощность рассеяного излучения(относительно мощности сигнала)
		public:
			float const get_pwr(
				time_band const&
				)const;	//!< возвращет относительная мощность реверберационной помехи в заданном временном интервале, ед;

			friend struct reverb;
		};

		float const get_pwr(
			time_band const&	//!< временной интервал для которого вычисляется усреднённая величина
			, space_band_v const&	//!< пространственный угол в ВП
			)const;	//!< возвращет относительная мощность реверберационной помехи в заданном пространственном угле, ед;

	private:
		
		std::map<float, field_reverb_slice> Q_angle_time;	//!< изменение поле реверберации по простраснтву (вертикальным углам)

		friend struct reverb;
	};
};

#endif /*HAC_FIELD_REVERB_H*/
