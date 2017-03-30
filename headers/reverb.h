#ifndef HAC_REVERB_H
#define HAC_REVERB_H

#include "../headers/main.h"
#include "../headers/bands.h"
#include "../headers/field_reverb.h"
#include "../headers/field_clutter.h"
#include "../headers/ray.h"
#include "../headers/bellhop_hydrology.h"
#include "../../GISDataBase/headers/det_hydro.h"

HAC_LIBRARY_NAMESPACE
{
	struct HAC_EXPORT reverb	//! реверберация
	{
		// при расчетах используем предположение, что углов много и угловое расстояние между ними строго меньше ширины ХН антенны, а сама антенна близка к квадратной
		enum class type_reverb	//! тип реверберации
		{
			surface			//!< поверхностная
			, volume	//!< объемная
			, seabed	//!< донная
		};		

		struct calc_param	//! структура для задания условий расчета самостоятельно
		{
			float t_sign;	//!< длительность зонлирующего сигнала, с 

			float t_lim;	//!< предельное время расчета реверберации, с

			std::list<type_reverb> type;	//!< тип рассчитываемой реверьерации (допускает расчёт нескольких типов одновременно)
		};

		bellhop_hydrology *local_hydrology = NULL;	//!< местная гидрология

		error const set_param(calc_param const&);		//!< запрос данных, установка начальных параметров и параметров по умолчанию
		error const get_clutter(
			field_reverb&, const std::map<float, std::vector<step_ray>> &
			)const;		//!< функция расчёта реверберации

		std::map<float, float> sum_spectrum(field_reverb&, float const &t);	//!< демонстрационная функция для ненаправленного приемника и излучателя, t - интервал суммирования
		// время, мощность
	private:
		calc_param local_calc_param;
 		bool surface_r;
 		bool volume_r;
 		bool seabed_r;
	};
};

#endif /*HAC_REVERB_H*/
