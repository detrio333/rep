#ifndef HAC_CLUTTER_H
#define HAC_CLUTTER_H

#include "../../GISDataBase/headers/det_hydro.h"

#include "main.h"
#include "bands.h"
#include "field_clutter.h"
#include "bellhop_hydrology.h"
#include "bellhop.h"

HAC_LIBRARY_NAMESPACE
{
	struct HAC_EXPORT clutter	//! структуры распределенной помехи
	{
		struct calc_param	//! структура для задания условий расчета самостоятельно
		{
			// производится расчёт поля сигнала на частоте 1000 Гц
			
			enum class clutter_types	//! сложность рассчитанной помехи (влияет на скорость расчетов)
			{
				hard,	//!< высокая
				normal,	//!< нормальная
				low,	//!< низкая (вполне хватает низкой, различия заметны при существенной неровности дна и сложной гидрологии)
			} clutter_type = clutter_types::hard;	//!< сложность рассчитанной помехи

			enum class attenuation_types	//! затухание звука
			{
				none,		//!< не считается
				thorpe,		//!< по формуле торпа
				database,	//!< значение приходит из базы
			}attenuation_type = attenuation_types::thorpe;	//!< затухание звука
		};

		// берутся из беллхопа и норм
		bellhop_hydrology *local_hydrology = NULL;	//!< местная гидрология

		error const set_param(calc_param const&);		//!< запрос данных, установка начальных параметров и параметров по умолчанию

		/*!
		\param[in] _depth_receiver глубина приёмника, м
		\param[out] _field_clutter выходная структура с описанием распределенной помехи
		*/
		error const get_clutter(
			const float &_depth_receiver,		// глубина приёмника, м
			field_clutter &_field_clutter		// помеха
			)const;		//!< функция расчёта распределенной помехи для произвольного дна	


	private:
		calc_param local_calc_param;	//!< параметры расчетов загнаны в одну структуру
	};
};

#endif /*HAC_CLUTTER_H*/
